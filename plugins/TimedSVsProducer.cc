// -*- C++ -*-
//
// Package:    KUCMSNtupleizer/KUCMSNtupleizer
// Class:      TimedSVsProducer
//
/**\class TimedSVsProducer TimedSVsProducer.cc KUCMSNtupleizer/KUCMSNtupleizer/plugins/TimedSVsProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 21 Oct 2023 20:06:20 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/HelixExtrapolatorToLine2Order.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

//
// class declaration
//

class TimedSVsProducer : public edm::stream::EDProducer<> {
public:
  explicit TimedSVsProducer(const edm::ParameterSet&);
  virtual ~TimedSVsProducer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  reco::Vertex ConvertFitVertex(const TransientVertex &vertex) const;
  reco::Track ConvertTrack(const reco::Track &track, const double t0, const double beta) const;
  double LinearDistance(const GlobalPoint &initialPosition, const GlobalPoint &destination) const;
  double PathLength(const reco::TransientTrack &track, const GlobalPoint &vertexPosition, const GlobalPoint &positionAtECAL) const;

  // ----------member data ---------------------------
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;

  std::string outputLabel_;
  //std::string gsfElectronLabel_, vertexLabel_;
};

TimedSVsProducer::TimedSVsProducer(const edm::ParameterSet& iConfig) {

  transientTrackBuilder_ = esConsumes(edm::ESInputTag("", "TransientTrackBuilder"));

  electronToken_ = consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electronSrc"));
  
  outputLabel_ = "timedSecondaryVertices";
  produces<reco::VertexCollection>(outputLabel_).setBranchAlias(outputLabel_);
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void TimedSVsProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  // Transient track builder for transforming GsfTracks into Transient Tracks
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);

  // Get input electron collection
  Handle<reco::ElectronCollection> electrons;
  iEvent.getByToken(electronToken_, electrons);
  const reco::ElectronCollection electronCollection(*electrons);

  // Smart pointers to containers of output collections
  std::unique_ptr<reco::VertexCollection> timedVertices = std::make_unique<reco::VertexCollection>();
									   
  // Create container for storing transient tracks (converted from reco::Tracks from input electrons)
  std::vector<reco::TransientTrack> transientTrackCollection;

  // Make sure there are enough electrons to construct a vertex from their tracks
  if(electronCollection.size() < 2) {
    iEvent.put(std::move(timedVertices), outputLabel_);
    return;
  }

  // Create map for matching gsfElectrons to the tracks after after the vertex fit.
  // Currently using the track chi2 and the track pT to perform the match, in that order. 
  std::map<std::pair<double, double>, const reco::Electron> trackToElectronMap;

  // Create transient track collection from gsfTracks for each of the input gsfElectrons
  reco::ElectronCollection::const_iterator eleIt;
  for (eleIt = electronCollection.begin(); eleIt != electronCollection.end(); eleIt++){
    const reco::Electron electron = *eleIt;
    const reco::TrackRef trackRef = electron.track();

    if(trackRef.isNull()) 
      continue;

    // Convert the electron gsfTrack to a transient track (needed by the AVR algorithm)
    const reco::TransientTrack track = ttBuilder->build(trackRef);
    
    // Fill map between transient track and electron 
    const double chi2 = track.track().chi2();
    const double pt  = track.track().pt();
    trackToElectronMap.emplace(std::make_pair(chi2, pt), electron);
    
    // Push converted track to the transient track collection
    transientTrackCollection.push_back(track);

  } // End electron loop
  
  // Ensure that there are at least 2 tracks to perform the vertex fit
  if(transientTrackCollection.size() < 2) {
    iEvent.put(std::move(timedVertices), outputLabel_);
    return;
  }

  // Create the vertex collection using the Adaptive Vertex Reconstructor (AVR)
  AdaptiveVertexReconstructor vertexConstructor;
  std::vector<TransientVertex> vertexCollection = vertexConstructor.vertices(transientTrackCollection);

  // Loop over the reconstructed vertices and their respective track collections
  for(auto const &vertex : vertexCollection) {

    if(!vertex.isValid()) 
      continue;

    // Construct a reco::Vertex using the same fields of the reconstructed TransientVertex
    reco::Vertex outputVertex = ConvertFitVertex(vertex);

    // Loop over tracks in current reconstructed vertex 
    for(auto const &track : vertex.originalTracks()) {
      
      const reco::Track recoTrack = track.track();

      // Match current track with its original electron and get hit position at ECAL 
      const reco::Electron trackElectron = trackToElectronMap[std::make_pair(recoTrack.chi2(), recoTrack.pt())];
      const math::XYZPoint trackDestination = trackElectron.superCluster()->position();//trackElectron.trackPositionAtCalo();
      const GlobalPoint trackPositionAtCalo = GlobalPoint(trackDestination.x(), trackDestination.y(), trackDestination.z());
      
      //Calculate path length with helical extrapolation 
      const double pathLength = PathLength(track, vertex.position(), trackPositionAtCalo);

      // Linear distance between the reconstructed vertex and the electron position at the ECAL
      const double linearDistance = LinearDistance(vertex.position(), trackPositionAtCalo);

      // Fill new track with path length and linear distance information
      const reco::Track updatedTrack = ConvertTrack(recoTrack, pathLength, linearDistance);

      // Convert to transient track (only way I found to save track using the reco::Vertex add method)
      const reco::TransientTrack ttrack = ttBuilder->build(updatedTrack);

      // Add tracks to output vertex (tracks are saved as refitted tracks)
      outputVertex.add(ttrack.trackBaseRef(), updatedTrack, 1.0); 

    }// Finish track loop

    // Collect the reconstructed vertices with attached tracks and path lengths
    timedVertices->push_back(outputVertex);

  }// Finish vertex loop

  iEvent.put(std::move(timedVertices), outputLabel_);
}// Producer end

reco::Vertex TimedSVsProducer::ConvertFitVertex(const TransientVertex &vertex) const {

  // Parameters needed for initializing new vertex object
  const math::XYZPoint vertexPosition = math::XYZPoint(vertex.position().x(), vertex.position().y(), vertex.position().z());
  const double vertexChi2 = vertex.totalChiSquared();
  const double vertexNdof = vertex.degreesOfFreedom();
  const size_t nTracks = vertex.originalTracks().size();
  const GlobalError error = vertex.positionError();
  const math::Error<3>::type vertexError(ROOT::Math::SVector<double, 6>(error.cxx(), error.cyx(), error.cyy(),
									error.czx(), error.czy(), error.czz()));

  // Save the new fitted vertex to a reco::Vertex and reserve space for the tracks
  reco::Vertex outputVertex(vertexPosition, vertexError, vertexChi2, vertexNdof, nTracks);
  outputVertex.reserve(nTracks);
  
  return outputVertex;
}

reco::Track TimedSVsProducer::ConvertTrack(const reco::Track &track, const double t0, const double beta) const {

  // Track parameters for creating new track instance with path length info
  const double chi2 = track.chi2();
  const int charge = track.charge();
  const double ndof = track.ndof();
  const math::XYZPoint trackRefPoint = track.referencePoint();
  const math::XYZVector trackMomentum = track.momentum();

  // Create track and save path length to t0 field (also saving the linear distance to the beta parameter)
  const reco::Track updatedTrack(chi2, ndof, trackRefPoint, trackMomentum, charge, track.covariance(),
				 reco::TrackBase::undefAlgorithm, reco::TrackBase::undefQuality,
				 t0, beta);

  return updatedTrack;
}

double TimedSVsProducer::LinearDistance(const GlobalPoint &initialPosition, const GlobalPoint &destination) const {

  double dx = destination.x() - initialPosition.x();
  double dy = destination.y() - initialPosition.y();
  double dz = destination.z() - initialPosition.z();
  double distance = sqrt(dx*dx + dy*dy + dz*dz);

  return distance;
}

double TimedSVsProducer::PathLength(const reco::TransientTrack &track, const GlobalPoint &vertexPosition, const GlobalPoint &positionAtECAL) const { 

  const FreeTrajectoryState initialState = track.initialFreeState();
  const GlobalPoint trackInitialPosition(initialState.position());
  const GlobalVector trackInitialMomentum(initialState.momentum());

  //Calculate path length with helical extrapolation
  //HelixExtrapolatorToLine2Order extrapolator(Basic3DVector<float>(trackInitialPosition), Basic3DVector<float>(trackInitialMomentum),
  //initialState.transverseCurvature(), PropagationDirection::oppositeToMomentum);
  //const double helix_length = extrapolator.pathLength(vertexPosition).second;

  //Define propagator tool for calculating path length between SV and position at Ecal
  SteppingHelixPropagator propagator(track.field());

  //Calculate path length using the stepping helix propagator
  double pathLength = propagator.propagateWithPath(initialState, positionAtECAL).second;// + fabs(helix_length);

  return pathLength;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void TimedSVsProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("electronSrc", edm::InputTag("displacedElectrons", "displacedElectrons"));
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TimedSVsProducer);
