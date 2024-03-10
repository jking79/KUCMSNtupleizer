// -*- C++ -*-
//
// Package:    KUNTupleFW/KUCMSNtupleizer
// Class:      ECALTracksProducer
//
/**\class ECALTracksProducer ECALTracksProducer.cc KUNTupleFW/KUCMSNtupleizer/plugins/ECALTracksProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Sat, 8 March 2024 10:55:00 GMT
//
//

// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/TrackExtrapolation.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/GeomPropagators/interface/HelixExtrapolatorToLine2Order.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

//Local includes
#include "KUNTupleFW/KUCMSNtupleizer/interface/Hungarian.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/MatchTracksToSC.h"

//
// class declaration
//

bool IsIndexMatched(const std::vector<int>& matchedIndexes, int index);
//template<typename T>
//std::vector<T> EraseIndices(const std::vector<T>& data, std::vector<size_t>& indicesToDelete);

typedef ROOT::Math::PtEtaPhiMVector LorentzVec;

class ECALTracksProducer : public edm::stream::EDProducer<> {
public:
  explicit ECALTracksProducer(const edm::ParameterSet&);
  ~ECALTracksProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef edm::Handle<reco::TrackCollection> TracksHandle;
  typedef edm::Handle<reco::GsfTrackCollection> GsfTracksHandle;
  typedef edm::Handle<reco::SuperClusterCollection> SCsHandle;

  template <typename T> using PropagatedTrack = std::pair<T, TrackDetMatchInfo>;
  template <typename T> using PropagatedTracks = std::vector<PropagatedTrack<T> >;
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  template <typename T, typename S>
  std::vector<T> RemoveOverlap(const std::vector<T> &toReduce, const std::vector<S> &toMatch, const double &deltaRcut); 

  template <typename T>
  PropagatedTracks<T> GetAcceptedTracks(const edm::Event &iEvent, const edm::EventSetup &iSetup,
					std::vector<T> &tracks, const MagneticField &magneticField); 

  TrackDetMatchInfo GetTrackDetMatchInfo(const edm::Event &iEvent, const edm::EventSetup &iSetup,
					 const reco::Track& track, const MagneticField &magneticField);

  // ----------member data ---------------------------
  edm::EDGetTokenT<reco::TrackCollection> generalTrackToken_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfElectronTrackToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> superClusterToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> ootSuperClusterToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  reco::GenParticleCollection signalGenElectrons_;
};

ECALTracksProducer::ECALTracksProducer(const edm::ParameterSet& iConfig) :
  generalTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTracksSrc")) ),
  gsfElectronTrackToken_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronTracksSrc")) ),
  superClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters")) ),
  ootSuperClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters")) ),
  caloGeometryToken_(esConsumes()),
  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()) 
{
  // TrackAssociator parameters
  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(parameters, iC);
  trackAssociator_.useDefaultPropagator();

  produces<reco::TrackCollection>("ecalGeneralTracks").setBranchAlias("ecalGeneralTracks");
  produces<reco::GsfTrackCollection>("ecalGsfTracks").setBranchAlias("ecalGsfTracks");
  //produces<PropagatedTracks<reco::Track> >("propagatedGeneralTracks").setBranchAlias("propagatedGeneralTracks");
  //produces<PropagatedTracks<reco::GsfTrack> >("propagatedGsfTracks").setBranchAlias("propagatedGsfTracks");
  produces<reco::SuperClusterCollection>("displacedElectronSCs").setBranchAlias("displacedElectronSCs");
}

ECALTracksProducer::~ECALTracksProducer() {}

//
// member functions
//
// ------------ method called to produce the data  ------------
void ECALTracksProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;

  typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<double> > Point;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<double> > Vector;

  auto magfield = iSetup.getTransientHandle(magneticFieldToken_);
  CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  edm::Handle<reco::TrackCollection> generalTracksHandle;
  iEvent.getByToken(generalTrackToken_, generalTracksHandle);

  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle;
  iEvent.getByToken(gsfElectronTrackToken_, gsfTracksHandle);

  edm::Handle<reco::SuperClusterCollection> superClusterHandle;
  iEvent.getByToken(superClusterToken_, superClusterHandle);

  edm::Handle<reco::SuperClusterCollection> ootSuperClusterHandle;
  iEvent.getByToken(ootSuperClusterToken_, ootSuperClusterHandle);

  // Smart pointers to containers of output collections
  std::unique_ptr<reco::TrackCollection> ecalGeneralTracks = std::make_unique<reco::TrackCollection>();
  std::unique_ptr<reco::GsfTrackCollection> ecalGsfTracks = std::make_unique<reco::GsfTrackCollection>();
  std::unique_ptr<PropagatedTracks<reco::Track> > propagatedGeneralTracks = std::make_unique<PropagatedTracks<reco::Track> >();
  std::unique_ptr<PropagatedTracks<reco::GsfTrack> > propagatedGsfTracks = std::make_unique<PropagatedTracks<reco::GsfTrack> >();
  std::unique_ptr<reco::SuperClusterCollection> displacedElectronSCs = std::make_unique<reco::SuperClusterCollection>();
  
  // Remove overlap between general tracks and gsf tracks. Gsf tracks are preferred.
  reco::TrackCollection generalTracks(*generalTracksHandle);
  reco::GsfTrackCollection gsfTracks(*gsfTracksHandle);

  // Retain only General tracks that make it to the ECAL
  generalTracks = RemoveOverlap<reco::Track, reco::GsfTrack>(generalTracks, gsfTracks, 0.1);
  PropagatedTracks<reco::Track> propGeneralTracks = GetAcceptedTracks(iEvent, iSetup, generalTracks, *(magfield.product()) );
  ecalGeneralTracks->insert(ecalGeneralTracks->end(), generalTracks.begin(), generalTracks.end());
  propagatedGeneralTracks->insert(propagatedGeneralTracks->end(), propGeneralTracks.begin(), propGeneralTracks.end());

  // Retain only Gsf tracks that make it to the ECAL
  PropagatedTracks<reco::GsfTrack> propGsfTracks = GetAcceptedTracks(iEvent, iSetup, gsfTracks, *(magfield.product()) );
  ecalGsfTracks->insert(ecalGsfTracks->end(), gsfTracks.begin(), gsfTracks.end());
  propagatedGsfTracks->insert(propagatedGsfTracks->end(), propGsfTracks.begin(), propGsfTracks.end());

  //if(ecalGsfTracks->size() != gsfTracksHandle->size())
  //std::cout << "out of " << gsfTracksHandle->size() << " gsf tracks only " << ecalGsfTracks->size() << " make it to the ECAL!" << std::endl; 

  // Remove overlap between superClusters and out-of-time superClusters
  reco::SuperClusterCollection superClusters(*superClusterHandle);
  reco::SuperClusterCollection ootSuperClusters(*ootSuperClusterHandle);

  superClusters = RemoveOverlap<reco::SuperCluster, reco::SuperCluster>(superClusters, ootSuperClusters, 0.1);
  displacedElectronSCs->insert(displacedElectronSCs->end(), superClusters.begin(), superClusters.end());

  iEvent.put(std::move(ecalGeneralTracks), "ecalGeneralTracks");
  iEvent.put(std::move(ecalGsfTracks), "ecalGsfTracks");
  //iEvent.put(std::move(propagatedGeneralTracks), "propagatedGeneralTracks");
  //iEvent.put(std::move(propagatedGsfTracks), "propagatedGsfTracks");
  iEvent.put(std::move(displacedElectronSCs), "displacedElectronSCs");
}// Producer end

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ECALTracksProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  /*
  desc.add<edm::InputTag>("generalTracksSrc", edm::InputTag("generalTracks"));
  desc.add<edm::InputTag>("gsfElectronTracksSrc", edm::InputTag("electronGsfTracks"));
  desc.add<edm::InputTag>("superClusters", edm::InputTag("particleFlowEGamma"));
  desc.add<edm::InputTag>("ootSuperClusters", edm::InputTag("particleFlowSuperClusterOOTECAL"));
  //desc.add<edm::InputTag>("tkAssocParamBlock");
  */
  descriptions.addDefault(desc);
}

template <typename T, typename S>
std::vector<T> ECALTracksProducer::RemoveOverlap(const std::vector<T> &toReduce, const std::vector<S> &toMatch, const double &deltaRcut) {

  DeltaRMatchHungarian<T, S> assigner(toReduce, toMatch);

  std::vector<size_t> matchedIndexes;

  for(auto const &pair : assigner.GetPairedObjects()) {
    if(pair.GetDeltaR() < deltaRcut) {
      matchedIndexes.push_back(size_t(pair.GetIndexA()));
    }
  }

  return EraseIndices<T>(toReduce, matchedIndexes);
}

template <typename T>
ECALTracksProducer::PropagatedTracks<T> ECALTracksProducer::GetAcceptedTracks(const edm::Event &iEvent,
									      const edm::EventSetup &iSetup,
									      std::vector<T> &tracks, 
									      const MagneticField &magneticField) {

  ECALTracksProducer::PropagatedTracks<T> propagatedTracks;
  std::vector<T> acceptedTracks;
  for(const auto &track : tracks) {
    
    TrackDetMatchInfo detInfo(GetTrackDetMatchInfo(iEvent, iSetup, track, magneticField));
    
    if(detInfo.crossedEcalIds.size() != 0 && track.pt() > 0.95) {
      propagatedTracks.emplace_back(std::make_pair(track, detInfo));
      acceptedTracks.emplace_back(track);
    }
  } 

  tracks = acceptedTracks;

  return propagatedTracks;
}

TrackDetMatchInfo ECALTracksProducer::GetTrackDetMatchInfo(const edm::Event &iEvent,
							   const edm::EventSetup &iSetup,
							   const reco::Track& track, 
							   const MagneticField &magneticField) {

  FreeTrajectoryState initialState = trajectoryStateTransform::initialFreeState(track, &magneticField);
  
  const TrackDetMatchInfo detInfo(trackAssociator_.associate(iEvent, iSetup, trackAssocParameters_, &initialState));

  return detInfo;

}

//define this as a plug-in
DEFINE_FWK_MODULE(ECALTracksProducer);
