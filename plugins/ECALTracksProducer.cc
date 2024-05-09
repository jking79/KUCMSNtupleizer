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
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackTools.h"


//
// class declaration
//

class ECALTracksProducer : public edm::stream::EDProducer<> {
public:
  explicit ECALTracksProducer(const edm::ParameterSet&);
  ~ECALTracksProducer() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  typedef edm::Handle<reco::TrackCollection> TracksHandle;
  typedef edm::Handle<reco::GsfTrackCollection> GsfTracksHandle;
  typedef edm::Handle<reco::SuperClusterCollection> SCsHandle;

  void produce(edm::Event&, const edm::EventSetup&) override;

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
  produces<reco::TrackCollection>("ecalTracks").setBranchAlias("ecalTracks");
  //produces<std::vector<GlobalPoint> >("generalTracksECALPosition").setBranchAlias("generalTracksECALPosition");
  //produces<std::vector<GlobalPoint> >("gsfTracksECALPosition").setBranchAlias("gsfTracksECALPosition");
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
  std::unique_ptr<reco::TrackCollection> ecalTracks = std::make_unique<reco::TrackCollection>();
  //std::unique_ptr<std::vector<GlobalPoint> > generalTracksECALPosition = std::make_unique<std::vector<GlobalPoint> >();
  //std::unique_ptr<std::vector<GlobalPoint> > gsfTracksECALPosition = std::make_unique<std::vector<GlobalPoint> >();
  std::unique_ptr<reco::SuperClusterCollection> displacedElectronSCs = std::make_unique<reco::SuperClusterCollection>();
  
  // Remove overlap between general tracks and gsf tracks. Gsf tracks are preferred.
  reco::TrackCollection generalTracks(*generalTracksHandle);
  reco::GsfTrackCollection gsfTracks(*gsfTracksHandle);

  // Retain only General tracks that make it to the ECAL
  generalTracks = RemoveOverlap<reco::Track, reco::GsfTrack>(generalTracks, gsfTracks, 0.1);
  TrackPropagator<reco::Track> generalTracksPropagator(iEvent, iSetup, magfield, trackAssocParameters_, generalTracks);
  generalTracks = generalTracksPropagator.GetTracks();
  ecalGeneralTracks->insert(ecalGeneralTracks->end(), generalTracks.begin(), generalTracks.end());

  // Retain only Gsf tracks that make it to the ECAL
  TrackPropagator<reco::GsfTrack> gsfTracksPropagator(iEvent, iSetup, magfield, trackAssocParameters_, gsfTracks);
  gsfTracks = gsfTracksPropagator.GetTracks();
  ecalGsfTracks->insert(ecalGsfTracks->end(), gsfTracks.begin(), gsfTracks.end());

  reco::TrackCollection ecalTracksTemp(generalTracks);
  for(const auto &track : gsfTracks)
    ecalTracksTemp.emplace_back(track);
  for(const auto &track : ecalTracksTemp)
    ecalTracks->emplace_back(track);

  // Remove overlap between superClusters and out-of-time superClusters
  reco::SuperClusterCollection superClusters(*superClusterHandle);
  reco::SuperClusterCollection ootSuperClusters(*ootSuperClusterHandle);

  superClusters = RemoveOverlap<reco::SuperCluster, reco::SuperCluster>(superClusters, ootSuperClusters, 0.1);
  displacedElectronSCs->insert(displacedElectronSCs->end(), superClusters.begin(), superClusters.end());

  iEvent.put(std::move(ecalGeneralTracks), "ecalGeneralTracks");
  iEvent.put(std::move(ecalGsfTracks), "ecalGsfTracks");
  iEvent.put(std::move(ecalTracks), "ecalTracks");
  //iEvent.put(std::move(generalTracksDetMatchInfo), "generalTracksDetMatchInfo");
  //iEvent.put(std::move(gsfTracksDetMatchInfo), "gsfTracksDetMatchInfo");
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

//define this as a plug-in
DEFINE_FWK_MODULE(ECALTracksProducer);
