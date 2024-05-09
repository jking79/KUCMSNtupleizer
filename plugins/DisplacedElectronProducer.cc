#ifndef DisplacedElectronProducer_cc
#define DisplacedElectronProducer_cc
// -*- C++ -*-
//
// Package:    KUCMSNtupleizer/KUCMSNtupleizer
// Class:      DisplacedElectronProducer
//
/**\class DisplacedElectronProducer DisplacedElectronProducer.cc KUCMSNtupleizer/KUCMSNtupleizer/plugins/DisplacedElectronProducer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andres Abreu
//         Created:  Wed, 10 Jan 2024 12:11:43 GMT
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

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "DataFormats/Math/interface/LorentzVector.h"

//Local includes
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/Hungarian.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchTracksToSC.h"

//
// class declaration
//

class DisplacedElectronProducer : public edm::stream::EDProducer<> {

public:
  explicit DisplacedElectronProducer(const edm::ParameterSet&);
  virtual ~DisplacedElectronProducer() = default;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  
  void produce(edm::Event&, const edm::EventSetup&) override;

  reco::ElectronCollection MakeElectrons(const MatchedTrackSCPairs<reco::Track> &candidates) const;

  // ----------member data ---------------------------  
  edm::EDGetTokenT<reco::TrackCollection> ecalTrackToken_;
  edm::EDGetTokenT<reco::TrackCollection> generalTrackToken_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfElectronTrackToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> superClusterToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> ootSuperClusterToken_;
  //edm::EDGetTokenT<reco::GenParticleCollection> genParticleToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  edm::Handle<reco::TrackCollection> ecalTracksHandle_;
  edm::Handle<reco::TrackCollection> generalTracksHandle_;
  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle_;
  edm::Handle<reco::SuperClusterCollection> superClusterHandle_;
  edm::Handle<reco::SuperClusterCollection> ootSuperClusterHandle_;
  //edm::Handle<reco::GenParticleCollection> genParticleHandle_;
};

DisplacedElectronProducer::DisplacedElectronProducer(const edm::ParameterSet& iConfig) :
  ecalTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("ecalTracksSrc")) ),
  generalTrackToken_(consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("generalTracksSrc")) ),
  gsfElectronTrackToken_(consumes<reco::GsfTrackCollection>(iConfig.getParameter<edm::InputTag>("gsfElectronTracksSrc")) ),
  superClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters")) ),
  ootSuperClusterToken_(consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters")) ),
//genParticleToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles")) ),
  caloGeometryToken_(esConsumes()),
  magneticFieldToken_(esConsumes<MagneticField, IdealMagneticFieldRecord>()) 
{

  // TrackAssociator parameters
  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters_.loadParameters(parameters, iC);
  trackAssociator_.useDefaultPropagator();

  produces<reco::TrackCollection>("displacedCandidateTracks").setBranchAlias("displacedCandidateTracks");
  produces<reco::SuperClusterCollection>("displacedCandidateSCs").setBranchAlias("displacedCandidateSCs");
  produces<reco::ElectronCollection>("displacedElectrons").setBranchAlias("displacedElectrons");
  //produces<reco::ElectronCollection>("signalDisplacedElectrons").setBranchAlias("signalDisplacedElectrons");
}

//
// member functions
//
// ------------ method called to produce the data  ------------
void DisplacedElectronProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace std;

  auto magfield = iSetup.getTransientHandle(magneticFieldToken_);
  CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  iEvent.getByToken(ecalTrackToken_, ecalTracksHandle_);
  iEvent.getByToken(generalTrackToken_, generalTracksHandle_);
  iEvent.getByToken(gsfElectronTrackToken_, gsfTracksHandle_);
  iEvent.getByToken(superClusterToken_, superClusterHandle_);
  iEvent.getByToken(ootSuperClusterToken_, ootSuperClusterHandle_);
  //iEvent.getByToken(genParticleToken_, genParticleHandle_);

  // Smart pointers to containers of output collections
  std::unique_ptr<reco::TrackCollection> displacedTracks = std::make_unique<reco::TrackCollection>();
  std::unique_ptr<reco::SuperClusterCollection> displacedElectronSCs = std::make_unique<reco::SuperClusterCollection>();
  std::unique_ptr<reco::ElectronCollection> displacedElectrons = std::make_unique<reco::ElectronCollection>();
  //std::unique_ptr<reco::ElectronCollection> signalDisplacedElectrons = std::make_unique<reco::ElectronCollection>();  

  // Remove overlap between superClusters and out-of-time superClusters
  reco::SuperClusterCollection superClusters(*superClusterHandle_);
  superClusters.insert(superClusters.end(), ootSuperClusterHandle_->begin(), ootSuperClusterHandle_->end());
 
 // Match ECAL tracks to super clusters
  MatchedTrackSCPairs<reco::Track> matchedTracksToSCs;
  if(ecalTracksHandle_->size() > 0) {
    MatchTracksToSC<reco::Track> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, *ecalTracksHandle_, superClusters);
    matchedTracksToSCs = assigner.GetMatchedTrackSCPairs();
  }

  // Fill the output containers with their appropriate objects
  for(const auto &pair : matchedTracksToSCs) {
    displacedTracks->emplace_back(pair.GetTrack());
    displacedElectronSCs->emplace_back(pair.GetSuperCluster());
  }

  // Construct reco::Electron collection from matched pairs
  reco::ElectronCollection electrons(MakeElectrons(matchedTracksToSCs));
  displacedElectrons->insert(displacedElectrons->end(), electrons.begin(), electrons.end());

  // Find signal electrons
  //reco::ElectronCollection signalElectrons = GetSignalElectrons(electrons, *genParticleHandle_);
  //signalDisplacedElectrons->insert(signalDisplacedElectrons->end(), signalElectrons.begin(), signalElectrons.end());

  iEvent.put(std::move(displacedTracks), "displacedCandidateTracks"); 
  iEvent.put(std::move(displacedElectronSCs), "displacedCandidateSCs");
  iEvent.put(std::move(displacedElectrons), "displacedElectrons");  
  //iEvent.put(std::move(signalDisplacedElectrons), "signalDisplacedElectrons");

}// Producer end

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void DisplacedElectronProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  /*
  desc.add<edm::InputTag>("generalTracksSrc", edm::InputTag("ecalTracks", "ecalGeneralTracks"));
  desc.add<edm::InputTag>("gsfElectronTracksSrc", edm::InputTag("ecalTracks", "ecalGsfTracks"));
  desc.add<edm::InputTag>("superClusters", edm::InputTag("ecalTracks", "displacedElectronSCs"));
  desc.add<edm::InputTag>("ootSuperClusters", edm::InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"));
  desc.add<edm::ESInputTag>("tkAssocParamBlock", edm::ESInputTag("", "TransientTrackBuilder"));
  */
  descriptions.addDefault(desc);
}

reco::ElectronCollection DisplacedElectronProducer::MakeElectrons(const MatchedTrackSCPairs<reco::Track> &candidates) const {

  const int nSCs = superClusterHandle_->size();

  reco::ElectronCollection electrons;

  for(const auto &pair : candidates) {
    const reco::Track track(pair.GetTrack());

    const ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > electron4Vec(track.px(), track.py(), track.pz(), pair.GetSuperCluster().rawEnergy());
    reco::Electron matchedElectron(track.charge(), electron4Vec, math::XYZPoint(track.vx(), track.vy(), track.vz()));

    reco::SuperClusterRef scRef = (pair.GetSCIndex() < nSCs) ? reco::SuperClusterRef(superClusterHandle_, pair.GetSCIndex()) :
                                                               reco::SuperClusterRef(ootSuperClusterHandle_, pair.GetSCIndex() - nSCs);    
    matchedElectron.setSuperCluster(scRef);
    matchedElectron.setTrack(reco::TrackRef(ecalTracksHandle_, pair.GetTrackIndex()));

    electrons.push_back(matchedElectron);
  }

  return electrons;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedElectronProducer);

#endif
