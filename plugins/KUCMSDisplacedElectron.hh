#ifndef KUCMSDisplacedElectron_hh
#define KUCMSDisplacedElectron_hh
// -*- C++ -*-
//
//
// Original Author:  Andres Abreu
// 
// 
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

// Includes for all object types

#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"

#include "TrackingTools/Records/interface/DetIdAssociatorRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

// Add includes for specfic collections
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchTracksToSC.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"

#include "TVector3.h"

//  KUCMS Object includes
#include "KUCMSGenObjects.hh"
#include "KUCMSBranchManager.hh"
#include "KUCMSObjectBase.hh"
#include "DisplacedElectronIsolation.hh"

//#define DEBUG true
#define DEBUG false

template class MatchTracksToSC<reco::TransientTrack>;

using namespace edm; 

class KUCMSDisplacedElectron : public KUCMSObjectBase {

public:

  // use base class constructor
  KUCMSDisplacedElectron( const edm::ParameterSet& iConfig );
  virtual ~KUCMSDisplacedElectron() = default;
  
  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadECALTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ ecalTracksToken_ = token; }
  void LoadGeneralTrackTokens( edm::EDGetTokenT<edm::View<reco::Track>> token ){ generalTracksToken_ = token; }
  void LoadGsfTrackTokens( edm::EDGetTokenT<edm::View<reco::GsfTrack>> token ){ gsfTracksToken_ = token; }
  void LoadSuperClusterTokens( edm::EDGetTokenT<edm::View<reco::SuperCluster>> token ){ superClusterToken_ = token; }
  void LoadOotSuperClusterTokens( edm::EDGetTokenT<edm::View<reco::SuperCluster>> token ){ ootSuperClusterToken_ = token; }
  void LoadDisplacedElectrons( edm::EDGetTokenT<reco::ElectronCollection> token ) {electronToken_ = token;}
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
  void LoadGenObject(KUCMSGenObject* genObjs){ genObjs_ = genObjs; };
  void LoadIsolation(DisplacedElectronIsolation* isoTool) { isoTool_ = isoTool; }
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }

  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree ); 
  
  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar ) {}
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(){} // do any jobs that need to be done after main event loop

private:

  int nGeneralTracks_;
  int nGsfTracks_;

  std::vector<reco::TransientTrack> transientTracks_;
  reco::SuperClusterCollection mergedSuperClusters_;
  reco::GenParticleCollection signalGenElectrons_;

  //ECAL tracks
  edm::EDGetTokenT<reco::TrackCollection> ecalTracksToken_;
  edm::Handle<reco::TrackCollection> ecalTracksHandle_;

  // General Tracks
  edm::EDGetTokenT<edm::View<reco::Track>> generalTracksToken_;
  edm::Handle<edm::View<reco::Track> > generalTracksHandle_;

  // Gsf Tracks
  edm::EDGetTokenT<edm::View<reco::GsfTrack>> gsfTracksToken_;
  edm::Handle<edm::View<reco::GsfTrack> > gsfTracksHandle_;

  // In-Time SuperClusters
  edm::EDGetTokenT<edm::View<reco::SuperCluster>> superClusterToken_;
  edm::Handle<edm::View<reco::SuperCluster> > superClusterHandle_;

  // OOT SuperClusters
  edm::EDGetTokenT<edm::View<reco::SuperCluster>> ootSuperClusterToken_;
  edm::Handle<edm::View<reco::SuperCluster> > ootSuperClusterHandle_;

  // Displaced Electrons
  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;
  edm::Handle<reco::ElectronCollection> electronHandle_;

  KUCMSGenObject* genObjs_;
  DisplacedElectronIsolation* isoTool_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;
  
  MatchedTrackSCPairs<reco::TransientTrack> matchedTracksToSCs_;

  double SaveGenInfo(const reco::GenParticleCollection &genCollection);

};//<<>>class KUCMSDisplacedElectron : public KUCMSObjectBase

KUCMSDisplacedElectron::KUCMSDisplacedElectron( const edm::ParameterSet& iConfig ) 
  : caloGeometryToken_(esConsumes()) {

  trackAssociator_.useDefaultPropagator();
  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );

}//<<>>KUCMSDisplacedElectron::KUCMSDisplacedElectron( const edm::ParameterSet& iConfig )

void KUCMSDisplacedElectron::InitObject( TTree* fOutTree ){

  Branches.makeBranch("DisplacedElectron_nSCs","DisplacedElectron_nSCs", INT);
  Branches.makeBranch("DisplacedElectron_nOOTSCs","DisplacedElectron_nOOTSCs", INT);
  Branches.makeBranch("DisplacedElectron_nGenSig","DisplacedElectron_nGenSig", INT); 
  Branches.makeBranch("DisplacedElectron_nTotal", "DisplacedElectron_nTotal", INT);
  Branches.makeBranch("DisplacedElectron_charge", "DisplacedElectron_charge", VINT);
  Branches.makeBranch("DisplacedElectron_p", "DisplacedElectron_p", VFLOAT);
  Branches.makeBranch("DisplacedElectron_pt", "DisplacedElectron_pt", VFLOAT);
  Branches.makeBranch("DisplacedElectron_eta", "DisplacedElectron_eta", VFLOAT);
  Branches.makeBranch("DisplacedElectron_phi", "DisplacedElectron_phi", VFLOAT);
  Branches.makeBranch("DisplacedElectron_energy", "DisplacedElectron_energy", VFLOAT);
  Branches.makeBranch("DisplacedElectron_x", "DisplacedElectron_x", VFLOAT);
  Branches.makeBranch("DisplacedElectron_y", "DisplacedElectron_y", VFLOAT);
  Branches.makeBranch("DisplacedElectron_z", "DisplacedElectron_z", VFLOAT);
  Branches.makeBranch("DisplacedElectron_etaECAL", "DisplacedElectron_etaECAL", VFLOAT);
  Branches.makeBranch("DisplacedElectron_phiECAL", "DisplacedElectron_phiECAL", VFLOAT);
  Branches.makeBranch("DisplacedElectron_xECAL", "DisplacedElectron_xECAL", VFLOAT);
  Branches.makeBranch("DisplacedElectron_yECAL", "DisplacedElectron_yECAL", VFLOAT);
  Branches.makeBranch("DisplacedElectron_zECAL", "DisplacedElectron_zECAL", VFLOAT);
  Branches.makeBranch("DisplacedElectron_dxy", "DisplacedElectron_dxy", VFLOAT);
  Branches.makeBranch("DisplacedElectron_dRSC", "DisplacedElectron_dRSC", VFLOAT);
  Branches.makeBranch("DisplacedElectron_dRGen", "DisplacedElectron_dRGen", VFLOAT);
  Branches.makeBranch("DisplacedElectron_costGen", "DisplacedElectron_costGen", FLOAT);
  Branches.makeBranch("DisplacedElectron_costSC", "DisplacedElectron_costSC", FLOAT);
  Branches.makeBranch("DisplacedElectron_nHits", "DisplacedElectron_nHits", VINT);
  Branches.makeBranch("DisplacedElectron_indexGen", "DisplacedElectron_indexGen", VINT);
  Branches.makeBranch("DisplacedElectron_indexTrack", "DisplacedElectron_indexTrack", VINT);
  Branches.makeBranch("DisplacedElectron_indexSC", "DisplacedElectron_indexSC", VINT);
  Branches.makeBranch("DisplacedElectron_isGeneralTrack","DisplacedElectron_isGeneralTrack", VBOOL);
  Branches.makeBranch("DisplacedElectron_isGsfTrack","DisplacedElectron_isGsfTrack", VBOOL);
  Branches.makeBranch("DisplacedElectron_isSC","DisplacedElectron_isSC", VBOOL);
  Branches.makeBranch("DisplacedElectron_isOOTSC","DisplacedElectron_isOOTSC", VBOOL);
  Branches.makeBranch("DisplacedElectron_isBarrel","DisplacedElectron_isBarrel", VBOOL);
  Branches.makeBranch("DisplacedElectron_isEndCap","DisplacedElectron_isEndCap", VBOOL);
  Branches.makeBranch("DisplacedElectron_isSignal","DisplacedElectron_isSignal", VBOOL);
  Branches.makeBranch("DisplacedElectron_isTau","DisplacedElectron_isTau", VBOOL);
  Branches.makeBranch("DisplacedElectron_isLight","DisplacedElectron_isLight", VBOOL);
  Branches.makeBranch("DisplacedElectron_isHeavy","DisplacedElectron_isHeavy", VBOOL);
  Branches.makeBranch("DisplacedElectron_isConversion","DisplacedElectron_isConversion", VBOOL);
  Branches.makeBranch("DisplacedElectron_isUnmatched","DisplacedElectron_isUnmatched", VBOOL);
  Branches.makeBranch("DisplacedElectron_genSigIndex", "DisplacedElectron_genSigIndex", VINT);
  Branches.makeBranch("DisplacedElectron_pdgId", "DisplacedElectron_pdgId", VINT);
  Branches.makeBranch("DisplacedElectron_trackIso0p4", "DisplacedElectron_trackIso0p4", VFLOAT);
  Branches.makeBranch("DisplacedElectron_pfIso0p4", "DisplacedElectron_pfIso0p4", VFLOAT);
  Branches.makeBranch("DisplacedElectron_scIso0p4", "DisplacedElectron_scIso0p4", VFLOAT);
  Branches.makeBranch("DisplacedElectron_ip3D", "DisplacedElectron_ip3D", VFLOAT);
  Branches.makeBranch("DisplacedElectron_ip2D", "DisplacedElectron_ip2D", VFLOAT);
  Branches.makeBranch("DisplacedElectron_sip3D", "DisplacedElectron_sip3D", VFLOAT);
  Branches.makeBranch("DisplacedElectron_sip2D", "DisplacedElectron_sip2D", VFLOAT);
  Branches.makeBranch("DisplacedElectron_ip3Dbs", "DisplacedElectron_ip3Dbs", VFLOAT);
  Branches.makeBranch("DisplacedElectron_ip2Dbs", "DisplacedElectron_ip2Dbs", VFLOAT);
  Branches.makeBranch("DisplacedElectron_sip3Dbs", "DisplacedElectron_sip3Dbs", VFLOAT);
  Branches.makeBranch("DisplacedElectron_sip2Dbs", "DisplacedElectron_sip2Dbs", VFLOAT);
  Branches.attachBranches(fOutTree);
  
}//<<>>void KUCMSDisplacedElectron::InitObject( TTree* fOutTree )

void KUCMSDisplacedElectron::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  if( DEBUG ) std::cout << "Collecting Electron Building Materials" << std::endl;

  // Get event track and super cluster information from AOD
  iEvent.getByToken( ecalTracksToken_, ecalTracksHandle_ );
  iEvent.getByToken( generalTracksToken_, generalTracksHandle_ );
  iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_ );
  iEvent.getByToken( superClusterToken_, superClusterHandle_ );
  iEvent.getByToken( ootSuperClusterToken_, ootSuperClusterHandle_ );
  iEvent.getByToken( electronToken_, electronHandle_ );

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  // Clear track and SCs collections  
  transientTracks_.clear();
  mergedSuperClusters_.clear();
  matchedTracksToSCs_.clear();

  // Collect general tracks
  nGeneralTracks_ = 0; 
  for(const auto &track : *generalTracksHandle_) {
    transientTracks_.emplace_back(ttBuilder->build(track));
    nGeneralTracks_++;
  }

  // Collect Gsf tracks
  nGsfTracks_ = 0;
  for(const auto &track : *gsfTracksHandle_) {
    transientTracks_.emplace_back(ttBuilder->build(track));
    nGsfTracks_++;
  }

  // Collect SuperClusters
  for(const auto &sc : *superClusterHandle_)
    mergedSuperClusters_.emplace_back(sc);

  // Merge SuperCluster collections
  mergedSuperClusters_.insert(mergedSuperClusters_.end(), ootSuperClusterHandle_->begin(), ootSuperClusterHandle_->end());

  // Match entire track collection to merged supercluster collection and find the electron candidates in the event
  if(transientTracks_.size() > 0 && mergedSuperClusters_.size() > 0) {
    MatchTracksToSC<reco::TransientTrack> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, transientTracks_, mergedSuperClusters_);
    matchedTracksToSCs_ = assigner.GetMatchedTrackSCPairs();
  }

}//<<>>void KUCMSDisplacedElectron::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSDisplacedElectron::PostProcessEvent( ItemManager<float>& geVar ){

  if( DEBUG ) std::cout << "Processing Electron Materials" << std::endl;

  Branches.clearBranches();
  signalGenElectrons_.clear();

  // Get Gen-Matching information, if available
  float genCost = -1;
  if(cfFlag("hasGenInfo")) {

    signalGenElectrons_ = genObjs_->GetSignalGenElectrons();

    int index(0);
    for(auto const &gen : genObjs_->GetGenParticles()) {
      if(abs(gen.pdgId()) == 11 && gen.status() == 1) {
        const LepMomType type(genObjs_->ClassifyGenElectron(gen));
        if(type == kZ || type == kSusy)
	  Branches.fillBranch("DisplacedElectron_genSigIndex", int(index));
      }
      index++;
    }

    genCost = SaveGenInfo(genObjs_->GetGenParticles());
  }

  if( DEBUG ) std::cout << " - Entering Electron Builder loop" << std::endl;

  Branches.fillBranch("DisplacedElectron_nGenSig", int(signalGenElectrons_.size()) );
  Branches.fillBranch("DisplacedElectron_nSCs",    int(superClusterHandle_->size()) );
  Branches.fillBranch("DisplacedElectron_nOOTSCs", int(ootSuperClusterHandle_->size()) );
  Branches.fillBranch("DisplacedElectron_nTotal",  int(matchedTracksToSCs_.size()) );
  Branches.fillBranch("DisplacedElectron_costSC",  float(matchedTracksToSCs_.GetCost()) );
  Branches.fillBranch("DisplacedElectron_costGen", float(genCost));
  
  isoTool_->CalculateIsolation(matchedTracksToSCs_);

  int index = 0;
  for(auto const &candidate : matchedTracksToSCs_) {

    const auto track(candidate.GetTrack().track());
    const reco::SuperCluster superCluster(candidate.GetSuperCluster());
    const GlobalPoint trackAtECAL(candidate.GetTrackLocationAtECAL());
    const int trackIndex(candidate.GetTrackIndex());
    const int scIndex(candidate.GetSCIndex()); 
    
    Branches.fillBranch("DisplacedElectron_charge",         int(track.charge()) );
    Branches.fillBranch("DisplacedElectron_p",              float(track.p()) );
    Branches.fillBranch("DisplacedElectron_pt",             float(track.pt()) );
    Branches.fillBranch("DisplacedElectron_eta",            float(track.eta()) );
    Branches.fillBranch("DisplacedElectron_phi",            float(track.phi()) );
    Branches.fillBranch("DisplacedElectron_energy",         float(superCluster.rawEnergy()) );
    Branches.fillBranch("DisplacedElectron_x",              float(track.vx()) );
    Branches.fillBranch("DisplacedElectron_y",              float(track.vy()) );
    Branches.fillBranch("DisplacedElectron_z",              float(track.vz()) );
    Branches.fillBranch("DisplacedElectron_etaECAL",        float(trackAtECAL.eta()) );
    Branches.fillBranch("DisplacedElectron_phiECAL",        float(trackAtECAL.phi()) );
    Branches.fillBranch("DisplacedElectron_xECAL",          float(trackAtECAL.x()) );
    Branches.fillBranch("DisplacedElectron_yECAL",          float(trackAtECAL.y()) );
    Branches.fillBranch("DisplacedElectron_zECAL",          float(trackAtECAL.z()) );
    Branches.fillBranch("DisplacedElectron_dxy",            float(TrackHelper::GetDXY(track)) );
    Branches.fillBranch("DisplacedElectron_dRSC",           float(candidate.GetDeltaR()) );
    Branches.fillBranch("DisplacedElectron_dRGen",          float(candidate.GetGenInfo().GetDeltaR()) );
    Branches.fillBranch("DisplacedElectron_nHits",          int(track.numberOfValidHits()) );
    Branches.fillBranch("DisplacedElectron_indexGen",       int(candidate.GetGenInfo().GetIndex()) );
    Branches.fillBranch("DisplacedElectron_isBarrel",       bool(fabs(trackAtECAL.eta()) <= 1.479) );
    Branches.fillBranch("DisplacedElectron_isEndCap",       bool(fabs(trackAtECAL.eta()) > 1.479) );
    Branches.fillBranch("DisplacedElectron_isGeneralTrack", bool(candidate.GetTrackIndex() < int(nGeneralTracks_)) );
    Branches.fillBranch("DisplacedElectron_isGsfTrack",     bool(candidate.GetTrackIndex() >= int(nGeneralTracks_)) );
    Branches.fillBranch("DisplacedElectron_isSC",           bool(candidate.GetSCIndex() < int(superClusterHandle_->size())) );
    Branches.fillBranch("DisplacedElectron_isOOTSC",        bool(candidate.GetSCIndex() >= int(superClusterHandle_->size())) );
    Branches.fillBranch("DisplacedElectron_indexTrack",     trackIndex);
    Branches.fillBranch("DisplacedElectron_indexSC",        scIndex);
    Branches.fillBranch("DisplacedElectron_trackIso0p4",    float(candidate.GetIsoInfo()["trackIso0p4"]) );
    Branches.fillBranch("DisplacedElectron_pfIso0p4",       float(candidate.GetIsoInfo()["pfIso0p4"]) );
    Branches.fillBranch("DisplacedElectron_scIso0p4",       float(candidate.GetIsoInfo()["scIso0p4"]) );
    Branches.fillBranch("DisplacedElectron_ip3D",           float(candidate.GetIsoInfo()["ip3D"]) );
    Branches.fillBranch("DisplacedElectron_ip2D",           float(candidate.GetIsoInfo()["ip2D"]) );
    Branches.fillBranch("DisplacedElectron_sip3D",          float(candidate.GetIsoInfo()["sip3D"]) );
    Branches.fillBranch("DisplacedElectron_sip2D",          float(candidate.GetIsoInfo()["sip2D"]) );
    Branches.fillBranch("DisplacedElectron_ip3Dbs",         float(candidate.GetIsoInfo()["ip3Dbs"]) );
    Branches.fillBranch("DisplacedElectron_ip2Dbs",         float(candidate.GetIsoInfo()["ip2Dbs"]) );
    Branches.fillBranch("DisplacedElectron_sip3Dbs",        float(candidate.GetIsoInfo()["sip3Dbs"]) );
    Branches.fillBranch("DisplacedElectron_sip2Dbs",        float(candidate.GetIsoInfo()["sip2Dbs"]) );
    
    if(cfFlag("hasGenInfo")) {
      GenLeptonInfo genInfo(candidate.GetGenInfo());
      Branches.fillBranch("DisplacedElectron_isSignal",     bool(genInfo.GetLepMomType() == kZ || genInfo.GetLepMomType() == kSusy) );
      Branches.fillBranch("DisplacedElectron_isTau",        bool(genInfo.GetLepMomType() == kTau) );
      Branches.fillBranch("DisplacedElectron_isLight",      bool(genInfo.GetLepMomType() == kLight) );
      Branches.fillBranch("DisplacedElectron_isHeavy",      bool(genInfo.GetLepMomType() == kHeavy) );
      Branches.fillBranch("DisplacedElectron_isConversion", bool(genInfo.GetLepMomType() == kConversion) );
      Branches.fillBranch("DisplacedElectron_isUnmatched",  bool(genInfo.GetLepMomType() == kUnmatched) );
      Branches.fillBranch("DisplacedElectron_pdgId",        int(genInfo.GetPdgID()) );

    }
    index++;
  }
}//<<>>void KUCMSDisplacedElectron::PostProcessEvent()

double KUCMSDisplacedElectron::SaveGenInfo(const reco::GenParticleCollection &genCollection) {

  if(matchedTracksToSCs_.size() == 0)
    return -1;

  DeltaRGenMatchHungarian<reco::Electron> genToTrackAssigner(*electronHandle_, genCollection);
  for(auto const &pair : genToTrackAssigner.GetPairedObjects()) {
    const reco::GenParticle genElectron(pair.GetObjectB());
    
    LepMomType type = kUnmatched;
    if(abs(genElectron.pdgId()) == 11)
      type = genObjs_->ClassifyGenElectron(genElectron);
    
    GenLeptonInfo genType(pair.GetIndexB(), genElectron, type, pair.GetDeltaR());
    matchedTracksToSCs_[pair.GetIndexA()].setGenMatchInfo(genType);
  }
  
  return genToTrackAssigner.GetCost();

}

#endif
