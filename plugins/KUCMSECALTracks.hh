#ifndef KUCMSECALTracks_hh
#define KUCMSECALTracks_hh
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
#include "KUNTupleFW/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/Hungarian.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/MatchTracksToSC.h"
#include "KUNTupleFW/KUCMSNtupleizer/interface/TrackTools.h"

#include "TVector3.h"

//  KUCMS Object includes
#include "KUCMSGenObjects.hh"
#include "KUCMSBranchManager.hh"
#include "KUCMSObjectBase.hh"

//#define DEBUG true
#define DEBUG false

using namespace edm; 

class KUCMSECALTracks : public KUCMSObjectBase {

public:

  // use base class constructor
  KUCMSECALTracks( const edm::ParameterSet& iConfig );
  ~KUCMSECALTracks() = default;
  
  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadGeneralTrackTokens( edm::EDGetTokenT<edm::View<reco::Track>> token ){ generalTracksToken_ = token; }
  void LoadGsfTrackTokens( edm::EDGetTokenT<edm::View<reco::GsfTrack>> token ){ gsfTracksToken_ = token; }
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }  
  void LoadGenObject(KUCMSGenObject* genObjs){ genObjs_ = genObjs; };

  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree ); 
  
  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar );
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(); // do any jobs that need to be done after main event loop
    
private:

  typedef std::vector<std::vector<GlobalPoint> > TrackRecHitLocations;

  PropagatedTracks<reco::Track> generalECALTracks_;
  PropagatedTracks<reco::GsfTrack> gsfECALTracks_;
  std::vector<reco::TransientTrack> transientTracks_;

  TrackRecHitLocations generalTrackRecHitLocations_;
  TrackRecHitLocations gsfTrackRecHitLocations_;

  // General Tracks
  edm::EDGetTokenT<edm::View<reco::Track>> generalTracksToken_;
  edm::Handle<edm::View<reco::Track> > generalTracksHandle_;

  // Gsf Tracks
  edm::EDGetTokenT<edm::View<reco::GsfTrack>> gsfTracksToken_;
  edm::Handle<edm::View<reco::GsfTrack> > gsfTracksHandle_;

  KUCMSGenObject* genObjs_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  template <typename T>
  void FillTrackBranches(const PropagatedTracks<T> &propagatedTracks);

  template <typename T>
  void FillDetIdBranches(const std::vector<DetId> &detIDs, const PropagatedTrack<T> &trackDet, const bool isECAL);

};//<<>>class KUCMSECALTracks : public KUCMSObjectBase

KUCMSECALTracks::KUCMSECALTracks( const edm::ParameterSet& iConfig ) 
  : caloGeometryToken_(esConsumes()) {

  trackAssociator_.useDefaultPropagator();
  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
}//<<>>KUCMSECALTracks::KUCMSECALTracks( const edm::ParameterSet& iConfig )

void KUCMSECALTracks::InitObject( TTree* fOutTree ){

  Branches.makeBranch("ECALTrack_nGeneral", "ECALTrack_nGeneral", INT);
  Branches.makeBranch("ECALTrack_nGsf", "ECALTrack_nGsf", INT); 
  Branches.makeBranch("ECALTrack_nTracks", "ECALTrack_nTracks", INT);
  Branches.makeBranch("ECALTrack_charge", "ECALTrack_charge", VINT);
  Branches.makeBranch("ECALTrack_p", "ECALTrack_p", VFLOAT);
  Branches.makeBranch("ECALTrack_px", "ECALTrack_px", VFLOAT);
  Branches.makeBranch("ECALTrack_py", "ECALTrack_py", VFLOAT);
  Branches.makeBranch("ECALTrack_pz", "ECALTrack_pz", VFLOAT);
  Branches.makeBranch("ECALTrack_pt", "ECALTrack_pt", VFLOAT);
  Branches.makeBranch("ECALTrack_eta", "ECALTrack_eta", VFLOAT);
  Branches.makeBranch("ECALTrack_phi", "ECALTrack_phi", VFLOAT);
  Branches.makeBranch("ECALTrack_theta", "ECALTrack_theta", VFLOAT);
  Branches.makeBranch("ECALTrack_lambda", "ECALTrack_lambda", VFLOAT);
  Branches.makeBranch("ECALTrack_vx", "ECALTrack_vx", VFLOAT);
  Branches.makeBranch("ECALTrack_vy", "ECALTrack_vy", VFLOAT);
  Branches.makeBranch("ECALTrack_vz", "ECALTrack_vz", VFLOAT);
  Branches.makeBranch("ECALTrack_dxy", "ECALTrack_dxy", VFLOAT);
  Branches.makeBranch("ECALTrack_ndof", "ECALTrack_ndof", VFLOAT);
  Branches.makeBranch("ECALTrack_chi2", "ECALTrack_chi2", VFLOAT);
  Branches.makeBranch("ECALTrack_normalizedChi2", "ECALTrack_normalizedChi2", VFLOAT);
  Branches.makeBranch("ECALTrack_qoverp", "ECALTrack_qoverp", VFLOAT);
  Branches.makeBranch("ECALTrack_ptError", "ECALTrack_ptError", VFLOAT);
  Branches.makeBranch("ECALTrack_etaError", "ECALTrack_etaError", VFLOAT);
  Branches.makeBranch("ECALTrack_phiError", "ECALTrack_phiError", VFLOAT);
  Branches.makeBranch("ECALTrack_qualityMask", "ECALTrack_qualityMask", VINT);
  Branches.makeBranch("ECALTrack_nValidHits", "ECALTrack_nValidHits", VINT);
  Branches.makeBranch("ECALTrack_nLostHits", "ECALTrack_nLostHits", VINT);
  Branches.makeBranch("ECALTrack_isGeneral", "ECALTrack_isGeneral", VBOOL);
  Branches.makeBranch("ECALTrack_isGsf", "ECALTrack_isGsf", VBOOL);
  Branches.makeBranch("ECALTrackDetID_detId", "ECALTrackDetID_detId", VUINT);
  Branches.makeBranch("ECALTrackDetID_trackIndex", "ECALTrackDetID_trackIndex", VUINT);
  Branches.makeBranch("ECALTrackDetID_genIndex", "ECALTrackDetID_genIndex", VINT);
  Branches.makeBranch("ECALTrackDetID_isECAL", "ECALTrackDetID_isECAL", VBOOL);
  Branches.makeBranch("ECALTrackDetID_isHCAL", "ECALTrackDetID_isHCAL", VBOOL);

  Branches.attachBranches(fOutTree);
  
}//<<>>void KUCMSECALTracks::InitObject( TTree* fOutTree )

void KUCMSECALTracks::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){
  
  if( DEBUG ) std::cout << "Collecting Tracks" << std::endl;
  
  // Get event track and super cluster information from AOD
  iEvent.getByToken( generalTracksToken_, generalTracksHandle_ );
  iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_ );

  reco::TrackCollection generalTracks;
  reco::GsfTrackCollection gsfTracks;
  
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);
  
  // Collect general tracks
  for(const auto &track : *generalTracksHandle_) {
    generalTracks.emplace_back(track);
    transientTracks_.emplace_back(ttBuilder->build(track));
  }
  // Collect Gsf tracks
  for(const auto &track : *gsfTracksHandle_) {
    gsfTracks.emplace_back(track);
    transientTracks_.emplace_back(ttBuilder->build(track));
  }

  // Repropagate gsf and general tracks (ideally this would be done in the producer once but it is not straightforward)
  // Keep only tracks that reach the ECAL.
  TrackPropagator<reco::Track> generalTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, generalTracks);
  generalECALTracks_ = generalTrackPropagator.GetPropagatedTracks();

  TrackPropagator<reco::GsfTrack> gsfTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, gsfTracks);
  gsfECALTracks_ = gsfTrackPropagator.GetPropagatedTracks();


}//<<>>void KUCMSECALTracks::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSECALTracks::ProcessEvent( ItemManager<float>& geVar ){}

void KUCMSECALTracks::PostProcessEvent( ItemManager<float>& geVar ){
  
  if( DEBUG ) std::cout << "Processing Tracks" << std::endl;
  
  Branches.clearBranches();
  
  if( DEBUG ) std::cout << " - Entering Track loop" << std::endl;
  
  Branches.fillBranch("ECALTrack_nGeneral", int(generalECALTracks_.size()) );
  Branches.fillBranch("ECALTrack_nGsf", int(gsfECALTracks_.size()) );
  Branches.fillBranch("ECALTrack_nTracks", int(generalECALTracks_.size() + gsfECALTracks_.size()) );
  
  FillTrackBranches<reco::Track>(generalECALTracks_);
  FillTrackBranches<reco::GsfTrack>(gsfECALTracks_);
  
}//<<>>void KUCMSECALTracks::PostProcessEvent()

void KUCMSECALTracks::EndJobs(){}

template <typename T>
void KUCMSECALTracks::FillTrackBranches(const PropagatedTracks<T> &propagatedTracks) {

  if(cfFlag("hasGenInfo")) {
    
  }

  for(const auto &trackDet : propagatedTracks) {
    const T track = trackDet.GetTrack();
    const TrackDetMatchInfo detInfo = trackDet.GetDetInfo();
    
    Branches.fillBranch("ECALTrack_charge", int(track.charge()) );
    Branches.fillBranch("ECALTrack_p", float(track.p()) );
    Branches.fillBranch("ECALTrack_px", float(track.px()) );
    Branches.fillBranch("ECALTrack_py", float(track.py()) );
    Branches.fillBranch("ECALTrack_pz", float(track.pz()) );
    Branches.fillBranch("ECALTrack_pt", float(track.pt()) );
    Branches.fillBranch("ECALTrack_eta", float(track.eta()) );
    Branches.fillBranch("ECALTrack_phi", float(track.phi()) );
    Branches.fillBranch("ECALTrack_theta", float(track.theta()) );
    Branches.fillBranch("ECALTrack_lambda", float(track.lambda()) );
    Branches.fillBranch("ECALTrack_vx", float(track.vx()) );
    Branches.fillBranch("ECALTrack_vy", float(track.vy()) );
    Branches.fillBranch("ECALTrack_vz", float(track.vz()) );
    Branches.fillBranch("ECALTrack_dxy", float(track.dxy()) );
    Branches.fillBranch("ECALTrack_ndof", float(track.ndof()) );
    Branches.fillBranch("ECALTrack_chi2", float(track.chi2()) );
    Branches.fillBranch("ECALTrack_normalizedChi2", float(track.normalizedChi2()) );
    Branches.fillBranch("ECALTrack_qoverp", float(track.qoverp()) );
    Branches.fillBranch("ECALTrack_ptError", float(track.ptError()) );
    Branches.fillBranch("ECALTrack_etaError", float(track.etaError()) );
    Branches.fillBranch("ECALTrack_phiError", float(track.phiError()) );
    Branches.fillBranch("ECALTrack_qualityMask", int(track.qualityMask()) );
    Branches.fillBranch("ECALTrack_nValidHits", int(track.numberOfValidHits()) );
    Branches.fillBranch("ECALTrack_nLostHits", int(track.numberOfLostHits()) );
    Branches.fillBranch("ECALTrack_isGeneral", bool(typeid(T) == typeid(reco::Track)) );
    Branches.fillBranch("ECALTrack_isGsf", bool(typeid(T) == typeid(reco::GsfTrack)) );
    
    FillDetIdBranches<T>(detInfo.crossedEcalIds, trackDet, true);
    FillDetIdBranches<T>(detInfo.crossedHcalIds, trackDet, false);
    //Branches.fillBranch("",);
  }
}

template <typename T>
void KUCMSECALTracks::FillDetIdBranches(const std::vector<DetId> &detIDs, const PropagatedTrack<T> &trackDet, const bool isECAL) {
  for(const auto &detID : detIDs) {
    Branches.fillBranch("ECALTrackDetID_detId", unsigned(detID.rawId()) );
    Branches.fillBranch("ECALTrackDetID_trackIndex", unsigned(trackDet.GetIndex()) );
    Branches.fillBranch("ECALTrackDetID_isECAL", isECAL);
    Branches.fillBranch("ECALTrackDetID_isHCAL", !isECAL);
  }
}

#endif
