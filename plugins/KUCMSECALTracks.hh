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
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"

// Add includes for interface collections
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackInfo.h"

#include "TVector3.h"

//  KUCMS Object includes
#include "KUCMSGenObjects.hh"
#include "KUCMSBranchManager.hh"
#include "KUCMSObjectBase.hh"

#define ETDEBUG true
//#define ETDEBUG false

using namespace edm; 

class trackSCMatch {

public:
  trackSCMatch(const reco::Track &track, const reco::SuperCluster &sc, double deltaR) :
    track_(track), sc_(sc), deltaR_(deltaR) {}

  reco::Track GetTrack() const {return track_;}
  reco::SuperCluster GetSC() const {return sc_;}
  double GetDeltaR() const {return deltaR_;}
private:

  reco::Track track_;
  reco::SuperCluster sc_;
  double deltaR_;  
};

class KUCMSECALTracks : public KUCMSObjectBase {

public:

  // use base class constructor
  KUCMSECALTracks( const edm::ParameterSet& iConfig );
  ~KUCMSECALTracks() = default;
  
  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadECALTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ ecalTracksToken_ = token; }
  void LoadGeneralTrackTokens( edm::EDGetTokenT<edm::View<reco::Track>> token ){ generalTracksToken_ = token; }
  void LoadMergedSCs( edm::EDGetTokenT<reco::SuperClusterCollection> token ) {mergedSCsToken_ = token;}
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
  void LoadGenObject(KUCMSGenObject* genObjs){ genObjs_ = genObjs; };
  void LoadBeamSpot(edm::EDGetTokenT<reco::BeamSpot> token) { beamspotToken_ = token; }
 
  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree ); 
  
  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar ){}
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs() {} // do any jobs that need to be done after main event loop
    
private:

  typedef std::vector<std::vector<GlobalPoint> > TrackRecHitLocations;

  reco::TrackCollection generalTracks_;
  reco::TrackCollection ecalTracks_;
  PropagatedTracks<reco::Track> generalECALTracks_;
  PairedObjectCollection<reco::Track, reco::GenParticle> genMatches_;

  //ECAL tracks 
  edm::EDGetTokenT<reco::TrackCollection> ecalTracksToken_;
  edm::Handle<reco::TrackCollection> ecalTracksHandle_;
  
  // General Tracks
  edm::EDGetTokenT<edm::View<reco::Track>> generalTracksToken_;
  edm::Handle<edm::View<reco::Track> > generalTracksHandle_;

  // Merged SCs
  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle_;
  
  KUCMSGenObject* genObjs_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  reco::BeamSpot beamSpot_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  std::vector<trackSCMatch> trackSCMatches_;

  template <typename T>
  void FillTrackBranches(PropagatedTracks<T> &propagatedTracks);

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
  Branches.makeBranch("ECALTrack_genDeltaR", "ECALTrack_genDeltaR", VFLOAT);
  Branches.makeBranch("ECALTrackDetID_detId", "ECALTrackDetID_detId", VUINT);
  Branches.makeBranch("ECALTrackDetID_trackIndex", "ECALTrackDetID_trackIndex", VUINT);
  Branches.makeBranch("ECALTrackDetID_isECAL", "ECALTrackDetID_isECAL", VBOOL);
  Branches.makeBranch("ECALTrackDetID_isHCAL", "ECALTrackDetID_isHCAL", VBOOL);
  Branches.makeBranch("ECALTrack_genIndex", "ECALTrack_genIndex", VINT);
  Branches.makeBranch("ECALTrack_pdgId", "ECALTrack_pdgId", VINT);

  Branches.makeBranch("TrackSCMatch_scEnergyToTrackPRatio", "TrackSCMatch_scEnergyToTrackPRatio", VFLOAT);
  Branches.makeBranch("TrackSCMatch_deltaR", "TrackSCMatch_deltaR", VFLOAT);
  Branches.attachBranches(fOutTree);
  
}//<<>>void KUCMSECALTracks::InitObject( TTree* fOutTree )

void KUCMSECALTracks::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){
  
  if( ETDEBUG ) std::cout << "Collecting Tracks" << std::endl;

  ecalTracks_.clear();
  
  // Get event track and super cluster information from AOD
  iEvent.getByToken( ecalTracksToken_, ecalTracksHandle_ );
  iEvent.getByToken( generalTracksToken_, generalTracksHandle_ );
  iEvent.getByToken( mergedSCsToken_, mergedSCsHandle_);
  beamSpot_ = iEvent.get(beamspotToken_);

  reco::TrackCollection generalTracks;
  
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);
  
  // Collect general tracks
  for(const auto &track : *generalTracksHandle_) 
    generalTracks.emplace_back(track);

  // Repropagate gsf and general tracks (ideally this would be done in the producer once but it is not straightforward)
  TrackPropagator<reco::Track> generalTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, generalTracks);
  generalECALTracks_ = generalTrackPropagator.GetPropagatedTracks();
  for(const auto &propTrack : generalECALTracks_) ecalTracks_.emplace_back(propTrack.GetTrack());

  //std::cout << "\nThis event: " << std::endl;
  for(const auto &sc : *mergedSCsHandle_) {
    const double scEta(sc.eta());
    const double scPhi(sc.phi());
    double minDeltaR(999.);

    //std::cout << "  SuperCluster Info: eta = " << scEta << ", phi = " << scPhi << ", energy = " << sc.correctedEnergy() << std::endl;

    reco::Track matchedTrack;
    //int matchedID(-999);
    for(const auto &propTrack : generalECALTracks_) {

      const TrackDetMatchInfo detInfo(propTrack.GetDetInfo());
      
      for(const auto & detID : detInfo.crossedEcalIds) {
	const GlobalPoint trackHitAtEcal(ecalGeometry.getGeometry(detID)->getPosition());

	const double ecalEta(trackHitAtEcal.eta());
	const double ecalPhi(trackHitAtEcal.phi());
	const double deltaR(sqrt(reco::deltaR2(scEta, scPhi, ecalEta, ecalPhi)));

	if(deltaR < minDeltaR) {
	  minDeltaR = deltaR;
	  matchedTrack = propTrack.GetTrack();
	  //matchedID = int(detID);
	}
      }
    }
    //std::cout << "  Matched track: eta " << matchedTrack.eta() << ", phi = " << matchedTrack.phi() << ", p = " << matchedTrack.p() << std::endl;
    //std::cout << "    Match details: deltaR = " << minDeltaR << ", matched detID = " << matchedID << std::endl;
    trackSCMatches_.emplace_back(trackSCMatch(matchedTrack, sc, minDeltaR));
  }
  
  
}//<<>>void KUCMSECALTracks::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSECALTracks::PostProcessEvent( ItemManager<float>& geVar ){
  
  if( ETDEBUG ) std::cout << "Processing Tracks" << std::endl;
  
  Branches.clearBranches();
  if( ETDEBUG ) std::cout << " - Entering Track loop" << std::endl;

  Branches.fillBranch("ECALTrack_nTracks", int(ecalTracks_.size()) );
  
  FillTrackBranches<reco::Track>(generalECALTracks_);

  for(const auto &match : trackSCMatches_) {
    Branches.fillBranch("TrackSCMatch_scEnergyToTrackPRatio", float(match.GetSC().correctedEnergy()/match.GetTrack().p()));
    Branches.fillBranch("TrackSCMatch_deltaR", float(match.GetDeltaR()));
  }
}//<<>>void KUCMSECALTracks::PostProcessEvent()

template <typename T>
void KUCMSECALTracks::FillTrackBranches(PropagatedTracks<T> &propagatedTracks) {

  int dbgcnt(0);
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
    
    FillDetIdBranches<T>(detInfo.crossedEcalIds, trackDet, true);
    FillDetIdBranches<T>(detInfo.crossedHcalIds, trackDet, false);
    //Branches.fillBranch("",);

    dbgcnt++;
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

