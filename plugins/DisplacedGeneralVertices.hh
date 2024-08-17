#ifndef DisplacedGeneralVertices_hh
#define DisplacedGeneralVertices_hh
// -*- C++ -*-
//
//
// Original Author:  Andres Abreu
// 
// DisplacedGeneralVertices class
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

// Includes for all object types

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

// Add includes for interface collections
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedTrackSCPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexAssembly.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"

// KUCMS Object includes
#include "KUCMSObjectBase.hh"
#include "KUCMSBranchManager.hh"

#define DEBUG false

class DisplacedGeneralVertices : public KUCMSObjectBase {

public:

  // use base class constructor
  DisplacedGeneralVertices( const edm::ParameterSet& iConfig );
  virtual ~DisplacedGeneralVertices() = default; 

  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadGeneralTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ generalTracksToken_ = token; } 
  void LoadGsfTracksToken( edm::EDGetTokenT<reco::GsfTrackCollection> token ){ gsfTracksToken_ = token; }
  void LoadECALTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ ecalTracksToken_ = token; }
  void LoadDisplacedElectrons( edm::EDGetTokenT<reco::ElectronCollection> token ) {electronToken_ = token;}
  void LoadMergedSCs( edm::EDGetTokenT<reco::SuperClusterCollection> token ) {mergedSCsToken_ = token;}
  void LoadPrimaryVertex( edm::EDGetTokenT<reco::VertexCollection> token ) {pvToken_ = token;}
  void LoadPFCandidatesToken( edm::EDGetTokenT<reco::PFCandidateCollection> token ){ pfCandidatesToken_ = token; }
  void LoadBeamSpot(edm::EDGetTokenT<reco::BeamSpot> token) { beamspotToken_ = token; }
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }  
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> token ){ genToken_ = token; }

  //void LoadDisplacedElectronObj(const KUCMSDisplacedElectron* electrons) { displacedElectronObj_ = electrons; }
  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree );

  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar ) {}
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs() {} // do any jobs that need to be done after main event loop
  
private:

  reco::TrackCollection ecalTracks_;
  reco::TrackCollection tracks_;
  reco::TrackCollection signalTracks_;
  std::vector<reco::TransientTrack> ttracks_;
  reco::VertexCollection generalVertices_;
  PropagatedTracks<reco::Track> propECALTracks_;
  PairedObjectCollection<reco::Track, reco::GenParticle> genMatchedTracks_;

  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
  edm::Handle<reco::TrackCollection> generalTracksHandle_;

  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle_;

  edm::EDGetTokenT<reco::TrackCollection> ecalTracksToken_;
  edm::Handle<reco::TrackCollection> ecalTracksHandle_;

  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;
  edm::Handle<reco::ElectronCollection> electronHandle_;

  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle_;

  MatchedTrackSCPairs<reco::Track> matchedTracksToSCs_;
  
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::Handle<reco::GenParticleCollection> genHandle_;

  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  reco::BeamSpot beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  reco::Vertex primaryVertex_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  VertexAssembly vertexBuilder_;

  bool VertexHasSignal(const reco::Vertex &vertex, int &counter) const;

};//<<>>class DisplacedGeneralVertices : public KUCMSObjectBase

DisplacedGeneralVertices::DisplacedGeneralVertices( const edm::ParameterSet& iConfig ) 
  : caloGeometryToken_(esConsumes()) {   
  // ---- end constructor initilizations  --------------------------
  trackAssociator_.useDefaultPropagator();
  
  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
  
}//<<>>DisplacedGeneralVertices::DisplacedGeneralVertices( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void DisplacedGeneralVertices::InitObject( TTree* fOutTree ) {

  //Beamspot
  Branches.makeBranch("Beamspot_x", "Beamspot_x", FLOAT);
  Branches.makeBranch("Beamspot_y", "Beamspot_y", FLOAT);
  Branches.makeBranch("Beamspot_z", "Beamspot_z", FLOAT);
  
  // Track variables
  Branches.makeBranch("GeneralTrack_nTotal","GeneralTrack_nTotal", UINT);
  Branches.makeBranch("Track_nTotal","Track_nTotal", UINT);
  Branches.makeBranch("Track_nECALTracks","Track_nECALTracks", UINT);
  Branches.makeBranch("Track_ecalTrackIndex","Track_ecalTrackIndex", VINT);
  Branches.makeBranch("Track_pt","Track_pt", VFLOAT);
  Branches.makeBranch("Track_eta","Track_eta", VFLOAT);
  Branches.makeBranch("Track_phi","Track_phi", VFLOAT);
  Branches.makeBranch("Track_charge", "Track_charge", VINT);
  Branches.makeBranch("Track_p", "Track_p", VFLOAT);
  Branches.makeBranch("Track_px", "Track_px", VFLOAT);
  Branches.makeBranch("Track_py", "Track_py", VFLOAT);
  Branches.makeBranch("Track_pz", "Track_pz", VFLOAT);
  Branches.makeBranch("Track_vx", "Track_vx", VFLOAT);
  Branches.makeBranch("Track_vy", "Track_vy", VFLOAT);
  Branches.makeBranch("Track_vz", "Track_vz", VFLOAT);
  Branches.makeBranch("Track_dxy", "Track_dxy", VFLOAT);
  Branches.makeBranch("Track_ndof", "Track_ndof", VFLOAT);
  Branches.makeBranch("Track_chi2", "Track_chi2", VFLOAT);
  Branches.makeBranch("Track_normalizedChi2", "Track_normalizedChi2", VFLOAT);
  Branches.makeBranch("Track_qoverp", "Track_qoverp", VFLOAT);
  Branches.makeBranch("Track_ptError", "Track_ptError", VFLOAT);
  Branches.makeBranch("Track_etaError", "Track_etaError", VFLOAT);
  Branches.makeBranch("Track_phiError", "Track_phiError", VFLOAT);
  Branches.makeBranch("Track_ip2D", "Track_ip2D", VFLOAT);
  Branches.makeBranch("Track_sip2D", "Track_sip2D", VFLOAT);
  Branches.makeBranch("Track_ip3D", "Track_ip3D", VFLOAT);
  Branches.makeBranch("Track_sip3D", "Track_sip3D", VFLOAT);
  Branches.makeBranch("Track_qualityMask", "Track_qualityMask", VINT);
  Branches.makeBranch("Track_nValidHits", "Track_nValidHits", VINT);
  Branches.makeBranch("Track_nLostHits", "Track_nLostHits", VINT);
  Branches.makeBranch("Track_isECAL", "Track_isECAL", VBOOL);
  Branches.makeBranch("Track_isSCmatch", "Track_isSCmatch", VBOOL);
  Branches.makeBranch("Track_SCDR", "Track_SCDR", VFLOAT);
  Branches.makeBranch("Track_energySC", "Track_energySC", VFLOAT);
  Branches.makeBranch("Track_ratioPtToEnergySC", "Track_ratioPtToEnergySC", VFLOAT);
  
  // Gen Particle Variables
  Branches.makeBranch("GenParticle_nMatches","GenParticle_nMatches", UINT);
  Branches.makeBranch("GenParticle_nSignal","GenParticle_nSignal", UINT);
  Branches.makeBranch("GenParticle_pdgId","GenParticle_pdgId", VINT);
  Branches.makeBranch("GenParticle_matchedTrackIndex","GenParticle_matchedTrackIndex", VINT);
  Branches.makeBranch("GenParticle_charge","GenParticle_charge", VINT);
  Branches.makeBranch("GenParticle_deltaR","GenParticle_deltaR", VFLOAT);
  Branches.makeBranch("GenParticle_pt","GenParticle_pt", VFLOAT);
  Branches.makeBranch("GenParticle_eta","GenParticle_eta", VFLOAT);
  Branches.makeBranch("GenParticle_phi","GenParticle_phi", VFLOAT);
  Branches.makeBranch("GenParticle_x","GenParticle_x", VFLOAT);
  Branches.makeBranch("GenParticle_y","GenParticle_y", VFLOAT);
  Branches.makeBranch("GenParticle_z","GenParticle_z", VFLOAT);
  Branches.makeBranch("GenParticle_dxy","GenParticle_dxy", VFLOAT);
  Branches.makeBranch("GenParticle_p","GenParticle_p", VFLOAT);
  Branches.makeBranch("GenParticle_isSignal","GenParticle_isSignal", VBOOL);
  
  // Primary Vertex
  Branches.makeBranch("PV_nTracks","PV_nTracks", UINT);
  Branches.makeBranch("PV_signalCount","PV_signalCount", UINT);
  Branches.makeBranch("PV_hasSignal","PV_hasSignal", BOOL);
  Branches.makeBranch("PV_x","PV_x", FLOAT);
  Branches.makeBranch("PV_y","PV_y", FLOAT);
  Branches.makeBranch("PV_z","PV_z", FLOAT);
  Branches.makeBranch("PV_sumPt","PV_sumPt", FLOAT);
  Branches.makeBranch("PV_weightedSumPt","PV_weightedSumPt", FLOAT);
  Branches.makeBranch("PV_chi2","PV_chi2", FLOAT);
  Branches.makeBranch("PV_ndof","PV_ndof", FLOAT);
  Branches.makeBranch("PV_normalizedChi2","PV_normalizedChi2", FLOAT);
  Branches.makeBranch("PV_ecalness","PV_ecalness", FLOAT);
  Branches.makeBranch("PV_trackIndex","PV_trackIndex", VUINT);
  Branches.makeBranch("PV_trackWeight","PV_trackWeight", VFLOAT);
  Branches.makeBranch("PVCollection_nTotal", "PVCollection_nTotal", UINT);
  Branches.makeBranch("PVCollection_nTracks", "PVCollection_nTracks", UINT, "Total unique tracks in PV");
  Branches.makeBranch("PVCollection_signalWeight", "PVCollection_signalWeight", VFLOAT, "weight of matched signal track in vertex");
  Branches.makeBranch("PVCollection_pvIndex", "PVCollection_pvIndex", VUINT);
  Branches.makeBranch("PVCollection_trackIndex", "PVCollection_trackIndex", VUINT);
  Branches.makeBranch("PVCollection_genIndex", "PVCollection_genIndex", VUINT);
  Branches.makeBranch("PVCollection_x", "PVCollection_x", VFLOAT);
  Branches.makeBranch("PVCollection_y", "PVCollection_y", VFLOAT);
  Branches.makeBranch("PVCollection_z", "PVCollection_z", VFLOAT);

  // Vertex Variables
  Branches.makeBranch("Vertex_nTotal","Vertex_nTotal", UINT);
  Branches.makeBranch("Vertex_signalCount","Vertex_signalCount", VUINT);
  Branches.makeBranch("Vertex_hasSignal","Vertex_hasSignal", VBOOL);
  Branches.makeBranch("Vertex_nTracks","Vertex_nTracks", VUINT);
  Branches.makeBranch("Vertex_mass","Vertex_mass", VFLOAT, "assigned to zero");
  Branches.makeBranch("Vertex_x","Vertex_x", VFLOAT);
  Branches.makeBranch("Vertex_y","Vertex_y", VFLOAT);
  Branches.makeBranch("Vertex_z","Vertex_z", VFLOAT);
  Branches.makeBranch("Vertex_p","Vertex_p", VFLOAT);
  Branches.makeBranch("Vertex_px", "Vertex_px", VFLOAT);
  Branches.makeBranch("Vertex_py", "Vertex_py", VFLOAT);
  Branches.makeBranch("Vertex_pz", "Vertex_pz", VFLOAT);
  Branches.makeBranch("Vertex_pt","Vertex_pt", VFLOAT);
  Branches.makeBranch("Vertex_eta","Vertex_eta", VFLOAT);
  Branches.makeBranch("Vertex_phi","Vertex_phi", VFLOAT);
  Branches.makeBranch("Vertex_dxy","Vertex_dxy", VFLOAT);
  Branches.makeBranch("Vertex_chi2","Vertex_chi2", VFLOAT);
  Branches.makeBranch("Vertex_normalizedChi2","Vertex_normalizedChi2", VFLOAT);
  Branches.makeBranch("Vertex_ndof","Vertex_ndof", VFLOAT);
  Branches.makeBranch("Vertex_ecalness","Vertex_ecalness", VFLOAT);
  Branches.makeBranch("Vertex_sumPt","Vertex_sumPt", VFLOAT);
  Branches.makeBranch("Vertex_vertexIndex","Vertex_vertexIndex", VUINT);
  Branches.makeBranch("Vertex_trackIndex","Vertex_trackIndex", VUINT);
  Branches.makeBranch("Vertex_trackWeight","Vertex_trackWeight", VFLOAT);
  Branches.makeBranch("Vertex_totalTrackWeight","Vertex_totalTrackWeight", VFLOAT);
  Branches.makeBranch("Vertex_trackRatioWithPV","Vertex_trackRatioWithPV", VFLOAT);

  //SignalSV
  Branches.makeBranch("SignalSV_nTotal","SignalSV_nTotal", INT);
  Branches.makeBranch("SignalSV_nTracks","SignalSV_nTracks", VINT);
  Branches.makeBranch("SignalSV_trackIndex","SignalSV_trackIndex", VUINT);
  Branches.makeBranch("SignalSV_vertexIndex","SignalSV_vertexIndex", VUINT);
  Branches.makeBranch("SignalSV_trackWeight","SignalSV_trackWeight", VFLOAT);
  Branches.makeBranch("SignalSV_genIndex","SignalSV_genIndex", VUINT);
  Branches.makeBranch("SignalSV_x","SignalSV_x", VFLOAT);
  Branches.makeBranch("SignalSV_y","SignalSV_y", VFLOAT);
  Branches.makeBranch("SignalSV_z","SignalSV_z", VFLOAT);
  Branches.makeBranch("SignalSV_dxy","SignalSV_dxy", VFLOAT);
  Branches.makeBranch("SignalSV_chi2","SignalSV_chi2", VFLOAT);
  Branches.makeBranch("SignalSV_normalizedChi2","SignalSV_normalizedChi2", VFLOAT);
  Branches.makeBranch("SignalSV_ndof","SignalSV_ndof", VFLOAT);

  //Branches.makeBranch("","", INT); 
  Branches.attachBranches(fOutTree);
}

void DisplacedGeneralVertices::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  // Get handles on the input edm collections
  iEvent.getByToken( generalTracksToken_, generalTracksHandle_);
  iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_);
  iEvent.getByToken( ecalTracksToken_, ecalTracksHandle_);
  iEvent.getByToken( electronToken_, electronHandle_);
  iEvent.getByToken( mergedSCsToken_, mergedSCsHandle_);
  iEvent.getByToken( pfCandidatesToken_, pfCandidatesHandle_);
  iEvent.getByToken(genToken_, genHandle_);
  iEvent.getByToken(pvToken_, pvHandle_);
  beamSpot_ = iEvent.get(beamspotToken_);
  primaryVertex_ = iEvent.get(pvToken_).at(0);

  tracks_.clear();
  ttracks_.clear();
  signalTracks_.clear();

  if( DEBUG ) std::cout << "Collecting Examples" << std::endl;

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  const double ptCut = 5;
  
  vertexBuilder_ = VertexAssembly(ttBuilder, generalTracksHandle_, beamSpot_);
  //generalVertices_ = vertexBuilder_.CreateVertexCollection(ptCut);

  int vertexCount(0);
  for(const auto &vtx : generalVertices_)
    vertexCount += int(vtx.tracksSize());

  reco::TrackCollection selectTracks;
  
  for(const auto &track : *generalTracksHandle_) {
    if(track.pt() > ptCut) {
      tracks_.emplace_back(track);
      ttracks_.emplace_back(ttBuilder->build(track));

      bool foundInPV(false);
      for(const auto &pv : *pvHandle_)
	for(const auto &pvTrackRef : pv.tracks())
	  if(TrackHelper::SameTrack(track, *pvTrackRef) && pv.trackWeight(pvTrackRef) > 0.5) {
	    foundInPV = true;
	    break;
	  }

      if(!foundInPV)
	selectTracks.emplace_back(track);
    }
  }

  generalVertices_ = vertexBuilder_.CreateVertexCollection(selectTracks, false); 
  
  TrackPropagator<reco::Track> generalTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, tracks_);
  propECALTracks_ = generalTrackPropagator.GetPropagatedTracks();
  ecalTracks_ = generalTrackPropagator.GetTracks();  

  if(ecalTracks_.size() > 0 && mergedSCsHandle_->size() > 0) {
    MatchTracksToSC<reco::Track> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, ecalTracks_, *mergedSCsHandle_);
    matchedTracksToSCs_ = assigner.GetMatchedTrackSCPairs();
  }
    
  if(cfFlag("hasGenInfo")) {
    int sigCount(0);
    for(const auto &gen : *genHandle_) {
      if(abs(gen.pdgId()) == 11 && gen.status() == 1)
        if(isSignalGenElectron(gen)) sigCount++;
    }

    DeltaRGenMatchHungarian<reco::Track> assigner(tracks_, *genHandle_);
    genMatchedTracks_ = assigner.GetPairedObjects();

    for(const auto &pair : genMatchedTracks_)
      if(isSignalGenElectron(pair.GetObjectB())) {
	signalTracks_.emplace_back(pair.GetObjectA());
      }
  }

}//<<>>void DisplacedGeneralVertices::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void DisplacedGeneralVertices::PostProcessEvent( ItemManager<float>& geVar ) {

  Branches.clearBranches();

  // Fill Beamspot branches
  Branches.fillBranch("Beamspot_x", float(beamSpot_.x0()));
  Branches.fillBranch("Beamspot_y", float(beamSpot_.y0()));
  Branches.fillBranch("Beamspot_z", float(beamSpot_.z0()));

  // Fill track branches
  Branches.fillBranch("Track_nTotal", unsigned(ttracks_.size()) );
  Branches.fillBranch("Track_nECALTracks", unsigned(ecalTracks_.size()) );
  Branches.fillBranch("GeneralTrack_nTotal", unsigned(generalTracksHandle_->size()) );
  for(const auto &ttrack : ttracks_) {
    
    reco::Track track(ttrack.track());
    GlobalVector direction(track.px(), track.py(), track.pz());
    const int ecalTrackIndex(TrackHelper::FindTrackIndex(track, ecalTracks_));
    const int scTrackIndex(TrackHelper::FindTrackIndex(track, matchedTracksToSCs_.GetTracks()));
    //const double scDeltaR(matchedTracksToSCs_[scTrackIndex].GetDeltaR());
    //reco::SuperCluster superCluster(matchedTracksToSCs_[scTrackIndex].GetSuperCluster());
    std::pair<bool, Measurement1D> ip3D = IPTools::signedImpactParameter3D(ttrack, direction, primaryVertex_);
    std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(ttrack, direction, primaryVertex_);

    double scDeltaR(-999.);
    double scEnergy(-999.);
    if(scTrackIndex >= 0) {
      const reco::SuperCluster superCluster(matchedTracksToSCs_[scTrackIndex].GetSuperCluster());
      scDeltaR = matchedTracksToSCs_[scTrackIndex].GetDeltaR();
      scEnergy = superCluster.energy();
    }
    
    Branches.fillBranch("Track_charge", int(track.charge()) );
    Branches.fillBranch("Track_p", float(track.p()) );
    Branches.fillBranch("Track_px", float(track.px()) );
    Branches.fillBranch("Track_py", float(track.py()) );
    Branches.fillBranch("Track_pz", float(track.pz()) );
    Branches.fillBranch("Track_pt", float(track.pt()) );
    Branches.fillBranch("Track_eta", float(track.eta()) );
    Branches.fillBranch("Track_phi", float(track.phi()) );
    Branches.fillBranch("Track_vx", float(track.vx()) );
    Branches.fillBranch("Track_vy", float(track.vy()) );
    Branches.fillBranch("Track_vz", float(track.vz()) );
    Branches.fillBranch("Track_dxy", float(track.dxy()) );
    Branches.fillBranch("Track_ndof", float(track.ndof()) );
    Branches.fillBranch("Track_chi2", float(track.chi2()) );
    Branches.fillBranch("Track_normalizedChi2", float(track.normalizedChi2()) );
    Branches.fillBranch("Track_qoverp", float(track.qoverp()) );
    Branches.fillBranch("Track_ptError", float(track.ptError()) );
    Branches.fillBranch("Track_etaError", float(track.etaError()) );
    Branches.fillBranch("Track_phiError", float(track.phiError()) );
    Branches.fillBranch("Track_ip2D", float(ip2D.second.value()) );
    Branches.fillBranch("Track_sip2D", float(ip2D.second.value() / ip2D.second.error()) );
    Branches.fillBranch("Track_ip3D", float(ip3D.second.value()) );
    Branches.fillBranch("Track_sip3D", float(ip3D.second.value() / ip3D.second.error()) );
    Branches.fillBranch("Track_qualityMask", int(track.qualityMask()) );
    Branches.fillBranch("Track_nValidHits", int(track.numberOfValidHits()) );
    Branches.fillBranch("Track_nLostHits", int(track.numberOfLostHits()) );
    Branches.fillBranch("Track_ecalTrackIndex", int(ecalTrackIndex));
    Branches.fillBranch("Track_isECAL", bool(ecalTrackIndex >= 0));
    Branches.fillBranch("Track_isSCmatch", bool(scTrackIndex >= 0 && scDeltaR < 0.1));
    Branches.fillBranch("Track_SCDR", float(scDeltaR) );
    Branches.fillBranch("Track_energySC", float(scEnergy) );
    Branches.fillBranch("Track_ratioPtToEnergySC", float( fabs(scEnergy - track.p())/ scEnergy) );
  }

  // Fill PV branches
  Branches.fillBranch("PV_nTracks", unsigned(primaryVertex_.tracksSize()) );
  Branches.fillBranch("PV_x", float(primaryVertex_.x()) );
  Branches.fillBranch("PV_y", float(primaryVertex_.y()) );
  Branches.fillBranch("PV_z", float(primaryVertex_.z()) );
  Branches.fillBranch("PV_ecalness", float(VertexHelper::CalculateEcalness(primaryVertex_, ecalTracks_)) );
  Branches.fillBranch("PV_sumPt", float(VertexHelper::CalculateTotalPt(primaryVertex_)) );
  Branches.fillBranch("PV_weightedSumPt", float(VertexHelper::CalculateTotalWeightedPt(primaryVertex_)) );
  Branches.fillBranch("PV_chi2", float(primaryVertex_.chi2()) );
  Branches.fillBranch("PV_ndof", float(primaryVertex_.ndof()) );
  Branches.fillBranch("PV_normalizedChi2", float(primaryVertex_.normalizedChi2()) );

  for(const auto &trackRef : primaryVertex_.tracks()) {
    const reco::Track track(*trackRef);
    const double trackWeight = primaryVertex_.trackWeight(trackRef);
    Branches.fillBranch("PV_trackIndex", unsigned(TrackHelper::FindTrackIndex(track, tracks_)) );
    Branches.fillBranch("PV_trackWeight", float(trackWeight) );
  }

  int uniqueTrackCnt(0);
  std::vector<double> ptVec({-1.});
  for(const auto &pv : *pvHandle_) {
    for(const auto &track : VertexHelper::GetTracks(pv)) {
      for(const auto &pt : ptVec)
	if(pt != track.pt()) {
	  ptVec.push_back(track.pt());
	  uniqueTrackCnt++;
	  break;
	}
	else break;
    }
  }
  Branches.fillBranch("PVCollection_nTracks", unsigned(uniqueTrackCnt));
  
  int sigCounter(0);
  if(cfFlag("hasGenInfo")) {

    Branches.fillBranch("PV_hasSignal", VertexHasSignal(primaryVertex_, sigCounter));
    Branches.fillBranch("PV_signalCount", unsigned(sigCounter));

    Branches.fillBranch("PVCollection_nTotal", unsigned(pvHandle_->size()));
    int pvIndex(0);
    for(const auto &pv : *pvHandle_) {
      for(const auto &matchPair : genMatchedTracks_) {
	const reco::Track track(matchPair.GetObjectA());
	const reco::GenParticle gen(matchPair.GetObjectB());
	if(VertexHelper::isInVertex(pv, track) && isSignalGenElectron(gen)) {
	  const int index = TrackHelper::FindTrackIndex(track, VertexHelper::GetTracks(pv));
	  const double trackWeight = pv.trackWeight(pv.tracks().at(index));
	  Branches.fillBranch("PVCollection_signalWeight", float(trackWeight));
	  Branches.fillBranch("PVCollection_pvIndex", unsigned(pvIndex));
	  Branches.fillBranch("PVCollection_trackIndex", unsigned(TrackHelper::FindTrackIndex(track, tracks_)) );
	  Branches.fillBranch("PVCollection_genIndex", unsigned(TrackHelper::FindIndex(gen, genMatchedTracks_.GetObjectBList())) );
	  Branches.fillBranch("PVCollection_x", float(pv.x()));
	  Branches.fillBranch("PVCollection_y", float(pv.y()));
	  Branches.fillBranch("PVCollection_z", float(pv.z()));
	}
      }
      pvIndex++;
    }
  }
  
  // Fill Vertex branches
  Branches.fillBranch("Vertex_nTotal", unsigned(generalVertices_.size()) );
  int vtxIndex(0);
  for(const auto &vertex : generalVertices_) {
    LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));
    
    Branches.fillBranch("Vertex_nTracks", unsigned(vertex.tracksSize()) );
    Branches.fillBranch("Vertex_hasSignal", VertexHasSignal(vertex, sigCounter));
    Branches.fillBranch("Vertex_signalCount", unsigned(sigCounter));
    Branches.fillBranch("Vertex_x", float(vertex.x()) );
    Branches.fillBranch("Vertex_y", float(vertex.y()) );
    Branches.fillBranch("Vertex_z", float(vertex.z()) );
    Branches.fillBranch("Vertex_p", float(vertex4Vec.P()) );
    Branches.fillBranch("Vertex_px", float(vertex4Vec.px()) );
    Branches.fillBranch("Vertex_py", float(vertex4Vec.py()) );
    Branches.fillBranch("Vertex_pz", float(vertex4Vec.pz()) );
    Branches.fillBranch("Vertex_pt", float(vertex4Vec.pt()) );
    Branches.fillBranch("Vertex_eta", float(vertex4Vec.eta()) );
    Branches.fillBranch("Vertex_phi", float(vertex4Vec.phi()) );
    Branches.fillBranch("Vertex_p", float(vertex4Vec.P()) );    
    Branches.fillBranch("Vertex_mass", float(vertex4Vec.M()) );
    Branches.fillBranch("Vertex_chi2", float(vertex.chi2()) );
    Branches.fillBranch("Vertex_normalizedChi2", float(vertex.normalizedChi2()) );
    Branches.fillBranch("Vertex_ndof", float(vertex.ndof()) );
    Branches.fillBranch("Vertex_dxy", float(sqrt(vertex.x()*vertex.x()+vertex.y()*vertex.y())) );
    Branches.fillBranch("Vertex_ecalness", float(VertexHelper::CalculateEcalness(vertex, ecalTracks_)) );
    Branches.fillBranch("Vertex_sumPt", float(VertexHelper::CalculateTotalPt(vertex)) );
    Branches.fillBranch("Vertex_trackRatioWithPV", float(VertexHelper::CalculateTrackOverlap(vertex, primaryVertex_)) );
    Branches.fillBranch("Vertex_totalTrackWeight", float(VertexHelper::CalculateTotalTrackWeight(vertex)));
    
    for(const auto &trackRef : vertex.tracks()) {
      const reco::Track track(*trackRef);
      const double trackWeight = vertex.trackWeight(trackRef);
      Branches.fillBranch("Vertex_vertexIndex", unsigned(vtxIndex) );
      Branches.fillBranch("Vertex_trackIndex", unsigned(TrackHelper::FindTrackIndex(track, tracks_)) );      
      Branches.fillBranch("Vertex_trackWeight", float(trackWeight) );
    }
    vtxIndex++;
  }

  // Fill branches that depend on gen information
  if(cfFlag("hasGenInfo")) {
  
    // Signal SVs
    reco::VertexCollection signalVertices = vertexBuilder_.CreateVertexCollection(signalTracks_, false);
    
    Branches.fillBranch("SignalSV_nTotal", int(signalVertices.size()) );

    vtxIndex = 0;
    for(reco::Vertex &vertex : signalVertices) {
      
      Branches.fillBranch("SignalSV_x", float(vertex.x()) );
      Branches.fillBranch("SignalSV_y", float(vertex.y()) );
      Branches.fillBranch("SignalSV_z", float(vertex.z()) );
      Branches.fillBranch("SignalSV_dxy", float(sqrt(vertex.x()*vertex.x()+vertex.y()*vertex.y())) );
      Branches.fillBranch("SignalSV_chi2", float(vertex.chi2()) );
      Branches.fillBranch("SignalSV_normalizedChi2", float(vertex.normalizedChi2()) );
      Branches.fillBranch("SignalSV_ndof", float(vertex.ndof()) );
      Branches.fillBranch("SignalSV_nTracks", int(vertex.tracksSize()) );
      
      for(const auto &trackRef : vertex.tracks()) {
	const reco::Track track(*trackRef);
	const int trackIndex = TrackHelper::FindTrackIndex(track, tracks_);
	const int genIndex = TrackHelper::FindIndex(genMatchedTracks_.FindObjectB(track), genMatchedTracks_.GetObjectBList());
	const double trackWeight = vertex.trackWeight(trackRef);
	
	Branches.fillBranch("SignalSV_trackIndex", unsigned(trackIndex) );
	Branches.fillBranch("SignalSV_vertexIndex", unsigned(vtxIndex) );
	Branches.fillBranch("SignalSV_trackWeight", float(trackWeight) );
	Branches.fillBranch("SignalSV_genIndex", unsigned(genIndex));
      }
      vtxIndex++;
    }
    
    // Gen Particle
    sigCounter = 0;
    Branches.fillBranch("GenParticle_nMatches", unsigned(genMatchedTracks_.size()) );

    for(const auto &pair : genMatchedTracks_) {
      const reco::Track track(pair.GetObjectA());
      const reco::GenParticle gen(pair.GetObjectB());

      const bool isSignal(abs(gen.pdgId()) == 11 && isSignalGenElectron(gen));

      Branches.fillBranch("GenParticle_matchedTrackIndex", int(TrackHelper::FindTrackIndex(track, tracks_)));
      Branches.fillBranch("GenParticle_pdgId", int(gen.pdgId()));
      Branches.fillBranch("GenParticle_charge", int(gen.charge()));
      Branches.fillBranch("GenParticle_deltaR", float(pair.GetDeltaR()));
      Branches.fillBranch("GenParticle_pt", float(gen.pt()));
      Branches.fillBranch("GenParticle_eta", float(gen.pt()));
      Branches.fillBranch("GenParticle_phi", float(gen.phi()));
      Branches.fillBranch("GenParticle_x", float(gen.vx()));
      Branches.fillBranch("GenParticle_y", float(gen.vy()));
      Branches.fillBranch("GenParticle_z", float(gen.vz()));
      Branches.fillBranch("GenParticle_dxy", float( sqrt(gen.vx()*gen.vx() + gen.vy()*gen.vy())) );
      Branches.fillBranch("GenParticle_p", float(gen.p()));
      Branches.fillBranch("GenParticle_isSignal", isSignal);
      if(isSignal) sigCounter++;
    }
    Branches.fillBranch("GenParticle_nSignal", unsigned(sigCounter) );
  }
}

bool DisplacedGeneralVertices::VertexHasSignal(const reco::Vertex &vertex, int &counter) const {
  counter = 0;
  bool hasSignal(false);
  for(const auto &track : signalTracks_) {
    if(VertexHelper::isInVertex(vertex, track)) {
      hasSignal = true;
      counter++;
    }
  }

  return hasSignal;
}

#endif
