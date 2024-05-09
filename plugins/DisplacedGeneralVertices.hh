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
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackTools.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/IsolationInfo.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedTrackSCPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexAssembly.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"

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
  std::vector<PairedObjects<reco::Track, reco::GenParticle> > genMatchedTracks_;

  edm::EDGetTokenT<reco::TrackCollection> generalTracksToken_;
  edm::Handle<reco::TrackCollection> generalTracksHandle_;

  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTracksToken_;
  edm::Handle<reco::GsfTrackCollection> gsfTracksHandle_;

  edm::EDGetTokenT<reco::TrackCollection> ecalTracksToken_;
  edm::Handle<reco::TrackCollection> ecalTracksHandle_;

  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;
  edm::Handle<reco::ElectronCollection> electronHandle_;

  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandidatesToken_;
  edm::Handle<reco::PFCandidateCollection> pfCandidatesHandle_;

  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::Handle<reco::GenParticleCollection> genHandle_;

  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  reco::BeamSpot beamSpot_;

  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
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
  
  // Track variables
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

  // Gen Particle Variables
  Branches.makeBranch("GenParticle_nMatches","GenParticle_nMatches", UINT);
  Branches.makeBranch("GenParticle_nSignal","GenParticle_nSignal", UINT);
  Branches.makeBranch("GenParticle_pdgId","GenParticle_pdgId", VINT);
  Branches.makeBranch("GenParticle_matchedTrackIndex","GenParticle_matchedTrackIndex", VINT);
  Branches.makeBranch("GenParticle_charge","GenParticle_charge", VINT);
  Branches.makeBranch("GenParticle_pt","GenParticle_pt", VFLOAT);
  Branches.makeBranch("GenParticle_eta","GenParticle_eta", VFLOAT);
  Branches.makeBranch("GenParticle_phi","GenParticle_phi", VFLOAT);
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
  Branches.makeBranch("PV_chi2","PV_chi2", FLOAT);
  Branches.makeBranch("PV_ndof","PV_ndof", FLOAT);
  Branches.makeBranch("PV_normalizedChi2","PV_normalizedChi2", FLOAT);
  Branches.makeBranch("PV_ecalness","PV_ecalness", FLOAT);

  // Vertex Variables
  Branches.makeBranch("Vertex_nTotal","Vertex_nTotal", UINT);
  Branches.makeBranch("Vertex_signalCount","Vertex_signalCount", VUINT);
  Branches.makeBranch("Vertex_hasSignal","Vertex_hasSignal", VBOOL);
  Branches.makeBranch("Vertex_nTracks","Vertex_nTracks", VUINT);
  Branches.makeBranch("Vertex_x","Vertex_x", VFLOAT);
  Branches.makeBranch("Vertex_y","Vertex_y", VFLOAT);
  Branches.makeBranch("Vertex_z","Vertex_z", VFLOAT);
  Branches.makeBranch("Vertex_dxy","Vertex_dxy", VFLOAT);
  Branches.makeBranch("Vertex_chi2","Vertex_chi2", VFLOAT);
  Branches.makeBranch("Vertex_normalizedChi2","Vertex_normalizedChi2", VFLOAT);
  Branches.makeBranch("Vertex_ndof","Vertex_ndof", VFLOAT);
  Branches.makeBranch("Vertex_ecalness","Vertex_ecalness", VFLOAT);
  Branches.makeBranch("Vertex_sumPt","Vertex_sumPt", VFLOAT);
  Branches.makeBranch("Vertex_vertexIndex","Vertex_vertexIndex", VUINT);
  Branches.makeBranch("Vertex_trackIndex","Vertex_trackIndex", VUINT);
  Branches.makeBranch("Vertex_trackWeight","Vertex_trackWeight", VFLOAT);
  Branches.makeBranch("Vertex_trackRatioWithPV","Vertex_trackRatioWithPV", VFLOAT);

  //SignalSV
  Branches.makeBranch("SignalSV_nTotal","SignalSV_nTotal", INT);
  Branches.makeBranch("SignalSV_nTracks","SignalSV_nTracks", VINT);
  Branches.makeBranch("SignalSV_trackIndex","SignalSV_trackIndex", VUINT);
  Branches.makeBranch("SignalSV_vertexIndex","SignalSV_vertexIndex", VUINT);
  Branches.makeBranch("SignalSV_trackWeight","SignalSV_trackWeight", VFLOAT);
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
  iEvent.getByToken( pfCandidatesToken_, pfCandidatesHandle_);
  iEvent.getByToken(genToken_, genHandle_);
  beamSpot_ = iEvent.get(beamspotToken_);
  primaryVertex_ = iEvent.get(pvToken_).at(0);

  tracks_.clear();
  ttracks_.clear();
  signalTracks_.clear();

  if( DEBUG ) std::cout << "Collecting Examples" << std::endl;

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  const double ptCut = 0.95;

  vertexBuilder_ = VertexAssembly(ttBuilder, generalTracksHandle_, beamSpot_);
  generalVertices_ = vertexBuilder_.CreateVertexCollection(ptCut);

  int vertexCount(0);
  for(const auto &vtx : generalVertices_)
    vertexCount += int(vtx.tracksSize());

  for(const auto &track : *generalTracksHandle_)
    if(track.pt() > 0.95) {
      tracks_.emplace_back(track);
      ttracks_.emplace_back(ttBuilder->build(track));
    }

  TrackPropagator<reco::Track> generalTrackPropagator(iEvent, iSetup, magfield, trackAssocParameters_, tracks_);
  propECALTracks_ = generalTrackPropagator.GetPropagatedTracks();
  ecalTracks_ = generalTrackPropagator.GetTracks();  

  if(cfFlag("hasGenInfo")) {
    int sigCount(0);
    for(const auto &gen : *genHandle_) {
      if(abs(gen.pdgId()) == 11 && gen.status() == 1)
        if(isSignalGenElectron(gen))
          sigCount++;
    }

    DeltaRGenMatchHungarian<reco::Track> assigner(ecalTracks_, *genHandle_);
    genMatchedTracks_ = assigner.GetPairedObjects();

    for(const auto &pair : genMatchedTracks_)
      if(isSignalGenElectron(pair.GetObjectB()))
	signalTracks_.emplace_back(pair.GetObjectA());
    
  }

}//<<>>void DisplacedGeneralVertices::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void DisplacedGeneralVertices::PostProcessEvent( ItemManager<float>& geVar ) {

  Branches.clearBranches();

  // Fill track branches
  Branches.fillBranch("Track_nTotal", unsigned(ttracks_.size()) );
  Branches.fillBranch("Track_nECALTracks", unsigned(ecalTracks_.size()) );
  for(const auto &ttrack : ttracks_) {
    
    reco::Track track(ttrack.track());
    GlobalVector direction(track.px(), track.py(), track.pz());
    std::pair<bool, Measurement1D> ip3D = IPTools::signedImpactParameter3D(ttrack, direction, primaryVertex_);
    std::pair<bool, Measurement1D> ip2D = IPTools::signedTransverseImpactParameter(ttrack, direction, primaryVertex_);

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
    Branches.fillBranch("Track_ecalTrackIndex", int(FindTrackIndex(track, ecalTracks_)) );
    Branches.fillBranch("Track_isECAL", bool(FindTrackIndex(track, ecalTracks_) >= 0) );
  }

  // Fill PV branches
  int sigCounter(0);
  Branches.fillBranch("PV_nTracks", unsigned(primaryVertex_.tracksSize()) );
  Branches.fillBranch("PV_hasSignal", VertexHasSignal(primaryVertex_, sigCounter));
  Branches.fillBranch("PV_signalCount", unsigned(sigCounter));
  Branches.fillBranch("PV_x", float(primaryVertex_.x()) );
  Branches.fillBranch("PV_y", float(primaryVertex_.y()) );
  Branches.fillBranch("PV_z", float(primaryVertex_.z()) );
  Branches.fillBranch("PV_ecalness", float(VertexHelper::CalculateEcalness(primaryVertex_, ecalTracks_)) );
  Branches.fillBranch("PV_sumPt", float(VertexHelper::CalculateTotalPt(primaryVertex_)) );
  Branches.fillBranch("PV_chi2", float(primaryVertex_.chi2()) );
  Branches.fillBranch("PV_ndof", float(primaryVertex_.ndof()) );
  Branches.fillBranch("PV_normalizedChi2", float(primaryVertex_.normalizedChi2()) );

  // Fill Vertex branches
  Branches.fillBranch("Vertex_nTotal", unsigned(generalVertices_.size()) );
  int vtxIndex(0);
  for(const auto &vertex : generalVertices_) {
    Branches.fillBranch("Vertex_nTracks", unsigned(vertex.tracksSize()) );
    Branches.fillBranch("Vertex_hasSignal", VertexHasSignal(vertex, sigCounter));
    Branches.fillBranch("Vertex_signalCount", unsigned(sigCounter));
    Branches.fillBranch("Vertex_x", float(vertex.x()) );
    Branches.fillBranch("Vertex_y", float(vertex.y()) );
    Branches.fillBranch("Vertex_z", float(vertex.z()) );
    Branches.fillBranch("Vertex_chi2", float(vertex.chi2()) );
    Branches.fillBranch("Vertex_normalizedChi2", float(vertex.normalizedChi2()) );
    Branches.fillBranch("Vertex_ndof", float(vertex.ndof()) );
    Branches.fillBranch("Vertex_dxy", float(sqrt(vertex.x()*vertex.x()+vertex.y()*vertex.y())) );
    Branches.fillBranch("Vertex_ecalness", float(VertexHelper::CalculateEcalness(vertex, ecalTracks_)) );
    Branches.fillBranch("Vertex_sumPt", float(VertexHelper::CalculateTotalPt(vertex)) );
    Branches.fillBranch("Vertex_trackRatioWithPV", float(VertexHelper::CalculateTrackOverlap(vertex, primaryVertex_)) );

    for(const auto &trackRef : vertex.tracks()) {
      const reco::Track track(*trackRef);
      const int trackIndex = FindTrackIndex(track, *generalTracksHandle_);
      const double trackWeight = vertex.trackWeight(reco::TrackRef(generalTracksHandle_, trackIndex));
      Branches.fillBranch("Vertex_vertexIndex", unsigned(vtxIndex) );
      Branches.fillBranch("Vertex_trackIndex", unsigned(FindTrackIndex(track, tracks_)) );      
      Branches.fillBranch("Vertex_trackWeight", float(trackWeight) );
    }
    vtxIndex++;
  }

  // Signal SVs
  reco::VertexCollection signalVertices = vertexBuilder_.CreateVertexCollection(signalTracks_);

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
      const int trackIndex = FindTrackIndex(track, *generalTracksHandle_);
      const double trackWeight = vertex.trackWeight(reco::TrackRef(generalTracksHandle_, trackIndex));
      Branches.fillBranch("SignalSV_trackIndex", unsigned(trackIndex) );
      Branches.fillBranch("SignalSV_vertexIndex", unsigned(vtxIndex) );
      Branches.fillBranch("SignalSV_trackWeight", float(trackWeight) );
    }
    vtxIndex++;
  }

  // Fill gen branches 
  if(cfFlag("hasGenInfo")) {
    
    sigCounter = 0;
    Branches.fillBranch("GenParticle_nMatches", unsigned(genMatchedTracks_.size()) );

    for(const auto &pair : genMatchedTracks_) {
      const reco::Track track(pair.GetObjectA());
      const reco::GenParticle gen(pair.GetObjectB());
      const bool isSignal(abs(gen.pdgId()) == 11 && isSignalGenElectron(gen));

      Branches.fillBranch("GenParticle_matchedTrackIndex", int(FindTrackIndex(track, tracks_)));
      Branches.fillBranch("GenParticle_pdgId", int(gen.pdgId()));
      Branches.fillBranch("GenParticle_charge", int(gen.charge()));
      Branches.fillBranch("GenParticle_pt", float(gen.pt()));
      Branches.fillBranch("GenParticle_eta", float(gen.pt()));
      Branches.fillBranch("GenParticle_phi", float(gen.phi()));
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
