#ifndef KUCMSDisplacedVertex_HH
#define KUCMSDisplacedVertex_HH
// -*- C++ -*-
//
//
// Original Author: Andres Abreu
//
// KUCMSDisplacedVertex Object
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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

// Add includes for specfic collections
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackPropagator.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchedTrackSCPair.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatch.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TTBuilderWrapper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/MatchTracksToSC.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DisplacedGenZ.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#include <iomanip>
#include <map>
#include <sstream>

using namespace edm;

class KUCMSDisplacedVertexMini : public KUCMSObjectBase {

public:

  // use base class constructor
  KUCMSDisplacedVertexMini( const edm::ParameterSet& iConfig );
  virtual ~KUCMSDisplacedVertexMini() = default;

  void LoadSip2DMuonEnhancedTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ muonEnhancedTracksToken_ = token; }
  void LoadMuonEnhancedTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ muonEnhancedTracksToken_ = token; }
  void LoadMuonTracks( edm::EDGetTokenT<reco::TrackCollection> token ) { muonTracksToken_ = token; }
  void LoadPrimaryVertex( edm::EDGetTokenT<reco::VertexCollection> token ) {pvToken_ = token;}
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadMergedSCs( edm::EDGetTokenT<reco::SuperClusterCollection> token ) {mergedSCsToken_ = token;}
  void LoadGenParticlesToken(edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> token ){ genToken_ = token; }
  void LoadPrunedGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> token){ prunedGenToken_ = token; }
  void LoadLeptonicSVsToken( edm::EDGetTokenT<reco::VertexCollection> token ){ leptonicSVsToken_ = token; }
  void LoadHadronicSVsToken( edm::EDGetTokenT<reco::VertexCollection> token ){ hadronicSVsToken_ = token; }

  void InitObject( TTree* fOutTree );
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  void ProcessEvent( ItemManager<float>& geVar );
  void PostProcessEvent( ItemManager<float>& geVar ) {}
  void EndJobs() {}

  template<class T>
  bool FoundLeptonMatch(const T &lepton) const;

  template<class T>
  bool IsPromptLepton(const T &lepton) const;

  template<class T>
  std::pair<double, double> BestMatch(const T &lepton) const;

private:

  reco::Vertex primaryVertex_;
  TrackVertexSetCollection generalVertices_;
  std::vector<DisplacedGenZ> signalZs_;
  reco::TrackCollection signalTracks_, electronTracks_;

  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  TTBuilderWrapper ttBuilder_;

  std::vector<MatchedTrackSCPair<reco::Track>> trackSCPairs_;

  edm::EDGetTokenT<reco::TrackCollection> muonTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> muonEnhancedTracksToken_;
  edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  edm::EDGetTokenT<reco::VertexCollection> leptonicSVsToken_;
  edm::EDGetTokenT<reco::VertexCollection> hadronicSVsToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> prunedGenToken_;

  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;

  edm::Handle<reco::TrackCollection> muonTracksHandle_;
  edm::Handle<reco::TrackCollection> muonEnhancedTracksHandle_;
  edm::Handle<std::vector<pat::PackedGenParticle>> genHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle_;
  edm::Handle<reco::VertexCollection> leptonicSVsHandle_;
  edm::Handle<reco::VertexCollection> hadronicSVsHandle_;
  edm::Handle<reco::GenParticleCollection> prunedGenHandle_;

  bool IsBronze(const reco::Vertex &vertex) const;
  bool IsSilver(const reco::Vertex &vertex) const;
  bool IsGold(const reco::Vertex &vertex) const;
  int FindGenVertexIndex(const reco::Vertex &vertex) const;
  int FindNearestGenVertexIndex(const reco::Vertex &vertex, double &distance) const;
  bool getSCMatch(const reco::Track &track, reco::SuperCluster &sc, double &deltaR) const;
  // Build a TrackVertexSet from a reco::Vertex using the sip2D track handle for refs
  TrackVertexSet buildTrackVertexSet(const reco::Vertex &vertex, const TransientTrackBuilder* ttBuilder) const;
  // Look up which Z mode a track came from (returns ZDecayMode::Unknown if not a signal track)
  ZDecayMode getZModeFromTrack(const reco::Track &track) const;
  // Print a per-event summary table (only when leptonic/hadronic Zs are present)
  void printEventSummaryTable() const;

};//<<>>class KUCMSDisplacedVertexMini : public KUCMSObjectBase

KUCMSDisplacedVertexMini::KUCMSDisplacedVertexMini( const edm::ParameterSet& iConfig ) :
  caloGeometryToken_(esConsumes()){
  // ---- end constructor initilizations  --------------------------

  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
  cfFlag.set( "verboseEventTable", iConfig.existsAs<bool>("verboseEventTable") ? iConfig.getParameter<bool>("verboseEventTable") : false );

  trackAssociator_.useDefaultPropagator();

}//<<>>KUCMSDisplacedVertexMini::KUCMSDisplacedVertexMini( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSDisplacedVertexMini::InitObject( TTree* fOutTree ){

  // Vertex Branches
  Branches.makeBranch("Vertex_nTotal","Vertex_nTotal", UINT);
  Branches.makeBranch("Vertex_nTracks","Vertex_nTracks", VUINT);
  Branches.makeBranch("Vertex_mass","Vertex_mass", VFLOAT);
  Branches.makeBranch("Vertex_cosTheta","Vertex_cosTheta", VFLOAT);
  Branches.makeBranch("Vertex_decayAngle","Vertex_decayAngle", VFLOAT);
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
  Branches.makeBranch("Vertex_cxx","Vertex_cxx", VFLOAT);
  Branches.makeBranch("Vertex_cyy","Vertex_cyy", VFLOAT);
  Branches.makeBranch("Vertex_czz","Vertex_czz", VFLOAT);
  Branches.makeBranch("Vertex_cxy","Vertex_cxy", VFLOAT);
  Branches.makeBranch("Vertex_cxz","Vertex_cxz", VFLOAT);
  Branches.makeBranch("Vertex_cyz","Vertex_cyz", VFLOAT);
  Branches.makeBranch("Vertex_dxy","Vertex_dxy", VFLOAT);
  Branches.makeBranch("Vertex_dxyError","Vertex_dxyError", VFLOAT);
  Branches.makeBranch("Vertex_chi2","Vertex_chi2", VFLOAT);
  Branches.makeBranch("Vertex_normalizedChi2","Vertex_normalizedChi2", VFLOAT);
  Branches.makeBranch("Vertex_ndof","Vertex_ndof", VFLOAT);
  Branches.makeBranch("Vertex_sumCharge","Vertex_sumCharge", VINT);
  Branches.makeBranch("Vertex_isUnique", "Vertex_isUnique", VBOOL);
  Branches.makeBranch("Vertex_scMatchRatio", "Vertex_scMatchRatio", VFLOAT);
  Branches.makeBranch("Vertex_passLooseMuonID", "Vertex_passLooseMuonID", VBOOL);
  Branches.makeBranch("Vertex_passLooseElectronID", "Vertex_passLooseElectronID", VBOOL);
  Branches.makeBranch("Vertex_passLooseID", "Vertex_passLooseID", VBOOL);

  // VertexTrack Branches
  Branches.makeBranch("VertexTrack_vertexIndex","VertexTrack_vertexIndex", VUINT);
  Branches.makeBranch("VertexTrack_trackIndex","VertexTrack_trackIndex", VUINT);
  Branches.makeBranch("VertexTrack_trackCosTheta","VertexTrack_trackCosTheta", VFLOAT);
  Branches.makeBranch("VertexTrack_trackCosThetaAtCM","VertexTrack_trackCosThetaAtCM", VFLOAT);
  Branches.makeBranch("VertexTrack_trackCompatibility","VertexTrack_trackCompatibility", VFLOAT);
  Branches.makeBranch("VertexTrack_shiftDzAfterTrackRemoval", "VertexTrack_shiftDzAfterTrackRemoval", VFLOAT);
  Branches.makeBranch("VertexTrack_shift3DAfterTrackRemoval", "VertexTrack_shift3DAfterTrackRemoval", VFLOAT);
  Branches.makeBranch("VertexTrack_SCDR", "VertexTrack_SCDR", VFLOAT);
  Branches.makeBranch("VertexTrack_energySC", "VertexTrack_energySC", VFLOAT);
  Branches.makeBranch("VertexTrack_ratioPToEnergySC", "VertexTrack_ratioPToEnergySC", VFLOAT);

  // Gen Vertices
  if(cfFlag("hasGenInfo")) {
    Branches.makeBranch("Vertex_genVertexIndex","Vertex_genVertexIndex", VINT);
    Branches.makeBranch("Vertex_nearestGenVertexIndex","Vertex_nearestGenVertexIndex", VINT);
    Branches.makeBranch("Vertex_min3D", "Vertex_min3D", VFLOAT);
    Branches.makeBranch("Vertex_matchRatio", "Vertex_matchRatio", VFLOAT);
    Branches.makeBranch("Vertex_isBronze","Vertex_isBronze", VBOOL);
    Branches.makeBranch("Vertex_isSilver","Vertex_isSilver", VBOOL);
    Branches.makeBranch("Vertex_isGold","Vertex_isGold", VBOOL);
    Branches.makeBranch("VertexTrack_isSignalTrack", "VertexTrack_isSignalTrack", VBOOL);
    Branches.makeBranch("VertexTrack_isSignalElectron", "VertexTrack_isSignalElectron", VBOOL);
    Branches.makeBranch("VertexTrack_isSignalMuon", "VertexTrack_isSignalMuon", VBOOL);

    Branches.makeBranch("GenVertex_nTotal", "GenVertex_nTotal", UINT);
    Branches.makeBranch("GenVertex_nTracks", "GenVertex_nTracks", VUINT);
    Branches.makeBranch("GenVertex_nElectron", "GenVertex_nElectron", UINT, "amount of electron gen vertices");
    Branches.makeBranch("GenVertex_nMuon", "GenVertex_nMuon", UINT, "amount of muon gen vertices");
    Branches.makeBranch("GenVertex_nHadronic", "GenVertex_nHadronic", UINT, "amount of hadronic gen vertices");
    Branches.makeBranch("GenVertex_mass","GenVertex_mass", VFLOAT);
    Branches.makeBranch("GenVertex_x","GenVertex_x", VFLOAT);
    Branches.makeBranch("GenVertex_y","GenVertex_y", VFLOAT);
    Branches.makeBranch("GenVertex_z","GenVertex_z", VFLOAT);
    Branches.makeBranch("GenVertex_p","GenVertex_p", VFLOAT);
    Branches.makeBranch("GenVertex_px", "GenVertex_px", VFLOAT);
    Branches.makeBranch("GenVertex_py", "GenVertex_py", VFLOAT);
    Branches.makeBranch("GenVertex_pz", "GenVertex_pz", VFLOAT);
    Branches.makeBranch("GenVertex_pt","GenVertex_pt", VFLOAT);
    Branches.makeBranch("GenVertex_eta","GenVertex_eta", VFLOAT);
    Branches.makeBranch("GenVertex_phi","GenVertex_phi", VFLOAT);
    Branches.makeBranch("GenVertex_dxy","GenVertex_dxy", VFLOAT);
    Branches.makeBranch("GenVertex_isElectron","GenVertex_isElectron", VBOOL);
    Branches.makeBranch("GenVertex_isMuon","GenVertex_isMuon", VBOOL);
    Branches.makeBranch("GenVertex_isHadronic","GenVertex_isHadronic", VBOOL);
    Branches.makeBranch("GenVertex_passSelection","GenVertex_passSelection", VBOOL);
    Branches.makeBranch("GenVertex_passSelectionAndCuts","GenVertex_passSelectionAndCuts", VBOOL);
  }

  Branches.attachBranches(fOutTree);

}//<<>>void KUCMSDisplacedVertexMini::InitObject( TTree* fOutTree )

void KUCMSDisplacedVertexMini::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  iEvent.getByToken( muonTracksToken_, muonTracksHandle_);
  iEvent.getByToken( muonEnhancedTracksToken_, muonEnhancedTracksHandle_);
  if(cfFlag("hasGenInfo")) iEvent.getByToken( genToken_, genHandle_);
  iEvent.getByToken( pvToken_, pvHandle_);
  iEvent.getByToken( mergedSCsToken_, mergedSCsHandle_);
  primaryVertex_ = iEvent.get(pvToken_).at(0);

  // Load HYDDRA pre-reconstructed SV collections
  iEvent.getByToken( leptonicSVsToken_, leptonicSVsHandle_);
  iEvent.getByToken( hadronicSVsToken_, hadronicSVsHandle_);

  signalTracks_.clear();
  electronTracks_.clear();
  signalZs_.clear();

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  ttBuilder_ = TTBuilderWrapper(ttBuilder);

  // Build TrackVertexSet collection from the pre-reconstructed HYDDRA SVs.
  generalVertices_.clear();

  for(const auto& vertex : *leptonicSVsHandle_) {
    TrackVertexSet tvs = buildTrackVertexSet(vertex, ttBuilder);
    if(tvs.isValid()) generalVertices_.add(tvs);
    else cout << "An already valid leptonic vertex failed refitting" << endl;
  }

  for(const auto& vertex : *hadronicSVsHandle_) {
    TrackVertexSet tvs = buildTrackVertexSet(vertex, ttBuilder);
    if(tvs.isValid()) generalVertices_.add(tvs);
    else cout << "An already valid hadronic vertex failed refitting" << endl;
  }

  MatchTracksToSC<reco::Track> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, generalVertices_.tracks(), *mergedSCsHandle_);
  trackSCPairs_ = assigner.GetMatchedTrackSCPairs();

  for(const auto &pair : trackSCPairs_) {
    if(pair.GetDeltaR() < 0.04)
      electronTracks_.emplace_back(pair.GetTrack());
  }

  if(cfFlag("hasGenInfo")) {
    iEvent.getByToken(prunedGenToken_, prunedGenHandle_);

    // Build transient tracks from the muon-enhanced track collection
    std::vector<reco::TransientTrack> ttracks;
    for(const auto &track : *muonEnhancedTracksHandle_)
      ttracks.emplace_back(ttBuilder->build(track));

    // Build DisplacedGenZ objects from the pruned+packed gen collections
    signalZs_ = DisplacedGenZ::build(prunedGenHandle_, genHandle_);

    // Match tracks to each Z boson's daughters independently (no cross-Z competition).
    // This mirrors HyddraSVAnalyzer's per-vertex matching and ensures that daughters
    // of one Z cannot steal good track matches from another Z's daughters.
    bool anyDaughters = false;
    for(int i = 0; i < (int)signalZs_.size(); i++) {
      if(!signalZs_[i].isLeptonic() && !signalZs_[i].isHadronic()) continue;
      auto daughters = signalZs_[i].getTrackableDaughters();
      if(daughters.empty()) continue;
      anyDaughters = true;
      DeltaRMatchHungarian<reco::TransientTrack, const pat::PackedGenParticle*>
          matcher(ttracks, daughters);
      signalZs_[i].setMatches(matcher.GetPairedObjects());
    }

    if(!anyDaughters) return;

    // Populate signalTracks_ from all Z match results
    for(const auto& genZ : signalZs_)
      for(const auto& track : genZ.getTracks())
        signalTracks_.emplace_back(track);
  }

}//<<>>void KUCMSDisplacedVertexMini::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSDisplacedVertexMini::ProcessEvent( ItemManager<float>& geVar ){

    Branches.clearBranches();

	auto nGVertices = generalVertices_.size();
    Branches.fillBranch("Vertex_nTotal", unsigned(nGVertices) );
	int nDisSV(0);
    int vtxIndex(0);
    for(const auto &tvertex : generalVertices_) {

      reco::Vertex vertex(tvertex);
      LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));
      const bool passLooseMuonID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, *muonTracksHandle_) == 2);
      const bool passLooseElectronID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, electronTracks_) == 2);
      const bool passLooseID(vertex.tracksSize() < 4? (passLooseMuonID || passLooseElectronID) : true);

	  bool lepsv = ( unsigned(vertex.tracksSize()) == 2 ) && ( passLooseMuonID || passLooseElectronID );
	  bool hadsv = ( unsigned(vertex.tracksSize()) >= 5 ) && ( (float(vertex4Vec.M())/unsigned(vertex.tracksSize()) ) > 1 );
	  if( lepsv || hadsv ) nDisSV++;

      Branches.fillBranch("Vertex_nTracks", unsigned(vertex.tracksSize()) );
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
      Branches.fillBranch("Vertex_mass", float(vertex4Vec.M()) );
      Branches.fillBranch("Vertex_cosTheta", float(tvertex.cosTheta(primaryVertex_)));
      Branches.fillBranch("Vertex_decayAngle", float(VertexHelper::CalculateDecayAngle(vertex)));
      Branches.fillBranch("Vertex_cxx", float(vertex.error().At(0,0)) );
      Branches.fillBranch("Vertex_cyy", float(vertex.error().At(1,1)) );
      Branches.fillBranch("Vertex_czz", float(vertex.error().At(2,2)) );
      Branches.fillBranch("Vertex_cxy", float(vertex.error().At(0,1)) );
      Branches.fillBranch("Vertex_cxz", float(vertex.error().At(0,2)) );
      Branches.fillBranch("Vertex_cyz", float(vertex.error().At(1,2)) );
      Branches.fillBranch("Vertex_chi2", float(vertex.chi2()) );
      Branches.fillBranch("Vertex_normalizedChi2", float(vertex.normalizedChi2()) );
      Branches.fillBranch("Vertex_ndof", float(vertex.ndof()) );
      Branches.fillBranch("Vertex_dxy", float(VertexHelper::CalculateDxy(vertex, primaryVertex_)));
      Branches.fillBranch("Vertex_dxyError", float(VertexHelper::CalculateDxyError(vertex, primaryVertex_)));
      Branches.fillBranch("Vertex_sumCharge", int(VertexHelper::CalculateTotalCharge(vertex)) );
      Branches.fillBranch("Vertex_isUnique", bool(generalVertices_.isVertexUnique(tvertex)));
      Branches.fillBranch("Vertex_scMatchRatio", float(VertexHelper::CountInstances(vertex, electronTracks_)/float(tvertex.size())));
      Branches.fillBranch("Vertex_passLooseMuonID", bool(passLooseMuonID));
      Branches.fillBranch("Vertex_passLooseElectronID", bool(passLooseElectronID));
      Branches.fillBranch("Vertex_passLooseID", bool(passLooseID));

      if(cfFlag("hasGenInfo")) {
	double minDistance;
	int nearestGenVertexIndex(FindNearestGenVertexIndex(vertex, minDistance));
	Branches.fillBranch("Vertex_isBronze", bool(IsBronze(vertex)));
	Branches.fillBranch("Vertex_isSilver", bool(IsSilver(vertex)));
	Branches.fillBranch("Vertex_isGold", bool(IsGold(vertex)));
	Branches.fillBranch("Vertex_genVertexIndex", int(FindGenVertexIndex(vertex)));
	Branches.fillBranch("Vertex_nearestGenVertexIndex", int(nearestGenVertexIndex));
	Branches.fillBranch("Vertex_min3D", float(minDistance));
	Branches.fillBranch("Vertex_matchRatio", float(nearestGenVertexIndex >= 0
	    ? signalZs_[nearestGenVertexIndex].matchRatio(vertex) : -1));
      }//<<>>if(cfFlag("hasGenInfo"))

      for(const auto &trackRef : vertex.tracks()) {
	const reco::TrackRef ref(trackRef.castTo<reco::TrackRef>());
	if(ref.isNull()) continue;
	const reco::Track track(*ref);

	double deltaR(-1.);
	reco::SuperCluster sc;
	const bool isSCMatched(getSCMatch(track, sc, deltaR));
	Branches.fillBranch("VertexTrack_vertexIndex", unsigned(vtxIndex) );
	Branches.fillBranch("VertexTrack_trackIndex", unsigned(ref.key()) );
	Branches.fillBranch("VertexTrack_trackCosTheta", float(TrackHelper::CalculateCosTheta(primaryVertex_, vertex, track)));
	Branches.fillBranch("VertexTrack_trackCosThetaAtCM", float(VertexHelper::CalculateCMCosTheta(vertex, track)));
	Branches.fillBranch("VertexTrack_trackCompatibility", float(tvertex.compatibility(ref)));
	Branches.fillBranch("VertexTrack_shiftDzAfterTrackRemoval", float(tvertex.shiftDzAfterTrackRemoval(ref)));
	Branches.fillBranch("VertexTrack_shift3DAfterTrackRemoval", float(tvertex.shift3DAfterTrackRemoval(ref)));
	Branches.fillBranch("VertexTrack_SCDR", float(deltaR));
	Branches.fillBranch("VertexTrack_energySC", float(isSCMatched? sc.correctedEnergy() : -1.));
	Branches.fillBranch("VertexTrack_ratioPToEnergySC", float(isSCMatched? track.p()/sc.correctedEnergy() : -1.));

	if(cfFlag("hasGenInfo")) {
	  const bool isSignal(TrackHelper::FindTrackIndex(track, signalTracks_) >= 0);
	  Branches.fillBranch("VertexTrack_isSignalTrack", isSignal);
	  const ZDecayMode zMode = isSignal ? getZModeFromTrack(track) : ZDecayMode::Unknown;
	  Branches.fillBranch("VertexTrack_isSignalElectron", zMode == ZDecayMode::Electron);
	  Branches.fillBranch("VertexTrack_isSignalMuon", zMode == ZDecayMode::Muon);
	}//<<>>if(cfFlag("hasGenInfo"))
      }//<<>>for(const auto &trackRef : vertex.tracks())
      vtxIndex++;
    }//<<>>for(const auto &tvertex : generalVertices_)
    geVar.set("nDisSVs",nDisSV);

    if(cfFlag("hasGenInfo")) {

      int nSigElectrons(0), nSigMuons(0), nSigHadrons(0);
      Branches.fillBranch("GenVertex_nTotal", unsigned(signalZs_.size()));
      for(const auto& genZ : signalZs_) {

        bool passSelectionAndCuts = false;
        if(genZ.hasTracks()) {
          auto checkVertex = [&](const reco::Vertex& recoVtx) {
            for(const auto& pair : genZ.getMatches())
              if(VertexHelper::isInVertex(recoVtx, pair.GetObjectA().track())) return true;
            return false;
          };
          for(const auto& vtx : *leptonicSVsHandle_)
            if(checkVertex(vtx)) { passSelectionAndCuts = true; break; }
          if(!passSelectionAndCuts)
            for(const auto& vtx : *hadronicSVsHandle_)
              if(checkVertex(vtx)) { passSelectionAndCuts = true; break; }
        }//<<>>if(genZ.hasTracks())

        Branches.fillBranch("GenVertex_nTracks", unsigned(genZ.getMatches().size()));
        Branches.fillBranch("GenVertex_mass", float(genZ.mass()));
        Branches.fillBranch("GenVertex_x", float(genZ.x()));
        Branches.fillBranch("GenVertex_y", float(genZ.y()));
        Branches.fillBranch("GenVertex_z", float(genZ.z()));
        Branches.fillBranch("GenVertex_p", float(genZ.p()));
        Branches.fillBranch("GenVertex_px", float(genZ.px()));
        Branches.fillBranch("GenVertex_py", float(genZ.py()));
        Branches.fillBranch("GenVertex_pz", float(genZ.pz()));
        Branches.fillBranch("GenVertex_pt", float(genZ.pt()));
        Branches.fillBranch("GenVertex_eta", float(genZ.eta()));
        Branches.fillBranch("GenVertex_phi", float(genZ.phi()));
        Branches.fillBranch("GenVertex_dxy", float(genZ.Lxy()));
        Branches.fillBranch("GenVertex_isElectron", bool(genZ.mode() == ZDecayMode::Electron));
        Branches.fillBranch("GenVertex_isMuon", bool(genZ.mode() == ZDecayMode::Muon));
        Branches.fillBranch("GenVertex_isHadronic", bool(genZ.isHadronic()));
        Branches.fillBranch("GenVertex_passSelection", bool(genZ.hasTracks()));
        Branches.fillBranch("GenVertex_passSelectionAndCuts", bool(passSelectionAndCuts));

        if(genZ.mode() == ZDecayMode::Electron) nSigElectrons++;
        if(genZ.mode() == ZDecayMode::Muon) nSigMuons++;
        if(genZ.isHadronic()) nSigHadrons++;

      }//<<>>for(const auto& genZ : signalZs_)
      Branches.fillBranch("GenVertex_nElectron", unsigned(nSigElectrons));
      Branches.fillBranch("GenVertex_nMuon", unsigned(nSigMuons));
      Branches.fillBranch("GenVertex_nHadronic", unsigned(nSigHadrons));

      if(cfFlag("verboseEventTable")) printEventSummaryTable();

    }//<<>>if(cfFlag("hasGenInfo"))

}//<<>>void KUCMSDisplacedVertexMini::ProcessEvent()


void KUCMSDisplacedVertexMini::printEventSummaryTable() const {

  // Tally signal Z types (only print when there is at least one leptonic or hadronic Z)
  int nElec(0), nMuon(0), nHad(0);
  for(const auto& z : signalZs_) {
    if(z.mode() == ZDecayMode::Electron)  nElec++;
    else if(z.mode() == ZDecayMode::Muon) nMuon++;
    else if(z.isHadronic())               nHad++;
  }
  if(nElec + nMuon + nHad == 0) return;

  const int nInputTracks = (int)muonEnhancedTracksHandle_->size();
  const int nLepSVs      = (int)leptonicSVsHandle_->size();
  const int nHadSVs      = (int)hadronicSVsHandle_->size();
  if(nLepSVs + nHadSVs == 0) return;

  // Collect per-Z row data
  struct ZRow {
    int         idx;
    std::string type;
    int         chargedDaus;
    int         inSV;
    std::string quality;
    double      lxy;
  };

  std::vector<ZRow> rows;
  int zIdx = 0;
  for(const auto& z : signalZs_) {
    if(!z.isLeptonic() && !z.isHadronic()) { zIdx++; continue; }

    ZRow r;
    r.idx        = zIdx++;
    r.chargedDaus = (int)z.getTrackableDaughters().size();

    if(z.mode() == ZDecayMode::Electron)  r.type = "Electron";
    else if(z.mode() == ZDecayMode::Muon) r.type = "Muon";
    else                                   r.type = "Hadronic";

    // Count how many of this Z's matched tracks appear in any reco SV
    r.inSV = 0;
    for(const auto& pair : z.getMatches()) {
      const reco::Track& trk = pair.GetObjectA().track();
      bool found = false;
      for(const auto& vtx : *leptonicSVsHandle_)
        if(VertexHelper::isInVertex(vtx, trk)) { found = true; break; }
      if(!found)
        for(const auto& vtx : *hadronicSVsHandle_)
          if(VertexHelper::isInVertex(vtx, trk)) { found = true; break; }
      if(found) r.inSV++;
    }

    // Quality: leptonic -> GOLD/SILVER/BRONZE, hadronic -> best match ratio
    if(z.isLeptonic()) {
      bool gold(false), silver(false), bronze(false);
      for(const auto& vtx : *leptonicSVsHandle_) {
        if(z.isGold(vtx))        { gold   = true; break; }
        if(z.isSilver(vtx))        silver = true;
        if(z.isBronze(vtx))        bronze = true;
      }
      if(gold)        r.quality = "GOLD";
      else if(silver) r.quality = "SILVER";
      else if(bronze) r.quality = "BRONZE";
      else            r.quality = "---";
    } else {
      double bestRatio = -1.;
      for(const auto& vtx : *hadronicSVsHandle_) {
        const double ratio = z.matchRatio(vtx);
        if(ratio > bestRatio) bestRatio = ratio;
      }
      if(bestRatio >= 0.) {
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3) << bestRatio;
        r.quality = "ratio: " + oss.str();
      } else {
        r.quality = "---";
      }
    }
    r.lxy = z.Lxy();
    rows.emplace_back(r);
  }

  // ── Layout ───────────────────────────────────────────────────────────────
  // Column content widths (no padding)
  const int cw0 = 4;   // Z#
  const int cw1 = 10;  // Type     ("Hadronic" = 8)
  const int cw2 = 13;  // Charged  ("Charged Daus" = 12)
  const int cw3 = 9;   // In SV    ("In SV" = 5, data "99/99" = 5)
  const int cw4 = 18;  // Quality  ("ratio: 0.750" = 12)
  const int cw5 = 9;   // Lxy      ("Lxy [cm]" = 8)
  // Total table width: 7 borders + 6*(col+2 padding)
  const int tw = 7 + (cw0+2) + (cw1+2) + (cw2+2) + (cw3+2) + (cw4+2) + (cw5+2); // = 83

  // Helpers
  auto centered = [](const std::string& s, int w) -> std::string {
    int sp = w - (int)s.size();
    if(sp <= 0) return s.substr(0, w);
    return std::string(sp/2, ' ') + s + std::string(sp - sp/2, ' ');
  };
  auto lpad = [](const std::string& s, int w) -> std::string {
    if((int)s.size() >= w) return s.substr(0, w);
    return s + std::string(w - (int)s.size(), ' ');
  };
  auto hline = [&](const std::string& l, const std::string& m, const std::string& r) -> std::string {
    return l + std::string(cw0+2,'-') + m + std::string(cw1+2,'-') + m
             + std::string(cw2+2,'-') + m + std::string(cw3+2,'-') + m
             + std::string(cw4+2,'-') + m + std::string(cw5+2,'-') + r;
  };
  auto dataRow = [&](const std::string& s0, const std::string& s1,
                     const std::string& s2, const std::string& s3,
                     const std::string& s4, const std::string& s5) -> std::string {
    return "| " + centered(s0,cw0) + " | " + lpad(s1,cw1) + " | "
               + centered(s2,cw2) + " | " + centered(s3,cw3) + " | "
               + lpad(s4,cw4) + " | " + centered(s5,cw5) + " |";
  };
  // Full-width banner line (spans all columns)
  auto banner = [&](const std::string& s) -> std::string {
    const int inner = tw - 2;
    if((int)s.size() >= inner) return "|" + s.substr(0, inner) + "|";
    return "|" + s + std::string(inner - (int)s.size(), ' ') + "|";
  };
  const std::string boxTop = "+" + std::string(tw-2, '-') + "+";

  // ── Type summary string ───────────────────────────────────────────────────
  std::string types;
  if(nElec) types += std::to_string(nElec) + " electron";
  if(nMuon) { if(!types.empty()) types += ", "; types += std::to_string(nMuon) + " muon"; }
  if(nHad)  { if(!types.empty()) types += ", "; types += std::to_string(nHad)  + " hadronic"; }

  // ── Print ─────────────────────────────────────────────────────────────────
  cout << "\n";
  cout << boxTop << "\n";
  cout << banner("  DISPLACED VERTEX EVENT SUMMARY") << "\n";
  cout << banner("") << "\n";
  cout << banner("  Input tracks : " + std::to_string(nInputTracks)) << "\n";
  cout << banner("  Z bosons found: " + std::to_string(nElec+nMuon+nHad)
                 + "  (" + types + ")") << "\n";
  cout << banner("  Reco SVs: " + std::to_string(nLepSVs+nHadSVs) + " total"
                 + "   [ " + std::to_string(nLepSVs) + " leptonic"
                 + "  |  " + std::to_string(nHadSVs) + " hadronic ]") << "\n";
  cout << hline("+", "+", "+") << "\n";
  cout << dataRow("Z #", "Type", "Charged Daus", "In SV", "Quality", "Lxy [cm]") << "\n";
  cout << hline("+", "+", "+") << "\n";
  for(const auto& r : rows) {
    const std::string inSvStr = std::to_string(r.inSV) + "/" + std::to_string(r.chargedDaus);
    std::ostringstream lxyOss;
    lxyOss << std::fixed << std::setprecision(2) << r.lxy;
    cout << dataRow(std::to_string(r.idx), r.type,
                    std::to_string(r.chargedDaus), inSvStr, r.quality, lxyOss.str()) << "\n";
  }
  cout << hline("+", "+", "+") << "\n\n";

}//<<>>void KUCMSDisplacedVertexMini::printEventSummaryTable()


TrackVertexSet KUCMSDisplacedVertexMini::buildTrackVertexSet(const reco::Vertex &vertex, const TransientTrackBuilder* ttBuilder) const {

  std::vector<reco::TrackRef> trackRefs;
  for(auto it = vertex.tracks_begin(); it != vertex.tracks_end(); ++it) {
    reco::TrackRef ref = it->castTo<reco::TrackRef>();
    if(ref.isNonnull()) trackRefs.emplace_back(ref);
    else
      cout << "This is fucked" << endl;
  }
  return TrackVertexSet(trackRefs, ttBuilder);
}

bool KUCMSDisplacedVertexMini::IsGold(const reco::Vertex &vertex) const {

  for(const auto& genZ : signalZs_)
    if(genZ.isGold(vertex)) return true;
  return false;
}

bool KUCMSDisplacedVertexMini::IsSilver(const reco::Vertex &vertex) const {

  for(const auto& genZ : signalZs_)
    if(genZ.isSilver(vertex)) return true;
  return false;
}

bool KUCMSDisplacedVertexMini::IsBronze(const reco::Vertex &vertex) const {

  for(const auto& genZ : signalZs_)
    if(genZ.isBronze(vertex)) return true;
  return false;
}

int KUCMSDisplacedVertexMini::FindGenVertexIndex(const reco::Vertex &vertex) const {

  int index(0);
  for(const auto& genZ : signalZs_) {
    if(genZ.hasTracks() && (genZ.isBronze(vertex) || genZ.isSilver(vertex) || genZ.isGold(vertex)))
      return index;
    index++;
  }
  return -1;
}

int KUCMSDisplacedVertexMini::FindNearestGenVertexIndex(const reco::Vertex &vertex, double &distance) const {

  int index(0), minIndex(-1);
  distance = std::numeric_limits<double>::max();
  for(const auto& genZ : signalZs_) {
    if(!genZ.isHadronic()) { index++; continue; }
    const double d(genZ.distance3D(vertex));
    if(d < distance) { distance = d; minIndex = index; }
    index++;
  }
  if(minIndex < 0) distance = -1.;
  return minIndex;
}

bool KUCMSDisplacedVertexMini::getSCMatch(const reco::Track &track, reco::SuperCluster &sc, double &deltaR) const {

  for(const auto &pair : trackSCPairs_) {
    if(TrackHelper::SameTrack(track, pair.GetTrack())) {
      sc = pair.GetSuperCluster();
      deltaR = pair.GetDeltaR();
      return true;
    }
  }
  return false;
}

ZDecayMode KUCMSDisplacedVertexMini::getZModeFromTrack(const reco::Track &track) const {

  for(const auto& genZ : signalZs_)
    for(const auto& pair : genZ.getMatches())
      if(TrackHelper::SameTrack(pair.GetObjectA().track(), track))
        return genZ.mode();
  return ZDecayMode::Unknown;
}

template<class T>
bool KUCMSDisplacedVertexMini::FoundLeptonMatch(const T &lepton) const {

  for(const auto &sv : generalVertices_) {
    if(sv.size() > 2)
      continue;

    for(const auto &track : sv.tracks()) {
      const double deltaR = sqrt(DeltaR2(lepton, *track));
      if(deltaR < 0.02) return true;
    }
  }
  return false;
}

template<class T>
bool KUCMSDisplacedVertexMini::IsPromptLepton(const T &lepton) const {
    double deltaRCut = 0.02;
    double relPtDiffCut = 0.01;

    if constexpr (std::is_same_v<T, reco::Electron>) {
        deltaRCut = 0.01;
        relPtDiffCut = 0.05;
    } else if constexpr (std::is_same_v<T, reco::Muon>) {
        deltaRCut = 0.001;
        relPtDiffCut = 0.005;
    }

    for(const auto &track : *muonEnhancedTracksHandle_) {
        const double deltaR = sqrt(DeltaR2(lepton, track));
        const double trackPt(track.pt()), leptonPt(lepton.pt());
        const double relPtDiffMin = 0.5*(trackPt-leptonPt)/(trackPt+leptonPt);
        if(deltaR < deltaRCut && fabs(relPtDiffMin) < relPtDiffCut) {
            return false;
        }
    }

    return true;
}

template<class T>
std::pair<double, double> KUCMSDisplacedVertexMini::BestMatch(const T &lepton) const {

  double bestDeltaRMatch(10), relPtDiffMin(-1);
  for(const auto &track : *muonEnhancedTracksHandle_) {
    const double deltaR = sqrt(DeltaR2(lepton, track));
    if(deltaR < bestDeltaRMatch) {
      const double trackPt(track.pt()), leptonPt(lepton.pt());
      relPtDiffMin = 0.5*(trackPt-leptonPt)/(trackPt+leptonPt);
      bestDeltaRMatch = deltaR;
    }
  }

  return std::make_pair(bestDeltaRMatch, relPtDiffMin);
}

#endif
