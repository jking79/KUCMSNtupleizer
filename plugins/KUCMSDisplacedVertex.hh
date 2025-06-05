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
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/NewVertexMerger.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/GenVertex.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackVertexSet.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TTBuilderWrapper.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

using namespace edm; 

class KUCMSDisplacedVertex : public KUCMSObjectBase {

public:
  
  // use base class constructor
  KUCMSDisplacedVertex( const edm::ParameterSet& iConfig );
  virtual ~KUCMSDisplacedVertex() = default;

  void LoadMuonEnhancedTracksToken( edm::EDGetTokenT<reco::TrackCollection> token ){ muonEnhancedTracksToken_ = token; }
  void LoadMuonTracks( edm::EDGetTokenT<reco::TrackCollection> token ) { muonTracksToken_ = token; }
  void LoadPrimaryVertex( edm::EDGetTokenT<reco::VertexCollection> token ) {pvToken_ = token;}
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
  void LoadMagneticField( edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> token){magneticFieldToken_ = token; }
  void LoadAssociationParameters(  TrackAssociatorParameters parameters){ trackAssocParameters_ = parameters;}
  void LoadMergedSCs( edm::EDGetTokenT<reco::SuperClusterCollection> token ) {mergedSCsToken_ = token;}
  void LoadGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> token ){ genToken_ = token; }

  void InitObject( TTree* fOutTree ); 
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  void ProcessEvent( ItemManager<float>& geVar );
  void PostProcessEvent( ItemManager<float>& geVar ) {}
  void EndJobs() {}

private:

  reco::Vertex primaryVertex_;
  TrackVertexSetCollection generalVertices_;
  GenVertices genElectronVertices_, genMuonVertices_, genVertices_;
  std::map<GenVertex, GenMatches> chargedMatches_;
  reco::TrackCollection signalTracks_, electronTracks_;
  
  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters trackAssocParameters_;

  TTBuilderWrapper ttBuilder_;

  std::vector<MatchedTrackSCPair<reco::Track>> trackSCPairs_;

  edm::EDGetTokenT<reco::TrackCollection> muonTracksToken_;
  edm::EDGetTokenT<reco::TrackCollection> muonEnhancedTracksToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::EDGetTokenT<reco::VertexCollection> pvToken_;
  edm::EDGetTokenT<reco::SuperClusterCollection> mergedSCsToken_;
  
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticFieldToken_;
  
  edm::Handle<reco::TrackCollection> muonTracksHandle_;
  edm::Handle<reco::TrackCollection> muonEnhancedTracksHandle_;
  edm::Handle<reco::GenParticleCollection> genHandle_;
  edm::Handle<reco::VertexCollection> pvHandle_;
  edm::Handle<reco::SuperClusterCollection> mergedSCsHandle_;

  bool IsBronze(const reco::Vertex &vertex) const;
  bool IsSilver(const reco::Vertex &vertex) const;
  bool IsGold(const reco::Vertex &vertex) const;
  int FindGenVertexIndex(const reco::Vertex &vertex) const;
  int FindNearestGenVertexIndex(const reco::Vertex &vertex, double &distance) const;
  std::vector<reco::TrackRef> ConvertTracksToRefs(const edm::Handle<reco::TrackCollection> &trackHandle) const;
  std::vector<reco::TrackRef> ConvertTracksToRefs(const reco::TrackCollection &tracks, const edm::Handle<reco::TrackCollection> &trackHandle) const;
  bool getSCMatch(const reco::Track &track, reco::SuperCluster &sc, double &deltaR) const;
  double matchRatio(const TrackVertexSet &vertex, const GenVertex &genVertex, const double threshold = 0.1) const;
  
};//<<>>class KUCMSDisplacedVertex : public KUCMSObjectBase

KUCMSDisplacedVertex::KUCMSDisplacedVertex( const edm::ParameterSet& iConfig ) :
  caloGeometryToken_(esConsumes()){   
  // ---- end constructor initilizations  --------------------------
  
  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );

  trackAssociator_.useDefaultPropagator();
  
}//<<>>KUCMSDisplacedVertex::KUCMSDisplacedVertex( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSDisplacedVertex::InitObject( TTree* fOutTree ){

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

}//<<>>void KUCMSDisplacedVertex::InitObject( TTree* fOutTree )

void KUCMSDisplacedVertex::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  iEvent.getByToken( muonTracksToken_, muonTracksHandle_);
  iEvent.getByToken( muonEnhancedTracksToken_, muonEnhancedTracksHandle_);
  if(cfFlag("hasGenInfo")) iEvent.getByToken( genToken_, genHandle_);
  iEvent.getByToken( pvToken_, pvHandle_);
  iEvent.getByToken( mergedSCsToken_, mergedSCsHandle_);
  primaryVertex_ = iEvent.get(pvToken_).at(0);

  signalTracks_.clear();
  genVertices_.clear();
  
  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  const edm::ESTransientHandle<MagneticField> magfield = iSetup.getTransientHandle(magneticFieldToken_);
  const CaloGeometry ecalGeometry = iSetup.getData(caloGeometryToken_);

  vector<reco::TransientTrack> ttracks;
  for(const auto &track : *muonEnhancedTracksHandle_)
    ttracks.emplace_back(ttBuilder->build(track));
  
  ttBuilder_ = TTBuilderWrapper(ttBuilder);
  
  NewVertexMerger diLeptonMerger(primaryVertex_, 0.75), hadronicMerger(primaryVertex_, 0.);
  if(cfFlag("hasGenInfo")) {
    GenVertices allSignalSVs(*genHandle_);
    DeltaRGenMatchHungarian<reco::TransientTrack> assigner(ttracks, allSignalSVs.getAllGenParticles());
    
    genVertices_ = GenVertices(assigner.GetPairedObjects().ConvertFromTTracks(), 0.02);
    allSignalSVs += genVertices_;
    genVertices_ = allSignalSVs;    
    diLeptonMerger.setGenVertices(genVertices_);
    hadronicMerger.setGenVertices(genVertices_);

    for(const auto &genVertex : genVertices_) {
      DeltaRGenMatchHungarian<reco::Track> chargedParticleAssigner(*muonEnhancedTracksHandle_, genVertex.getStableChargedDaughters(*genHandle_));
      chargedMatches_[genVertex] = chargedParticleAssigner.GetPairedObjects();
      
      if(!genVertex.hasTracks()) continue;

      for(const auto &pair : genVertex.genMatches())
        signalTracks_.emplace_back(pair.GetObjectA());
    }
  }

  std::vector<reco::TrackRef> trackRefs(ConvertTracksToRefs(muonEnhancedTracksHandle_));
  TrackVertexSetCollection diLeptonCandidates(diLeptonMerger.trackVertexSets(trackRefs, ttBuilder).fullSelection(primaryVertex_));
  TrackVertexSetCollection hadronicCandidates(TrackVertexSetHelper::hadronicVertexCandidates(hadronicMerger.trackVertexSets(trackRefs, ttBuilder), primaryVertex_));

  generalVertices_ = diLeptonCandidates;
  generalVertices_.insert(hadronicCandidates.begin(), hadronicCandidates.end());

  MatchTracksToSC<reco::Track> assigner(iEvent, iSetup, magfield, ecalGeometry, trackAssocParameters_, generalVertices_.tracks(), *mergedSCsHandle_);
  trackSCPairs_ = assigner.GetMatchedTrackSCPairs();
  
  for(const auto &pair : trackSCPairs_) {
    if(pair.GetDeltaR() < 0.04)
      electronTracks_.emplace_back(pair.GetTrack());
  }

}//<<>>void KUCMSDisplacedVertex::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSDisplacedVertex::ProcessEvent( ItemManager<float>& geVar ){

    Branches.clearBranches();

    map<string, TrackVertexSetCollection> vertexMap = {{"passSignificanceCut", TrackVertexSetHelper::applySignificanceCut(generalVertices_, primaryVertex_)},
						       {"passTrackCuts", TrackVertexSetHelper::applyPostDisambiguationTrackCuts(generalVertices_, primaryVertex_)}};

    Branches.fillBranch("Vertex_nTotal", unsigned(generalVertices_.size()) );
    int vtxIndex(0);
    for(const auto &tvertex : generalVertices_) {
      reco::Vertex vertex(tvertex);
      LorentzVec vertex4Vec(VertexHelper::GetVertex4Vector(vertex));
      LorentzVec vertexWeighted4Vec(VertexHelper::GetVertexWeighted4Vector(vertex));
      const bool passLooseMuonID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, *muonTracksHandle_) == 2);
      const bool passLooseElectronID(vertex.tracksSize() == 2 && VertexHelper::CountInstances(vertex, electronTracks_) == 2);
      const bool passLooseID(vertex.tracksSize() < 4? (passLooseMuonID || passLooseElectronID) : true);

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
	Branches.fillBranch("Vertex_matchRatio", float(nearestGenVertexIndex>=0? matchRatio(tvertex, genVertices_[nearestGenVertexIndex]) : -1));
      }
      
      for(const auto &trackRef : vertex.tracks()) {
	const reco::Track track(*trackRef);
	const reco::TrackRef ref(TrackHelper::GetTrackRef(*trackRef, muonEnhancedTracksHandle_));
	
	double deltaR(-1.);
	reco::SuperCluster sc;
	const bool isSCMatched(getSCMatch(track, sc, deltaR));
	Branches.fillBranch("VertexTrack_vertexIndex", unsigned(vtxIndex) );
	Branches.fillBranch("VertexTrack_trackIndex", unsigned(TrackHelper::FindTrackIndex(track, *muonEnhancedTracksHandle_)) );
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
	  Branches.fillBranch("VertexTrack_isSignalElectron", isSignal? genVertices_.getGenVertexFromTrack(track).isGenElectron() : false);
	  Branches.fillBranch("VertexTrack_isSignalMuon", isSignal? genVertices_.getGenVertexFromTrack(track).isGenMuon() : false);
	}
      }
      vtxIndex++;
    }

    if(cfFlag("hasGenInfo")) {
      // Gen Vertex
      GenVertices allSignalSVs(*genHandle_);
      allSignalSVs += genVertices_;

      int nSigElectrons(0), nSigMuons(0), nSigHadrons(0);
      Branches.fillBranch("GenVertex_nTotal", unsigned(allSignalSVs.size()));
      for(const auto &genVertex : allSignalSVs) {
	
	TrackVertexSet tvertex(ConvertTracksToRefs(genVertex.tracks(), muonEnhancedTracksHandle_), ttBuilder_.getBuilder());
	TrackVertexSetCollection vertexSetCollection;
	if(tvertex.isValid()) {
	  vertexSetCollection.add(tvertex);
	  vertexSetCollection = vertexSetCollection.fullSelection(primaryVertex_);
	}

	Branches.fillBranch("GenVertex_nTracks", unsigned(genVertex.tracks().size()));
	Branches.fillBranch("GenVertex_mass", float(genVertex.mass()));
	Branches.fillBranch("GenVertex_x", float(genVertex.x()));
	Branches.fillBranch("GenVertex_y", float(genVertex.y()));
	Branches.fillBranch("GenVertex_z", float(genVertex.z()));
	Branches.fillBranch("GenVertex_p", float(genVertex.p()));
	Branches.fillBranch("GenVertex_px", float(genVertex.px()));
	Branches.fillBranch("GenVertex_py", float(genVertex.py()));
	Branches.fillBranch("GenVertex_pz", float(genVertex.pz()));
	Branches.fillBranch("GenVertex_pt", float(genVertex.pt()));
	Branches.fillBranch("GenVertex_eta", float(genVertex.eta()));
	Branches.fillBranch("GenVertex_phi", float(genVertex.phi()));
	Branches.fillBranch("GenVertex_dxy", float(genVertex.dxy()));
	Branches.fillBranch("GenVertex_isElectron", bool(genVertex.isGenElectron()));
	Branches.fillBranch("GenVertex_isMuon", bool(genVertex.isGenMuon()));
	Branches.fillBranch("GenVertex_isHadronic", bool(genVertex.isGenHadronic()));
	Branches.fillBranch("GenVertex_passSelection", bool(genVertex.hasTracks()));
	Branches.fillBranch("GenVertex_passSelectionAndCuts", bool(genVertex.hasTracks() && !vertexSetCollection.empty()));
	
	if(genVertex.isGenElectron()) nSigElectrons++;
	if(genVertex.isGenMuon()) nSigMuons++;
	if(genVertex.isGenHadronic()) nSigHadrons++;
      }
      Branches.fillBranch("GenVertex_nElectron", unsigned(nSigElectrons));
      Branches.fillBranch("GenVertex_nMuon", unsigned(nSigMuons));
      Branches.fillBranch("GenVertex_nHadronic", unsigned(nSigHadrons));
    }

}//<<>>void KUCMSDisplacedVertex::ProcessEvent()


std::vector<reco::TrackRef> KUCMSDisplacedVertex::ConvertTracksToRefs(const edm::Handle<reco::TrackCollection> &trackHandle) const {
  std::vector<reco::TrackRef> trackRefs;
  for(const auto &track : *trackHandle)
    trackRefs.emplace_back(TrackHelper::GetTrackRef(track, trackHandle));

  return trackRefs;
}

std::vector<reco::TrackRef> KUCMSDisplacedVertex::ConvertTracksToRefs(const reco::TrackCollection &tracks, const edm::Handle<reco::TrackCollection> &trackHandle) const {
  std::vector<reco::TrackRef> trackRefs;
  for(const auto &track : tracks)
    trackRefs.emplace_back(TrackHelper::GetTrackRef(track, trackHandle));

  return trackRefs;
}

bool KUCMSDisplacedVertex::IsGold(const reco::Vertex &vertex) const {

  bool isGold(false);
  for(const auto &genVertex : genVertices_) {
    if(isGold) break;
    isGold = genVertex.isGold(vertex);
  }

  return isGold;
}

bool KUCMSDisplacedVertex::IsSilver(const reco::Vertex &vertex) const {

  bool isSilver(false);
  for(const auto &genVertex : genVertices_) {
    if(isSilver) break;
    isSilver = genVertex.isSilver(vertex);
  }

  return isSilver;
}

bool KUCMSDisplacedVertex::IsBronze(const reco::Vertex &vertex) const {

  bool isBronze(false);
  for(const auto &genVertex : genVertices_) {
    if(isBronze) break;
    isBronze = genVertex.isBronze(vertex);
  }

  return isBronze;
}

int KUCMSDisplacedVertex::FindGenVertexIndex(const reco::Vertex &vertex) const {

  int index(0);
  for(const auto &genVertex : genVertices_) {
    if(genVertex.hasTracks() && (genVertex.isBronze(vertex) || genVertex.isSilver(vertex) || genVertex.isGold(vertex)))
      return index;
    index++;
  }

  return -1;
}

int KUCMSDisplacedVertex::FindNearestGenVertexIndex(const reco::Vertex &vertex, double &distance) const {
  int index(0), minIndex(-1);
  distance = std::numeric_limits<double>::max();
  for(const auto &genVertex : genVertices_) {
    if(!genVertex.isGenHadronic()) continue;
    
    const double tempMin(genVertex.distance3D(vertex));
    if(tempMin < distance) {
      distance = tempMin;
      minIndex = index;
    }
    index++;
  }

  if(minIndex < 0)
    distance = -1.;
  
  return minIndex;
}

bool KUCMSDisplacedVertex::getSCMatch(const reco::Track &track, reco::SuperCluster &sc, double &deltaR) const {

  bool isMatched(false);
  for(const auto &pair : trackSCPairs_) {

    if(TrackHelper::SameTrack(track, pair.GetTrack())) {
      isMatched = true;
      sc = pair.GetSuperCluster();
      deltaR = pair.GetDeltaR();
      break;
    }
  }

  return isMatched;
}

double KUCMSDisplacedVertex::matchRatio(const TrackVertexSet &vertex, const GenVertex &genVertex, const double threshold) const {

  int count(0);
  std::vector<reco::TrackRef> refVec;
  std::set<reco::TrackRef> refSet;
  for(const auto &pair : chargedMatches_.at(genVertex)) {
    reco::TrackRef ref(TrackHelper::GetTrackRef(pair.GetObjectA(), muonEnhancedTracksHandle_));
    refVec.emplace_back(ref);
    refSet.insert(ref);

    const double genPt(pair.GetObjectB().pt());
    const double trackPt(pair.GetObjectA().pt());
    const double relPtDiff((genPt-trackPt)/genPt);
    const double deltaR(pair.GetDeltaR());

    if(vertex.contains(ref) && ((deltaR < 0.01 && fabs(relPtDiff) < (-10*deltaR + 0.15)) ||
                                (deltaR > 0.01 && deltaR < threshold && fabs(relPtDiff) < (-deltaR/2 + 0.055)))) {
      count++;
    }
  }
  return double(count)/vertex.size();
}

#endif
