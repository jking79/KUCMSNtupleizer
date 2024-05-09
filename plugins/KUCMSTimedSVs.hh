// -*- C++ -*-
//
//
// Original Author:  Andres Abreu
// 
// Kucms Timed Secondary Vertex Object
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

#ifndef KUCMSTimedSVs_hh
#define KUCMSTimedSVs_hh

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

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackTools.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/VertexAssembly.h"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/DeltaRMatchApplications.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#define DEBUG false

class KUCMSTimedSVs : public KUCMSObjectBase {
  
public:
  
  KUCMSTimedSVs(const edm::ParameterSet& iConfig) {}
  virtual ~KUCMSTimedSVs() = default;

  void LoadSVTokens(edm::EDGetTokenT<edm::View<reco::Vertex> > svToken){ svToken_ = svToken; } 
  void LoadElectronToken(edm::EDGetTokenT<reco::ElectronCollection> electronToken ){ electronToken_ = electronToken; }
  void LoadElectronTracksToken(edm::EDGetTokenT<reco::TrackCollection> token ){ trackToken_ = token; }
  void LoadGenParticlesToken(edm::EDGetTokenT<reco::GenParticleCollection> token ){ genToken_ = token; }
  void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }

  // sets up branches, do preloop jobs 
  void InitObject(TTree* fOutTree); 

  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar);

  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent(ItemManager<float>& geVar) {}
  void PostProcessEvent(ItemManager<float>& geVar);

  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(){}

private:
  
  reco::VertexCollection timedVertices_;
  reco::ElectronCollection signalElectrons_;
  reco::TrackCollection tracks_;
 
  // SV collection from timed SVs producer
  edm::EDGetTokenT<edm::View<reco::Vertex> > svToken_;
  edm::Handle<edm::View<reco::Vertex> > svHandle_;

  // DisplacedElectron collection 
  edm::EDGetTokenT<reco::ElectronCollection> electronToken_;
  edm::Handle<reco::ElectronCollection> electronHandle_;

  // Displaced Electron Tracks
  edm::EDGetTokenT<reco::TrackCollection> trackToken_;
  edm::Handle<reco::TrackCollection> trackHandle_;

  //Gen particles
  edm::EDGetTokenT<reco::GenParticleCollection> genToken_;
  edm::Handle<reco::GenParticleCollection> genHandle_;

  // Transient track builder
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

  VertexAssembly vertexBuilder_;

  void SortByPt(reco::TrackCollection &tracks) const;
  void SortByChi2(reco::VertexCollection &vertices) const;
  int FindIndex(const reco::Track &track) const;  

  reco::TrackCollection ExtractTracks(const reco::ElectronCollection &electrons) const;
  
};

void KUCMSTimedSVs::InitObject(TTree* fOutTree) {  

  Branches.makeBranch("TimedSV_nTotal","TimedSV_nTotal", INT);
  Branches.makeBranch("TimedSV_nTracks","TimedSV_nTracks", VINT);
  Branches.makeBranch("TimedSV_vertexIndex","TimedSV_vertexIndex", VUINT);
  Branches.makeBranch("TimedSV_displacedElectronIndex","TimedSV_displacedElectronIndex", VINT);
  Branches.makeBranch("TimedSV_x","TimedSV_x", VFLOAT);
  Branches.makeBranch("TimedSV_y","TimedSV_y", VFLOAT);
  Branches.makeBranch("TimedSV_z","TimedSV_z", VFLOAT);
  Branches.makeBranch("TimedSV_dxy","TimedSV_dxy", VFLOAT);
  Branches.makeBranch("TimedSV_chi2","TimedSV_chi2", VFLOAT);
  Branches.makeBranch("TimedSV_normalizedChi2","TimedSV_normalizedChi2", VFLOAT);
  Branches.makeBranch("TimedSV_ndof","TimedSV_ndof", VFLOAT);
  Branches.makeBranch("TimedSV_pathLength","TimedSV_pathLength", VFLOAT);
  Branches.makeBranch("TimedSV_linearDistance","TimedSV_linearDistance", VFLOAT);
  Branches.makeBranch("TimedSV_trackPt","TimedSV_trackPt", VFLOAT);

  //SignalSV
  Branches.makeBranch("SignalSV_nTotal","SignalSV_nTotal", INT);
  Branches.makeBranch("SignalSV_nTracks","SignalSV_nTracks", VINT);
  Branches.makeBranch("SignalSV_vertexIndex","SignalSV_vertexIndex", VUINT);
  Branches.makeBranch("SignalSV_displacedElectronIndex","SignalSV_displacedElectronIndex", VINT);
  Branches.makeBranch("SignalSV_trackWeight","SignalSV_trackWeight", VFLOAT);
  Branches.makeBranch("SignalSV_x","SignalSV_x", VFLOAT);
  Branches.makeBranch("SignalSV_y","SignalSV_y", VFLOAT);
  Branches.makeBranch("SignalSV_z","SignalSV_z", VFLOAT);
  Branches.makeBranch("SignalSV_dxy","SignalSV_dxy", VFLOAT);
  Branches.makeBranch("SignalSV_chi2","SignalSV_chi2", VFLOAT);
  Branches.makeBranch("SignalSV_normalizedChi2","SignalSV_normalizedChi2", VFLOAT);
  Branches.makeBranch("SignalSV_ndof","SignalSV_ndof", VFLOAT);
  //Branches.makeBranch("","", );
  Branches.attachBranches(fOutTree);
}

void KUCMSTimedSVs::LoadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar) {

  iEvent.getByToken(svToken_, svHandle_);
  iEvent.getByToken(electronToken_, electronHandle_);
  iEvent.getByToken(trackToken_, trackHandle_);
  iEvent.getByToken(genToken_, genHandle_);

  const TransientTrackBuilder* ttBuilder = &iSetup.getData(transientTrackBuilder_);
  vertexBuilder_ = VertexAssembly(ttBuilder, trackHandle_);

  if(DEBUG) std::cout << "Collecting Timed SVs" << std::endl;
  timedVertices_.clear();
  tracks_.clear();

  // Fill vertex container for the event 
  for(const auto &vertex : *svHandle_) 
    timedVertices_.push_back(vertex);

  tracks_ = ExtractTracks(*electronHandle_);
  signalElectrons_ = GetSignalElectrons(*electronHandle_, *genHandle_);
}

void KUCMSTimedSVs::PostProcessEvent(ItemManager<float>& geVar) {

  if(DEBUG) std::cout << "Processing timed secondary vertices" << std::endl;

  Branches.clearBranches();
  
  SortByChi2(timedVertices_);
  Branches.fillBranch("TimedSV_nTotal", int(timedVertices_.size()) );
  
  // Loop over vertex collection
  int vertexIdx = 0;
  for(reco::Vertex const &vertex : timedVertices_) {

    Branches.fillBranch("TimedSV_x", float(vertex.x()) ); 
    Branches.fillBranch("TimedSV_y", float(vertex.y()) );
    Branches.fillBranch("TimedSV_z", float(vertex.z()) );
    Branches.fillBranch("TimedSV_chi2", float(vertex.chi2()) );
    Branches.fillBranch("TimedSV_normalizedChi2", float(vertex.normalizedChi2()) );
    Branches.fillBranch("TimedSV_ndof", float(vertex.ndof()) );
    Branches.fillBranch("TimedSV_dxy", float(sqrt(vertex.x()*vertex.x()+vertex.y()*vertex.y())) );

    const int nTracks = vertex.refittedTracks().size();
    Branches.fillBranch("TimedSV_nTracks", nTracks);
    
    // Get tracks and sort them by pt
    reco::TrackCollection tracks = vertex.refittedTracks();
    SortByPt(tracks);
  
    // Loop over tracks in vertex
    for(reco::Track const &track : tracks) { 
      Branches.fillBranch("TimedSV_displacedElectronIndex", FindIndex(track) );
      Branches.fillBranch("TimedSV_vertexIndex", unsigned(vertexIdx) );
      Branches.fillBranch("TimedSV_pathLength", float(track.t0()) );
      Branches.fillBranch("TimedSV_linearDistance", float(track.beta()) );
      Branches.fillBranch("TimedSV_trackPt", float(track.pt()) );
    }
    vertexIdx++;
  }

  //Signal Vertices
  reco::VertexCollection signalVertices = vertexBuilder_.CreateVertexCollection(signalElectrons_);

  Branches.fillBranch("SignalSV_nTotal", int(signalVertices.size()));

  vertexIdx = 0;
  for(reco::Vertex &vertex : signalVertices) {

    Branches.fillBranch("SignalSV_x", float(vertex.x()) );
    Branches.fillBranch("SignalSV_y", float(vertex.y()) );
    Branches.fillBranch("SignalSV_z", float(vertex.z()) );
    Branches.fillBranch("SignalSV_dxy", float(sqrt(vertex.x()*vertex.x()+vertex.y()*vertex.y())) );
    Branches.fillBranch("SignalSV_chi2", float(vertex.chi2()) );
    Branches.fillBranch("SignalSV_normalizedChi2", float(vertex.normalizedChi2()) );
    Branches.fillBranch("SignalSV_ndof", float(vertex.ndof()) );

    reco::TrackCollection tracks;
    for(const auto &track : vertex.tracks()) tracks.emplace_back(*(track.get()));

    Branches.fillBranch("SignalSV_nTracks", int(tracks.size()) );
    for(reco::Track const &track : tracks) {
      Branches.fillBranch("SignalSV_trackWeight", float(vertex.trackWeight(reco::TrackRef(trackHandle_, FindIndex(track)))));
      Branches.fillBranch("SignalSV_displacedElectronIndex", FindIndex(track) );
      Branches.fillBranch("SignalSV_vertexIndex", unsigned(vertexIdx) );
    }
    vertexIdx++;
  } 
}

// Private methods
void KUCMSTimedSVs::SortByPt(reco::TrackCollection &tracks) const {

  // Lambda function for sorting in descending order based on pt
  auto compareByPt = [](const reco::Track& a, const reco::Track& b) {
    return a.pt() > b.pt();
  };

  // Sorting the vector using the lambda function
  std::sort(tracks.begin(), tracks.end(), compareByPt);
}

void KUCMSTimedSVs::SortByChi2(reco::VertexCollection &vertices) const {

  // Lambda function for sorting in descending order based on pt
  auto compareByChi2 = [](const reco::Vertex& a, const reco::Vertex& b) {
    return a.chi2() < b.chi2();
  };

  // Sorting the vector using the lambda function
  std::sort(vertices.begin(), vertices.end(), compareByChi2);
}

int KUCMSTimedSVs::FindIndex(const reco::Track &track) const {

  int index(-1);  
  for(int t = 0; t < int(tracks_.size()); t++) {
    if(track.pt() == tracks_[t].pt() && track.chi2() == tracks_[t].chi2())
      index = t;
  } 
  
  return index;
}

reco::TrackCollection KUCMSTimedSVs::ExtractTracks(const reco::ElectronCollection &electrons) const {

  reco::TrackCollection tracks;

  for(const auto &electron : electrons) {
    tracks.emplace_back(*(electron.track().get()) );
  }
  return tracks;
}

#endif
