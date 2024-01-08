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
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"
#include "KUCMSGenObjects.hh"

#define DEBUG false

class KUCMSTimedSVs : public KUCMSObjectBase {
  
public:
  
  KUCMSTimedSVs(const edm::ParameterSet& iConfig);
  ~KUCMSTimedSVs(){}

  void LoadSVTokens(edm::EDGetTokenT<edm::View<reco::Vertex> > svToken){ svToken_ = svToken; } 
  void LoadElectronTokens(edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronToken ){ electronToken_ = electronToken; }
  void LoadGenObject(KUCMSGenObject* genObjs){ genObjs_ = genObjs; };

  // sets up branches, do preloop jobs 
  void InitObject(TTree* fOutTree); 

  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar);

  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent(ItemManager<float>& geVar);
  void PostProcessEvent(ItemManager<float>& geVar){}

  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(){}

private:
  
  typedef std::map<std::pair<float, float>, reco::GsfElectron> TrackToElectronMap;

  reco::VertexCollection timedVertices_;
  reco::GsfElectronCollection electrons_;
  TrackToElectronMap trackToElectronMap_;

  // SV collection from timed SVs producer
  edm::EDGetTokenT<edm::View<reco::Vertex> > svToken_;
  edm::Handle<edm::View<reco::Vertex> > svHandle_;

  // Electron collection (currently using low pt electrons)
  edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronToken_;
  edm::Handle<edm::View<reco::GsfElectron> > electronHandle_;

  // Gen Particles collection
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticlesToken_;
  edm::Handle<std::vector<reco::GenParticle>> genParticles_;

  KUCMSGenObject* genObjs_;

  void SortByPt(reco::TrackCollection &tracks) const;
  void SortByChi2(reco::VertexCollection &vertices) const;
};

KUCMSTimedSVs::KUCMSTimedSVs(const edm::ParameterSet& iConfig) {
  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );

}

void KUCMSTimedSVs::InitObject(TTree* fOutTree) {  

  // Vertex Branches
  Branches.makeBranch("nVertex","nVertex", INT);
  Branches.makeBranch("Vertex_x","Vertex_x", VFLOAT);
  Branches.makeBranch("Vertex_y","Vertex_y", VFLOAT);
  Branches.makeBranch("Vertex_z","Vertex_z", VFLOAT);
  Branches.makeBranch("Vertex_chi2","Vertex_chi2", VFLOAT);
  Branches.makeBranch("Vertex_normalizedChi2","Vertex_normalizedChi2", VFLOAT);
  Branches.makeBranch("Vertex_ndof","Vertex_ndof", VFLOAT);

  // Track Branches
  Branches.makeBranch("nTracks","nTracks",VINT);
  Branches.makeBranch("LinearDist","LinearDist", VFLOAT);
  Branches.makeBranch("PathLength","PathLength",VFLOAT);
  Branches.makeBranch("Track_pt", "Track_pt", VFLOAT);
  Branches.makeBranch("Track_eta", "Track_eta", VFLOAT);
  Branches.makeBranch("Track_phi", "Track_phi", VFLOAT);
  Branches.makeBranch("Track_x", "Track_x", VFLOAT);
  Branches.makeBranch("Track_y", "Track_y", VFLOAT);
  Branches.makeBranch("Track_z", "Track_z", VFLOAT);
  Branches.makeBranch("Track0_pt", "Track0_pt", VFLOAT);
  Branches.makeBranch("Track1_pt", "Track1_pt", VFLOAT);
  Branches.makeBranch("Track0_pathLength", "Track0_pathLength", VFLOAT);
  Branches.makeBranch("Track1_pathLength", "Track1_pathLength", VFLOAT);

  // Matched Electron Branches
  Branches.makeBranch("VertexID","VertexID", VINT);
  Branches.makeBranch("PositionAtECAL_x","PositionAtECAL_x", VFLOAT);
  Branches.makeBranch("PositionAtECAL_y","PositionAtECAL_y", VFLOAT);
  Branches.makeBranch("PositionAtECAL_z","PositionAtECAL_z", VFLOAT);
  Branches.makeBranch("Electron0_pt", "Electron0_pt", VFLOAT);
  Branches.makeBranch("Electron1_pt", "Electron1_pt", VFLOAT);
  Branches.makeBranch("LowPtElectron_mass", "LowPtElectron_mass", VFLOAT);
  Branches.makeBranch("LowPtElectron_energy", "LowPtElectron_energy", VFLOAT);
  Branches.makeBranch("LowPtElectron_et", "LowPtElectron_et", VFLOAT);
  Branches.makeBranch("LowPtElectron_pt", "LowPtElectron_pt", VFLOAT);
  Branches.makeBranch("LowPtElectron_eta", "LowPtElectron_eta", VFLOAT);
  Branches.makeBranch("LowPtElectron_theta", "LowPtElectron_theta", VFLOAT);
  Branches.makeBranch("LowPtElectron_phi", "LowPtElectron_phi", VFLOAT);
  Branches.makeBranch("LowPtElectron_charge", "LowPtElectron_charge", VFLOAT);
  Branches.makeBranch("isEB","isEB", VBOOL);
  Branches.makeBranch("isEE","isEE", VBOOL);
  Branches.makeBranch("SC_energy","SC_energy", VFLOAT); 
  
  // Gen matching Branches
  Branches.makeBranch("PdgID", "PdgID", VINT);
  Branches.makeBranch("SPdgID", "SPdgID", VINT);
  Branches.makeBranch("LowPtElectronGenIdx", "LowPtElectronGenIdx", VINT);
  Branches.makeBranch("LowPtElectronGenDr", "LowPtElectronGenDr", VFLOAT);
  Branches.makeBranch("LowPtElectronGenDp", "LowPtElectronGenDp", VFLOAT);
  Branches.makeBranch("LowPtElectronGenSIdx", "LowPtElectronGenSIdx", VINT);
  Branches.makeBranch("LowPtElectronGenSDr", "LowPtElectronGenSDr", VFLOAT);
  Branches.makeBranch("LowPtElectronGenSDr", "LowPtElectronGenSDr", VFLOAT);
  Branches.makeBranch("DiffVertex_x","DiffVertex_x", VFLOAT, "Difference between gen and reco vertex in x coordinate");
  Branches.makeBranch("DiffVertex_y","DiffVertex_y", VFLOAT, "Difference between gen and reco vertex in y coordinate");
  Branches.makeBranch("DiffVertex_z","DiffVertex_z", VFLOAT, "Difference between gen and reco vertex in z coordinate");
  Branches.makeBranch("GenVertex_x","GenVertex_x", VFLOAT, "Gen particle vertex x coordinate");
  Branches.makeBranch("GenVertex_y","GenVertex_y", VFLOAT, "Gen particle vertex y coordinate");
  Branches.makeBranch("GenVertex_z","GenVertex_z", VFLOAT, "Gen particle vertex z coordinate");
  Branches.makeBranch("GenVertexPV_x","GenVertexPV_x", FLOAT, "Gen particle primary vertex x coordinate");
  Branches.makeBranch("GenVertexPV_y","GenVertexPV_y", FLOAT, "Gen particle primary vertex y coordinate");
  Branches.makeBranch("GenVertexPV_z","GenVertexPV_z", FLOAT, "Gen particle primary vertex z coordinate");
  Branches.makeBranch("isW","isW", VBOOL);
  Branches.makeBranch("isZ","isZ", VBOOL);
  Branches.makeBranch("isTau","isTau", VBOOL);
  Branches.makeBranch("isLight","isLight", VBOOL);
  Branches.makeBranch("isHeavy","isHeavy", VBOOL);
  Branches.makeBranch("isConversion","isConversion", VBOOL);
  Branches.makeBranch("isUnmatched","isUnmatched", VBOOL);

  //Branches.makeBranch("","", VINT);
  //Branches.makeBranch("","", VFLOAT);
  //Branches.makeBranch("","", VBOOL);
  Branches.attachBranches(fOutTree);
}

void KUCMSTimedSVs::LoadEvent(const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar) {

  iEvent.getByToken(svToken_, svHandle_);
  iEvent.getByToken(electronToken_, electronHandle_);

  if(DEBUG) std::cout << "Collecting Timed SVs" << std::endl;
  timedVertices_.clear();

  // Fill vertex container for the event 
  for(const auto &vertex : *svHandle_) {
    timedVertices_.push_back(vertex);
  }

  for(const auto &electron : *electronHandle_) {
    electrons_.push_back(electron);
    const reco::GsfTrack *track = electron.gsfTrack().get();
    const double pt = track->pt();
    const double chi2 = track->chi2();

    trackToElectronMap_.emplace(std::make_pair(chi2, pt), electron);
  }
}

void KUCMSTimedSVs::ProcessEvent(ItemManager<float>& geVar) {

  if(DEBUG) std::cout << "Processing timed secondary vertices" << std::endl;

  Branches.clearBranches();
  
  SortByChi2(timedVertices_);
  Branches.fillBranch("nVertex", int(timedVertices_.size()) );

  if(cfFlag("hasGenInfo")) {
    Branches.fillBranch("GenVertexPV_x", float(genObjs_->GenVertex().x()) );
    Branches.fillBranch("GenVertexPV_y", float(genObjs_->GenVertex().y()) );
    Branches.fillBranch("GenVertexPV_z", float(genObjs_->GenVertex().z()) );
  }

  // Loop over vertex collection
  int vertexIdx = 0;
  for(reco::Vertex const &vertex : timedVertices_) {

    Branches.fillBranch("Vertex_x", float(vertex.x()) ); 
    Branches.fillBranch("Vertex_y", float(vertex.y()) );
    Branches.fillBranch("Vertex_z", float(vertex.z()) );
    Branches.fillBranch("Vertex_chi2", float(vertex.chi2()) );
    Branches.fillBranch("Vertex_normalizedChi2", float(vertex.normalizedChi2()) );
    Branches.fillBranch("Vertex_ndof", float(vertex.ndof()) );

    const int nTracks = vertex.refittedTracks().size();
    Branches.fillBranch("nTracks", nTracks);

    // Get tracks and sort them by pt
    reco::TrackCollection tracks = vertex.refittedTracks();
    SortByPt(tracks);

    // Get leading and subleading electron and track information
    const reco::GsfElectron leadingElectron = trackToElectronMap_[std::make_pair(tracks[0].chi2(), tracks[0].pt())];
    const reco::GsfElectron subLeadingElectron = trackToElectronMap_[std::make_pair(tracks[1].chi2(), tracks[1].pt())];

    Branches.fillBranch("Track0_pt", float(tracks[0].pt()) );
    Branches.fillBranch("Track1_pt", float(tracks[1].pt()) );
    Branches.fillBranch("Track0_pathLength", float(tracks[0].t0()) );
    Branches.fillBranch("Track1_pathLength", float(tracks[1].t0()) );
    Branches.fillBranch("Electron0_pt", float(leadingElectron.pt()) );
    Branches.fillBranch("Electron1_pt", float(subLeadingElectron.pt()) );

    // Loop over tracks in vertex
    for(reco::Track const &track : tracks){
      
      const reco::GsfElectron electron = trackToElectronMap_[std::make_pair(track.chi2(), track.pt())];
      Branches.fillBranch("VertexID", vertexIdx);
      Branches.fillBranch("PositionAtECAL_x", float(electron.trackPositionAtCalo().x()) );
      Branches.fillBranch("PositionAtECAL_y", float(electron.trackPositionAtCalo().y()) );
      Branches.fillBranch("PositionAtECAL_z", float(electron.trackPositionAtCalo().z()) ); 
      Branches.fillBranch("PathLength", float(track.t0()) );
      Branches.fillBranch("LinearDist", float(track.beta()) );
      Branches.fillBranch("Track_x", float(track.vx()) );
      Branches.fillBranch("Track_y", float(track.vy()) );
      Branches.fillBranch("Track_z", float(track.vz()) );
      Branches.fillBranch("Track_pt", float(track.pt()) );
      Branches.fillBranch("Track_eta", float(track.eta()) );
      Branches.fillBranch("Track_phi", float(track.phi()) );
      Branches.fillBranch("isEB", bool(electron.isEB()) );
      Branches.fillBranch("isEE", bool(electron.isEE()) );
      Branches.fillBranch("LowPtElectron_mass", float(electron.mass()) );
      Branches.fillBranch("LowPtElectron_energy", float(electron.energy()) );
      Branches.fillBranch("LowPtElectron_et", float(electron.et()) );
      Branches.fillBranch("LowPtElectron_pt", float(electron.pt()) );
      Branches.fillBranch("LowPtElectron_eta", float(electron.eta()) );
      Branches.fillBranch("LowPtElectron_theta", float(electron.theta()) );
      Branches.fillBranch("LowPtElectron_phi", float(electron.phi()) );
      Branches.fillBranch("LowPtElectron_charge", float(electron.charge()) );
      Branches.fillBranch("SC_energy", float(electron.superCluster()->correctedEnergy()) ); 
      //Branches.fillBranch("", );

      if(!cfFlag("hasGenInfo"))
	continue;

      std::map<string,float> genInfo = genObjs_->getGenPartMatch(*(electron.superCluster().get()), electron.pt() );
      std::map<string, bool> momType = genObjs_->MotherID(genInfo["matchedIdx"]);

      Branches.fillBranch("PdgID", int(genInfo["matchedId"]));
      Branches.fillBranch("SPdgID", int(genInfo["matchedSId"]));
      Branches.fillBranch("isW", momType["isW"]);
      Branches.fillBranch("isZ", momType["isZ"]);
      Branches.fillBranch("isTau", momType["isTau"]);
      Branches.fillBranch("isLight", momType["isLight"]);
      Branches.fillBranch("isHeavy", momType["isHeavy"]);
      Branches.fillBranch("isConversion", momType["isConversion"]);
      Branches.fillBranch("isUnmatched", momType["isUnmatched"]);
      Branches.fillBranch("LowPtElectronGenIdx", int(genInfo["matchedIdx"]));
      Branches.fillBranch("LowPtElectronGenDr", genInfo["GenDR"]);
      Branches.fillBranch("LowPtElectronGenDp", genInfo["GenDpt"]);
      Branches.fillBranch("LowPtElectronGenSIdx", int(genInfo["matchedSIdx"]));
      Branches.fillBranch("LowPtElectronGenSDr", genInfo["GenSDR"]);
      Branches.fillBranch("LowPtElectronGenSDr", genInfo["GenSDpt"]);
      Branches.fillBranch("GenVertex_x", genInfo["GenVertex_x"]);
      Branches.fillBranch("GenVertex_y", genInfo["GenVertex_y"]);
      Branches.fillBranch("GenVertex_z", genInfo["GenVertex_z"]);
      Branches.fillBranch("DiffVertex_x", float(genInfo["GenVertex_x"] - vertex.x()) );
      Branches.fillBranch("DiffVertex_y", float(genInfo["GenVertex_y"] - vertex.y()) );
      Branches.fillBranch("DiffVertex_z", float(genInfo["GenVertex_z"] - vertex.z()) );
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

#endif
