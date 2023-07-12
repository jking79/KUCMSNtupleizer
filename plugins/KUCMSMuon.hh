// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Example Object
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

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSMuon_hh
#define KUCMSMuon_hh

//#define ExampleEBUG true
#define DEBUG false

using namespace edm; 

class KUCMSMuon : public KUCMSObjectBase {

public:
  
  // use base class constructor
  KUCMSMuon(const edm::ParameterSet& iConfig);
  ~KUCMSMuon(){};
  
  // load tokens for event based collections
  void LoadMuonTokens(edm::EDGetTokenT<std::vector<reco::Muon>> muonsToken ){ muonsToken_ = muonsToken; }; 

  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree ); 

  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );

  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar );
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(); // do any jobs that need to be done after main event loop
  
private:

  std::vector<reco::Muon> fmuons;  

  const edm::InputTag muonsTag;
  edm::EDGetTokenT<std::vector<reco::Muon> > muonsToken_;
  edm::Handle<std::vector<reco::Muon> > muons_;
  
};

KUCMSMuon::KUCMSMuon(const edm::ParameterSet& iConfig) :
  muonsTag(iConfig.getParameter<edm::InputTag>("muons")) {   
  cfPrm.set( "minMuonE", iConfig.existsAs<double>("minMuonE") ? iConfig.getParameter<double>("minMuonE") : 1.0 );    
}

void KUCMSMuon::InitObject( TTree* fOutTree ){
  
  Branches.makeBranch("muonCharge","Muon_charge",VFLOAT);
  Branches.makeBranch("muonPt","Muon_pt",VFLOAT);
  Branches.makeBranch("muonPhi","Muon_phi",VFLOAT);
  Branches.makeBranch("muonEta","Muon_eta",VFLOAT);  
  Branches.makeBranch("muonMass","Muon_mass",VFLOAT);
  Branches.attachBranches(fOutTree);
}

void KUCMSMuon::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){
  
  iEvent.getByToken(muonsToken_, muons_);
  
  if(DEBUG) std::cout << "Collecting Muons" << std::endl;
  fmuons.clear();
  
  for(const auto &muon : *muons_) { 
    if(muon.energy() > cfPrm("minMuonE")) 
      fmuons.push_back(muon); 
  }
  
}

void KUCMSMuon::ProcessEvent( ItemManager<float>& geVar ){
  
  if(DEBUG) std::cout << "Processing Muons" << std::endl;
  Branches.clearBranches();

  if(DEBUG) std::cout << " - entering Muon loop" << std::endl;

  for(auto const &muon : fmuons){
    const int charge = muon.charge();
    const float pt = muon.pt();
    const float phi = muon.phi();
    const float eta = muon.eta();
    const float mass = muon.mass();

    Branches.fillBranch("muonCharge", charge);
    Branches.fillBranch("muonPt", pt);
    Branches.fillBranch("muonPhi", phi);
    Branches.fillBranch("muonEta", eta);
    Branches.fillBranch("muonMass", mass);
    
    //Branches.attachBranches(outTree);

  }  
}

void KUCMSMuon::PostProcessEvent( ItemManager<float>& geVar ){
  
  // process anything that must be done after all crosstalk with other objects
  
}//<<>>void KUCMSExampleObject::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSMuon::EndJobs(){}

#endif
