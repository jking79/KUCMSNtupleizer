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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"

// Add includes for specfic collections
//#include "DataFormats/Example/interface/Example.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSSecondaryVertex_hh
#define KUCMSSecondaryVertex_hh

//#define DEBUG true
#define DEBUG false

using namespace edm; 

class KUCMSSecondaryVertex : public KUCMSObjectBase {

public:

  // use base class constructor
  KUCMSSecondaryVertex(const edm::ParameterSet& iConfig);
  ~KUCMSSecondaryVertex(){};
  
  // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
  // load tokens for eventt based collections
  void LoadSVTokens(edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svToken) {svToken_ = svToken; }; 
  // sets up branches, do preloop jobs 
  void InitObject( TTree* fOutTree ); 
  // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
  // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
  
  // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
  // get collections, do initial processing
  void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
  // do cross talk jobs with other objects, do event processing, and load branches
  void ProcessEvent( ItemManager<float>& geVar );
  void PostProcessEvent( ItemManager<float>& geVar );
  
  // if there are any final tasks be to done after the event loop via objectManager
  void EndJobs(); // do any jobs that need to be done after main event loop
  
  // New functions specific to this collection
  // void answerCrossTalk(); // define functions that will be called in another object - this is an example
  // ect ...
  
private:
  
  std::vector<reco::VertexCompositePtrCandidate> fsvs;
  
  const edm::InputTag svTag_;
  edm::EDGetTokenT<edm::View<reco::VertexCompositePtrCandidate>> svToken_;
  edm::Handle<edm::View<reco::VertexCompositePtrCandidate>> svs_;
  
  // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
  // exampleObject* otherObjectPtr;
  
};//<<>>class KUCMSExample : public KUCMSObjectBase

KUCMSSecondaryVertex::KUCMSSecondaryVertex(const edm::ParameterSet& iConfig) :
  svTag_(iConfig.getParameter<edm::InputTag>("svSrc")) {   

}//<<>>KUCMSExample::KUCMSExample( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSSecondaryVertex::InitObject( TTree* fOutTree ){

    Branches.makeBranch("nSV","nSV",VINT);
    Branches.makeBranch("svCharge","SV_charge",VINT);
    Branches.makeBranch("svPt","SV_pt",VFLOAT);
    Branches.makeBranch("svEta","SV_eta",VFLOAT);
    Branches.makeBranch("svPhi","SV_phi",VFLOAT);
    Branches.makeBranch("svMass","SV_mass",VFLOAT);
    Branches.makeBranch("svX","SV_x",VFLOAT);
    Branches.makeBranch("svY","SV_y",VFLOAT);
    Branches.makeBranch("svZ","SV_z",VFLOAT);
    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSExample::InitObject( TTree* fOutTree )

void KUCMSSecondaryVertex::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

  iEvent.getByToken(svToken_, svs_);

  if(DEBUG) std::cout << "Collecting SVs" << std::endl;
  fsvs.clear();

  for(const auto &sv : *svs_) {
    fsvs.push_back(sv);
  }


}//<<>>void KUCMSExample::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSSecondaryVertex::ProcessEvent( ItemManager<float>& geVar ){

    if(DEBUG) std::cout << "Processing SVs" << std::endl;

    Branches.clearBranches();

    if(DEBUG) std::cout << " - entering SV loop" << std::endl;

    VertexDistance3D vdist;
    VertexDistanceXY vdistXY;

    const int nSV = fsvs.size();
    for (auto const &sv : fsvs) {

      int charge = 0;
      const float pt = sv.pt();
      const float eta = sv.eta();
      const float phi = sv.phi();
      const float mass = sv.mass();
      const float x = sv.vx();
      const float y = sv.vy();
      const float z = sv.vz();

      // Add total charge in SV
      for (unsigned int id = 0; id < sv.numberOfDaughters(); ++id) {
	const reco::Candidate* daughter = sv.daughter(id);
	charge += daughter->charge();
      }

      Branches.fillBranch("nSV", nSV);
      Branches.fillBranch("svCharge", charge);
      Branches.fillBranch("svPt", pt);
      Branches.fillBranch("svEta", eta);
      Branches.fillBranch("svPhi", phi);
      Branches.fillBranch("svMass", mass);
      Branches.fillBranch("svX", x);
      Branches.fillBranch("svY", y);
      Branches.fillBranch("svZ", z);
    }
    
}//<<>>void KUCMSExample::ProcessEvent()

void KUCMSSecondaryVertex::PostProcessEvent( ItemManager<float>& geVar ){

    // process anything that must be done after all crosstalk with other objects

}//<<>>void KUCMSExampleObject::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSSecondaryVertex::EndJobs(){}

#endif
