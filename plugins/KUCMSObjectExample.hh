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

// Add includes for specfic collections
//#include "DataFormats/Example/interface/Example.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSExampleObjectHeader
#define KUCMSExampleObjectHeader

//#define ExampleEBUG true
#define ExampleDEBUG false

using namespace edm; 

class KUCMSExampleObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSExampleObject( const edm::ParameterSet& iConfig );
    ~KUCMSExampleObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadExampleTokens( edm::EDGetTokenT<example> exampleToken_ ); 
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar );

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    private:

    std::vector<examnple> fexamples;

    edm::EDGetTokenT<example> exampleToken_;
    edm::Handle<example> example_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSExample : public KUCMSObjectBase

KUCMSExampleObject::KUCMSExampleObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );

}//<<>>KUCMSExample::KUCMSExample( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSExampleObject::LoadExampleTokens( edm::EDGetTokenT<ExampleCol> ExampleEBToken, edm::EDGetTokenT<ExampleCol> ExampleEEToken ){

    ExampleEBToken_ = ExampleEBToken;

}//<<>>void LoadExampleTokens( edm::EDGetTokenT<ExampleCol>* ExampleEBToken_, edm::EDGetTokenT<ExampleCol>* ExampleEEToken_ )

void KUCMSExampleObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("rhx","Example_rhx",VFLOAT);
    Branches.makeBranch("rhy","Example_rhy",VFLOAT);
    Branches.makeBranch("rhz","Example_rhz",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSExample::InitObject( TTree* fOutTree )

void KUCMSExampleObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS
    iEvent.getByToken( ExampleEBToken_, ExampleEB_);

    if( ExampleDEBUG ) std::cout << "Collecting Examples" << std::endl;
    fexamples.clear();
    for( const auto &Example : *ExampleEB_ ){ if(Example.energy() > cfPrm("minRHEi")) fexamples.push_back(Example);}

}//<<>>void KUCMSExample::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSExampleObject::ProcessEvent( ItemManager<float>& geVar ){

    if( ExampleDEBUG ) std::cout << "Processing Examples" << std::endl;

    Branches.clearBranches();

    if( ExampleDEBUG ) std::cout << " - enetering Example loop" << std::endl;

    int nExamples = fexamples.size();
    for ( int it = 0; it < nExamples; it++ ){

	//,..
        const float phi = ExamplePos.phi();
        const auto rhX = ExamplePos.x();
        const auto rhY = ExamplePos.y();
	//...
        Branches.fillBranch("rhy",rhY);
        Branches.fillBranch("rhz",rhZ);

    }//<<>>for (const auto Example : *ExampleEB_ ) 

}//<<>>void KUCMSExample::ProcessEvent()

void KUCMSExampleObject::EndJobs(){}

#endif
