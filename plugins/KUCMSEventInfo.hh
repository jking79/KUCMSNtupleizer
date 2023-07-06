// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS EventInfo Object
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
//#include "DataFormats/EventInfo/interface/EventInfo.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSEventInfoObjectHeader
#define KUCMSEventInfoObjectHeader

//#define EventInfoEBUG true
#define eiDEBUG false

using namespace edm; 

class KUCMSEventInfoObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSEventInfoObject( const edm::ParameterSet& iConfig );
    ~KUCMSEventInfoObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadVertexTokens( edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_ ){ verticesToken = verticesToken_; }; 
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 


    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar );
    void PostProcessEvent( ItemManager<float>& geVar ){};

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    private:

    ItemManager<uInt> eventVar;

    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken;
    edm::Handle<std::vector<reco::Vertex>> vertices_;   

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSEventInfo : public KUCMSObjectBase

KUCMSEventInfoObject::KUCMSEventInfoObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );

}//<<>>KUCMSEventInfo::KUCMSEventInfo( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSEventInfoObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("erun","Evt_run",UINT,"Run number of event");
    Branches.makeBranch("elumi","Evt_luminosityBlock",UINT,"Lumi section of event");
    Branches.makeBranch("event","Evt_event",UINT);
    Branches.makeBranch("vnVtx","PV_npvs",UINT);
    Branches.makeBranch("vtxX","PV_x",FLOAT);
    Branches.makeBranch("vtxY","PV_y",FLOAT);
    Branches.makeBranch("vtxZ","PV_z",FLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSEventInfo::InitObject( TTree* fOutTree )

void KUCMSEventInfoObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS
    iEvent.getByToken( verticesToken, vertices_ );

    //if( eiDEBUG ) std::cout << "Collecting EventInfos" << std::endl;

    eventVar.set("run",iEvent.id().run());
    eventVar.set("lumi", iEvent.luminosityBlock());
    eventVar.set("event", iEvent.id().event());

    if( eiDEBUG ) std::cout << "Processing event: " << eventVar("event");
    if( eiDEBUG ) std::cout << " in run: " << eventVar("run") << " and lumiblock: " << eventVar("lumi") << std::endl;

    const auto & primevtx = vertices_->front();
    const float vtxX = primevtx.position().x();
    const float vtxY = primevtx.position().y();
    const float vtxZ = primevtx.position().z();

    // set global event varibles for vertix
    geVar.set("vtxX",vtxX);
    geVar.set("vtxY",vtxY);
    geVar.set("vtxZ",vtxZ);

}//<<>>void KUCMSEventInfo::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSEventInfoObject::ProcessEvent( ItemManager<float>& geVar ){

    if( eiDEBUG ) std::cout << "Processing EventInfo" << std::endl;

    Branches.clearBranches();

    Branches.fillBranch("erun",eventVar("run"));
    Branches.fillBranch("elumi",eventVar("lumi"));
    Branches.fillBranch("event",eventVar("event"));

    // do not need to pass in collectioon handel -> class varible   --  Fix This
    //edm::Handle<std::vector<reco::Vertex>> vertices_
    uInt nvtx = vertices_->size();

    Branches.fillBranch("vnVtx",nvtx);
    Branches.fillBranch("vtxX",geVar("vtxX"));
    Branches.fillBranch("vtxY",geVar("vtxY"));
    Branches.fillBranch("vtxZ",geVar("vtxZ"));

}//<<>>void KUCMSEventInfo::ProcessEvent()

void KUCMSEventInfoObject::EndJobs(){}

#endif
