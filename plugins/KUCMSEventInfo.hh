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
#define EventInfoDEBUG false

using namespace edm; 

class KUCMSEventInfoObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSEventInfoObject( const edm::ParameterSet& iConfig );
    ~KUCMSEventInfoObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadVertexTokens( edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_ ){ verticesToken = verticesToken_; }; 
    void LoadTriggerTokens( edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_ , edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_ ){ triggerResultsToken = triggerResultsToken_; triggerEventToken = triggerEventToken_; }; 
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

    ItemManager<uInt> eventVar;

    edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
    edm::Handle<edm::TriggerResults> triggerResults;
    edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken;
    edm::Handle<trigger::TriggerEvent> triggerEvent;
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken;
    edm::Handle<std::vector<reco::Vertex>> vertices_;   

	std::map<std::string,bool> trigFlags;
	std::vector<std::string> triggerList;
	std::vector<std::string> fullTriggerList;

    std::map<std::string,bool> flags;
	void fillFlagBranch( std::string flag ){ if( flags.find(flag) != flags.end() ) Branches.fillBranch( flag, flags[flag] ); }
	// for fillFlagBranch to work the branch reffrence name and the cms process path name for the flag must be the same !!!!!!
	// and the flags map must be use the branch/path name as its key !!!!!!!

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

};//<<>>class KUCMSEventInfo : public KUCMSObjectBase

KUCMSEventInfoObject::KUCMSEventInfoObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
	cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );
	cfFlag.set( "makeTriggerList", iConfig.existsAs<bool>("makeTriggerList") ? iConfig.getParameter<bool>("makeTriggerList") : false );

	if( iConfig.existsAs<std::vector<std::string>>("triggerList") ){ 
		triggerList = iConfig.getParameter<std::vector<std::string>>("triggerList");
	}//<<>>if( iConfig.existsAs<std::vector<std::string>>("triggerList") )

	fullTriggerList.clear();

}//<<>>KUCMSEventInfo::KUCMSEventInfo( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSEventInfoObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("erun","Evt_run",UINT,"Run number of event");
    Branches.makeBranch("elumi","Evt_luminosityBlock",UINT,"Lumi section of event");
    Branches.makeBranch("event","Evt_event",UINT);
    Branches.makeBranch("vnVtx","PV_npvs",UINT);
    Branches.makeBranch("vtxX","PV_x",FLOAT);
    Branches.makeBranch("vtxY","PV_y",FLOAT);
    Branches.makeBranch("vtxZ","PV_z",FLOAT);

    Branches.makeBranch("Flag_goodVertices","Flag_goodVertices",BOOL);
    Branches.makeBranch("Flag_globalSuperTightHalo2016Filter","Flag_globalSuperTightHalo2016Filter",BOOL);
    Branches.makeBranch("Flag_HBHENoiseFilter","Flag_HBHENoiseFilter",BOOL);
    Branches.makeBranch("Flag_HBHENoiseIsoFilter","Flag_HBHENoiseIsoFilter",BOOL);
    Branches.makeBranch("Flag_EcalDeadCellTriggerPrimitiveFilter","Flag_EcalDeadCellTriggerPrimitiveFilter",BOOL);
    Branches.makeBranch("Flag_BadPFMuonFilter","Flag_BadPFMuonFilter",BOOL);
    Branches.makeBranch("Flag_BadPFMuonDzFilter","Flag_BadPFMuonDzFilter",BOOL);
    Branches.makeBranch("Flag_hfNoisyHitsFilter","Flag_hfNoisyHitsFilter",BOOL);
    Branches.makeBranch("Flag_BadChargedCandidateFilter","Flag_BadChargedCandidateFilter",BOOL);
    Branches.makeBranch("Flag_eeBadScFilter","Flag_eeBadScFilter",BOOL);
    Branches.makeBranch("Flag_ecalBadCalibFilter","Flag_ecalBadCalibFilter",BOOL);

	for( auto trigName : triggerList ){

		std::string branchName = "Trigger_" + trigName;
		Branches.makeBranch( trigName, branchName, BOOL );
		//trigFlags[trigName] = false;

	}//<<>>for( auto trigName : triggerList )

    //Branches.makeBranch("triggerList","Evt_triggerList",VSTR,"Name of all triggers the event passed");

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSEventInfo::InitObject( TTree* fOutTree )

void KUCMSEventInfoObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // Load handels from tokens
    iEvent.getByToken( verticesToken, vertices_ );
    iEvent.getByToken( triggerResultsToken, triggerResults );
    iEvent.getByToken( triggerEventToken, triggerEvent );

    const edm::TriggerNames& triggerNames = iEvent.triggerNames( *triggerResults );
    const uInt nTriggerNames = triggerNames.size();
	flags.clear();
    //std::cout << " ---- Trigger Results :" << std::endl;
	for (auto itrig = 0U; itrig < nTriggerNames; itrig++){ 
		auto name = triggerNames.triggerName(itrig);
		auto result = triggerResults->accept(itrig);	
		flags[name] = result;
		//std::cout << " " << name << " = " << result << std::endl; 
	}//for (auto itrig = 0U; itrig < nTriggerNames; itrig++)
	
	trigFlags.clear();
	for( auto trigName : triggerList ){ trigFlags[trigName] = false; }

    //if( EventInfoDEBUG ) std::cout << "Collecting EventInfos" << std::endl;

	//edm::Hadndle<bool> goodVerticesFilter;
    //iEvent.getByLabel("metfilters_path",goodVerticesFilter);
	//std::cout << " -- goodVerticesFilter : " << *goodVerticesFilter << std::endl;

    eventVar.set("run",iEvent.id().run());
    eventVar.set("lumi", iEvent.luminosityBlock());
    eventVar.set("event", iEvent.id().event());

    if( EventInfoDEBUG ) std::cout << "Processing event: " << eventVar("event");
    if( EventInfoDEBUG ) std::cout << " in run: " << eventVar("run") << " and lumiblock: " << eventVar("lumi") << std::endl;

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

    if( EventInfoDEBUG ) std::cout << "Processing EventInfo" << std::endl;

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

    //if( flags.find("metfilters_path") != flags.end() ) std::cout << "metfilters_path = " << flags["metfilters_path"] << std::endl;
	//if( flags.find("Flag_goodVertices") != flags.end() ) Branches.fillBranch( "Flag_goodVertices", flags["Flag_goodVertices"] );
	fillFlagBranch("Flag_goodVertices");
    fillFlagBranch("Flag_globalSuperTightHalo2016Filter");
    fillFlagBranch("Flag_HBHENoiseFilter");
    fillFlagBranch("Flag_HBHENoiseIsoFilter");
    fillFlagBranch("Flag_EcalDeadCellTriggerPrimitiveFilter");
    fillFlagBranch("Flag_BadPFMuonFilter");
    fillFlagBranch("Flag_BadPFMuonDzFilter");
    fillFlagBranch("Flag_hfNoisyHitsFilter");
    fillFlagBranch("Flag_BadChargedCandidateFilter");
    fillFlagBranch("Flag_eeBadScFilter");
    fillFlagBranch("Flag_ecalBadCalibFilter");


	//std::cout << " -------------------------- Trigger Event :" << std::endl;
	std::vector<std::string> triggetList;
	const unsigned sizeFilters(triggerEvent->sizeFilters());
    for (size_t iF = 0; iF < sizeFilters; ++iF){
    
		const std::string nameFilter( triggerEvent->filterLabel(iF) );
		//triggetList.push_back(nameFilter);
		//Branches.fillBranch("triggerList",nameFilter);
		if( cfFlag("makeTriggerList") ){
			if( std::find( fullTriggerList.begin(), fullTriggerList.end(), nameFilter ) == fullTriggerList.end() ){
				fullTriggerList.push_back(nameFilter);
			}//<<>>if( std::find( fullTriggerList.begin(), fullTriggerList.end(), nameFilter ) == fullTriggerList.end() )
		}//<<>>if( cfFlag("makeTriggerList") )

		if( std::find( triggerList.begin(), triggerList.end(), nameFilter ) != triggerList.end() ){
			trigFlags[nameFilter] = true;
		}//<<>>if( std::find( triggerList.begin(), triggerList.end(), nameFilter ) == triggerList.end() )
		//std::cout << " : " << nameFilter << std::endl;

		//const trigger::Keys& keys = triggerEvent->filterKeys(iF);
		//const trigger::Vids& types = triggerEvent->filterIds(iF);
		//const unsigned nTriggers = types.size();
		//for (unsigned iTrig = 0; iTrig < nTriggers; ++iTrig){
		//	trigger::TriggerObject trigObj = trigObjs[keys[iTrig]];
		//	std::cout << " " << nameFilter << std::endl;
		//}//<<>>for (unsigned iTrig = 0; iTrig < nTriggers; ++iTrig)
    
	}//<<>>for (size_t iF = 0; iF < sizeFilters; ++iF)
	// this section access the actual trigger flags - move to load event section simular to flags with triggerResults

	for( auto trigName : triggerList ){ Branches.fillBranch( trigName, trigFlags[trigName] ); }
	//Branches.fillBranch("triggerList",triggetList);

}//<<>>void KUCMSEventInfo::ProcessEvent()

void KUCMSEventInfoObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSEventInfoObject::EndJobs(){

	if( cfFlag("makeTriggerList") ){
		
		std::cout << "List of Triggers appearing in the events processed in this dataset : " << std::endl;
		for( auto trigName : fullTriggerList ) std::cout << trigName << std::endl;

	}//<<>>if( cfFlag("makeTriggerList") )

}//<<>>void KUCMSEventInfoObject::EndJobs()

#endif
