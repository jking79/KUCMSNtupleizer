// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// Object Manager and
// KUCMS Base Object for Object Manager
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// basic C++ types
#include <cstdlib>
#include <vector>
#include <map>
#include <string>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// DataFormats
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

// TOOLS
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TFormula.h"
#include "TF1.h"
#include "TTree.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

// KUCMS
#include "KUCMSItemManager.hh"
#include "KUCMSBranchManager2.hh"

#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

#ifndef KUCMSObjectBaseHeader
#define KUCMSObjectBaseHeader

// This file contains both ObjectManager and KUCMSObjectBase

//.............................................................................................
// KUCMSObjectBase----------------------------------------------------------------------------
//.............................................................................................

typedef edm::View<reco::Candidate> CandidateView;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;

#define BUNCHES 3564
#define SAMPLES 10

class KUCMSObjectBase : public edm::one::EDAnalyzer<edm::one::SharedResources> {

    public:

    KUCMSObjectBase(){}; 
    virtual ~KUCMSObjectBase(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    virtual void InitObject( TTree* fOutTree ){}; // sets up branches, do preloop jobs 
    // Load tokens for collections derived from iEvent 
	// void LoadSCToken( edm::EDGetTokenT<reco::SuperClusterCollection> scToken_ ); <- example code
    // crossload functions are to be declared in derived classes
    // void LoadObject( KUCMSObjectBase* otherObject ){ otherObjectPtr = otherObject}; // overload with specific KUCMS object need here

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    virtual void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){}; 
    // do cross talk jobs with other objects, do event processing, and load branches
    virtual void ProcessEvent( ItemManager<float>& geVar ){}; 

    // if there are any final tasks be to done after the event loop via objectManager
    virtual void EndJobs(){}; // do any jobs that need to be done after main event loop

    void setName( std::string name ){ title = name; };
    std::string getName(){ return title; }; 

    protected:

    // input parameters
    ItemManager<double> cfPrm; // confguration paramters loaded in object constructor step from iConfig
	ItemManager<bool> cfFlag;    

    // Branch Manager
    KUCMSBranchManager Branches;

	// other object pointers loaded with void LoadObject( KUCMSObjectBase* otherObject )
	// KUCMSObjectBase* otherObjectPtr

    // Object Name
    std::string title;

	private:

    void analyze( const edm::Event& iEvent, const edm::EventSetup& iSetup)
		{ std::cout << "Error!!!! analyze called in KUCMSObjectBase" << std::endl; }	

};

// ---------------------------------------------------------------------------------------
// ObjectManager -------------------------------------------------------------------------
// ---------------------------------------------------------------------------------------

class ObjectManager {

    public:

    ~ObjectManager();

    void Load( std::string key, std::string name, KUCMSObjectBase* object );
    void Load( std::string name, KUCMSObjectBase* object );
	void Init( TTree* fOutTree ); 
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar );
    void ProcessEvent( ItemManager<float>& geVar );
    void EndJobs();

    private:

    std::map< std::string, KUCMSObjectBase* > objects;
    bool valid( std::string key );

};

ObjectManager::~ObjectManager(){

    for( auto & object : objects ){ delete object.second; }

}//<<>>ObjectManager::~ObjectManager()

void ObjectManager::Load( std::string key, std::string name, KUCMSObjectBase* object  ){

    if( objects.find(key) == objects.end() ){
        object->setName(name);
        objects.insert( {key,object} );
    } else std::cout << " -- Error: " << key << "  already exists !!!!!!!!!!!!!!!!!!!" << std::endl;

}//>><<void ObjectManager::set( std::string key, std::string name, std::string doc )

void ObjectManager::Load( std::string name, KUCMSObjectBase* object ){ Load( name, name, object ); }

void ObjectManager::Init( TTree* fOutTree ){

    for( auto & object : objects ){ object.second->InitObject( fOutTree ); }

}//<<>>void ObjectManager::Init( ItemManager<float>& geVar, TTree* fOutTree )

void ObjectManager::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    for( auto & object : objects ){ object.second->LoadEvent( iEvent, iSetup, geVar ); }

}//<<>>void ObjectManager::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void ObjectManager::ProcessEvent( ItemManager<float>& geVar ){ for( auto & object : objects ){ object.second->ProcessEvent( geVar ); } }

void ObjectManager::EndJobs(){ for( auto & object : objects ){ object.second->EndJobs(); } }

bool ObjectManager::valid( std::string key ){

    if( objects.find(key) == objects.end() ){
        std::cout << " -- Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return false;
    } else return true;

}//<<>>bool ObjectManager::valid( std::string key )

#endif
