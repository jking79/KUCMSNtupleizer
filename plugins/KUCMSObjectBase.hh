// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Object for Object Manager
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

//.............................................................................................

typedef edm::View<reco::Candidate> CandidateView;

typedef unsigned int uInt;
typedef unsigned long int ulInt;
typedef const float CFlt;
typedef const double CDbl;
typedef const vector<float> CVFlt;
typedef const vector<double> CVDbl;

typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;

#define nEBEEMaps 0
#define nHists 4
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ...  
#define BUNCHES 3564
#define SAMPLES 10

const auto sortByPt = [](auto & obj1, auto & obj2) {return obj1.pt() > obj2.pt();};

class KUCMSObjectBase {

    public:

    KUCMSObjectBase( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag ) : iConfig(iConfig), cfFlag(cfFlag) {}; 
    virtual ~KUCMSObjectBase(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    virtual void InitObject( ItemManager<float>& geVar, TTree* fOutTree ); // load geVars and sets up branches, do preloop jobs 
    // crossload functions are to be declared in derived classes
    // void LoadObject( KUCMSObjectBase& otherObject ){ otherObjectRef = otherObject}; // overload with specific KUCMS object need here

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    virtual void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup ); // get collections, do initial processing
    virtual void ProcessEvent(); // do cross talk jobs with other objects, do event processing, and load branches

    // if there are any final tasks be to done after the event loop via objectManager
    virtual void EndJobs(); // do any jobs that need to be done after main event loop

    void setName( std::string name ){ title = name; };
    std::string getName(){ return title; }; 

    private:

	const edm::ParameterSet iConfig;
    const ItemManager<bool> cfFlag

    // input parameters
    ItemManager<double> cfPrm; // confguration paramters loaded in object constructor step from iConfig
    
    // global event varibles
    std::reference_wrapper<ItemManager<float>> eventvars; // used to pass event level inforamtion

    // Branch Manager
    KUCMSBranchManager branches;

    // Object Name
    std::string title;

};

class exampleObject : public KUCMSObjectBase {};

class DummyObject : public KUCMSObjectBase {

    // use base class constructor
    using KUCMSObjectBase::KUCMSObjectBase;

    public:

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    void InitObject( ItemManager<float>& geVar, TTree* fOutTree ); // load geVars and sets up branches, do preloop jobs
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    void LoadObject( exampleObject& otherObject ){ otherObjectRef = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup ); // get collections, do initial processing
    void ProcessEvent(); // do cross talk jobs with other objects, do event processing, and load branches

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    private: // example from photons

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    std::reference_wrapper<exampleObject> otherObjectRef;

};

class ObjectManager {

    public:

    void Load( std::string key, std::string name, KUCMSObjectBase& object );
    void Load( std::string name, KUCMSObjectBase& object );
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup );
    void ProcessEvent();
    void EndJobs();

    private:

    std::map< std::string, KUCMSObjectBase& > objects;
    bool valid( std::string key );

};

void ObjectManager::Load( std::string key, std::string name, KUCMSObjectBase& object  ){

    if( objects.find(key) == objects.end() ){
        object.setName(name);
        objects.insert( {key,object} );
    } else std::cout << " -- Error: " << key << "  already exists !!!!!!!!!!!!!!!!!!!" << std::endl;

}//>><<void ObjectManager::set( std::string key, std::string name, std::string doc )

void ObjectManager::Load( std::string name, KUCMSObjectBase& object ){ Load( name, name, object ); }

void ObjectManager::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup ){

    for( auto & object : objects ){ object.second.LoadEvent( iEvent, iSetup ); }

}//<<>>void ObjectManager::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void ObjectManager::ProcessEvent(){ for( auto & object : objects ){ object.second.ProcessEvent(); } }

void ObjectManager::EndJobs(){ for( auto & object : objects ){ object.second.EndJobs(); } }

bool ObjectManager::valid( std::string key ){

    if( objects.find(key) == objects.end() ){
        std::cout << " -- Error: No Such Key : " << key << " !!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
        return false;
    } else return true;

}//<<>>bool ObjectManager::valid( std::string key )

#endif
