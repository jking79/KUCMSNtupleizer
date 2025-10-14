// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS MET Filters Object
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

// Includes for all object types

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/Exception.h"

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
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
//#include "KUCMSEcalRechit.hh"
//#include "KUCMSGenObjects.hh"

#ifndef KUCMSEventFilterObjectHeader
#define KUCMSEventFilterObjectHeader

#define EventFilterDEBUG true
//#define EventFilterDEBUG false

using namespace edm; 

class KUCMSEventFilterObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSEventFilterObject( const edm::ParameterSet& iConfig );
    ~KUCMSEventFilterObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadHcalNoiseToken( edm::EDGetTokenT<HcalNoiseSummary> noiseToken){ noiseToken_ = noiseToken; }
//    void LoadElectronTokens( edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken ){ electronsToken_ = electronsToken; }; 
   void LoadBunchSpacing( edm::EDGetTokenT<unsigned int> bunchSpacing ) { bunchSpacing_ = bunchSpacing; } 

    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    //void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
   // void LoadGenObject( KUCMSGenObject* genObjs_ ){ genObjs = genObjs_; };

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
    //bool getElectronVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ); 
    //int getIndex( float kideta, float kidphi );

    private:

    //std::vector<examnple> fexamples;
    
    //HBHENoiseFilter
    bool HBHENoiseFilterResult_;
    edm::EDGetTokenT<HcalNoiseSummary> noiseToken_;
    edm::Handle<HcalNoiseSummary> summary_h;
    edm::EDGetTokenT<unsigned int> bunchSpacing_; 
    edm::Handle<unsigned int> bunchSpacingH_;    

    std::map<std::string, bool> decisionMap_;
  


    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    
    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;
    //KUCMSEcalRecHitObject* rhObj;
    //KUCMSGenObject* genObjs;

};//<<>>class KUCMSEventFilter : public KUCMSObjectBase

KUCMSEventFilterObject::KUCMSEventFilterObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------
   //default values are user set - these are from the HBHENoiseFilterResultProducer_cfi.py config file in CMSSW
    cfInt.set("minHPDHits",iConfig.existsAs<int>("minHPDHits") ? iConfig.getParameter<int>("minHPDHits") : 17);
    cfInt.set("minHPDNoOtherHits",iConfig.existsAs<int>("minHPDNoOtherHits") ? iConfig.getParameter<int>("minHPDNoOtherHits") : 10);
    cfInt.set("minZeros",iConfig.existsAs<int>("minZeros") ? iConfig.getParameter<int>("minZeros") : 9999);

    cfFlag.set( "IgnoreTS4TS5ifJetInLowBVRegion_", iConfig.existsAs<bool>("IgnoreTS4TS5ifJetInLowBVRegion") ? iConfig.getParameter<bool>("IgnoreTS4TS5ifJetInLowBVRegion") : true); 
    cfStr.set( "defaultHBHEDecision", iConfig.existsAs<std::string>("defaultDecision") ? iConfig.getParameter<std::string>("defaultDecision") : "HBHENoiseFilterResultRun2Loose");

    cfInt.set("minNumIsolatedNoiseChannels",iConfig.existsAs<int>("minNumIsolatedNoiseChannels") ? iConfig.getParameter<int>("minNumIsolatedNoiseChannels") : 10);
    cfPrm.set("minIsolatedNoiseSumE",iConfig.existsAs<double>("minIsolatedNoiseSumE") ? iConfig.getParameter<double>("minIsolatedNoiseSumE") : 50.0);
    cfPrm.set("minIsolatedNoiseSumEt",iConfig.existsAs<double>("minIsolatedNoiseSumEt") ? iConfig.getParameter<double>("minIsolatedNoiseSumEt") : 25.0);

    cfFlag.set( "useBranchSpacingProducer", iConfig.existsAs<bool>("useBunchSpacingProducer") ? iConfig.getParameter<bool>("useBunchSpacingProducer") : true);



}//<<>>KUCMSEventFilter::KUCMSEventFilter( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSEventFilterObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("HBHENoiseFilter","HBHENoiseFilter",VBOOL);
    //Branches.makeBranch("Pt","Electron_pt",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSEventFilter::InitObject( TTree* fOutTree )

void KUCMSEventFilterObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    if( EventFilterDEBUG ) std::cout << "Collecting EventFilters" << std::endl;
   // get tokens for event
    iEvent.getByToken(noiseToken_, summary_h);
 if (!summary_h.isValid()) {
    throw edm::Exception(edm::errors::ProductNotFound) << " could not find HcalNoiseSummary.\n";
    return;
  }
   //set decision map based on HcalNoiseSummary
    const HcalNoiseSummary& summary(*summary_h);
    bool goodJetFoundInLowBVRegion = false;
     if (cfFlag("IgnoreTS4TS5ifJetInLowBVRegion"))
       goodJetFoundInLowBVRegion = summary.goodJetFoundInLowBVRegion();

    const bool failCommon = summary.maxHPDHits() >= cfInt("minHPDHits") || summary.maxHPDNoOtherHits() >= cfInt("minHPDNoOtherHits") || summary.maxZeros() >= cfInt("minZeros");

    const bool failRun1 = failCommon || (summary.HasBadRBXTS4TS5() && !goodJetFoundInLowBVRegion);
    decisionMap_["HBHENoiseFilterResultRun1"] = failRun1;

    const bool failRun2Loose = failCommon || (summary.HasBadRBXRechitR45Loose() && !goodJetFoundInLowBVRegion);
    decisionMap_["HBHENoiseFilterResultRun2Loose"] = failRun2Loose;
    
    const bool failRun2Tight = failCommon || (summary.HasBadRBXRechitR45Tight() && !goodJetFoundInLowBVRegion);
    decisionMap_["HBHENoiseFilterResultRun2Tight"] = failRun2Tight;


    // Overwrite defaultDecision_ dynamically based on bunchSpacingProducer
    if (cfFlag("useBunchSpacingProducer")) {
      iEvent.getByToken(bunchSpacing_, bunchSpacingH_);
      unsigned int bunchspacing = 0;
      if (bunchSpacingH_.isValid()) {
        bunchspacing = *bunchSpacingH_;
        if (bunchspacing == 50) {
          cfStr.set("defaultHBHEDecision","HBHENoiseFilterResultRun1");
        } else {
          cfStr.set("defaultHBHEDecision","HBHENoiseFilterResultRun2Loose");
        }
      }
    }

 
}//<<>>void KUCMSEventFilter::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSEventFilterObject::ProcessEvent( ItemManager<float>& geVar ){

    if( EventFilterDEBUG ) std::cout << "Processing EventFilters" << std::endl;

    Branches.clearBranches();

    std::map<std::string, bool>::const_iterator it = decisionMap_.find(cfStr("defaultHBHEDecision"));
    if( EventFilterDEBUG) std::cout << "HBHENoiseFilterResultProducer parameter: " << cfStr("defaultHBHEDecision") << "\n" << std::endl;
    if(it == decisionMap_.end()){
      throw cms::Exception("Invalid HBHENoiseFilterResultProducer parameter"); 
    }
    const bool HBHENoiseFilter = !it->second;
    Branches.fillBranch("HBHENoiseFilter",HBHENoiseFilter);

}//<<>>void KUCMSEventFilter::ProcessEvent()

void KUCMSEventFilterObject::EndJobs(){}

#endif
