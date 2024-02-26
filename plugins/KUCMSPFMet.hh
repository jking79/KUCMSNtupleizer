// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS PFMet Object
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
//#include "DataFormats/PFMet/interface/PFMet.h"

#include "DataFormats/PatCandidates/interface/MET.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

#ifndef KUCMSPFMetObjectHeader
#define KUCMSPFMetObjectHeader

//#define PFMetDEBUG true
#define PFMetDEBUG false

using namespace edm; 

class KUCMSPFMetObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSPFMetObject( const edm::ParameterSet& iConfig );
    ~KUCMSPFMetObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadPFMetTokens( edm::EDGetTokenT<std::vector<reco::PFMET>> pfmetsToken ){ pfmetsToken_ = pfmetsToken; }; 
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
	void LoadPhotonObject( KUCMSPhotonObject* otherObject ){ phoObj = otherObject; };

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

    //std::vector<examnple> fexamples;
    std::vector<reco::PFMET> fpfmet;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    //const edm::InputTag pfmetsTag;
    edm::EDGetTokenT<std::vector<reco::PFMET>> pfmetsToken_;
    edm::Handle<std::vector<reco::PFMET>> pfmets_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;
    KUCMSPhotonObject* phoObj;

};//<<>>class KUCMSPFMet : public KUCMSObjectBase

KUCMSPFMetObject::KUCMSPFMetObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );

}//<<>>KUCMSPFMet::KUCMSPFMet( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSPFMetObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("SumEt","Met_sumEt",FLOAT,"scalar sum of Et");
    Branches.makeBranch("Pt","Met_pt",FLOAT,"pt");
    Branches.makeBranch("Eta","Met_eta",FLOAT);
    Branches.makeBranch("Phi","Met_phi",FLOAT);
    Branches.makeBranch("Px","Met_px",FLOAT);
    Branches.makeBranch("Py","Met_py",FLOAT);

    //Branches.makeBranch("signif","Met_significance",FLOAT,"MET significance");
    //Branches.makeBranch("covXX","Met_covXX",FLOAT,"xx element of met covariance matrix");
    //Branches.makeBranch("covXY","Met_covXY",FLOAT,"xy element of met covariance matrix");
    //Branches.makeBranch("covYY","Met_covYY",FLOAT,"yy element of met covariance matrix");

    Branches.makeBranch("CPt","Met_CPt",FLOAT,"Met Pt corrected for oootPhotons");
    Branches.makeBranch("CSumEt","Met_CsumEt",FLOAT,"SumEt corrected for ootPhotons");
    Branches.makeBranch("CPx","Met_Cpx",FLOAT,"Met Px corrected for oootPhotons");
    Branches.makeBranch("CPy","Met_Cpy",FLOAT,"Met Py corrected for oootPhotons");

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSPFMet::InitObject( TTree* fOutTree )

void KUCMSPFMetObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS
    iEvent.getByToken( pfmetsToken_, pfmets_ );

    if( PFMetDEBUG ) std::cout << "Collecting PFMets" << std::endl;

    fpfmet.clear();
    //for( const auto met : *pfmets_ ){ if( met.energy() > minMetE ) 
    fpfmet.push_back((*pfmets_)[0]);

}//<<>>void KUCMSPFMet::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSPFMetObject::ProcessEvent( ItemManager<float>& geVar ){}

void KUCMSPFMetObject::PostProcessEvent( ItemManager<float>& geVar ){

    if( PFMetDEBUG ) std::cout << "Processing PFMets" << std::endl;

    Branches.clearBranches();

    if( PFMetDEBUG ) std::cout << " - enetering PFMet loop" << std::endl;

    auto t1pfMET = fpfmet[0];
    const float SumEt = t1pfMET.sumEt();
    const float Pt = t1pfMET.pt();
    const float Px = t1pfMET.px();
    const float Py = t1pfMET.py();
    const float Phi = t1pfMET.phi();

    //const float signif = t1pfMET.significance();
    //const float covXX = t1pfMET.getSignificanceMatrix().At(0,0);
    //const float covXY = t1pfMET.getSignificanceMatrix().At(0,1);
    //const float covYY = t1pfMET.getSignificanceMatrix().At(1,1);

    float CSumEt = SumEt;
    float CPx = Px;
    float CPy = Py;

    phoObj->correctMet( CSumEt, CPx, CPy );

    auto evtMET = std::sqrt( sq2(CPx) + sq2(CPy));
    geVar.set( "evtMET", evtMET );

    Branches.fillBranch("SumEt",SumEt);
    Branches.fillBranch("Pt",Pt);
    Branches.fillBranch("Phi",Phi);
    Branches.fillBranch("Px",Px);
    Branches.fillBranch("Py",Py);

    //Branches.fillBranch("signif",signif);
    //Branches.fillBranch("covXX",covXX);
    //Branches.fillBranch("covXY",covXY);
    //Branches.fillBranch("covYY",covYY);

    Branches.fillBranch("CPt",evtMET);;
    Branches.fillBranch("CSumEt",CSumEt);
    Branches.fillBranch("CPx",CPx);
    Branches.fillBranch("CPy",CPy);

}//<<>>void KUCMSPFMet::ProcessEvent()

void KUCMSPFMetObject::EndJobs(){}

#endif
