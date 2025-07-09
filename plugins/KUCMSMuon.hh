// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Muon Object
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
#include "DataFormats/PatCandidates/interface/Muon.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSDisplacedVertex.hh"

#ifndef KUCMSMuonObjectHeader
#define KUCMSMuonObjectHeader

//#define MuonDEBUG true
#define MuonDEBUG false

using namespace edm; 

class KUCMSMuonObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSMuonObject( const edm::ParameterSet& iConfig );
    ~KUCMSMuonObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadMuonTokens( edm::EDGetTokenT<edm::View<reco::Muon>> muonsToken ){ muonsToken_ = muonsToken; }; 
    ////void LoadConversionTokens( edm::EDGetTokenT<reco::ConversionCollection> conversionsToken ){ conversionsToken_ = conversionsToken; };
    ////void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObject* genObjs_ ){ genObjs = genObjs_; };
    void LoadDisplacedVertexObject( KUCMSDisplacedVertex* svObj_ ){ svObj = svObj_; };

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
    //bool getMuonVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ); 
	//float getMuTrackZMatch( const reco::Photon & photon ); 
    int getIndex( float kideta, float kidphi );

    private:

    //std::vector<examnple> fexamples;
    std::vector<reco::Muon> fmuons;
    std::vector<int> muIds;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    //const edm::InputTag muonsTag;
    //edm::EDGetTokenT<std::vector<reco::GsfMuon> > muonsToken_;
    edm::EDGetTokenT<edm::View<reco::Muon>> muonsToken_;
    //edm::Handle<std::vector<reco::GsfMuon> > muons_;
    edm::Handle<edm::View<reco::Muon> > muons_;
    //std::vector<reco::GsfMuon> muons;

    //const edm::InputTag muMVAIDLooseMapTag;
    //edm::EDGetTokenT<edm::ValueMap<bool>> muMVAIDLooseMapToken_;
    //edm::Handle<edm::ValueMap<bool>> muMVAIDLooseMap_;

    // Conversions ( for muon veto ) ---------------------------------
    //const edm::InputTag conversionsTag;
    //edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    //edm::Handle<reco::ConversionCollection> conversions_;

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    //edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    //edm::Handle<reco::BeamSpot> beamSpot_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;
    KUCMSEcalRecHitObject* rhObj;
    KUCMSGenObject* genObjs;
    KUCMSDisplacedVertex* svObj;

};//<<>>class KUCMSMuon : public KUCMSObjectBase

KUCMSMuonObject::KUCMSMuonObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    //cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minMuE", iConfig.existsAs<double>("minMuE") ? iConfig.getParameter<double>("minMuE") : 2.0 );
    //cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );

}//<<>>KUCMSMuon::KUCMSMuon( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSMuonObject::InitObject( TTree* fOutTree ){

    //Branches.makeBranch("RhIds","Muon_rhIds",VVUINT);
    Branches.makeBranch("nMu","Muon_nMuons",INT);
    //Branches.makeBranch("scIndex","Muon_scIndex",VINT,"index of supercluster");
    //Branches.makeBranch("SeedTOFTime","Muon_seedTOFTime",VFLOAT);
    Branches.makeBranch("Pt","Muon_pt",VFLOAT);
    Branches.makeBranch("Energy","Muon_energy",VFLOAT);
    Branches.makeBranch("Phi","Muon_phi",VFLOAT);
    Branches.makeBranch("Eta","Muon_eta",VFLOAT);
    Branches.makeBranch("Px","Muon_px",VFLOAT);
    Branches.makeBranch("Py","Muon_py",VFLOAT);
    Branches.makeBranch("Pz","Muon_pz",VFLOAT);
    //Branches.makeBranch("Time","Muon_time",VFLOAT);
    //Branches.makeBranch("timeValid","Muon_timeValid",VBOOL);
    //Branches.makeBranch("TrackZ","Muon_trackz",VFLOAT);
    Branches.makeBranch("looseID","Muon_isLoose",VBOOL);
    Branches.makeBranch("mediumID","Muon_isMedium",VBOOL);
    Branches.makeBranch("nSMu","Muon_nSelMuons",INT);
    Branches.makeBranch("GenIdx","Muon_genIdx",VINT);
    //Branches.makeBranch("GenXMomIdx","Muon_genSigXMomId",VINT);
    //Branches.makeBranch("GenWZIdx","Muon_genSigWZId",VINT);
    //Branches.makeBranch("Sieie","Muon_Sieie",VFLOAT);
    //Branches.makeBranch("DetaSCTV","Muon_DetaSCTV",VFLOAT);
    //Branches.makeBranch("DphiSCTV","Muon_DphiSCTV",VFLOAT);
    //Branches.makeBranch("HOE","Muon_HOE",VFLOAT);
    Branches.makeBranch("svMatch","Muon_hasSVMatch",VBOOL);
    Branches.makeBranch("nSVMatch","Muon_nSVMatched",INT);
    Branches.makeBranch("IsPrompt","Muon_isPrompt",VBOOL);
    Branches.makeBranch("nPrompt","Muon_nPrompt",INT);

    //Branches.makeBranch("GenDr","Muon_genDr",VFLOAT);
    //Branches.makeBranch("GenDp","Muon_genDp",VFLOAT);
    //Branches.makeBranch("GenSIdx","Muon_genSIdx",VINT);
    //Branches.makeBranch("GenSDr","Muon_genSDr",VFLOAT);
    //Branches.makeBranch("GenSDr","Muon_genSDp",VFLOAT);
    //Branches.makeBranch("GenLlpId","Muon_genLlpId",VFLOAT);
    //Branches.makeBranch("GenSLlpId","Muon_genSLlpId",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSMuon::InitObject( TTree* fOutTree )

void KUCMSMuonObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){


	iEvent.getByToken(muonsToken_, muons_);
    //iEvent.getByToken(conversionsToken_, conversions_);
    //iEvent.getByToken(beamLineToken_, beamSpot_);

    if( MuonDEBUG ) std::cout << "Collecting Muons" << std::endl;

    fmuons.clear();
    muIds.clear();// indexed by mu index ( 0,1,2 ) * number of ids ( 1 current, 6? possible ) + index of ID wanted
    for (edm::View<reco::Muon>::const_iterator itMuon = muons_->begin(); itMuon != muons_->end(); itMuon++) {
        //auto idx = itMuon - muons_->begin();//unsigned int
        //auto muonRef = muons_->refAt(idx);//edm::RefToBase<reco::GsfMuon> 
        auto &muon = (*itMuon);
        //auto passIdCut = true; //muon.muonID(muCutLoose);// pat muon ( miniAOD ) method
        //muIdBools.push_back((*muMVAIDLooseMap_)[muonRef]);// order is important, track how this vector is loaded
        auto passEnergyCut = muon.energy() > cfPrm("minMuE");
        if( passEnergyCut ) fmuons.push_back(muon);
    }//<<>>for( const auto muon : *muons_ )


}//<<>>void KUCMSMuon::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSMuonObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSMuonObject::ProcessEvent( ItemManager<float>& geVar ){

    if( MuonDEBUG ) std::cout << "Processing Muons" << std::endl;

    Branches.clearBranches();

    if( MuonDEBUG ) std::cout << " - enetering Muon loop" << std::endl;

	int muIndx = 0;
	int nSelMu = 0;
    int nSVMatched = 0;
	int nPromptMuon = 0;
    //scGroup scptrs;
    //std::vector<float> scptres;
    if( MuonDEBUG ) std::cout << "Processing Muons" << std::endl;
    for( const auto &muon : fmuons ){

        const float muPt = muon.pt();
        const float muEnergy = muon.energy();
        const float muPhi = muon.phi();
        const float muEta = muon.eta();
        const float muPx = muon.px();
        const float muPy = muon.py();
        const float muPz = muon.pz();
		//const float time = muon.time();
		//const bool timeValid = muon.isTimeValid();
		const bool isLoose = muon::isLooseMuon( muon );
		const bool isMedium = muon::isMediumMuon( muon );

		if( MuonDEBUG ) std::cout << " -- " << muPt << " - " << muEta << " - " << isLoose << " - " << isMedium << std::endl;	
        ////const float muTrackZ = muon.trackPositionAtVtx().Z();

        Branches.fillBranch("Pt",muPt);
        Branches.fillBranch("Energy",muEnergy);
        Branches.fillBranch("Phi",muPhi);
        Branches.fillBranch("Eta",muEta);
        Branches.fillBranch("Px",muPx);
        Branches.fillBranch("Py",muPy);
        Branches.fillBranch("Pz",muPz);
    	//Branches.fillBranch("Time",time);
		Branches.fillBranch("looseID",isLoose);
        Branches.fillBranch("mediumID",isMedium);
		////Branches.fillBranch("TrackZ",muTrackZ);
    	//Branches.fillBranch("timeValid","Muon_timeValid",VBOOL);

        if( MuonDEBUG ) std::cout << " --- Proccesssing : " << muon << std::endl;

        bool hasSVMatch = false;
        if( cfFlag("doSVModule") ) hasSVMatch = svObj->FoundLeptonMatch( muon );
        if( hasSVMatch ) nSVMatched++;
        Branches.fillBranch("svMatch",hasSVMatch);

        bool isPrompt = false;
        if( cfFlag("doSVModule") ) isPrompt = svObj->IsPromptLepton( muon );

        bool mptc = muPt >= 20 ;
        if( mptc && isLoose ) nSelMu++;

        bool isPromptMuon = false;
        if( isPrompt && isLoose && mptc ){ isPromptMuon = true; nPromptMuon++; }
        Branches.fillBranch("IsPrompt",isPromptMuon);

		int genIndex = -9;
		if( cfFlag("hasGenInfo") ) genIndex = genObjs->getGenParticleIndex( muon, 13 ); 
		Branches.fillBranch("GenIdx", genIndex );

	    muIndx++;

    }//<<>>for( const auto muon : *muons_ )
	Branches.fillBranch("nMu",muIndx);
    Branches.fillBranch("nSMu",nSelMu);
    Branches.fillBranch("nSVMatch",nSVMatched);
	Branches.fillBranch("nPrompt",nPromptMuon);
	geVar.set("nSelMu",nSelMu);


}//<<>>void KUCMSMuon::ProcessEvent()

void KUCMSMuonObject::EndJobs(){}

//float KUCMSMuonObject::getMuTrackZMatch( const reco::Photon & photon ){
//
//	float bestmatch(0.1);
//	float elTrackZ(1000.0);
//	for (const auto & muon : fmuons){
//		auto match = reco::deltaR(muon,photon);
//		if ( match < bestmatch ){ elTrackZ = muon.trackPositionAtVtx().Z(); bestmatch = match; }
//	}//<<>>for (const auto muon : fmuons)
//	return elTrackZ;
//
//}//<<>>float KUCMSMuonObject::getMuTrackZMatch( const reco::Photon & photon )

//bool KUCMSMuonObject::getMuonVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ){
//
//	return conversions_.isValid() ?
//            not ConversionTools::hasMatchedPromptMuon( phosc, fmuons, *conversions_, beamSpot_->position() ) : false;

//}//<<>> bool KUCMSMuonObject::getMuonVeto( const reco::SuperCluster & phosc )

int KUCMSMuonObject::getIndex( float kideta, float kidphi ){

    auto bestdr = 0.01;
    int idx = -1;
    int cnt = 0;
    for( const auto &muon : fmuons  ){
        auto mueta = muon.eta();
        auto muphi = muon.phi();
        auto kidpfdr = std::sqrt( reco::deltaR2( mueta, muphi, kideta, kidphi ));
        if( kidpfdr < bestdr ){ bestdr = kidpfdr; idx = cnt; }
        cnt++;
    }//<<>>for( const auto muon : fmuons  )

	return idx;

}//<<>>uInt KUCMSMuonObject::getIndex( float kideta, float kidphi )

#endif
