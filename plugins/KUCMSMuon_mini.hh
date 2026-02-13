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
#include "KUCMSEcalRechit_mini.hh"
#include "KUCMSGenObjects_mini.hh"
#include "KUCMSDisplacedVertex_mini.hh"

#ifndef KUCMSMuonObjectMiniHeader
#define KUCMSMuonObjectMiniHeader

//#define MuonDEBUG true
#define MuonDEBUG false

using namespace edm; 

class KUCMSMuonObjectMini : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSMuonObjectMini( const edm::ParameterSet& iConfig );
    ~KUCMSMuonObjectMini(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadMuonTokens( edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken ){ muonsToken_ = muonsToken; }; 
    void LoadVertexTokens( edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_ ){ verticesToken = verticesToken_; };

    ////void LoadConversionTokens( edm::EDGetTokenT<reco::ConversionCollection> conversionsToken ){ conversionsToken_ = conversionsToken; };
    ////void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    void LoadRecHitObject( KUCMSEcalRecHitObjectMini* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObjectMini* genObjs_ ){ genObjs = genObjs_; };
    void LoadDisplacedVertexObject( KUCMSDisplacedVertexMini* svObj_ ){ svObj = svObj_; };
    void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }

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
    std::vector<pat::Muon> fmuons;
    std::vector<int> muIds;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    //const edm::InputTag muonsTag;
    //edm::EDGetTokenT<std::vector<reco::GsfMuon> > muonsToken_;
    edm::EDGetTokenT<edm::View<pat::Muon>> muonsToken_;
    //edm::Handle<std::vector<reco::GsfMuon> > muons_;
    edm::Handle<edm::View<pat::Muon> > muons_;
    //std::vector<reco::GsfMuon> muons;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken;
    edm::Handle<std::vector<reco::Vertex>> vertices_;

	//Tracks
    //TransientTrackBuilder* ttBuilder;
	TTBuilderWrapper ttBuilder;

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
    KUCMSEcalRecHitObjectMini* rhObj;
    KUCMSGenObjectMini* genObjs;
    KUCMSDisplacedVertexMini* svObj;

};//<<>>class KUCMSMuon : public KUCMSObjectBase

KUCMSMuonObjectMini::KUCMSMuonObjectMini( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    //cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minMuE", iConfig.existsAs<double>("minMuE") ? iConfig.getParameter<double>("minMuE") : 2.0 );
    cfPrm.set( "maxEta",iConfig.existsAs<double>("maxEta")? iConfig.getParameter<double>("maxEta") : 2.44 );
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );

}//<<>>KUCMSMuon::KUCMSMuon( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSMuonObjectMini::InitObject( TTree* fOutTree ){

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
    Branches.makeBranch("tightID","Muon_isTight",VBOOL);
    Branches.makeBranch("nSMu","Muon_nSelMuons",INT);

	if( cfFlag("hasGenInfo") ){
	Branches.makeBranch("isLWZP", "Muon_isLWZP", VBOOL );
	Branches.makeBranch("GenIdx","Muon_genIdx",VINT);
	Branches.makeBranch("GenWZIdx","Muon_genSigWZId",VINT);
	Branches.makeBranch("GenXMomIdx","Muon_genSigXMomId",VINT);
    Branches.makeBranch("hasGenIdx","Muon_hasGenMatch",VBOOL);
	}//<<>>if( cfFlag("hasGenInfo") )
    Branches.makeBranch("nGenMatch","Muon_nGenMatch",INT);

    //Branches.makeBranch("Sieie","Muon_Sieie",VFLOAT);
    //Branches.makeBranch("DetaSCTV","Muon_DetaSCTV",VFLOAT);
    //Branches.makeBranch("DphiSCTV","Muon_DphiSCTV",VFLOAT);
    //Branches.makeBranch("HOE","Muon_HOE",VFLOAT);
    Branches.makeBranch("svMatch","Muon_hasSVMatch",VBOOL);
    Branches.makeBranch("nSVMatch","Muon_nSVMatched",INT);
    Branches.makeBranch("passPrompt","Muon_passPrompt",VBOOL);
    //Branches.makeBranch("hasMinPt","Muon_hasMinPt",VBOOL);
    Branches.makeBranch("IsPromptMuon","Muon_isPromptMuon",VBOOL);
    Branches.makeBranch("nPrompt","Muon_nPrompt",INT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSMuon::InitObject( TTree* fOutTree )

void KUCMSMuonObjectMini::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){


	iEvent.getByToken(muonsToken_, muons_);
    //iEvent.getByToken(conversionsToken_, conversions_);
    //iEvent.getByToken(beamLineToken_, beamSpot_);
    iEvent.getByToken( verticesToken, vertices_ );

	ttBuilder = TTBuilderWrapper( &iSetup.getData(transientTrackBuilder_) );

    if( MuonDEBUG ) std::cout << "Collecting Muons" << std::endl;

    fmuons.clear();
    muIds.clear();// indexed by mu index ( 0,1,2 ) * number of ids ( 1 current, 6? possible ) + index of ID wanted
    for (edm::View<pat::Muon>::const_iterator itMuon = muons_->begin(); itMuon != muons_->end(); itMuon++) {
        //auto idx = itMuon - muons_->begin();//unsigned int
        //auto muonRef = muons_->refAt(idx);//edm::RefToBase<reco::GsfMuon> 
        auto &muon = (*itMuon);
        //auto passIdCut = true; //muon.muonID(muCutLoose);// pat muon ( miniAOD ) method
        //muIdBools.push_back((*muMVAIDLooseMap_)[muonRef]);// order is important, track how this vector is loaded
        auto passEnergyCut = muon.energy() > cfPrm("minMuE");
		auto passEtaCut = std::abs( muon.eta() ) < cfPrm("maxEta");
        if( passEnergyCut && passEtaCut ) fmuons.push_back(muon);
    }//<<>>for( const auto muon : *muons_ )


}//<<>>void KUCMSMuon::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSMuonObjectMini::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSMuonObjectMini::ProcessEvent( ItemManager<float>& geVar ){

    if( MuonDEBUG ) std::cout << "Processing Muons" << std::endl;

    Branches.clearBranches();

    if( MuonDEBUG ) std::cout << " - enetering Muon loop" << std::endl;

	int muIndx = 0;
	int nSelMu = 0;
    int nSVMatched = 0;
	int nPromptMuon = 0;
    std::vector<v3fPoint> scvertex;
    std::vector<float> scptres;
	std::vector<reco::TransientTrack> sctrks;
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
		//const bool isLoose = muon::isLooseMuon( muon );
		const bool isLoose = muon.passed( reco::Muon::CutBasedIdLoose );
		//const bool isMedium = muon::isMediumMuon( muon );
        const bool isMedium = muon.passed( reco::Muon::CutBasedIdMedium );
		//const auto & primevtx = vertices_->front();
		//const bool isTight = muon::isTightMuon( muon, primevtx );
        //const bool isTight = muon.passed( reco::Muon::TkIsoTight );
		const bool isTight = muon.passed( reco::Muon::PFIsoVeryVeryTight );
        //const bool isTight = muon.passed( reco::Muon::CutBasedIdTight );
		//const bool isGlobal = muon.isGlobalMuon();

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
        Branches.fillBranch("tightID",isTight);
		////Branches.fillBranch("TrackZ",muTrackZ);
    	//Branches.fillBranch("timeValid","Muon_timeValid",VBOOL);

        if( MuonDEBUG ) std::cout << " --- Proccesssing : " << muon << std::endl;

        bool hasSVMatch = false;
        if( cfFlag("doSVModule") ) hasSVMatch = svObj->FoundLeptonMatch( muon );
        if( hasSVMatch ) nSVMatched++;
        Branches.fillBranch("svMatch",hasSVMatch);

        bool vmptc = ( muPt >= 20 );
        bool smptc = ( muPt >= 10 );
        if( smptc && isMedium ) nSelMu++;
        //Branches.fillBranch("hasMinPt",mptc);

        bool passPrompt = false;
        if( cfFlag("doSVModule") ) passPrompt = svObj->IsPromptLepton( muon );
        Branches.fillBranch("passPrompt",passPrompt);

        bool isPromptMuon = false;
        //if( passPrompt && isLoose && vmptc ){ isPromptMuon = true; nPromptMuon++; }
        if( isLoose && vmptc ){ isPromptMuon = true; nPromptMuon++; }
        Branches.fillBranch("IsPromptMuon",isPromptMuon);

	    muIndx++;

		//std::cout << " -- Muon getting sc : " << std::endl;
        const auto btptr = muon.bestTrack();
        auto ttrack = ttBuilder.build( *btptr );

        // GenParticle Info for muon  -------------------------------------------------------------------
        if( cfFlag("hasGenInfo") ){
            //scptrs.push_back(*scptr);
            //std::cout << " -- Muon sc : " << itptr->vx() << " - " << itptr->vy() << " - " << itptr->vz() << std::endl;
            //v3fPoint itv( itloc.x(), itloc.y(), itloc.z() );
			sctrks.push_back( ttrack );
			v3fPoint itv( muEta, muPhi, 0.f );
            scvertex.push_back( itv );
            scptres.push_back( muEnergy );
        }//<<>>if( hasGenInfo )


    }//<<>>for( const auto muon : *muons_ )
	Branches.fillBranch("nMu",muIndx);
    Branches.fillBranch("nSMu",nSelMu);
    Branches.fillBranch("nSVMatch",nSVMatched);
	Branches.fillBranch("nPrompt",nPromptMuon);
	geVar.set("nSelMu",nSelMu);

    int nGenMatched = 0;
    if( cfFlag("hasGenInfo") ){
        auto genInfo = genObjs->getGenMuonMatch( scvertex, scptres, sctrks );
        for( auto genidx : genInfo ){
            Branches.fillBranch("GenIdx",genidx);
			bool hasGenMatch = ( genidx > -1 );
			Branches.fillBranch("hasGenIdx", hasGenMatch );
			if( hasGenMatch ) nGenMatched++;
			int genWZMomIndx = hasGenMatch ? genObjs->getGenSigEleXMother( genidx ) : -5;
			Branches.fillBranch("GenWZIdx",genWZMomIndx);
			int genGMomIndx = ( genWZMomIndx > -1 ) ? genObjs->getGenSigEleXMother( genWZMomIndx ) : -5;
			Branches.fillBranch("GenXMomIdx",genGMomIndx);
			int momIdx = genObjs->getGenMomIdx( genidx ); 	
			int genMomPdgId = genObjs->getGenPdgId( momIdx ); 
			int gmomIdx = genObjs->getGenMomIdx( momIdx );	
			int genGMomPdgId = genObjs->getGenPdgId( gmomIdx );
			bool isLWZ =  ( genMomPdgId == 23 || genMomPdgId == 24 );		
			bool isLWZP = isLWZ && ( ( genGMomPdgId > 0 && genGMomPdgId < 10 ) || ( genGMomPdgId == 1000023 ) );	
			Branches.fillBranch("isLWZP", isLWZP );
        }//<<>>for( auto genidx : genInfo )
    }//<<>>if( cfFlag("hasGenInfo") )
	Branches.fillBranch("nGenMatch",nGenMatched);

}//<<>>void KUCMSMuon::ProcessEvent()

void KUCMSMuonObjectMini::EndJobs(){}

//float KUCMSMuonObjectMini::getMuTrackZMatch( const reco::Photon & photon ){
//
//	float bestmatch(0.1);
//	float elTrackZ(1000.0);
//	for (const auto & muon : fmuons){
//		auto match = reco::deltaR(muon,photon);
//		if ( match < bestmatch ){ elTrackZ = muon.trackPositionAtVtx().Z(); bestmatch = match; }
//	}//<<>>for (const auto muon : fmuons)
//	return elTrackZ;
//
//}//<<>>float KUCMSMuonObjectMini::getMuTrackZMatch( const reco::Photon & photon )

//bool KUCMSMuonObjectMini::getMuonVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ){
//
//	return conversions_.isValid() ?
//            not ConversionTools::hasMatchedPromptMuon( phosc, fmuons, *conversions_, beamSpot_->position() ) : false;

//}//<<>> bool KUCMSMuonObjectMini::getMuonVeto( const reco::SuperCluster & phosc )

int KUCMSMuonObjectMini::getIndex( float kideta, float kidphi ){

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

}//<<>>uInt KUCMSMuonObjectMini::getIndex( float kideta, float kidphi )

#endif
