// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Electron Object
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
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "CommonTools/Egamma/interface/EffectiveAreas.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSDisplacedVertex.hh"

#ifndef KUCMSElectronObjectHeader
#define KUCMSElectronObjectHeader

//#define ElectronDEBUG true
#define ElectronDEBUG false

using namespace edm; 

class KUCMSElectronObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSElectronObject( const edm::ParameterSet& iConfig );
    ~KUCMSElectronObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadElectronTokens( edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken ){ electronsToken_ = electronsToken; }; 
    void LoadConversionTokens( edm::EDGetTokenT<reco::ConversionCollection> conversionsToken ){ conversionsToken_ = conversionsToken; };
    void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    void LoadRhoTokens( edm::EDGetTokenT<double> RhoToken ){ RhoToken_ = RhoToken; };
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
    bool getElectronVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ); 
	float getEleTrackZMatch( const reco::Photon & photon ); 
    int getIndex( float kideta, float kidphi );

    private:

    //std::vector<examnple> fexamples;
    std::vector<reco::GsfElectron> felectrons;
    std::vector<int> eleIds;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    //const edm::InputTag electronsTag;
    //edm::EDGetTokenT<std::vector<reco::GsfElectron> > electronsToken_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken_;
    //edm::Handle<std::vector<reco::GsfElectron> > electrons_;
    edm::Handle<edm::View<reco::GsfElectron> > electrons_;
    //std::vector<reco::GsfElectron> electrons;

    //const edm::InputTag eleMVAIDLooseMapTag;
    //edm::EDGetTokenT<edm::ValueMap<bool>> eleMVAIDLooseMapToken_;
    //edm::Handle<edm::ValueMap<bool>> eleMVAIDLooseMap_;
	edm::EDGetTokenT<double> RhoToken_;
	edm::Handle<double>  rho_;

    // Conversions ( for electron veto ) ---------------------------------
    //const edm::InputTag conversionsTag;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::Handle<reco::ConversionCollection> conversions_;

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    edm::Handle<reco::BeamSpot> beamSpot_;

	// for cut based IDs
	EffectiveAreas effectiveAreas_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;
    KUCMSEcalRecHitObject* rhObj;
    KUCMSGenObject* genObjs;
	KUCMSDisplacedVertex* svObj;

};//<<>>class KUCMSElectron : public KUCMSObjectBase

KUCMSElectronObject::KUCMSElectronObject( const edm::ParameterSet& iConfig ) :
	effectiveAreas_((iConfig.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath())
{   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minEleE", iConfig.existsAs<double>("minEleE") ? iConfig.getParameter<double>("minEleE") : 2.0 );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );

}//<<>>KUCMSElectron::KUCMSElectron( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSElectronObject::InitObject( TTree* fOutTree ){

    //Branches.makeBranch("RhIds","Electron_rhIds",VVUINT);
    Branches.makeBranch("nEle","Electron_nElectrons",INT);
    Branches.makeBranch("scIndex","Electron_scIndex",VINT,"index of supercluster");
    Branches.makeBranch("SeedTOFTime","Electron_seedTOFTime",VFLOAT);
    Branches.makeBranch("Pt","Electron_pt",VFLOAT);
    Branches.makeBranch("Energy","Electron_energy",VFLOAT);
    Branches.makeBranch("Phi","Electron_phi",VFLOAT);
    Branches.makeBranch("Eta","Electron_eta",VFLOAT);
    Branches.makeBranch("Px","Electron_px",VFLOAT);
    Branches.makeBranch("Py","Electron_py",VFLOAT);
    Branches.makeBranch("Pz","Electron_pz",VFLOAT);
    Branches.makeBranch("TrackZ","Electron_trackz",VFLOAT);
    Branches.makeBranch("nSele","Electron_nSelElectrons",INT);
    Branches.makeBranch("GenIdx","Electron_genIdx",VINT);
    /////Branches.makeBranch("LepGenIdx","Electron_lepGenIdx",VINT);
    Branches.makeBranch("GenXMomIdx","Electron_genSigXMomId",VINT);
    Branches.makeBranch("GenWZIdx","Electron_genSigWZId",VINT);
    Branches.makeBranch("Sieie","Electron_Sieie",VFLOAT);
    Branches.makeBranch("DetaSCTV","Electron_DetaSCTV",VFLOAT);
    Branches.makeBranch("DphiSCTV","Electron_DphiSCTV",VFLOAT);
    Branches.makeBranch("HOE","Electron_HOE",VFLOAT);

    Branches.makeBranch("IsMedium","Electron_isMedium",VBOOL);
    Branches.makeBranch("IsLoose","Electron_isLoose",VBOOL);
    Branches.makeBranch("IsVeto","Electron_isVeto",VBOOL);
    Branches.makeBranch("svMatch","Electron_hasSVMatch",VBOOL);
    Branches.makeBranch("nSVMatch","Electron_nSVMatched",INT);
    Branches.makeBranch("IsPromptEle","Electron_isPrompt",VBOOL);
    Branches.makeBranch("passPrompt","Electron_passPrompt",VBOOL);
    Branches.makeBranch("nPrompt","Electron_nPrompt",INT);
    Branches.makeBranch("dRGenMatch","Electron_dRGenMatch",VFLOAT);
    Branches.makeBranch("hasGenIdx","Electron_hasGenMatch",VBOOL);
    Branches.makeBranch("nGenMatch","Electron_nGenMatch",INT);
    Branches.makeBranch("IsLVeto","Electron_isLVeto",VBOOL);
    Branches.makeBranch("Vsieiec","Electron_Vsieiec",VBOOL);
    Branches.makeBranch("Vdetaseed","Electron_Vdetaseed",VBOOL);
    Branches.makeBranch("Vdpscvc","Electron_Vdpscvc",VBOOL);
    Branches.makeBranch("Vhoec","Electron_Vhoec",VBOOL);
    Branches.makeBranch("VrelIsoWithEA","Electron_VrelIsoWithEA",VBOOL);
    Branches.makeBranch("Vooemoop","Electron_Vooemoop",VBOOL);


    //Branches.makeBranch("GenDr","Electron_genDr",VFLOAT);
    //Branches.makeBranch("GenDp","Electron_genDp",VFLOAT);
    //Branches.makeBranch("GenSIdx","Electron_genSIdx",VINT);
    //Branches.makeBranch("GenSDr","Electron_genSDr",VFLOAT);
    //Branches.makeBranch("GenSDr","Electron_genSDp",VFLOAT);
    //Branches.makeBranch("GenLlpId","Electron_genLlpId",VFLOAT);
    //Branches.makeBranch("GenSLlpId","Electron_genSLlpId",VFLOAT);
    //Branches.makeBranch("ninovx","Electron_genSigMomVx",VFLOAT);
    //Branches.makeBranch("ninovy","Electron_genSigMomVy",VFLOAT);
    //Branches.makeBranch("ninovz","Electron_genSigMomVz",VFLOAT);
    //Branches.makeBranch("ninopx","Electron_genSigMomPx",VFLOAT);
    //Branches.makeBranch("ninopy","Electron_genSigMomPy",VFLOAT);
    //Branches.makeBranch("ninopz","Electron_genSigMomPz",VFLOAT);
    //Branches.makeBranch("ninoeta","Electron_genSigMomEta",VFLOAT);
    //Branches.makeBranch("ninophi","Electron_genSigMomPhi",VFLOAT);
    //Branches.makeBranch("ninoe","Electron_genSigMomEnergy",VFLOAT);
    //Branches.makeBranch("ninom","Electron_genSigMomMass",VFLOAT);
    //Branches.makeBranch("ninopt","Electron_genSigMomPt",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSElectron::InitObject( TTree* fOutTree )

void KUCMSElectronObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){


	iEvent.getByToken(electronsToken_, electrons_);
    iEvent.getByToken(conversionsToken_, conversions_);
    iEvent.getByToken(beamLineToken_, beamSpot_);
    iEvent.getByToken(RhoToken_, rho_);

    if( ElectronDEBUG ) std::cout << "Collecting Electrons" << std::endl;

    felectrons.clear();
    eleIds.clear();// indexed by ele index ( 0,1,2 ) * number of ids ( 1 current, 6? possible ) + index of ID wanted
    for (edm::View<reco::GsfElectron>::const_iterator itElectron = electrons_->begin(); itElectron != electrons_->end(); itElectron++) {
        //auto idx = itElectron - electrons_->begin();//unsigned int
        //auto electronRef = electrons_->refAt(idx);//edm::RefToBase<reco::GsfElectron> 
        auto &electron = (*itElectron);
        //if ( cfFlag("onlyEB") && ( std::abs(electron.eta()) > cfPrm("ebMaxEta") ) ) continue;
        //auto passIdCut = true; //electron.electronID(eleCutLoose);// pat electron ( miniAOD ) method
        //eleIdBools.push_back((*eleMVAIDLooseMap_)[electronRef]);// order is important, track how this vector is loaded
        bool passEnergyCut = electron.energy() > cfPrm("minEleE");
        bool passEtaCut = true; //std::abs( energy.eta() ) < cfPrm("ebMaxEta");
        if( passEnergyCut && passEtaCut ) felectrons.push_back(electron);
    }//<<>>for( const auto electron : *electrons_ )


}//<<>>void KUCMSElectron::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSElectronObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSElectronObject::ProcessEvent( ItemManager<float>& geVar ){

    if( ElectronDEBUG ) std::cout << "Processing Electrons" << std::endl;

    Branches.clearBranches();

    if( ElectronDEBUG ) std::cout << " - enetering Electron loop" << std::endl;

	int eleIndx = 0;
	int nSelEle = 0;
	int nSVMatched = 0;
    int nPromptEle = 0;
    int nGenMatched = 0;
    scGroup scptrs;
    std::vector<float> scptres;
    if( ElectronDEBUG ) std::cout << "Processing Electrons" << std::endl;
    for( const auto &electron : felectrons ){

        const float elePt = electron.pt();
        const float eleEnergy = electron.energy();
        const float elePhi = electron.phi();
        const float eleEta = electron.eta();
        const float elePx = electron.px();
        const float elePy = electron.py();
        const float elePz = electron.pz();
        const float eleTrackZ = electron.trackPositionAtVtx().Z();

        Branches.fillBranch("Pt",elePt);
        Branches.fillBranch("Energy",eleEnergy);
        Branches.fillBranch("Phi",elePhi);
        Branches.fillBranch("Eta",eleEta);
        Branches.fillBranch("Px",elePx);
        Branches.fillBranch("Py",elePy);
        Branches.fillBranch("Pz",elePz);
		Branches.fillBranch("TrackZ",eleTrackZ);

        if( ElectronDEBUG ) std::cout << " --- Proccesssing : " << electron << std::endl;
        const auto &elesc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
        const auto scptr = elesc.get();
		//int encIndx = eleIndx+1;
        const auto scIndex = rhObj->getSuperClusterIndex(scptr,11,eleIndx);
        Branches.fillBranch("scIndex",scIndex);
        //const scGroup eleSCGroup{*scptr};
        //const auto eleRhGroup = rhObj->getRHGroup( eleSCGroup, 0.2 );
        //const auto eleRhIdsGroup = rhObj->getRhGrpIDs( eleRhGroup );
        //Branches.fillBranch("RhIds",eleRhIdsGroup );
        //rhObj->setRecHitUsed( eleRhIdsGroup );
        //if( ElectronDEBUG ) std::cout << " -- Electrons : " << scptr << " #: " << eleRhGroup.size() << std::endl;
        //const auto tofTimes = rhObj->getLeadTofRhTime( eleRhGroup, geVar("vtxX"), geVar("vtxY"), geVar("vtxZ") );
        //auto timeStats = getTimeDistStats( tofTimes, eleRhGroup );
        const float seedTOFTime = rhObj->getSeedTofTime( *scptr, geVar("vtxX"), geVar("vtxY"), geVar("vtxZ") );

		float dEtaSeed = std::abs( electron.deltaEtaSuperClusterTrackAtVtx() - scptr->eta() + scptr->seed()->eta() );
		float rho = rho_.isValid() ? (float)(*rho_) : 0;
		float epinv = std::abs( 1.0 - electron.eSuperClusterOverP() )/electron.ecalEnergy();
		constexpr auto missingHitType = reco::HitPattern::MISSING_INNER_HITS;
		auto eleMisTrack = electron.gsfTrack()->hitPattern().numberOfLostHits(missingHitType);
		const reco::GsfElectron::PflowIsolationVariables& pfIso = electron.pfIsolationVariables();
		const float chad = pfIso.sumChargedHadronPt;
		const float nhad = pfIso.sumNeutralHadronEt;
		const float pho = pfIso.sumPhotonEt;
		double absEta = std::abs(electron.superCluster()->eta());
		float  eA = effectiveAreas_.getEffectiveArea( absEta );
		float iso = (chad + std::max(0.0f, nhad + pho - rho*eA))/elePt;

		bool passConversionVeto = !ConversionTools::hasMatchedConversion( electron, *conversions_, beamSpot_->position());
        bool mistrack = ( eleMisTrack <= 1 );
		bool ebee = ( eleEta <= 1.479 );

        bool Lsieiec = ebee ? electron.full5x5_sigmaIetaIeta() < 0.0112 : electron.full5x5_sigmaIetaIeta() < 0.0425;
        bool Ldetaseed = ebee ? dEtaSeed < 0.00377 : dEtaSeed < 0.00674;
        bool Ldpscvc = ebee ? electron.deltaPhiSuperClusterTrackAtVtx() < 0.0884 : electron.deltaPhiSuperClusterTrackAtVtx() < 0.169;
        bool Lhoec = ebee ? electron.hadronicOverEm() < ( 0.05 + 1.16/scptr->energy() + 0.0324*rho/scptr->energy() ) 
						  : electron.hadronicOverEm() < ( 0.0441 + 2.54/scptr->energy() + 0.183*rho/scptr->energy() );
        bool Looemoop = ebee ? epinv < 0.193 : epinv < 0.111;
        bool LrelIsoWithEA = ebee ? iso < ( 0.112 + 0.506/elePt ) : iso < ( 0.108 + 0.963/elePt );

        bool Msieiec = ebee ? electron.full5x5_sigmaIetaIeta() < 0.0106 : electron.full5x5_sigmaIetaIeta() < 0.0387;
        bool Mdetaseed = ebee ? dEtaSeed < 0.0032 : dEtaSeed < 0.00632;
        bool Mdpscvc = ebee ? electron.deltaPhiSuperClusterTrackAtVtx() < 0.0547 : electron.deltaPhiSuperClusterTrackAtVtx() < 0.0394;
        bool Mhoec = ebee ? electron.hadronicOverEm() < ( 0.046 + 1.16/scptr->energy() + 0.0324*rho/scptr->energy() ) 
						  : electron.hadronicOverEm() < ( 0.0275 + 2.52/scptr->energy() + 0.183*rho/scptr->energy() );
        bool MrelIsoWithEA = ebee ? iso < ( 0.0478 + 0.506/elePt ) : iso < ( 0.0658 + 0.963/elePt );
        bool Mooemoop = ebee ? epinv < 0.184 : epinv < 0.0721;

		bool Vsieiec = ebee ? electron.full5x5_sigmaIetaIeta() < 0.0126 : electron.full5x5_sigmaIetaIeta() < 0.0457;
        bool Vdetaseed = ebee ? dEtaSeed < 0.00463 : dEtaSeed < 0.00814;
        bool Vdpscvc = ebee ? electron.deltaPhiSuperClusterTrackAtVtx() < 0.148 : electron.deltaPhiSuperClusterTrackAtVtx() < 0.19;
        bool Vhoec = ebee ? electron.hadronicOverEm() < ( 0.05 + 1.16/scptr->energy() + 0.0324*rho/scptr->energy() ) 
                          : electron.hadronicOverEm() < ( 0.05 + 2.54/scptr->energy() + 0.183*rho/scptr->energy() );
		bool VrelIsoWithEA = ebee ? iso < ( 0.198 + 0.506/elePt ) : iso < ( 0.203 + 0.963/elePt );
        bool Vooemoop = ebee ? epinv < 0.209 : epinv < 0.132; 
        bool Vmistrack = ebee ? ( eleMisTrack <= 2 ) : ( eleMisTrack <= 3 );

		Branches.fillBranch("Vsieiec",Vsieiec);
		Branches.fillBranch("Vdetaseed",Vdetaseed);
		Branches.fillBranch("Vdpscvc",Vdpscvc);
		Branches.fillBranch("Vhoec",Vhoec);
		Branches.fillBranch("VrelIsoWithEA",VrelIsoWithEA);
		Branches.fillBranch("Vooemoop",Vooemoop);

        bool septc = ( elePt >= 20 );
		bool veptc = ( elePt >= 10 );

        bool isMedium = Msieiec && Mdetaseed && Mdpscvc && Mhoec && Mooemoop && mistrack && MrelIsoWithEA && passConversionVeto;
		bool isLoose = Lsieiec && Ldetaseed && Ldpscvc && Lhoec && Looemoop && mistrack && LrelIsoWithEA && passConversionVeto;
		bool isVeto = Vsieiec && Vdetaseed && Vdpscvc && Vhoec && Vooemoop && Vmistrack && VrelIsoWithEA && passConversionVeto;
        bool isLVeto = Vdpscvc && Vooemoop && Vmistrack && passConversionVeto;
		if( septc && isLoose ) nSelEle++;

		Branches.fillBranch("IsLVeto",isLVeto);
        Branches.fillBranch("IsVeto",isVeto);
        Branches.fillBranch("IsLoose",isLoose);
        Branches.fillBranch("IsMedium",isMedium);

        Branches.fillBranch("SeedTOFTime",seedTOFTime);
    	Branches.fillBranch("Sieie",electron.full5x5_sigmaIetaIeta());
    	Branches.fillBranch("DetaSCTV",electron.deltaEtaSuperClusterTrackAtVtx());
    	Branches.fillBranch("DphiSCTV",electron.deltaPhiSuperClusterTrackAtVtx());
    	Branches.fillBranch("HOE",electron.hadronicOverEm());
		
		bool hasSVMatch = false;
		if( cfFlag("doSVModule") ) hasSVMatch = svObj->FoundLeptonMatch( electron );
		if( hasSVMatch ) nSVMatched++;
        Branches.fillBranch("svMatch",hasSVMatch);

        float dRGenMatch = -9;
        if( cfFlag("hasGenInfo") ) dRGenMatch = genObjs->getGenParticleIndex( electron, 11 );
        //Branches.fillBranch("LepGenIdx", genIndex );
        bool hasGenMatch = cfFlag("hasGenInfo")  ? ( dRGenMatch > -1 ) : false;
		if( hasGenMatch ) nGenMatched++;
		Branches.fillBranch("hasGenIdx", hasGenMatch );
    	Branches.fillBranch("dRGenMatch", dRGenMatch);

        bool passPrompt = false;
        if( cfFlag("doSVModule") ) passPrompt = svObj->IsPromptLepton( electron );
        Branches.fillBranch("passPrompt",passPrompt);

        bool isPromptEle = false;
        if( passPrompt && veptc && isLVeto && hasGenMatch ){ isPromptEle = true; nPromptEle++; }
        //if( isPrompt && isVeto && hasGenMatch ){ isPromptEle = true; nPromptEle++; }
        //if( isPrompt && isVeto && eptc ){ isPromptEle = true; nPromptEle++; }
        Branches.fillBranch("IsPromptEle",isPromptEle);


	    eleIndx++;
        // GenParticle Info for electron  -------------------------------------------------------------------
        if( cfFlag("hasGenInfo") ){
            scptrs.push_back(*scptr);
            scptres.push_back(eleEnergy);
        }//<<>>if( hasGenInfo )

    }//<<>>for( const auto electron : *electrons_ )
	Branches.fillBranch("nEle",eleIndx);
    Branches.fillBranch("nSele",nSelEle);
    Branches.fillBranch("nSVMatch",nSVMatched);
	Branches.fillBranch("nPrompt",nPromptEle);
    Branches.fillBranch("nGenMatch",nGenMatched);
	geVar.set("nSelEle",nSelEle);

    if( cfFlag("hasGenInfo") ){
        auto genInfo = genObjs->getGenEleMatch( scptrs, scptres );
        //if( ElectronDEBUG) 
		//std::cout << " Electron Gen Match ------------------------- " << std::endl;
        for( auto genidx : genInfo ){
			//std::cout << " -- Filling genindex : " << genidx << std::endl;
            Branches.fillBranch("GenIdx",genidx);
            if( genidx > -1 ){
                int genWZMomIndx = genObjs->getGenSigEleXMother( genidx );
                //std::cout << " -- Electron Match : level 1 : genidx: " << genidx << " mom: " << genMomIndx  << std::endl;
                Branches.fillBranch("GenWZIdx",genWZMomIndx);
				if( genWZMomIndx > -1 ){
					//uInt wzIndx = genMomIndx;
                	int genGMomIndx = genObjs->getGenSigEleXMother( genWZMomIndx );
                    //std::cout << "  -- Electron Match : level 2 : mom: " << genMomIndx << " gmom: " << genGMomIndx  << std::endl;
					Branches.fillBranch("GenXMomIdx",genGMomIndx);
				} else { Branches.fillBranch("GenXMomIdx",-5); }//<<>>if( genMomIndx > -1 )
            } else { Branches.fillBranch("GenWZIdx",-5); Branches.fillBranch("GenXMomIdx",-5); }//<<>>if( genidx > -1 )
			//std::cout << " --- next electron -------------------------- " << std::endl;
        }//<<>>for( auto genidx : genInfo )
		//std::cout << " Electron Gen Match Finished ------------------------- " << std::endl;
    }//<<>>if( cfFlag("hasGenInfo") )

}//<<>>void KUCMSElectron::ProcessEvent()

void KUCMSElectronObject::EndJobs(){}

float KUCMSElectronObject::getEleTrackZMatch( const reco::Photon & photon ){

	float bestmatch(0.1);
	float elTrackZ(1000.0);
	for (const auto & electron : felectrons){
		auto match = reco::deltaR(electron,photon);
		if ( match < bestmatch ){ elTrackZ = electron.trackPositionAtVtx().Z(); bestmatch = match; }
	}//<<>>for (const auto electron : felectrons)
	return elTrackZ;

}//<<>>float KUCMSElectronObject::getEleTrackZMatch( const reco::Photon & photon )

bool KUCMSElectronObject::getElectronVeto( const edm::Ref<std::vector<reco::SuperCluster>> phosc ){

	return conversions_.isValid() ?
            not ConversionTools::hasMatchedPromptElectron( phosc, felectrons, *conversions_, beamSpot_->position() ) : false;

}//<<>> bool KUCMSElectronObject::getElectronVeto( const reco::SuperCluster & phosc )

int KUCMSElectronObject::getIndex( float kideta, float kidphi ){

    auto bestdr = 0.01;
    int idx = -1;
    int cnt = 0;
    for( const auto &electron : felectrons  ){
        auto eleeta = electron.eta();
        auto elephi = electron.phi();
        auto kidpfdr = std::sqrt( reco::deltaR2( eleeta, elephi, kideta, kidphi ));
        if( kidpfdr < bestdr ){ bestdr = kidpfdr; idx = cnt; }
        cnt++;
    }//<<>>for( const auto electron : felectrons  )

	return idx;

}//<<>>uInt KUCMSElectronObject::getIndex( float kideta, float kidphi )

#endif
