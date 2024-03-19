// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Photon Object
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
//#include "DataFormats/Photon/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

// Supercluster info
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSElectron.hh"
#include "KUCMSGenObjects.hh"

#ifndef KUCMSPhotonObjectHeader
#define KUCMSPhotonObjectHeader

//#define PhotonDEBUG true
#define PhotonDEBUG false

using namespace edm; 

typedef edm::View<reco::Photon> phoColView;

class KUCMSPhotonObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSPhotonObject( const edm::ParameterSet& iConfig );
    ~KUCMSPhotonObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadPhotonTokens( edm::EDGetTokenT<phoColView> gedPhotonsToken, edm::EDGetTokenT<phoColView> ootPhotonsToken ); 
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObject*  genObj_ ){ genObj = genObj_; };
    void LoadElectronObject( KUCMSElectronObject* electronObj_ ){ electronObj = electronObj_; };

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
    int getIndex( float kideta, float kidphi );
    float getPhotonSeedTime( pat::Photon photon );
    float getPhotonSeedTime( reco::Photon photon );
	void correctMet( float & CSumEt, float & CPx, float & CPy );

    private:

    //std::vector<examnple> fexamples;
    std::vector<reco::Photon> fphotons;
    std::vector<bool> phoExcluded;
    std::vector<bool> phoIsOotPho;
    std::vector<int> phoIds;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;

    // gedPhotons
    //const edm::InputTag gedPhotonsTag;
    edm::EDGetTokenT<phoColView> gedPhotonsToken_;
    edm::Handle<phoColView> gedPhotons_;

    // ootPhotons
    //const edm::InputTag ootPhotonsTag;
    edm::EDGetTokenT<phoColView> ootPhotonsToken_;
    edm::Handle<phoColView> ootPhotons_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    KUCMSEcalRecHitObject* rhObj;
    KUCMSGenObject*  genObj;
    KUCMSElectronObject* electronObj;

};//<<>>class KUCMSPhoton : public KUCMSObjectBase

KUCMSPhotonObject::KUCMSPhotonObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfPrm.set( "minPhoE", iConfig.existsAs<double>("minPhoE") ? iConfig.getParameter<double>("minPhoE") : 2.0 );
    cfPrm.set( "phoMinPt", iConfig.existsAs<double>("phoMinPt") ? iConfig.getParameter<double>("phoMinPt") : 0.0 );
    cfPrm.set( "phoMinSeedTime", iConfig.existsAs<double>("phoMinSeedTime") ? iConfig.getParameter<double>("phoMinSeedTime") : 30.0 );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : false );

}//<<>>KUCMSPhoton::KUCMSPhoton( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSPhotonObject::LoadPhotonTokens( edm::EDGetTokenT<phoColView> gedPhotonsToken, edm::EDGetTokenT<phoColView> ootPhotonsToken ){

    gedPhotonsToken_ = gedPhotonsToken;
    ootPhotonsToken_ = ootPhotonsToken;

}//<<>>void LoadPhotonTokens( edm::EDGetTokenT<PhotonCol>* PhotonEBToken_, edm::EDGetTokenT<PhotonCol>* PhotonEEToken_ )

void KUCMSPhotonObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("IsOotPho","Photon_isOot",VBOOL);
    Branches.makeBranch("Excluded","Photon_excluded",VBOOL);
    //Branches.makeBranch("hasSC","Photon_hasSC",VBOOL);

    Branches.makeBranch("Pt","Photon_pt",VFLOAT);
    Branches.makeBranch("Energy","Photon_energy",VFLOAT);
    Branches.makeBranch("Phi","Photon_phi",VFLOAT);
    Branches.makeBranch("Eta","Photon_eta",VFLOAT);
    Branches.makeBranch("Px","Photon_px",VFLOAT);
    Branches.makeBranch("Py","Photon_py",VFLOAT);
    Branches.makeBranch("Pz","Photon_pz",VFLOAT);
    Branches.makeBranch("nPho","Photon_nPhoton",INT);

    Branches.makeBranch("EnergyErr","Photon_energyErr",VFLOAT,"energy error of the cluster from regression");//
    //Branches.makeBranch("EnergyRaw","Photon_energyRaw",VFLOAT,"raw energy of photon supercluster");//

    //Branches.makeBranch("JetIdx",VUINT);// index of matching jet/ele -> can do at skimmer ..
    //Branches.makeBranch("eleIdx",VUINT);// index of matching jet/ele -> can do at skimmer ..
    // superclusterEta  = Var("superCluster().eta()",float,doc="supercluster eta",precision=10),

    Branches.makeBranch("SeedTOFTime","Photon_seedTOFTime",VFLOAT,"time of flight from PV to photon seed crystal");
    //Branches.makeBranch("RhIds","Photon_rhIds",VVUINT,"list of rechit raw ids in hits and fractions list from supercluster");
    Branches.makeBranch("scIndex","Photon_scIndex",VINT,"index of supercluster");
    Branches.makeBranch("hasPixelSeed","Photon_pixelSeed",VBOOL,"has pixel seed");
    Branches.makeBranch("eleVeto","Photon_electronVeto",VBOOL,"pass electron veto");//
    //Branches.makeBranch("isEB","Photon_seedIsEB",VBOOL,"photon supercluster seed crystal is in ecal barrel");
    Branches.makeBranch("hasConversionTracks", VBOOL, "Variable specifying if photon has associated conversion tracks (one-legged or two-legged)");
	// hasConversionTracks = Var("hasConversionTracks()",bool,doc="Variable specifying if photon has associated conversion tracks (one-legged or two-legged)"

    Branches.makeBranch("R9","Photon_r9",VFLOAT,"R9 of the supercluster, calculated with full 5x5 region");
    Branches.makeBranch("SigmaIEtaIEta","Photon_SigmaIEtaIEta",VFLOAT);
    Branches.makeBranch("sieie","Photon_sieie",VFLOAT,"sigma_IetaIeta of supercluster, calc w/ full 5x5 region");
    Branches.makeBranch("sieip","Photon_sieip",VFLOAT,"sigmaIphiIphi of supercluster");//
    Branches.makeBranch("sipip","Photon_sipip",VFLOAT,"sigma_IetaIphi of supercluster, calc w/ full 5x5 region");//
    Branches.makeBranch("S4","Photon_s4",VFLOAT,"e2x2/e5x5 of supercluster, calc w/ full 5x5 region");//

    Branches.makeBranch("HadTowOverEM","Photon_hadTowOverEM",VFLOAT);
    Branches.makeBranch("HadOverEM","Photon_hadOverEM",VFLOAT);//
    Branches.makeBranch("EcalRHSumEtConeDR04","Photon_ecalRHSumEtConeDR04",VFLOAT);
    Branches.makeBranch("HcalTowerSumEtConeDR04","Photon_hcalTowerSumEtConeDR04",VFLOAT);
    Branches.makeBranch("HcalTowerSumEtBcConeDR04","Photon_hcalTowerSumEtBcConeDR04",VFLOAT);
    Branches.makeBranch("TrkSumPtSolidConeDR04","Photon_trkSumPtSolidConeDR04",VFLOAT);
    Branches.makeBranch("TrkSumPtHollowConeDR04","Photon_trkSumPtHollowConeDR04",VFLOAT);
    Branches.makeBranch("TrkSumPtHollowConeDR03","Photon_trkSumPtHollowConeDR03",VFLOAT,"Sum of track pT in a hollow cone of outer radius, inner radius");// nano -> DR03?
    Branches.makeBranch("NTrkSolidConeDR04","Photon_nTrkSolidConeDR04",VFLOAT);
    Branches.makeBranch("NTrkHollowConeDR04","Photon_nTrkHollowConeDR04",VFLOAT);

    Branches.makeBranch("ecalPFClusterIso","Photon_ecalPFClusterIso",VFLOAT);//
    Branches.makeBranch("hcalPFClusterIso","Photon_hcalPFClusterIso",VFLOAT);//
	// ecalPFClusterIso = Var("ecalPFClusterIso()",float,doc="sum pt of ecal clusters, vetoing clusters part of photon", precision=8),
	// hcalPFClusterIso = Var("hcalPFClusterIso()",float,doc="sum pt of hcal clusters, vetoing clusters part of photon", precision=8), 
	// two above are : //backwards compat functions for pat::Photon not filled in AOD

	Branches.makeBranch("pfChargedIso", "Photon_pfChargedIso", VFLOAT, "PF absolute isolation dR=0.3, charged component with dxy,dz match to PV");
	// pfChargedIso = Var("chargedHadronIso()",float,doc="PF absolute isolation dR=0.3, charged component with dxy,dz match to PV", precision=8),  
    Branches.makeBranch("pfPhoIso03","Photon_pfPhoIso03",VFLOAT,"PF abs iso dR=0.3, photon component (uncorrected)");//
    Branches.makeBranch("pfChargedIsoPFPV","Photon_pfChargedIsoPFPV",VFLOAT,"PF abs iso dR=0.3, charged component (PF PV only)");//
    Branches.makeBranch("pfChargedIsoWorstVtx","Photon_pfChargedIsoWorstVtx",VFLOAT,"PF abs iso dR=0.3, charged component (Vertex w/ largest iso)");//
    Branches.makeBranch("pfRelIso03_chg_quadratic",VFLOAT);//
    Branches.makeBranch("pfRelIso03_all_quadratic",VFLOAT);//
    Branches.makeBranch("hoe_PUcorr","Photon_Hoe_PUcorr",VFLOAT, "PU crt H/E (cone-based w/ quadraticEA*rho*rho + linearEA*rho Winter22V1 corts)");//UF

    //Branches.makeBranch("isScEtaEB","Photon_isScEtaEB",VBOOL,"is supercluster eta within barrel acceptance");//
    //Branches.makeBranch("isScEtaEE","Photon_isScEtaEE",VBOOL,"is supercluster eta within endcap acceptance");//

    //Branches.makeBranch("seedGain",VINT);// ? can find through rh collection if needed
    //Branches.makeBranch("seediEtaOriX","Photon_seediEtaOriX",VINT,"iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100.");//
    //Branches.makeBranch("seediPhiOriY","Photon_seediPhiOriY",VINT,"iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100.");//
    //Branches.makeBranch("x_calo","Photon_x_calo",VFLOAT,"photon supercluster position on calorimeter, x coordinate (cm)");//
    //Branches.makeBranch("y_calo","Photon_y_calo",VFLOAT,"photon supercluster position on calorimeter, y coordinate (cm)");//
    //Branches.makeBranch("z_calo","Photon_z_calo",VFLOAT,"photon supercluster position on calorimeter, z coordinate (cm)");//

    //Branches.makeBranch("esEffSigmaRR","Photon_esEffSigmaRR",VFLOAT,"preshower sigmaRR");//
    //Branches.makeBranch("esEnergyOverRawE","Photon_esEnergyOverRawE",VFLOAT,"ratio of preshower energy to raw supercluster energy");
    //Branches.makeBranch("haloTaggerMVAVal","Photon_haloTaggerMVAVal",VFLOAT,"Value of MVA based beam halo tagger in the Ecal endcap (valid for pT > 200 GeV)");//

    Branches.makeBranch("GenIdx","Photon_genIdx",VINT);
    //Branches.makeBranch("GenDr","Photon_genDr",VFLOAT);
    //Branches.makeBranch("GenDp","Photon_genDp",VFLOAT);
    //Branches.makeBranch("GenSIdx","Photon_genSIdx",VINT);
    //Branches.makeBranch("GenSDr","Photon_genSDr",VFLOAT);
    //Branches.makeBranch("GenSDp","Photon_genSDp",VFLOAT);
    //Branches.makeBranch("GenLlpId","Photon_genLlpId",VFLOAT);
    //Branches.makeBranch("GenSLlpId","Photon_genSLlpId",VFLOAT);
    Branches.makeBranch("GenMomIdx","Photon_genSigMomId",VINT);
    //Branches.makeBranch("ninovx","Photon_genSigMomVx",VFLOAT);
    //Branches.makeBranch("ninovy","Photon_genSigMomVy",VFLOAT);
    //Branches.makeBranch("ninovz","Photon_genSigMomVz",VFLOAT);
    //Branches.makeBranch("ninopx","Photon_genSigMomPx",VFLOAT);
    //Branches.makeBranch("ninopy","Photon_genSigMomPy",VFLOAT);
    //Branches.makeBranch("ninopz","Photon_genSigMomPz",VFLOAT);
    //Branches.makeBranch("ninoeta","Photon_genSigMomEta",VFLOAT);
    //Branches.makeBranch("ninophi","Photon_genSigMomPhi",VFLOAT);
    //Branches.makeBranch("ninoe","Photon_genSigMomEnergy",VFLOAT);
    //Branches.makeBranch("ninom","Photon_genSigMomMass",VFLOAT);
    //Branches.makeBranch("ninopt","Photon_genSigMomPt",VFLOAT);


    //Branches.makeBranch("etaWidth","Photon_etaWidth",VFLOAT,"Width of the photon supercluster in eta");//
    //Branches.makeBranch("phiWidth","Photon_phiWidth",VFLOAT,"Width of the photon supercluster in phi");//
    //Branches.makeBranch("SMaj","Photon_smaj",VFLOAT);
    //Branches.makeBranch("SMin","Photon_smin",VFLOAT);
    //Branches.makeBranch("SAlp","Photon_salp",VFLOAT);
    //Branches.makeBranch("CovEtaEta","Photon_covEtaEta",VFLOAT);
    //Branches.makeBranch("CovEtaPhi","Photon_covEtaPhi",VFLOAT);
    //Branches.makeBranch("CovPhiPhi","Photon_covPhiPhi",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSPhoton::InitObject( TTree* fOutTree )

void KUCMSPhotonObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // PHOTONS
    iEvent.getByToken(gedPhotonsToken_, gedPhotons_);
    //iEvent.getByToken(phoCBIDLooseMapToken_, phoCBIDLooseMap_);
    iEvent.getByToken(ootPhotonsToken_, ootPhotons_);

    fphotons.clear();
    phoExcluded.clear();
    phoIsOotPho.clear();
	phoIds.clear();

    if( PhotonDEBUG ) std::cout << "Collecting Photons/OOTPhotons" << std::endl;
    std::vector<reco::Photon> fphotons_temp;
    std::vector<bool> phoExcluded_temp;
    std::vector<bool> phoIsOotPho_temp;
    //std::vector<int> phoIds_temp;// indexed by pho index ( 0,1,2 ) * number of ids ( 1 current, 6? possible ) + index of ID wanted
    for (edm::View<reco::Photon>::const_iterator itPhoton = ootPhotons_->begin(); itPhoton != ootPhotons_->end(); itPhoton++) {
        //auto idx = itPhoton - ootPhotons_->begin();//unsigned int
        //auto ootPhoRef = ootPhotons_->refAt(idx);//edm::RefToBase<reco::GsfElectron> 
        auto &ootPho = (*itPhoton);
        //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
        auto minPhoPt = ootPho.pt() < 0.0; //cfPrm("phoMinPt");
        auto phoSeedTime = getPhotonSeedTime(ootPho);
        auto minTime = std::abs(phoSeedTime) > cfPrm("phoMinSeedTime");
        auto minEnergy = ootPho.energy() < cfPrm("minPhoE");
        if( minPhoPt || minTime || minEnergy ) continue;
        fphotons_temp.push_back(ootPho);
        phoIsOotPho_temp.push_back(true);
        //phoIdBools.push_back((*phoCBIDLooseMap_)[ootPhoRef]);// not implimented 
/*
		// Set photon excluded via dr and pt critera
		//
        double minDr(1.0);
        double dRmatch(10.0);
        float matchpt(0);
        auto oEta = ootPho.eta();
        auto oPhi = ootPho.phi();
        auto oPt = ootPho.pt();
        for( const auto &gedPho : *gedPhotons_ ){
            //if( cfFlag("onlyEB") && gedPho.isEE() ) continue;
            auto pEta = gedPho.eta();
            auto pPhi = gedPho.phi();
            auto pPt = gedPho.pt();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
			Branches.fillBranch("zscnOExDr",dRmatch);
			Branches.fillBranch("zscnOExDr",dRmatch);
			Branches.fillBranch("zscnOExDr",dRmatch);
			Branches.fillBranch("zscnOExDr",dRmatch);
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = pPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( dRmatch < 0.3 && oPt < matchpt ) phoExcluded_temp.push_back(true);
        else phoExcluded_temp.push_back(false);
		////////////////////////////////////
*/
		// Set photon excluded false always
		//
		phoExcluded_temp.push_back(false);
		///////////////////////////////////
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )
    for (edm::View<reco::Photon>::const_iterator itPhoton = gedPhotons_->begin(); itPhoton != gedPhotons_->end(); itPhoton++) {
        //auto idx = itPhoton - gedPhotons_->begin();//unsigned int
        //auto gedPhoRef = gedPhotons_->refAt(idx);//edm::RefToBase<reco::GsfElectron> 
        auto &gedPho = (*itPhoton);
        if( cfFlag("onlyEB") && gedPho.isEE() ) continue;
        auto minPt = gedPho.pt() < 0; //cfPrm("phoMinPt");
        auto phoSeedTime = getPhotonSeedTime(gedPho);
        auto minTime = std::abs(phoSeedTime) > cfPrm("phoMinSeedTime");
        auto minEnergy = gedPho.energy() < cfPrm("minPhoE");
        if( minPt || minTime || minEnergy ) continue;
        fphotons_temp.push_back(gedPho);
        phoIsOotPho_temp.push_back(false);
        //phoIdBools.push_back((*phoCBIDLooseMap_)[gedPhoRef]);
/*
		// Set photon excluded via dr and pt critera
		//
        double minDr(1.0);
        double dRmatch(10.0);
        float matchpt(0);
        auto pEta = gedPho.eta();
        auto pPhi = gedPho.phi();
        auto pPt = gedPho.pt();
        for( const auto &ootPho : *ootPhotons_ ){
            //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
            auto oEta = ootPho.eta();
            auto oPhi = ootPho.phi();
            auto oPt = ootPho.pt();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = oPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( dRmatch < 0.3 && pPt < matchpt ) phoExcluded_temp.push_back(true);
        else phoExcluded_temp.push_back(false);
		/////////////////////////////////////////
*/
        // Set photon excluded if dr match to oot photon
        //
        double minDr(1.0);
        double dRmatch(10.0);
        //float matchpt(0);
        auto pEta = gedPho.eta();
        auto pPhi = gedPho.phi();
        for( const auto &ootPho : *ootPhotons_ ){
            //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
            auto oEta = ootPho.eta();
            auto oPhi = ootPho.phi();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
            if( dRmatch < minDr ){ minDr = dRmatch; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( dRmatch < 0.3 ) phoExcluded_temp.push_back(true);
        else phoExcluded_temp.push_back(false);
		///////////////////////////////////////

    }//<<>>for( int io = 0; io < nOotPhotons; io++ )

	std::map< float,int > phoOrderIndx;
	int it = 0;
	for( auto pho : fphotons_temp ){ 
		float ordpt = pho.pt();
		float iordpt = ordpt;
		float nordpt = ordpt - 0.00002;
		while( phoOrderIndx.count(ordpt) > 0 ){ nordpt = ( iordpt + nordpt ) / 2.0; ordpt = nordpt; }
		phoOrderIndx[ordpt] = it;
		it++;
	}//<<>>for( auto pho : fphotons_temp ){

	for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ ){

		fphotons.push_back(fphotons_temp[phoptit->second]);
        phoExcluded.push_back(phoExcluded_temp[phoptit->second]);
        phoIsOotPho.push_back(phoIsOotPho_temp[phoptit->second]);

	}//<<>>for( auto phoptit = phoOrderIndx.crbegin(); phoptit != phoOrderIndx.crend(); phoptit++ )


}//<<>>void KUCMSPhoton::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSPhotonObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSPhotonObject::ProcessEvent( ItemManager<float>& geVar ){

    if( PhotonDEBUG ) std::cout << "Processing Photons" << std::endl;

    Branches.clearBranches();

    if( PhotonDEBUG ) std::cout << " - enetering Photon loop" << std::endl;

    int phoIdx = 0;
	scGroup scptrs;
	std::vector<float> scptres;
	//std::vector<int> scmatched;
    for( const auto &photon : fphotons ){

        Branches.fillBranch("IsOotPho",phoIsOotPho[phoIdx]);
        Branches.fillBranch("Excluded",phoExcluded[phoIdx]);

        const float phoPt = photon.pt();
		if( PhotonDEBUG ) std::cout << " -- Pho# " << phoIdx << " Pt: " << phoPt << " oot: " << phoIsOotPho[phoIdx]  << std::endl;
        const float phoEnergy = photon.energy();
        const float phoPhi = photon.phi();
        const float phoEta = photon.eta();
        const float phoPx = photon.px();
        const float phoPy = photon.py();
        const float phoPz = photon.pz();

        const float phoEnergyErr = photon.getCorrectedEnergyError(reco::Photon::regression2);
        //const float haloTaggerMVAVal = photon.haloTaggerMVAVal();
        const bool phoHasPixelSeed = photon.hasPixelSeed();
        const bool hasConversionTracks = photon.hasConversionTracks();

        const float phoHadOverEM = photon.hadronicOverEm();
        //const float phoHadOverEMVaid = photon.hadronicOverEmValid();
        const float phoHadTowOverEM = photon.hadTowOverEm();
        //const float phoHadTowOverEMValid = photon.hadTowOverEmValid();
        const float hoe_PUcorr = 0;//UF needs PhoHoverEValueMapProducer slimmedPhotons, hasConversionTracks ( not in AOD ), 
        //const float phoMaxEnergyXtal = photon.maxEnergyXtal();
        const float phoSigmaIEtaIEta = photon.sigmaIetaIeta();
        const float sieie = photon.showerShapeVariables().sigmaIetaIeta;
        const float sieip = photon.showerShapeVariables().sigmaIetaIphi;
        const float sipip = photon.showerShapeVariables().sigmaIphiIphi;
        const float s4 = photon.full5x5_showerShapeVariables().e2x2/photon.full5x5_showerShapeVariables().e5x5;
    	//const float smaj = photon.showerShapeVariables().smMajor;
        //const float smin = photon.showerShapeVariables().smMinor;
        //const float salpha = photon.showerShapeVariables().smAlpha;
    	//const float esEffSigmaRR = photon.showerShapeVariables().effSigmaRR;

        //const float phoR1x5 = photon.r1x5();
        //const float phoR2x5 = photon.r2x5();
        const float phoR9 = photon.r9();

        const float phoEcalRHSumEtConeDR04 = photon.ecalRecHitSumEtConeDR04();
        const float phoHcalTwrSumEtConeDR04 = photon.hcalTowerSumEtConeDR04();
        const float phoHcalTowerSumEtBcConeDR04 = photon.hcalTowerSumEtBcConeDR04();
        const float phoTrkSumPtSolidConeDR04 = photon.trkSumPtSolidConeDR04();
        const float phoTrkSumPtHollowConeDR04 = photon.trkSumPtHollowConeDR04();
        const float phoNTrkSolidConeDR04 = photon.nTrkSolidConeDR04();
        const float phoNTrkHollowConeDR04 = photon.nTrkHollowConeDR04();
        const float phoTrkSumPtHollowConeDR03 = photon.trkSumPtHollowConeDR03();

        const float pfPhoIso = photon.chargedHadronIso();
        const float pfPhoIso03 = photon.photonIso();
        const float pfChargedIsoPFPV = photon.chargedHadronPFPVIso();
        const float pfChargedIsoWorstVtx = photon.chargedHadronWorstVtxIso();
        const float pfRelIso03_chg_quadratic = 0;//UF needs  isoForPho producer : effAreaPhotons_cone03_pfChargedHadrons_90percentBased_V2.txt
        const float pfRelIso03_all_quadratic = 0;//UF needs  isoForPho producer : w/ _pfNeutralHadrons_ && _pfPhotons_ variates for both
        const float ecalPFClusterIso = photon.ecalPFClusterIso();
        const float hcalPFClusterIso = photon.hcalPFClusterIso();
        //const float hoe_PUcorr = 0;//UF needs PhoHoverEValueMapProducer slimmedPhotons, hasConversionTracks ( not in AOD ),  QuadraticEAFile_HoverE 

        Branches.fillBranch("Pt",phoPt);
        Branches.fillBranch("Energy",phoEnergy);
        Branches.fillBranch("Phi",phoPhi);
        Branches.fillBranch("Eta",phoEta);
        Branches.fillBranch("Px",phoPx);
        Branches.fillBranch("Py",phoPy);
        Branches.fillBranch("Pz",phoPz);

        Branches.fillBranch("S4",s4);
        //Branches.fillBranch("esEffSigmaRR",esEffSigmaRR);
        Branches.fillBranch("SigmaIEtaIEta",phoSigmaIEtaIEta);
        Branches.fillBranch("sieie",sieie);
        Branches.fillBranch("sieip",sieip);
        Branches.fillBranch("sipip",sipip);
        //Branches.fillBranch("SMaj",smaj);
        //Branches.fillBranch("SMin",smin);
        //Branches.fillBranch("SAlp",salpha);

        Branches.fillBranch("EnergyErr",phoEnergyErr);
        //Branches.fillBranch("haloTaggerMVAVal",haloTaggerMVAVal);
        Branches.fillBranch("hasPixelSeed",phoHasPixelSeed);
        Branches.fillBranch("hasConversionTracks",hasConversionTracks);

        Branches.fillBranch("HadOverEM",phoHadOverEM);
        //Branches.fillBranch("HadOverEMVaid",phoHadOverEmValid);
        Branches.fillBranch("HadTowOverEM",phoHadTowOverEM);
        //Branches.fillBranch("hadTowOverEMValid",phoHadTowOverEmValid);
        //Branches.fillBranch("MaxEnergyXtal",phoMaxEnergyXtal);
        Branches.fillBranch("hoe_PUcorr",hoe_PUcorr);

        //Branches.fillBranch("R1x5",phoR1x5);
        //Branches.fillBranch("R2x5",phoR2x5);
        Branches.fillBranch("R9",phoR9);

        Branches.fillBranch("EcalRHSumEtConeDR04",phoEcalRHSumEtConeDR04);
        Branches.fillBranch("HcalTowerSumEtConeDR04",phoHcalTwrSumEtConeDR04);
        Branches.fillBranch("HcalTowerSumEtBcConeDR04",phoHcalTowerSumEtBcConeDR04);
        Branches.fillBranch("TrkSumPtSolidConeDR04",phoTrkSumPtSolidConeDR04);
        Branches.fillBranch("TrkSumPtHollowConeDR04",phoTrkSumPtHollowConeDR04);
        Branches.fillBranch("TrkSumPtHollowConeDR03",phoTrkSumPtHollowConeDR03);
        Branches.fillBranch("NTrkSolidConeDR04",phoNTrkSolidConeDR04);
        Branches.fillBranch("NTrkHollowConeDR04",phoNTrkHollowConeDR04);

        Branches.fillBranch("ecalPFClusterIso",ecalPFClusterIso);
        Branches.fillBranch("hcalPFClusterIso",hcalPFClusterIso);

        Branches.fillBranch("pfChargedIso",pfPhoIso);
        Branches.fillBranch("pfPhoIso03",pfPhoIso03);
        Branches.fillBranch("pfChargedIsoPFPV",pfChargedIsoPFPV);
        Branches.fillBranch("pfChargedIsoWorstVtx",pfChargedIsoWorstVtx);
        Branches.fillBranch("pfRelIso03_chg_quadratic",pfRelIso03_chg_quadratic);
        Branches.fillBranch("pfRelIso03_all_quadratic",pfRelIso03_all_quadratic);


        if( PhotonDEBUG ) std::cout << " --- Proccesssing : " << photon << std::endl;
        const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
        const auto scptr = phosc.get();
/*
        const float phoEnergyRaw = scptr->rawEnergy();
        const bool isScEtaEB = abs(scptr->eta()) < 1.4442;
        const bool isScEtaEE = abs(scptr->eta()) > 1.566 && abs(scptr->eta()) < 2.5;
        const int seediEtaOriX = scptr->seedCrysIEtaOrIx();
        const int seediPhiOriY = scptr->seedCrysIPhiOrIy();
        const float x_calo = scptr->seed()->position().x();
        const float y_calo = scptr->seed()->position().y();
        const float z_calo = scptr->seed()->position().z();
        //const float esEnergyOverRawE = scptr->preshowerEnergy()/phoEnergyRaw;
        const float etaWidth = scptr->etaWidth();
        const float phiWidth = scptr->phiWidth();

        if( PhotonDEBUG ) std::cout << " --- Finding Moments & Covariences : " << scptr << std::endl;
    	const auto &seedDetId = scptr->seed()->seed();// seed detid
    	const auto isEB = (seedDetId.subdetId() == EcalBarrel);// which subdet
        const auto ph2ndMoments = rhObj->getCluster2ndMoments( scptr );
        const auto lCov = rhObj->getCovariances( scptr );

        if( PhotonDEBUG ) std::cout << " --- Storing Moments & Covariences : " << scptr << std::endl;
        const float phoSMaj = ph2ndMoments.sMaj;
        const float phoSMin = ph2ndMoments.sMin;
        const float phoSAlp = ph2ndMoments.alpha;
        const float phoCovEtaEta = lCov[0];
        const float phoCovEtaPhi = lCov[1];
        const float phoCovPhiPhi = lCov[2];

        Branches.fillBranch("isEB",isEB);
        Branches.fillBranch("EnergyRaw",phoEnergyRaw);
        Branches.fillBranch("isScEtaEB",isScEtaEB);
        Branches.fillBranch("isScEtaEE",isScEtaEE);
        Branches.fillBranch("seediEtaOriX",seediEtaOriX);
        Branches.fillBranch("seediPhiOriY",seediPhiOriY);
        Branches.fillBranch("x_calo",x_calo);
        Branches.fillBranch("y_calo",y_calo);
        Branches.fillBranch("z_calo",z_calo);
        //Branches.fillBranch("esEnergyOverRawE",esEnergyOverRawE);
        Branches.fillBranch("etaWidth",etaWidth);
        Branches.fillBranch("phiWidth",phiWidth);

        Branches.fillBranch("SMaj",phoSMaj);
        Branches.fillBranch("SMin",phoSMin);
        Branches.fillBranch("SAlp",phoSAlp);
        Branches.fillBranch("CovEtaEta",phoCovEtaEta);
        Branches.fillBranch("CovEtaPhi",phoCovEtaPhi);
        Branches.fillBranch("CovPhiPhi",phoCovPhiPhi);
*/

		auto passelectronveto = electronObj->getElectronVeto( phosc ); 
        Branches.fillBranch("eleVeto",passelectronveto);

        if( PhotonDEBUG ) std::cout << " --- Gathering SC info : " << scptr << std::endl;
		//int encIndx = int(phoIdx)*100;
		const auto scIndex = rhObj->getSuperClusterIndex(scptr,22,phoIdx);
		Branches.fillBranch("scIndex",scIndex);
		//scmatched.push_back(scIndex);
		//if( scIndex == -1 ) std::cout << " - phoSC: " << scIndex << " - " << phoEnergyRaw << " - " << " - " << seedDetId.rawId() << std::endl;   
        //const scGroup phoSCGroup{*scptr};
        //const auto phoRhGroup = rhObj->getRHGroup( phoSCGroup, 0.2 );
        //const auto phoRhIdsGroup = rhObj->getRhGrpIDs( phoRhGroup );
        //Branches.fillBranch("RhIds",phoRhIdsGroup);
        //rhObj->setRecHitUsed(phoRhIdsGroup);
        //if( PhotonDEBUG ) std::cout << " -- gedPhotons : " << scptr << " #: " << phoRhGroup.size() << std::endl;
        //auto tofTimes = rhObj->getLeadTofRhTime( phoRhGroup, geVar("vtxX"), geVar("vtxY"), geVar("vtxZ") );
        //auto timeStats = getTimeDistStats( tofTimes, phoRhGroup );
        const auto seedTOFTime = rhObj->getSeedTofTime( *scptr, geVar("vtxX"), geVar("vtxY"), geVar("vtxZ") );
        //auto phoLeadTOFTime =  rhObj->getLeadTofTime( phoRhGroup, geVar("vtxX"), geVar("vtxY"), geVar("vtxZ") );
        Branches.fillBranch("SeedTOFTime",seedTOFTime);
        //Branches.fillBranch("CMeanTime",timeStats[6]);

        // GenParticle Info for photon  -------------------------------------------------------------------
        if( cfFlag("hasGenInfo") ){
			scptrs.push_back(*scptr);
			scptres.push_back(phoEnergy);
        }//<<>>if( hasGenInfo )

        phoIdx++;
    }//<<>>for( const auto &photon : fPhotons 
	Branches.fillBranch("nPho",phoIdx);

	if( cfFlag("hasGenInfo") ){
		auto genInfo = genObj->getGenPhoMatch( scptrs, scptres );
		//std::cout << " - genmatched - scmatched comp : " << std::endl;
		//std::cout << " Photon Gen Match ------------------------- " << std::endl;
		//for( uInt i = 0; i < genInfo.size(); i++ ){ std::cout << " -- gen: " << genInfo[i] << " sc: " << scmatched[i] << std::endl; }
		//if( PhotonDEBUG) std::cout << " Photon Match ------------------------- " << std::endl;
		for( auto genidx : genInfo ){
			//std::cout << " -- Filling genindex : " << genidx << std::endl;  
			Branches.fillBranch("GenIdx",genidx);
			if( genidx > -1 ){ 
			 	int genMomIndx = genObj->getGenSigPhoInfo( genidx );
				//std::cout << " -- Photon Match : level 1 : genidx: " << genidx << " mom: " << genMomIndx  << std::endl;		
				Branches.fillBranch("GenMomIdx",genMomIndx);
			} else { Branches.fillBranch("GenMomIdx",-5); }//<<>>if( genidx > -1 )
			//std::cout << " --- next photon -------------------------- " << std::endl;
		}//<<>>for( auto genidx : genInfo )
		//std::cout << " Photon Gen Match Finished ------------------------- " << std::endl;
	}//<<>>if( cfFlag("hasGenInfo") )

}//<<>>void KUCMSPhoton::ProcessEvent()

void KUCMSPhotonObject::EndJobs(){}

float KUCMSPhotonObject::getPhotonSeedTime( pat::Photon photon ){

    const auto & phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
    return rhObj->getSuperClusterSeedTime( phosc );

}//<<>>float KUCMSPhotonObject::getPhotonSeedTime( pat::Photon photon )

float KUCMSPhotonObject::getPhotonSeedTime( reco::Photon photon ){

    const auto & phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
    return rhObj->getSuperClusterSeedTime( phosc );

}//<<>>float KUCMSPhotonObject::getPhotonSeedTime( reco::Photon photon )

int KUCMSPhotonObject::getIndex( float kideta, float kidphi ){

    auto bestdr = 0.01;
    int idx = -1;
    int cnt = 0;
    for( const auto &photon : fphotons  ){
        auto phoeta = photon.eta();
        auto phophi = photon.phi();
        auto kidpfdr = std::sqrt( reco::deltaR2( phoeta, phophi, kideta, kidphi ));
        if( kidpfdr < bestdr ){ bestdr = kidpfdr; idx = cnt; }
        cnt++;
    }//<<>>for( const auto photon : fphotons  )

	return idx;

}//<<>>uInt KUCMSPhotonObject::getIndex( float kideta, float kidphi )

void KUCMSPhotonObject::correctMet( float & CSumEt, float & CPx, float & CPy ){

    int nphidx = fphotons.size();
    for( auto phidx = 0; phidx < nphidx; phidx++ ){

        auto photon = fphotons[phidx];
        if( phoIsOotPho[phidx] && not phoExcluded[phidx] ){
            CPx -= ((photon.pt())*std::cos(photon.phi()));
            CPy -= ((photon.pt())*std::sin(photon.phi()));
            CSumEt -= photon.et();
        }//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )
        if( not phoIsOotPho[phidx] && phoExcluded[phidx] ){
            CPx += ((photon.pt())*std::cos(photon.phi()));
            CPy += ((photon.pt())*std::sin(photon.phi()));
            CSumEt += photon.et();
        }//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )   

    }//<<>>for( auto phidx = 0; phidx < fphotons.size(); phidx++ )

}//<<>>void KUCMSPhotonObject::correctMet( float & CSumEt, float & CPx, float & CPy )

#endif
>>>>>>> 2f0e277fb0524b1d969b7ed35559585a2c17f0ba
