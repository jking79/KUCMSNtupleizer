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

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"

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
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObject* genObjs_ ){ genObjs = genObjs_; };

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

    // Conversions ( for electron veto ) ---------------------------------
    //const edm::InputTag conversionsTag;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::Handle<reco::ConversionCollection> conversions_;

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    edm::Handle<reco::BeamSpot> beamSpot_;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;
    KUCMSEcalRecHitObject* rhObj;
    KUCMSGenObject* genObjs;

};//<<>>class KUCMSElectron : public KUCMSObjectBase

KUCMSElectronObject::KUCMSElectronObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "minEleE", iConfig.existsAs<double>("minEleE") ? iConfig.getParameter<double>("minEleE") : 2.0 );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );

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
    Branches.makeBranch("GenIdx","Electron_genIdx",VINT);
    Branches.makeBranch("GenMomIdx","Electron_genSigMomId",VINT);
    Branches.makeBranch("GenWZIdx","Electron_genSigWZId",VINT);
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

    if( ElectronDEBUG ) std::cout << "Collecting Electrons" << std::endl;

    felectrons.clear();
    eleIds.clear();// indexed by ele index ( 0,1,2 ) * number of ids ( 1 current, 6? possible ) + index of ID wanted
    for (edm::View<reco::GsfElectron>::const_iterator itElectron = electrons_->begin(); itElectron != electrons_->end(); itElectron++) {
        //auto idx = itElectron - electrons_->begin();//unsigned int
        //auto electronRef = electrons_->refAt(idx);//edm::RefToBase<reco::GsfElectron> 
        auto &electron = (*itElectron);
        if ( cfFlag("onlyEB") && std::abs(electron.eta()) > cfPrm("ebMaxEta") ) continue;
        //auto passIdCut = true; //electron.electronID(eleCutLoose);// pat electron ( miniAOD ) method
        //eleIdBools.push_back((*eleMVAIDLooseMap_)[electronRef]);// order is important, track how this vector is loaded
        auto passEnergyCut = electron.energy() > cfPrm("minEleE");
        if( passEnergyCut ) felectrons.push_back(electron);
    }//<<>>for( const auto electron : *electrons_ )


}//<<>>void KUCMSElectron::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSElectronObject::PostProcessEvent( ItemManager<float>& geVar ){}

void KUCMSElectronObject::ProcessEvent( ItemManager<float>& geVar ){

    if( ElectronDEBUG ) std::cout << "Processing Electrons" << std::endl;

    Branches.clearBranches();

    if( ElectronDEBUG ) std::cout << " - enetering Electron loop" << std::endl;

	int eleIndx = 0;
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

        Branches.fillBranch("SeedTOFTime",seedTOFTime);

	    eleIndx++;
        // GenParticle Info for electron  -------------------------------------------------------------------
        if( cfFlag("hasGenInfo") ){
            scptrs.push_back(*scptr);
            scptres.push_back(eleEnergy);
        }//<<>>if( hasGenInfo )

    }//<<>>for( const auto electron : *electrons_ )
	Branches.fillBranch("nEle",eleIndx);

    if( cfFlag("hasGenInfo") ){
        auto genInfo = genObjs->getGenEleMatch( scptrs, scptres );
        //if( ElectronDEBUG) 
		//std::cout << " Electron Gen Match ------------------------- " << std::endl;
        for( auto genidx : genInfo ){
			//std::cout << " -- Filling genindex : " << genidx << std::endl;
            Branches.fillBranch("GenIdx",genidx);
            if( genidx > -1 ){
                int genMomIndx = genObjs->getGenSigEleInfo( genidx );
                //std::cout << " -- Electron Match : level 1 : genidx: " << genidx << " mom: " << genMomIndx  << std::endl;
                Branches.fillBranch("GenWZIdx",genMomIndx);
				if( genMomIndx > -1 ){
					//uInt wzIndx = genMomIndx;
                	int genGMomIndx = genObjs->getGenSigEleInfo( genMomIndx );
                    //std::cout << "  -- Electron Match : level 2 : mom: " << genMomIndx << " gmom: " << genGMomIndx  << std::endl;
					Branches.fillBranch("GenMomIdx",genGMomIndx);
				} else { Branches.fillBranch("GenMomIdx",-5); }//<<>>if( genMomIndx > -1 )
            } else { Branches.fillBranch("GenWZIdx",-5); Branches.fillBranch("GenMomIdx",-5); }//<<>>if( genidx > -1 )
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
