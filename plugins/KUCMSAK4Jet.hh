// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS AK4Jet Object
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
//#include "DataFormats/AK4Jet/interface/AK4Jet.h"

#include "DataFormats/PatCandidates/interface/Jet.h"

// JECS
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JERs
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSElectron.hh"
#include "KUCMSPhoton.hh"

#ifndef KUCMSAK4JetObjectHeader
#define KUCMSAK4JetObjectHeader

//#define AK4JetEBUG true
#define AK4JetDEBUG false

using namespace edm; 

class KUCMSAK4JetObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSAK4JetObject( const edm::ParameterSet& iConfig );
    ~KUCMSAK4JetObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadAK4JetTokens( edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken ){ jetsToken_ = jetsToken; }; 
    void LoadPfcandTokens( edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken ){ pfcandToken_ = pfcandToken; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObject* genObjs_ ){ genObjs = genObjs_; };
    void LoadElectronObject( KUCMSElectronObject* electronObj_ ){ electronObj = electronObj_; };
    void LoadPhotonObject( KUCMSPhotonObject* photonObj_ ){ photonObj = photonObj_; };

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

    private:

    //std::vector<examnple> fexamples;
    std::vector<reco::PFJet> fjets;
    std::vector<int> jetIds;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    // Jets ----------------------------------------------------------

    //const edm::InputTag jetsTag;
    edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken_;
    edm::Handle<std::vector<reco::PFJet>> jets_;

    // PF Candidates View 
    //const edm::InputTag pfcandTag;
    edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken_;
    edm::Handle<edm::View<reco::Candidate>> pfcands_;

    // Jet corrections  ----------------------------------------------------------

    // JEC & JER  -> From miniAOD code not converted
    // JECs
    edm::ESHandle<JetCorrectorParametersCollection> jetCorrH;
    // JERs
    JME::JetResolution jetRes;
    JME::JetResolutionScaleFactor jetRes_sf;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    // exampleObject* otherObjectPtr;

    KUCMSEcalRecHitObject* rhObj;
    KUCMSGenObject* genObjs;
    KUCMSElectronObject* electronObj;
    KUCMSPhotonObject* photonObj;

};//<<>>class KUCMSAK4Jet : public KUCMSObjectBase

KUCMSAK4JetObject::KUCMSAK4JetObject( const edm::ParameterSet& iConfig ){   
// ---- end constructor initilizations  --------------------------

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );
    cfPrm.set( "jetPTmin", iConfig.existsAs<double>("jetPTmin") ? iConfig.getParameter<double>("jetPTmin") : 15.0 );
    cfPrm.set( "jetEtaMax", iConfig.existsAs<double>("jetEtaMax") ? iConfig.getParameter<double>("jetEtaMax") : 3.0 );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );

}//<<>>KUCMSAK4Jet::KUCMSAK4Jet( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSAK4JetObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("E","Jet_energy",VFLOAT);
    Branches.makeBranch("M","Jet_mass",VFLOAT);
    Branches.makeBranch("Pt","Jet_pt",VFLOAT);
    Branches.makeBranch("Eta","Jet_eta",VFLOAT);
    Branches.makeBranch("Phi","Jet_phi",VFLOAT);
    Branches.makeBranch("Area","Jet_area",VFLOAT,"jet catchment area, for JECs");

    Branches.makeBranch("nKids","Jet_nConstituents",VINT,"Number of particles in the jet");
    Branches.makeBranch("NHF","Jet_neHEF",VFLOAT,"neutral Hadron Energy Fraction");
    Branches.makeBranch("NEMF","Jet_neEmEF",VFLOAT,"neutral Electromagnetic Energy Fraction");
    Branches.makeBranch("CHF","Jet_chHEF",VFLOAT,"charged Hadron Energy Fraction");
    Branches.makeBranch("CEMF","Jet_chEmEF",VFLOAT,"charged Electromagnetic Energy Fraction");
    Branches.makeBranch("MUF","Jet_muEF",VFLOAT,"muon Energy Fraction");
    Branches.makeBranch("NHM","Jet_neHM",VFLOAT,"neutral hadron multiplicity");
    Branches.makeBranch("CHM","Jet_chHM",VFLOAT,"charged hadron multiplicity");

    Branches.makeBranch("Parts","Jet_egIndxs",VVUINT);
    Branches.makeBranch("DrRhIds","Jet_drRhIds",VVUINT);

    Branches.makeBranch("GenImpactAngle","Jet_genImpactAngle",VFLOAT);
    Branches.makeBranch("GenTime","Jet_genTime",VFLOAT);
    Branches.makeBranch("GenPt","Jet_genPt",VFLOAT);
    Branches.makeBranch("GenEta","Jet_genEta",VFLOAT);
    Branches.makeBranch("GenPhi","Jet_genPhi",VFLOAT);
    Branches.makeBranch("GenEnergy","Jet_genEnergy",VFLOAT);
    Branches.makeBranch("GenDrMatch","Jet_genDrMatch",VFLOAT);
    Branches.makeBranch("GenDptMatch","Jet_genDptMatch",VFLOAT);
    Branches.makeBranch("GenTimeLLP","Jet_genTimeLLP",VFLOAT);
    Branches.makeBranch("GenTOF","Jet_genTOF",VFLOAT);

    Branches.attachBranches(fOutTree);

}//<<>>void KUCMSAK4Jet::InitObject( TTree* fOutTree )

void KUCMSAK4JetObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // JETS
	iEvent.getByToken(jetsToken_, jets_);
    iEvent.getByToken(pfcandToken_, pfcands_);

    if( AK4JetDEBUG ) std::cout << "Collecting AK4Jets" << std::endl;
    fjets.clear();
    jetIds.clear();
    for(const auto &jet : *jets_ ){ // Filters jet collection & sorts by pt

        if( jet.pt() < cfPrm("jetPTmin") ) continue;
          if( cfFlag("onlyEB") && std::abs(jet.eta()) > cfPrm("ebMaxEta") ) continue;
        if( std::abs(jet.eta()) > cfPrm("jetEtaMax") ) continue;
          // save the jets, and then store the ID
          fjets.emplace_back(jet);
        auto jetID = 0;
        jetIds.push_back(jetID);

    }//<<>>for(const auto& jet : *jets_ )
	//std::sort(fjets.begin(),fjets.end(),sortByPt);


}//<<>>void KUCMSAK4Jet::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSAK4JetObject::ProcessEvent( ItemManager<float>& geVar ){

    if( AK4JetDEBUG ) std::cout << "Processing AK4Jets" << std::endl;

    Branches.clearBranches();

    if( AK4JetDEBUG ) std::cout << " - enetering AK4Jet loop" << std::endl;

	vector<uInt> noIdGrp(1,0);
    // **** const reco::TrackRefVector& associatedTracks() const;
    // **** const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
    // **** reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }

    auto nJets = fjets.size();
    float jetHt = 0;
    if(AK4JetDEBUG ) std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl;
    for ( uInt ijet(0); ijet < nJets; ijet++ ){

        const auto & jet = fjets[ijet];
        //const auto jetid = fjetsID[ijet];
        if(AK4JetDEBUG ) std::cout << " --- Proccesssing : " << jet  << std::endl;

        jetHt += jet.pt();

        const float jetE = jet.energy();
        const float jetM = jet.mass();
        const float jetPt = jet.pt();
        const float jetEta = jet.eta();
        const float jetPhi = jet.phi();
        const float jetArea = jet.jetArea();

        const float jetNHF = jet.neutralHadronEnergyFraction();
        const float jetNEMF = jet.neutralEmEnergyFraction();
        const float jetCHF = jet.chargedHadronEnergyFraction();
        const float jetCEMF = jet.chargedEmEnergyFraction();
        const float jetMUF = jet.muonEnergyFraction();
        const float jetNHM = jet.neutralMultiplicity();
        const float jetCHM = jet.chargedMultiplicity();

        Branches.fillBranch("E",jetE);
        Branches.fillBranch("M",jetM);
        Branches.fillBranch("Pt",jetPt);
        Branches.fillBranch("Eta",jetEta);
        Branches.fillBranch("Phi",jetPhi);
        Branches.fillBranch("Area",jetArea);

        Branches.fillBranch("NHF",jetNHF);
        Branches.fillBranch("NEMF",jetNEMF);
        Branches.fillBranch("CHF",jetCHF);
        Branches.fillBranch("CEMF",jetCEMF);
        Branches.fillBranch("MUF",jetMUF);
        Branches.fillBranch("NHM",jetNHM);
        Branches.fillBranch("CHM",jetCHM);

        // dR matched RH group  ----------------------------------------------- 
        // --------------------------------------------------------------------

        if(AK4JetDEBUG ) std::cout << "Getting jet dR rechit group" << std::endl;
        auto deltaRmin = 0.8;
        auto minRhE = 0.5;
        const auto jetDrRhGroup = rhObj->getRHGroup( jetEta, jetPhi, deltaRmin, minRhE );
        const auto jetDrRhIdsGroup = rhObj->getRhGrpIDs( jetDrRhGroup );
        Branches.fillBranch("DrRhIds",jetDrRhIdsGroup);
        rhObj->setRecHitUsed( jetDrRhIdsGroup );

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////
        // !!!!!!!!  create list of associated photons and electrons with jet via daughter particles of jet !!!!!!!!!!
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////

        const int nKids = jet.nConstituents();
        Branches.fillBranch("nKids",nKids);

        std::vector<uInt> jetParts;
        for( const auto &kid : jet.daughterPtrVector() ){

            auto kidcand = (pfcands_->ptrAt(kid.key())).get();
            auto kideta = kidcand->eta();
            auto kidphi = kidcand->phi();
            auto idx = photonObj->getIndex( kideta, kidphi );
            if( idx >= 0 ){ jetParts.push_back(idx); }
            else {
				idx = electronObj->getIndex( kideta, kidphi );
                if( idx >= 0 ) jetParts.push_back(1000+idx);
            }//<<>>if( phoidx >= 0 )

        }//<<>.for( const auto kid : jet.daughterPtrVector() )

        Branches.fillBranch("Parts",jetParts);

        // GenJet Info for MC  -------------------------------------------------------------------
        // ---------------------------------------------------------------------------------------

        if(AK4JetDEBUG ) std::cout << "Getting jetGenParton Information" << std::endl;


        if( cfFlag("hasGenInfo") ){

			auto genJetInfo = genObjs->getGenJetInfo( jetEta, jetPhi, jetPt );
            // load event level vectors for this jet with gen info

            Branches.fillBranch("GenImpactAngle",genJetInfo[0]);
            Branches.fillBranch("GenTime",genJetInfo[1]);
            Branches.fillBranch("GenPt",genJetInfo[2]);
            Branches.fillBranch("GenEta",genJetInfo[3]);
            Branches.fillBranch("GenPhi",genJetInfo[4]);
            Branches.fillBranch("GenEnergy",genJetInfo[5]);
            Branches.fillBranch("GenDrMatch",genJetInfo[7]);
            Branches.fillBranch("GenDptMatch",genJetInfo[8]);
            Branches.fillBranch("GenTimeLLP",genJetInfo[9]);
            Branches.fillBranch("GenTOF",genJetInfo[10]);

        }//<<>>if( hasGenInfo )
        if(AK4JetDEBUG ) std::cout << "Next Jet .......................... " << std::endl;

    }//<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
    // ** end of jets   ***************************************************************************************************

}//<<>>void KUCMSAK4Jet::ProcessEvent()

void KUCMSAK4JetObject::EndJobs(){}

#endif
