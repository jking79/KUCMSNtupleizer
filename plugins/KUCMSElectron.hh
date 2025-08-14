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
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"

//  KUCMS Object includes
#include "KUCMSObjectBase.hh"

// includes for other KUCMSObjects
#include "KUCMSEcalRechit.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSDisplacedVertex.hh"
#include "KUCMSNtupleizer/KUCMSNtupleizer/interface/TrackHelper.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#ifndef KUCMSElectronObjectHeader
#define KUCMSElectronObjectHeader

//#define ElectronDEBUG true
#define ElectronDEBUG false

using namespace edm; 

class KUCMSElectronObject : public KUCMSObjectBase {

    public:

	typedef pat::IsolatedTrack::PolarLorentzVector PolarLorentzVector;

    // use base class constructor
    KUCMSElectronObject( const edm::ParameterSet& iConfig );
    //KUCMSElectronObject( const edm::ParameterSet& iConfig );
    ~KUCMSElectronObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    // load tokens for eventt based collections
    void LoadElectronTokens( edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken ){ electronsToken_ = electronsToken; }; 
    void LoadConversionTokens( edm::EDGetTokenT<reco::ConversionCollection> conversionsToken ){ conversionsToken_ = conversionsToken; };
    void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    void LoadRhoTokens( edm::EDGetTokenT<double> RhoToken ){ RhoToken_ = RhoToken; };
    void LoadGeneralTrackTokens( edm::EDGetTokenT<edm::View<reco::Track>> token ){ generalTracksToken_ = token; }
    void LoadGsfTrackTokens( edm::EDGetTokenT<edm::View<reco::GsfTrack>> token ){ gsfTracksToken_ = token; }
    void LoadVertexTokens( edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_ ){ verticesToken = verticesToken_; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 
    void LoadRecHitObject( KUCMSEcalRecHitObject* rhObj_ ){ rhObj = rhObj_; }; // define with specific KUCMS object(s) needed 
    void LoadGenObject( KUCMSGenObject* genObjs_ ){ genObjs = genObjs_; };
	void LoadDisplacedVertexObject( KUCMSDisplacedVertex* svObj_ ){ svObj = svObj_; };
    void LoadTTrackBuilder(edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbuilder) {transientTrackBuilder_ = ttbuilder; }
    void LoadPfcandTokens( edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken ){ pfcandToken_ = pfcandToken; };

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
	void getTrkIso(const PolarLorentzVector& p4, int pc_idx, pat::PFIsolation& iso, pat::PFIsolation& miniiso) const;

    private:

    //std::vector<examnple> fexamples;
    std::vector<reco::GsfElectron> felectrons;
    std::vector<int> eleIds;

  	// General Tracks
  	edm::EDGetTokenT<edm::View<reco::Track>> generalTracksToken_;
  	edm::Handle<edm::View<reco::Track> > generalTracksHandle_;

    // GsfTracks
    edm::EDGetTokenT<edm::View<reco::GsfTrack>> gsfTracksToken_;
    edm::Handle<edm::View<reco::GsfTrack> > gsfTracksHandle_;

    //edm::EDGetTokenT<example> exampleToken_;
    //edm::Handle<example> example_;
    //const edm::InputTag electronsTag;
    //edm::EDGetTokenT<std::vector<reco::GsfElectron> > electronsToken_;
    edm::EDGetTokenT<edm::View<reco::GsfElectron>> electronsToken_;
    //edm::Handle<std::vector<reco::GsfElectron> > electrons_;
    edm::Handle<edm::View<reco::GsfElectron> > electrons_;
    //std::vector<reco::GsfElectron> electrons;
    edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> transientTrackBuilder_;

    //const edm::InputTag eleMVAIDLooseMapTag;
    //edm::EDGetTokenT<edm::ValueMap<bool>> eleMVAIDLooseMapToken_;
    //edm::Handle<edm::ValueMap<bool>> eleMVAIDLooseMap_;
	edm::EDGetTokenT<double> RhoToken_;
	edm::Handle<double>  rho_;

    // Tracks
    TTBuilderWrapper ttBuilder;

    // Conversions ( for electron veto ) ---------------------------------
    //const edm::InputTag conversionsTag;
    edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
    edm::Handle<reco::ConversionCollection> conversions_;

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    edm::Handle<reco::BeamSpot> beamSpot_;

	// PV
    edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken;
    edm::Handle<std::vector<reco::Vertex>> vertices_;

    // PF Candidates View 
    //const edm::InputTag pfcandTag;
    edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken_;
    edm::Handle<edm::View<reco::Candidate>> pfcands_;

	// for cut based IDs
	EffectiveAreas effectiveAreas_;
    //         run                   cut                    id                 eb?    value 
    std::map<std::string,std::map<std::string,std::map<std::string,std::map<bool,float>>>> cutValMap;

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
    prmStr.set( "RunEra", iConfig.existsAs<std::string>("RunEra") ? iConfig.getParameter<std::string>("RunEra") : "Run2" );

    //         run                   cut                    id                 eb?    value 
    //std::map<std::string,std::map<std::string,std::map<std::string,std::map<bool,float>>>> cutValMap;

	cutValMap["Run3"]["sigmaIeta"]["V"][true] = 0.0117; cutValMap["Run3"]["sigmaIeta"]["V"][false] = 0.0298; 
    cutValMap["Run3"]["dEtaSeed"]["V"][true] = 0.0071; cutValMap["Run3"]["dEtaSeed"]["V"][false] = 0.0173;
    cutValMap["Run3"]["deltaPhi"]["V"][true] = 0.208; cutValMap["Run3"]["deltaPhi"]["V"][false] = 0.234;
    cutValMap["Run3"]["HOE1"]["V"][true] = 0.05; cutValMap["Run3"]["HOE1"]["V"][false] = 0.05;
    cutValMap["Run3"]["HOE2"]["V"][true] = 1.28; cutValMap["Run3"]["HOE2"]["V"][false] = 2.3;
    cutValMap["Run3"]["HOE3"]["V"][true] = 0.0422; cutValMap["Run3"]["HOE3"]["V"][false] = 0.262;
    cutValMap["Run3"]["ooemoop"]["V"][true] = 0.178; cutValMap["Run3"]["ooemoop"]["V"][false] = 0.137;
    cutValMap["Run3"]["EA1"]["V"][true] = 0.406; cutValMap["Run3"]["EA1"]["V"][false] = 0.342;
    cutValMap["Run3"]["EA2"]["V"][true] = 0.535; cutValMap["Run3"]["EA2"]["V"][false] = 0.519;

    cutValMap["Run3"]["sigmaIeta"]["L"][true] = 0.0107; cutValMap["Run3"]["sigmaIeta"]["L"][false] = 0.0275;
    cutValMap["Run3"]["dEtaSeed"]["L"][true] = 0.00691; cutValMap["Run3"]["dEtaSeed"]["L"][false] = 0.0121;
    cutValMap["Run3"]["deltaPhi"]["L"][true] = 0.175; cutValMap["Run3"]["deltaPhi"]["L"][false] = 0.228;
    cutValMap["Run3"]["HOE1"]["L"][true] = 0.05; cutValMap["Run3"]["HOE1"]["L"][false] = 0.05;
    cutValMap["Run3"]["HOE2"]["L"][true] = 1.28; cutValMap["Run3"]["HOE2"]["L"][false] = 2.3;
    cutValMap["Run3"]["HOE3"]["L"][true] = 0.0422; cutValMap["Run3"]["HOE3"]["L"][false] = 0.262;
    cutValMap["Run3"]["ooemoop"]["L"][true] = 0.138; cutValMap["Run3"]["ooemoop"]["L"][false] = 0.127;
    cutValMap["Run3"]["EA1"]["L"][true] = 0.194; cutValMap["Run3"]["EA1"]["L"][false] = 0.184;
    cutValMap["Run3"]["EA2"]["L"][true] = 0.535; cutValMap["Run3"]["EA2"]["L"][false] = 0.519;

    cutValMap["Run3"]["sigmaIeta"]["M"][true] = 0.0103; cutValMap["Run3"]["sigmaIeta"]["M"][false] = 0.0272;
    cutValMap["Run3"]["dEtaSeed"]["M"][true] = 0.00481; cutValMap["Run3"]["dEtaSeed"]["M"][false] = 0.00951;
    cutValMap["Run3"]["deltaPhi"]["M"][true] = 0.127; cutValMap["Run3"]["deltaPhi"]["M"][false] = 0.221;
    cutValMap["Run3"]["HOE1"]["M"][true] = 0.0241; cutValMap["Run3"]["HOE1"]["M"][false] = 0.05;
    cutValMap["Run3"]["HOE2"]["M"][true] = 1.28; cutValMap["Run3"]["HOE2"]["M"][false] = 2.3;
    cutValMap["Run3"]["HOE3"]["M"][true] = 0.0422; cutValMap["Run3"]["HOE3"]["M"][false] = 0.262;
    cutValMap["Run3"]["ooemoop"]["M"][true] = 0.0966; cutValMap["Run3"]["ooemoop"]["M"][false] = 0.111;
    cutValMap["Run3"]["EA1"]["M"][true] = 0.0837; cutValMap["Run3"]["EA1"]["M"][false] = 0.0741;
    cutValMap["Run3"]["EA2"]["M"][true] = 0.535; cutValMap["Run3"]["EA2"]["M"][false] = 0.519;

    cutValMap["Run2"]["sigmaIeta"]["V"][true] = 0.0126; cutValMap["Run2"]["sigmaIeta"]["V"][false] = 0.0457;
    cutValMap["Run2"]["dEtaSeed"]["V"][true] = 0.00463; cutValMap["Run2"]["dEtaSeed"]["V"][false] = 0.00814;
    cutValMap["Run2"]["deltaPhi"]["V"][true] = 0.148; cutValMap["Run2"]["deltaPhi"]["V"][false] = 0.19;
    cutValMap["Run2"]["HOE1"]["V"][true] = 0.05; cutValMap["Run2"]["HOE1"]["V"][false] = 0.05;
    cutValMap["Run2"]["HOE2"]["V"][true] = 1.16; cutValMap["Run2"]["HOE2"]["V"][false] = 2.54;
    cutValMap["Run2"]["HOE3"]["V"][true] = 0.0324; cutValMap["Run2"]["HOE3"]["V"][false] = 0.183;
    cutValMap["Run2"]["ooemoop"]["V"][true] = 0.209; cutValMap["Run2"]["ooemoop"]["V"][false] = 0.132;
    cutValMap["Run2"]["EA1"]["V"][true] = 0.198; cutValMap["Run2"]["EA1"]["V"][false] = 0.203;
    cutValMap["Run2"]["EA2"]["V"][true] = 0.506; cutValMap["Run2"]["EA2"]["V"][false] = 0.963;

    cutValMap["Run2"]["sigmaIeta"]["L"][true] = 0.0112; cutValMap["Run2"]["sigmaIeta"]["L"][false] = 0.0425;
    cutValMap["Run2"]["dEtaSeed"]["L"][true] = 0.00377; cutValMap["Run2"]["dEtaSeed"]["L"][false] = 0.00674;
    cutValMap["Run2"]["deltaPhi"]["L"][true] = 0.0884; cutValMap["Run2"]["deltaPhi"]["L"][false] = 0.169;
    cutValMap["Run2"]["HOE1"]["L"][true] = 0.05; cutValMap["Run2"]["HOE1"]["L"][false] = 0.0441;
    cutValMap["Run2"]["HOE2"]["L"][true] = 1.16; cutValMap["Run2"]["HOE2"]["L"][false] = 2.54;
    cutValMap["Run2"]["HOE3"]["L"][true] = 0.0324; cutValMap["Run2"]["HOE3"]["L"][false] = 0.183;
    cutValMap["Run2"]["ooemoop"]["L"][true] = 0.193; cutValMap["Run2"]["ooemoop"]["L"][false] = 0.111;
    cutValMap["Run2"]["EA1"]["L"][true] = 0.112; cutValMap["Run2"]["EA1"]["L"][false] = 0.108;
    cutValMap["Run2"]["EA2"]["L"][true] = 0.506; cutValMap["Run2"]["EA2"]["L"][false] = 0.963;

    cutValMap["Run2"]["sigmaIeta"]["M"][true] = 0.0106; cutValMap["Run2"]["sigmaIeta"]["M"][false] = 0.0387;
    cutValMap["Run2"]["dEtaSeed"]["M"][true] = 0.0032; cutValMap["Run2"]["dEtaSeed"]["M"][false] = 0.00632;
    cutValMap["Run2"]["deltaPhi"]["M"][true] = 0.0547; cutValMap["Run2"]["deltaPhi"]["M"][false] = 0.0394;
    cutValMap["Run2"]["HOE1"]["M"][true] = 0.046; cutValMap["Run2"]["HOE1"]["M"][false] = 0.0275;
    cutValMap["Run2"]["HOE2"]["M"][true] = 1.16; cutValMap["Run2"]["HOE2"]["M"][false] = 2.52;
    cutValMap["Run2"]["HOE3"]["M"][true] = 0.0324; cutValMap["Run2"]["HOE3"]["M"][false] = 0.183;
    cutValMap["Run2"]["ooemoop"]["M"][true] = 0.184; cutValMap["Run2"]["ooemoop"]["M"][false] = 0.0721;
    cutValMap["Run2"]["EA1"]["M"][true] = 0.0478; cutValMap["Run2"]["EA1"]["M"][false] = 0.0658;
    cutValMap["Run2"]["EA2"]["M"][true] = 0.506; cutValMap["Run2"]["EA2"]["M"][false] = 0.963;

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

    //Branches.makeBranch("nIsoTracks","PF_nIsoTracks",INT);
    //Branches.makeBranch("absiso","PF_absIso",VFLOAT);
    //Branches.makeBranch("reliso","PF_relIso",VFLOAT);
    //Branches.makeBranch("minireliso","PF_miniRelIso",VFLOAT);


    if( cfFlag("hasGenInfo") ){
    Branches.makeBranch("GenIdx","Electron_genIdx",VINT);
    Branches.makeBranch("GenXMomIdx","Electron_genSigXMomId",VINT);
    Branches.makeBranch("GenWZIdx","Electron_genSigWZId",VINT);
    Branches.makeBranch("hasGenIdx","Electron_hasGenMatch",VBOOL);
    Branches.makeBranch("isLWZP", "Electron_isLWZP", VBOOL );
	}//<<>>if( cfFlag("hasGenInfo") )
    Branches.makeBranch("nGenMatch","Electron_nGenMatch",INT);

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

    Branches.makeBranch("IsLVeto","Electron_isLVeto",VBOOL);
    Branches.makeBranch("Vsieiec","Electron_Vsieiec",VBOOL);
    Branches.makeBranch("Vdetaseed","Electron_Vdetaseed",VBOOL);
    Branches.makeBranch("Vdpscvc","Electron_Vdpscvc",VBOOL);
    Branches.makeBranch("Vhoec","Electron_Vhoec",VBOOL);
    Branches.makeBranch("VrelIsoWithEA","Electron_VrelIsoWithEA",VBOOL);
    Branches.makeBranch("Vooemoop","Electron_Vooemoop",VBOOL);
    Branches.makeBranch("passConversionVeto","Electron_passConversionVeto",VBOOL);
    Branches.makeBranch("Vmistrack","Electron_Vmistrack",VBOOL);

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


	iEvent.getByToken( electronsToken_, electrons_ );
    iEvent.getByToken( conversionsToken_, conversions_ );
    iEvent.getByToken( beamLineToken_, beamSpot_ );
    iEvent.getByToken( RhoToken_, rho_ );
    iEvent.getByToken( pfcandToken_, pfcands_ );
  	iEvent.getByToken( generalTracksToken_, generalTracksHandle_ );
    iEvent.getByToken( gsfTracksToken_, gsfTracksHandle_ );
    iEvent.getByToken( verticesToken, vertices_ );

    ttBuilder = TTBuilderWrapper( &iSetup.getData(transientTrackBuilder_) );

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
		float eta = std::abs( electron.eta() );
        bool passEnergyCut = electron.energy() > cfPrm("minEleE");
        bool passEtaCut = ( eta < 1.444 ) || ( ( eta > 1.566 ) && ( eta < 2.5 ) );
        if( passEnergyCut && passEtaCut ) felectrons.push_back(electron);
    }//<<>>for( const auto electron : *electrons_ )

}//<<>>void KUCMSElectron::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSElectronObject::PostProcessEvent( ItemManager<float>& geVar ){

}//<<>>void KUCMSElectronObject::PostProcessEvent( ItemManager<float>& geVar )

void KUCMSElectronObject::ProcessEvent( ItemManager<float>& geVar ){

    if( ElectronDEBUG ) std::cout << "Processing Electrons" << std::endl;

    Branches.clearBranches();

    if( ElectronDEBUG ) std::cout << " - enetering Electron loop" << std::endl;

	int eleIndx = 0;
	int nSelEle = 0;
	int nSVMatched = 0;
    int nPromptEle = 0;
    //scGroup scptrs;
	std::vector<v3fPoint> scvertex;
    std::vector<float> scptres;
    std::vector<reco::TransientTrack> sctrks;
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

		std::string re = prmStr("RunEra");
		std::map<std::string,std::map<std::string,bool>> id;

		float dEtaSeed = std::abs( electron.deltaEtaSuperClusterTrackAtVtx() - scptr->eta() + scptr->seed()->eta() );
		float epinv = std::abs( 1.0 - electron.eSuperClusterOverP() )/electron.ecalEnergy();
		auto eleMisTrack = electron.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);

		const reco::GsfElectron::PflowIsolationVariables& pfIso = electron.pfIsolationVariables();
		const float chad = pfIso.sumChargedHadronPt;
		const float nhad = pfIso.sumNeutralHadronEt;
		const float pho = pfIso.sumPhotonEt;
		double absEta = std::abs(electron.superCluster()->eta());
		float  eA = effectiveAreas_.getEffectiveArea( absEta );
        float rho = rho_.isValid() ? (float)(*rho_) : 0;
		//std::cout << " -- EA : rho " << rho << " pho " << pho << " chad " << chad << " nhad " << nhad << " eA " << eA << std::endl;
		float relIsoWithEA = (chad + std::max(0.0f, nhad + pho - rho*eA))/elePt;
        bool ebee = ( eleEta <= 1.479 );

		bool passConversionVeto = !ConversionTools::hasMatchedConversion( electron, *conversions_, beamSpot_->position());
        bool mistrack = ( eleMisTrack <= 1 );
		bool Vmistrack = ebee ? ( eleMisTrack <= 2 ) : ( eleMisTrack <= 3 );
		bool CvTr = mistrack && passConversionVeto;
        bool VCvTr = Vmistrack && passConversionVeto;

		std::vector<std::string> idLev = { "V", "L", "M" };
		std::vector<std::string> idVar = { "sigmaIeta", "dEtaSeed", "deltaPhi", "HOE", "EA", "ooemoop" };

		float sce = scptr->energy();
		for( auto idl : idLev ){

			id["sigmaIeta"][idl] = electron.full5x5_sigmaIetaIeta() < cutValMap[re]["sigmaIeta"][idl][ebee];
        	id["dEtaSeed"][idl] = dEtaSeed < cutValMap[re]["dEtaSeed"][idl][ebee];
        	id["deltaPhi"][idl] = electron.deltaPhiSuperClusterTrackAtVtx() < cutValMap[re]["deltaPhi"][idl][ebee];
			float hoev =  cutValMap[re]["HOE1"][idl][ebee] + cutValMap[re]["HOE2"][idl][ebee]/sce + cutValMap[re]["HOE3"][idl][ebee]*rho/sce;
			id["HOE"][idl] = electron.hadronicOverEm() < hoev;
        	id["EA"][idl] = relIsoWithEA < ( cutValMap[re]["EA1"][idl][ebee] + cutValMap[re]["EA2"][idl][ebee]/elePt );
        	id["ooemoop"][idl] = epinv < cutValMap[re]["ooemoop"][idl][ebee];

		}//<<>>for( auto idl : idLev )

		Branches.fillBranch("Vsieiec",id["sigmaIeta"]["V"]);
		Branches.fillBranch("Vdetaseed",id["dEtaSeed"]["V"]);
		Branches.fillBranch("Vdpscvc",id["deltaPhi"]["V"]);
		Branches.fillBranch("Vhoec",id["HOE"]["V"]);
		Branches.fillBranch("VrelIsoWithEA",id["EA"]["V"] );
		Branches.fillBranch("Vooemoop",id["ooemoop"]["V"] );
        Branches.fillBranch("passConversionVeto",passConversionVeto);
        Branches.fillBranch("Vmistrack",Vmistrack);

		//if( id["EA"]["V"] ) std::cout << " -- EA : rho " << rho << " pho " << pho << " chad " << chad << " nhad " << nhad << " eA " << eA << std::endl;

        bool septc = ( elePt >= 10 );
		bool veptc = ( elePt >= 20 );

		std::map<std::string,bool> idBool;
		for( auto lev : idLev ){ idBool[lev] = true; for( auto var : idVar ){ idBool[lev] = idBool[lev] && id[var][lev]; }}
		idBool["L"] = idBool["L"] && CvTr;
        idBool["M"] = idBool["M"] && CvTr;
        idBool["V"] = idBool["V"] && VCvTr;
        bool isLVeto = id["dEtaSeed"]["V"] && id["deltaPhi"]["V"] && id["ooemoop"]["V"] && id["HOE"]["V"] && VCvTr;

		if( septc && idBool["L"] ) nSelEle++;

		Branches.fillBranch("IsLVeto",isLVeto);
        Branches.fillBranch("IsVeto",idBool["V"]);
        Branches.fillBranch("IsLoose",idBool["L"]);
        Branches.fillBranch("IsMedium",idBool["M"]);

        Branches.fillBranch("SeedTOFTime",seedTOFTime);
    	Branches.fillBranch("Sieie",electron.full5x5_sigmaIetaIeta());
    	Branches.fillBranch("DetaSCTV",electron.deltaEtaSuperClusterTrackAtVtx());
    	Branches.fillBranch("DphiSCTV",electron.deltaPhiSuperClusterTrackAtVtx());
    	Branches.fillBranch("HOE",electron.hadronicOverEm());
		
		bool hasSVMatch = false;
		if( cfFlag("doSVModule") ) hasSVMatch = svObj->FoundLeptonMatch( electron );
		if( hasSVMatch ) nSVMatched++;
        Branches.fillBranch("svMatch",hasSVMatch);

        bool passPrompt = false;
        if( cfFlag("doSVModule") ) passPrompt = svObj->IsPromptLepton( electron );
        Branches.fillBranch("passPrompt",passPrompt);

        bool isPromptEle = false;
        //if( passPrompt && veptc && isVeto ){ isPromptEle = true; nPromptEle++; }
        if( veptc && isLVeto ){ isPromptEle = true; nPromptEle++; }
        Branches.fillBranch("IsPromptEle",isPromptEle);

	    eleIndx++;

        //std::cout << " -- Muon getting sc : " << std::endl;
        const auto btptr =electron.gsfTrack();
        auto ttrack = ttBuilder.build( *btptr );

        // GenParticle Info for electron  -------------------------------------------------------------------
        if( cfFlag("hasGenInfo") ){
            //scptrs.push_back(*scptr);
            sctrks.push_back( ttrack );
			//v3fPoint scv( scptr->x(), scptr->y(), scptr->z() );
            v3fPoint scv( eleEta, elePhi, 0.f );
			scvertex.push_back( scv );
            scptres.push_back(eleEnergy);
        }//<<>>if( hasGenInfo )

    }//<<>>for( const auto electron : *electrons_ )
	Branches.fillBranch("nEle",eleIndx);
    Branches.fillBranch("nSele",nSelEle);
    Branches.fillBranch("nSVMatch",nSVMatched);
	Branches.fillBranch("nPrompt",nPromptEle);
	geVar.set("nSelEle",nSelEle);

    int nGenMatched = 0;
    if( cfFlag("hasGenInfo") ){
        auto genInfo = genObjs->getGenEleMatch( scvertex, scptres, sctrks );
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
