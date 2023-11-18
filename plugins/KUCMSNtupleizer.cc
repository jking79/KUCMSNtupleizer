// -*- C++ -*-
//
// Package:    KUCMSNtupilizer
// Class:      KUCMSNtupilizer
//
/**\class KUCMSNtupilizer KUCMSNtupilizer.cc LLPgammaAnalyzer/plugins/KUCMSNtupilizer.cc

         Description: [one line class summary]

         Implementation: [Notes on implementation]

*/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//
//
//----------------------------------------  cc file   --------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include "KUCMSNtupleizer.hh"
#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"

#include "KUCMSEventInfo.hh"
#include "KUCMSEcalRechit.hh"
#include "KUCMSPhoton.hh"
#include "KUCMSAK4Jet.hh"
#include "KUCMSPFMet.hh"
#include "KUCMSElectron.hh"
#include "KUCMSMuon.hh"
#include "KUCMSSecondaryVertex.hh"
#include "KUCMSTimedSVs.hh"
//#include "KUCMSConversions.hh"
#include "KUCMSGenObjects.hh"

using namespace std;

//#define DEBUG true
#define DEBUG false

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig):

    // ECAL RECORDS
    caloGeometryToken_(esConsumes()),
    caloTopologyToken_{esConsumes()},
    ecalLaserDbServiceToken_(esConsumes()),
    ecalIntercalibConstantsToken_(esConsumes()),
    ecalADCToGeVConstantToken_(esConsumes()),
    EcalPedestalsToken_(esConsumes())

{ 

    usesResource();
    usesResource("TFileService");

    // -- set input flags and parmeters  ----------------------------------------------------------
    // accuire input values form python config master file ( the one we use cmsRun with interactivilly )

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );
    cfPrm.set( "minEvtMet", iConfig.existsAs<double>("minEvtMet") ? iConfig.getParameter<double>("minEvtMet") : 150.0 );

    // -- consume tags ------------------------------------------------------------
    // creats "token" for all collections that we wish to process
    
    //if( DEBUG ) std::cout << "In constructor for KUCMSNtupilizer - tag and tokens" << std::endl;

    // Triggers
    //triggerResultsToken_ = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
    //triggerObjectsToken_ = consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"));
    // tracks 
    //tracksToken_ = consumes<std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"));
    // rho
    //rhoToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));

    if( DEBUG ) std::cout << "Create Object Classes" << std::endl;

    auto eventInfoObj = new KUCMSEventInfoObject(  iConfig ); 
    auto vertexToken = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
    eventInfoObj->LoadVertexTokens( vertexToken );

    auto recHitsObj = new KUCMSEcalRecHitObject( iConfig );
    auto rhEBtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEB"));
    auto rhEEtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEE"));
    recHitsObj->LoadRecHitTokens( rhEBtoken, rhEEtoken );
    auto sctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters"));
    auto ootsctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters"));
    recHitsObj->LoadSCTokens( sctoken, ootsctoken );
    auto ccltoken = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("caloClusters"));
    recHitsObj->LoadClusterTokens( ccltoken );

    auto electronsObj = new KUCMSElectronObject( iConfig );
    auto electronToken = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"));
    electronsObj->LoadElectronTokens( electronToken );
    auto conversionsToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
    electronsObj->LoadConversionTokens( conversionsToken );
    auto beamLineToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    electronsObj->LoadBeamSpotTokens( beamLineToken );
    electronsObj->LoadRecHitObject( recHitsObj );

    //muons
    auto muonsObj = new KUCMSMuon(iConfig);
    auto muonToken = consumes<std::vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
    muonsObj->LoadMuonTokens(muonToken);

    auto photonsObj = new KUCMSPhotonObject( iConfig );
    auto photonToken = consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("gedPhotons"));
    auto ootPhotonToken = consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("ootPhotons"));
    photonsObj->LoadPhotonTokens( photonToken, ootPhotonToken );
    photonsObj->LoadRecHitObject( recHitsObj );
    photonsObj->LoadElectronObject( electronsObj );
    
    auto ak4jetObj = new KUCMSAK4JetObject( iConfig );
    auto jetsToken = consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"));
    ak4jetObj->LoadAK4JetTokens( jetsToken );
    edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pfcandidates"));
    ak4jetObj->LoadPfcandTokens( pfcandToken );
    ak4jetObj->LoadRecHitObject( recHitsObj );
    ak4jetObj->LoadPhotonObject( photonsObj );
    ak4jetObj->LoadElectronObject( electronsObj );

    // SVs
    auto svsObj = new KUCMSSecondaryVertex(iConfig);
    auto svsToken = consumes<edm::View<reco::VertexCompositePtrCandidate>>(iConfig.getParameter<edm::InputTag>("svSrc"));
    svsObj->LoadSVTokens(svsToken);
    // Timed SVs
    auto timedSVsObj = new KUCMSTimedSVs(iConfig);
    auto timedSVsToken = consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("timedVertices"));
    auto lowPtElectronToken = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("lowPtElectrons"));
    timedSVsObj->LoadSVTokens(timedSVsToken);
    timedSVsObj->LoadElectronTokens(lowPtElectronToken);

    //Conversions
    /*
    auto conversionsObj = new KUCMSConversions(iConfig);
    auto conversionToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
    conversionsObj->LoadConversionTokens(conversionToken);
    conversionsObj->LoadElectronObject(electronsObj);
    conversionsObj->LoadRecHitObject(recHitsObj);
    */
    auto pfmetObj = new KUCMSPFMetObject( iConfig );
    auto pfmetToken = consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("mets")); 
    pfmetObj->LoadPFMetTokens( pfmetToken );
    pfmetObj->LoadPhotonObject( photonsObj );
    
    if( DEBUG ) std::cout << "Loading Object Manager" << std::endl;

    ObjMan.Load( "EventInfo", eventInfoObj );
    ObjMan.Load( "Electrons", electronsObj );
    ObjMan.Load( "Muons", muonsObj );
    ObjMan.Load( "Photons", photonsObj );
    ObjMan.Load( "JetsAK4", ak4jetObj );
    ObjMan.Load( "SVs", svsObj );
    ObjMan.Load( "TimedSVs", timedSVsObj);
    //ObjMan.Load( "Conversions", conversionsObj );
    ObjMan.Load( "PFMet", pfmetObj );
    ObjMan.Load( "ECALRecHits", recHitsObj );

    if( cfFlag("hasGenInfo") ){

        if( DEBUG ) std::cout << "Create & Load Gen Information" << std::endl;

        auto genObjs = new KUCMSGenObject( iConfig );
        auto genEvtInfoToken = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvt"));
        genObjs->LoadGenEvtInfoTokens( genEvtInfoToken );
        auto gent0Token = consumes<float>(iConfig.getParameter<edm::InputTag>("gent0"));
        genObjs->LoadGenT0Tokens( gent0Token );
        auto genxyz0Token = consumes<Point3D>(iConfig.getParameter<edm::InputTag>("genxyz0"));
        genObjs->LoadGenXYZ0Tokens( genxyz0Token );
        auto pileupInfosToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileups"));
        genObjs->LoadGenPUInfosTokens( pileupInfosToken );
        auto genPartToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
        genObjs->LoadgTokens( genPartToken );
        auto genJetsToken = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
        genObjs->LoadGenJetsTokens( genJetsToken );

        // Load gen object into other collections
	electronsObj->LoadGenObject( genObjs );
        photonsObj->LoadGenObject( genObjs );
        ak4jetObj->LoadGenObject( genObjs );

        // Load gen object into objman last, should be no dependence with other objects
        ObjMan.Load( "GenObjects", genObjs );

    }//<<>>if( hasGenInfo )

// ---------------------------------------------------------------------------------
}//>>>>KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ Event Selection   -----------------------------
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool KUCMSNtupilizer::selectedEvent(){

    bool hasMinMet = geVar("evtMET") > cfPrm("minEvtMet");

    bool selected = hasMinMet;
    return selected;

}//<<>>bool KUCMSNtupilizer::selectedEvent()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   -------------------  Normally nothing needs to be modified below this point  --------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ Destructor -----------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////i

KUCMSNtupilizer::~KUCMSNtupilizer(){

    ///////////////////////////////////////////////////////////////////
    // do anything here that needs to be done at desctruction time   //
    // (e.g. close files, deallocate resources etc.)                 //
    ///////////////////////////////////////////////////////////////////

}//>>>>KUCMSNtupilizer::~KUCMSNtupilizer()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ Analyzer Inherited Class Functions ------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void KUCMSNtupilizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    using namespace edm;

    // -- Consume Tokens --------------------------------------------
    // gets pointer to the collections from cmssw using the "token" for that collection

    //if( DEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

    // TRIGGER
    //iEvent.getByToken(triggerResultsToken_,triggerResults_);
    //iEvent.getByToken(triggerObjectsToken_,triggerObjects_);
    // TRACKS
    //iEvent.getByToken(tracksToken_, tracks_);
    // RHO
    //iEvent.getByToken(rhoToken_, rho_);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ---  Collecting objests ( preprocessing object pruning ) ---------------------------------------
    // -------------------------------------------------------------------------------------------------
    // -- Process Event  ---------------------------------------    
    // ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
    // ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << "***********************************************************************************************" << std::endl;

    // clear global event varibles 
    geVar.clear(); // floats

    // -------------------------------------------------------------------------------------------------
    // ---  Collecting objects ( preprocessing object pruning & fill global object vectors )------------
    // -------------------------------------------------------------------------------------------------

    if( DEBUG ) std::cout << "LoadEvent ObjMan" << std::endl;
	ObjMan.LoadEvent( iEvent, iSetup, geVar );

    //------------------------------------------------------------------------------------
    // ----   Object processing ----------------------------------------------------------
    //------------------------------------------------------------------------------------
    // call functions to process collections and fill tree varibles to be saved
    // varibles to be saved to ttree are declared in the header
    // use LoadEvent() for any processing that must be done before crosstalk 
    // use PostProcessEvent() for any processing that must be done after crosstalk

    if( DEBUG ) std::cout << "ProcessEvent ObjMan" << std::endl;
	ObjMan.ProcessEvent( geVar );
    ObjMan.PostProcessEvent( geVar );

    //------------------------------------------------------------------------------------
    //---- Object processing Completed ----------------------------------------------------------
    //------------------------------------------------------------------------------------

    // -- Fill output trees ------------------------------------------

    if( DEBUG ) std::cout << "Select Event and Fill Tree" << std::endl;
    //if( selectedEvent() ) {
    //std::cout << "filled events" << std::endl;
    outTree->Fill();
    //}
    //else 
    //std::cout << "events not filled" << std::endl;

    // -- EOFun ------------------------------------------------------
    //     #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //     ESHandle<SetupData> pSetup;
    //     iSetup.get<SetupRecord>().get(pSetup);
    //     #endif

}//>>>>void KUCMSNtupilizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ beginJob/endJob methods called once each job just before/after starting event loop    ------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void KUCMSNtupilizer::beginJob(){

    // Global Varibles

    // Book output files and trees
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree>("llpgtree","KUCMSNtuple");

    // Book //histograms ( if any )
    
    std::cout << "Services Booked" << std::endl;

    if( DEBUG ) std::cout << "Init ObjMan" << std::endl;
	ObjMan.Init( outTree );

}//>>>>void KUCMSNtupilizer::beginJob()


// ------------ method called once each job just after ending the event loop    ------------
void KUCMSNtupilizer::endJob(){ 

    if( DEBUG ) std::cout << "ObjMan EndJobs" << std::endl;
    ObjMan.EndJobs(); 

}//>>>>void KUCMSNtupilizer::endJob()

// ------------ method fills 'descriptions' with the allowed parameters for the module    ------------
void KUCMSNtupilizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);

    //Specify that only 'tracks' is allowed
    //To use, remove the default given above and uncomment below
    //ParameterSetDescription desc;
    //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
    //descriptions.addDefault(desc);

}//>>>>void KUCMSNtupilizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)

////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------- CMSSW Ana Helper Functions ---------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

/*  Not need at this level
vector<float> KUCMSNtupilizer::getTimeDistStats( vector<float> times, rhGroup rechits ){

    //if( rechits.size() == 0 ){ vector<float> result{-99.0}; return result; }

    // N 3.64, C 0.3000  s^2 = (N/(rhe))^2 + 2C^2
    float N(3.64);
    float C(0.3000);

    vector<float> wts;
    for( uInt it(0); it < rechits.size(); it++ ){
    //auto wt = 1/std::sqrt(sq2(N/rechits[it].energy())+2*sq2(C)); 
        auto wt = 1/(sq2(N/rechits[it].energy())+2*sq2(C));
        wts.push_back(wt);
    }//<<>>for( uInt it(0); it < rechits.size(); it++ )

    return getDistStats( times, wts );

}//>>>>vector<float> KUCMSNtupilizer::getTimeDistStats( vector<float> times, rhGroup rechits ){
*/

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ define this as a plug-in --------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(KUCMSNtupilizer);

