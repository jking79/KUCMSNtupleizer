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

#include "KUCMSNtupleizer_mini.hh"
//#include "KUCMSHelperFunctions.hh"
//#include "KUCMSRootHelperFunctions.hh"

#include "KUCMSEventInfo_mini.hh"
#include "KUCMSEcalRechit_mini.hh"
#include "KUCMSPhoton_mini.hh"
#include "KUCMSAK4Jet_mini.hh"
#include "KUCMSPFMet_mini.hh"
//#include "KUCMSECALTracks.hh"
#include "KUCMSElectron_mini.hh"
//#include "KUCMSMuon_mini.hh"
//#include "KUCMSDisplacedElectron.hh"
#include "KUCMSGenObjects_mini.hh"
#include "KUCMSDisplacedVertex.hh"
#include "KUCMSTrack_mini.hh"

using namespace std;

//#define NTDEBUG true
#define NTDEBUG false

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// constructor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KUCMSNtupilizerMini::KUCMSNtupilizerMini(const edm::ParameterSet& iConfig):

    // ECAL RECORDS
    caloGeometryToken_(esConsumes()),
    caloTopologyToken_{esConsumes()},
    ecalLaserDbServiceToken_(esConsumes()),
    ecalIntercalibConstantsToken_(esConsumes()),
    ecalADCToGeVConstantToken_(esConsumes()),
    EcalPedestalsToken_(esConsumes())

{ 

    //usesResource();
    usesResource("TFileService");

    // -- set input flags and parmeters  ----------------------------------------------------------
    // accuire input values form python config master file ( the one we use cmsRun with interactivilly )

    cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );
    cfFlag.set( "doDisEleModule", iConfig.existsAs<bool>("doDisEleModule") ? iConfig.getParameter<bool>("doDisEleModule") : false );
    cfFlag.set( "doECALTrackOnly", iConfig.existsAs<bool>("doECALTrackOnly") ? iConfig.getParameter<bool>("doECALTrackOnly") : false );

    // -- consume tags ------------------------------------------------------------
    // creats "token" for all collections that we wish to process
    
    //if( NTDEBUG ) 
    std::cout << "In constructor for KUCMSNtupilizer - tag and tokens" << std::endl;

	// -----  set event skim selection ----------------------------------------------------------------
	std::string filterSelect(iConfig.existsAs<std::string>("fltrSelection")? iConfig.getParameter<std::string>("fltrSelection") : "None" );
	ntupleFilter.setEventSelectionTag( filterSelect );

	if( NTDEBUG ) std::cout << "Create Object Classes" << std::endl;
	// note : objects stored in a map - the order in which the objects are called is set by the alphanumeric order of the object names


	//Event Info
    auto eventInfoObj = new KUCMSEventInfoObjectMini(  iConfig ); 
    auto vertexToken = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
    auto triggerFlagResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerFlagResults"));
    auto triggerHltResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerHLTResults"));
    eventInfoObj->LoadVertexTokens( vertexToken );
    eventInfoObj->LoadTriggerTokens( triggerFlagResultsToken, triggerHltResultsToken );
    ObjMan.Load( "EventInfo", eventInfoObj );

    //Rechits ECAL
    auto recHitsObj = new KUCMSEcalRecHitObjectMini( iConfig );
    auto rhEBtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEB"));
    auto rhEEtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEE"));
    recHitsObj->LoadRecHitTokens( rhEBtoken, rhEEtoken );
    auto sctoken  = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters"));
    auto ootsctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters"));
    auto othersctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("otherSuperClusters"));
	recHitsObj->LoadSCTokens( sctoken, ootsctoken );
    recHitsObj->LoadSCTokens( othersctoken );
    auto ccltoken = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("caloClusters"));
    recHitsObj->LoadClusterTokens( ccltoken );
    auto beamLineToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    recHitsObj->LoadBeamSpotTokens( beamLineToken );

	/*
    //Track General + GSF
    auto TrackObj = new KUCMSTrackObjectMini( iConfig );
    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcandToken = consumes<edm::View<reco::Candidate>>(iConfig.getParameter<edm::InputTag>("pfcandidates"));
    TrackObj->LoadPfcandTokens( pfcandToken );
    auto parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
	TrackAssociatorParameters trackAssocParameters;
	edm::ConsumesCollector iC = consumesCollector();
	trackAssocParameters.loadParameters(parameters, iC);
	TrackObj->LoadAssociationParameters(trackAssocParameters);
    auto transientTrackBuilderToken = esConsumes<TransientTrackBuilder,TransientTrackRecord>(edm::ESInputTag("","TransientTrackBuilder"));
    TrackObj->LoadTTrackBuilder(transientTrackBuilderToken);
    auto ogGeneralTracksToken = consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("ogGeneralTracks"));
    TrackObj->LoadGeneralTrackTokens(ogGeneralTracksToken);
    auto ogGsfTracksToken = consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("ogGsfTracks"));
    TrackObj->LoadGsfTrackTokens(ogGsfTracksToken);
	auto magneticFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
	TrackObj->LoadMagneticField(magneticFieldToken);
    TrackObj->LoadBeamSpotTokens( beamLineToken );
    TrackObj->LoadVertexTokens( vertexToken );
    TrackObj->LoadRecHitObject( recHitsObj );
    ObjMan.Load( "Tracks", TrackObj );
	*/

    edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcandToken = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfcandidates"));
    auto transientTrackBuilderToken = esConsumes<TransientTrackBuilder,TransientTrackRecord>(edm::ESInputTag("","TransientTrackBuilder"));
    auto magneticFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
    auto parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
    TrackAssociatorParameters trackAssocParameters;
    edm::ConsumesCollector iC = consumesCollector();
    trackAssocParameters.loadParameters(parameters, iC);

    //Electrons 
    auto electronsObj = new KUCMSElectronObjectMini( iConfig );
    auto rhoToken = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
    electronsObj->LoadRhoTokens( rhoToken );
    auto electronToken = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));
    electronsObj->LoadElectronTokens( electronToken );
    auto conversionsToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
    electronsObj->LoadConversionTokens( conversionsToken );
    electronsObj->LoadPfcandTokens( pfcandToken );
	electronsObj->LoadTTrackBuilder(transientTrackBuilderToken);
    //electronsObj->LoadGeneralTrackTokens(ogGeneralTracksToken);
    auto ogGsfTracksToken = consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("ogGsfTracks"));
    electronsObj->LoadGsfTrackTokens(ogGsfTracksToken);
    electronsObj->LoadBeamSpotTokens( beamLineToken );
    electronsObj->LoadVertexTokens( vertexToken );
    electronsObj->LoadRecHitObject( recHitsObj );
    ObjMan.Load( "Electrons", electronsObj );

/*
    issue with edm::View<reco::Muon>::const_iterator 
    //Muons    
    auto muonObj = new KUCMSMuonObject( iConfig );
    auto muonToken = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
    muonObj->LoadVertexTokens( vertexToken );
    muonObj->LoadMuonTokens( muonToken );
	muonObj->LoadTTrackBuilder(transientTrackBuilderToken);
    ObjMan.Load( "Muons", muonObj );
*/


    KUCMSDisplacedVertexMini* displacedVertexObj = NULL;
    if( cfFlag("doSVModule") ){

		auto mergedSCToken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("displacedSCs"));
			
	    // Displaced Vertices
	    auto combinedMuonTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("combinedMuonTracks"));
	    auto muonEnhancedTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("sip2DMuonEnhancedTracks"));

	    displacedVertexObj = new KUCMSDisplacedVertexMini(iConfig);
	    displacedVertexObj->LoadMuonTracks(combinedMuonTracks);
	    displacedVertexObj->LoadMuonEnhancedTracksToken(muonEnhancedTracks);
	    displacedVertexObj->LoadPrimaryVertex(vertexToken);
	    displacedVertexObj->LoadTTrackBuilder(transientTrackBuilderToken);
	    displacedVertexObj->LoadMagneticField(magneticFieldToken);
	    displacedVertexObj->LoadAssociationParameters(trackAssocParameters);
	    displacedVertexObj->LoadMergedSCs(mergedSCToken);
  	
		ObjMan.Load( "DisplacedVertex", displacedVertexObj );

		electronsObj->LoadDisplacedVertexObject( displacedVertexObj );
		//muonObj->LoadDisplacedVertexObject( displacedVertexObj );

    }//<<>>if( cfFlag("doSVModule") )
	if( not cfFlag("doSVModule") ) geVar.set("nDisSVs",0.f);

	//Photons
    auto photonsObj = new KUCMSPhotonObjectMini( iConfig );
    auto photonToken = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("gedPhotons"));
    auto ootPhotonToken = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("ootPhotons"));
    photonsObj->LoadPhotonTokens( photonToken, ootPhotonToken );
    photonsObj->LoadRecHitObject( recHitsObj );
    photonsObj->LoadElectronObject( electronsObj );
    ObjMan.Load( "Photons", photonsObj );

	//Jets
    auto ak4jetObj = new KUCMSAK4JetObjectMini( iConfig );
    auto jetsToken = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
    ak4jetObj->LoadAK4JetTokens( jetsToken );
    ak4jetObj->LoadPfcandTokens( pfcandToken );
    ak4jetObj->LoadRecHitObject( recHitsObj );
    ak4jetObj->LoadPhotonObject( photonsObj );
    ak4jetObj->LoadElectronObject( electronsObj );
    ObjMan.Load( "JetsAK4", ak4jetObj );

	//MET
    auto pfmetObj = new KUCMSPFMetObjectMini( iConfig );
    auto pfmetToken = consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets")); 
    pfmetObj->LoadPFMetTokens( pfmetToken );
    pfmetObj->LoadPhotonObject( photonsObj );
    ObjMan.Load( "PFMet", pfmetObj );

	if( NTDEBUG ) std::cout << "Loading Final Objects" << std::endl;
    ObjMan.Load( "ECALRecHits", recHitsObj );// loaded last to process feedback from other objects

    if( cfFlag("hasGenInfo") ){

        if( NTDEBUG ) std::cout << "Create & Load Gen Information" << std::endl;

        auto genObjs = new KUCMSGenObjectMini( iConfig );
        auto genEvtInfoToken = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("genEvt"));
        genObjs->LoadGenEvtInfoTokens( genEvtInfoToken );
        auto gent0Token = consumes<float>(iConfig.getParameter<edm::InputTag>("gent0"));
        genObjs->LoadGenT0Tokens( gent0Token );
        auto genxyz0Token = consumes<Point3D>(iConfig.getParameter<edm::InputTag>("genxyz0"));
        genObjs->LoadGenXYZ0Tokens( genxyz0Token );
        auto pileupInfosToken = consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileups"));
        genObjs->LoadGenPUInfosTokens( pileupInfosToken );
        auto genPartToken = consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParticles"));
        genObjs->LoadgTokens( genPartToken );
        auto genJetsToken = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
        genObjs->LoadGenJetsTokens( genJetsToken );
		auto genPrunToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartPruned"));
        genObjs->LoadgpTokens( genPrunToken );

        // Load gen object into other collections
	    electronsObj->LoadGenObject( genObjs );
		//muonObj->LoadGenObject( genObjs );
        photonsObj->LoadGenObject( genObjs );
        ak4jetObj->LoadGenObject( genObjs );
		//if( cfFlag("doSVModule") || cfFlag("doDisEleModule") ){ ecalTracksObj->LoadGenObject( genObjs ); }
		if( cfFlag("doSVModule") ){ displacedVertexObj->LoadGenParticlesToken(genPartToken); }
	    //if( cfFlag("doDisEleModule") ){	displacedElectronObj->LoadGenObject( genObjs ); }

        // Load gen object into objman last, should be no dependence with other objects
        ObjMan.Load( "GenObjects", genObjs );

    }//<<>>if( hasGenInfo )

// ---------------------------------------------------------------------------------
}//>>>>KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ define this as a plug-in --------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(KUCMSNtupilizerMini);

