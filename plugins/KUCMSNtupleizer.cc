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
//#include "KUCMSHelperFunctions.hh"
//#include "KUCMSRootHelperFunctions.hh"

#include "KUCMSEventInfo.hh"
#include "KUCMSEcalRechit.hh"
#include "KUCMSPhoton.hh"
#include "KUCMSAK4Jet.hh"
#include "KUCMSPFMet.hh"
#include "KUCMSECALTracks.hh"
#include "KUCMSElectron.hh"
#include "KUCMSMuon.hh"
#include "KUCMSDisplacedElectron.hh"
#include "KUCMSGenObjects.hh"
#include "KUCMSDisplacedVertex.hh"

using namespace std;

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
    cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );
    cfFlag.set( "doDisEleModule", iConfig.existsAs<bool>("doDisEleModule") ? iConfig.getParameter<bool>("doDisEleModule") : false );
    //cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    //cfFlag.set( "motherChase", iConfig.existsAs<bool>("doGenMotherChase") ? iConfig.getParameter<bool>("doGenMotherChase") : false );
    //cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );

    // -- consume tags ------------------------------------------------------------
    // creats "token" for all collections that we wish to process
    
    //if( DEBUG ) 
    std::cout << "In constructor for KUCMSNtupilizer - tag and tokens" << std::endl;

	// -----  set event skim selection ----------------------------------------------------------------
	std::string filterSelect( iConfig.existsAs<std::string>("fltrSelection") ? iConfig.getParameter<std::string>("fltrSelection") : "None" );
	ntupleFilter.setEventSelectionTag( filterSelect );

    // Triggers
    //auto triggerObjectsToken = consumes<std::vector<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"));

	if( DEBUG ) std::cout << "Create Object Classes" << std::endl;

	//Event Info
    auto eventInfoObj = new KUCMSEventInfoObject(  iConfig ); 
    auto vertexToken = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"));
    auto triggerResultsToken = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"));
    auto triggerEventToken = consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
    //auto rhoToken_ = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
    eventInfoObj->LoadVertexTokens( vertexToken );
    eventInfoObj->LoadTriggerTokens( triggerResultsToken, triggerEventToken );
    ObjMan.Load( "EventInfo", eventInfoObj );

    //Rechits ECAL
    auto recHitsObj = new KUCMSEcalRecHitObject( iConfig );
    auto rhEBtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEB"));
    auto rhEEtoken = consumes<recHitCol>(iConfig.getParameter<edm::InputTag>("recHitsEE"));
    recHitsObj->LoadRecHitTokens( rhEBtoken, rhEEtoken );
    auto sctoken  = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("superClusters"));
    auto ootsctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("ootSuperClusters"));
    auto othersctoken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("otherSuperClusters"));
    if( cfFlag("doSVModule") ) recHitsObj->LoadSCTokens( sctoken, ootsctoken ); else recHitsObj->LoadSCTokens( othersctoken, ootsctoken );
    recHitsObj->LoadSCTokens( othersctoken );
    auto ccltoken = consumes<std::vector<reco::CaloCluster>>(iConfig.getParameter<edm::InputTag>("caloClusters"));
    recHitsObj->LoadClusterTokens( ccltoken );
    auto beamLineToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    recHitsObj->LoadBeamSpotTokens( beamLineToken );
    //ObjMan.Load( "ECALRecHits", recHitsObj );// loaded last to process feedback from other objects

    //Electrons 
    auto electronsObj = new KUCMSElectronObject( iConfig );
    auto electronToken = consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("electrons"));
    electronsObj->LoadElectronTokens( electronToken );
    auto conversionsToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
    electronsObj->LoadConversionTokens( conversionsToken );
    //auto beamLineToken = consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    electronsObj->LoadBeamSpotTokens( beamLineToken );
    electronsObj->LoadRecHitObject( recHitsObj );
    ObjMan.Load( "Electrons", electronsObj );

	KUCMSECALTracks* ecalTracksObj;
	KUCMSDisplacedElectron* displacedElectronObj;
    KUCMSDisplacedVertex* displacedVertexObj;
	if( cfFlag("doSVModule") || cfFlag("doDisEleModule") ){

	    //ECAL Tracks
	    ecalTracksObj = new KUCMSECALTracks(iConfig);
	    auto ecalTracksToken = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("ecalTracks"));
	    auto generalTracksToken = consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("tracks"));
	    //auto gsfTracksToken = consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("gsfTracksSrc"));
	    auto parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
	    auto magneticFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
	auto transientTrackBuilderToken = esConsumes<TransientTrackBuilder,TransientTrackRecord>(edm::ESInputTag("","TransientTrackBuilder"));
		auto mergedSCToken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("displacedSCs"));
			
    
	    // Setup the track associator
	    TrackAssociatorParameters trackAssocParameters;
	    edm::ConsumesCollector iC = consumesCollector();
	    trackAssocParameters.loadParameters(parameters, iC);
	
	    ecalTracksObj->LoadECALTracksToken(ecalTracksToken);
	    ecalTracksObj->LoadGeneralTrackTokens(generalTracksToken);
	    //ecalTracksObj->LoadGsfTrackTokens(gsfTracksToken);
	    ecalTracksObj->LoadMergedSCs(mergedSCToken);
	    ecalTracksObj->LoadAssociationParameters(trackAssocParameters);
	    ecalTracksObj->LoadMagneticField(magneticFieldToken);
	    ecalTracksObj->LoadBeamSpot(beamLineToken);

        ObjMan.Load( "ECALTracks", ecalTracksObj );

		if( cfFlag("doSVModule") ){

	    	// Displaced Vertices
	    	auto combinedMuonTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("combinedMuonTracks"));
	    	auto muonEnhancedTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("sip2DMuonEnhancedTracks"));

	    	displacedVertexObj = new KUCMSDisplacedVertex(iConfig);
	    	displacedVertexObj->LoadMuonTracks(combinedMuonTracks);
	    	displacedVertexObj->LoadMuonEnhancedTracksToken(muonEnhancedTracks);
	    	displacedVertexObj->LoadPrimaryVertex(vertexToken);
	    	displacedVertexObj->LoadTTrackBuilder(transientTrackBuilderToken);
	    	displacedVertexObj->LoadMagneticField(magneticFieldToken);
	    	displacedVertexObj->LoadAssociationParameters(trackAssocParameters);
	    	displacedVertexObj->LoadMergedSCs(mergedSCToken);
  	
		ObjMan.Load( "DisplacedVertex", displacedVertexObj );
    	}//<<>>if( cfFlag("doSVModule") )

		if( cfFlag("doDisEleModule") ){

        	//Displaced Electrons
        	displacedElectronObj = new KUCMSDisplacedElectron(iConfig);
        	auto gsfTracksToken = consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("gsfTracksSrc"));
        	auto displacedElectronToken = consumes<reco::ElectronCollection>(iConfig.getParameter<edm::InputTag>("displacedElectrons"));
			auto displacedSCToken = consumes<edm::View<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("displacedSCs"));
			auto ootSuperClusterToken = consumes<edm::View<reco::SuperCluster>>(iConfig.getParameter<edm::InputTag>("ootSuperClusters"));

        	displacedElectronObj->LoadECALTracksToken(ecalTracksToken);
        	displacedElectronObj->LoadGeneralTrackTokens(generalTracksToken);
        	displacedElectronObj->LoadGsfTrackTokens(gsfTracksToken);
        	displacedElectronObj->LoadSuperClusterTokens(displacedSCToken);
        	displacedElectronObj->LoadOotSuperClusterTokens(ootSuperClusterToken);
        	displacedElectronObj->LoadDisplacedElectrons(displacedElectronToken);
        	displacedElectronObj->LoadAssociationParameters(trackAssocParameters);
        	displacedElectronObj->LoadMagneticField(magneticFieldToken);
        	displacedElectronObj->LoadTTrackBuilder(transientTrackBuilderToken);

        	ObjMan.Load( "DisplacedElectrons", displacedElectronObj );

		}//<<>>if( cfFlag("doDisEleModule") )

    }//<<>>if( cfFlag("doSVModule") )
	if( not cfFlag("doSVModule") ) geVar.set("nDisSVs",0.f);

	//Muons    
	auto muonObj = new KUCMSMuonObject( iConfig );
    auto muonToken = consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
    muonObj->LoadMuonTokens( muonToken );
    ObjMan.Load( "Muons", muonObj );


	//Photons
    auto photonsObj = new KUCMSPhotonObject( iConfig );
    auto photonToken = consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("gedPhotons"));
    auto ootPhotonToken = consumes<edm::View<reco::Photon>>(iConfig.getParameter<edm::InputTag>("ootPhotons"));
    photonsObj->LoadPhotonTokens( photonToken, ootPhotonToken );
    photonsObj->LoadRecHitObject( recHitsObj );
    photonsObj->LoadElectronObject( electronsObj );
    ObjMan.Load( "Photons", photonsObj );

	//Jets
    auto ak4jetObj = new KUCMSAK4JetObject( iConfig );
    auto jetsToken = consumes<std::vector<reco::PFJet>>(iConfig.getParameter<edm::InputTag>("jets"));
    ak4jetObj->LoadAK4JetTokens( jetsToken );
    edm::EDGetTokenT<edm::View<reco::Candidate>> pfcandToken = consumes<CandidateView>(iConfig.getParameter<edm::InputTag>("pfcandidates"));
    ak4jetObj->LoadPfcandTokens( pfcandToken );
    ak4jetObj->LoadRecHitObject( recHitsObj );
    ak4jetObj->LoadPhotonObject( photonsObj );
    ak4jetObj->LoadElectronObject( electronsObj );
    ObjMan.Load( "JetsAK4", ak4jetObj );

	//MET
    auto pfmetObj = new KUCMSPFMetObject( iConfig );
    auto pfmetToken = consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("mets")); 
    pfmetObj->LoadPFMetTokens( pfmetToken );
    pfmetObj->LoadPhotonObject( photonsObj );
    ObjMan.Load( "PFMet", pfmetObj );

	if( DEBUG ) std::cout << "Loading Final Objects" << std::endl;
    ObjMan.Load( "ECALRecHits", recHitsObj );// loaded last to process feedback from other objects

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
		if( cfFlag("doSVModule") || cfFlag("doDisEleModule") ){ ecalTracksObj->LoadGenObject( genObjs ); }
		if( cfFlag("doSVModule") ){ displacedVertexObj->LoadGenParticlesToken(genPartToken); }
	    if( cfFlag("doDisEleModule") ){	displacedElectronObj->LoadGenObject( genObjs ); }

        // Load gen object into objman last, should be no dependence with other objects
        ObjMan.Load( "GenObjects", genObjs );

    }//<<>>if( hasGenInfo )

// ---------------------------------------------------------------------------------
}//>>>>KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ define this as a plug-in --------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(KUCMSNtupilizer);

