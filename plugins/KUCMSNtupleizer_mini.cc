// -*- C++ -*-
//
// Package:    KUCMSNtupleizer
// Class:      KUCMSNtupilizerMini
//
/**\class KUCMSNtupilizerMini KUCMSNtupleizer_mini.cc

MiniAOD analyzer that assembles event, detector-object, optional
displaced-vertex, and generator information into ntuple objects.
The constructor establishes CMSSW tokens and the dependency order
used by ObjectManager during event processing.

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

#include "KUCMSEventInfo_mini.hh"
#include "KUCMSEcalRechit_mini.hh"
#include "KUCMSPhoton_mini.hh"
#include "KUCMSAK4Jet_mini.hh"
#include "KUCMSPFMet_mini.hh"
#include "KUCMSElectron_mini.hh"
#include "KUCMSMuon_mini.hh"
#include "KUCMSGenObjects_mini.hh"
#include "KUCMSDisplacedVertex_mini.hh"
#include "KUCMSTrack_mini.hh"

using namespace std;

//#define NTDEBUG true
#define NTDEBUG false

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Constructor: read configuration flags, create input tokens, and wire the
// object modules together in the order required by their data dependencies.
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

  // Read optional feature flags from the cmsRun configuration.  Defaults
  // preserve the historical MiniAOD behavior when a flag is omitted.

  cfFlag.set( "hasGenInfo", iConfig.existsAs<bool>("hasGenInfo") ? iConfig.getParameter<bool>("hasGenInfo") : true );
  cfFlag.set( "doSVModule", iConfig.existsAs<bool>("doSVModule") ? iConfig.getParameter<bool>("doSVModule") : true );
  cfFlag.set( "doDisEleModule", iConfig.existsAs<bool>("doDisEleModule") ? iConfig.getParameter<bool>("doDisEleModule") : false );
  cfFlag.set( "doECALTrackOnly", iConfig.existsAs<bool>("doECALTrackOnly") ? iConfig.getParameter<bool>("doECALTrackOnly") : false );

  // Create tokens for the event collections and EventSetup records consumed
  // by the object modules below.

  //if( NTDEBUG )
  std::cout << "In constructor for KUCMSNtupilizer - tag and tokens" << std::endl;

  // Configure the event-level selection helper, if a selection tag was supplied.
  std::string filterSelect(iConfig.existsAs<std::string>("fltrSelection")? iConfig.getParameter<std::string>("fltrSelection") : "None" );
  ntupleFilter.setEventSelectionTag( filterSelect );

  if( NTDEBUG ) std::cout << "Create Object Classes" << std::endl;
  // ObjectManager stores modules by name; its name ordering determines the
  // processing order, so dependent objects are named and loaded accordingly.


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

  // Track inputs and shared propagation/association services.
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> pfcandToken = consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("pfcandidates"));
  auto parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  TrackAssociatorParameters trackAssocParameters;
  edm::ConsumesCollector iC = consumesCollector();
  trackAssocParameters.loadParameters(parameters, iC);
  auto transientTrackBuilderToken = esConsumes<TransientTrackBuilder,TransientTrackRecord>(edm::ESInputTag("","TransientTrackBuilder"));
  auto magneticFieldToken = esConsumes<MagneticField, IdealMagneticFieldRecord>();
  auto ogGeneralTracksToken = consumes<edm::View<reco::Track>>(iConfig.getParameter<edm::InputTag>("ogGeneralTracks"));
  auto ogGsfTracksToken = consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("ogGsfTracks"));

  if( cfFlag("doSVModule") ){
    //Track General + GSF
    auto TrackObj = new KUCMSTrackObjectMini( iConfig );
    TrackObj->LoadPfcandTokens( pfcandToken );
    TrackObj->LoadAssociationParameters(trackAssocParameters);
    TrackObj->LoadTTrackBuilder(transientTrackBuilderToken);
    TrackObj->LoadGeneralTrackTokens(ogGeneralTracksToken);
    TrackObj->LoadGsfTrackTokens(ogGsfTracksToken);
    TrackObj->LoadMagneticField(magneticFieldToken);
    TrackObj->LoadBeamSpotTokens( beamLineToken );
    TrackObj->LoadVertexTokens( vertexToken );
    TrackObj->LoadRecHitObject( recHitsObj );
    ObjMan.Load( "Tracks", TrackObj );
  }//<<>>if( cfFlag("doSVModule") )

  //Electrons
  auto electronsObj = new KUCMSElectronObjectMini( iConfig );
  auto rhoToken = consumes<double>(iConfig.getParameter<edm::InputTag>("rho"));
  electronsObj->LoadRhoTokens( rhoToken );
  auto electronToken = consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"));
  electronsObj->LoadElectronTokens( electronToken );
  // Conversion-token wiring is intentionally disabled; the current producer
  // configuration does not provide a supported conversion collection here.
  //auto conversionsToken = consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"));
  //electronsObj->LoadConversionTokens( conversionsToken );
  // Other optional electron associations remain disabled until their input
  // contracts are enabled in the producer configuration.
  electronsObj->LoadTTrackBuilder(transientTrackBuilderToken);
  electronsObj->LoadBeamSpotTokens( beamLineToken );
  electronsObj->LoadRecHitObject( recHitsObj );
  ObjMan.Load( "Electrons", electronsObj );

  //Muons
  auto muonObj = new KUCMSMuonObjectMini( iConfig );
  auto muonToken = consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"));
  muonObj->LoadVertexTokens( vertexToken );
  muonObj->LoadMuonTokens( muonToken );
  muonObj->LoadTTrackBuilder(transientTrackBuilderToken);
  ObjMan.Load( "Muons", muonObj );

  KUCMSDisplacedVertexMini* displacedVertexObj = NULL;
  if( cfFlag("doSVModule") ){

    auto mergedSCToken = consumes<reco::SuperClusterCollection>(iConfig.getParameter<edm::InputTag>("displacedSCs"));

    // Displaced Vertices
    auto muonEnhancedTracks = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("sip2DMuonEnhancedTracks"));
    auto leptonicSVsToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("leptonicSVs"));
    auto hadronicSVsToken = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("hadronicSVs"));

    displacedVertexObj = new KUCMSDisplacedVertexMini(iConfig);
    displacedVertexObj->LoadMuonToken(muonToken);
    displacedVertexObj->LoadMuonEnhancedTracksToken(muonEnhancedTracks);
    displacedVertexObj->LoadPrimaryVertex(vertexToken);
    displacedVertexObj->LoadTTrackBuilder(transientTrackBuilderToken);
    displacedVertexObj->LoadMagneticField(magneticFieldToken);
    displacedVertexObj->LoadAssociationParameters(trackAssocParameters);
    displacedVertexObj->LoadMergedSCs(mergedSCToken);
    displacedVertexObj->LoadLeptonicSVsToken(leptonicSVsToken);
    displacedVertexObj->LoadHadronicSVsToken(hadronicSVsToken);
    displacedVertexObj->LoadRecHitObject(recHitsObj);

    ObjMan.Load( "DisplacedVertex", displacedVertexObj );

    electronsObj->LoadDisplacedVertexObject( displacedVertexObj );
    muonObj->LoadDisplacedVertexObject( displacedVertexObj );

  }//<<>>if( cfFlag("doSVModule") )
  if( not cfFlag("doSVModule") ) geVar.set("nDisSVs",0.f);

  // Photon object: depends on ECAL rechits and the electron object.
  auto photonsObj = new KUCMSPhotonObjectMini( iConfig );
  auto photonToken = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("gedPhotons"));
  auto ootPhotonToken = consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("ootPhotons"));
  photonsObj->LoadPhotonTokens( photonToken, ootPhotonToken );
  photonsObj->LoadRecHitObject( recHitsObj );
  photonsObj->LoadElectronObject( electronsObj );
  ObjMan.Load( "Photons", photonsObj );

  // AK4 jets: depends on particle-flow candidates and selected objects.
  auto ak4jetObj = new KUCMSAK4JetObjectMini( iConfig );
  auto jetsToken = consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"));
  ak4jetObj->LoadAK4JetTokens( jetsToken );
  ak4jetObj->LoadPfcandTokens( pfcandToken );
  ak4jetObj->LoadRecHitObject( recHitsObj );
  ak4jetObj->LoadPhotonObject( photonsObj );
  ak4jetObj->LoadElectronObject( electronsObj );
  ObjMan.Load( "JetsAK4", ak4jetObj );

  // PF MET: retain the photon association used by the output object.
  auto pfmetObj = new KUCMSPFMetObjectMini( iConfig );
  auto pfmetToken = consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"));
  pfmetObj->LoadPFMetTokens( pfmetToken );
  pfmetObj->LoadPhotonObject( photonsObj );
  ObjMan.Load( "PFMet", pfmetObj );

  if( NTDEBUG ) std::cout << "Loading Final Objects" << std::endl;
  // Load rechits last so the module can receive feedback from the objects
  // that depend on its shared rechit state.
  ObjMan.Load( "ECALRecHits", recHitsObj );

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
    auto genJetsToken = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genjets"));
    genObjs->LoadGenJetsTokens( genJetsToken );
    auto genPrunToken = consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genPartPruned"));
    genObjs->LoadgpTokens( genPrunToken );

    // Make generator information available to reconstructed-object modules
    // before registering the generator module itself.

    electronsObj->LoadGenObject( genObjs );
    muonObj->LoadGenObject( genObjs );
    photonsObj->LoadGenObject( genObjs );
    ak4jetObj->LoadGenObject( genObjs );

    if( cfFlag("doSVModule") ){
      displacedVertexObj->LoadGenParticlesToken(genPartToken);
      displacedVertexObj->LoadPrunedGenParticlesToken(genPrunToken);
    }//<<>>if( cfFlag( ) )

    // Generator information is registered last because no later module
    // should be required to populate it.
    ObjMan.Load( "GenObjects", genObjs );

  }//<<>>if( hasGenInfo )

  // ---------------------------------------------------------------------------------
}//>>>>KUCMSNtupilizerMini::KUCMSNtupilizerMini(const edm::ParameterSet& iConfig)

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ define this as a plug-in --------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

DEFINE_FWK_MODULE(KUCMSNtupilizerMini);
