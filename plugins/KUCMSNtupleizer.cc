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

#include "./KUCMSNtupilizer.hh"
#include "./KUCMSHelperFunctions.hh"
#include "./KUCMSRootHelperFunctions.hh"

using namespace std;

//#define DEBUG true
#define DEBUG false

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// constructors and destructor
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig) :

	// -- declare tags and flags ----------------------------------------------------------
	// accuire input values form python config master file ( the one we use cmsRun with interactivilly )

	// flags
	hasGenInfo(iConfig.existsAs<bool>("hasGenInfo")  ? iConfig.getParameter<bool>("hasGenInfo")  : true),

	// triggers
	//triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
	//triggerObjectsTag(iConfig.getParameter<edm::InputTag>("triggerObjects")),

	// tracks
	tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),
	
	// pfcands
	pfcandTag(iConfig.getParameter<edm::InputTag>("pfcandidates")),
    pfCanTag(iConfig.getParameter<edm::InputTag>("particleflow")),
	pfCanPhoMapTag(iConfig.getParameter<edm::InputTag>("pfcanphomap")),
    pfCanOOTPhoMapTag(iConfig.getParameter<edm::InputTag>("pfcanootphomap")),
    pfCanEleMapTag(iConfig.getParameter<edm::InputTag>("pfcanelemap")),
	
	// vertices
	verticesTag(iConfig.getParameter<edm::InputTag>("vertices")),

	// rho
	//rhoTag(iConfig.getParameter<edm::InputTag>("rho")),

	// mets
	pfmetsTag(iConfig.getParameter<edm::InputTag>("mets")),  

	// supercluster
	superClusterCollectionTag(iConfig.getParameter<edm::InputTag>("superClusters")),
	ootSuperClusterCollectionTag(iConfig.getParameter<edm::InputTag>("ootSuperClusters")),

	// caloclusters
	caloClusterTag(iConfig.getParameter<edm::InputTag>("caloClusters")),

	// jets
	jetsTag(iConfig.getParameter<edm::InputTag>("jets")), 

    // calojets
    caloJetsTag(iConfig.getParameter<edm::InputTag>("calojets")),

    // genJets
    genJetsTag(iConfig.getParameter<edm::InputTag>("genjets")),

	// electrons
	electronsTag(iConfig.getParameter<edm::InputTag>("electrons")),  

	// muons
	muonsTag(iConfig.getParameter<edm::InputTag>("muons")),  

	// recHits
	recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),  
	recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

	// gedphotons
	gedPhotonsTag(iConfig.getParameter<edm::InputTag>("gedPhotons")),
    phoCBIDLooseMapTag(iConfig.getParameter<edm::InputTag>("phoCBIDLooseMap")),

	// ootPhotons
	ootPhotonsTag(iConfig.getParameter<edm::InputTag>("ootPhotons")),

	// pfcand ref
	//reco2pfTag(iConfig.getParameter<edm::InputTag>("recoToPFMap")),
	//reco2pf_(mayConsume<edm::ValueMap<std::vector<reco::PFCandidateRef> > >(iConfig.getParameter<edm::InputTag>("recoToPFMap"))),

	// genParticles
	genEvtInfoTag(iConfig.getParameter<edm::InputTag>("genEvt")),
  	gent0Tag(iConfig.getParameter<edm::InputTag>("gent0")),
  	genxyz0Tag(iConfig.getParameter<edm::InputTag>("genxyz0")),
  	pileupInfosTag(iConfig.getParameter<edm::InputTag>("pileups")),
  	genParticlesTag(iConfig.getParameter<edm::InputTag>("genParticles"))

	// -- end of tag declarations ---------------------------------------
{ //<<<< KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig) :

	usesResource();
	usesResource("TFileService");

	// -- consume tags ------------------------------------------------------------
	// creats "token" for all collections that we wish to process
	
	if( DEBUG ) std::cout << "In constructor for KUCMSNtupilizer - tag and tokens" << std::endl;

	// Triggers
	//triggerResultsToken_ 	= consumes<edm::TriggerResults>(triggerResultsTag);
	//triggerObjectsToken_ 	= consumes<std::vector<pat::TriggerObjectStandAlone>>(triggerObjectsTag);

	// tracks 
	tracksToken_				= consumes<std::vector<reco::Track>>(tracksTag);

	// genparticles
	// genpart_token_              = consumes<std::vector<pat::PackedGenParticle>>(genpartTag);

	// pfcandidates
	pfCan_token_        		= consumes<std::vector<reco::PFCandidate>>(pfCanTag);
	pfCanPhotonMap_token_		= consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanPhoMapTag);
    pfCanOOTPhotonMap_token_    = consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanOOTPhoMapTag);
    pfCanElectronMap_token_     = consumes<edm::ValueMap<edm::Ptr<reco::PFCandidate>>>(pfCanEleMapTag);

	// pfcandidates view
    pfcand_token_               = consumes<CandidateView>(pfcandTag);

	// vertices
	verticesToken_				= consumes<std::vector<reco::Vertex>>(verticesTag);

	// rho
	//rhoToken_					= consumes<double>(rhoTag);

	// mets
	pfmetsToken_				= consumes<std::vector<reco::PFMET>>(pfmetsTag);

	// supercluster
	scToken_             		= consumes<reco::SuperClusterCollection>(superClusterCollectionTag);
	ootScToken_          		= consumes<reco::SuperClusterCollection>(ootSuperClusterCollectionTag); 

	// caloClusters
	ccToken_			     	= consumes<std::vector<reco::CaloCluster>>(caloClusterTag);	

	// jets
	jetsToken_					= consumes<std::vector<reco::PFJet>>(jetsTag);
    caloJetsToken_              = consumes<std::vector<reco::CaloJet>>(caloJetsTag);

	// leptons
	electronsToken_				= consumes<std::vector<reco::GsfElectron>>(electronsTag);

	// muons
	muonsToken_					= consumes<std::vector<reco::Muon>>(muonsTag);

	// rechits
	recHitsEBToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
	recHitsEEToken_				= consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

	// photons
	gedPhotonsToken_ 			= consumes<std::vector<reco::Photon>>(gedPhotonsTag);
	phoCBIDLooseMapToken_		= consumes<edm::ValueMap<bool>>(phoCBIDLooseMapTag); 
	ootPhotonsToken_ 			= consumes<std::vector<reco::Photon>>(ootPhotonsTag);

	// pfref
	//reco2pf_					= consumes<edm::ValueMap<std::vector<reco::PFCandidateRef>>>(reco2pfTag);

	//genParticles & genJets
	if( hasGenInfo ){

		genJetsToken_ 			= consumes<std::vector<reco::GenJet>>(genJetsTag);
		genEvtInfoToken_   		= consumes<GenEventInfoProduct>(genEvtInfoTag);
    	gent0Token_        		= consumes<float>(gent0Tag);
    	genxyz0Token_      		= consumes<Point3D>(genxyz0Tag);
    	pileupInfosToken_  		= consumes<std::vector<PileupSummaryInfo>>(pileupInfosTag);
    	genParticlesToken_ 		= consumes<std::vector<reco::GenParticle>>(genParticlesTag);

	}//<<>>if( hasGenInfo )

// ---------------------------------------------------------------------------------
}//>>>>KUCMSNtupilizer::KUCMSNtupilizer(const edm::ParameterSet& iConfig)

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

	if( DEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

	// TRIGGER
	//iEvent.getByToken(triggerResultsToken_,triggerResults_);
	//iEvent.getByToken(triggerObjectsToken_,triggerObjects_);

	// TRACKS
	iEvent.getByToken(tracksToken_, tracks_);

	// PFCANDIDATES
	iEvent.getByToken(pfcand_token_, pfcands_);
	iEvent.getByToken(pfCan_token_, pfCans_);
	iEvent.getByToken(pfCanPhotonMap_token_, pfCanPhotonMap_);
    iEvent.getByToken(pfCanOOTPhotonMap_token_, pfCanOOTPhotonMap_);
    iEvent.getByToken(pfCanElectronMap_token_, pfCanElectronMap_);

	// VERTICES
	iEvent.getByToken(verticesToken_, vertices_);

	// RHO
	//iEvent.getByToken(rhoToken_, rho_);

	// METS
	iEvent.getByToken(pfmetsToken_, pfmets_);

	// SUPERCLUSTERS
	iEvent.getByToken(scToken_, superCluster_);  
	iEvent.getByToken(ootScToken_, ootSuperCluster_);

	// CALOCLUSTERS
	iEvent.getByToken(ccToken_, caloCluster_);

	// JETS
	iEvent.getByToken(jetsToken_, jets_);

    // CALOJETS
    iEvent.getByToken(caloJetsToken_, caloJets_);

	// GENPARTICLES & GENJETS
	if( hasGenInfo ){

		iEvent.getByToken(genEvtInfoToken_, genEvtInfo_);
    	iEvent.getByToken(gent0Token_, gent0_);
    	iEvent.getByToken(genxyz0Token_, genxyz0_);
    	iEvent.getByToken(pileupInfosToken_, pileupInfos_);
    	iEvent.getByToken(genParticlesToken_, genParticles_);
    	iEvent.getByToken(genJetsToken_, genJets_);

	}//<<>>if( hasGenInfo )

	// ELECTRONS
	iEvent.getByToken(electronsToken_, electrons_);

	// PHOTONS
	iEvent.getByToken(gedPhotonsToken_, gedPhotons_);
    iEvent.getByToken(phoCBIDLooseMapToken_, phoCBIDLooseMap_);
	iEvent.getByToken(ootPhotonsToken_, ootPhotons_);

	// MUONS
    iEvent.getByToken(muonsToken_, muons_);

	// ECAL RECHITS
	iEvent.getByToken(recHitsEBToken_, recHitsEB_);
	iEvent.getByToken(recHitsEEToken_, recHitsEE_);

	// GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
	iSetup.get<CaloGeometryRecord>().get(caloGeo_); 
	barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
	endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap); 
	ecalGeometry = caloGeo_.product();

    // CaloTopology
    iSetup.get<CaloTopologyRecord>().get(caloTopo_);
	barrelTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalBarrel);
    endcapTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalEndcap);
	ecalTopology = caloTopo_.product();

  	// Laser constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalLaserDbRecord>().get(laserH);
	evTime = iEvent.time();

  	// Intercalibration constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalIntercalibConstantsRcd>().get(interCalibH);
	interCalibMap = &interCalibH->getMap();  

  	// ADCToGeV : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
  	iSetup.get<EcalADCToGeVConstantRcd>().get(adcToGeVH);
	adcToGeVEB = adcToGeVH->getEBValue();
	adcToGeVEE = adcToGeVH->getEEValue();

  	// Pedestals : https://github.com/ferriff/usercode/blob/master/DBDump/plugins/DBDump.cc
  	iSetup.get<EcalPedestalsRcd>().get(pedestalsH);

	// PFCand refs
	//iEvent.getByToken(reco2pf_, reco2pf);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////
    // ---  Collecting objests ( preprocessing object pruning ) ---------------------------------------
    // -------------------------------------------------------------------------------------------------
	// -- Process Event  ---------------------------------------    
	// ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
	// ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 
	//////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << "***********************************************************************************************" << std::endl;
    if( DEBUG ) std::cout << "Processing event: " << event << " in run: " << run << " and lumiblock: " << lumi << std::endl;

	// declare local vector for filtered set of collections to process

	std::vector<reco::PFJet> 		fjets;
  	std::vector<int> 				fjetsID;
  	std::vector<reco::CaloCluster> 	fbclusts;
  	std::vector<EcalRecHit> 		frechits;
  	std::vector<reco::Photon> 		fphotons;
  	//std::vector<reco::Photon> 		footphotons;
  	std::vector<reco::GsfElectron> 	felectrons;
  	std::vector<reco::CaloJet> 		fcalojets;
	std::vector<reco::Muon> 		fmuons;
	std::vector<reco::PFMET> 		fpfmet;
	std::vector<reco::GenParticle>  fgenparts;
  	std::vector<int>				fgenpartllp;
	std::vector<reco::GenJet>		fgenjets;

	// Selection criteria declartation for filtering -------------------------------
	// ------  move these to config input values via tag system --------------------------

    float   jetPTmin        = 1.0;// for energy/time comp
	int 	jetIDmin		= 2; //3;
	float 	jetETAmax		= 1.5; //1.5;

	float 	minRHE			= 0.5;

    float   minMuonE       	= 1.0;

   	float 	minMetE			= 1.0;

    float   minCaloCltrE    = 1.0;
    float   minCaloJetE     = 1.0;

    float   minPhoE         = 1.0;
    float 	phoMinPt 		= 5.0;
    float 	phoMinSeedTime 	= -25.0;
    string  phoMvaWp80("mvaPhoID-RunIIFall17-v1-wp80");//2018
    //string phoMvaWp80("mvaPhoID-RunIIFall17-v2-wp80");//r3
    string  phoCutTight("cutBasedPhotonID-Fall17-94X-V2-tight");//2022
    string  phoCutLoose("cutBasedPhotonID-Fall17-94X-V2-loose");//2022
    //string phoCutLoose("cutBasedPhotonID-Fall17-94X-V1-loose");//2018i

    float   minEleE         = 1.0;
    string eleMvaWp80("mvaEleID-Fall17-noIso-V2-wp80");
    string eleMvaWpLoose("mvaEleID-Fall17-noIso-V2-wpLoose");//r3
    string eleCutTight("cutBasedElectronID-Fall17-94X-V2-tight");
    string eleCutLoose("cutBasedElectronID-Fall17-94X-V2-loose");//2022
    //string eleCutLoose("cutBasedElectronID-Fall17-94X-V1-loose");//2018

    float   minGenE         = 0.0;

	bool 	onlyEB 			= true;

	// -----------------------------------------------------------------------------------------------
	// ---  Collecting objests ( preprocessing object pruning ) ---------------------------------------
	// -------------------------------------------------------------------------------------------------
	
	// fill filtered collection vectors

    if( DEBUG ) std::cout << "Collecting ECAL RecHits" << std::endl;
    for( const auto recHit : *recHitsEB_ ){ if( recHit.energy() > minRHenr ) frechits.push_back(recHit); }
    for( const auto recHit : *recHitsEE_ ){ if( recHit.energy() > minRHenr ) frechits.push_back(recHit); }

    if( DEBUG ) std::cout << "Collecting Muons" << std::endl;
    for( const auto muon : *muons_ ){ if( muon.energy() > minMuonE ) fmuons.push_back(muon); }

    if( DEBUG ) std::cout << "Collecting Met" << std::endl;
    //for( const auto met : *pfmets_ ){ if( met.energy() > minMetE ) 
	fpfmet.push_back((*pfmets_)[0]);

	if( DEBUG ) std::cout << "Filter Jets" << std::endl;
	for(const auto &jet : *jets_ ){ // Filters jet collection & sorts by pt

		if (jet.pt() < jetPTmin) continue;
	  	if (std::abs(jet.eta()) > jetETAmax) continue;
	  	const auto jetID = getPFJetID(jet);
	  	if (jetID < jetIDmin) continue;
	  	// save the jets, and then store the ID
	  	fjets.emplace_back(jet);
	  	//fjets.back().addUserInt("jetID",jetID);
		fjetsID.push_back(jetID);	  

	  	std::sort(fjets.begin(),fjets.end(),sortByPt);

	}//<<>>for(const auto& jet : *jets_ )

	if( DEBUG ) std::cout << "Collecting GenJets" << std::endl;	
	for(const auto& genjet : *genJets_ ){ fgenjets.push_back(genjet); }

	//calojets
	if( DEBUG ) std::cout << "Collecting CaloJets" << std::endl;
	for(const auto &calojet : *caloJets_ ){ if( calojet.energy() > minCaloJetE ) fcalojets.push_back(calojet);  }

	//if( DEBUG ) std::cout << "Collecting Calo Clusters" << std::endl;
	//for( const auto &bclust : *caloCluster_ ){ if( bclust.energy() > minCaloCltrE ) fbclusts.push_back(bclust); }

    if( DEBUG ) std::cout << "Collecting Photons/OOTPhotons" << std::endl;

    phoExcluded.clear();
    phoIsOotPho.clear();
    for( const auto &ootPho : *ootPhotons_ ){
		auto passIdCut = ootPho.photonID(phoCutLoose);
        auto minPhoPt = ootPho.pt() < phoMinPt;
        auto phoSeedTime = getPhotonSeedTime(ootPho);
        auto timecut = phoSeedTime < phoMinSeedTime;
        auto energycut = ootPho.energy() < minPhoE;
		if( not passIdCut || minPhoPt || timecut || energycut ) continue;
        double minDr(0.5);
        double dRmatch(10.0);
        float matchpt(0);
        auto oEta = ootPho.eta();
        auto oPhi = ootPho.phi();
		auto oPt = ootPho.pt();
        for( const auto gedPho : *gedPhotons_ ){
            auto pEta = gedPho.eta();
            auto pPhi = gedPho.phi();
            auto pPt = gedPho.pt();
            dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = pPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
		fphotons.push_back(ootPho);
		phoIsOotPho.push_back(true);
        if( dRmatch < 0.1 && oPt < matchpt ) phoExcluded.push_back(true); 
		else phoExcluded.push_back(false);
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )
    for( const auto gedPho : *gedPhotons_ ){
        auto passIdCut = gedPho.photonID(phoCutLoose);
        auto ptCut = gedPho.pt() < phoMinPt;
        auto phoSeedTime = getPhotonSeedTime(gedPho);
        auto timeCut = phoSeedTime < phoMinSeedTime;
        auto energyCut = gedPho.energy() < minPhoE;
        if( not passIdCut || ptCut || timeCut || energyCut ) continue;
        double minDr(0.5);
        double dRmatch(10.0);
        float matchpt(0);
		auto pEta = gedPho.eta();
		auto pPhi = gedPho.phi();
		auto pPt = gedPho.pt();
        for( const auto &ootPho : *ootPhotons_ ){
			auto oEta = ootPho.eta();
			auto oPhi = ootPho.phi();
			auto oPt = ootPho.pt(); 
            dRmatch = deltaR( pEta, oEta, pPhi, oPhi );
            if( dRmatch < minDr ){ minDr = dRmatch; matchpt = oPt; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        fphotons.push_back(gedPho);
        phoIsOotPho.push_back(false);
        if( dRmatch < 0.1 && pPt < matchpt ) phoExcluded.push_back(true);
        else phoExcluded.push_back(false); 
    }//<<>>for( int io = 0; io < nOotPhotons; io++ )

    if( DEBUG ) std::cout << "Collecting Electrons" << std::endl;
    for( const auto electron : *electrons_ ){ 
		auto passIdCut = electron.electronID(eleCutLoose);
		auto passEnergyCut = electron.energy() > minEleE;
		if( passIdCut && passEnergyCut ) felectrons.push_back(electron); }
	}//<<>>for( const auto electron : *electrons_ )

	if( DEBUG ) std::cout << "Collecting Gen Particles" << std::endl;
    int nGenPart(0);
    nGenPart = 0;
	for( const auto genPart : *genParticles_ ){ 

		if( genPart.energy() < minGenE ) continue;
		nGenPart++;
		if( genPart.status() == 1 ){
			fgenparts.push_back(genPart);
			auto genPartID = llpGenChaseP(genPart,0);
			fgenpartllp.push_back(genPartID);
		}//<<>>if( genPart.isLastCopy() )

	}//<<>>for(const auto& genPart : *genParticles_ )
    if( DEBUG ) std::cout << "Finished Collecting Gen Particles" << std::endl;

    //-------------`-----------------------------------------------------------------------
	// ----   Object processing ----------------------------------------------------------
    //------------------------------------------------------------------------------------
	// call functions to process collections and fill tree varibles to be saved
	// varibles to be saved to ttree are declared in the header

	processEvent( iEvent );
	processVtx( vertices_ );
	processGenPart( fgenparts );
    processRecHits( frechits );
	processMet( fpfmet, fphotons );
    processClJet( fcalojets );
	processPhotons( fphotons, fgenparts );
    processElectrons( felectrons, fgenparts );
	processJets( fjets, fjetsID, fgenparts );

    //------------------------------------------------------------------------------------
    //---- Object processing Completed ----------------------------------------------------------
    //------------------------------------------------------------------------------------

	// -- Fill output trees ------------------------------------------
	if( DEBUG ) std::cout << "---------- Next Event -----" << std::endl;
	outTree->Fill();

	// -- EOFun ------------------------------------------------------
	//	 #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	//	 ESHandle<SetupData> pSetup;
	//	 iSetup.get<SetupRecord>().get(pSetup);
	//	 #endif

}//>>>>void KUCMSNtupilizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ beginJob/endJob methods called once each job just before/after starting event loop	------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void KUCMSNtupilizer::beginJob(){

	// Global Varibles
	nGoodJetEvents = 0;

	// Book output files and trees
	edm::Service<TFileService> fs;
	outTree = fs->make<TTree>("llpgtree","llpgtree");

	// Book //histograms ( if any )
	
	std::cout << "Histograms Booked" << std::endl;

	// Create output Tree branches -----------------------------
	// via set branch functions that call : outTree->Branch("",&); 
	// for the varibles to be saved in the tree for that particular collection

	setBranchesEvent();
	setBranchesVtx();
	setBranchesMet();
	setBranchesJets();
    setBranchesClusterJets();
	setBranchesPhotons();
	setBranchesElectrons();
	setBranchesGenParts();
	setBranchesRecHits();

}//>>>>void KUCMSNtupilizer::beginJob()


// ------------ method called once each job just after ending the event loop	------------
void KUCMSNtupilizer::endJob(){

}//>>>>void KUCMSNtupilizer::endJob()


// ------------ method fills 'descriptions' with the allowed parameters for the module	------------
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

//////////////////////////////////////////////////////////////////////////////////////////////////
//  --------------  collection processing and branch declaration functions ----------------------
//////////////////////////////////////////////////////////////////////////////////////////////////

//---------------------------------------------
void KUCMSNtupilizer::setBranchesEvent(){

    // Run, Lumi, Event info
    //    ulInt               event; // technically unsigned long long in Event.h...
    //    uInt                run, lumi;
    outTree->Branch("run", &run);
    outTree->Branch("lumi", &lumi);
    outTree->Branch("event", &event, "event/l");

}//<<>>setBranchesEvent()

void KUCMSNtupilizer::processEvent( const edm::Event& iEvent ){

    run   = iEvent.id().run();
    lumi  = iEvent.luminosityBlock();
    event = iEvent.id().event();

}//<<>>processEvent()

//----------------------------------------------
void KUCMSNtupilizer::setBranchesVtx(){

    //    int                 nVtx;
    //    float               vtxX, vtxY, vtxZ;
    //outTree->Branch("nVtx", &nVtx);
    outTree->Branch("vtxX", &vtxX);
    outTree->Branch("vtxY", &vtxY);
    outTree->Branch("vtxZ", &vtxZ);

}//<<>>setBranchesVtx()

void KUCMSNtupilizer::processVtx( edm::Handle<std::vector<reco::Vertex>> vertices ){

	// do not need to pass in collectioon handel -> class varible   --  Fix This
	//edm::Handle<std::vector<reco::Vertex>> vertices_
    const auto & primevtx = vertices->front();
    vtxX = primevtx.position().x();
    vtxY = primevtx.position().y();
    vtxZ = primevtx.position().z();

}//<<>>processVtx()

//-------------------------------------------------
void KUCMSNtupilizer::setBranchesMet(){

    //    float  metSumEt, metPt, metPx, metPy, metPhi;
    outTree->Branch("metSumEt", &metSumEt);
    outTree->Branch("metPt", &metPt);
    outTree->Branch("metPx", &metPx);
    outTree->Branch("metPy", &metPy);
    outTree->Branch("metPhi", &metPhi);
    outTree->Branch("metCSumEt", &metCSumEt);
    outTree->Branch("metCPx", &metCPx);
    outTree->Branch("metCPy", &metCPy);

}//<<>>void setBranchesMet()

void KUCMSNtupilizer::processMet( std::vector<reco::PFMET> fpfmet, std::vector<reco::Photon> fphotons ){

    if( DEBUG ) std::cout << "Processing MET" << std::endl;

    auto t1pfMET = fpfmet[0];
    metSumEt = t1pfMET.sumEt();
    metPt = t1pfMET.pt();
    metPx = t1pfMET.px();
    metPy = t1pfMET.py();
    metPhi = t1pfMET.phi();

    metCSumEt = metSumEt;
    metCPx = metPx;
    metCPy = metPy;

    // ???  need corrected values for photon Et and Pt
    int nphidx = fphotons.size();
    for( auto phidx = 0; phidx < nphidx; phidx++ ){

        auto photon = fphotons[phidx];
        if( phoIsOotPho[phidx] && not phoExcluded[phidx] ){
            metCPx -= ((photon.pt())*std::cos(photon.phi()));
            metCPy -= ((photon.pt())*std::sin(photon.phi()));
            metCSumEt -= photon.et();
        }//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )
        if( not phoIsOotPho[phidx] && phoExcluded[phidx] ){
            metCPx += ((photon.pt())*std::cos(photon.phi()));
            metCPy += ((photon.pt())*std::sin(photon.phi()));
            metCSumEt += photon.et();
        }//<<>>if( phoIsOotPho[phidx] && not phoExcluded[phidx] )   

    }//<<>>for( auto phidx = 0; phidx < fphotons.size(); phidx++ )

}//<<>>processMet( std::vector<reco::PFMET> fpfmet, std::vector<reco::Photon> fphotons )

//---------------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesJets(){

    //    uInt                nJets, nGoodDrJets, nGoodScJets, nGoodBcJets, nUnJets;
    //    float               jetHt;
    // xxx   uInt                nGoodJetEvents;
    outTree->Branch("jetHt", &jetHt);
    outTree->Branch("nJets", &nJets);
    outTree->Branch("nGoodDrJets", &nGoodDrJets);
    outTree->Branch("nGoodScJets", &nGoodScJets);
    outTree->Branch("nGoodBcJets", &nGoodBcJets);
    outTree->Branch("nUnJets", &nUnJets);

    // Jet info
    //    std::vector<float>  jetE, jetPt, jetPhi, jetEta;
    // xxx   std::vector<float>  jetEtaetaMmt, jetPhiphiMnt, jetEtaphiMnt, jetMaxD, jetConPtDis, jetConEtaPhiSprd, jetArea;
    // xxx   std::vector<uInt>   jetNCarry, jetNConst;
    //    std::vector<float>  jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF, jetNHM, jetCHM, jetC, jetPHE, jetPHEF;
    //    std::vector<float>  jetELE, jetELEF, jetMUE, jetCharge;
  	outTree->Branch("jetE", &jetE);
  	outTree->Branch("jetPt", &jetPt);
  	outTree->Branch("jetEta", &jetEta);
  	outTree->Branch("jetPhi", &jetPhi);
  	outTree->Branch("jetNHF", &jetNHF);
  	outTree->Branch("jetNEMF", &jetNEMF);  
  	outTree->Branch("jetCHF", &jetCHF);
  	outTree->Branch("jetCEMF", &jetCEMF);
  	outTree->Branch("jetMUF", &jetMUF);
  	outTree->Branch("jetNHM", &jetNHM);
  	outTree->Branch("jetCHM", &jetCHM);
 	outTree->Branch("jetPHM", &jetPHM);
  	outTree->Branch("jetELM", &jetELM);
    outTree->Branch("jetC", &jetC);
  	outTree->Branch("jetPHE", &jetPHE);
  	outTree->Branch("jetPHEF", &jetPHEF);
  	outTree->Branch("jetELE", &jetELE);
  	outTree->Branch("jetELEF", &jetELEF);
  	outTree->Branch("jetMUE", &jetMUE);
    outTree->Branch("jetCharge", &jetCharge);
    outTree->Branch("jetID", &jetID);
    //    std::vector<int>    jetID, njetKids, jetKidOfJet, njetSubs, njetRecHits, jetRecHitOfJet;
    //    std::vector<int>    jetKidPdgID, jetKidCharge, jetKid3Charge, jetPHM, jetELM;
    //    std::vector<uInt>   jetRecHitId;
    //    std::vector<bool>   jetKidLLP;
    //    std::vector<double> jetKidMass, jetKidVx, jetKidVy, jetKidVz;
    //    std::vector<float>  jetKidE, jetKidPt, jetKidPhi, jetKidEta, jetKidTime, jetKidMedTime;

    // xxx    std::vector<float>  jetPhMuTime, jetOOTPhMuTime, jetEleMuTime;

    //    std::vector<float>  jetSumEPFrac, jetEPEnergy, jetEMEnergy, jetEMEnrFrac, jetEPEnrFrac;

    outTree->Branch("jetSumEPFrac", &jetSumEPFrac);
    outTree->Branch("jetEPEnergy", &jetEPEnergy);
    outTree->Branch("jetEMEnergy", &jetEMEnergy);
    outTree->Branch("jetEMEnrFrac", &jetEMEnrFrac);
    outTree->Branch("jetEPEnrFrac", &jetEPEnrFrac);

    //    std::vector<float>  jetDrLeadEta, jetDrLeadPhi, jetDrLeadEnr;
    //    std::vector<uInt>   jetDrRhCnt;
    //    std::vector<float>  sJetDrRHEnergy, jetDrEMF;

    //    std::vector<uInt>   nJetScMatch, jetScRhCnt;
    //    std::vector<float>  sJetScEnergy, sJetScPhEnergy, sJetScRhEnergy, jetScEMF;

    //    std::vector<uInt>   jetBcTimesCnt, jetBcRhCnt, jetBcGrpCnt;
    //    std::vector<float>  jetBcSumRHEnr, jetBcEMFr;

    //    std::vector<float>  jetMuTime, jetTimeError, jetTimeRMS, jetMedTime, jetCMuTime, jetCMedTime;
    //    std::vector<float>  jetSCMuTime, jetSCMedTime, jetCSCMuTime, jetCSCMedTime, jetCBCMuTime, jetCBCMedTime;

    outTree->Branch("jetDrLeadEta", &jetDrLeadEta);
    outTree->Branch("jetDrLeadPhi", &jetDrLeadPhi);
    outTree->Branch("jetDrLeadEnr", &jetDrLeadEnr);

    outTree->Branch("sJetDrRHEnergy", &sJetDrRHEnergy);
    outTree->Branch("jetDrEMF", &jetDrEMF);
    outTree->Branch("jetDrRhCnt", &jetDrRhCnt);
    //outTree->Branch("jetDrRhIds", &jetDrRhIds);

    //    std::vector<float>  jetGenImpactAngle, jetGenTime, jetGenPt, jetGenEta, jetGenEnergy, jetGenEMFrac, jetGenDrMatch;
    //    std::vector<float>  jetGenTimeVar, jetGenTimeLLP, jetGenLLPPurity, jetGenNextBX, jetGenNKids, jetGenTOF;

    outTree->Branch("jetGenImpactAngle", &jetGenImpactAngle);
    outTree->Branch("jetGenTime", &jetGenTime);
    outTree->Branch("jetGenPt", &jetGenPt);
    outTree->Branch("jetGenEta", &jetGenEnergy);
    outTree->Branch("jetGenEnergy", &jetGenEnergy);
    outTree->Branch("jetGenEMFrac", &jetGenEMFrac);
    outTree->Branch("jetGenDrMatch", &jetGenDrMatch);
    outTree->Branch("jetGenTimeVar", &jetGenTimeVar);
    outTree->Branch("jetGenTimeLLP", &jetGenTimeLLP);
    outTree->Branch("jetGenLLPPurity", &jetGenLLPPurity);
    outTree->Branch("jetGenNextBX", &jetGenNextBX);
    outTree->Branch("jetGenNKids", &jetGenNKids);
    outTree->Branch("jetGenTOF", &jetGenTOF);

}//<<>> setBranchesJets(){

void KUCMSNtupilizer::processJets( std::vector<reco::PFJet> fjets, std::vector<int> fjetsID, std::vector<reco::GenJet> fgenjets ){

    nUnJets = (*jets_).size();
    nJets = fjets.size();
    if( DEBUG ) std::cout << "Init for Jet Loop with " << nJets << " jets"<< std::endl;

    nGoodDrJets = 0;
    nGoodScJets = 0;
    nGoodBcJets = 0;
    vector<uInt> noIdGrp(1,0);

    //set the number of leading jets to skim ( = nJets for all )
    //auto stJets = nJets; 

    jetHt = 0.0;
    jetE.clear();
    jetPt.clear();
    jetPhi.clear();
    jetEta.clear();
    jetEtaetaMmt.clear();
    jetPhiphiMnt.clear();
    jetEtaphiMnt.clear();
    jetMaxD.clear();
    jetConPtDis.clear();
    jetConEtaPhiSprd.clear();
    jetArea.clear();
    jetNCarry.clear();
    jetNConst.clear();

    jetID.clear();
    jetNHF.clear();
    jetNEMF.clear();
    jetCHF.clear();
    jetCEMF.clear();
    jetMUF.clear();
    jetNHM.clear();
    jetCHM.clear();

    jetMuTime.clear();
    jetTimeError.clear();
    jetTimeRMS.clear();
    jetMedTime.clear();
    jetSCMuTime.clear();
    jetSCMedTime.clear();
    jetCMuTime.clear();
    jetCMedTime.clear();
    jetCSCMuTime.clear();
    jetCSCMedTime.clear();
    jetCBCMuTime.clear();
    jetCBCMedTime.clear();
    jetPhMuTime.clear();
    jetOOTPhMuTime.clear();
    jetEleMuTime.clear();

    //jetDrRhIds.clear();

    jetPHM.clear();
    jetELM.clear();
    //jetC.clear();
    jetPHE.clear();
    jetPHEF.clear();
    jetELE.clear();
    jetELEF.clear();
    jetMUE.clear();
    //jetCharge.clear();

    njetRecHits.clear();
    jetRecHitOfJet.clear();
    jetRecHitId.clear();

    jetSumEPFrac.clear();
    jetEPEnergy.clear();
    jetEMEnergy.clear();
    jetEMEnrFrac.clear();
    jetEPEnrFrac.clear();

    jetDrLeadEta.clear();
    jetDrLeadPhi.clear();
    jetDrLeadEnr.clear();

    sJetDrRHEnergy.clear();
    jetDrEMF.clear();
    jetDrRhCnt.clear();

    jetGenImpactAngle.clear();
    jetGenTime.clear();
    jetGenPt.clear();
    jetGenEta.clear();
    jetGenEnergy.clear();
    jetGenEMFrac.clear();
    jetGenDrMatch.clear();
    jetGenTimeVar.clear();
    jetGenTimeLLP.clear();
    jetGenLLPPurity.clear();
    jetGenNextBX.clear();
    jetGenNKids.clear();
    jetGenTOF.clear();

	// **** const reco::TrackRefVector& associatedTracks() const;
	// **** const reco::PixelClusterTagInfo* tagInfoPixelCluster(const std::string& label = "") const; // from b tagging info methods
	// **** reco::PFCandidateFwdPtrVector const& pfCandidatesFwdPtr() const { return pfCandidatesFwdPtr_; }
    bool goodJetEvent(false);
    if( DEBUG ) std::cout << "-----------------------------------------------------------" << std::endl;
    if( DEBUG ) std::cout << "Starting Jet Loop for " << nJets << " jets " << std::endl;
    for ( uInt ijet(0); ijet < nJets; ijet++ ){

        const auto & jet = fjets[ijet];
        const auto jetid = fjetsID[ijet];
        if( DEBUG ) std::cout << " --- Proccesssing : " << jet  << std::endl;

        jetHt += jet.pt();

        jetE.push_back(jet.energy());
        jetPt.push_back(jet.pt());
        jetPhi.push_back(jet.phi());
        jetEta.push_back(jet.eta());

        jetEtaetaMmt.push_back(jet.etaetaMoment());
        jetPhiphiMnt.push_back(jet.phiphiMoment());
        jetEtaphiMnt.push_back(jet.etaphiMoment());
        jetMaxD.push_back(jet.maxDistance());
        jetConPtDis.push_back(jet.constituentPtDistribution());
        jetConEtaPhiSprd.push_back(jet.constituentEtaPhiSpread());
        jetArea.push_back(jet.jetArea());
        jetNCarry.push_back(jet.nCarrying(0.1));
        jetNConst.push_back(jet.nConstituents());

        jetID.push_back(jetid);
        jetNHF.push_back(jet.neutralHadronEnergyFraction());
        jetNEMF.push_back(jet.neutralEmEnergyFraction());
        jetCHF.push_back(jet.chargedHadronEnergyFraction());
        jetCEMF.push_back(jet.chargedEmEnergyFraction());
        jetMUF.push_back(jet.muonEnergyFraction());
        jetNHM.push_back(jet.neutralMultiplicity());
        jetCHM.push_back(jet.chargedMultiplicity());
        //jetCharge.push_back(jet.jetCharge());

        jetPHE.push_back(jet.photonEnergy());
        jetPHEF.push_back(jet.photonEnergyFraction());
        jetELE.push_back(jet.electronEnergy());
        jetELEF.push_back(jet.electronEnergyFraction());
        jetMUE.push_back(jet.muonEnergy());
        jetPHM.push_back(jet.photonMultiplicity());
        jetELM.push_back(jet.electronMultiplicity());

        const auto jetepafrac   = jet.photonEnergyFraction() + jet.electronEnergyFraction();
        const auto jetepe       = jet.photonEnergy() + jet.electronEnergy();
        const auto jeteme       = jet.chargedEmEnergy() + jet.neutralEmEnergy();
        const auto jetemfrac    = jeteme/jet.energy();
        const auto jetepfrac    = jetepe/jet.energy();

        jetSumEPFrac.push_back(jetepafrac);
        jetEPEnergy.push_back(jetepe);
        jetEMEnergy.push_back(jeteme);
        jetEMEnrFrac.push_back(jetemfrac);
        jetEPEnrFrac.push_back(jetepfrac);

        // dR matched RH group  -----------------------------------------------
        // --------------------------------------------------------------------

        if( DEBUG ) std::cout << "Getting jet dR rechit group" << std::endl;
        auto jetDrRhGroup = getRHGroup( jet.eta(), jet.phi(), deltaRminJet, 0.5 );
        auto jetDrRhIdsGroup = getRhGrpIDs( jetDrRhGroup );
        jetDrRhIds.push_back( jetDrRhIdsGroup );

        auto rhCount = jetDrRhGroup.size();
        auto sumdrrhe = getRhGrpEnr( jetDrRhGroup );
        auto dremf = sumdrrhe/jet.energy();
        sJetDrRHEnergy.push_back(sumdrrhe);
        jetDrEMF.push_back(dremf);
        jetDrRhCnt.push_back(rhCount);

        if( DEBUG ) std::cout << " - get jetDRtofTimes " << std::endl;
        auto tofTimes = getLeadTofRhTime( jetDrRhGroup, vtxX, vtxY, vtxZ );
        auto leadJetRh = getLeadRh( jetDrRhGroup );

        if( DEBUG ) std::cout << "Starting RecHit Loop" << std::endl;
        for ( uInt irhg = 0; irhg < rhCount; irhg++){

            jetRecHitOfJet.push_back(ijet);
            auto detid = (jetDrRhGroup[irhg]).detid();
            jetRecHitId.push_back(detid.rawId());

        }//<<>>for ( uInt irhg = 0; irhg < rhCount; irhg++)

        const auto leadJetRhId = leadJetRh.detid();
        const auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
        auto sc_eta = leadJetRhIdPos.eta();
        auto sc_phi = leadJetRhIdPos.phi();
        auto sc_enr = leadJetRh.energy();
        //  make jettime varible
        auto jetTimeStats = getTimeDistStats( tofTimes, jetDrRhGroup );
        auto jmutime = jetTimeStats[0];
        auto jterr = jetTimeStats[1];
        auto jtrms = jetTimeStats[4];
        auto jmedtime = jetTimeStats[2];
        //auto mederr = jetTimeStats[3];
        auto jcmutime = jetTimeStats[6];
        auto jcmedtime = jetTimeStats[10];

        jetMuTime.push_back(jmutime);
        jetTimeError.push_back(jterr);
        jetTimeRMS.push_back(jtrms);
        jetMedTime.push_back(jmedtime);
        jetCMuTime.push_back(jcmutime);
        jetCMedTime.push_back(jcmedtime);

        jetDrLeadEta.push_back(sc_eta);
        jetDrLeadPhi.push_back(sc_phi);
        jetDrLeadEnr.push_back(sc_enr);

        if( jcmutime > -28.9 ) nGoodDrJets++;

        // GenJet Info for MC  -------------------------------------------------------------------
        // ---------------------------------------------------------------------------------------

        if( DEBUG ) std::cout << "Getting jetGenParton Information" << std::endl;

        //const reco::GenParticle * jetGenParton(0);
        float genImpactAngle(7.0);
        float genTime(-99.9);
        float genPt(-1.0);
        float genEta(10.0);
        float genEnergy(-1.0);
        float genEMFrac(-1.0);
        float genDrMatch(-1.0);
        float genTimeVar(-1.0);
        float genTimeLLP(0.0);
        float genLLPPurity(-1.0);
        float genNextBX(-1.0);
        float genNKids(-1.0);
        float genTOF(-99.9);

        if( hasGenInfo ){

            float goodDr(0.1);
            int matchedIdx(-1);
            int index(0);
            for(const auto& genJet : fgenjets ){

                auto gjeta = genJet.eta();
                auto gjphi = genJet.phi();
                auto jtgjdr = std::sqrt(reco::deltaR2(gjeta, gjphi, jet.eta(), jet.phi() ));
                if( jtgjdr <= goodDr ){
                    goodDr = jtgjdr;
                    matchedIdx = index;
                }//<<>>if( jtgjdr <= goodDr )
                index++;

            }//<<>>for(const auto& genJet : fgenjets ) 

            if( matchedIdx >= 0 ){

                auto genJet = fgenjets[matchedIdx];
                if( DEBUG ) std::cout << " --- Jet-GenJet dR match : " << goodDr << std::endl;

                auto nSources = genJet.numberOfSourceCandidatePtrs();
                if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
                if( DEBUG ) std::cout << " - genJet srcs : " << nSources << " PV(" << vtxX << "," << vtxY << "," << vtxZ << ")" << std::endl;
                auto kids = genJet.daughterPtrVector();
				if( 
                auto leadJetRh = getLeadRh( jetDrRhGroup );
                auto leadJetRhId = leadJetRh.detid();
                auto leadJetRhIdPos = barrelGeometry->getGeometry(leadJetRhId)->getPosition();
                auto cx = leadJetRhIdPos.x();
                auto cy = leadJetRhIdPos.y();
                auto cz = leadJetRhIdPos.z();
                auto tofcor = hypo( cx, cy, cz )/SOL;
                //if( DEBUG ) kidChase( kids, vtxX, vtxY, vtxZ );
                auto genKidInfo = kidTOFChain( kids, cx, cy, cz );
                if( DEBUG ) std::cout << " - genJet GenTime noTOF : " << genKidInfo[0] << " rhPos: " << cx;
				if( DEBUG ) std::cout << "," << cy << "," << cz << std::endl;
                genEta = genJet.eta();
                if( genKidInfo[0] > 25.0 ) genTime = -28.0;
                else if( genKidInfo[0] > -25.0 ) genTime = genKidInfo[0]-tofcor;
                else genTime = -27.0;
                genImpactAngle = genKidInfo[1];
                if( DEBUG ) std::cout << " - genJet GenTime : " << genTime << " Angle: " << genImpactAngle << std::endl;
                genPt = genJet.pt();
                genEnergy = genJet.energy();
                genEMFrac = (genJet.chargedEmEnergy() + genJet.neutralEmEnergy())/genEnergy;
                genDrMatch = goodDr; //std::sqrt(reco::deltaR2(jet.eta(), jet.phi(), genJet.eta(), genJet.phi()));
                genTimeVar = genKidInfo[2];
                genNextBX = genKidInfo[3];
                genTimeLLP = genKidInfo[4];
                genLLPPurity = genKidInfo[5];
                genNKids = genKidInfo[6];
                genTOF = tofcor;
                if( DEBUG ) std::cout << " -- Energy : " << genEnergy << " Pt : " << genPt << " EMfrac : " << genEMFrac << std::endl;

            }//<<>>if( matchedIdx >= 0 ) 

            // load event level vectors for this jet with gen info
            jetGenImpactAngle.push_back(genImpactAngle);
            jetGenTime.push_back(genTime);
            jetGenPt.push_back(genPt);
            jetGenEta.push_back(genEta);
            jetGenEnergy.push_back(genEnergy);
            jetGenEMFrac.push_back(genEMFrac);
            jetGenDrMatch.push_back(genDrMatch);
            jetGenTimeVar.push_back(genTimeVar);
            jetGenTimeLLP.push_back(genTimeLLP);
            jetGenLLPPurity.push_back(genLLPPurity);
            jetGenNextBX.push_back(genNextBX);
            jetGenNKids.push_back(genNKids);
            jetGenTOF.push_back(genTOF);

            if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;

        }//<<>>if( hasGenInfo )

        if( DEBUG ) std::cout << "Next Jet .......................... " << std::endl;
    }//<<>>for ( uInt ijet = 0; ijet < nJets; ijet++ )
    // ** end of jets   ***************************************************************************************************
    if( goodJetEvent ) nGoodJetEvents++;

}//<<>>processJets( fjets, fjetsID, fgenparts )

//-----------------------------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesClusterJets(){

    //    uInt                nCaloJets;
    //    std::vector<float>  cljSeedTOFTimes, cljCMeanTimes;
    //    std::vector<double> cljPt, cljEnergy, cljPhi, cljEta, cljPx, cljPy, cljPz;
    outTree->Branch("cljBcCnt", &cljBcCnt);
    outTree->Branch("nCaloJets", &nCaloJets);
    //outTree->Branch("cljRhIds", &cljRhIds);
    outTree->Branch("cljSeedTOFTime", &cljSeedTOFTime);
    outTree->Branch("cljCMeanTime", &cljCMeanTime);
    outTree->Branch("cljCDrMeanTime", &cljCDrMeanTime);
    outTree->Branch("cljPt", &cljPt);
    outTree->Branch("cljEnergy", &cljEnergy);
    outTree->Branch("cljPhi", &cljPhi);
    outTree->Branch("cljEta", &cljEta);
    outTree->Branch("cljPx", &cljPx);
    outTree->Branch("cljPy", &cljPy);
    outTree->Branch("cljPz", &cljPz);

}//<<>>setBranchesClusterJets()

void KUCMSNtupilizer::processClJet( std::vector<reco::CaloJet> fcalojets );
    
	if( DEBUG ) std::cout << "Processing CaloJets" << std::endl;
    nCaloJets = 0;
    //cljSeedTOFTime.clear();
    cljBcCnt.clear();
    cljCMeanTime.clear();
    cljBc3dEx.clear();
    cljBc3dEy.clear();
    cljBc3dEz.clear();
    cljBc3dEv.clear();
    cljBc3dEslope.clear();
    cljBc3dEchisp.clear();
    cljBc2dEx.clear();
    cljBc2dEy.clear();
    cljBc2dEv.clear();
    cljBc2dEslope.clear();
    cljBc2dEchisp.clear();
    cljCDrMeanTime.clear();
    cljPt.clear();
    cljEnergy.clear();
    cljPhi.clear();
    cljEta.clear();
    cljPx.clear();
    cljPy.clear();
    cljPz.clear();
    //cljRhIds.clear();
    cljEMFrac.clear();

    int iCaloJet(0);
    for (const auto caloJet : fcalojets ){

        cljPt.push_back(caloJet.pt());
        cljEnergy.push_back(caloJet.energy());
        cljPhi.push_back(caloJet.phi());
        cljEta.push_back(caloJet.eta());
        cljPx.push_back(caloJet.px());
        cljPy.push_back(caloJet.py());
        cljPz.push_back(caloJet.pz());
        cljEMFrac.push_back(caloJet.emEnergyFraction());
        iCaloJet++;

        rhGroup cljRhGroup;
        int nCljBcs(0);
        bool first = true;
        auto cjeta = caloJet.eta();
        auto cjphi = caloJet.phi();
        for( const auto cluster : fbclusts ){
            auto cleta = cluster.eta();
            auto clphi = cluster.phi();
            auto cjcldr = std::sqrt(reco::deltaR2( cjeta, cjphi, cleta, clphi ));
            if( cjcldr < 0.4 ){ //std::cout << "CJ - BC Match : " << cjcldr << std::endl;
                nCljBcs++;
                if( first ){ cljRhGroup = getRHGroup( cluster, 0.5 ); first = false; }
                else { auto rhgrp = getRHGroup( cluster, 0.5 ); mrgRhGrp( cljRhGroup, rhgrp ); }
            }//<<>>if( cjcldr < 0.1 )
        }//<<>>for( const auto cluster : fbclusts )
        auto cljRhIdsGroup = getRhGrpIDs( cljRhGroup );
        cljRhIds.push_back(cljRhIdsGroup);
        if( cljRhGroup.size() < 3 ){ 
            cljCMeanTime.push_back( -29.75 ); 
            cljSeedTOFTime.push_back( -29.75 );
        } else { //<<>>if( cljRhGroup.size() < minRHcnt ) ***** IF cljRhGroup.size() > minRHcnt BELOW THIS POINT IN LOOP *****
            auto tofTimes = getLeadTofRhTime( cljRhGroup, vtxX, vtxY, vtxZ );
            auto timeStats = getTimeDistStats( tofTimes, cljRhGroup );
            auto leadJetRh = getLeadRh( cljRhGroup );
            auto seedTOFTime = leadJetRh.time();
            cljSeedTOFTime.push_back(seedTOFTime);
            cljCMeanTime.push_back(timeStats[6]);
            cljBcCnt.push_back(nCljBcs);
        }//<<>>if( cljRhGroup.size() < minRHcnt )

        const auto cljDrGroup = getRHGroup( cjeta, cjphi, 0.4, 2.0 );
        if( cljDrGroup.size() < 13 ){ cljCDrMeanTime.push_back( -29.75 ); }
        else { // <<>> if( cljRhGroup.size() < 3 )
            auto tofTimes = getLeadTofRhTime( cljRhGroup, vtxX, vtxY, vtxZ );
            auto timeStats = getTimeDistStats( tofTimes, cljRhGroup );
            cljCDrMeanTime.push_back(timeStats[6]);
        }//<<>>if( cljRhGroup.size() < 3 )


    }//<<>>for (const auto caloJet : fcalojets )
    nCaloJets = iCaloJet;

}//<<>>processClJet( fcalojets )

//----------------------------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesPhotons(){

    //    uInt                nPhotons;
    //    std::vector<float>  phoSeedTOFTimes, phoCMeanTimes;
    //    std::vector<float>  phoSc3dEx, phoSc3dEy, phoSc3dEz, phoSc3dEv, phoSc3dEslope, phoSc3dEchisp;
    //    std::vector<float>  phoSc2dEx, phoSc2dEy, phoSc2dEv, phoSc2dEslope, phoSc2dEchisp;
    //    std::vector<double> phoPt, phoEnergy, phoPhi, phoEta, phoPx, phoPy, phoPz;
    outTree->Branch("nPhotons", &nPhotons);
    //outTree->Branch("phoIsOotPho", &phoIsOotPho);
    //outTree->Branch("phoExcluded", &phoExcluded);
    outTree->Branch("phoSeedTOFTime", &phoSeedTOFTime);
    outTree->Branch("phoCMeanTime", &phoCMeanTime);
    outTree->Branch("phoPt", &phoPt);
    outTree->Branch("phoEnergy", &phoEnergy);
    outTree->Branch("phoPhi", &phoPhi);
    outTree->Branch("phoEta", &phoEta);
    outTree->Branch("phoPx", &phoPx);
    outTree->Branch("phoPy", &phoPy);
    outTree->Branch("phoPz", &phoPz);
    outTree->Branch("phoRhIds", &phoRhIds);
	//    outTree->Branch("phoIsPFPhoton", &phoIsPFPhoton);
	//    outTree->Branch("phoIsStdPhoton", &phoIsStdPhoton);
    outTree->Branch("phoHasConTracks", &phoHasConTracks);
    outTree->Branch("phoIsPixelSeed", &phoIsPixelSeed);
	//    outTree->Branch("phoIsPhoton", &phoIsPhoton);
    outTree->Branch("phoIsEB", &phoIsEB);
    outTree->Branch("phoIsEE", &phoIsEE);
    //outTree->Branch("phoIsEBGap", &phoIsEBGap);
    //outTree->Branch("phoIsEBEtaGap", &phoIsEBEtaGap);
    //outTree->Branch("phoIsEBPhiGap", &phoIsEBPhiGap);
    //outTree->Branch("phoIsEEGap", &phoIsEEGap);
    //outTree->Branch("phoIsEERingGap", &phoIsEERingGap);
    //outTree->Branch("phoIsEEDeeGap", &phoIsEEDeeGap);
    //outTree->Branch("phoIsEBEEGap", &phoIsEBEEGap);

    outTree->Branch("phoHadOverEM", &phoHadOverEM);
    outTree->Branch("phoHadD1OverEM", &phoHadD1OverEM);
    outTree->Branch("phoHadD2OverEM", &phoHadD2OverEM);
    outTree->Branch("phoHadOverEMVaid", &phoHadOverEMVaid);
    outTree->Branch("phohadTowOverEM", &phohadTowOverEM);
    outTree->Branch("phohadTowD10OverEM", &phohadTowD10OverEM);
    outTree->Branch("phohadTowD20OverEM", &phohadTowD20OverEM);
    outTree->Branch("phohadTowOverEMValid", &phohadTowOverEMValid);
    outTree->Branch("phoE1x5", &phoE1x5);
    outTree->Branch("phoE2x5", &phoE2x5);
    outTree->Branch("phoE3x3", &phoE3x3);
    outTree->Branch("phoE5x5", &phoE5x5);
    outTree->Branch("phoMaxEnergyXtal", &phoMaxEnergyXtal);
    outTree->Branch("phoSigmaEtaEta", &phoSigmaEtaEta);
    outTree->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta);

    outTree->Branch("phoR1x5", &phoR1x5);
    outTree->Branch("phoR2x5", &phoR2x5);
    outTree->Branch("phoR9", &phoR9);
    outTree->Branch("phoFull5x5_e1x5", &phoFull5x5_e1x5);
    outTree->Branch("phoFull5x5_e2x5", &phoFull5x5_e2x5);
    outTree->Branch("phoFull5x5_e3x3", &phoFull5x5_e3x3);
    outTree->Branch("phoFull5x5_e5x5", &phoFull5x5_e5x5);
    outTree->Branch("phoFull5x5_maxEnergyXtal", &phoFull5x5_maxEnergyXtal);
    outTree->Branch("phoFull5x5_sigmaEtaEta", &phoFull5x5_sigmaEtaEta);
    outTree->Branch("phoFull5x5_sigmaIEtaIEta", &phoFull5x5_sigmaIEtaIEta);
    //outTree->Branch("phoFull5x5_r1x5 ", &phoFull5x5_r1x5);
    //outTree->Branch("phoFull5x5_r2x5 ", &phoFull5x5_r2x5);
    outTree->Branch("phoFull5x5_r9", &phoFull5x5_r9);

    //outTree->Branch("phoNSatXtals", &phoNSatXtals);
    //outTree->Branch("phoIsSeedSat", &phoIsSeedSat);

    //outTree->Branch("phoMipChi2", &phoMipChi2);
    //outTree->Branch("phoMipTotEnergy", &phoMipTotEnergy);
    //outTree->Branch("phoMipSlope", &phoMipSlope);
    //outTree->Branch("phoMipInter", &phoMipInter);

    //outTree->Branch("phoMipNHitCone", &phoMipNHitCone);
    //outTree->Branch("phoMipIsHalo", &phoMipIsHalo);

    outTree->Branch("phoEcalRHSumEtConeDR04", &phoEcalRHSumEtConeDR04);
    outTree->Branch("phoHcalTwrSumEtConeDR04", &phoHcalTwrSumEtConeDR04);
    outTree->Branch("phoHcalDepth1TowerSumEtConeDR04", &phoHcalDepth1TowerSumEtConeDR04);
    outTree->Branch("phoCalDepth2TowerSumEtConeDR04", &phoCalDepth2TowerSumEtConeDR04);
    outTree->Branch("phoHcalTowerSumEtBcConeDR04", &phoHcalTowerSumEtBcConeDR04);
    outTree->Branch("phoHcalDepth1TowerSumEtBcConeDR04", &phoHcalDepth1TowerSumEtBcConeDR04);
    outTree->Branch("phoHcalDepth2TowerSumEtBcConeDR04", &phoHcalDepth2TowerSumEtBcConeDR04);
    outTree->Branch("phoTrkSumPtSolidConeDR04", &phoTrkSumPtSolidConeDR04);
    outTree->Branch("phoTrkSumPtHollowConeDR04", &phoTrkSumPtHollowConeDR04);
    outTree->Branch("phoNTrkSolidConeDR04", &phoNTrkSolidConeDR04);
    outTree->Branch("phoNTrkHollowConeDR04", &phoNTrkHollowConeDR04);

    //outTree->Branch("genPhoPt", &genPhoPt);
    //outTree->Branch("genPhoEnergy", &genPhoEnergy);
    //outTree->Branch("genPhoPhi", &genPhoPhi);
    //outTree->Branch("genPhoEta", &genPhoEta);
    //outTree->Branch("genPhoPx", &genPhoPx);
    //outTree->Branch("genPhoPy", &genPhoPy);
    //outTree->Branch("genPhoPz", &genPhoPz);
    //outTree->Branch("genPhoPdgId", &genPhoPdgId);
    //outTree->Branch("genPhoLLP", &genPhoLlp);
    outTree->Branch("phoGenIdx", &phoGenIdx);
    outTree->Branch("phoGenDr", &phoGenDr);

    outTree->Branch("phoSMaj", &phoSMaj);
    outTree->Branch("phoSMin", &phoSMin);
    outTree->Branch("phoSAlp", &phoSAlp);
    outTree->Branch("phoCovEtaEta", &phoCovEtaEta);
    outTree->Branch("phoCovEtaPhi", &phoCovEtaPhi);
    outTree->Branch("phoCovPhiPhi", &phoCovPhiPhi);

}//<<>>setBranchesPhotons()

void KUCMSNtupilizer::processPhotons( std::vector<reco::Photon> fphotons, std::vector<reco::GenParticle> fgenparts ){

    if( DEBUG ) std::cout << "Processing Photons" << std::endl;
    int iGedPhos(0);
    nPhotons = 0;
    phoSeedTOFTime.clear();
    phoCMeanTime.clear();
    phoSc3dEx.clear();
    phoSc3dEy.clear();
    phoSc3dEz.clear();
    phoSc3dEv.clear();
    phoSc3dEslope.clear();
    phoSc3dEchisp.clear();
    phoSc2dEx.clear();
    phoSc2dEy.clear();
    phoSc2dEv.clear();
    phoSc2dEslope.clear();
    phoSc2dEchisp.clear();
    phoPt.clear();
    phoEnergy.clear();
    phoPhi.clear();
    phoEta.clear();
    phoPx.clear();
    phoPy.clear();
    phoPz.clear();
    phoRhIds.clear();

    phoIsPFPhoton.clear();
    phoIsStdPhoton.clear();
    phoHasConTracks.clear();
    phoIsPixelSeed.clear();
    phoIsPhoton.clear();
    phoIsEB.clear();
    phoIsEE.clear();
    phoIsEBGap.clear();
    phoIsEBEtaGap.clear();
    phoIsEBPhiGap.clear();
    phoIsEEGap.clear();
    phoIsEERingGap.clear();
    phoIsEEDeeGap.clear();
    phoIsEBEEGap.clear();

    phoHadOverEM.clear();
    phoHadD1OverEM.clear();
    phoHadD2OverEM.clear();
    phoHadOverEMVaid.clear();
    phohadTowOverEM.clear();
    phohadTowD10OverEM.clear();
    phohadTowD20OverEM.clear();
    phohadTowOverEMValid.clear();
    phoE1x5.clear();
    phoE2x5.clear();
    phoE3x3.clear();
    phoE5x5.clear();
    phoMaxEnergyXtal.clear();
    phoSigmaEtaEta.clear();
    phoSigmaIEtaIEta.clear();

    phoR1x5.clear();
    phoR2x5.clear();
    phoR9.clear();
    phoFull5x5_e1x5.clear();
    phoFull5x5_e2x5.clear();
    phoFull5x5_e3x3.clear();
    phoFull5x5_e5x5.clear();
    phoFull5x5_maxEnergyXtal.clear();
    phoFull5x5_sigmaEtaEta.clear();
    phoFull5x5_sigmaIEtaIEta.clear();
    phoFull5x5_r1x5.clear();
    phoFull5x5_r2x5.clear();
    phoFull5x5_r9.clear();
    phoNSatXtals.clear();
    phoIsSeedSat.clear();

    phoMipChi2.clear();
    phoMipTotEnergy.clear();
    phoMipSlope.clear();
    phoMipInter.clear();
    phoMipNHitCone.clear();
    phoMipIsHalo.clear();

    phoEcalRHSumEtConeDR04.clear();
    phoHcalTwrSumEtConeDR04.clear();
    phoHcalDepth1TowerSumEtConeDR04.clear();
    phoCalDepth2TowerSumEtConeDR04.clear();
    phoHcalTowerSumEtBcConeDR04.clear();
    phoHcalDepth1TowerSumEtBcConeDR04.clear();
    phoHcalDepth2TowerSumEtBcConeDR04.clear();
    phoTrkSumPtSolidConeDR04.clear();
    phoTrkSumPtHollowConeDR04.clear();
    phoNTrkSolidConeDR04.clear();
    phoNTrkHollowConeDR04.clear();

    phoGenIdx.clear();
    phoGenDr.clear();

    phoSMaj.clear();
    phoSMin.clear();
    phoSAlp.clear();

    phoCovEtaEta.clear();
    phoCovEtaPhi.clear();
    phoCovPhiPhi.clear();

    for( const auto photon : fphotons ){

        if( onlyEB && photon.isEE() ) continue;

        phoPt.push_back(photon.pt());
        phoEnergy.push_back(photon.energy());
        phoPhi.push_back(photon.phi());
        phoEta.push_back(photon.eta());
        phoPx.push_back(photon.px());
        phoPy.push_back(photon.py());
        phoPz.push_back(photon.pz());

        phoIsPFPhoton.push_back(photon.isPFlowPhoton()) ;
        phoIsStdPhoton.push_back(photon.isStandardPhoton());
        phoHasConTracks.push_back(photon.hasConversionTracks());
        phoIsPixelSeed.push_back(photon.hasPixelSeed());
        phoIsPhoton.push_back(photon.isPhoton());
        phoIsEB.push_back(photon.isEB());
        phoIsEE.push_back(photon.isEE());
        phoIsEBGap.push_back(photon.isEBGap());
        phoIsEBEtaGap.push_back(photon.isEBEtaGap());
        phoIsEBPhiGap.push_back(photon.isEBPhiGap());
        phoIsEEGap.push_back(photon.isEEGap());
        phoIsEERingGap.push_back(photon.isEERingGap());
        phoIsEEDeeGap.push_back(photon.isEEDeeGap());
        phoIsEBEEGap.push_back(photon.isEBEEGap());

        phoHadOverEM.push_back(photon.hadronicOverEm());
        phoHadD1OverEM.push_back(photon.hadronicDepth1OverEm() );
        phoHadD2OverEM.push_back(photon.hadronicDepth2OverEm() );
        phoHadOverEMVaid.push_back(photon.hadronicOverEmValid());
        phohadTowOverEM.push_back(photon.hadTowOverEm());
        phohadTowD10OverEM.push_back(photon.hadTowDepth1OverEm());
        phohadTowD20OverEM.push_back(photon.hadTowDepth2OverEm() );
        phohadTowOverEMValid.push_back(photon.hadTowOverEmValid());
        phoE1x5.push_back(photon.e1x5());
        phoE2x5.push_back(photon.e2x5());
        phoE3x3.push_back(photon.e3x3());
        phoE5x5.push_back(photon.e5x5());
        phoMaxEnergyXtal.push_back(photon.maxEnergyXtal());
        phoSigmaEtaEta.push_back(photon.sigmaEtaEta() );
        phoSigmaIEtaIEta.push_back(photon.sigmaIetaIeta() );

        phoR1x5.push_back(photon.r1x5());
        phoR2x5.push_back(photon.r2x5());
        phoR9.push_back(photon.r9());
        phoFull5x5_e1x5.push_back(photon.full5x5_e1x5());
        phoFull5x5_e2x5.push_back(photon.full5x5_e2x5());
        phoFull5x5_e3x3.push_back(photon.full5x5_e3x3());
        phoFull5x5_e5x5.push_back(photon.full5x5_e5x5());
        phoFull5x5_maxEnergyXtal.push_back(photon.full5x5_maxEnergyXtal());
        phoFull5x5_sigmaEtaEta.push_back(photon.full5x5_sigmaEtaEta());
        phoFull5x5_sigmaIEtaIEta.push_back(photon.full5x5_sigmaIetaIeta());
        //phoFull5x5_r1x5.push_back(photon.full5x5_r1x5());
        //phoFull5x5_r2x5.push_back(photon.full5x5_r2x5());
        phoFull5x5_r9.push_back(photon.full5x5_r9());

        phoNSatXtals.push_back(photon.nSaturatedXtals());
        phoIsSeedSat.push_back(photon.isSeedSaturated());
        phoMipChi2.push_back(photon.mipChi2());
        phoMipTotEnergy.push_back(photon.mipTotEnergy());
        phoMipSlope.push_back(photon.mipSlope());
        phoMipInter.push_back(photon.mipIntercept());
        phoMipNHitCone.push_back(photon.mipNhitCone());
        phoMipIsHalo.push_back(photon.mipIsHalo());

        phoEcalRHSumEtConeDR04.push_back(photon.ecalRecHitSumEtConeDR04());
        phoHcalTwrSumEtConeDR04.push_back(photon.hcalTowerSumEtConeDR04());
        phoHcalDepth1TowerSumEtConeDR04.push_back(photon.hcalDepth1TowerSumEtConeDR04());
        phoCalDepth2TowerSumEtConeDR04.push_back(photon.hcalDepth2TowerSumEtConeDR04());
        phoHcalTowerSumEtBcConeDR04.push_back(photon.hcalTowerSumEtBcConeDR04() );
        phoHcalDepth1TowerSumEtBcConeDR04.push_back(photon.hcalDepth1TowerSumEtBcConeDR04());
        phoHcalDepth2TowerSumEtBcConeDR04.push_back(photon.hcalDepth2TowerSumEtBcConeDR04());
        phoTrkSumPtSolidConeDR04.push_back(photon.trkSumPtSolidConeDR04());
        phoTrkSumPtHollowConeDR04.push_back(photon.trkSumPtHollowConeDR04());
        phoNTrkSolidConeDR04.push_back(photon.nTrkSolidConeDR04());
        phoNTrkHollowConeDR04.push_back(photon.nTrkHollowConeDR04());

        if( DEBUG ) std::cout << " --- Proccesssing : " << photon << std::endl;
        const auto &phosc = photon.superCluster().isNonnull() ? photon.superCluster() : photon.parentSuperCluster();
        const auto scptr = phosc.get();
        const auto &seedDetId = scptr->seed()->seed(); // seed detid
        const auto isEB = (seedDetId.subdetId() == EcalBarrel); // which subdet
        const auto recHits = ( isEB ? recHitsEB_ : recHitsEE_ );
        const auto & ph2ndMoments = noZS::EcalClusterTools::cluster2ndMoments( *scptr, *recHits);
        const auto  lCov = EcalClusterTools::covariances( *(phosc->seed()), &(*recHits), ecalTopology, ecalGeometry);
        phoSMaj.push_back(ph2ndMoments.sMaj);
        phoSMin.push_back(ph2ndMoments.sMin);
        phoSAlp.push_back(ph2ndMoments.alpha);
        phoCovEtaEta.push_back(lCov[0]);
        phoCovEtaPhi.push_back(lCov[1]);
        phoCovPhiPhi.push_back(lCov[2]);

        scGroup phoSCGroup{*scptr};
        auto phoRhGroup = getRHGroup( phoSCGroup, 0.5 );
        auto phoRhIdsGroup = getRhGrpIDs( phoRhGroup );
        phoRhIds.push_back(phoRhIdsGroup);
        if( DEBUG ) std::cout << " -- gedPhotons : " << scptr << " #: " << phoRhGroup.size() << std::endl;
        auto tofTimes = getLeadTofRhTime( phoRhGroup, vtxX, vtxY, vtxZ );
        auto timeStats = getTimeDistStats( tofTimes, phoRhGroup );
        auto seedTOFTime = getSeedTofTime( *scptr, vtxX, vtxY, vtxZ );
        //auto phoLeadTOFTime =  getLeadTofTime( phoRhGroup, vtxX, vtxY, vtxZ );
        phoSeedTOFTime.push_back(seedTOFTime);
        phoCMeanTime.push_back(timeStats[6]);

        // GenParticle Info for photon  -------------------------------------------------------------------
        if( hasGenInfo ){

            auto genInfo = getGenPartMatch( *scptr, fgenparts );
            if( DEBUG) std::cout << " Photon Match ------------------------- " << std::endl;
            if( DEBUG) std::cout << " Matched Eta: " <<  electron.eta() << " gen: " << genInfo[4] << " match: " << genInfo[2] << std::endl;
            if( DEBUG) std::cout << " Matched Phi: " <<  electron.phi() << " gen: " << genInfo[5] << " match: " << genInfo[3] << std::endl;
            if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
            phoGenIdx.push_back(genInfo[1]);
            phoGenDr.push_back(genInfo[0]);

        }//<<>>if( hasGenInfo )


        }//<<>>if( hasGenInfo )
        iGedPhos++;
    }//<<>>for( const auto photon : *gedPhotons_ )
    nPhotons = iGedPhos;

}//<<>>processPhotons( fphotons, fgenparts )

//-------------------------------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesElectrons(){

    //    uInt                nElectrons;
    //    std::vector<float>  eleSeedTOFTimes, eleCMeanTimes;
    //    std::vector<float>  eleSc3dEx, eleSc3dEy, eleSc3dEz, eleSc3dEv, eleSc3dEslope, eleSc3dEchisp;
    //    std::vector<float>  eleSc2dEx, eleSc2dEy, eleSc2dEv, eleSc2dEslope, eleSc2dEchisp;
    //    std::vector<double> elePt, eleEnergy, elePhi, eleEta, elePx, elePy, elePz;
    outTree->Branch("nElectrons", &nElectrons);
    //outTree->Branch("eleRhIds", &eleRhIds); 
    outTree->Branch("eleSeedTOFTime", &eleSeedTOFTime);
    outTree->Branch("eleCMeanTime", &eleCMeanTime);
    outTree->Branch("elePt", &elePt);
    outTree->Branch("eleEnergy", &eleEnergy);
    outTree->Branch("elePhi", &elePhi);
    outTree->Branch("eleEta", &eleEta);
    outTree->Branch("elePx", &elePx);
    outTree->Branch("elePy", &elePy);
    outTree->Branch("elePz", &elePz);
    outTree->Branch("eleGenIdx", &phoGenIdx);
    outTree->Branch("eleGenDr", &phoGenDr);

}//<<>>setBranchesElectrons()

void KUCMSNtupilizer::processElectrons( std::vector<reco::GsfElectron> felectrons, std::vector<reco::GenParticle> fgenparts ){

    nElectrons = 0;
    eleSeedTOFTime.clear();
    eleCMeanTime.clear();
    elePt.clear();
    eleEnergy.clear();
    elePhi.clear();
    eleEta.clear();
    elePx.clear();
    elePy.clear();
    elePz.clear();
    //eleRhIdsGroup.clear();
    eleGenIdx.clear();
    eleGenDr.clear();

    int iElectros(0);
    scGroup jetEleSCGroup;
    if( DEBUG ) std::cout << "Processing Electrons" << std::endl;
    for( const auto electron : felectrons ){

        elePt.push_back(electron.pt());
        eleEnergy.push_back(electron.energy());
        elePhi.push_back(electron.phi());
        eleEta.push_back(electron.eta());
        elePx.push_back(electron.px());
        elePy.push_back(electron.py());
        elePz.push_back(electron.pz());

        iElectros++;

        if( DEBUG ) std::cout << " --- Proccesssing : " << electron << std::endl;
        const auto &elesc = electron.superCluster().isNonnull() ? electron.superCluster() : electron.parentSuperCluster();
        const auto scptr = elesc.get();
        scGroup eleSCGroup{*scptr};
        auto eleRhGroup = getRHGroup( eleSCGroup, 0.5 );
        //auto eleRhIdsGroup = getRhGrpIDs( eleRhGroup );
        //eleRhIds.push_back(eleRhIdsGroup);
        if( DEBUG ) std::cout << " -- Electrons : " << scptr << " #: " << eleRhGroup.size() << std::endl;
        auto tofTimes = getLeadTofRhTime( eleRhGroup, vtxX, vtxY, vtxZ );
        auto timeStats = getTimeDistStats( tofTimes, eleRhGroup );
        auto seedTOFTime = getSeedTofTime( *scptr, vtxX, vtxY, vtxZ );

        eleSeedTOFTime.push_back(seedTOFTime);
        eleCMeanTime.push_back(timeStats[6]);

        // GenParticle Info for electron  -------------------------------------------------------------------
        if( hasGenInfo ){

			auto genInfo = getGenPartMatch( *scptr, fgenparts );
            if( DEBUG) std::cout << " Electron Match ------------------------- " << std::endl;
            if( DEBUG) std::cout << " Matched Eta: " <<  electron.eta() << " gen: " << genInfo[4] << " match: " << genInfo[2] << std::endl;
            if( DEBUG) std::cout << " Matched Phi: " <<  electron.phi() << " gen: " << genInfo[5] << " match: " << genInfo[3] << std::endl;
            if( DEBUG ) std::cout << " ---------------------------------------------------- " << std::endl;
            eleGenIdx.push_back(genInfo[1]);
            eleGenDr.push_back(genInfo[0]);

        }//<<>>if( hasGenInfo )

    }//<<>>for( const auto electron : *electrons_ )
    nElectrons = iElectros;

}//<<>>processElectrons( felectrons, fgenparts )

//----------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesGenParts(){

    //    uInt nGenParts;
    //    std::vector<double> genPt, genEnergy, genPhi, genEta, genPx, genPy, genPz;
    //    std::vector<int>    genPdgId, genLLP;
    outTree->Branch("nGenParts", &nGenParts);
    outTree->Branch("genPt", &genPt);
    outTree->Branch("genEnergy", &genEnergy);
    outTree->Branch("genPhi", &genPhi);
    outTree->Branch("genEta", &genEta);
    outTree->Branch("genPx", &genPx);
    outTree->Branch("genPy", &genPy);
    outTree->Branch("genPz", &genPz);
    outTree->Branch("genPdgId", &genPdgId);
    outTree->Branch("genLLP", &genLLP);

}//<<>>setBranchesGenParts()

void KUCMSNtupilizer::processGenPart( std::vector<reco::GenParticle> fgenparts ){

    if( DEBUG ) std::cout << "Processing GenParticles" << std::endl;

    nGenParts = 0;
    genPt.clear();
    genEnergy.clear();
    genPhi.clear();
    genEta.clear();
    genPx.clear();
    genPy.clear();
    genPz.clear();
    genPdgId.clear();
    genLLP.clear();
    if( DEBUG ) std::cout << " - enetering Genpart loop" << std::endl;
    auto ngens(0);
    for (const auto genpart : fgenparts ){

        ngens++;
        genPt.push_back(genpart.pt());
        genEnergy.push_back(genpart.energy());
        genPhi.push_back(genpart.phi());
        genEta.push_back(genpart.eta());
        genPdgId.push_back(genpart.pdgId());
        genLLP.push_back(fgenpartllp[ngens-1]);
        genPx.push_back(genpart.px());
        genPy.push_back(genpart.py());
        genPz.push_back(genpart.pz());

    }//<<>> for (const auto genpart : fgenparts )
    nGenParts = ngens;

}//<<>>processGenPart( std::vector<reco::GenParticle>  fgenparts )

//-----------------------------------------------------------------------
void KUCMSNtupilizer::setBranchesRecHits(){

    //    int nRecHits;
    //    std::vector<float>  rhPosX, rhPosY, rhPosZ, rhPosEta, rhPosPhi;
    //    std::vector<uInt>   rhID, rhXtalI1, rhXtalI2, rhSubdet;
    //    std::vector<float>  rhEnergy, rhTime, rhTimeErr, rhTOF;
    //    std::vector<bool>   rhisOOT, rhisGS6, rhisGS1;
    //    std::vector<float>  rhadcToGeV;
    //    std::vector<float>  rhped12, rhped6, rhped1;
    //    std::vector<float>  rhpedrms12, rhpedrms6, rhpedrms1;
    outTree->Branch("nRecHits", &nRecHits);
    outTree->Branch("rhPosX", &rhPosX);
    outTree->Branch("rhPosY", &rhPosY);
    outTree->Branch("rhPosZ", &rhPosZ);
    outTree->Branch("rhPosEta", &rhPosEta);
    outTree->Branch("rhPosPhi", &rhPosPhi);
    outTree->Branch("rhEnergy", &rhEnergy);
    outTree->Branch("rhTime", &rhTime);
    //outTree->Branch("rhTimeErr", &rhTimeErr);
    outTree->Branch("rhTOF", &rhTOF);
    outTree->Branch("rhID", &rhID);
    //outTree->Branch("rhXtalI1", &rhXtalI1);
    //outTree->Branch("rhXtalI2", &rhXtalI2);
    //outTree->Branch("rhSubdet", &rhSubdet);
    outTree->Branch("rhisOOT", &rhisOOT);
    //outTree->Branch("rhAmplitude", &rhAmplitude); // <<<< ---------   calc and save this todo !!!!1
    //outTree->Branch("rhisGS6", &rhisGS6);
    //outTree->Branch("rhisGS1", &rhisGS1);
    //outTree->Branch("rhadcToGeV", &rhadcToGeV);
    //outTree->Branch("rhped12", &rhped12);
    //outTree->Branch("rhped6", &rhped6);
    //outTree->Branch("rhped1", &rhped1);
    //outTree->Branch("rhpedrms12", &rhpedrms12);
    //outTree->Branch("rhpedrms6", &rhpedrms6);
    //outTree->Branch("rhpedrms1", &rhpedrms1);

}//<<>>setBranchesRecHits()

void KUCMSNtupilizer::processRecHits( std::vector<EcalRecHit> frechits ){

    if( DEBUG ) std::cout << "Processing RecHits" << std::endl;

    auto nRecHitCnt(0);
    nRecHits = 0;
    rhPosX.clear();
    rhPosY.clear();
    rhPosZ.clear();
    rhPosEta.clear();
    rhPosPhi.clear();
    rhID.clear();
    rhXtalI1.clear();
    rhXtalI2.clear();
    rhSubdet.clear();
    rhEnergy.clear();
    rhTime.clear();
    rhTimeErr.clear();
    rhTOF.clear();
    rhisOOT.clear();
    rhisGS6.clear();
    rhisGS1.clear();
    rhisWeird.clear();
    rhisDiWeird.clear();
    rhadcToGeV.clear();
    rhSwCross.clear();
    rhped12.clear();
    rhped6.clear();
    rhped1.clear();
    rhpedrms12.clear();
    rhpedrms6.clear();
    rhpedrms1.clear();

    if( DEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    for (const auto recHit : frechits ){

        const auto recHitID = getRawID(recHit);
        auto isEB = getIsEB(recHit); // which subdet
        if( onlyEB && not isEB ) continue;
        const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        auto recHitPos = geometry->getGeometry(recHit.detid())->getPosition();
        const auto rhX = recHitPos.x();
        const auto rhY = recHitPos.y();
        const auto rhZ = recHitPos.z();
        const auto d_rh = hypo(rhX,rhY,rhZ)/SOL;
        const auto d_pv = hypo(rhX-vtxX,rhY-vtxY,rhZ-vtxZ)/SOL;
        float swisscross(0.0);

        if( isEB ) swisscross = EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true);
        else swisscross = EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);

        const auto laser = laserH->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        const float adcToGeV( isEB ? adcToGeVEB : adcToGeVEE );
        const auto & pediter = pedestalsH->find(recHitID);

        rhID.push_back(recHitID);
        rhPosX.push_back(rhX);
        rhPosY.push_back(rhY);
        rhPosZ.push_back(rhZ);
        rhTOF.push_back(d_pv);
        rhPosEta.push_back(recHitPos.eta());
        rhPosPhi.push_back(recHitPos.phi());
        rhTime.push_back(recHit.time()-d_rh); // Note : Margret adds d_rh to the time in her code & subtracts d_pv ( or TOF )
        rhTimeErr.push_back(recHit.timeError());
        rhisOOT.push_back(recHit.checkFlag(EcalRecHit::kOutOfTime));
        rhEnergy.push_back(recHit.energy());
        rhSwCross.push_back(swisscross);
        rhisWeird.push_back(recHit.checkFlag(EcalRecHit::kWeird));
        rhisDiWeird.push_back(recHit.checkFlag(EcalRecHit::kDiWeird));
        rhisGS6.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain6));
        rhisGS1.push_back(recHit.checkFlag(EcalRecHit::kHasSwitchToGain1));
        rhadcToGeV.push_back(laser*interCalib*adcToGeV);

        if (pediter != pedestalsH->end()){
            const auto & ped = (*pediter);
            rhped12.push_back(ped.mean(1));
            rhped6.push_back(ped.mean(2));
            rhped1.push_back(ped.mean(3));
            rhpedrms12.push_back(ped.rms(1));
            rhpedrms6.push_back(ped.rms(2));
            rhpedrms1.push_back(ped.rms(3));
        } else {
            rhped12.push_back(0.0);
            rhped6.push_back(0.0);
            rhped1.push_back(0.0);
            rhpedrms12.push_back(0.0);
            rhpedrms6.push_back(0.0);
            rhpedrms1.push_back(0.0);
        }//<<>>if (pediter != pedestalsH->end())
        nRecHitCnt++;

    }//<<>>for (const auto recHit : *recHitsEB_ )   
    nRecHits = nRecHitCnt;

}//<<>>processRecHits( std::vector<EcalRecHit> frechits )


////////////////////////////////////////////////////////////////////////////////////////////////////
// ---------------- CMSSW Ana Helper Functions ---------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////

int KUCMSNtupilizer::getPFJetID(const reco::PFJet & jet){

    const auto eta  = std::abs(jet.eta());
    const auto NHF  = jet.neutralHadronEnergyFraction();
    const auto NEMF = jet.neutralEmEnergyFraction();
    const auto CHF  = jet.chargedHadronEnergyFraction();
    const auto CEMF = jet.chargedEmEnergyFraction();
    const auto NHM  = jet.neutralMultiplicity();
    const auto CHM  = jet.chargedMultiplicity();
    const auto SHM  = jet.chargedMultiplicity()+jet.neutralMultiplicity();
    const auto MUF  = jet.muonEnergyFraction();

    int tighter = 3;
    int tightLepVeto = 0;
    int tight = 2;
    int loose = 1;

    bool nhfup  = NHF  < 0.90;
    bool nhflw  = NHF  > 0.2;

    bool nemfup1 = NEMF < 0.90;
    bool nemfup2 = NEMF < 0.99;
    bool nemf80 = NEMF < 0.80;
    bool nemflw = NEMF > 0.01;
    bool nemf10 = NEMF > 0.1;

    bool shm1  = SHM  > 1;
    bool muf8  = MUF  < 0.80;
    bool chf0  = CHF  > 0;
    bool chf10  = CHF  > 0.10;
    bool chm0  = CHM  > 0;
    bool cemf8 = CEMF > 0.80;
    bool nhm2  = NHM  > 1;
    bool nhm10 = NHM  > 10;

    bool eta1 = eta <= 2.6;
    bool eta2 = eta <= 2.7;
    bool eta3 = eta <= 3.0;

    if (eta1){
        if      (nhfup && nemfup1 && shm1 && muf8 && chf0 && chm0 && cemf8) return tightLepVeto;
        else if (nhfup && nemf80 && shm1 && chf10 && chm0) return tighter;
        else if (nhfup && nemfup1 && shm1 && chf0 && chm0) return tight;
        else    return loose;
    } else if (!eta1 && eta2 ){ //<<>>if (eta1)
        if      (nhfup && nemfup2 && chm0 && muf8 && cemf8) return tightLepVeto;
        else if (nhfup && nemf80 && chm0) return tighter;
        else if (nhfup && nemfup2 && chm0) return tight;
        else    return loose;
    } else if (!eta2 && eta3){ //<<>>if (eta1) : else if
        if      (nemf10 && nemf80 && nhm2) return tighter;
        else if (nemflw && nemfup2 && nhm2) return tight;
        else    return loose;
    } else { //<<>>if (eta1) : else if : else if
        if      (nhflw && nemfup1 && nhm10) return tight;
        else    return loose;
    }//<<>>if (eta1) : else if : else if : else

    return -1; // should not happen

}//>>>>int KUCMSNtupilizer::getPFJetID(const reco::PFJet & jet)

rhGroup KUCMSNtupilizer::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 ){

    rhGroup result;
    for (const auto recHit : *recHitsEB_ ){
        if( recHit.checkFlag(EcalRecHit::kHasSwitchToGain6) ) continue;
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId(recHit.detid());
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        const auto dr = std::sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
        if( dr > drmin ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        if( recHit.checkFlag(EcalRecHit::kHasSwitchToGain6) ) continue;
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId(recHit.detid());
        const auto recHitPos = endcapGeometry->getGeometry(recHitId)->getPosition();
        const auto dr = std::sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
        if( dr > drmin ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup KUCMSNtupilizer::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 )

rhGroup KUCMSNtupilizer::getRHGroup( uInt detid ){

    rhGroup result;
    for (const auto recHit : *recHitsEB_ ){
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( detid != rawId ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( detid != rawId ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;
}//>>>>rhGroup KUCMSNtupilizer::getRHGroup( uInt detid )

rhGroup KUCMSNtupilizer::getRHGroup(){

    rhGroup result;
    for (const auto recHit : *recHitsEB_ ) result.push_back(recHit);
    for (const auto recHit : *recHitsEE_ ) result.push_back(recHit);

    return result;

}//>>>>rhGroup KUCMSNtupilizer::getRHGroup()


rhGroup KUCMSNtupilizer::getRHGroup( const scGroup superClusterGroup, float minenr ){

    rhGroup result;
    vector<uInt> rawIds;
    for ( const auto superCluster : superClusterGroup ){
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )  
    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup KUCMSNtupilizer::getRHGroup( const scGroup superClusterGroup, float minenr = 0.0 )

rhGroup KUCMSNtupilizer::getRHGroup( const scGroup superClstrGrp, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax ){

    rhGroup result;
    vector<uInt> rawIds;

    float dRcut = 0.05;

    int iter = -1;
    for ( const auto superCluster : superClstrGrp ){
        iter++;
        if( phDr[iter] > dRcut ) continue;
        if( phEnergy[iter] != phEnMax ) continue;
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClstrGrp )   
    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup KUCMSNtupilizer::getRHGroup( const scGroup superClstrGrp, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax )

rhGroup KUCMSNtupilizer::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 ){

    rhGroup result;
    vector<uInt> rawIds;
    auto & hitsAndFractions = basicCluster.hitsAndFractions();
    const auto nHAF = hitsAndFractions.size();
    //std::cout << " --- getRHGroup for basic cluster with " << nHAF << " rechits " << std::endl;
    for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
        const auto detId = hitsAndFractions[iHAF].first;
        const auto rawId = detId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
    }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    for (const auto recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )       

    return result;

}////rhGroup KUCMSNtupilizer::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 )

vector<float> KUCMSNtupilizer::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    vector<float> result;
	if( rechits.size() == 0 ){ result.push_back(-99.0); return result;
    for (const auto recHit : recHits ){
        const auto rht = recHit.time();
        //std::cout << " ----- Get TOF Time rh time: " << rht << std::endl;
        const auto recHitId(recHit.detid());
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        const auto rhPosX = recHitPos.x();
        const auto rhPosY = recHitPos.y();
        const auto rhPosZ = recHitPos.z();
        //std::cout << " ----- Get TOF Time rh POS: " << rhPosX << " " <<   rhPosY << " " << rhPosZ << std::endl;
        const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
        const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
        const auto tof = (d_rh-d_pv)/SOL;
        //std::cout << " ----- Get TOF Time rh tof: " << tof << std::endl;
        result.push_back(rht-tof);
    }//<<>>for (const auto recHit : recHits )   

    return result;

}//>>>>vector<float>  KUCMSNtupilizer::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float KUCMSNtupilizer::getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ ){

    const auto recHitId(rechit.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    return tof;

}//>>>>>float getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ )

EcalRecHit KUCMSNtupilizer::getLeadRh( rhGroup recHits ){

    EcalRecHit result;
    float enr(0.0);
    for (const auto recHit : recHits ){
        auto rhenr = recHit.energy();
        if( rhenr < enr ) continue;
        enr = rhenr;
        result = recHit;
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit KUCMSNtupilizer::getLeadRh( rhGroup recHits )

vector<float> KUCMSNtupilizer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    vector<float> result;
	if( rechits.size() == 0 ){ result.push_back(-99.0); return result;
    auto lrh = getLeadRh(recHits);
    const auto recHitId(lrh.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for (const auto recHit : recHits ){result.push_back(recHit.time()-tof);}
    return result;

}//>>>>vector<float>  KUCMSNtupilizer::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float KUCMSNtupilizer::getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ ){

    float rhTime(-100.0);
    const auto seedDetId = sprclstr.seed().get()->seed(); // seed detid
    const auto seedRawId = seedDetId.rawId(); // crystal number
    const auto isEB = (seedDetId.subdetId() == EcalSubdetector::EcalBarrel); // which subdet
    const auto recHits = (isEB ? *recHitsEB_ : *recHitsEE_ );
    for( const auto rechit : recHits ){
        const auto recHitId = rechit.detid();
        const auto rawId = recHitId.rawId();
        if( rawId == seedRawId ){ rhTime = rechit.time(); continue; }
    }//<<>>for (const auto recHit : *recHitsEE_ ) 

    const auto recHitPos = barrelGeometry->getGeometry(seedDetId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    const auto seedTofTime = rhTime - tof;
    return seedTofTime;

}//>>>>float  KUCMSNtupilizer::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float KUCMSNtupilizer::getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

	if( rechits.size() == 0 ) return -99.0;
    auto lrh = getLeadRh(recHits);
    const auto recHitId(lrh.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    return lrh.time()-tof;

}//>>>>float  KUCMSNtupilizer::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

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

void KUCMSNtupilizer::mrgRhGrp( rhGroup & x, rhGroup & y ){

    //std::cout << " --- In mrgRhGrp " << std::endl;
    bool matched(false);
    for(const auto rhy : y ){
        for(const auto rhx : x ){ if( rhMatch(rhx,rhy) ) matched = true; }
        if( not matched ) x.push_back(rhy);
    }//<<>>for(const auto rhx : y )
   //std::cout << " ---- Finished mrgRhGrp " << std::endl;

    return;

}//>>>>rhGroup KUCMSNtupilizer::mrgRhGrp(const rhGroup x, const rhGroup y)

bool KUCMSNtupilizer::reduceRhGrps( vector<rhGroup> & x ){

    //std::cout << " --- In reduceRhGrps " << std::endl;
    uInt s = x.size();
    uInt a(0), b(0);
    bool match(false);
    for( uInt ita(0); ita < s; ita++ ){
        for( auto itb = ita+1; itb < s; itb++ ){
            if( dupRhFnd(x[ita], x[itb]) ){ a = ita; b = itb; match = true; }
            if(match) break;
    }//<<>>for( int itb = ita+1; itb < s; itb++ ) 
        if(match) break;
    }//<<>>for( int ita = 0; ita < nBcRhGroups; ita++ )
    if(match){
        mrgRhGrp( x[a], x[b] );
        x.erase(x.begin()+b);
    } //<<>>if(match) : else
    //std::cout << " ---- Finished reduceRhGrps " << std::endl;

    return match;

}//>>>>vector<rhGroup> KUCMSNtupilizer::reduceRhGrps(const vector<rhGroup> x)

string KUCMSNtupilizer::bigKidChase( std::vector<reco::CandidatePtr> kids, float vx ){

    std::string result("");
    float maxpt(0.0);
    const reco::Candidate* bigKid(NULL);
    for( auto kid : kids ){
        if( kid->pt() > maxpt ){
            maxpt = kid->pt();
            bigKid = kid.get();
        }//<<>>if( kid->pt() > maxpt )
    }//<<>>for( auto kid : kids )
    result += std::string(" -- BigKid : pdgID : ") + bigKid->pdgId() + std::string(" pt : ") + bigKid->pt();
    result += std::string(" vtx : ") + bigKid->mother(0)->vx() + std::string(" / ") + vx ;
    return result;

}//<<>>string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx )

void KUCMSNtupilizer::motherChase( const reco::Candidate* kid, string & depth ){

    depth += "-";
    for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        std::cout <<  depth  << " gMomID : " << kid->mother(gmit)->pdgId() << " pt : " << kid->mother(gmit)->pt();
        std::cout << " Vertix (" << kid->mother(gmit)->vx() << "," << kid->mother(gmit)->vy() << "," << kid->mother(gmit)->vz() << ")";
        std::cout << " nGMothers " << kid->mother(gmit)->numberOfMothers() << std::endl;
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    if( not depth.empty() ) depth = depth.substr (0,depth.length()-1);

}//<<>> void MotherChase( Candidate* kid, string depth  )

int KUCMSNtupilizer::llpGenChase( const reco::Candidate* kid ){

    auto match = ( std::abs(kid->pdgId()) > 1000020 ) && ( std::abs(kid->pdgId()) < 1000040 );
    if( match ) return kid->pdgId();
    for( long unsigned int gmit(0); gmit < kid->numberOfMothers(); gmit++ ){
        if( kid->mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid->mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

int KUCMSNtupilizer::llpGenChase( const reco::GenParticle kid ){

    auto match = ( std::abs(kid.pdgId()) > 1000020 ) && ( std::abs(kid.pdgId()) < 1000040 );
    auto hard = kid.isLastCopy();
    if( match && hard ) return kid.pdgId();
    for( long unsigned int gmit(0); gmit < kid.numberOfMothers(); gmit++ ){
        //if( kid->pdgId() == 6000113 ) return true;
        if( kid.mother(gmit)->numberOfMothers() > 0 ){ auto r = llpGenChase(kid.mother(gmit)); if(r != 0) return r; }
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return 0;

}//<<>> void llpChase( Candidate* kid )

int KUCMSNtupilizer::llpGenChaseP( const reco::Candidate* kid, int depth ){

    int nMoms = kid->numberOfMothers();
    if( ! nMoms ) return 999;
    int genPartID = 1000;
    auto kidPdgID = std::abs(kid->pdgId());
    int wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? 1 : 0;
    for( int gmit(0); gmit < nMoms; gmit++ ){
        auto genmom = kid->mother(gmit);
        auto momPdgID = std::abs(genmom->pdgId());
        if( momPdgID == 1000022 && kidPdgID == 22 && depth == 0 ) return 1;
        else if( momPdgID == 1000022 ) { if( wzFlag ) return 13; else return 12; }
        else if( momPdgID == 1000023 ){ if( wzFlag ) return 23; else return 22; }
        else if( momPdgID == 1000025 ){ if( wzFlag ) return 33; else return 32; }
        else if( momPdgID == 1000035 ){ if( wzFlag ) return 43; else return 42; }
        else if( momPdgID == 1000024 && kidPdgID == 22 && depth == 0 ) return 111;
        else if( momPdgID == 1000024 ) { if( wzFlag ) return 113; else return 112; }
        else if( momPdgID == 1000037 ){ if( wzFlag ) return 123; else return 122; }
        else if( (momPdgID >= 1000001 && momPdgID <= 1000006) || momPdgID == 1000021 ) return 500;

        auto mGenPartID = llpGenChaseP(kid->mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return genPartID;

}//<<>> void llpChase( const reco::Candidate* kid, int wzFlag ){)

int KUCMSNtupilizer::llpGenChaseP( const reco::GenParticle kid, int depth ){

    int nMoms = kid.numberOfMothers();
    if( ! nMoms ) return 999;
    int genPartID = 1000;
    auto kidPdgID = std::abs(kid.pdgId());
    int wzFlag = ( kidPdgID == 23 || kidPdgID == 24 ) ? 1 : 0;
    for( int gmit(0); gmit < nMoms; gmit++ ){
        auto genmom = kid.mother(gmit);
        auto momPdgID = std::abs(genmom->pdgId());

        if( momPdgID == 1000022 && kidPdgID == 22 && depth == 0 ) return 1;
        else if( momPdgID == 1000022 ) { if( wzFlag ) return 13; else return 12; }
        else if( momPdgID == 1000023 ){ if( wzFlag ) return 23; else return 22; }
        else if( momPdgID == 1000025 ){ if( wzFlag ) return 33; else return 32; }
        else if( momPdgID == 1000035 ){ if( wzFlag ) return 43; else return 42; }
        else if( momPdgID == 1000024 && kidPdgID == 22 && depth == 0 ) return 111;
        else if( momPdgID == 1000024 ) { if( wzFlag ) return 113; else return 112; }
        else if( momPdgID == 1000037 ){ if( wzFlag ) return 123; else return 122; }
        else if( (momPdgID >= 1000001 && momPdgID <= 1000006) || momPdgID == 1000021 ) return 500;

        auto mGenPartID = llpGenChaseP(kid.mother(gmit),depth+1);
        if( mGenPartID < genPartID ) genPartID = mGenPartID;
    }//<<>>for( long unsigned int gmit(0); gmit < nKMother; gmit++ )
    return genPartID;

}//<<>> void llpChase( Candidate* kid )

void KUCMSNtupilizer::kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz ){

    for( auto kid : kids ){

        bool done(false);
        auto mom = kid->mother();
        std::cout << " -- Kid pdgID : " << kid->pdgId() << " p : " << kid->p() << " e : " << kid->energy() << std::endl;
        bool top(true);
        while( not done ){
            if( mom->pt() == 0 ) done = true;
            else {
                auto gmom = mom->mother();
                if( top ){
                    std::cout << " ---- KidStep pdgID : " << mom->pdgId() << " p : " << mom->p() << " e : " << mom->energy();
                    std::cout << " vert : (" << mom->vx() << "," << mom->vy() << "," << mom->vz() << ")" << std::endl;
                }//<<>>if( top )
                if( gmom->pt() == 0 ) done = true;
                else {
                    if( mom->vx() == gmom->vx() ) top = false;
                    else top = true;
                    mom = gmom;
                }//<<>>if( gmom->pt() == 0 )
            }//<<>>if( mom->pt() == 0 )                 
        }//<<>>while( not done )
    }//<<>>for( auto kid : kids )

}//<<>>string kidChase( std::vector<reco::CandidatePtr> kids, float vx )

vector<float> KUCMSNtupilizer::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  ){
// redo this function to give tof and impact angle for input gen particle

    vector<float> result;
    vector<float> kidtime;
    vector<float> kide;
    vector<float> impact;
    bool llpcut(false);
    float nextBX(0.0);
    float llpkid(0.0);
    for( auto kid : kids ){

        bool done(false);
        bool first(true);
        bool llp(false);
        bool bquark(false);
        auto mom = kid->mother();
        vector<double> stepp{kid->p()};
        vector<double> stepe{kid->energy()};
        vector<double> orignvx;
        vector<double> orignvy;
        vector<double> orignvz;
        bool top(true);
        int steps(0);
        while( not done ){
            auto pdgid = std::abs(mom->pdgId());
            if( pdgid == 6000113 ) llp = true;
            if( pdgid == 5 ) bquark = true;
            if( std::floor(pdgid/100) == 5 ) bquark = true;
            if( mom->pt() == 0 ) done = true;
            else {
                auto gmom = mom->mother();
                if( top ){
                    if( first ) first = false;
                    else {
                        stepp.push_back(mom->p());
                        stepe.push_back(mom->energy());
                    }//<<>>if( first )
                    steps++;
                    //if( llp || bquark ) llpstep.push_back(1.0); else llpstep.push_back(0.0);
                    orignvx.push_back(mom->vx());
                    orignvy.push_back(mom->vy());
                    orignvz.push_back(mom->vz());
                }//<<>>if( top )
                if( gmom->pt() == 0 ) done = true;
                else {
                    if( mom->vx() == gmom->vx() ) top = false;
                    else top = true;
                    mom = gmom;
                }//<<>>if( gmom->pt() == 0 )
            }//<<>>if( mom->pt() == 0 )                 
        }//<<>>while( not done )
        first = true;

        if( DEBUG && false ) std::cout << " ---- jetGenTime Calc Steps : " << steps << std::endl;
        bool stepcut( steps == 0 );
        if( DEBUG && false ) std::cout << " ---- jetGenTime Calc llpcut : " << llpcut << std::endl;
        if( not stepcut ){
            float maxe(0.0);
            float totaltime(0.0);
            auto destvx(cx);
            auto destvy(cy);
            auto destvz(cz);
            for( int it(0); it < steps; it++ ){
                if( stepe[it] > maxe ) maxe = stepe[it];
                auto beta = stepp[it]/stepe[it];
                auto legtime = hypo( (destvx-orignvx[it]), (destvy-orignvy[it]), (destvz-orignvz[it]) )/(SOL*beta);
                totaltime += legtime;
                //std::cout << "Sum path legs beta: " << beta << " x: " << (destvx-orignvx[it]) << " y: " << (destvy-orignvy[it]);
                //std::cout << " z: " << (destvz-orignvz[it]) << " legtime: " << legtime << " totaltime: " << totaltime << std::endl;
                //std::cout << " -- dest : " << destvx << ", " << destvy << ", " << destvz; 
                //std::cout << " orig : " << orignvx[it] << ", " << orignvy[it] << ", " << orignvz[it] << std::endl;
                destvx = orignvx[it];
                destvy = orignvy[it];
                destvz = orignvz[it];
            }//<<>>for( int it(0); it < steps; it++ )
            auto r = hypo((cx-orignvx[0]), (cy-orignvy[0]));
            auto iangle = getAngle((cz-orignvz[0]),r);
            //std::cout << "getAngle: " << iangle << " r: " << r << " z: " << cz-orignvz[0] << " x: " << cx-orignvx[0];
            //std::cout << " y: " << cy-orignvy[0] << std::endl; 
            if( totaltime > 25.0 ) {
                nextBX = 1.0;
            } else {
                impact.push_back(iangle);
                kidtime.push_back(totaltime);
                kide.push_back(maxe);
                llpcut = llpcut || llp || bquark;
                if( llp || bquark ) llpkid++;
            }//<<>>if( totaltime > 25.0 )
        }//<<>>if( steps == 0 )

    }//<<>>for( auto kid : kids )

    auto mtime(-30.0);
    auto mangle(7.0);
    auto tvar(-2.0);
    auto ktsize = kidtime.size();
    //if( llpcut && ( ktsize > 0 ) ){
    if( ktsize > 0 ){
        if( ktsize > 1 ){
            mtime = mean(kidtime,kide);
            mangle = mean(impact,kide);
            tvar = var(kidtime,mtime,kide);
        } else {
            mtime = kidtime[0];
            mangle = impact[0];
            tvar = -1.0;
        }//<<>>if( kidtime.size() > 1 )
    }//<<>>if( llpcut )
    result.push_back(mtime);//0
    result.push_back(mangle);//1
    result.push_back(tvar);//2
    result.push_back(nextBX);//3
    if( llpcut ) result.push_back(1.0); else result.push_back(0.0);//4
    result.push_back(llpkid/ktsize);//5
    result.push_back(ktsize);//6
    return result;

}//>>>>vector<float> KUCMSNtupilizer::kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  )


vector<float>  KUCMSNtupilizer::getGenPartMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts ){

    // GenParticle Info   -------------------------------------------------------------------
    if( DEBUG ) std::cout << "Getting phoGenParton Information" << std::endl;
    // set defaults for no match here
    vector<float> results;
    auto mgenDr = -1.0;// 0
    bool matchfound(false);
    float goodDr(0.4);
    int matchedIdx(-1);// 1
    int index(0);
    float meta(9);// 2
    float mphi(9);// 3
    float geta(9);// 4
    float gphi(9);// 5
    for(const auto& genPart : fgenparts ){

        const auto rhX = scptr->x();
        const auto rhY = scptr->y();
        const auto rhZ = scptr->z();
        const auto gnX = genPart.vx();
        const auto gnY = genPart.vy();
        const auto gnZ = genPart.vz();
        auto cphoEta = std::asinh((rhZ-gnZ)/hypo(rhX-gnX,rhY-gnY));
        auto cphoPhi = std::atan2(rhY-gnY,rhX-gnX);
        auto dr = std::sqrt(reco::deltaR2(genPart.eta(), genPart.phi(), cphoEta, cphoPhi ));
        if( dr < goodDr ){
            matchfound = true;
            goodDr = dr;
            matchedIdx = index;
            meta = cphoEta;
            mphi = cphoPhi;
            geta = genPart.eta();
            gphi = genPart.phi();
        }//<<>>if( jtgjdr <= goodDr )
        index++;

    }//<<>>for(const auto& genPart : fgenparts  ) 
    if( matchedIdx >= 0 ){ mgenDr = goodDr; }//<<>>if( matchedIdx >= 0 ) 

	results.push_back( mgenDr );
    results.push_back( matchedIdx );
    results.push_back( meta );
    results.push_back( mphi );
    results.push_back( geta );
    results.push_back( gphi );
	return results;

}//<<>>getGenPartMatch( reco::SuperClusterCollection *scptr, std::vector<reco::GenParticle> fgenparts )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ define this as a plug-in --------------------------------------------------------------------------------
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
DEFINE_FWK_MODULE(KUCMSNtupilizer);

