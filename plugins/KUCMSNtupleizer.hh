// -*- C++ -*-
//
// Package:    KUCMSNtupilizer
// Class:      KUCMSNtupilizer
//
/**\class KUCMSNtupilizer LLPgammaAnalyzer.cc LLPgammaAnalyzer/plugins/KUCMSNtupilizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

// basic C++ types
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// HLT + Trigger info
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// Gen Info
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

// DataFormats
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

// DetIds 
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

// Ecal RecHits
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"

// Supercluster info
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

// EGamma Tools
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

// Geometry
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

// Topology 
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

// ECAL Record info (Laser Constants)
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"

// ECAL Record info (Intercalibration Constants)
#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"

// ECAL Record info (ADCToGeV)
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

// ECAL Record info (ADCToGeV)
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"

// ECAL Record info (Pedestals)
#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

// JECS
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

// JERs
#include "CondFormats/JetMETObjects/interface/JetResolutionObject.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"

// TOOLS
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/RefToPtr.h"

// ROOT
#include "TH1.h"
#include "TH2.h"
#include "TFormula.h"
#include "TF1.h"
#include "TTree.h"
#include "Math/PositionVector3D.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TMatrixDSymEigen.h"
#include "TGraph.h"
#include "TMathBase.h"

#ifndef KUCMSNtupilizerHeader
#define KUCMSNtupilizerHeader

using namespace std;
using namespace edm;

//
// In class declaration :
// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.
//

//
// sort functions
//

const auto sortByPt = [](auto & obj1, auto & obj2) {return obj1.pt() > obj2.pt();};

//
//  constants, enums and typedefs
//

using reco::TrackCollection;

typedef edm::View<reco::Candidate> CandidateView;
typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> recHitCol;
typedef vector<EcalRecHit> rhGroup;
typedef vector<reco::SuperCluster> scGroup;
typedef vector<reco::CaloCluster> bcGroup;

typedef unsigned int uInt;
typedef	unsigned long int ulInt;
typedef const float CFlt;
typedef const double CDbl;
typedef const vector<float> CVFlt;
typedef const vector<double> CVDbl;

typedef vector<uInt> rhIdGroup;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;

#define nEBEEMaps 0
#define nHists 4
//const float sol = 29.9792458; // speed of light in cm/ns
#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ...  
#define BUNCHES 3564
#define SAMPLES 10

enum class ECAL {EB, EM, EP, EE, NONE};

//
//  Class Declaration
//

class KUCMSNtupilizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

	public:

    	explicit KUCMSNtupilizer(const edm::ParameterSet&);
      	~KUCMSNtupilizer();

      	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:

    	virtual void beginJob() override;
      	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      	virtual void endJob() override;

        void setBranchesEvent();
		void processEvent( const edm::Event& iEvent );
        void setBranchesVtx();
		void processVtx();

        void setBranchesGenParts();
        void processGenPart();

        void setBranchesMet();
		void processMet();
        void setBranchesPhotons();
		void processPhotons();
        void setBranchesElectrons();
		void processElectrons();
        void setBranchesJets();
		void processJets();
        void setBranchesClusterJets();
        void processClJet();	

        void setBranchesRecHits();
        void processRecHits( float minRHEf );

		////////////////////////////////////////////////////
      	// ----------member data ---------------------------
		////////////////////////////////////////////////////

		// filtered collection vectors
    	std::vector<reco::PFJet>        fjets;
    	std::vector<int>                fjetsID;
    	std::vector<reco::CaloCluster>  fbclusts;
    	std::vector<EcalRecHit>         frechits;
    	std::vector<bool>               frhused;
    	std::vector<reco::Photon>       fphotons;
    	//std::vector<reco::Photon>     footphotons;
        std::vector<int>                fphotonID;
    	std::vector<reco::GsfElectron>  felectrons;
        std::vector<int>                felectronID;
    	std::vector<reco::CaloJet>      fcalojets;
    	std::vector<reco::Muon>         fmuons;
    	std::vector<reco::PFMET>        fpfmet;
    	std::vector<reco::GenParticle>  fgenparts;
    	std::vector<int>                fgenpartllp;
    	std::vector<reco::GenJet>       fgenjets;

		// oputput tree
      	TTree *outTree;

		// histograms
		TH1D *hist1d[nHists];
      	TH2D *hist2d[nHists];

		// Flags
		const bool hasGenInfo;

		////////////////////////////////////////////////////////////
		// ----------- collection and branch declaration ----------
		////////////////////////////////////////////////////////////

      	// Event
      	ulInt 	event; // technically unsigned long long in Event.h...
      	uInt 	run, lumi; 

      	// Tracks
      	const edm::InputTag tracksTag;
      	edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
      	edm::Handle<std::vector<reco::Track>> tracks_;

		// Trigger Info -> From miniAOD code not converted
		// trigger inputs
        const std::string inputPaths;
        std::vector<std::string> pathNames;
        //strBitMap triggerBitMap;
        const std::string inputFilters;
        std::vector<std::string> filterNames;
		// trigger results
        const edm::InputTag triggerResultsTag;
        edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;
        edm::Handle<edm::TriggerResults> triggerResultsH;
        // trigger objects
        const edm::InputTag triggerObjectsTag;
        edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsToken;
        edm::Handle<std::vector<pat::TriggerObjectStandAlone> > triggerObjectsH;
        // output triggers
        std::map<std::string,std::vector<pat::TriggerObjectStandAlone> > triggerObjectsByFilterMap; 
			// first index is filter label, second is trigger objects

      	// PF Candidates View
      	const edm::InputTag pfcandTag;
      	edm::EDGetTokenT<edm::View<reco::Candidate>> pfcand_token_;
      	edm::Handle<edm::View<reco::Candidate>> pfcands_;

		// PF Candidates
        const edm::InputTag pfCanTag;
        edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCan_token_;
		edm::Handle<std::vector<reco::PFCandidate>> pfCans_; 
        const edm::InputTag pfCanPhoMapTag;
		edm::EDGetTokenT<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanPhotonMap_token_;
        edm::Handle<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanPhotonMap_;
        const edm::InputTag pfCanOOTPhoMapTag;
        edm::EDGetTokenT<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanOOTPhotonMap_token_;
        edm::Handle<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanOOTPhotonMap_;
        const edm::InputTag pfCanEleMapTag;
        edm::EDGetTokenT<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanElectronMap_token_;
        edm::Handle<edm::ValueMap<edm::Ptr<reco::PFCandidate>>> pfCanElectronMap_;

      	// vertices
      	const edm::InputTag verticesTag;
      	edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
      	edm::Handle<std::vector<reco::Vertex>> vertices_;

      	int		nVtx;
      	float 	vtxX, vtxY, vtxZ;

      	// mets
      	const edm::InputTag pfmetsTag;
      	edm::EDGetTokenT<std::vector<reco::PFMET>> pfmetsToken_;
      	edm::Handle<std::vector<reco::PFMET>> pfmets_;
		//reco::PFMET t1pfMET;	

        float  metSumEt, metPt, metPx, metPy, metPhi;
        float  metCSumEt, metCPx, metCPy;

		// Met Filters -> From miniAOD code not converted
  		// met filter inputs
  		const std::string inputFlags;
  		std::vector<std::string> flagNames;
  		//strBitMap triggerFlagMap;
  		// met filters
  		const edm::InputTag triggerFlagsTag;
  		edm::EDGetTokenT<edm::TriggerResults> triggerFlagsToken;
  		edm::Handle<edm::TriggerResults> triggerFlagsH;

	   	// supercluster
      	const edm::InputTag superClusterCollectionTag;
      	edm::EDGetTokenT<reco::SuperClusterCollection> scToken_;
      	edm::Handle<reco::SuperClusterCollection> superCluster_;

      	const edm::InputTag ootSuperClusterCollectionTag;
     	edm::EDGetTokenT<reco::SuperClusterCollection> ootScToken_;
      	edm::Handle<reco::SuperClusterCollection> ootSuperCluster_;

		// calocluster
      	const edm::InputTag caloClusterTag;
      	edm::EDGetTokenT<vector<reco::CaloCluster>> ccToken_;
      	edm::Handle<vector<reco::CaloCluster>> caloCluster_;

      	// jets
      	const edm::InputTag jetsTag;
      	edm::EDGetTokenT<std::vector<reco::PFJet>> jetsToken_;
      	edm::Handle<std::vector<reco::PFJet>> jets_;

		float 				jetHt;
      	std::vector<float>  jetE, jetPt, jetPhi, jetEta, jetC, jetCharge; 
        std::vector<float>  jetEtaetaMmt, jetPhiphiMnt, jetEtaphiMnt, jetMaxD, jetConPtDis, jetConEtaPhiSprd, jetArea;
        std::vector<uInt>   jetNCarry, jetNConst;
      	std::vector<float>  jetMuTime, jetTimeError, jetTimeRMS, jetMedTime, jetCMuTime, jetCMedTime;
      	std::vector<float>  jetSCMuTime, jetSCMedTime, jetCSCMuTime, jetCSCMedTime, jetCBCMuTime, jetCBCMedTime;
      	std::vector<float>  jetPhMuTime, jetOOTPhMuTime, jetEleMuTime;
      	std::vector<int>    jetID, njetKids, jetKidOfJet, njetSubs, njetRecHits, jetRecHitOfJet;
      	std::vector<int>    jetKidPdgID, jetKidCharge, jetKid3Charge, jetPHM, jetELM;
      	std::vector<uInt>   jetRecHitId;
      	std::vector<bool>   jetKidLLP;
      	std::vector<float>  jetKidMass, jetKidVx, jetKidVy, jetKidVz;
      	std::vector<float>  jetKidE, jetKidPt, jetKidPhi, jetKidEta, jetKidTime, jetKidMedTime;
      	std::vector<float>  jetNHF, jetNEMF, jetCHF, jetCEMF, jetMUF, jetNHM, jetCHM, jetPHE, jetPHEF;
      	std::vector<float>  jetELE, jetELEF, jetMUE;

        std::vector<float>  jetSumEPFrac, jetEPEnergy, jetEMEnergy, jetEMEnrFrac, jetEPEnrFrac;
        std::vector<float>  jetDrLeadEta, jetDrLeadPhi, jetDrLeadEnr;
        std::vector<uInt>   jetDrRhCnt;
        std::vector<float>  sJetDrRHEnergy, jetDrEMF;
        std::vector<rhIdGroup> jetDrRhIds;

        std::vector<uInt>   nJetScMatch, jetScRhCnt;
        std::vector<float>  sJetScEnergy, sJetScPhEnergy, sJetScRhEnergy, jetScEMF;
        std::vector<rhIdGroup> jetScRhIds;

        std::vector<uInt>   jetBcTimesCnt, jetBcRhCnt, jetBcGrpCnt;
        std::vector<float>  jetBcSumRHEnr, jetBcEMFr; 
        std::vector<rhIdGroup> jetBcRhIds;

         // calojets
        const edm::InputTag caloJetsTag;
        edm::EDGetTokenT<std::vector<reco::CaloJet>> caloJetsToken_;
        edm::Handle<std::vector<reco::CaloJet>> caloJets_;

		std::vector<uInt>   cljBcCnt;
        std::vector<float>  cljPt, cljEnergy, cljPhi, cljEta, cljPx, cljPy, cljPz, cljEMFrac;
        std::vector<float>  cljSeedTOFTime, cljCMeanTime, cljCDrMeanTime;
        std::vector<float>  cljBc3dEx, cljBc3dEy, cljBc3dEz, cljBc3dEv, cljBc3dEslope, cljBc3dEchisp;
        std::vector<float>  cljBc2dEx, cljBc2dEy, cljBc2dEv, cljBc2dEslope, cljBc2dEchisp;
        std::vector<rhIdGroup> cljRhIds;

         // genJets
        const edm::InputTag genJetsTag;
        edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
        edm::Handle<std::vector<reco::GenJet>> genJets_;

        std::vector<float>  jetGenImpactAngle, jetGenTime, jetGenPt, jetGenEta, jetGenEnergy, jetGenEMFrac, jetGenDrMatch;
        std::vector<float>  jetGenTimeVar, jetGenTimeLLP, jetGenLLPPurity, jetGenNextBX, jetGenNKids, jetGenTOF;

		// JEC & JER  -> From miniAOD code not converted
  		// JECs
  		edm::ESHandle<JetCorrectorParametersCollection> jetCorrH;
  		// JERs
  		JME::JetResolution jetRes;
  		JME::JetResolutionScaleFactor jetRes_sf;

      	// electrons
      	const edm::InputTag electronsTag;
      	edm::EDGetTokenT<std::vector<reco::GsfElectron> > electronsToken_;
      	edm::Handle<std::vector<reco::GsfElectron> > electrons_;
      	std::vector<reco::GsfElectron> electrons;
     
        std::vector<int> 	eleCharge;
        std::vector<float>  elePt, eleEnergy, elePhi, eleEta, elePx, elePy, elePz;
        std::vector<float>  eleSeedTOFTime, eleCMeanTime;
        std::vector<float>  eleSc3dEx, eleSc3dEy, eleSc3dEz, eleSc3dEv, eleSc3dEslope, eleSc3dEchisp;
        std::vector<float>  eleSc2dEx, eleSc2dEy, eleSc2dEv, eleSc2dEslope, eleSc2dEchisp; 
        std::vector<rhIdGroup> eleRhIds;
        std::vector<float>  eleGenDr; 
        std::vector<int>    eleGenIdx;

      	// muons
      	const edm::InputTag muonsTag;
      	edm::EDGetTokenT<std::vector<reco::Muon> > muonsToken_;
      	edm::Handle<std::vector<reco::Muon> > muons_;
      	std::vector<reco::Muon> muons;
      
      	// RecHits
      	const edm::InputTag recHitsEBTag;
      	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEBToken_;
      	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEB_;
      	const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEB;

      	const edm::InputTag recHitsEETag;
      	edm::EDGetTokenT<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEEToken_;
      	edm::Handle<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > > recHitsEE_;
      	const edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit> > * recHitsEE;

        const edm::InputTag ecalBadCalibFlagTag;
        edm::EDGetTokenT<bool> ecalBadCalibFlagToken;
        edm::Handle<bool> ecalBadCalibFlagH;

      	std::vector<float> 	rhPosX, rhPosY, rhPosZ, rhPosEta, rhPosPhi;
      	std::vector<uInt> 	rhID;
		std::vector<int>    rhXtalI1, rhXtalI2, rhSubdet;
        std::vector<float>  rhEnergy, rhTime, rhTimeErr, rhTOF;  
      	std::vector<bool> 	rhisOOT, rhisGS6, rhisGS1, rhisWeird, rhisDiWeird;
      	std::vector<float> 	rhadcToGeV, rhSwCross;
      	std::vector<float> 	rhped12, rhped6, rhped1;
      	std::vector<float> 	rhpedrms12, rhpedrms6, rhpedrms1;

      	// gedPhotons
      	const edm::InputTag gedPhotonsTag;
      	edm::EDGetTokenT<std::vector<reco::Photon> > gedPhotonsToken_;
      	edm::Handle<std::vector<reco::Photon> > gedPhotons_;

        const edm::InputTag phoCBIDLooseMapTag;
        edm::EDGetTokenT<edm::ValueMap<bool>> phoCBIDLooseMapToken_;
        edm::Handle<edm::ValueMap<bool>> phoCBIDLooseMap_;

        std::vector<std::vector<uInt>> phoRhIds;

		std::vector<float>  phoPt, phoEnergy, phoPhi, phoEta, phoPx, phoPy, phoPz;
        std::vector<bool>   phoIsOotPho, phoExcluded;

        std::vector<float>  phoSeedTOFTime, phoCMeanTime;

        std::vector<bool>	phoIsPFPhoton, phoIsStdPhoton, phoHasConTracks, phoIsPixelSeed, phoIsPhoton, phoIsEB, phoIsEE;
        std::vector<bool>   phoIsEBGap, phoIsEBEtaGap, phoIsEBPhiGap, phoIsEEGap, phoIsEERingGap, phoIsEEDeeGap, phoIsEBEEGap;

        std::vector<float>	phoHadOverEM, phoHadOverEMVaid, phohadTowOverEM;
        std::vector<float>  phohadTowOverEMValid, phoE1x5, phoE2x5, phoE3x3, phoE5x5, phoMaxEnergyXtal; 
        std::vector<float>  phoSigmaEtaEta, phoSigmaIEtaIEta;
        std::vector<float>  phoR1x5, phoR2x5, phoR9, phoFull5x5_e1x5, phoFull5x5_e2x5, phoFull5x5_e3x3, phoFull5x5_e5x5; 
        std::vector<float>  phoFull5x5_maxEnergyXtal;
        std::vector<float>  phoFull5x5_sigmaEtaEta, phoFull5x5_sigmaIEtaIEta, phoFull5x5_r1x5, phoFull5x5_r2x5, phoFull5x5_r9;
        std::vector<int>    phoNSatXtals, phoMipNHitCone;
        std::vector<bool>   phoIsSeedSat, phoMipIsHalo;
        std::vector<float>  phoMipChi2, phoMipTotEnergy, phoMipSlope, phoMipInter;

        std::vector<float>  phoEcalRHSumEtConeDR04, phoHcalTwrSumEtConeDR04;
        std::vector<float>  phoHcalTowerSumEtBcConeDR04;
        std::vector<float>  phoTrkSumPtSolidConeDR04, phoTrkSumPtHollowConeDR04;
        std::vector<int>    phoNTrkSolidConeDR04, phoNTrkHollowConeDR04;

		std::vector<float>  phoSMaj, phoSMin, phoSAlp, phoCovEtaEta, phoCovEtaPhi, phoCovPhiPhi;
        std::vector<float>  phoGenDr; 
        std::vector<int> 	phoGenIdx;

      	// ootPhotons
      	const edm::InputTag ootPhotonsTag;
      	edm::EDGetTokenT<std::vector<reco::Photon> > ootPhotonsToken_;
      	edm::Handle<std::vector<reco::Photon> > ootPhotons_;
     
        // genEvtInfo
        const edm::InputTag genEvtInfoTag;
        edm::EDGetTokenT<GenEventInfoProduct> genEvtInfoToken_;
        edm::Handle<GenEventInfoProduct> genEvtInfo_;

        // gen time
        const edm::InputTag gent0Tag;
        edm::EDGetTokenT<float> gent0Token_;
        edm::Handle<float> gent0_;

        // gen vertex
        const edm::InputTag genxyz0Tag;
        edm::EDGetTokenT<Point3D> genxyz0Token_;
        edm::Handle<Point3D> genxyz0_;

        // pileups
        const edm::InputTag pileupInfosTag;
        edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupInfosToken_;
        edm::Handle<std::vector<PileupSummaryInfo> > pileupInfos_;

        // genParticles
        const edm::InputTag genParticlesTag;
        edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticlesToken_;
        edm::Handle<std::vector<reco::GenParticle> > genParticles_;

		uInt nGenParts;
        std::vector<float> genPt, genEnergy, genPhi, genEta, genPx, genPy, genPz;
        std::vector<int> genPdgId, genLLP;

      	// geometry CaloSubdetectorGeometry
        edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
      	edm::ESHandle<CaloGeometry> caloGeo_;
      	const CaloSubdetectorGeometry * barrelGeometry;
      	const CaloSubdetectorGeometry * endcapGeometry;  
		const CaloGeometry *ecalGeometry;

        // CaloTopology
        const edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
        edm::ESHandle<CaloTopology> caloTopo_;
        const CaloSubdetectorTopology * barrelTopology;
        const CaloSubdetectorTopology * endcapTopology;
		const CaloTopology * ecalTopology; 

        // lasers
        edm::ESGetToken<EcalLaserDbService, EcalLaserDbRecord> ecalLaserDbServiceToken_;
        edm::ESHandle<EcalLaserDbService> laser_;
        edm::Timestamp evTime;

        // inter calibration
        edm::ESGetToken<EcalIntercalibConstants, EcalIntercalibConstantsRcd> ecalIntercalibConstantsToken_;
        edm::ESHandle<EcalIntercalibConstants> interCalib_;
        const EcalIntercalibConstantMap * interCalibMap;

        // ADCToGeV
        edm::ESGetToken<EcalADCToGeVConstant, EcalADCToGeVConstantRcd> ecalADCToGeVConstantToken_;
        edm::ESHandle<EcalADCToGeVConstant> adcToGeV_;
        float adcToGeVEB;
        float adcToGeVEE;

        // pedestals
        edm::ESGetToken<EcalPedestals, EcalPedestalsRcd> EcalPedestalsToken_;
        edm::ESHandle<EcalPedestals> pedestals_;

	    // pf cand refs
		const edm::InputTag reco2pfTag;	
		const edm::EDGetTokenT<edm::ValueMap<std::vector<reco::PFCandidateRef> > > reco2pf_;	
		edm::Handle<edm::ValueMap<std::vector<reco::PFCandidateRef>>> reco2pf;

		///////////////////////////////////////////////////////////////////////////////////////
		// ------------- CMSSW Ana Class Helper Functions ------------------------------------ 
		///////////////////////////////////////////////////////////////////////////////////////

		// jet functions
        int getPFJetID(const reco::PFJet & jet);

		// sc functions
		float getSuperClusterSeedTime( reco::SuperClusterRef supercluster );
		float getPhotonSeedTime( pat::Photon photon );
		float getPhotonSeedTime( reco::Photon photon );

		// rechit group functions
        rhGroup getRHGroup( float eta, float phi, float drmin, float minenr );
        rhGroup getRHGroup( const scGroup superClusterGroup, float minenr );
        rhGroup getRHGroup( const scGroup superClusterGroup, float minenr, vector<float> phEnergy, vector<float> phDr, float phEnMax );
        rhGroup getRHGroup( const reco::CaloCluster basicCluster, float minenr );
        rhGroup getRHGroup( uInt detid );
        rhGroup getRHGroup();

        EcalRecHit getLeadRh( rhGroup recHits );
        vector<float> getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
        vector<float> getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );

        float getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ );
        float getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
        float getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ );
        void  mrgRhGrp( rhGroup & x, rhGroup & y );
        bool  reduceRhGrps( vector<rhGroup> & x );
		float setRecHitUsed( rhIdGroup idgroup );

		const auto getRawID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.rawId();}
		const auto getIsEB(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return (recHitId.subdetId() == EcalBarrel)?1:0;}
		const auto getSubDetID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.subdetId();}
		const auto rhMatch(const EcalRecHit rhx, const EcalRecHit rhy){ return getRawID(rhx) == getRawID(rhy);}
		const auto dupRhFnd(const rhGroup x, const rhGroup y){
	            		for(auto rhx : x ){ for(auto rhy : y ){if(rhMatch(rhx,rhy)){ return true;}}} return false;}
		const auto isRhGrpEx(const rhGroup x){ int s=x.size();for( int i=0;i<s;i++){
	                    for( int j=i+1;j<s;j++){if(rhMatch(x[i],x[j])) return false;}} return true;}
		const auto getRhGrpEnr(const rhGroup x){float e(0.0);for( auto ix : x ){e+=ix.energy();} return e;}
		const auto getDupCnt(const vector<rhGroup> x){ int c=0; int s=x.size();
                        for( int a=0;a<s;a++){for( int b=a+1;b<s;b++){if(dupRhFnd(x[a],x[b]))c++;}} return c;}
        const auto getRhGrpIDs(const rhGroup rhs ){ rhIdGroup rt; if(rhs.empty()){ rt.push_back(0);} else{ 
						for(const auto rh : rhs ){ rt.push_back(getRawID(rh));}} return rt;}


		// gen particle functions
		string bigKidChase( std::vector<reco::CandidatePtr> kids, float vx );
		void motherChase( const reco::Candidate* kid, string & depth );
		int llpGenChase( const reco::Candidate* kid );
		int llpGenChase( const reco::GenParticle kid );
		int llpGenChaseP( const reco::Candidate* kid, int depth );
		int llpGenChaseP( const reco::GenParticle kid, int depth );
		void kidChase( std::vector<reco::CandidatePtr> kids, float vx, float vy, float vz );
        vector<float> kidTOFChain( std::vector<reco::CandidatePtr> kids, float cx, float cy, float cz  );		
        vector<float>  getGenPartMatch( const reco::SuperCluster* scptr, std::vector<reco::GenParticle>  fgenparts );
    
};//<<>>class KUCMSNtupilizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    
#endif

//-------------------------------------------------------------------------------------------------------------------
