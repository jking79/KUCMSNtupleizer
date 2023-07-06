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

// Supercluster info
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

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

#include "KUCMSItemManager.hh"
#include "KUCMSBranchManager.hh"
#include "KUCMSObjectBase.hh"

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

//const auto sortByPt = [](auto & obj1, auto & obj2) {return obj1.pt() > obj2.pt();};

//
//  constants, enums and typedefs
//

using reco::TrackCollection;

typedef edm::View<reco::Candidate> CandidateView;
typedef vector<reco::SuperCluster> scGroup;
typedef vector<reco::CaloCluster> bcGroup;

typedef unsigned int uInt;
typedef unsigned long int ulInt;
typedef const float CFlt;
typedef const double CDbl;
typedef const vector<float> CVFlt;
typedef const vector<double> CVDbl;

typedef vector<uInt> rhIdGroup;
typedef ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<float>,ROOT::Math::DefaultCoordinateSystemTag> Point3D;

#define SOL 29.9792458 // speed of light in cm/ns
#define PI 3.14159265358979323846 // pie ...  
#define BUNCHES 3564
#define SAMPLES 10
#define nHists 10

//
//  Class Declaration
//

class KUCMSNtupilizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

    public:

        explicit KUCMSNtupilizer(const edm::ParameterSet&);
        ~KUCMSNtupilizer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

		friend class KUCMSEcalRechit;

    private:

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

        bool selectedEvent();

        void setBranches();
        void processEvent( const edm::Event& iEvent );

        ////////////////////////////////////////////////////
        // ----------member data ---------------------------
        ////////////////////////////////////////////////////
        //
        // input parameters
        ItemManager<double> cfPrm;
        ItemManager<int> cfCnt;
        ItemManager<std::string> cfStr;
        ItemManager<bool> cfFlag;

		//KUCMSEcalRecHit recHitsObj; 

		ObjectManager ObjMan;

        // global event varibles
        ItemManager<float> geVar;

        // oputput tree
        TTree *outTree;

        // histograms
        TH1D *hist1d[nHists];
        TH2D *hist2d[nHists];

        ////////////////////////////////////////////////////////////
        // ----------- collection and branch declaration ----------
        ////////////////////////////////////////////////////////////

        // Event ----------------------------------------------------------

        KUCMSBranchManager EventInfoBM;

        // Tracks ----------------------------------------------------------
        
        //const edm::InputTag tracksTag;
        //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
        //edm::Handle<std::vector<reco::Track>> tracks_;

        // Trigger  ----------------------------------------------------------

//        // Trigger Info -> From miniAOD code not converted
//        // trigger inputs
//        const std::string inputPaths;
//        std::vector<std::string> pathNames;
//        //strBitMap triggerBitMap;
//        const std::string inputFilters;
//        std::vector<std::string> filterNames;
//        // trigger results
//        const edm::InputTag triggerResultsTag;
//        edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
//        edm::Handle<edm::TriggerResults> triggerResults_;
//        // trigger objects
//        const edm::InputTag triggerObjectsTag;
//        edm::EDGetTokenT<std::vector<pat::TriggerObjectStandAlone>> triggerObjectsToken_;
//        edm::Handle<std::vector<pat::TriggerObjectStandAlone>> triggerObjects_;
//        // output triggers
//        std::map<std::string,std::vector<pat::TriggerObjectStandAlone>> triggerObjectsByFilterMap; 
//        // first index is filter label, second is trigger objects

        // Vertices ----------------------------------------------------------

        //const edm::InputTag verticesTag;
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
        edm::Handle<std::vector<reco::Vertex>> vertices_;

//        // Met Filters  ----------------------------------------------------------
//
//        // Met Filters -> From miniAOD code not converted
//        // met filter inputs
//        const std::string inputFlags;
//        std::vector<std::string> flagNames;
//        //strBitMap triggerFlagMap;
//        // met filters
//        const edm::InputTag triggerFlagsTag;
//        edm::EDGetTokenT<edm::TriggerResults> triggerFlagsToken;
//        edm::Handle<edm::TriggerResults> triggerFlagsH;

//        // Jet corrections  ----------------------------------------------------------
//
//        // JEC & JER  -> From miniAOD code not converted
//        // JECs
//        edm::ESHandle<JetCorrectorParametersCollection> jetCorrH;
//        // JERs
//        JME::JetResolution jetRes;
//        JME::JetResolutionScaleFactor jetRes_sf;

        // Setup info  ----------------------------------------------------------

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

};//<<>>class KUCMSNtupilizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
    
#endif

//-------------------------------------------------------------------------------------------------------------------
