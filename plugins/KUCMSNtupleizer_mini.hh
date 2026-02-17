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
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

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
#include "KUCMSEventSelection.hh"

#ifndef KUCMSNtupilizerHeader
#define KUCMSNtupilizerHeader

using namespace std;
using namespace edm;

//#define NTHDEBUG true
#define NTHDEBUG false

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

typedef edm::View<pat::PackedCandidate> CandidateView;
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

class KUCMSNtupilizerMini : public edm::one::EDAnalyzer<edm::one::SharedResources>  {

    public:

        explicit KUCMSNtupilizerMini(const edm::ParameterSet&);
        ~KUCMSNtupilizerMini();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


    private:

        virtual void beginJob() override;
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endJob() override;

		// setup config tree
		void InitConfigTree( TTree* fConfigTree );

        //bool selectedEvent();
		//ntple event selection
		KUCMSEventSelection ntupleFilter;

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

		// config - event count and wgt info
		KUCMSBranchManager ConfigBranches;
    	int nTotEvts, nFltrdEvts, metFltrdEvts, phoFltrdEvts;
    	float sumEvtWgt, sumFltrdEvtWgt;

        // oputput tree
        TTree *outTree;
        TTree *configTree;

        // histograms
        TH1D *hist1d[nHists];
        TH2D *hist2d[nHists];

        ////////////////////////////////////////////////////////////
        // ----------- collection and branch declaration ----------
        ////////////////////////////////////////////////////////////

//        // Tracks ----------------------------------------------------------
//        
//       const edm::InputTag tracksTag;
//        edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
//        edm::Handle<std::vector<reco::Track>> tracks_;
//
//        // Trigger  ----------------------------------------------------------
//
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
//
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
//
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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   -------------------  Normally nothing needs to be modified below this point  --------------------------------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ Destructor -----------------------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////i

KUCMSNtupilizerMini::~KUCMSNtupilizerMini(){
 
    ///////////////////////////////////////////////////////////////////
    // do anything here that needs to be done at desctruction time   //
    // (e.g. close files, deallocate resources etc.)                 //
    ///////////////////////////////////////////////////////////////////
 
}//>>>>KUCMSNtupilizer::~KUCMSNtupilizer()

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ Analyzer Inherited Class Functions ------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void KUCMSNtupilizerMini::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){

    using namespace edm;

    // -- Consume Tokens --------------------------------------------
    // gets pointer to the collections from cmssw using the "token" for that collection

    if( NTHDEBUG ) std::cout << "Consume Tokens -------------------------------------------- " << std::endl;

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

    if( NTHDEBUG ) std::cout << "********************************************************************************************" << std::endl;

    // clear global event varibles 
    geVar.clear(); // floats

	// increment event counts and wgts
	int single = 1;
	geVar.set("genWgt",single);
    nTotEvts++;;

    // -------------------------------------------------------------------------------------------------
    // ---  Collecting objects ( preprocessing object pruning & fill global object vectors )------------
    // -------------------------------------------------------------------------------------------------

    if( NTHDEBUG ) std::cout << "LoadEvent ObjMan" << std::endl;
    ObjMan.LoadEvent( iEvent, iSetup, geVar );

    //------------------------------------------------------------------------------------
    // ----   Object processing ----------------------------------------------------------
    //------------------------------------------------------------------------------------
    // call functions to process collections and fill tree varibles to be saved
    // varibles to be saved to ttree are declared in the header
    // use LoadEvent() for any processing that must be done before crosstalk 
    // use PostProcessEvent() for any processing that must be done after crosstalk <<<  Most work should be done now.


    if( NTHDEBUG ) std::cout << "ProcessEvent ObjMan" << std::endl;
    ObjMan.ProcessEvent( geVar );
    ObjMan.PostProcessEvent( geVar );

    //------------------------------------------------------------------------------------
    //---- Object processing Completed ----------------------------------------------------------
    //------------------------------------------------------------------------------------

    // -- Fill output trees ------------------------------------------

    if( NTHDEBUG ) std::cout << "Select Event and Fill Tree" << std::endl;
    if( ntupleFilter.selectEvent( geVar ) ){
 
		outTree->Fill();
		nFltrdEvts++;
		sumFltrdEvtWgt += geVar("genWgt");

	}//<<>>if( ntupleSkim.selectEvent( geVar ) )

	metFltrdEvts += ntupleFilter.trackMetFilter() ? 1 : 0;
	phoFltrdEvts += ntupleFilter.trackPhoFilter() ? 1 : 0;
	sumEvtWgt += geVar("genWgt");

    // -- EOFun ------------------------------------------------------
    //     #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //     ESHandle<SetupData> pSetup;
    //     iSetup.get<SetupRecord>().get(pSetup);
    //     #endif

}//>>>>void KUCMSNtupilizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ------------ beginJob/endJob methods called once each job just before/after starting event loop    ------------
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void KUCMSNtupilizerMini::beginJob(){
 
    // Global Varibles
 
    // Book output files and trees
    edm::Service<TFileService> fs;
    outTree = fs->make<TTree>("llpgtree","KUCMSNtuple");
    configTree = fs->make<TTree>("configtree","KUCMSNtuple");
 
	InitConfigTree( configTree );

    // Book //histograms ( if any )
 
    std::cout << "Services Booked" << std::endl;
 
    if( NTHDEBUG ) std::cout << "Init ObjMan" << std::endl;
    ObjMan.Init( outTree );

    nTotEvts = 0; 
	nFltrdEvts = 0;
	sumEvtWgt = 0;
    sumFltrdEvtWgt = 0;
	metFltrdEvts = 0;
	phoFltrdEvts = 0;
 
}//>>>>void KUCMSNtupilizer::beginJob()

void KUCMSNtupilizerMini::InitConfigTree( TTree* fConfigTree ){

    ConfigBranches.makeBranch("nTotEvts","nTotEvts",INT);
    ConfigBranches.makeBranch("nFltrdEvts","nFltrdEvts",INT);
    ConfigBranches.makeBranch("sumEvtWgt","sumEvtWgt",FLOAT);
    ConfigBranches.makeBranch("sumFltrdEvtWgt","sumFltrdEvtWgt",FLOAT);
    ConfigBranches.makeBranch("nMetFltrdEvts","nMetFltrdEvts",INT);
    ConfigBranches.makeBranch("nPhoFltrdEvts","nPhoFltrdEvts",INT);

    ConfigBranches.attachBranches(fConfigTree);

}//<<>>void KUCMSNtupilizer::InitConfigTree( TTree* fConfigTree )

// ------------ method called once each job just after ending the event loop    ------------
void KUCMSNtupilizerMini::endJob(){
 
    if( NTHDEBUG ) std::cout << "ObjMan EndJobs" << std::endl;
    ObjMan.EndJobs();

	ConfigBranches.clearBranches();
	ConfigBranches.fillBranch("nTotEvts",nTotEvts);
	ConfigBranches.fillBranch("nFltrdEvts",nFltrdEvts);
	ConfigBranches.fillBranch("sumEvtWgt",sumEvtWgt);
    ConfigBranches.fillBranch("sumFltrdEvtWgt",sumFltrdEvtWgt);
    ConfigBranches.fillBranch("nMetFltrdEvts",metFltrdEvts);
    ConfigBranches.fillBranch("nPhoFltrdEvts",phoFltrdEvts);
	configTree->Fill();
 
}//>>>>void KUCMSNtupilizer::endJob()

// ------------ method fills 'descriptions' with the allowed parameters for the module    ------------
void KUCMSNtupilizerMini::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
 
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
    
#endif

//-------------------------------------------------------------------------------------------------------------------
