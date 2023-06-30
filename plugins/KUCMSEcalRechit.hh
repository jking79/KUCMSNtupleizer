// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Object for Object Manager
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

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

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

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
#include "CommonTools/Egamma/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/EcalRegressionData.h"

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


// KUCMS
#include "KUCMSObjectBase.hh"

#ifndef KUCMSEcalRecHitObjectHeader
#define KUCMSEcalRecHitObjectHeader

//#define ERHODEBUG true
#define ERHODEBUG false

using namespace edm;

typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> recHitCol;
typedef std::vector<EcalRecHit> rhGroup;
typedef std::vector<reco::SuperCluster> scGroup;
typedef std::vector<reco::CaloCluster> bcGroup;
typedef std::vector<uInt> rhIdGroup;

enum class ECAL {EB, EM, EP, EE, NONE};

//class exampleObject : public KUCMSObjectBase {};

class KUCMSEcalRecHitObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    //KUCMSEcalRecHitObject();
    KUCMSEcalRecHitObject( const edm::ParameterSet& iConfig );
	~KUCMSEcalRecHitObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
	//void ConfigObject( const edm::ParameterSet& iConfig );
	void LoadRecHitTokens( edm::EDGetTokenT<recHitCol>* recHitsEBToken_, edm::EDGetTokenT<recHitCol>* recHitsEEToken_ );
  void LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection>* scToken_, edm::EDGetTokenT<reco::SuperClusterCollection>* ootScToken_ );
    void LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>>* ccToken_ );
    void InitObject( TTree* fOutTree ); // sets up branches, do preloop jobs
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ); 
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar ); 

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...
/*
    // sc functions
    float getSuperClusterSeedTime( reco::SuperClusterRef supercluster );
    //float getPhotonSeedTime( pat::Photon photon );
    //float getPhotonSeedTime( reco::Photon photon );

    // rechit group functions
    rhGroup getRHGroup( float eta, float phi, float drmin, float minenr );
    rhGroup getRHGroup( const scGroup superClusterGroup, float minenr );
    rhGroup getRHGroup( const scGroup superClusterGroup, float minenr, std::vector<float> phEnergy, std::vector<float> phDr, float phEnMax );
    rhGroup getRHGroup( const reco::CaloCluster basicCluster, float minenr );
    rhGroup getRHGroup( uInt detid );
    rhGroup getRHGroup();

    EcalRecHit getLeadRh( rhGroup recHits );
    std::vector<float> getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
    std::vector<float> getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );

    float getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ );
    float getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
    float getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ );
    void  mrgRhGrp( rhGroup & x, rhGroup & y );
    bool  reduceRhGrps( std::vector<rhGroup> & x );
    void  setRecHitUsed( rhIdGroup idgroup );

*/
    const auto getRawID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.rawId();}
    const auto getIsEB(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return (recHitId.subdetId() == EcalBarrel)?1:0;}
    const auto getSubDetID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.subdetId();}
    const auto rhMatch(const EcalRecHit rhx, const EcalRecHit rhy){ return getRawID(rhx) == getRawID(rhy);}
    const auto dupRhFnd(const rhGroup x, const rhGroup y){
               for(auto rhx : x ){ for(auto rhy : y ){if(rhMatch(rhx,rhy)){ return true;}}} return false;}
    const auto isRhGrpEx(const rhGroup x){ int s=x.size();for( int i=0;i<s;i++){
               for( int j=i+1;j<s;j++){if(rhMatch(x[i],x[j])) return false;}} return true;}
    const auto getRhGrpEnr(const rhGroup x){float e(0.0);for( auto ix : x ){e+=ix.energy();} return e;}
    const auto getDupCnt(const std::vector<rhGroup> x){ int c=0; int s=x.size();
               for( int a=0;a<s;a++){for( int b=a+1;b<s;b++){if(dupRhFnd(x[a],x[b]))c++;}} return c;}
    const auto getRhGrpIDs(const rhGroup rhs ){ rhIdGroup rt; if(rhs.empty()){ rt.push_back(0);} else{
               for(const auto rh : rhs ){ rt.push_back(getRawID(rh));}} return rt;}

    private:

    std::vector<EcalRecHit> frechits;
    std::vector<bool> frhused;

    //const edm::InputTag recHitsEBTag;
    edm::EDGetTokenT<recHitCol>* recHitsEBToken_;
    edm::Handle<recHitCol> recHitsEB_;
    //const recHitCol* recHitsEB;

    //const edm::InputTag recHitsEETag;
    edm::EDGetTokenT<recHitCol>* recHitsEEToken_;
    edm::Handle<recHitCol> recHitsEE_;
    //const recHitCol* recHitsEE;

    const edm::InputTag ecalBadCalibFlagTag;
    edm::EDGetTokenT<bool>* ecalBadCalibFlagToken;
    edm::Handle<bool> ecalBadCalibFlagH;

    // supercluster
    const edm::InputTag superClusterCollectionTag;
    edm::EDGetTokenT<reco::SuperClusterCollection>* scToken_;
    edm::Handle<reco::SuperClusterCollection> superCluster_;

    const edm::InputTag ootSuperClusterCollectionTag;
    edm::EDGetTokenT<reco::SuperClusterCollection>* ootScToken_;
    edm::Handle<reco::SuperClusterCollection> ootSuperCluster_;

    // calocluster
    const edm::InputTag caloClusterTag;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>>* ccToken_;
    edm::Handle<std::vector<reco::CaloCluster>> caloCluster_;

    // geometry CaloSubdetectorGeometry
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeo_;
    const CaloSubdetectorGeometry * barrelGeometry;
    const CaloSubdetectorGeometry * endcapGeometry;
    const CaloGeometry *ecalGeometry;

    // CaloTopology
    edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
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

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    //std::reference_wrapper<exampleObject> otherObjectRef;

};//<<>>class KUCMSEcalRecHit : public KUCMSObjectBase

//KUCMSEcalRecHitObject::KUCMSEcalRecHitObject( const edm::ParameterSet& iConfig ){ ConfigObject( iConfig ); } 

//void KUCMSEcalRecHitObject::ConfigObject( const edm::ParameterSet& iConfig ){

KUCMSEcalRecHitObject::KUCMSEcalRecHitObject( const edm::ParameterSet& iConfig ):

    // recHits
    //recHitsEBTag(iConfig.getParameter<edm::InputTag>("recHitsEB")),
    //recHitsEETag(iConfig.getParameter<edm::InputTag>("recHitsEE")),

    // supercluster
    //superClusterCollectionTag(iConfig.getParameter<edm::InputTag>("superClusters")),
    //ootSuperClusterCollectionTag(iConfig.getParameter<edm::InputTag>("ootSuperClusters")),

    // caloclusters
    //caloClusterTag(iConfig.getParameter<edm::InputTag>("caloClusters")),

    // ECAL RECORDS 
    caloGeometryToken_(esConsumes()),
    caloTopologyToken_{esConsumes()},
    ecalLaserDbServiceToken_(esConsumes()),
    ecalIntercalibConstantsToken_(esConsumes()),
    ecalADCToGeVConstantToken_(esConsumes()),
    EcalPedestalsToken_(esConsumes())

{   // ---- end constructor initilizations  --------------------------

    // supercluster
    //if( ERHODEBUG ) std::cout << "Getting supercluster token" << std::endl;
    //scToken_ = IC->consumes<reco::SuperClusterCollection>(superClusterCollectionTag);
    //ootScToken_ = IC->consumes<reco::SuperClusterCollection>(ootSuperClusterCollectionTag);

    // caloClusters
    //ccToken_ = IC->consumes<std::vector<reco::CaloCluster>>(caloClusterTag);

    // rechits
    //if( ERHODEBUG ) std::cout << "Getting rechit token" << std::endl;
    //recHitsEBToken_ = IC->consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEBTag);
    //recHitsEEToken_ = IC->consumes<edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>>>(recHitsEETag);

    cfPrm.set( "minRHEi", iConfig.existsAs<double>("minRHEi") ? iConfig.getParameter<double>("minRHEi") : 0.0 );
    cfPrm.set( "minRHEf", iConfig.existsAs<double>("minRHEf") ? iConfig.getParameter<double>("minRHEf") : 2.0 );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : true );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );


}//<<>>KUCMSEcalRecHit::KUCMSEcalRecHit( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSEcalRecHitObject::LoadRecHitTokens( edm::EDGetTokenT<recHitCol>* recHitsEBToken, edm::EDGetTokenT<recHitCol>* recHitsEEToken ){

	recHitsEBToken_ = recHitsEBToken;
	recHitsEEToken_ = recHitsEEToken;

}//<<>>void LoadRecHitTokens( edm::EDGetTokenT<recHitCol>* recHitsEBToken_, edm::EDGetTokenT<recHitCol>* recHitsEEToken_ )

void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection>* scToken, edm::EDGetTokenT<reco::SuperClusterCollection>* ootScToken ){

	scToken_ = scToken;
	ootScToken_ = ootScToken;

}//<<>>void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection>*

void KUCMSEcalRecHitObject::LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>>* ccToken ){

	ccToken_ = ccToken;

}//<<>>void KUCMSEcalRecHitObject::LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>>* ccToken )

void KUCMSEcalRecHitObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("Energy","omERH_energy",VFLOAT);
    Branches.makeBranch("Time","omERH_time",VFLOAT);
    Branches.makeBranch("TOF","omERH_TOF",VFLOAT);
    Branches.makeBranch("ID","omERH_ID",VUINT);
    Branches.makeBranch("isOOT","omERH_isOOT",VBOOL);
    Branches.makeBranch("SwCross","omERH_swCross",VFLOAT);
    Branches.makeBranch("eta","omERH_eta",VFLOAT);
    Branches.makeBranch("phi","omERH_phi",VFLOAT);
    Branches.makeBranch("rhx","omERH_rhx",VFLOAT);
    Branches.makeBranch("rhy","omERH_rhy",VFLOAT);
    Branches.makeBranch("rhz","omERH_rhz",VFLOAT);

    Branches.attachBranches(fOutTree);	

}//<<>>void KUCMSEcalRecHit::InitObject( TTree* fOutTree )

void KUCMSEcalRecHitObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS
    iEvent.getByToken( *recHitsEBToken_, recHitsEB_);
    iEvent.getByToken( *recHitsEEToken_, recHitsEE_);

    // SUPERCLUSTERS
    iEvent.getByToken( *scToken_, superCluster_);
    iEvent.getByToken( *ootScToken_, ootSuperCluster_);

    // CALOCLUSTERS
    iEvent.getByToken( *ccToken_, caloCluster_);

    // GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
    caloGeo_ = iSetup.getHandle(caloGeometryToken_);
    barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
    endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap);
    ecalGeometry = &iSetup.getData(caloGeometryToken_);

    // CaloTopology
    caloTopo_ = iSetup.getHandle(caloTopologyToken_);
    ecalTopology = &iSetup.getData(caloTopologyToken_);
    barrelTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalBarrel);
    endcapTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalEndcap);

    // Laser constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
    laser_ = iSetup.getHandle(ecalLaserDbServiceToken_);
    evTime = iEvent.time();

    // Intercalibration constants : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
    interCalib_ = iSetup.getHandle(ecalIntercalibConstantsToken_);
    interCalibMap = &interCalib_->getMap();

    // ADCToGeV : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h
    adcToGeV_ = iSetup.getHandle(ecalADCToGeVConstantToken_);
    adcToGeVEB = adcToGeV_->getEBValue();
    adcToGeVEE = adcToGeV_->getEEValue();

    // Pedestals : https://github.com/ferriff/usercode/blob/master/DBDump/plugins/DBDump.cc
    pedestals_ = iSetup.getHandle(EcalPedestalsToken_);

    if( ERHODEBUG ) std::cout << "Collecting ECAL RecHits" << std::endl;
    frechits.clear();
    frhused.clear();
    for( const auto &recHit : *recHitsEB_ ){ if(recHit.energy() > cfPrm("minRHEi")) frechits.push_back(recHit); frhused.push_back(false);}
    for( const auto &recHit : *recHitsEE_ ){ if(recHit.energy() > cfPrm("minRHEi")) frechits.push_back(recHit); frhused.push_back(false);}

}//<<>>void KUCMSEcalRecHit::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )


void KUCMSEcalRecHitObject::ProcessEvent( ItemManager<float>& geVar ){

    if( ERHODEBUG ) std::cout << "Processing RecHits" << std::endl;

    Branches.clearBranches();

    if( ERHODEBUG ) std::cout << " - enetering RecHit loop" << std::endl;

    int nRecHits = frechits.size();
    for ( int it = 0; it < nRecHits; it++ ){

        auto recHit = frechits[it];
        auto used = frhused[it];
        //if( not used ) continue;
        if( not used && recHit.energy() < cfPrm("minRHEf") ) continue;

        const auto recHitID = getRawID(recHit);
        const bool isEB = getIsEB(recHit); // which subdet
        const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        const auto recHitPos = geometry->getGeometry(recHit.detid())->getPosition();
        const float eta = recHitPos.eta();
        const float phi = recHitPos.phi();
        const auto rhX = recHitPos.x();
        const auto rhY = recHitPos.y();
        const auto rhZ = recHitPos.z();
        const float d_rh = hypo(rhX,rhY,rhZ)/SOL;
        const float d_pv = hypo(rhX-geVar("vtxX"),rhY-geVar("vtxY"),rhZ-geVar("vtxZ"))/SOL;
        //float swisscross(0.0);
        //if( isEB ) swisscross = EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true);
        //else swisscross = EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);
        const float swisscross = (isEB) ?
                EcalTools::swissCross(recHitID, *recHitsEB_, 0.0, true) : EcalTools::swissCross(recHitID, *recHitsEE_, 0.0, true);


        //const auto laser = laser_->getLaserCorrection(recHitID,evTime);
        //const auto interCalibIter = interCalibMap->find(recHitID);
        //const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        //const float adcToGeV( isEB ? adcToGeVEB : adcToGeVEE );
        //const auto & pediter = pedestals_->find(recHitID);

        const float rhTime = recHit.time();
        const bool rhIsOOT = recHit.checkFlag(EcalRecHit::kOutOfTime);
        const float rhEnergy = recHit.energy();
        const float rhAdjTime = rhTime-d_rh;  // Note : Margret adds d_rh to the time in her code & subtracts d_pv ( or TOF )

        Branches.fillBranch("ID",recHitID);
        Branches.fillBranch("TOF",d_pv);
        Branches.fillBranch("Time",rhAdjTime);
        Branches.fillBranch("isOOT",rhIsOOT);
        Branches.fillBranch("Energy",rhEnergy);
        Branches.fillBranch("SwCross",swisscross);
        Branches.fillBranch("eta",eta);
        Branches.fillBranch("phi",phi);
        Branches.fillBranch("rhx",rhX);
        Branches.fillBranch("rhy",rhY);
        Branches.fillBranch("rhz",rhZ);
        //Branches.fillBranch("rhisWeird",recHit.checkFlag(EcalRecHit::kWeird));
        //Branches.fillBranch("rhisDiWeird",recHit.checkFlag(EcalRecHit::kDiWeird));

    }//<<>>for (const auto recHit : *recHitsEB_ ) 

}//<<>>void KUCMSEcalRecHit::ProcessEvent()

void KUCMSEcalRecHitObject::EndJobs(){}

#endif
