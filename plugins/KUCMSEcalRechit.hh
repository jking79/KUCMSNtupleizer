// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS ECAL RecHit Object
// This object process information for ECAL rechits and SuperClusters   
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

// Specific includes for this Object

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


// KUCMS Object includes

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

typedef reco::SuperClusterCollection supClusterCol;

class KUCMSEcalRecHitObject : public KUCMSObjectBase {

    public:

    // use base class constructor
    KUCMSEcalRecHitObject( const edm::ParameterSet& iConfig );
	~KUCMSEcalRecHitObject(){};

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
	// load tokens for eventt based collections
	void LoadRecHitTokens( edm::EDGetTokenT<recHitCol> recHitsEBToken_, edm::EDGetTokenT<recHitCol> recHitsEEToken_ );
    void LoadSCTokens( edm::EDGetTokenT<supClusterCol> scToken_, edm::EDGetTokenT<supClusterCol> ootScToken_ );
    void LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>> ccToken_ );
    void LoadSCTokens( edm::EDGetTokenT<supClusterCol> otherscToken_ );
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ); 
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar ){}; 
	void PostProcessEvent( ItemManager<float>& geVar );

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    // sc functions
    float getSuperClusterSeedTime( reco::SuperClusterRef supercluster );
	int getSuperClusterIndex( const reco::SuperCluster* supercluster );

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

	// Helper functions for rechit infomation & IDs
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
    const auto getOverLapCnt(const rhIdGroup x, const rhIdGroup y){ int c=0; int sx=x.size(); int sy=y.size();
               for( int a=0;a<sx;a++){for( int b=0;b<sy;b++){if(x[a]==y[b]){c++;break;}}} return c;}

	// Wrapper functions for EcalClusterTools functions 
	Cluster2ndMoments getCluster2ndMoments( const reco::SuperCluster* scptr );
    std::array<float, 3> getCovariances( const reco::SuperCluster* scptr );

    private:

    std::vector<EcalRecHit> frechits;
    std::vector<bool> frhused;

    std::vector<reco::SuperCluster> fsupclstrs;
    std::vector<bool> fscExclude;
    std::vector<bool> fscIsOOT;

    //const edm::InputTag recHitsEBTag;
    edm::EDGetTokenT<recHitCol> recHitsEBToken_;
    edm::Handle<recHitCol> recHitsEB_;
    //const recHitCol* recHitsEB;

    //const edm::InputTag recHitsEETag;
    edm::EDGetTokenT<recHitCol> recHitsEEToken_;
    edm::Handle<recHitCol> recHitsEE_;
    //const recHitCol* recHitsEE;

    const edm::InputTag ecalBadCalibFlagTag;
    edm::EDGetTokenT<bool> ecalBadCalibFlagToken;
    edm::Handle<bool> ecalBadCalibFlagH;

    // supercluster
    const edm::InputTag superClusterCollectionTag;
    edm::EDGetTokenT<reco::SuperClusterCollection> scToken_;
    edm::Handle<reco::SuperClusterCollection> superCluster_;

    const edm::InputTag ootSuperClusterCollectionTag;
    edm::EDGetTokenT<reco::SuperClusterCollection> ootScToken_;
    edm::Handle<reco::SuperClusterCollection> ootSuperCluster_;

	int nOtherSC;
    int nOtherSCEx;
    const edm::InputTag otherSuperClusterCollectionTag;
    edm::EDGetTokenT<reco::SuperClusterCollection> otherScToken_;
    edm::Handle<reco::SuperClusterCollection> otherSuperCluster_;

    // calocluster
    const edm::InputTag caloClusterTag;
    edm::EDGetTokenT<std::vector<reco::CaloCluster>> ccToken_;
    edm::Handle<std::vector<reco::CaloCluster>> caloCluster_;

    // geometry CaloSubdetectorGeometry
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeo_;
    const CaloSubdetectorGeometry * barrelGeometry;
    const CaloSubdetectorGeometry * endcapGeometry;
    const CaloGeometry * ecalGeometry;

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

KUCMSEcalRecHitObject::KUCMSEcalRecHitObject( const edm::ParameterSet& iConfig ):

    // ECAL RECORDS 
    caloGeometryToken_(esConsumes()),
    caloTopologyToken_{esConsumes()},
    ecalLaserDbServiceToken_(esConsumes()),
    ecalIntercalibConstantsToken_(esConsumes()),
    ecalADCToGeVConstantToken_(esConsumes()),
    EcalPedestalsToken_(esConsumes())

{   // ---- end constructor initilizations  --------------------------

    cfPrm.set( "minRHEi", iConfig.existsAs<double>("minRHEi") ? iConfig.getParameter<double>("minRHEi") : 0.0 );
    cfPrm.set( "minRHEf", iConfig.existsAs<double>("minRHEf") ? iConfig.getParameter<double>("minRHEf") : 0.0 );
    cfPrm.set( "minSCE", iConfig.existsAs<double>("minSCE") ? iConfig.getParameter<double>("minSCE") : 0.0 );
    cfFlag.set( "onlyEB", iConfig.existsAs<bool>("onlyEB") ? iConfig.getParameter<bool>("onlyEB") : false );
    cfPrm.set( "ebMaxEta",iConfig.existsAs<double>("ebMaxEta")? iConfig.getParameter<double>("ebMaxEta") : 1.479 );


}//<<>>KUCMSEcalRecHit::KUCMSEcalRecHit( const edm::ParameterSet& iConfig, const ItemManager<bool>& cfFlag )

void KUCMSEcalRecHitObject::LoadRecHitTokens( edm::EDGetTokenT<recHitCol> recHitsEBToken, edm::EDGetTokenT<recHitCol> recHitsEEToken ){

	recHitsEBToken_ = recHitsEBToken;
	recHitsEEToken_ = recHitsEEToken;

}//<<>>void LoadRecHitTokens( edm::EDGetTokenT<recHitCol>* recHitsEBToken_, edm::EDGetTokenT<recHitCol>* recHitsEEToken_ )

void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection> scToken, edm::EDGetTokenT<reco::SuperClusterCollection> ootScToken ){

	scToken_ = scToken;
	ootScToken_ = ootScToken;

}//<<>>void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection>*

void KUCMSEcalRecHitObject::LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>> ccToken ){

	ccToken_ = ccToken;

}//<<>>void KUCMSEcalRecHitObject::LoadClusterTokens( edm::EDGetTokenT<std::vector<reco::CaloCluster>>* ccToken )

void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection> otherScToken ){

    otherScToken_ = otherScToken;

}//<<>>void KUCMSEcalRecHitObject::LoadSCTokens( edm::EDGetTokenT<reco::SuperClusterCollection>*

void KUCMSEcalRecHitObject::InitObject( TTree* fOutTree ){

    Branches.makeBranch("Energy","ECALRecHit_energy",VFLOAT);
    Branches.makeBranch("Time","ECALRecHit_time",VFLOAT);
    Branches.makeBranch("TOFpv","ECALRecHit_pvTOF",VFLOAT);
    Branches.makeBranch("TOF0","ECALRecHit_0TOF",VFLOAT);
    Branches.makeBranch("ID","ECALRecHit_ID",VUINT);
    Branches.makeBranch("isOOT","ECALRecHit_isOOT",VBOOL);
    Branches.makeBranch("SwCross","ECALRecHit_swCross",VFLOAT);
    Branches.makeBranch("eta","ECALRecHit_eta",VFLOAT);
    Branches.makeBranch("phi","ECALRecHit_phi",VFLOAT);
    Branches.makeBranch("rhx","ECALRecHit_rhx",VFLOAT);
    Branches.makeBranch("rhy","ECALRecHit_rhy",VFLOAT);
    Branches.makeBranch("rhz","ECALRecHit_rhz",VFLOAT);
    Branches.makeBranch("amplitude","ECALRecHit_amplitude",VFLOAT);
    Branches.makeBranch("ampres","ECALRecHit_ampres",VFLOAT);

    Branches.makeBranch("zscnSC","SuperCluster_nSuperCluster",INT);
    Branches.makeBranch("zscOOT","SuperCluster_isOot",VBOOL);
    Branches.makeBranch("zscExcluded","SuperCluster_excluded",VBOOL);
    Branches.makeBranch("zscsid","SuperCluster_XtalSeedID",VUINT);
    Branches.makeBranch("zsclcx","SuperCluster_clcx",VFLOAT,"x coordinate of cluster centroid");
    Branches.makeBranch("zsclcy","SuperCluster_clcy",VFLOAT,"y coordinate of cluster centroid");
    Branches.makeBranch("zsclcz","SuperCluster_clcz",VFLOAT,"z coordinate of cluster centroid");
    Branches.makeBranch("zsceta","SuperCluster_eta",VFLOAT,"pseudorapidity of cluster centroid");
    Branches.makeBranch("zscphi","SuperCluster_phi",VFLOAT,"azimuthal angle of cluster centroid");
    Branches.makeBranch("zscenergy","SuperCluster_energy",VFLOAT);
    Branches.makeBranch("zscsize","SuperCluster_nXtals",VUINT);
    Branches.makeBranch("zscrhids","SuperCluster_rhIds",VVUINT,"list of rechit raw ids in hits and fractions list from supercluster");

    Branches.makeBranch("zscEnergyRaw","SuperCluster_energyRaw",VFLOAT,"raw energy of photon supercluster");//
    Branches.makeBranch("zscIsEB","SuperCluster_seedIsEB",VBOOL,"photon supercluster seed crystal is in ecal barrel");

    Branches.makeBranch("zscIsScEtaEB","SuperCluster_isScEtaEB",VBOOL,"is supercluster eta within barrel acceptance");//
    Branches.makeBranch("zscIsScEtaEE","SuperCluster_isScEtaEE",VBOOL,"is supercluster eta within endcap acceptance");//
    Branches.makeBranch("zscSeediEtaOriX","SuperCluster_seediEtaOriX",VINT,"iEta or iX of seed crystal. iEta is barrel-only, iX is endcap-only. iEta runs from -85 to +85, with no crystal at iEta=0. iX runs from 1 to 100.");//
    Branches.makeBranch("zscSeediPhiOriY","SuperCluster_seediPhiOriY",VINT,"iPhi or iY of seed crystal. iPhi is barrel-only, iY is endcap-only. iPhi runs from 1 to 360. iY runs from 1 to 100.");//
    Branches.makeBranch("zscx_calo","SuperCluster_x_calo",VFLOAT,"photon supercluster position on calorimeter, x coordinate (cm)");//
    Branches.makeBranch("zscy_calo","SuperCluster_y_calo",VFLOAT,"photon supercluster position on calorimeter, y coordinate (cm)");//
    Branches.makeBranch("zscz_calo","SuperCluster_z_calo",VFLOAT,"photon supercluster position on calorimeter, z coordinate (cm)");//

    Branches.makeBranch("zscEtaWidth","SuperCluster_etaWidth",VFLOAT,"Width of the photon supercluster in eta");//
    Branches.makeBranch("zscPhiWidth","SuperCluster_phiWidth",VFLOAT,"Width of the photon supercluster in phi");//
    Branches.makeBranch("zscSMaj","SuperCluster_smaj",VFLOAT);
    Branches.makeBranch("zscSMin","SuperCluster_smin",VFLOAT);
    Branches.makeBranch("zscSAlp","SuperCluster_salp",VFLOAT);
    Branches.makeBranch("zscCovEtaEta","SuperCluster_covEtaEta",VFLOAT);
    Branches.makeBranch("zscCovEtaPhi","SuperCluster_covEtaPhi",VFLOAT);
    Branches.makeBranch("zscCovPhiPhi","SuperCluster_covPhiPhi",VFLOAT);

    Branches.makeBranch("zscnOther","SuperCluster_nOther",INT);
    Branches.makeBranch("zscnOtherEx","SuperCluster_nOtherEx",INT);
    Branches.makeBranch("zscnOtherIn","SuperCluster_nOtherIn",INT);
    Branches.makeBranch("zscnOExDr","SuperCluster_nOExDr",VFLOAT);
    Branches.makeBranch("zscnOver","SuperCluster_nXtalOverlap",VINT);
    Branches.makeBranch("zscnOtherSID","SuperCluster_otherSeedID",VUINT);
    Branches.makeBranch("zscnOMatchSID","SuperCluster_otherMatchSeedID",VUINT);

    Branches.attachBranches(fOutTree);	

}//<<>>void KUCMSEcalRecHit::InitObject( TTree* fOutTree )

void KUCMSEcalRecHitObject::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ){

    // ECAL RECHITS
    iEvent.getByToken( recHitsEBToken_, recHitsEB_);
    iEvent.getByToken( recHitsEEToken_, recHitsEE_);

    // SUPERCLUSTERS
    iEvent.getByToken( scToken_, superCluster_);
    iEvent.getByToken( ootScToken_, ootSuperCluster_);
    iEvent.getByToken( otherScToken_, otherSuperCluster_);

    // CALOCLUSTERS
    iEvent.getByToken( ccToken_, caloCluster_);

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

    if( ERHODEBUG ) std::cout << "Collecting SuperClusters" << std::endl;
    fsupclstrs.clear();
    fscExclude.clear();
    fscIsOOT.clear();
    for( const auto &supclstr : *ootSuperCluster_ ){ 
		if(supclstr.energy() < cfPrm("minSCE"))  continue;
		fsupclstrs.push_back(supclstr);
		fscExclude.push_back(false);
		fscIsOOT.push_back(true);
	}//<<>>for( const auto &supclstr : *ootSuperCluster_ )
    for( const auto &supclstr : *superCluster_ ){
		if(supclstr.energy() < cfPrm("minSCE")) continue;
        fsupclstrs.push_back(supclstr);
        fscIsOOT.push_back(false);
        double minDr(10.0);
        double dRmatch(10.0);
        //float matchpt(0);
        auto pEta = supclstr.eta();
        auto pPhi = supclstr.phi();
        for( const auto &ootsc : *ootSuperCluster_ ){
            //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
            auto oEta = ootsc.eta();
            auto oPhi = ootsc.phi();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
            if( dRmatch < minDr ){ minDr = dRmatch; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( minDr < 0.05 ) fscExclude.push_back(true);
        else fscExclude.push_back(false);
	}//<<>>for( const auto &supclstr : *superCluster_ )

/*
    for( const auto &osupclstr : *otherSuperCluster_ ){
        if(osupclstr.energy() < cfPrm("minSCE")) continue;
        //fsupclstrs.push_back(osupclstr);
        //fscIsOOT.push_back(false);
        //const auto & oscSeedDetId = osupclstr.seed()->seed(); // get seed detid 
        //const auto oscSeedRawID = oscSeedDetId.rawId();
		//bool matched = false;
        double minDr(10.0);
        double dRmatch(10.0);
        //float matchpt(0);
        auto pEta = osupclstr.eta();
        auto pPhi = osupclstr.phi();
        for( const auto &fsupclstr : fsupclstrs ){
            //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
            auto oEta = fsupclstr.eta();
            auto oPhi = fsupclstr.phi();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
            if( dRmatch < minDr ){ minDr = dRmatch; }
            //const auto & fscSeedDetId = fsupclstr.seed()->seed(); // get seed detid 
			//const auto fscSeedRawID = fscSeedDetId.rawId();
			//if( oscSeedRawID == fscSeedRawID ){ matched = true; break; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( minDr < 0.03 ){ fscExclude.push_back(true); } //std::cout << " -- OtherSC is EX cluded !!!!! " << std::endl; }
		//if( matched ){ fscExclude.push_back(true); } //std::cout << " -- OtherSC is EX cluded !!!!! " << std::endl; }
        else { fscExclude.push_back(false); } //std::cout << " -- OtherSC is IN cluded !!!!! " << std::endl; }
    }//<<>>for( const auto &supclstr : *superCluster_ )
	//std::cout << " -- OtherSC is IN cluded : " << nOtherSC << std::endl;
*/

}//<<>>void KUCMSEcalRecHit::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )


void KUCMSEcalRecHitObject::PostProcessEvent( ItemManager<float>& geVar ){

    if( ERHODEBUG ) std::cout << "Processing RecHits" << std::endl;

    Branches.clearBranches();

    if( ERHODEBUG ) std::cout << " - enetering SuperCluster loop" << std::endl;
    int nSupClstrs = fsupclstrs.size();
    Branches.fillBranch("zscnSC",nSupClstrs);
    for ( int it = 0; it < nSupClstrs; it++ ){

        auto supclstr = fsupclstrs[it];
        bool excluded = fscExclude[it];
        bool isOOT = fscIsOOT[it];
        DetId xsDetId = supclstr.seed()->seed();
        uInt xseed = xsDetId.rawId();
        float clcx = supclstr.x();
        float clcy = supclstr.y();
        float clcz = supclstr.z();
        float sceta = supclstr.eta();
        float scphi = supclstr.phi();
        float scenergy = supclstr.rawEnergy();
        //uInt bcsize = supclstr.clustersSize();

        Branches.fillBranch("zscOOT",isOOT);
        Branches.fillBranch("zscExcluded",excluded);
        Branches.fillBranch("zscsid",xseed);
        Branches.fillBranch("zsclcx",clcx);
        Branches.fillBranch("zsclcy",clcy);
        Branches.fillBranch("zsclcz",clcz);
        Branches.fillBranch("zsceta",sceta);
        Branches.fillBranch("zscphi",scphi);
        Branches.fillBranch("zscenergy",scenergy);

        const scGroup phoSCGroup{supclstr};
        const auto phoRhGroup = getRHGroup( phoSCGroup, cfPrm("minRHEi") );
        const auto phoRhIdsGroup = getRhGrpIDs( phoRhGroup );
        Branches.fillBranch("zscrhids",phoRhIdsGroup);
        setRecHitUsed(phoRhIdsGroup);
        uInt scsize = phoRhIdsGroup.size();
        Branches.fillBranch("zscsize",scsize);

        const float phoEnergyRaw = supclstr.rawEnergy();
        const bool isScEtaEB = abs(supclstr.eta()) < 1.4442;
        const bool isScEtaEE = abs(supclstr.eta()) > 1.566 && abs(supclstr.eta()) < 2.5;
        const int seediEtaOriX = supclstr.seedCrysIEtaOrIx();
        const int seediPhiOriY = supclstr.seedCrysIPhiOrIy();
        const float x_calo = supclstr.seed()->position().x();
        const float y_calo = supclstr.seed()->position().y();
        const float z_calo = supclstr.seed()->position().z();
        //const float esEnergyOverRawE = supclstr.preshowerEnergy()/phoEnergyRaw;
        const float etaWidth = supclstr.etaWidth();
        const float phiWidth = supclstr.phiWidth();

        if( ERHODEBUG ) std::cout << " --- Finding Moments & Covariences : " << supclstr << std::endl;
        const auto &seedDetId = supclstr.seed()->seed();// seed detid
        const auto isEB = (seedDetId.subdetId() == EcalBarrel);// which subdet
        const auto ph2ndMoments = getCluster2ndMoments( &supclstr );
        const auto lCov = getCovariances( &supclstr );

        if( ERHODEBUG ) std::cout << " --- Storing Moments & Covariences : " << supclstr << std::endl;
        const float phoSMaj = ph2ndMoments.sMaj;
        const float phoSMin = ph2ndMoments.sMin;
        const float phoSAlp = ph2ndMoments.alpha;
        const float phoCovEtaEta = lCov[0];
        const float phoCovEtaPhi = lCov[1];
        const float phoCovPhiPhi = lCov[2];

        Branches.fillBranch("zscIsEB",isEB);
        Branches.fillBranch("zscEnergyRaw",phoEnergyRaw);

        Branches.fillBranch("zscIsScEtaEB",isScEtaEB);
        Branches.fillBranch("zscIsScEtaEE",isScEtaEE);
        Branches.fillBranch("zscSeediEtaOriX",seediEtaOriX);
        Branches.fillBranch("zscSeediPhiOriY",seediPhiOriY);
        Branches.fillBranch("zscx_calo",x_calo);
        Branches.fillBranch("zscy_calo",y_calo);
        Branches.fillBranch("zscz_calo",z_calo);

        Branches.fillBranch("zscEtaWidth",etaWidth);
        Branches.fillBranch("zscPhiWidth",phiWidth);
        Branches.fillBranch("zscSMaj",phoSMaj);
        Branches.fillBranch("zscSMin",phoSMin);
        Branches.fillBranch("zscSAlp",phoSAlp);
        Branches.fillBranch("zscCovEtaEta",phoCovEtaEta);
        Branches.fillBranch("zscCovEtaPhi",phoCovEtaPhi);
        Branches.fillBranch("zscCovPhiPhi",phoCovPhiPhi);

    }//<<>>for ( int it = 0; it < nSupClstrs; it++ )

    if( ERHODEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    int nRecHits = frechits.size();
    for ( int it = 0; it < nRecHits; it++ ){

        auto recHit = frechits[it];
        auto used = frhused[it];
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

        //if( DEBUG ) std::cout << " -- proccesing LASER info" << std::endl;
        // adcToGeVInfo : http://cmslxr.fnal.gov/source/RecoEcal/EgammaCoreTools/src/EcalClusterLazyTools.cc#0204
        const auto laser = laser_->getLaserCorrection(recHitID,evTime);
        const auto interCalibIter = interCalibMap->find(recHitID);
        const auto interCalib = ((interCalibIter != interCalibMap->end()) ? (*interCalibIter) : - 1.f);
        //if( DEBUG ) std::cout << " -- proccesing ADC info" << std::endl;
        //if ((laser > 0.f) && (interCalib > 0.f) && (adcToGeV > 0.f)) rhadcToGeV[pos] = (laser*interCalib*adcToGeV);
        const float adcToGeV0( isEB ? adcToGeVEB : adcToGeVEE );
        //if( DEBUG ) std::cout << " -- proccesing PED info" << std::endl;
        // pedestal info
        const auto & pediter = pedestals_->find(recHitID);
        const auto pedrms12 = (pediter != pedestals_->end()) ? pediter->rms(1) : 0.0;
        const float adcToGeV = laser*interCalib*adcToGeV0;
		const float amplitude = ( adcToGeV != 0 ) ? recHit.energy()/adcToGeV : 0;
        const float ampres = ( pedrms12 != 0 ) ? amplitude/pedrms12 : 0;


        const float rhTime = recHit.time();
        const bool rhIsOOT = recHit.checkFlag(EcalRecHit::kOutOfTime);
        const float rhEnergy = recHit.energy();
        //const float rhAdjTime = rhTime-d_rh;  // Note : Margret adds d_rh to the time in her code & subtracts d_pv ( or TOF )

        Branches.fillBranch("ID",recHitID);
        Branches.fillBranch("TOFpv",d_pv);
        Branches.fillBranch("TOF0",d_rh);
        Branches.fillBranch("Time",rhTime);
        Branches.fillBranch("isOOT",rhIsOOT);
        Branches.fillBranch("Energy",rhEnergy);
        Branches.fillBranch("SwCross",swisscross);
        Branches.fillBranch("eta",eta);
        Branches.fillBranch("phi",phi);
        Branches.fillBranch("rhx",rhX);
        Branches.fillBranch("rhy",rhY);
        Branches.fillBranch("rhz",rhZ);
		Branches.fillBranch("amplitude",amplitude);
		Branches.fillBranch("ampres",ampres);
        //Branches.fillBranch("rhisWeird",recHit.checkFlag(EcalRecHit::kWeird));
        //Branches.fillBranch("rhisDiWeird",recHit.checkFlag(EcalRecHit::kDiWeird));

    }//<<>>for (const auto recHit : *recHitsEB_ ) 


    nOtherSC = otherSuperCluster_->size();
    nOtherSCEx = 0;
	int nOtherSCIn = 0;
    for( const auto &osupclstr : *otherSuperCluster_ ){
	//for( const auto &osupclstr : fsupclstrs ){
        if(osupclstr.energy() < cfPrm("minSCE")) continue;
        const auto & oscSeedDetId = osupclstr.seed()->seed(); // get seed detid 
        const auto oscSeedRawID = oscSeedDetId.rawId();
		Branches.fillBranch("zscnOtherSID",oscSeedRawID);
        const scGroup oscgroup{osupclstr};
        const auto orhgroup = getRHGroup( oscgroup, 0 );
        const auto orhidsgroup = getRhGrpIDs( orhgroup );
        //bool matched = false;
        float minDr(10.0);
        float dRmatch(10.0);
        auto pEta = osupclstr.eta();
        auto pPhi = osupclstr.phi();
		int nOverlap(-10);
		uInt oMatchSID(0);
		int iter(-1);
        for( const auto &fsupclstr : fsupclstrs ){
			iter++;
			if( fscIsOOT[iter] == false ) continue;
            auto oEta = fsupclstr.eta();
            auto oPhi = fsupclstr.phi();
            dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, oEta, oPhi ));
			const scGroup fscgroup{fsupclstr};
        	const auto frhgroup = getRHGroup( fscgroup, 0 );
        	const auto frhidsgroup = getRhGrpIDs( frhgroup );
            const auto & fscSeedDetId = fsupclstr.seed()->seed(); // get seed detid 
            const auto fscSeedRawID = fscSeedDetId.rawId();
            if( dRmatch < minDr ){ 
				minDr = dRmatch; 
				nOverlap = getOverLapCnt(orhidsgroup,frhidsgroup); 
				oMatchSID = fscSeedRawID;
				//std::cout << " -- OtherSC new min : " << minDr << " With overlap: " << nOverlap;
				//std::cout << " from: " << orhidsgroup.size() << ", " << frhidsgroup.size() << std::endl;
			}//<<>>if( dRmatch < minDr )
        	//if( oscSeedRawID == fscSeedRawID ){ matched = true; break; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
		Branches.fillBranch("zscnOExDr",minDr);
        Branches.fillBranch("zscnOver",nOverlap);
		Branches.fillBranch("zscnOMatchSID",oMatchSID);
        if( minDr < 0.03 ){ nOtherSCIn++; } //s"
        //if( matched ){ nOtherSCIn++; } //s"
        //{ fscExclude.push_back(true); } //std::cout << " -- OtherSC is EX cluded !!!!! " << std::endl; }
        else { nOtherSCEx++; } //Branches.fillBranch("zscnOExDr",dRmatch); } //std::cout << " -- OtherSC is IN cluded !!!!! " << std::endl; }
    }//<<>>for( const auto &supclstr : *superCluster_ )
    //std::cout << " -- OtherSC is IN cluded : " << nOtherSC << std::endl;
    Branches.fillBranch("zscnOther",nOtherSC);
    Branches.fillBranch("zscnOtherEx",nOtherSCEx);
    Branches.fillBranch("zscnOtherIn",nOtherSCIn);

}//<<>>void KUCMSEcalRecHit::ProcessEvent()

void KUCMSEcalRecHitObject::EndJobs(){}

void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup){

    int nRecHits = frechits.size();
    for( auto rhid : idgroup ){
        for( int it = 0; it < nRecHits; it++ ){ if( getRawID( frechits[it] ) == rhid ){ frhused[it] = true; break; } }
    }//<<>>for( auto rhid : phoRhIdsGroup )

};//<<>>void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup)

float KUCMSEcalRecHitObject::getSuperClusterSeedTime( reco::SuperClusterRef supercluster ){

    const auto & seedDetId = supercluster->seed()->seed(); // get seed detid
    const auto recHits = ((seedDetId.subdetId() == EcalSubdetector::EcalBarrel) ? recHitsEB_ : recHitsEE_); // which recHits to use
    const auto seedHit = recHits->find(seedDetId); // get the underlying rechit
    const auto seedTime = ((seedHit != recHits->end()) ? seedHit->time() : -9999.f);
    //if( DEBUG && seedTime == -9999.f ) std::cout << "Bad Photon seed time !!!! " << std::endl;
    return seedTime;

}//<<>>getSuperClusterSeedTime( reco::SuperClusterRef supercluster )

int KUCMSEcalRecHitObject::getSuperClusterIndex( const reco::SuperCluster* supercluster ){

	const auto & scSeedDetId = supercluster->seed()->seed(); // get seed detid	
	const auto scSeedRawID = scSeedDetId.rawId();
	//const auto osceta = supercluster->eta();
    //const auto oscphi = supercluster->phi();
	int fscIter = 0;
	int retIndex = -1;
	//float mindr = 0.4;
    //std::cout << " -- SC matching : " << scSeedRawID << " eta: " << osceta << std::endl;
	for( const auto & fsclsrt : fsupclstrs ){ 
		
		const auto & fscSeedDetId = fsclsrt.seed()->seed(); // get seed detid 
		const auto fscSeedRawID = fscSeedDetId.rawId();
		//std::cout << " --- SC compared to : " << scSeedRawID << " ?= " << fscSeedRawID << std::endl;
		if( scSeedRawID == fscSeedRawID ){ retIndex = fscIter; break; }
		//const auto sceta = fsclsrt.eta();
		//const auto scphi = fsclsrt.phi();
		//auto dr = std::sqrt(reco::deltaR2( osceta, oscphi, sceta, scphi ));
		//if( dr < mindr ){ mindr = dr; retIndex = fscIter; std::cout << " --- SC matched: " << retIndex << " dr: " << mindr << std::endl; }
		fscIter++;		

	}//<<>>for( sclsrt : fsupclstrs )		
	//std::cout << " --- SC matched: " << retIndex << std::endl;
	return retIndex;

}//<<>>float KUCMSEcalRecHitObject::getSuperClusterIndex( reco::SuperCluster supercluster )

bool KUCMSEcalRecHitObject::reduceRhGrps( std::vector<rhGroup> & x ){

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

}//>>>>std::vector<rhGroup> KUCMSEcalRecHitObject::reduceRhGrps(const std::vector<rhGroup> x)

void KUCMSEcalRecHitObject::mrgRhGrp( rhGroup & x, rhGroup & y ){

    //std::cout << " --- In mrgRhGrp " << std::endl;
    bool matched(false);
    for(const auto rhy : y ){
        for(const auto rhx : x ){ if( rhMatch(rhx,rhy) ) matched = true; }
        if( not matched ) x.push_back(rhy);
    }//<<>>for(const auto rhx : y )
   //std::cout << " ---- Finished mrgRhGrp " << std::endl;

    return;

}//>>>>rhGroup KUCMSEcalRecHitObject::mrgRhGrp(const rhGroup x, const rhGroup y)

float KUCMSEcalRecHitObject::getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    if( recHits.size() == 0 ) return -99.0;
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

}//>>>>float  KUCMSEcalRecHitObject::getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

float KUCMSEcalRecHitObject::getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ ){

    float rhTime(-100.0);
    const auto seedDetId = sprclstr.seed().get()->seed(); // seed detid
    const auto seedRawId = seedDetId.rawId(); // crystal number
    const auto isEB = (seedDetId.subdetId() == EcalSubdetector::EcalBarrel); // which subdet
    const auto recHits = (isEB ? *recHitsEB_ : *recHitsEE_ );
    for( const auto &rechit : recHits ){
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

}//>>>>float  KUCMSEcalRecHitObject::getSeedTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

std::vector<float> KUCMSEcalRecHitObject::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    std::vector<float> result;
    if( recHits.size() == 0 ){ result.push_back(-99.0); return result; }
    auto lrh = getLeadRh(recHits);
    const auto recHitId(lrh.detid());
    const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
    const auto rhPosX = recHitPos.x();
    const auto rhPosY = recHitPos.y();
    const auto rhPosZ = recHitPos.z();
    const auto d_rh = hypo(rhPosX,rhPosY,rhPosZ);
    const auto d_pv = hypo(rhPosX-vtxX,rhPosY-vtxY,rhPosZ-vtxZ);
    const auto tof = (d_rh-d_pv)/SOL;
    for (const auto &recHit : recHits ){result.push_back(recHit.time()-tof);}
    return result;

}//>>>>std::vector<float>  KUCMSEcalRecHitObject::getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

EcalRecHit KUCMSEcalRecHitObject::getLeadRh( rhGroup recHits ){

    EcalRecHit result;
    float enr(0.0);
    for (const auto &recHit : recHits ){
        auto rhenr = recHit.energy();
        if( rhenr < enr ) continue;
        enr = rhenr;
        result = recHit;
    }//<<>>for (const auto recHit : recHits )

    return result;

}//>>>>EcalRecHit KUCMSEcalRecHitObject::getLeadRh( rhGroup recHits )

float KUCMSEcalRecHitObject::getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ ){

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

std::vector<float> KUCMSEcalRecHitObject::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ ){

    std::vector<float> result;
    if( recHits.size() == 0 ){ result.push_back(-99.0); return result; }
    for (const auto &recHit : recHits ){
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

}//>>>>std::vector<float>  KUCMSEcalRecHitObject::getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ )

rhGroup KUCMSEcalRecHitObject::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 ){

    rhGroup result;
    std::vector<uInt> rawIds;
    auto & hitsAndFractions = basicCluster.hitsAndFractions();
    const auto nHAF = hitsAndFractions.size();
    //std::cout << " --- getRHGroup for basic cluster with " << nHAF << " rechits " << std::endl;
    for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
        const auto detId = hitsAndFractions[iHAF].first;
        const auto rawId = detId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
    }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    for (const auto &recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto &recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )       

    return result;

}////rhGroup KUCMSEcalRecHitObject::getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.0 )

rhGroup KUCMSEcalRecHitObject::getRHGroup( const scGroup superClstrGrp, float minenr, std::vector<float> phEnergy, std::vector<float> phDr, float phEnMax ){

    rhGroup result;
    std::vector<uInt> rawIds;

    float dRcut = 0.05;

    int iter = -1;
    for ( const auto &superCluster : superClstrGrp ){
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
    for (const auto &recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto &recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup( const scGroup superClstrGrp, float minenr, std::vector<float> phEnergy, std::vector<float> phDr, float phEnMax )

rhGroup KUCMSEcalRecHitObject::getRHGroup( const scGroup superClusterGroup, float minenr ){

    rhGroup result;
    std::vector<uInt> rawIds;
    for ( const auto &superCluster : superClusterGroup ){
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
            if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ) rawIds.push_back(rawId);
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )  
    for (const auto &recHit : *recHitsEB_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto &recHit : *recHitsEE_ ){
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( std::find( rawIds.begin(), rawIds.end(), rawId ) != rawIds.end() ) result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup( const scGroup superClusterGroup, float minenr = 0.0 )

rhGroup KUCMSEcalRecHitObject::getRHGroup( uInt detid ){

    rhGroup result;
    for (const auto &recHit : *recHitsEB_ ){
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( detid != rawId ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto &recHit : *recHitsEE_ ){
        const auto recHitId = recHit.detid();
        const auto rawId = recHitId.rawId();
        if( detid != rawId ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEE_ )

    return result;
}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup( uInt detid )

rhGroup KUCMSEcalRecHitObject::getRHGroup(){

    rhGroup result;
    for (const auto &recHit : *recHitsEB_ ) result.push_back(recHit);
    for (const auto &recHit : *recHitsEE_ ) result.push_back(recHit);

    return result;

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup()

rhGroup KUCMSEcalRecHitObject::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 ){

    rhGroup result;
    for (const auto &recHit : *recHitsEB_ ){
        if( recHit.checkFlag(EcalRecHit::kHasSwitchToGain6) ) continue;
        auto enr = recHit.energy();
        if( enr <= minenr ) continue;
        const auto recHitId(recHit.detid());
        const auto recHitPos = barrelGeometry->getGeometry(recHitId)->getPosition();
        const auto dr = std::sqrt(reco::deltaR2(eta, phi, recHitPos.eta(), recHitPos.phi()));
        if( dr > drmin ) continue;
        result.push_back(recHit);
    }//<<>>for (const auto recHit : *recHitsEB_ )
    for (const auto &recHit : *recHitsEE_ ){
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

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup( float eta, float phi, float drmin, float minenr = 0.0 )

Cluster2ndMoments KUCMSEcalRecHitObject::getCluster2ndMoments( const reco::SuperCluster* scptr ){

    const auto &seedDetId = scptr->seed()->seed();// seed detid
	const auto isEB = (seedDetId.subdetId() == EcalBarrel);// which subdet
    const auto recHits = ( isEB ? recHitsEB_.product() : recHitsEE_.product() );
	return noZS::EcalClusterTools::cluster2ndMoments( *(scptr->seed()), *recHits );

}//<<>>Cluster2ndMoments KUCMSEcalRecHitObject::getCluster2ndMoments( reco::SuperCluster* scptr )

std::array<float, 3> KUCMSEcalRecHitObject::getCovariances( const reco::SuperCluster* scptr ){

    const auto &seedDetId = scptr->seed()->seed();// seed detid
    const auto isEB = (seedDetId.subdetId() == EcalBarrel);// which subdet
    const auto recHits = ( isEB ? recHitsEB_.product() : recHitsEE_.product() );
    return EcalClusterTools::covariances( *(scptr->seed()), recHits, ecalTopology, ecalGeometry );

}//<<>>std::array<float, 3> KUCMSEcalRecHitObject::getCovariances( reco::SuperCluster* scptr )

#endif
