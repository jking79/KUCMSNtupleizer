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
    void LoadBeamSpotTokens( edm::EDGetTokenT<reco::BeamSpot> beamLineToken ){ beamLineToken_ = beamLineToken; };
    // sets up branches, do preloop jobs 
    void InitObject( TTree* fOutTree ); 
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject* otherObject ){ otherObjectPtr = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    // get collections, do initial processing
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup, ItemManager<float>& geVar ); 
    // do cross talk jobs with other objects, do event processing, and load branches
    void ProcessEvent( ItemManager<float>& geVar ); 
	void PostProcessEvent( ItemManager<float>& geVar );

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    // sc functions
    float getSuperClusterSeedTime( reco::SuperClusterRef supercluster );
	int getSuperClusterIndex( const reco::SuperCluster* supercluster, int objectPdgId, int objectIdx );
    std::vector<int> getSuperClusterIndex( const rhIdGroup crystals, int objectIdx );

    // rechit group functions
    rhGroup getRHGroup( float eta, float phi, float drmin, float minenr = 0.f );
    rhGroup getRHGroup( const scGroup superClusterGroup, float minenr = 0.f );
    rhGroup getRHGroup( const scGroup superClusterGroup, float minenr, std::vector<float> phEnergy, std::vector<float> phDr, float phEnMax );
    rhGroup getRHGroup( const reco::CaloCluster basicCluster, float minenr = 0.f );
    rhGroup getRHGroup( uInt detid );
    rhGroup getRHGroup();


    EcalRecHit getLeadRh( rhGroup recHits );
    std::vector<float> getRhTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
    std::vector<float> getLeadTofRhTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
	rhIdGroup getRhRawIdList( const scGroup superClusterGroup );
    std::vector<float> getFractionList( const scGroup superClusterGroup, rhGroup recHits );
    std::vector<float> getMissFractionList( const scGroup superClusterGroup );
	std::tuple<uInt,float> getHFList( const scGroup superClusterGroup );
	float getRecHitEnergy( const uInt id );
    float getRecHitTime( const uInt id );

    float getRhTOF( EcalRecHit rechit, double vtxX, double vtxY, double vtxZ );
    float getLeadTofTime( rhGroup recHits, double vtxX, double vtxY, double vtxZ );
    float getSeedTofTime( reco::SuperCluster sprclstr, double vtxX, double vtxY, double vtxZ );
    void  mrgRhGrp( rhGroup & x, rhGroup & y );
    bool  reduceRhGrps( std::vector<rhGroup> & x );
    void  setRecHitUsed( rhIdGroup idgroup );
    void  setRecHitUsed( rhIdGroup idgroup, bool isOOT, float energy, std::vector<float> fracs );

	// Helper functions for rechit infomation & IDs
    const uInt getRawID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.rawId();}
    const bool getIsEB(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return (recHitId.subdetId() == EcalBarrel)?1:0;}
    const auto getSubDetID(const EcalRecHit recHit){ auto recHitId = recHit.detid(); return recHitId.subdetId();}
    const bool rhMatch(const EcalRecHit rhx, const EcalRecHit rhy){ return getRawID(rhx) == getRawID(rhy);}
    const bool dupRhFnd(const rhGroup x, const rhGroup y){
               for(auto rhx : x ){ for(auto rhy : y ){if(rhMatch(rhx,rhy)){ return true;}}} return false;}
    const bool isRhGrpEx(const rhGroup x){ int s=x.size();for( int i=0;i<s;i++){
               for( int j=i+1;j<s;j++){if(rhMatch(x[i],x[j])) return false;}} return true;}
    const float getRhGrpEnr(const rhGroup x){float e(0.0);for( auto ix : x ){e+=ix.energy();} return e;}
    const int getDupCnt(const std::vector<rhGroup> x){ int c=0; int s=x.size();
               for( int a=0;a<s;a++){for( int b=a+1;b<s;b++){if(dupRhFnd(x[a],x[b]))c++;}} return c;}
    const rhIdGroup getRhGrpIDs(const rhGroup rhs ){ rhIdGroup rt; if(rhs.empty()){ rt.push_back(0);} else{
               for(const auto rh : rhs ){ rt.push_back(getRawID(rh));}} return rt;}
    const int getOverLapCnt(const rhIdGroup x, const rhIdGroup y){ int c=0; int sx=x.size(); int sy=y.size();
               for( int a=0;a<sx;a++){for( int b=0;b<sy;b++){if(x[a]==y[b]){c++;break;}}} return c;}

	// Wrapper functions for EcalClusterTools functions 
	Cluster2ndMoments getCluster2ndMoments( const reco::SuperCluster* scptr );
    std::array<float, 3> getCovariances( const reco::SuperCluster* scptr );

    private:

    std::vector<EcalRecHit> frechits;
    std::vector<bool> frhused;

    std::vector<reco::SuperCluster> fsupclstrs;
	std::vector<rhIdGroup> fscRhIdGrpList;
    std::vector<bool> fscExclude;
    std::vector<bool> fscIsOOT;
    std::vector<int> fscOType;
    std::vector<int> fscPhoIndx;
    std::vector<int> fscEleIndx;
    std::vector<std::vector<int>> fscTrackIndx;
    std::vector<bool> fscOriginal;
	std::map<uInt,int> scRechitMap;

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

    // BeamSpot ----------------------------------------------------------
    //const edm::InputTag beamSpotTag;
    edm::EDGetTokenT<reco::BeamSpot> beamLineToken_;
    edm::Handle<reco::BeamSpot> beamSpot_;

    // geometry CaloSubdetectorGeometry
    edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
    edm::ESHandle<CaloGeometry> caloGeo_;
    const CaloSubdetectorGeometry * barrelGeometry;
    const CaloSubdetectorGeometry * endcapGeometry;
    //const CaloGeometry * ecalGeometry;

    // CaloTopology
    edm::ESGetToken<CaloTopology, CaloTopologyRecord> caloTopologyToken_;
    edm::ESHandle<CaloTopology> caloTopo_;
    //const CaloSubdetectorTopology * barrelTopology;
    //const CaloSubdetectorTopology * endcapTopology;
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
	cfFlag.set( "doProbeOut", iConfig.existsAs<bool>("doProbeOut") ? iConfig.getParameter<bool>("doProbeOut") : false );

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
    Branches.makeBranch("TimeError","ECALRecHit_timeError",VFLOAT);
    Branches.makeBranch("isTimeValid","ECALRecHit_isTimeValid",VBOOL);
    Branches.makeBranch("TOFpv","ECALRecHit_pvTOF",VFLOAT);
    Branches.makeBranch("TOF0","ECALRecHit_0TOF",VFLOAT);
    Branches.makeBranch("ID","ECALRecHit_ID",VUINT);
    Branches.makeBranch("isOOT","ECALRecHit_isOOT",VBOOL);
    Branches.makeBranch("isWrd","ECALRecHit_isWrd",VBOOL);
    Branches.makeBranch("isDiWrd","ECALRecHit_isDiWrd",VBOOL);
    Branches.makeBranch("isRecov","ECALRecHit_isRecov",VBOOL);
    Branches.makeBranch("isPoor","ECALRecHit_isPoor",VBOOL);
    Branches.makeBranch("isDead","ECALRecHit_isDead",VBOOL);
    Branches.makeBranch("isOther","ECALRecHit_isOther",VBOOL);
    Branches.makeBranch("hasGS1","ECALRecHit_hasGS1",VBOOL);
    Branches.makeBranch("hasGS6","ECALRecHit_hasGS6",VBOOL);
    Branches.makeBranch("SwCross","ECALRecHit_swCross",VFLOAT);
    Branches.makeBranch("eta","ECALRecHit_eta",VFLOAT);
    Branches.makeBranch("phi","ECALRecHit_phi",VFLOAT);
    Branches.makeBranch("rhx","ECALRecHit_rhx",VFLOAT);
    Branches.makeBranch("rhy","ECALRecHit_rhy",VFLOAT);
    Branches.makeBranch("rhz","ECALRecHit_rhz",VFLOAT);
    Branches.makeBranch("amplitude","ECALRecHit_amplitude",VFLOAT);
    Branches.makeBranch("ampres","ECALRecHit_ampres",VFLOAT);
    Branches.makeBranch("pused","ECALRecHit_precentUsed",VFLOAT);

    Branches.makeBranch("zscnSC","SuperCluster_nSuperCluster",INT);
    Branches.makeBranch("zscOOT","SuperCluster_isOot",VBOOL,"This SC is from the out of time SC collection");
    Branches.makeBranch("zscExcluded","SuperCluster_excluded",VBOOL,"This SC is superseeded by an overlapping OOT SC");
    Branches.makeBranch("zscOriginal","SuperCluster_original",VBOOL,"This SC is only found in the orignal propmt or oot collection");
    Branches.makeBranch("zscsid","SuperCluster_XtalSeedID",VUINT);
    Branches.makeBranch("zscstime","SuperCluster_SeedTime",VFLOAT);
    Branches.makeBranch("zscotype","SuperCluster_ObjectPdgId",VINT,"Currently is 0( no match ), 11, 22, or 33( matches an ele & pho )");
    Branches.makeBranch("zscphoindx","SuperCluster_PhotonIndx",VINT,"Index of matching photon");
    Branches.makeBranch("zsceleindx","SuperCluster_ElectronIndx",VINT,"index of matching electron");
    Branches.makeBranch("zsctrackindx","SuperCluster_TrackIndx",VVINT,"index of intersecting track");
    Branches.makeBranch("zsclcx","SuperCluster_clcx",VFLOAT,"x coordinate of cluster centroid");
    Branches.makeBranch("zsclcy","SuperCluster_clcy",VFLOAT,"y coordinate of cluster centroid");
    Branches.makeBranch("zsclcz","SuperCluster_clcz",VFLOAT,"z coordinate of cluster centroid");
    Branches.makeBranch("zsceta","SuperCluster_eta",VFLOAT,"pseudorapidity of cluster centroid");
    Branches.makeBranch("zscphi","SuperCluster_phi",VFLOAT,"azimuthal angle of cluster centroid");
    Branches.makeBranch("zscenergy","SuperCluster_energy",VFLOAT);
    Branches.makeBranch("zscnrh","SuperCluster_nRHXtals",VUINT, "size of the rechitgroup found from the H&F list for the supercluster");
    Branches.makeBranch("zscsize","SuperCluster_nHFXtals",VUINT,"size of the H&F list for the supercluster");
    //Branches.makeBranch("zscsizedif","SuperCluster_diffXtrals",VUINT,"diff of nRHXtals and nHFXtals");
    Branches.makeBranch("zscrhids","SuperCluster_rhIds",VVUINT,"list of rechit raw ids in hits and fractions list from supercluster");
    Branches.makeBranch("zscrhfracs","SuperCluster_rhFracs",VVFLOAT,"list of rechit energy fractions in hits and fractions list from supercluster");
    //Branches.makeBranch("zscmissrhfracs","SuperCluster_MissingRhFracs",VVFLOAT,"list of rechit fractions in H&F list from supercluster not in rechits");
    Branches.makeBranch("zscbcsize","SuperCluster_nBasicClusters",VUINT,"size of the BasicCluuster list for the supercluster");

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

    //Branches.makeBranch("zscnOther","OSuperCluster_nOther",INT);
    //Branches.makeBranch("zscnOtherEx","OSuperCluster_nOtherEx",INT);
    //Branches.makeBranch("zscnOtherIn","OSuperCluster_nOtherIn",INT);

    //Branches.makeBranch("zscnOExDr","OSuperCluster_dR",VFLOAT);
    //Branches.makeBranch("zscnOExDiffE","OSuperCluster_diffEnergy",VFLOAT);
    //Branches.makeBranch("zscnOExPerDiffE","OSuperCluster_perDiffEnergy",VFLOAT);
    //Branches.makeBranch("zscnOver","OSuperCluster_nXtalOverlap",VINT);
    //Branches.makeBranch("zscnOtherSID","OSuperCluster_otherSeedID",VUINT);
    //Branches.makeBranch("zscnOMatchSID","OSuperCluster_otherMatchSeedID",VUINT);
    //Branches.makeBranch("zscpOver","OSuperCluster_precentXtalOverlap",VFLOAT);

    //Branches.makeBranch("zscnOExDrX","OSuperCluster_XdR",VFLOAT);
    //Branches.makeBranch("zscnOExDiffEX","OSuperCluster_XdiffEnergy",VFLOAT);
    //Branches.makeBranch("zscnOExPerDiffEX","OSuperCluster_XPerDiffEnergy",VFLOAT);
    //Branches.makeBranch("zscnOverX","OSuperCluster_XnXtalOverlap",VINT);
    //Branches.makeBranch("zscnOMatchSIDX","OSuperCluster_XotherMatchSeedID",VUINT);
    //Branches.makeBranch("zscpOverX","OSuperCluster_XprecentXtalOverlap",VFLOAT);


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

	// BeamSpot
    iEvent.getByToken(beamLineToken_,beamSpot_);

    // GEOMETRY : https://gitlab.cern.ch/shervin/ECALELF
    caloGeo_ = iSetup.getHandle(caloGeometryToken_);
    barrelGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalBarrel);
    endcapGeometry = caloGeo_->getSubdetectorGeometry(DetId::Ecal, EcalSubdetector::EcalEndcap);
    //ecalGeometry = &iSetup.getData(caloGeometryToken_);

    // CaloTopology
    caloTopo_ = iSetup.getHandle(caloTopologyToken_);
    ecalTopology = &iSetup.getData(caloTopologyToken_);
    //barrelTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalBarrel);
    //endcapTopology = caloTopo_->getSubdetectorTopology(DetId::Ecal, EcalSubdetector::EcalEndcap);

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
	//std::cout << "Collecting ECAL RecHits w " << cfPrm("minRHEi") << std::endl;
    frechits.clear();
    frhused.clear();
    //for( const auto &recHit : *recHitsEB_ ){ if(recHit.energy() > cfPrm("minRHEi")){ frechits.push_back(recHit); frhused.push_back(false);}}
    //for( const auto &recHit : *recHitsEE_ ){ if(recHit.energy() > cfPrm("minRHEi")){ frechits.push_back(recHit); frhused.push_back(false);}}
    for( const auto &recHit : *recHitsEB_ ){ frechits.push_back(recHit); frhused.push_back(false);}
    for( const auto &recHit : *recHitsEE_ ){ frechits.push_back(recHit); frhused.push_back(false);}
/*
	std::string ConfigEB("/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/macros/ecal_config/fullinfo_detids_EE.txt");
	std::string newConfigEB("/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/macros/ecal_config/fullinfo_v2_detids_EE.txt");
	std::ifstream infile( ConfigEB, std::ios::in);
    std::ofstream outfile( newConfigEB, std::ios::out);	

    //UInt_t cmsswId, dbID;
    //Int_t hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    //std::string pos;
    UInt_t cmsswId, dbID;
    Int_t hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    std::string EE;
	//for( uInt i = 838904321; i <= 838970216; i++ ){ 
	//while (infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM){
    while (infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant){

		DetId thedetid(cmsswId);
		//auto geometry( isEB ? barrelGeometry : endcapGeometry );
        auto recHitPos = endcapGeometry->getGeometry(thedetid)->getPosition();
        const float eta = recHitPos.eta();
        const float phi = recHitPos.phi();
		outfile << cmsswId << " " << dbID << " " << hashedId << " " << side << " " << ix << " " << iy << " " << SC;
        outfile << " " << iSC << " " << Fed << " " << EE << " " << TTCCU << " " << strip << " " << Xtal << " " << quadrant;
        outfile << " " << phi << " " << eta << std::endl;

	}//<<>>for( uInt i = 838904321; i <= 838970216; i++ )
	infile.close();
	outfile.close();
*/
    if( ERHODEBUG ) std::cout << "Collecting SuperClusters" << std::endl;

    fsupclstrs.clear();
    fscExclude.clear();
    fscIsOOT.clear();
    fscOType.clear();
	fscPhoIndx.clear();
	fscEleIndx.clear();
	fscTrackIndx.clear();
    fscOriginal.clear();
	scRechitMap.clear();

    for( const auto &ootSupclstr : *ootSuperCluster_ ){
		if(ootSupclstr.energy() < cfPrm("minSCE"))  continue;
		fsupclstrs.push_back(ootSupclstr);
		fscIsOOT.push_back(true);
		fscOriginal.push_back(true);
	}//<<>>for( const auto &supclstr : *ootSuperCluster_ )

    for( const auto &supclstr : *superCluster_ ){
        if( supclstr.energy() < cfPrm("minSCE")) continue;
        bool found = false;
        auto pEta = supclstr.eta();
        auto pPhi = supclstr.phi();
        for( const auto &fsc : fsupclstrs ){
            auto fEta = fsc.eta();
            auto fPhi = fsc.phi();
            float dRmatch = std::sqrt(reco::deltaR2( pEta, pPhi, fEta, fPhi ));
			if( dRmatch < 0.1 ){ found = true; break; }
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( not found ){
            fsupclstrs.push_back(supclstr);
            fscIsOOT.push_back(false);
            fscOriginal.push_back(true);
        }//<<>>if( found == -1 )
		else {
            scGroup scSCGroup{supclstr};
            rhIdGroup scRhIdsGroup = getRhRawIdList( scSCGroup );
            setRecHitUsed( scRhIdsGroup );
		}//<<>>if( not found )
    }//<<>>for( const auto &supclstr : *superCluster_ )

    for( const auto &oSupclstr : *otherSuperCluster_ ){
        if(oSupclstr.energy() < cfPrm("minSCE")) continue;
        bool skip = false;
        const auto & pSeedDetId = oSupclstr.seed()->seed(); // get seed detid 
        const auto pSeedRawID = pSeedDetId.rawId();
        for( const auto &fsc : fsupclstrs ){
            const auto & fSeedDetId = fsc.seed()->seed(); // get seed detid 
            const auto fSeedRawID = fSeedDetId.rawId();
            if( pSeedRawID == fSeedRawID ){ skip = true;  continue; }
            //if( cfFlag("onlyEB") && ootPho.isEE() ) continue;
        }//<<>>for( int ip; ip < nPhotons; ip++ )
        if( not skip ){
            fsupclstrs.push_back(oSupclstr);
            fscIsOOT.push_back(false);
        	fscOriginal.push_back(false);
        }//<<>>if( not skip )
    }//<<>>for( const auto &supclstr : *superCluster_ )

	int scIndx = -1;
	int nSCRechits = 0;
	//int nFSCs = fsupclstrs.size();
	for( const auto &sc : fsupclstrs ){
	//for( int iter = nFSCs; iter < 0; iter-- ){	

			//auto sc = fsupclstrs[iter];
			scIndx++;
            fscOType.push_back(0);
            fscPhoIndx.push_back(-1);
            fscEleIndx.push_back(-1);
			std::vector<int> trackIndexList{-1};
            fscTrackIndx.push_back(trackIndexList);
            fscExclude.push_back(false);
            scGroup scSCGroup{sc};
            rhIdGroup scRhIdsGroup = getRhRawIdList( scSCGroup );
			setRecHitUsed( scRhIdsGroup );
			nSCRechits += scRhIdsGroup.size();
            fscRhIdGrpList.push_back(scRhIdsGroup);
			for( const auto rhid : scRhIdsGroup ){
				if( scRechitMap.find(rhid) == scRechitMap.end() ) scRechitMap[rhid] = scIndx;
				//scRechitMap[rhid] = scIndx;
			}//<<>>for( const auto rhid : scRhIdsGroup )

	}//<<>>for( const auto sc : fsupclstrs )

	//std::cout << "SC - Rechit Map ------------------------" << std::endl;
	//std::cout << " # SC  " << fsupclstrs.size() << " -- nSCRecHits : " << nSCRechits << " v " <<  scRechitMap.size() << std::endl;
	//for( auto & ent : scRechitMap ){ std::cout << " -- " << ent.first << " : " << ent.second << std::endl; }
	//std::cout << "SC - Rechit Map ------------------------" << std::endl;

}//<<>>void KUCMSEcalRecHit::LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup )

void KUCMSEcalRecHitObject::ProcessEvent( ItemManager<float>& geVar ){}

void KUCMSEcalRecHitObject::PostProcessEvent( ItemManager<float>& geVar ){

    if( ERHODEBUG ) std::cout << "Processing RecHits && SuperClusters" << std::endl;

    Branches.clearBranches();

    if( ERHODEBUG ) std::cout << " - enetering SuperCluster loop" << std::endl;

	float nNpSupClstrs( 0 );//# of non-prompt super clusters
    int nSupClstrs = fsupclstrs.size();
    Branches.fillBranch("zscnSC",nSupClstrs);
    for ( int it = 0; it < nSupClstrs; it++ ){

        auto supclstr = fsupclstrs[it];
        bool excluded = fscExclude[it];
        bool original = fscOriginal[it];
        bool isOOT = fscIsOOT[it];
        DetId xsDetId = supclstr.seed()->seed();
        uInt xseed = xsDetId.rawId();
		float time = getRecHitTime( xseed );
        float clcx = supclstr.x();
        float clcy = supclstr.y();
        float clcz = supclstr.z();
        float sceta = supclstr.eta();
        float scphi = supclstr.phi();
        float scenergy = supclstr.rawEnergy();
        uInt bcsize = supclstr.clustersSize();
		int scOType = fscOType[it];
        int scPhoIdx = fscPhoIndx[it];
        int scEleIdx = fscEleIndx[it];

        Branches.fillBranch("zscOOT",isOOT);
        Branches.fillBranch("zscExcluded",excluded);
        Branches.fillBranch("zscOriginal",original);
        Branches.fillBranch("zscsid",xseed);
        Branches.fillBranch("zscstime",time);
        Branches.fillBranch("zscotype",scOType);
        Branches.fillBranch("zscphoindx",scPhoIdx);
        Branches.fillBranch("zsceleindx",scEleIdx);
        Branches.fillBranch("zsclcx",clcx);
        Branches.fillBranch("zsclcy",clcy);
        Branches.fillBranch("zsclcz",clcz);
        Branches.fillBranch("zsceta",sceta);
        Branches.fillBranch("zscphi",scphi);
        Branches.fillBranch("zscenergy",scenergy);

        const scGroup scSCGroup{supclstr};
        const rhGroup scRhGroup = getRHGroup( scSCGroup );
        //const rhIdGroup scRhIdsGroup = getRhGrpIDs( scRhGroup );
		const rhIdGroup scRhIdsGroup = fscRhIdGrpList[it];
		std::vector<float> fracList = getFractionList( scSCGroup, scRhGroup );
        //std::vector<float> missFracList = getMissFractionList( scSCGroup );
        Branches.fillBranch("zscrhids",scRhIdsGroup);
        Branches.fillBranch("zscrhfracs",fracList);
        //Branches.fillBranch("zscmissrhfracs",missFracList);
        //setRecHitUsed(scRhIdsGroup); // now done above
        //setRecHitUsed(scRhIdsGroup, isOOT, scenergy, fracList );
        uInt sc_rhg_size = scRhIdsGroup.size();
        uInt sc_hfl_size = (supclstr.hitsAndFractions()).size();
		//uInt sc_dups_size = missFracList.size();
		//auto basicClusters = supclstr.clusters();
        Branches.fillBranch("zscsize",sc_hfl_size);
        Branches.fillBranch("zscnrh",sc_rhg_size);
        //Branches.fillBranch("zscsizedif",sc_dups_size);
        Branches.fillBranch("zscbcsize",bcsize);

 		Branches.fillBranch("zsctrackindx",fscTrackIndx[it]);

        const float scEnergyRaw = supclstr.rawEnergy();
        const bool isScEtaEB = abs(supclstr.eta()) < 1.4442;
        const bool isScEtaEE = abs(supclstr.eta()) > 1.566 && abs(supclstr.eta()) < 2.5;
        const int seediEtaOriX = supclstr.seedCrysIEtaOrIx();
        const int seediPhiOriY = supclstr.seedCrysIPhiOrIy();
        const float x_calo = supclstr.seed()->position().x();
        const float y_calo = supclstr.seed()->position().y();
        const float z_calo = supclstr.seed()->position().z();
        //const float esEnergyOverRawE = supclstr.preshowerEnergy()/scEnergyRaw;
        const float etaWidth = supclstr.etaWidth();
        const float phiWidth = supclstr.phiWidth();

        if( ERHODEBUG ) std::cout << " --- Finding Moments & Covariences : " << supclstr << std::endl;
        const auto &seedDetId = supclstr.seed()->seed();// seed detid
        const auto isEB = (seedDetId.subdetId() == EcalBarrel);// which subdet
        const auto ph2ndMoments = getCluster2ndMoments( &supclstr );
        const auto lCov = getCovariances( &supclstr );

        if( ERHODEBUG ) std::cout << " --- Storing Moments & Covariences : " << supclstr << std::endl;
        const float scSMaj = ph2ndMoments.sMaj;
        const float scSMin = ph2ndMoments.sMin;
        const float scSAlp = ph2ndMoments.alpha;
        const float scCovEtaEta = lCov[0];
        const float scCovEtaPhi = lCov[1];
        const float scCovPhiPhi = lCov[2];

        Branches.fillBranch("zscIsEB",isEB);
        Branches.fillBranch("zscEnergyRaw",scEnergyRaw);

        Branches.fillBranch("zscIsScEtaEB",isScEtaEB);
        Branches.fillBranch("zscIsScEtaEE",isScEtaEE);
        Branches.fillBranch("zscSeediEtaOriX",seediEtaOriX);
        Branches.fillBranch("zscSeediPhiOriY",seediPhiOriY);
        Branches.fillBranch("zscx_calo",x_calo);
        Branches.fillBranch("zscy_calo",y_calo);
        Branches.fillBranch("zscz_calo",z_calo);

        Branches.fillBranch("zscEtaWidth",etaWidth);
        Branches.fillBranch("zscPhiWidth",phiWidth);
        Branches.fillBranch("zscSMaj",scSMaj);
        Branches.fillBranch("zscSMin",scSMin);
        Branches.fillBranch("zscSAlp",scSAlp);
        Branches.fillBranch("zscCovEtaEta",scCovEtaEta);
        Branches.fillBranch("zscCovEtaPhi",scCovEtaPhi);
        Branches.fillBranch("zscCovPhiPhi",scCovPhiPhi);

		if( ( isOOT || ( time < -1.0 ) ) && scenergy > 5.0 ) nNpSupClstrs++;

    }//<<>>for ( int it = 0; it < nSupClstrs; it++ )
	//std::cout << " -- nNpSC : " << nNpSupClstrs << std::endl;
	geVar.set("nNpSC",nNpSupClstrs);

    if( ERHODEBUG ) std::cout << " - enetering RecHit loop" << std::endl;
    int nRecHits = frechits.size();
	float nUsed = 0;
	//std::cout << " - enetering RecHit loop with " << cfPrm("minRHEf") << std::endl;
    for ( int it = 0; it < nRecHits; it++ ){

        auto recHit = frechits[it];
        bool used = frhused[it];
		if( used ) nUsed++;
		//if( recHit.energy() < 1.0 ) std::cout << " -- checking rh " << getRawID(recHit) << " w/ " << recHit.energy() << " - " << used << std::endl;
		//if( not ( not used && recHit.energy() < cfPrm("minRHEf") ) ) if( recHit.energy() < 1.0 ) std::cout << " --- cut passed old " << std::endl;
		if( recHit.energy() < cfPrm("minRHEi") ) continue;
        if( ( not used ) &&  ( recHit.energy() < cfPrm("minRHEf") ) ) continue;
        //if( recHit.energy() < 1.0 ) std::cout << " --- cut passed  rh " << getRawID(recHit) << " w/ " << recHit.energy() << std::endl;

        const auto recHitID = getRawID(recHit);
        const bool isEB = getIsEB(recHit); // which subdet
        const auto geometry( isEB ? barrelGeometry : endcapGeometry );
        const auto recHitPos = geometry->getGeometry(recHit.detid())->getPosition();
        const float eta = recHitPos.eta();
        const float phi = recHitPos.phi();
        const float rhX = recHitPos.x();
        const float rhY = recHitPos.y();
        const float rhZ = recHitPos.z();
        const float bsX = beamSpot_->x0();
        const float bsY = beamSpot_->y0();
        const float bsZ = beamSpot_->z0();
		//std::cout << " -- BS: " << beamSpot_->sigmaZ() << " : " << beamSpot_->BeamWidthX() << " : " << beamSpot_->BeamWidthY() << std::endl; 
        const float d_rh = hypo(rhX-bsX,rhY-bsY,rhZ-bsZ)/SOL;
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

		const bool hasGS1 = recHit.checkFlag(EcalRecHit::kHasSwitchToGain1);
        const bool hasGS6 = recHit.checkFlag(EcalRecHit::kHasSwitchToGain6);
        const float rhTime = recHit.time();
        const bool rhIsOOT = recHit.checkFlag(EcalRecHit::kOutOfTime);
        const bool rhIsWrd = recHit.checkFlag(EcalRecHit::kWeird);
        const bool rhIsDiWrd = recHit.checkFlag(EcalRecHit::kDiWeird);
        const bool rhIsRecov = recHit.checkFlag(EcalRecHit::kLeadingEdgeRecovered) || recHit.checkFlag(EcalRecHit::kNeighboursRecovered) || recHit.checkFlag(EcalRecHit::kTowerRecovered);
        const bool rhIsPoor = recHit.checkFlag(EcalRecHit::kPoorReco) || recHit.checkFlag(EcalRecHit::kFaultyHardware) || recHit.checkFlag(EcalRecHit::kNoisy) || recHit.checkFlag(EcalRecHit::kPoorCalib) || recHit.checkFlag(EcalRecHit::kSaturated);
        const bool rhIsDead = recHit.checkFlag(EcalRecHit::kDead) || recHit.checkFlag(EcalRecHit::kKilled);
        const bool rhIsOther = recHit.checkFlag(EcalRecHit::kTPSaturated) || recHit.checkFlag(EcalRecHit::kL1SpikeFlag);
        //const bool rhIsNotGood = not recHit.checkFlag(EcalRecHit::kGood);

        const float rhEnergy = recHit.energy();
		float rhJErr = recHit.timeError();
        float rhChi2 = recHit.chi2();
		bool isTimeValid = rhJErr > 0.61;

		if( cfFlag("doProbeOut") ){
			//bool flagged = rhIsRecov || rhIsPoor || rhIsDead || rhIsOther;
			if( rhIsRecov ) 
			std::cout << "Recovered " << " time " << rhTime << " err " << rhJErr << " chi2 " << rhChi2 << " e " << rhEnergy << std::endl;
        	//if( rhIsPoor ) 
        	//std::cout << "RecHitIsPoorCali : " << " time " << rhTime << " err " << rhJErr << " chi2 " << rhChi2 << std::endl;
			//if( flagged ) 
			//std::cout << "Rec " << rhIsRecov << " Poor " << rhIsPoor << " Dead " << rhIsDead << " Other " << rhIsOther;
			//if( rhIsNotGood ) 
			//std::cout << "RecHitIsNotGood : " << " time " << rhTime << " err " << rhJErr << " chi2 " << rhChi2 << std::endl;
			//if( hasGS1 || hasGS6 ) 
			//std::cout << "GS1 " << hasGS1 << " GS6 " << hasGS6 << " err " << rhJErr << " chi2 " << rhChi2 << std::endl;
			//if( rhJErr == 10000 || rhJErr == 0 ) 
			//if( rhJErr == 0 || rhJErr == 10000 )
			//	std::cout << "JerrS " << " time " << rhTime << " err " << rhJErr << " chi2 " << rhChi2 << std::endl;
            //if( rhJErr < 1 )
            //    std::cout << "JerrC " << " time " << rhTime << " err " << rhJErr << " chi2 " << rhChi2 << std::endl;
		}//<<>>if( cfFlag("doProbeOut")

        Branches.fillBranch("ID",recHitID);
        Branches.fillBranch("TOFpv",d_pv);
        Branches.fillBranch("TOF0",d_rh);
        Branches.fillBranch("Time",rhTime);
        Branches.fillBranch("TimeError",rhJErr);
        Branches.fillBranch("isTimeValid",isTimeValid);
        Branches.fillBranch("isOOT",rhIsOOT);
        Branches.fillBranch("isWrd",rhIsWrd);
        Branches.fillBranch("isDiWrd",rhIsDiWrd);
        Branches.fillBranch("isRecov",rhIsRecov);
        Branches.fillBranch("isPoor",rhIsPoor);
        Branches.fillBranch("isDead",rhIsDead);
        Branches.fillBranch("isOther",rhIsOther);
        Branches.fillBranch("hasGS1",hasGS1);
        Branches.fillBranch("hasGS6",hasGS6);
        Branches.fillBranch("Energy",rhEnergy);
        Branches.fillBranch("SwCross",swisscross);
        Branches.fillBranch("eta",eta);
        Branches.fillBranch("phi",phi);
        Branches.fillBranch("rhx",rhX);
        Branches.fillBranch("rhy",rhY);
        Branches.fillBranch("rhz",rhZ);
		Branches.fillBranch("amplitude",amplitude);
		Branches.fillBranch("ampres",ampres);

    }//<<>>for (const auto recHit : *recHitsEB_ )
	float pUsed = nUsed/float(nRecHits); 
    Branches.fillBranch("pused",pUsed);
	//std::cout << " -- % rechits used : " << pUsed << std::endl;

}//<<>>void KUCMSEcalRecHit::ProcessEvent()

void KUCMSEcalRecHitObject::EndJobs(){}

void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup){

    int nRecHits = frechits.size();
    for( auto rhid : idgroup ){
        for( int it = 0; it < nRecHits; it++ ){ if( getRawID( frechits[it] ) == rhid ){ frhused[it] = true; break; } }
    }//<<>>for( auto rhid : phoRhIdsGroup )

};//<<>>void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup)

void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup, bool isOOT, float sce, std::vector<float> fracs ){

	//std::cout << " - enetering setRecHitUsed loop with sce " << sce << " isOOT:  " <<  isOOT << std::endl;
    int nRecHits = frechits.size();
	int idxRhid = 0;
	int nRhids = idgroup.size();
    for( auto rhid : idgroup ){
        for( int it = 0; it < nRecHits; it++ ){ 
			if( getRawID( frechits[it] ) == rhid ){ 
				frhused[it] = true;
				if( frechits[it].energy() < 0.2 ){ 
					std::cout << " - In setRecHitUsed loop with sce " << sce << " isOOT:  " << isOOT << " #rhs: " << nRhids;
					std::cout << " --- rh energy: " << frechits[it].energy() << " frac: " << fracs[idxRhid] << std::endl;  
				}//<<>>if( frechits[it].energy() < 0.2 )
				break; 
			}//<<>>if( getRawID( frechits[it] ) == rhid ) 
		}//<<>>for( int it = 0; it < nRecHits; it++ ){
		idxRhid++;
    }//<<>>for( auto rhid : idgroup )

};//<<>>void KUCMSEcalRecHitObject::setRecHitUsed( rhIdGroup idgroup)

float KUCMSEcalRecHitObject::getSuperClusterSeedTime( reco::SuperClusterRef supercluster ){

    const auto & seedDetId = supercluster->seed()->seed(); // get seed detid
    const auto recHits = ((seedDetId.subdetId() == EcalSubdetector::EcalBarrel) ? recHitsEB_ : recHitsEE_); // which recHits to use
    const auto seedHit = recHits->find(seedDetId); // get the underlying rechit
    const auto seedTime = ((seedHit != recHits->end()) ? seedHit->time() : -9999.f);
    //if( DEBUG && seedTime == -9999.f ) std::cout << "Bad Photon seed time !!!! " << std::endl;
    return seedTime;

}//<<>>getSuperClusterSeedTime( reco::SuperClusterRef supercluster )

std::vector<int> KUCMSEcalRecHitObject::getSuperClusterIndex( const rhIdGroup crystals, int objectIdx ){

    int fscIter = -1;
    std::vector<int> retIndex{-1};
	std::vector<int> prev;
    //std::cout << " -- SC Track matching Trk Indx : " << objectIdx << " with " << crystals.size() << " crystals" << std::endl;
	fscIter++;
	for( const auto cryst : crystals ){

		//std::cout << " ---- Crystal  : " << cryst << std::endl;;
		if( scRechitMap.find(cryst) != scRechitMap.end() ){
	
			int scIndex = scRechitMap[cryst];
            //std::cout << " ----- Matched to : " << scIndex << std::endl;
			if( std::find( prev.begin(), prev.end(), scIndex ) != prev.end() ) continue;
			prev.push_back(scIndex);
			if( fscTrackIndx[scIndex][0] == -1 ) fscTrackIndx[scIndex][0] = objectIdx;
			else fscTrackIndx[scIndex].push_back(objectIdx);
			if( retIndex[0] == -1 ) retIndex[0] = scIndex;
			else retIndex.push_back(scIndex);

		}//<<>>if( loc != scRechitMap.end() )

	}//<<>>for( const auto cryst : crystals )
    //std::cout << " ------------------------------------------ SC matched Idx: "; 
	//for( auto & x : retIndex ){ std::cout << x << " "; }  
	//std::cout << " in: " << objectIdx << std::endl;
    return retIndex;

}//<<>>float KUCMSEcalRecHitObject::getSuperClusterIndex( reco::SuperCluster supercluster )

int KUCMSEcalRecHitObject::getSuperClusterIndex( const reco::SuperCluster* supercluster, int objectID, int objectIdx ){

	const auto & scSeedDetId = supercluster->seed()->seed(); // get seed detid	
	const auto scSeedRawID = scSeedDetId.rawId();
	int fscIter = 0;
	int retIndex = -1;
    //std::cout << " -- SC matching : " << scSeedRawID << " eta: " << osceta << std::endl;
	for( const auto & fsclsrt : fsupclstrs ){ 
		
		const auto & fscSeedDetId = fsclsrt.seed()->seed(); // get seed detid 
		const auto fscSeedRawID = fscSeedDetId.rawId();
		//std::cout << " --- SC compared to : " << scSeedRawID << " ?= " << fscSeedRawID << std::endl;
		if( scSeedRawID == fscSeedRawID ){ 
			retIndex = fscIter; 
			fscOType[fscIter] += objectID;
			if( objectID == 22 ) fscPhoIndx[fscIter] = objectIdx;
            if( objectID == 11 ) fscEleIndx[fscIter] = objectIdx;
			//std::cout << " ----- SC Match Found : " << scSeedRawID << " for: " << objectID << " givin: " << fscOType[fscIter]; 
			//std::cout << " pho: " << fscPhoIndx[fscIter] << " ele: " << fscEleIndx[fscIter] << std::endl; 
			break; 
		}//<<>>if( scSeedRawID == fscSeedRawID )
		fscIter++;		

	}//<<>>for( sclsrt : fsupclstrs )		
	//std::cout << " ------------------------------------------ SC matched Idx: " << retIndex  << " in: " << objectID << std::endl;
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

rhGroup KUCMSEcalRecHitObject::getRHGroup( const reco::CaloCluster basicCluster, float minenr ){

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
    //std::vector<uInt> rawIds;
    for ( const auto &superCluster : superClusterGroup ){
		//std::vector<uInt> rawIds;
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
			//const auto frac = hitsAndFractions[iHAF].second; 
            const auto rawId = detId.rawId();
            //if( ERHODEBUG ) std::cout << " h&f : rawid " << rawId << " frac : " << frac << std::endl;

            for (const auto &recHit : frechits ){
                const auto recHitId = recHit.detid();
                const uInt rhRawId = recHitId.rawId();
                if( rhRawId == rawId ){ result.push_back(recHit); break; } //std::cout << " ---- Found " <<  std::endl; break; }
            }//<<>>for (const auto &recHit : recHits ){
        }//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
    }//<<>>for ( const auto superCluster : superClusterGroup )

	//std::cout << " -------------------------------------------------------------- " << std::endl;
    return result;

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup( const scGroup superClusterGroup, float minenr = 0.0 )

std::vector<float> KUCMSEcalRecHitObject::getFractionList( const scGroup superClusterGroup, rhGroup recHits ){

    std::vector<float> fracs;
	for (const auto &recHit : recHits ){
		const auto recHitId = recHit.detid();
		const uInt rhRawId = recHitId.rawId();
    	for ( const auto &superCluster : superClusterGroup ){
        	auto & hitsAndFractions = superCluster.hitsAndFractions();
        	const uInt nHAF = hitsAndFractions.size();
        	for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            	const auto scDetId = hitsAndFractions[iHAF].first;
				const float frac = hitsAndFractions[iHAF].second;	
            	const uInt scRawId = scDetId.rawId();
				if( rhRawId == scRawId ){ fracs.push_back(frac); break; }
			}//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ) 
		}//<<>>for ( const auto &superCluster : superClusterGroup )
    }//<<>>for (const auto &recHit : recHits )
	return fracs;

}//<<>>std::vector<float> KUCMSEcalRecHitObject::getFractionList( const scGroup superClusterGroup, rhGroup recHits )

rhIdGroup KUCMSEcalRecHitObject::getRhRawIdList( const scGroup superClusterGroup ){

    rhIdGroup list;
    for ( const auto &superCluster : superClusterGroup ){
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
			list.push_back(rawId);
		}//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
	}//<<>>for ( const auto &superCluster : superClusterGroup )
    return list;

}//<<>>std::vector<float> KUCMSEcalRecHitObject::getRhRawIdList( const scGroup superClusterGroup )

std::vector<float> KUCMSEcalRecHitObject::getMissFractionList( const scGroup superClusterGroup ){

	std::vector<uInt> rawIds;
    std::vector<uInt> rawIdsDup;
    std::vector<float> dupcnt;
    for ( const auto &superCluster : superClusterGroup ){
        auto & hitsAndFractions = superCluster.hitsAndFractions();
        const auto nHAF = hitsAndFractions.size();
        for( uInt iHAF = 0; iHAF < nHAF; iHAF++ ){
            const auto detId = hitsAndFractions[iHAF].first;
            const auto rawId = detId.rawId();
			if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() ){
				rawIds.push_back(rawId);
			} else {
				if( std::find( rawIdsDup.begin(), rawIdsDup.end(), rawId ) == rawIdsDup.end() ){
					rawIdsDup.push_back(rawId);
					float cnt = 0;
					for (const auto &recHit : frechits ){

                		const auto recHitId = recHit.detid();
                		const uInt rhRawId = recHitId.rawId();
						if( rhRawId == rawId ){
							dupcnt.push_back(recHit.energy());
							cnt += 1;
						}//<<>>if( rhRawId == rawId )
					}//<<>>for (const auto &recHit : frechits )	
					if( cnt > 1 ) std::cout << " More than 1 matching rechit ? " << std::endl;

				}//<<>>if( std::find( rawIdsDup.begin(), rawIdsDup.end(), rawId ) == rawIdsDup.end() )
			}//<<>>if( std::find( rawIds.begin(), rawIds.end(), rawId ) == rawIds.end() )
		}//<<>>for( uInt iHAF = 0; iHAF < nHAF; iHAF++ )
	}//<<>>for ( const auto &superCluster : superClusterGroup )
    return dupcnt;

}//<<>>std::vector<float> KUCMSEcalRecHitObject::getMissFractionList( const scGroup superClusterGroup, rhGroup recHits )

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

float KUCMSEcalRecHitObject::getRecHitEnergy( const uInt id ){

    if ( id != 0 ){
		for( auto rh : frechits ){
			if( getRawID(rh) == id ) return rh.energy();
    	}//<<>> for( auto rh : frechits )
	}//if ( id != 0 )
    return -999.f;

}//<<>>float EcalTools::recHitE( const DetId id, const EcalRecHitCollection &recHits )

float KUCMSEcalRecHitObject::getRecHitTime( const uInt id ){

    if ( id != 0 ){
        for( auto rh : frechits ){
            if( getRawID(rh) == id ) return rh.time();
        }//<<>> for( auto rh : frechits )
    }//if ( id != 0 )
    return -999.f;

}//<<>>float EcalTools::recHitE( const DetId id, const EcalRecHitCollection &recHits )

rhGroup KUCMSEcalRecHitObject::getRHGroup(){

    rhGroup result;
    for (const auto &recHit : *recHitsEB_ ) result.push_back(recHit);
    for (const auto &recHit : *recHitsEE_ ) result.push_back(recHit);

    return result;

}//>>>>rhGroup KUCMSEcalRecHitObject::getRHGroup()

rhGroup KUCMSEcalRecHitObject::getRHGroup( float eta, float phi, float drmin, float minenr ){

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
    //return EcalClusterTools::covariances( *(scptr->seed()), recHits, ecalTopology, ecalGeometry );
	return EcalClusterTools::localCovariances( *(scptr->seed()), recHits, &(*ecalTopology));

}//<<>>std::array<float, 3> KUCMSEcalRecHitObject::getCovariances( reco::SuperCluster* scptr )

#endif
