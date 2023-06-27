// -*- C++ -*-
//
//
// Original Author:  Jack W King III
// 
// KUCMS Object for Object Manager
//

//--------------------   hh file -------------------------------------------------------------
//--------------------------------------------------------------------------------------------

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

#ifndef KUCMSObjectBaseHeader
#define KUCMSObjectBaseHeader

typedef edm::SortedCollection<EcalRecHit,edm::StrictWeakOrdering<EcalRecHit>> recHitCol;
typedef vector<EcalRecHit> rhGroup;
typedef vector<reco::SuperCluster> scGroup;
typedef vector<reco::CaloCluster> bcGroup;
typedef vector<uInt> rhIdGroup;

enum class ECAL {EB, EM, EP, EE, NONE};

//class exampleObject : public KUCMSObjectBase {};

class KUCMSEcalRecHit : public KUCMSObjectBase {

    // use base class constructor
    using KUCMSObjectBase::KUCMSObjectBase;

    public:

    // object setup : 1) construct object 2) InitObject 3) CrossLoad 4) load into Object Manager
    void InitObject( ItemManager<float>& geVar, TTree* fOutTree ); // load geVars and sets up branches, do preloop jobs
    // new function needed for crosstalk - EXAMPLE CLASS USED HERE FOR REFRENCE ONLY -
    // void LoadObject( exampleObject& otherObject ){ otherObjectRef = otherObject; }; // define with specific KUCMS object(s) needed 

    // object processing : 1) LoadEvent prior to event loop 2) ProcessEvent during event loop via objectManager
    void LoadEvent( const edm::Event& iEvent, const edm::EventSetup& iSetup ); // get collections, do initial processing
    void ProcessEvent(); // do cross talk jobs with other objects, do event processing, and load branches

    // if there are any final tasks be to done after the event loop via objectManager
    void EndJobs(); // do any jobs that need to be done after main event loop

    // New functions specific to this collection
    // void answerCrossTalk(); // define functions that will be called in another object - this is an example
    // ect ...

    private: // example from photons

        std::vector<EcalRecHit>         frechits;
        std::vector<bool>               frhused;

    // Other object(s) need by this object - BASE CLASS USED HERE FOR REFRENCE ONLY -
    //std::reference_wrapper<exampleObject> otherObjectRef;

};

void KUCMSEcalRecHit::InitObject( ItemManager<float>& geVar, TTree* fOutTree ){

	

}

#endif
