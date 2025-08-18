//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

#include "KUCMSItemManager.hh"
#include "KUCMSBranchManager2.hh"

#include "RestFrames/RestFrames.hh"
//#include "BayesianClustering/BayesianClustering.hh"

#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"
#include "ecal_config/KUCMSEcalDetIDFunctions.hh"
//#include "ecal_config/KUCMSTimeCalibration.hh"

#include "KUCMSNtuple_base_v30.h"
//#include "KUCMSNtuple_base_noSV_v30.h"

#ifndef KUCMSAodSkimmer_header
#define KUCMSAodSkimmer_header
//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSAodSkimmer class ---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

using namespace RestFrames;

#define n1dHists 512
#define n2dHists 512 
#define n3dHists 16
#define nEBEEMaps 36

struct globalEventVars {

    ItemManager<std::vector<float>> fvects;
    ItemManager<uInt> cnts;
    ItemManager<float> fvars;

};

//class KUCMSAodSkimmer : public root_base_v16 {
class KUCMSAodSkimmer : public llpgtree {

    public:

    KUCMSAodSkimmer();
	~KUCMSAodSkimmer();

	// tchian processing functions
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo );
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect );
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect, int skipCnt );
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect, bool doSVs, int skipCnt );
    void kucmsAodSkimmer( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename, bool hasGenInfo, bool genSigPerfect, bool doSVs, int skipCnt, bool useEvtWgts );
    void initHists();
    bool eventLoop( Long64_t entry );
	void startJobs();
    void endJobs();
	void fillConfigTree( TTree* fOutTree );
    void setOutputBranches( TTree* fOutTree );

	// object processing & selection
	void processEvntVars();
	void processRechits();
    void processRechitsHist();
	void processGenParticles();
	void processCalojets();
	void processPhotons();
	void processElectrons();
	void processMuons();
    void processJets();
	void processMet();
	void processRJR( int type, bool newEvent );

	int getPhoQuality( int it );
    int getJetQuality( int it );	

	// event processing and selection
	bool eventSelection();

	// histogram and recdhit map varibles & functions
    TH1D *hist1d[n1dHists];
    TH2D *hist2d[n2dHists];
    TH3D *hist3d[n3dHists];

    int nMaps;
    bool fMap[nEBEEMaps];
    TH2D *ebeeMapP[nEBEEMaps], *ebeeMapT[nEBEEMaps], *ebeeMapR[nEBEEMaps];
    void makeEBEEMaps( int phoit );
    void makeEBEEMaps( std::vector<unsigned int> rhcol );
	bool useEvtWgt;

    // aod skimmer helper functions & varibles

    std::map<UInt_t,DetIDStruct> DetIDMap;
	//KUCMSTimeCalibration* timeCali;

    int getRhIdx( uInt rhDetID );
    uInt getLeadRhID( std::vector<uInt> recHitIds );
    float clstrR9( std::vector<uInt> recHitIds );
	std::vector<float> getRhGrpTimes( std::vector<uInt> rechitids );
    std::vector<float> getRhGrpEnergies( std::vector<uInt> rechitids );
	std::vector<float> getRhGrpEigenFromAngles( std::vector<uInt> rechitids );
    std::vector<float> getLeadTofRhTime( std::vector<uInt> recHitIds, double vtxX, double vtxY, double vtxZ );
    std::vector<float> getRhGrpEigen_sph( std::vector<float> times, std::vector<uInt> rechitids );

  	// RestFrames frames and friends

  	LabRecoFrame* LAB;
  	DecayRecoFrame* S;
	DecayRecoFrame*	X2a; 
	DecayRecoFrame* X2b;
  	DecayRecoFrame* Ja; 
    DecayRecoFrame* Jb;
    VisibleRecoFrame* J1a;
    VisibleRecoFrame* J1b;
    VisibleRecoFrame* J2a;
    VisibleRecoFrame* J2b;
  	InvisibleRecoFrame* X1a; 
    InvisibleRecoFrame* X1b;

  	InvisibleGroup* INV;
  	SetMassInvJigsaw* InvM;
  	SetRapidityInvJigsaw* InvEta;
  	MinMassesSqInvJigsaw* InvSplit;

  	CombinatoricGroup* COMB_J;
  	MinMassesSqCombJigsaw* CombSplit_J;
    MinMassesCombJigsaw* CombSplit_Ja;
    MinMassesCombJigsaw* CombSplit_Jb;

    // config vars
    std::string dataSetKey;
    float xsctn;
    float gmass;
	float xmass;
	float mcwgt;
	int mctype;

	// event varibles

    bool doGenInfo;
	globalEventVars gEvtVar;
    ItemManager<std::vector<float>> geVects;	
    ItemManager<uInt> geCnts;
	ItemManager<float> geVars;
	std::map<std::string,int> cutflow;
    std::map<std::string,int> configCnts;
    std::map<std::string,float> configWgts;
	uInt nEvents, nSelectedEvents;
    float sumEvtGenWgt;
    std::vector<bool> phoJetVeto;
	std::vector<bool> rhused;
    std::vector<float> rhispho;
    std::vector<float> rhisele;

    // Output Branch variables

	KUCMSBranchManager selEvtVars;
	KUCMSBranchManager selMet;
    KUCMSBranchManager selElectrons;
	KUCMSBranchManager selPhotons;
    KUCMSBranchManager selMuons;
	KUCMSBranchManager selJets;
	KUCMSBranchManager selRjrVars;
    KUCMSBranchManager selGenPart;


   std::vector<float>   *OSuperCluster_diffEnergy;
   std::vector<float>   *OSuperCluster_dR;
   std::vector<unsigned int> *OSuperCluster_otherMatchSeedID;
   Int_t           OSuperCluster_nOther;
   Int_t           OSuperCluster_nOtherEx;
   Int_t           OSuperCluster_nOtherIn;
   std::vector<unsigned int> *OSuperCluster_otherSeedID;
   std::vector<int>     *OSuperCluster_nXtalOverlap;


   TBranch        *b_OSuperCluster_diffEnergy;   //!
   TBranch        *b_OSuperCluster_dR;   //!
   TBranch        *b_OSuperCluster_otherMatchSeedID;   //!
   TBranch        *b_OSuperCluster_nOther;   //!
   TBranch        *b_OSuperCluster_nOtherEx;   //!
   TBranch        *b_OSuperCluster_nOtherIn;   //!
   TBranch        *b_OSuperCluster_otherSeedID;   //!
   TBranch        *b_OSuperCluster_nXtalOverlap;   //!


};

#endif
