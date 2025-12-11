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

//#include "KUCMSItemManager.hh"
#include "KUCMSBranchManager2.hh"

#include "RestFrames/RestFrames.hh"
#include "BayesianClustering/BayesianClustering.hh"

//#include "KUCMSHelperBaseClass.hh"
//#include "KUCMSRootHelperBaseClass.hh"
#include "KUCMSHelperFunctions.hh"
#include "KUCMSRootHelperFunctions.hh"


//#include "KUCMSEcalDetIDFunctions.hh"
#include "KUCMS_TimeCalibration.hh"

#include "KUCMSNtuple_base_v30.h"
//#include "KUCMSNtuple_base_v30_noSV.h"

#include "MaterialVeto.h"

#include<Math/Vector3D.h>
#include<Math/Vector3Dfwd.h>

#ifndef KUCMSAodSkimmer_header
#define KUCMSAodSkimmer_header

//--------------------------------------------------------------------------------------------------------------------------------------
// KUCMSAodSkimmer class ---------------------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------------------------------------------

//#define DEBUG true
#define DEBUG false


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

    void kucmsAodSkimmer( std::string infile, std::string outfilename);
    void kucmsAodSkimmer_Filelist( std::string eosdir, std::string infilelist, std::string outfilename);
    void kucmsAodSkimmer_listsOfLists( std::string eosdir, std::string infilelist, std::string outfilename );
    void kucmsAodSkimmer_local( std::string listdir, std::string eosdir, std::string infilelist, std::string outfilename );

    int ProcessFilelistOfLists(string eosdir, vector<string> processed_strings, TChain*& fInTree, TChain*& fInConfigTree);
    int ProcessFilelist(string eosdir, string infilename, TChain*& fInTree, TChain*& fInConfigTree);
    int ProcessFile(string infilename, TChain*& fInTree, TChain*& fInConfigTree);

    void ProcessMainLoop( TChain* fInTree, TChain* fInConfigTree );
    void ProcessConfigTree( TChain* fInConfigTree );
    void fillConfigTree( TTree* fOutTree );

    void initHists();
    bool eventLoop( Long64_t entry );
	void startJobs();
    void endJobs();

   // Class varibles accessors

    void SetEventRange(int evti, int evtj){ _evti = evti; _evtj = evtj;}
    void SetDataSetKey(string key){ dataSetKey = key; }
    void SetCrossSection(float xsec){ xsctn = xsec; }
    void SetGluinoMass(float glumass){ gmass = glumass; }
    void SetN2Mass(float n2mass){ xmass = n2mass; }
    void SetMCType(int type){ mctype = type; }
    void SetTimeCalibrationTag(string ttag){ tctag = ttag; }
    void SetMCWeight(float w){ mcwgt = w; }

	void SetGenInfoFlag( bool f ){ hasGenInfoFlag = f; } 
    void SetNoSVorPhoFlag( bool f ){ noSVorPhoFlag = f; }
    void SetUseEvtGenWgtFlag( bool f ){ useEvtGenWgtFlag = f; }
    void SetGenSigPerfectFlag( bool f ){ genSigPerfectFlag = f; }
    void SetDoBHC( bool f ){ doBHC = f; }
    void SetOutFileName( std::string  ofn ){ outFileName = ofn; }
    void SetDoSV( bool f ){ doSVs = f; }

	// set branches

    void setOutputBranches( TTree* fOutTree );
	void setEvtVarMetBranches( TTree* fOutTree );
	void setTrackBranches( TTree* fOutTree );
	void setEcalBranches( TTree* fOutTree );
	void setGenBranches( TTree* fOutTree );
	void setPhotonBranches( TTree* fOutTree );
	void setRJRBranches( TTree* fOutTree );
	void setRJRISRBranches( TTree* fOutTree );
	void setElectronBranches( TTree* fOutTree );
	void setJetsBranches( TTree* fOutTree );
	void setMuonsBranches( TTree* fOutTree );
	void setSVBranches( TTree* fOutTree );
	void setBCBranches( TTree* fOutTree );
	void setTSBranches( TTree* fOutTree );

	// object processing & selection

    void processTracks();
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
    void processSV();

	void processRJR( int type, bool newEvent );
	void processRJRISR();
	void processBHCPhotons();
    void processBHCJets();
	void processTimeSig();

	// TSig functions
	//float getTimeSig( int scIndex );
	float getTimeSig( int scIndex, float& num, float& denom, const map<unsigned int, float>& rhIdToBHCw  = {});
	float getTimeSig( vector<vector<unsigned int>> rhids );

	void MakePhotonIsoMap(int phoidx, map<string, double>& isomap);
	// object quality functions

	int getPhoQuality( int it );
    int getJetQuality( int it );	

	// good Lumisection check
	std::map<int, std::vector<std::map<int, int>>> lumiData;
	void loadLumiJson( std::string infile, bool verbose = false );
	bool isValidLumisection( int run, int lumi );

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

    //std::map<UInt_t,DetIDStruct> DetIDMap;
	KUCMS_TimeCalibration* timeCali;

    bool inHEMRegion( float eta, float phi );
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
    //VisibleRecoFrame* Ja;
    //VisibleRecoFrame* Jb;
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

    // ISR RestFrames frames and friends

    LabRecoFrame* LAB_c;
    DecayRecoFrame* CM_c;
    DecayRecoFrame* S_c;
    DecayRecoFrame* X2a_c;
    DecayRecoFrame* X2b_c;
    DecayRecoFrame* Ja_c;
    DecayRecoFrame* Jb_c;
    VisibleRecoFrame* ISR_c;
    VisibleRecoFrame* J1a_c;
    VisibleRecoFrame* J1b_c;
    VisibleRecoFrame* J2a_c;
    VisibleRecoFrame* J2b_c;
    InvisibleRecoFrame* X1a_c;
    InvisibleRecoFrame* X1b_c;

    InvisibleGroup* INV_c;
    SetMassInvJigsaw* InvM_c;
    SetRapidityInvJigsaw* InvEta_c;
    MinMassesSqInvJigsaw* InvSplit_c;

    CombinatoricGroup* COMB_J_c;
    MinMassesCombJigsaw* CombSplit_ISR_c; 	
    MinMassesSqCombJigsaw* CombSplit_J_c;
    MinMassesCombJigsaw* CombSplit_Ja_c;
    MinMassesCombJigsaw* CombSplit_Jb_c;

    // BHCish

    ClusterAnalyzer _ca;

    // config vars

    std::string dataSetKey;
    float xsctn;
    float gmass;
	float xmass;
	float mcwgt;
	int mctype;
	std::string tctag;

	bool genSigPerfectFlag;
	bool noSVorPhoFlag;
	bool useEvtGenWgtFlag; 
	bool hasGenInfoFlag;
    bool doBHC;
	bool doSVs;

	// input tree, paths, and file names

	std::string disphoTreeName;
	std::string configTreeName;
    std::string inFilePath;
    std::string inFileName;
    std::string listDirPath;
	std::string eosDirPath;
    std::string outFileName;

	// event varibles

	enum objID { Base, Loose, Tight };
	enum tSigID { Early, Prompt, Late };

	globalEventVars gEvtVar;
    ItemManager<std::vector<float>> geVects;	
    ItemManager<uInt> geCnts;
	ItemManager<float> geVars;

	std::map<std::string,int> cutflow;
    std::map<std::string,int> configCnts;
    std::map<std::string,float> configWgts;

    std::map<UInt_t, pair<int,int>> _detidmap;
	uInt nEvents, nSelectedEvents;
	int _evti, _evtj;
    float sumEvtGenWgt;

	// hem veto 
    bool hasHemObj;
	//std::map< std::string, bool > hemBits;
	ItemManager<bool> hemBits;

    std::vector<bool> phoJetVeto;
	std::vector<bool> rhused;
    std::vector<float> rhispho;
    std::vector<float> rhisele;

	// BC informaiton passing 
    std::vector<float> erh_corTime;
    std::vector<float> erh_timeRes;
	std::vector<bool> isSelPho;
    std::vector<bool> isSelJet;
	std::map<unsigned int, int> rhIDtoIterMap;

    // Output Branch variables

	KUCMSBranchManager selEvtVars;
	KUCMSBranchManager selMet;
    KUCMSBranchManager selElectrons;
	KUCMSBranchManager selPhotons;
    KUCMSBranchManager selMuons;
	KUCMSBranchManager selJets;
	KUCMSBranchManager selRjrVars;
    KUCMSBranchManager selRjrIsrVars;
    KUCMSBranchManager selGenPart;
    KUCMSBranchManager selSV;
    KUCMSBranchManager selTracks;
    KUCMSBranchManager BHCPhoInfo;
    KUCMSBranchManager BHCJetInfo;
    KUCMSBranchManager TSPhoInfo;

   // Other SC information ? ( this is a hack ) depreciated - or make so ?   
/*
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
*/

};

#endif
