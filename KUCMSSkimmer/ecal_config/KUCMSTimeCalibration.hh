//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSRootHelperBaseClass.hh"
#include <TRandom.h>
#include "TChain.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"

#ifndef KUCMSTimeCalibrationClass
#define KUCMSTimeCalibrationClass

//--------------------   hh file -------------------------------------------------------------
//---------------------------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////////////////
//
//  this code assumes that the ecal_config directory is included beside this class
//  ecal_config should contain :
//
//  	fullinfo_detids_EB.txt
//  	fullinfo_detids_EE.txt
//  	fullinfo_v2_detids_EB.txt
//  	fullinfo_v2_detids_EE.txt
//  	howto.txt
//  	reducedinfo_detids.txt
//  	rhid_i12_list.txt
//  	rhid_info_list.txt 
//
/////////////////////////////////////////////////////////////////////////////////////////

struct DetIDStruct {

    DetIDStruct(){}
    DetIDStruct( const int ni1, const int ni2, const int nTT, const int & necal, const float nphi, const float neta )
        : i1(ni1), i2(ni2), TT(nTT), ecal(necal), phi(nphi), eta(neta) {}

    int i1; // EB: iphi, EE: ix
    int i2; // EB: ieta, EE: iy
    int TT; // trigger tower
    int ecal; // EB, EM, EP
    float phi; // xtal phi
    float eta; // xtal eta

};//<<>>struct DetIDStruct

//---------------------------------------------------------------------------------------------------------------------

struct lumiRunStruct {

    lumiRunStruct(){}
    lumiRunStruct( int trun, int tfill, float tlumi )
        : run(trun), fill(tfill), lumi(tlumi) {}

    int run;
    int fill;
    float lumi;

};//<<>>lumiRunStruct

//---------------------------------------------------------------------------------------------------------------------

struct smearTagStruct {

    smearTagStruct(){}
    smearTagStruct( float tnoise, float tstoch, float tstant )
        : noise(tnoise), stoch(tstoch), stant(tstant) {}

    float noise;
    float stoch;
    float stant;

};//<<>>smearTagStruct

//---------------------------------------------------------------------------------------------------------------------

struct caliHistStruct {
    
    caliHistStruct(){}
    caliHistStruct( TH2F* thist, std::string histname, bool isnew, bool reshist )
        : h2f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h1f = NULL; }
    caliHistStruct( TH1F* thist, std::string histname, bool isnew, bool reshist )
        : h1f(thist), histName(histname), isNew(isnew), isResHist(reshist) { h2f = NULL; }

    TH2F* h2f;
    TH1F* h1f;
    std::string histName;
	bool isNew;
	bool isResHist;

};//<<>>struct caliHistStruct

//------------------------------------------------------------------------------------------------------------

struct sumCnt {

    sumCnt(){}
	sumCnt( float s, float s2, int c ) : sum(s), sumsqr(s2), cnt(c) {} 

    float sum;
    float sumsqr;
	int cnt;

};//<<>>struct sumCnt

class CaliRunClass : public KUCMSRootHelperBaseClass {
    
	public: 

    CaliRunClass(){}
    CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi );
	//~CaliRunClass(); 
 
    std::string histMapName; 
    int startRun;
    int endRun; 
    int lastRun;
    float lumi;

    float noise;
    float stoch;
    float stant;

	bool isNew;
	bool updated;
	bool has2DResMap;
	bool hasResParams;
	bool isExternal;

    std::map<uInt,sumCnt> sumCntMap;
    std::map<uInt,float> meanMap;
    std::map<uInt,float> errMap;
	std::map<uInt,TH1F*> detIdHists;
	//detIdHists[cmsswId] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",500,-5,5);	

	void fillSumCnt( uInt detID, float val, int cnt = 1 );
	void makeMeanMap( bool filter = false );

};//<<>>TimeCaliTagStruct

inline CaliRunClass::CaliRunClass( std::string tmpxtalmap, int tstart, int tend, int last, float tlumi )
	: histMapName(tmpxtalmap), startRun(tstart), endRun(tend), lastRun(last), lumi(tlumi) 
	{ isNew = true; isExternal = false; updated = false; has2DResMap = false; hasResParams = false; noise = 0; stoch = 0; stant = 0; }

inline void CaliRunClass::makeMeanMap( bool filter ){

	meanMap.clear();
	errMap.clear();
	//std::cout << " make mean maps: " << std::endl;
	for( auto& entry : sumCntMap ){

		float cnt = float( entry.second.cnt );
		if( cnt == 0 ) cnt = 1;
		float mean = entry.second.sum / cnt; 
		float err = sqrt( (entry.second.sumsqr/cnt - mean*mean)/cnt );
		if( err == 0 ){ err = 9; mean = 0; }
		//std::cout << " - calc: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl;
		if( filter && cnt > 100 ){
			auto fitFunc  = new TF1("gfit","gaus",-3.0,3.0);
			auto thefit = detIdHists[entry.first]->Fit("gfit","QNRL");
			float fmean = fitFunc->GetParameter(1);
			float ferr = fitFunc->GetParError(1);
			//std::cout << " --- fit: result " << thefit << " mean " << fmean << " +/- " << ferr << std::endl;
			if( thefit != 4 && ferr != 0 && ferr < err ){ mean = fmean; err = ferr; }
			delete fitFunc;
		}//<<>>if( filter )
		delete detIdHists[entry.first];
		meanMap[entry.first] = mean;
        errMap[entry.first] = err; 
		//std::cout << " - mmm: " << entry.first << " = " << mean << " +/- " << err << " occ: " << cnt << std::endl;
	}//<<>>for( auto& entry : sumCntMap )

}//<<>>void CaliRunClass::makeMeanMap()

inline void CaliRunClass::fillSumCnt( uInt detid, float val, int cnt ){

	if( endRun == lastRun ) return;
	//std::cout << "Filling " << detid << " with " << val << " " << cnt << std::endl;
	updated = true;
	if( sumCntMap.find(detid) != sumCntMap.end() ){ 
		sumCntMap[detid].sum += val;
        sumCntMap[detid].sumsqr += val*val; 
		sumCntMap[detid].cnt += cnt; 
	} else {
		sumCntMap[detid] = { val, val*val, cnt }; 
	}//<<>>if( sumCntMap.find(detid) != sumCntMap.end() )

	if( detIdHists.find(detid) != detIdHists.end() ){ detIdHists[detid]->Fill(val); } //std::cout << " -- hist Fill " << detid << " with " << val << std::endl; }
	else {
		//std::cout << " -- hist Make " << detid << " with " << val << std::endl;
		std::string histname = histMapName + "_SumCntHist_" + std::to_string( startRun ) + "_" + std::to_string(detid) ;
		detIdHists[detid] = new TH1F(histname.c_str(),"AveXtalTimeDist;XtalTime [ns]",200,-5,5);
		detIdHists[detid]->Sumw2();
		detIdHists[detid]->Fill(val);
	}//<<>>if( detIdHists.find(detid) != detIdHists.end() )

}//<<>>void CaliRunClass::fillSumCnt( uInt detID, float sum, int cnt )

//-------------------------------------------------------------------------------------------------------------

// fit result struct
struct TimeFitResult {

    float mu;
    float emu;
    float chi2ndf;
    float chi2prob;
    float sigma;
    float sigmaHigh;
    float sigmaLow;
    float esigma;
    float occ;
    float rms;
    float std;

};//<<>>struct TimeFitResult

struct SigmaFitResult {

    float noise;//noise
    float enoise;
    float stant;//constant
    float estant;
    float stoch;//stochastic
    float estoch;

};//<<>>struct SigmaFitResult

class ProfileTimeFit {

	public:
	
	ProfileTimeFit(){ profileHist = NULL; form = NULL; fit = NULL; };
  	ProfileTimeFit( TH1F* hist ) : profileHist(hist) { form = NULL; fit = NULL; };
	//~ProfileTimeFit(){ if( profileHist ) delete profileHist; if( form ) delete form; if( fit ) delete fit; };  

  	// internal data members
  	TH1F* profileHist;
  	TFormula* form;
  	TF1* fit;
    TF1* fitHigh;
    TF1* fitLow;
	TimeFitResult results;

    // helper functions for making fits to variables
    //bool isEmpty(){ std::cout << " -- ?  is empty " << profileHist << std::endl; return ( profileHist ) ? profileHist->GetEntries() == 0 : true; };
    void DoFit();
    TimeFitResult GetFitResult(){ return results; };
	void deleteHists(){ if( profileHist ) delete profileHist; if( form ) delete form; if( fit ) delete fit; };

};//<<>>class ProfileTimeFit

inline void ProfileTimeFit::DoFit(){

	std::cout << " -- ProfileTimeFit DoFit " << std::endl;
	//if( isEmpty() ){ std::cout << " --- No such profile hist !!!!!!! " << std::endl; return; }
	//std::cout << " -- yes profile exists " << std::endl;

  	float Sqrt2PI   = 2.506628274631001f;

  	// set tmp init vals
  	auto hsum = profileHist->Integral();
	if( hsum == 0 ){
		std::cout << " -- !!!!!! Intergral is Zero !!!!!!!! " << std::endl;
    	results.mu = 0;
    	results.emu = 0;
    	results.chi2ndf = 0;
    	results.chi2prob = 0;
    	results.sigma = 10;
    	results.sigmaHigh = 0;
    	results.sigmaLow = 0;
    	results.esigma = 10;
    	results.occ = 0;
    	results.rms = 0;
    	results.std = 0;
		return;
	}//<<>>if( hsum == 0 )

  	float mu = profileHist->GetMean();
  	float sigma = profileHist->GetStdDev();
  	float rms = profileHist->GetRMS();
    //std::cout << " -- getting intgral, mean, std 5" << std::endl;
  	//float norm = hsum/(sigma*Sqrt2PI);
    float norm = profileHist->GetBinContent(profileHist->GetMaximumBin());
    std::cout << " -- Params: sum " << hsum << " mu " << mu << " s " << sigma << " rms " << rms << " norm " << norm << std::endl;

    std::cout << " -- setting paramters " << std::endl;
  	// range vars   ? names ?
  	//float sigrange = 2.0;
  	float lowerBound = mu - 2*sigma;
  	float upperBound = mu + 2*sigma;
	float lowerHBound = mu - 4*sigma;
    float upperHBound = mu + 4*sigma;
    float lowerLBound = mu - 1*sigma;
    float upperLBound = mu + 1*sigma;

  	// names for fits and formulas
  	std::string histname = profileHist->GetName();
  	std::string formname = histname+"_formula";
  	std::string fitname  = histname+"_fit";
    std::string fitnameH  = histname+"_fitHigh";
    std::string fitnameL  = histname+"_fitLow";

    //if( lowerBound < -1.0 ) lowerBound = -0.5;
    //if( upperBound > 1.0 ) upperBound = 0.5;
    //if( abs(mu) > 0.03 ) mu = -0.02;
    //if( sigma > 1  ) sigma = 0.3;

    std::cout << " -- setting form and fit " << std::endl;
    form = new TFormula(formname.c_str(),"[0]*exp(-0.5*((x-[1])/[2])**2)");
    fit  = new TF1(fitname.c_str(),form->GetName(),lowerBound,upperBound);
    fitHigh  = new TF1(fitnameH.c_str(),form->GetName(),lowerHBound,upperHBound);
    fitLow  = new TF1(fitnameL.c_str(),form->GetName(),lowerLBound,upperLBound);

    fit->SetParName(0,"N");      
	fit->SetParameter(0,norm); 
	fit->SetParLimits(0,norm/10,norm*10);
    fit->SetParName(1,"#mu");    
	fit->SetParameter(1,mu); 
	//fit->SetParLimits(1,-0.03,0.03);
    fit->SetParLimits(1,lowerBound,upperBound);
    fit->SetParName(2,"#sigma"); 
	fit->SetParameter(2,sigma); 
	//fit->SetParLimits(2,0,1);
    fit->SetParLimits(2,sigma/10,sigma*10);

    fitHigh->SetParName(0,"N");
    fitHigh->SetParameter(0,norm);
    fitHigh->SetParLimits(0,norm/10,norm*10);
    fitHigh->SetParName(1,"#mu");
    fitHigh->SetParameter(1,mu);
    //fitHigh->SetParLimits(1,-0.03,0.03);
    fitHigh->SetParLimits(1,lowerHBound,upperHBound);
    fitHigh->SetParName(2,"#sigma");
    fitHigh->SetParameter(2,sigma);
    //fitHigh->SetParLimits(2,0,1);
    fitHigh->SetParLimits(2,sigma/10,sigma*10);

    fitLow->SetParName(0,"N");
    fitLow->SetParameter(0,norm);
    fitLow->SetParLimits(0,norm/10,norm*10);
    fitLow->SetParName(1,"#mu");
    fitLow->SetParameter(1,mu);
    //fitLow->SetParLimits(1,-0.03,0.03);
    fitLow->SetParLimits(1,lowerLBound,upperLBound);
    fitLow->SetParName(2,"#sigma");
    fitLow->SetParameter(2,sigma);
    //fitLow->SetParLimits(2,0,1);
    fitLow->SetParLimits(2,sigma/10,sigma*10);

    std::cout << " - Fit > Mu: " << mu << " ( " << lowerBound << " to " << upperBound << " )";
    std::cout << " Norm: " << norm << "( " << norm/10 << " to " << norm*10 << " )";
    std::cout  << " Sigma: " << sigma << "( " << sigma/10 << " to " << sigma*10 << " )" << std::endl;

    std::cout << " -- fiting " << std::endl;
  	profileHist->Fit(fit->GetName(),"RBQO");
    profileHist->Fit(fitHigh->GetName(),"RBQO");
    profileHist->Fit(fitLow->GetName(),"RBQO");

    std::cout << " -- setting results " << std::endl;
    results.mu = fit->GetParameter(1);
    results.emu = fit->GetParError (1);
    results.chi2ndf = fit->GetChisquare();
    results.chi2prob = fit->GetProb();
    results.sigma = fit->GetParameter(2);
    results.sigmaHigh = fitHigh->GetParameter(2);
    results.sigmaLow = fitLow->GetParameter(2);
    results.esigma = fit->GetParError (2);
    results.occ = hsum;
    results.rms = rms;
    results.std = sigma;

    std::cout << " - Fit Sigma: " << results.sigma << " +/- " << results.esigma;
    std::cout << " ( low " << results.sigmaLow << " high " << results.sigmaHigh << " ) " << std::endl;

}//<<>>void ProfileTimeFit::DoFit()

//-----------------------------------------------------------------------------------------------------------------------
//-----------------------------------------------------------------------------------------------------------------------

class KUCMSTimeCalibration : public KUCMSRootHelperBaseClass {

	public:

    KUCMSTimeCalibration( bool stayOpen = false, bool makeNew = false  );
	~KUCMSTimeCalibration();


	private:

	enum ECAL {EB, EM, EP, NONE};

    std::string detIDConfigEB;
    std::string detIDConfigEE;

    std::string caliFileDir;
    std::string caliRunConfig;
    std::string caliTTConfig;
    std::string caliTFileName;
    std::string cali2DResPlotsTFileName;
    std::string caliSmearConfig;

    std::string xtalHistMapName;
    std::string ttHistMapName;
    std::string lochist;
    std::string globhist;

    TFile* caliTFile;
    TFile* cali2DResTFile;

	std::string eosDir;
	std::string inDir;

    std::string curTTIov;
    std::string curXIov;
	std::string curTag; // used to change rereco version or campaian calibrations are for
    std::string curLumiTag;

    std::string smearTag;

	std::string xBinStr;
    std::string yBinStr;

	std::vector<std::string> maptypes;

    std::map<UInt_t,DetIDStruct> DetIDMap; // map of information by detid for each crystal
    std::map<int,std::map<int,std::map<int,uInt>>> InvDetIDMap;

	// do we need TT and X verions ?  yes for set and no for hist maps ?
	std::map<std::string,std::map<int,CaliRunClass>> CaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map
    std::map<std::string,std::map<int,CaliRunClass>> TTCaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map

	// interface between TFile and RunMapSet ?  off load to meanMap?
	std::map<std::string,caliHistStruct> CaliHists; // str is label for specific calibration histogram - stored as ref in calirunmap
    std::map<std::string,std::map<int,lumiRunStruct>> lumiRunMaps;// str is settag, int is run - ref for maps - give lumi by run
    std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>

    std::map<std::string,smearTagStruct> SmearTagSet;

	bool updated;
	bool lowEnergy;
	bool useEffEnergy;
    bool externalCali;
    bool useGSwitch;

    TRandom* getRandom;

	public:

	void SetupDetIDsEB();
	void SetupDetIDsEE();
	void SetupIovMaps();

    void SetupIovMap( std::string tag, float lumiMin );

	void ReadCaliRunFile();
    void SaveCaliRunFile();
    void ReadTTRunFile();
    void SaveTTRunFile();
	void ReadSmearFile();
    void SaveSmearFile();
	void ReadLumiFile( std::string lumifile, std::string tag );// only need to create caliRunMap entries - made with brilcalc

    void LoadCaliHists( bool stayOpen = false, bool makeNew = false );// loads up all starting info
    void SaveCaliHists();// save existing calimaphists to root TFile -- used if CaliHits left open
	void LoadExtCali( std::string calihist, std::string mapname, std::string tag, int startr, int endr );

	// use to create calibration files and add to DB
	// work in progress - still thinking this area trhough
	// 
	// ok - load hists from tfile then make maps from hists ( mean, count, sum )
	// 		use &/or add to maps 
	// 		when finshed update hists from maps then save hist to tfile
	//			will delete old tfile and recreate tfile every pass
	//
	void makeCaliMaps(); // make maps from hists
	void makeCaliHists(); // make hists from updated maps
	// ? void makeSmearMaps();

	// accessors to utilize calibration information
	float getCalibration( uInt rhid, int run, std::string tag  ); // tag indicates which calibration set to use
    float getCalibration( uInt rhid, int run )
			{ return getCalibration( rhid, run, curTag ); };
	float getCalibTime( float rhtime, uInt rhid, int run )
			{ return rhtime - getCalibration( rhid, run ); };
	//tag indicates which smear to use
	float getSmearedTime(  float rhtime, float rhamp, std::string stag );
    float getSmearedTime(  float rhtime, float rhamp )
            { return getSmearedTime( rhtime, rhamp, smearTag ); };
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run, std::string ctag, std::string stag );
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run )
			{ return getSmrdCalibTime( rhtime, rhamp, rhid, run, curTag, smearTag ); }; 

	float getTTCali( uInt rhid, int run, std::string tag );

	DetIDStruct& getDetIdInfo( uInt rhid );
	uInt getDetIdInfo( int i1, int i2, int ecal );
	uInt getTTId( uInt detId );
	std::pair<int,int> getTTInfo( uInt ttid );
    uInt getInvTTId( int i1, int i2 );

	void setTag( std::string tag ){ curTag = tag; }; // used to change rereco version or campaian calibrations are for
    void setXIov( std::string tag ){ curXIov = tag; };
    void setTTIov( std::string tag ){ curTTIov = tag; };
    void setLumiTag( std::string tag ){ curLumiTag = tag; };

    void setSmearTag( std::string tag ){ smearTag = tag; };

    void SetXBinStr( std::string xbins ){ xBinStr = xbins; };
    void SetYBinStr( std::string ybins ){ yBinStr = ybins; };
	void SetEosDir( std::string eosdir ){ eosDir = eosdir; }
    void SetInDir( std::string indir ){ inDir = indir; };

	void Set2DResTFileName( std::string name ){ cali2DResPlotsTFileName = name; };
    void SetCaliTFileName( std::string name ){ caliTFileName = name; };
    void SetCaliFileDir( std::string name ){ caliFileDir = name; };

	// function to create calibration and smear information 
	// make TT maps from egamma res ntuple from lpc/jwk space
	void makeCaliMapsEGR( std::string inputFileName, bool doTT, bool small = false, bool doCali = true );
	void makeTTCaliMapEGR( std::string inputFileName, bool small = false, bool doCali = true ) 
			{ makeCaliMapsEGR( inputFileName, true, small, doCali ); }; 
	void makeXCaliMapEGR( std::string inputFileName, bool small = false, bool doCali = true )
			{ makeCaliMapsEGR( inputFileName, false, small, doCali ); };

	void plot2dResolutionEGR( std::string inputFileName, bool scale = true, bool usecali = true, bool smear = false, std::string ext = "" );
    void plot2dResbyIovForEGR( std::string inputFileName, bool scale = true, bool usecali = true, bool smear = false, std::string ext = "" );
    SigmaFitResult runTimeFitter( TH2F* hist2D );
	void doResTimeFits( bool doLocal = false );
    void doResTimeFit( std::string histName );
	void load2DResHist( std::string histName );

    void makeSmearTag( std::string sourceName, std::string destName, std::string smearTag );

	void plotMeanRunTimeEGR( std::string inputFileName, int srun, int erun, bool usecali = true );
    void makeTTDiffMaps( int srun, int erun );
    void makeTTDriftMaps( std::string tag, int srun, int erun );

	void setUseEffEnergy( bool setting ){ useEffEnergy = setting; };
	void setLowEnergy( bool setting ){ lowEnergy = setting; }; 
	void useExtCali( bool ext ){ externalCali = ext; };
    void useGainSwitch( bool gs ){ useGSwitch = gs; };

};//<<>>class KUCMSTimeCalibration : KUCMSRootHelperBaseClass

///////////////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code  --  yes yes this is easier for me, im weird, will divide into hh/cc at end 
///////////////////////////////////////////////////////////////////////////////////////////////////

inline KUCMSTimeCalibration::KUCMSTimeCalibration( bool stayOpen, bool makeNew ){

    std::cout << "Initiating KUCMSTimeCalibrationClass" << std::endl;
	// parameter and setup hardcoded in constructer - some of this needs moved to "run" code

    caliFileDir = "ecal_config/";// move
	detIDConfigEB = caliFileDir + "fullinfo_v2_detids_EB.txt";
    detIDConfigEE = caliFileDir + "fullinfo_v2_detids_EE.txt";
	caliRunConfig = caliFileDir + "caliRunConfig.txt";
    caliTTConfig = caliFileDir + "caliTTConfig.txt";
    caliSmearConfig = caliFileDir + "caliSmearConfig.txt";

    xtalHistMapName = "_X_";
    ttHistMapName = "_TT_";
    lochist = "_SRO_Data_Hist";
    globhist = "_ZEE_Data_Hist";

    //std::cout << " - opening caliHistsTFile " << std::endl;
    caliTFileName = caliFileDir + "caliHistsTFile.root";// name configurable ?
	cali2DResPlotsTFileName = caliFileDir + "cali2dResPlotsTFile.root";
	caliTFile = NULL;
	cali2DResTFile = NULL;

    curTag = "EG_EOY_MINI";
	// PD: EG_EOY_MINI    MC: RunIIFall17DRPremix 
	curTTIov = "r2ul";
    curXIov = "prompt";
	curLumiTag = "r2ul";

	smearTag = "EG300202_DYF17";

	eosDir = "root://cmseos.fnal.gov//store/user/jaking/";
	inDir = "";

	lowEnergy = false;
	useEffEnergy = false;
	xBinStr = "VARIABLE 25 50 75 100 125 150 175 200 225 250 300 400 600 1000 1800"; 
    yBinStr = "CONSTANT 600 -3 3";

    //lowEnergy = true;
	//xBinStr = "VARIABLE 0.2 0.5 1 2 5 10 15 20 25 50";
	//yBinStr = "CONSTANT 1200 -3 3"; 

	maptypes = {"_MeanMap","_ErrMap","_SumMap","_Sum2Map","_OccMap"};

	updated = false;
	externalCali = false;
    useGSwitch = false;

    std::cout << " - setup DetID & IOV Maps " << std::endl;
    SetupDetIDsEB();
    SetupDetIDsEE();

    getRandom = new TRandom();
    getRandom->SetSeed(0);

    std::cout << " - reading in TTrun, Xrun, & smear cali files " << std::endl;
	ReadCaliRunFile();
	ReadTTRunFile();
	ReadSmearFile();

    //std::cout << " - loading cali hists " << std::endl;
	LoadCaliHists( stayOpen, makeNew );

}//<<>>KUCMSTimeCalibration()   

// move save functions external - iopen close TFile in read/save functions
inline KUCMSTimeCalibration::~KUCMSTimeCalibration(){

    std::cout << "Wrapping KUCMSTimeCalibrationClass" << std::endl;

    delete getRandom;

    //std::cout << "Finished Wrapping KUCMSTimeCalibrationClass" << std::endl;

}//<<>>KUCMSTimeCalibration::~KUCMSTimeCalibration()

inline void KUCMSTimeCalibration::SetupDetIDsEB(){

    std::ifstream infile( detIDConfigEB, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    float phi, eta;
    std::string pos;

    while( infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos
                >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM >> phi >> eta ){

        DetIDMap[cmsswId] = {iphi,ieta,TT25,0,phi,eta};
		InvDetIDMap[iphi][ieta][0] = cmsswId;

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )

inline void KUCMSTimeCalibration::SetupDetIDsEE(){

    std::ifstream infile( detIDConfigEE, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    float phi, eta;
    std::string EE;

    while( infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC
                >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant >> phi >> eta ){

        int ec = 1;
        if( side > 0 ) ec = 2;
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec,phi,eta};
        InvDetIDMap[ix][iy][ec] = cmsswId;

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )

inline void KUCMSTimeCalibration::SetupIovMaps(){

	//Prompt 
	std::map<int,int> promptIovMap;
	//promptIovMap[1] = 189158;//03/02/2012
	//promptIovMap[189159] = 203827;//26/03/2012
    //promptIovMap[203828] = 208838;//28/09/2012
    //promptIovMap[208839] = 253983;//11/12/2012

    //promptIovMap[253984] = 276811;//11/08/2015
    //promptIovMap[276812] = 293998;//14/07/2016
    promptIovMap[253984] = 273157;
    promptIovMap[273158] = 284044;

//2016C  	275657	275836
//2016D  	276315	276775
//2016E  	276831	277096
//2016F  	277992	278239 
//2016G  	278820	278969
//2016H 	281613	283830
// iov  253984 284044
// lumi 273158 284044

    //promptIovMap[293999] = 296398;
    promptIovMap[296399] = 297726;
    promptIovMap[297727] = 300201;
    promptIovMap[300202] = 301486;
    promptIovMap[301487] = 304475;
    promptIovMap[304476] = 306460;//05/10/2017
  
//2017B  	297114	297296
//2017C  	299368	301417
//2017D  	302031	302393
//2017E   	303832	304616
//2017F  	305044	305081
// iov  296399 306460
// lumi 297050 306460

    //promptIovMap[307555] = 314766;
    promptIovMap[314767] = 315343;
    promptIovMap[315344] = 316360;
    promptIovMap[316361] = 316569;
    promptIovMap[316570] = 318744;
    promptIovMap[318745] = 319329;
    promptIovMap[319330] = 321506;
    promptIovMap[321507] = 322717;
    promptIovMap[322718] = 323412;
    promptIovMap[323413] = 324305;
    promptIovMap[324306] = 325172;
    promptIovMap[327239] = 356513;//25/11/2018

//2018A		315257	316993
//2018B		317435	317435
//2018C 	319756	320038
//2018D	    320673	322356
// iov  314767 315343
// lumi 315257 325172

//Run 3

    promptIovMap[356514] = 357289;
    promptIovMap[357290] = 358883;
    promptIovMap[358884] = 359420;
    promptIovMap[359421] = 360089;
    promptIovMap[360090] = 360981;
    promptIovMap[360982] = 361416;
    promptIovMap[361417] = 362522;
    promptIovMap[362523] = 367094;//24/11/2022
    promptIovMap[367095] = 367190;//EOY 2022

// 2022

    promptIovMap[367191] = 367515;
    promptIovMap[367516] = 367882;
    promptIovMap[367884] = 368825;
    promptIovMap[368826] = 369801;
    promptIovMap[369802] = 369912;
    promptIovMap[369913] = 370496;
    //promptIovMap[370497] = 372308;// CC period only?
	promptIovMap[370497] = 373577;
	promptIovMap[373578] = 378745;// EOY 2023

/// 2023

    promptIovMap[378746] = 381306;//CC? 
    promptIovMap[381307] = 382007;//ratio timing algorithm 28/05/2024
    promptIovMap[382008] = 382298;
    promptIovMap[382299] = 382721;
	promptIovMap[382722] = 383362;
    promptIovMap[383363] = 383755;
    promptIovMap[383756] = 384187;
    promptIovMap[384188] = 384412;
    promptIovMap[384413] = 384755;
    promptIovMap[384756] = 385152;
    promptIovMap[385153] = 385727;
    promptIovMap[385728] = 386662;
    promptIovMap[386663] = 387741;
    promptIovMap[387742] = 390755;
    promptIovMap[390756] = 999999;// EOY 2024

// 2024

/*
 	// CC calibration periods in 2024
    promptIovMap[372309] = 372750;//"29/08/2023";//cc start
    promptIovMap[372751] = 373539;
    promptIovMap[373540] = 374564;
    promptIovMap[374565] = 375351;
    promptIovMap[375352] = 378673;//"25/10/2023";
    promptIovMap[378674] = 378793;
    promptIovMap[378794] = 379190;
    promptIovMap[379191] = 379323;
    promptIovMap[379324] = 379366;
    promptIovMap[379367] = 379433;
    promptIovMap[379434] = 380065;
    promptIovMap[380066] = 381383;
    promptIovMap[381384] = 382007;
    promptIovMap[382008] = 382208;
    promptIovMap[382209] = 382959;
    promptIovMap[382960] = 385285;
    promptIovMap[385286] = 386951;//"06/09/2024";
    promptIovMap[386952] = 999999;
*/

	iovMaps["prompt"] = promptIovMap;

	std::map<int,int> mcIovMap;
	mcIovMap[0] = 999999;
	iovMaps["mc"] = mcIovMap;

	//float ttlumi = 0.36; // in /fb
	//float xlumi = 9.0; // in /fb
	//float fb_to_ub = 1000000000;
    //float minttlumi =  fb_to_ub * ttlumi; // in /ub
    //float minxlumi = fb_to_ub * xlumi; // in /ub
    float minttlumi =  0.36; // in /fb
    float minxlumi = 9.0; // in /fb

    std::string r2ulTagTT( "r2ultt" );
    std::string r2ulTagX( "r2ulx" );
    std::string r2ulLumiTag( "r2ul" );
    //std::string lumiUL16Config = "UL2016_runlumi.txt";
    //std::string lumiUL17Config = "UL2017_runlumi.txt";
    //std::string lumiUL18Config = "UL2018_runlumi.txt";
    std::string lumiUL16Config = "brilcalc_lumi_2016.txt";
    std::string lumiUL17Config = "brilcalc_lumi_2017.txt";
    std::string lumiUL18Config = "brilcalc_lumi_2018.txt";
    std::string lumi22Config = "brilcalc_lumi_2022.txt";
    std::string lumi23Config = "brilcalc_lumi_2023.txt";
    std::string lumi24Config = "brilcalc_lumi_2024.txt";

    ReadLumiFile( caliFileDir+lumiUL16Config, r2ulLumiTag );
    ReadLumiFile( caliFileDir+lumiUL17Config, r2ulLumiTag );
    ReadLumiFile( caliFileDir+lumiUL18Config, r2ulLumiTag );
	setLumiTag(r2ulLumiTag);
    SetupIovMap( r2ulTagTT, minttlumi );
    SetupIovMap( r2ulTagX, minxlumi );
    setXIov("r2ulx");

    std::string r3TagTT( "r3tt" );
    std::string r3TagX( "r3x" );
    std::string r3LumiTag( "r3lumi" );

    ReadLumiFile( caliFileDir+lumi22Config, r3LumiTag );
    ReadLumiFile( caliFileDir+lumi23Config, r3LumiTag );
    ReadLumiFile( caliFileDir+lumi24Config, r3LumiTag );

    setLumiTag(r3LumiTag);
    SetupIovMap( r3TagTT, minttlumi );
    SetupIovMap( r3TagX, minxlumi );
    //setXIov("r3x");



}//<<>>void KUCMSTimeCalibration::SetupEraIovMap()

inline void KUCMSTimeCalibration::SetupIovMap( std::string tag, float minLumi ){

	bool newrange( true );
	int start = 0;
	int end = 0;
	int prev = 0;
	int oldrun = 0;
    int curxiovend = 0;
	float lumisum = 0;
	int currentrun = 0;
	std::map<int,int> theIovMap;
	auto promptIovMap = iovMaps["prompt"];
	// go through runs by lumi and add up lumi > check for end conditions
	for( auto& lumirun : lumiRunMaps[curLumiTag] ){
		currentrun = lumirun.second.run;
		if( newrange ){ // if we are starting a new range - find curent iov period for this run and set end/start of periods
			prev = start;
			start = ( prev == 0 ) ? currentrun : theIovMap[prev] + 1;
			int iovbound = ( start < oldrun ) ? oldrun : start;
			oldrun = 0;
			newrange = false;
			for( auto& xiov : promptIovMap ){ 
				if( iovbound >= xiov.first && iovbound <= xiov.second ){ 
                    curxiovend = xiov.second; 
					break; 
				}//<<>>if( start >= xiov.first && start <= xiov.second )
			}//<<>>for( auto& xiov : promptIovMap ) 
		}//<<>>if( newrange )
		if( currentrun > curxiovend ){ // tack small portion to previous and start new run period with this run
			if( prev > 0 ){
				end = curxiovend;
				if( end <= start ) end = start;
				theIovMap[prev] = end;
				start = prev;
				oldrun = currentrun;
				lumisum = lumirun.second.lumi;
				newrange = true;
			}//<<>>if( prev > 0 )
			else { 
				for( auto& xiov : promptIovMap ){
                	if( currentrun >= xiov.first && currentrun <= xiov.second ){
                    	curxiovend = xiov.second;
                    	break;
                	}//<<>>if( start >= xiov.first && start <= xiov.second )
            	}//<<>>for( auto& xiov : promptIovMap ) 
            }//<<>>else { if( prev > 0 ){ 
		} else {
			lumisum += lumirun.second.lumi;
			if( lumisum > minLumi ){ 
				end = currentrun;
				if( end <= start ) end = start;
				theIovMap[start] = end;
				oldrun = 0;
				lumisum = 0; 
				newrange = true;
			}//<<>>if( lumisum > minTTLumi )
		}//<<>>if( lumirun.second.run > promptIovMap[curxiov] )
	}//<<>>for( auto& lumirun : lumiRunMap )
	if( not newrange ) theIovMap[prev] = currentrun;
	iovMaps[tag] = theIovMap;		

	for( auto& iov : theIovMap ){ std::cout << "Iov map " << tag << " : " << iov.first << " " << iov.second << std::endl; }

}//<<>>void KUCMSTimeCalibration::makeTTIovMap()

inline void KUCMSTimeCalibration::ReadLumiFile( std::string lumifile, std::string tag ){

	std::cout << " - Read LumiFile : " << lumifile << std::endl;
    std::ifstream infile( lumifile, std::ios::in);
	if( not infile.is_open() ){ std::cout << " -- " << lumifile << " not opened. "  << std::endl; return; }
    std::string infilestr;
	while( std::getline( infile, infilestr ) ){

		if( infilestr[0] != '#' ){

			std::stringstream ss(infilestr);
    		std::string runfill, date, time;
    		int run, fill, nls, ncms;
    		float delivered, recorded;

			ss >> runfill >> date >> time >> nls >> ncms >> delivered >> recorded;
			auto splitRunFill = splitString( runfill, ":" );
			std::stringstream srun( splitRunFill[0] );
			srun >> run;
			std::stringstream sfil( splitRunFill[1] );
			sfil >> fill;		
			//std::cout << " --- " << run << " " << fill  << std::endl;

 			lumiRunMaps[tag][run] = { run, fill, recorded };

		}//<<>>if( infilestr[0] != '#' )

	}//<<>>while(std::getline(iflstream,infilestr))	
	infile.close();

}//<<>>void KUCMSTimeCalibration::upLoadLumiFile( std::string lumifile )

inline void KUCMSTimeCalibration::ReadCaliRunFile(){

	std::cout << " - Reading CaliRunFile " << caliRunConfig << std::endl;
    std::ifstream infile( caliRunConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliRunConfig << " not opened. "  << std::endl; return; }
    std::string caliMapName, tag;
	int srun, erun, lrun;
    float lumi, noise, stoch, stant;
	bool resmap, resparms;
    while( infile >> tag >> caliMapName >> srun >> erun >> lrun >> lumi >> resmap >> resparms >> noise >> stoch >> stant ){
        std::cout << " -- " << tag << " " << caliMapName << " " << srun << " " << erun << " " << lrun << " " << lumi;
		std::cout << " " << resmap << " " << resparms << " " << noise << " " << stoch <<  " " << stant << std::endl;
		CaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lrun, lumi };
        CaliRunMapSet[tag][srun].isNew = false;
		CaliRunMapSet[tag][srun].has2DResMap = resmap;
        CaliRunMapSet[tag][srun].hasResParams = resparms;
        CaliRunMapSet[tag][srun].noise = noise;
        CaliRunMapSet[tag][srun].stoch = stoch;
        CaliRunMapSet[tag][srun].stant = stant;
    }//<<>>while (infile >>
	infile.close();

}//<<>>void ReadTimeCaliTagFile()

inline void KUCMSTimeCalibration::SaveCaliRunFile(){

    std::cout << " - Saving CaliRunFile " << std::endl;
    std::ofstream outfile( caliRunConfig, std::ios::out | std::ios::trunc );
    if( not outfile.is_open() ){ std::cout << " -- " << caliRunConfig << " not opened. "  << std::endl; return; }
	for( auto& calirunmap : CaliRunMapSet ){

		std::string tag = calirunmap.first;
		for( auto& calirunsct : calirunmap.second ){
			if( calirunsct.second.isExternal ) continue;

            std::string mapName = calirunsct.second.histMapName;		
			int srun = calirunsct.second.startRun;
        	int erun = calirunsct.second.endRun;
            int lrun = calirunsct.second.lastRun;
        	float lumi = calirunsct.second.lumi;
			bool resmap = calirunsct.second.has2DResMap;
			bool resparms = calirunsct.second.hasResParams;
			float noise = calirunsct.second.noise;
			float stoch = calirunsct.second.stoch;
			float stant = calirunsct.second.stant;

            //std::cout << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
            outfile << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << " " << resmap; 
			outfile << " " << resparms << " " << noise << " " << stoch << " " << stant << std::endl;

		}//<<>>for( auto& calirunsct : calirunmap[tag] )

	}//<<>>for( auto& calirunmap : CaliRunMapSet )
	outfile.close();

}//<<>>void ReadTimeCaliTagFile()

inline void KUCMSTimeCalibration::ReadTTRunFile(){

    std::cout << " - Reading TTCaliRunFile " << caliTTConfig << std::endl;
    std::ifstream infile( caliTTConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliTTConfig << " not opened. "  << std::endl; return; }
    std::string caliMapName, tag;
    int srun, erun, lrun;
    float lumi;
    while( infile >> tag >> caliMapName >> srun >> erun >> lrun >> lumi ){
        //std::cout << " -- " << tag << " " << caliMapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
        TTCaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lrun, lumi };
        TTCaliRunMapSet[tag][srun].isNew = false;
    }//<<>>while (infile >>
	infile.close();

}//<<>>void KUCMSTimeCalibration::ReadTTRunFile()

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

inline void KUCMSTimeCalibration::SaveTTRunFile(){
    
	std::cout << " - Saving TTCaliRunFile " << std::endl;
    std::ofstream outfile( caliTTConfig, std::ios::out | std::ios::trunc );
    if( not outfile.is_open() ){ std::cout << "Failed to open file !!! " << std::endl; return; }
    for( auto& calirunmap : TTCaliRunMapSet ){
        
        std::string tag = calirunmap.first;
		for( auto& calirunsct : calirunmap.second ){

        	int srun = calirunsct.second.startRun;
        	int erun = calirunsct.second.endRun;
        	std::string mapName = calirunsct.second.histMapName;
        	float lumi = calirunsct.second.lumi;
			int lrun = calirunsct.second.lastRun;        
			//std::cout << " -- " << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
        	outfile << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;

		}//<<>>for( auto& calirunsct : calirunmap )

    }//<<>>for( auto& calirunmap : CaliRunMapSet )
	outfile.close();

}//<<>>void ReadTimeCaliTagFile()

inline void KUCMSTimeCalibration::ReadSmearFile(){

    std::cout << " - Reading SmearFile " << caliSmearConfig << std::endl;
    std::ifstream infile( caliSmearConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliSmearConfig << " not opened. "  << std::endl; return; }
    std::string tag;
    float noise, stoch, stant;
    while( infile >> tag >> noise >> stoch >> stant ){
        std::cout << " -- " << tag << " " << noise << " " << stoch <<  " " << stant << std::endl;
        SmearTagSet[tag] = { noise, stoch, stant };
    }//<<>>while (infile >>
    infile.close();

}//<<>>void KUCMSTimeCalibration::ReadSmearFile()

inline void KUCMSTimeCalibration::SaveSmearFile(){

    std::cout << " - Saving SmearFile " << std::endl;
    std::ofstream outfile( caliSmearConfig, std::ios::out | std::ios::trunc );
    if( not outfile.is_open() ){ std::cout << "Failed to open file !!! " << std::endl; return; }
    for( auto& smeartag : SmearTagSet ){

        std::string tag = smeartag.first;
		float noise = smeartag.second.noise;
        float stoch = smeartag.second.stoch;
        float stant = smeartag.second.stant;

        outfile << tag << " " << noise << " " << stoch << " " << stant << std::endl;

    }//<<>>for( auto& calirunmap : SmearTagSet )
    outfile.close();

}//<<>>void SaveSmearFile()

//std::string ttfilename = ttHistMapName+name+std::to_string(run);
//		MeanMap + pd/camp/tag name + start run  + TT/X ?
//		ErrMap + pd/camp/tag name + start run
//		SumMap + pd/camp/tag name + start run
//		Sum2Map + pd/camp/tag name + start run
//      OccMap + pd/camp/tag name + start run
//     maptypes = {"_MeanMap","_ErrMap","_SumMap","_Sum2Map","_OccMap"};
// std::string tfilename = calirunmap.first + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
//            //calirunsct.second.histMapName = tfilename;
//            std::string sfilename = tfilename + "SumMap";
//            std::string s2filename = tfilename + "Sum2Map";
//            std::string ofilename = tfilename + "OccMap";
//            std::string mfilename = tfilename + "MeanMap";
//            std::string efilename = tfilename + "ErrMap";

inline void KUCMSTimeCalibration::LoadCaliHists( bool stayOpen, bool makeNew ){

    std::cout << " - Loading CaliHists " << std::endl;

    std::cout << " - opening caliHistsTFile " << std::endl;
    //caliTFileName = caliFileDir + "caliHistsTFile.root";// name configurable ?
    if( makeNew ){ caliTFile = TFile::Open( caliTFileName.c_str(), "RECREATE" ); }
    if( stayOpen && not makeNew ){ caliTFile = TFile::Open( caliTFileName.c_str(), "UPDATE" ); }
	if( not stayOpen && not makeNew ){ caliTFile = TFile::Open( caliTFileName.c_str(), "READ" ); }
    caliTFile->cd();

	for( auto& calirunmap : TTCaliRunMapSet ){
		for( auto& calirunsct : calirunmap.second ){
			std::string ttfilename( calirunsct.second.histMapName );
			for( auto mapname : maptypes ){
				std::string filename = ttfilename + mapname;
				//std::cout << " -- loading : " << filename << " ( " << ttfilename << " ) ";
				if( ttfilename != "none" && CaliHists.find(filename) == CaliHists.end() ){
					TH2F* hist = (TH2F*)caliTFile->Get(filename.c_str());
                    //std::cout << " Found " << hist << std::endl;
					if( hist ) CaliHists[filename] = { hist, filename, false, false }; // histfile histname isnew isreshist
				}//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
				//else { std::cout << " Not Found " << std::endl; }
			}//<<>>for( auto mapname : mtype )
		}//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

    for( auto& calirunmap : CaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
        	std::string xtfilename( calirunsct.second.histMapName );
            for( auto mapname : maptypes ){
                std::string filename = xtfilename + mapname;
                //std::cout << " -- loading : " << filename << " ( " << xtfilename << " ) ";
                if( xtfilename != "none" && CaliHists.find(filename) == CaliHists.end() ){
                    TH2F* hist = (TH2F*)caliTFile->Get(filename.c_str());
                    //std::cout << " Found " << hist << std::endl;
                    if( hist ) CaliHists[filename] = { hist, filename, false, false }; // histfile histname isnew isreshist
                }//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
            }//<<>>for( auto mapname : mtype )
/*
			if( calirunsct.second.has2DResMap ){
				//std::cout << " -- Loading : reshists " << " ( " << xtfilename << " ) " << std::endl;	
                std::string lsfname = xtfilename+lochist;
                std::string gbfname = xtfilename+globhist;	
				TH2F* lhist = (TH2F*)caliTFile->Get(lsfname.c_str());
				//std::cout << " --- Found local " << lhist << std::endl;
				if( lhist ) CaliHists[lsfname] = { lhist, lsfname, false, true }; // histfile histname isnew isreshist
				TH2F* ghist = (TH2F*)caliTFile->Get(gbfname.c_str());
                //std::cout << " --- Found global " << ghist << std::endl;
				if( ghist ) CaliHists[gbfname] = { ghist, gbfname, false, true }; // histfile histname isnew isreshist
			}//<<>>if( calirunsct.second.has2DResMap )
*/
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

	makeCaliMaps();

	if( not stayOpen ){

    	for( auto& calihist : CaliHists ){
        	if( calihist.second.h2f != NULL ){ delete calihist.second.h2f; }
        	//if( calihist.second.h1f != NULL ) delete calihist.second.h1f;
    	}//<<>>for( auto& calimapsct : TTCaliMaps )	
		CaliHists.clear();
    	caliTFile->Close();

	}//<<>>if( not stayOpen )

}//<<>>void KUCMSTimeCalibration::SetupCaliHists()

//    TH2F* h2f;
//    std::string histName;
//    bool isNew;

// for updating and/or adding new info to califiles - do last - saving closes TCalifile and clears CaliHists
inline void KUCMSTimeCalibration::SaveCaliHists(){

    std::cout << " - Saving CaliHists " << std::endl;
	caliTFile->cd();
    //caliTFile->Delete("*;1");
	for( auto& calihist : CaliHists ){
		//std::cout << " -- saving : " << calihist.second.h2f << " " << calihist.second.histName << " "; 
        //std::cout << calihist.second.isNew << " " << calihist.second.isResHist << std::endl;
		//std::string command = calihist.second.histName + ";1";
		//bool notNew = not calihist.second.isNew;
        bool notResHist = not calihist.second.isResHist;
		if( calihist.second.h2f != NULL ){
			//std::cout << " -- tfile deleting : " << command << std::endl; 
			//if( notNew && notResHist ) caliTFile->Delete( command.c_str() );
			//std::cout << " -- writing : " << calihist.second.h2f << std::endl;
			if( notResHist ) calihist.second.h2f->Write( calihist.second.h2f->GetName(), TObject::kOverwrite );
            //if( notResHist ) calihist.second.h2f->Write();
			//std::cout << " -- deleting : " << calihist.second.h2f << std::endl;	
		}//<<>>if( calihist.second.h2f != NULL )
        if( calihist.second.h1f != NULL ){
			//if( notNew ) caliTFile->Delete( command.c_str() );
            calihist.second.h1f->Write( calihist.second.h1f->GetName(), TObject::kOverwrite );
            //calihist.second.h1f->Write();
        }//<<>>if( calihist.second.h2f != NULL )
		//std::cout << " -- saved " << std::endl;
	}//<<>>for( auto& calimapsct : TTCaliMaps )
	std::cout << " - Closing caliTFile " << std::endl;
    for( auto& calihist : CaliHists ){
		if( calihist.second.h2f != NULL ){ delete calihist.second.h2f; }
		if( calihist.second.h1f != NULL ){ delete calihist.second.h1f; }
	}//<<>>for( auto& calihist : CaliHists )
	CaliHists.clear();
    caliTFile->Close();

}///<<>>void KUCMSTimeCalibration::SaveCaliMaps()

//std::string ttfilename = ttHistMapName+name+std::to_string(run);
//      MeanMap + pd/camp/tag name + start run  + TT/X ?
//      ErrMap + pd/camp/tag name + start run
//      SumMap + pd/camp/tag name + start run
//      Sum2Map + pd/camp/tag name + start run
//      OccMap + pd/camp/tag name + start run
// maptypes = {"_MeanMap","_ErrMap","_SumMap","_Sum2Map","_OccMap"};
//    std::map<uInt,sumCnt> sumCntMap;
//    std::map<uInt,float> meanMap;
//    std::map<uInt,float> errMap;
// DetIDStruct idinfo = DetIDMap[];


inline void KUCMSTimeCalibration::makeCaliHists(){

	// assumed : all current hists are cleared and we remake everything ( maps in runsets are primary )
	// only process x & TT hists, ? only process updated ( new ) hists ?
	//std::cout << "Deleting Exisiting CaliHists " << std::endl;
	//std::vector<std::string> keys;
	//for( auto& entry : CaliHists ){ 
	//	if( entry.second.isResHist == false ){ 
	//		delete entry.second.h2f;
	//		entry.second.h2f = NULL;
	//		//keys.push_back(entry.first); 
	//	}//<<>>if( not entry.second.isResHist )
	//}//<<>>for( auto& entry : CaliHists )
    //std::cout << "Clearing CaliHists " << std::endl;
	//for( auto& key : keys ){ CaliHists.erase(key); }
	//keys.clear();
	//CaliHists.clear();

    //std::cout << "Recreate TFile " << std::endl;
    caliTFile->cd();
	//caliTFile->Delete("*;1");
    //caliTFile->Close();
    //caliTFile = TFile::Open( caliTFileName.c_str(), "RECREATE" );
    //caliTFile->cd();
	

    std::cout << "Making CaliHists from CaliRunMapSet" << std::endl;
    for( auto& calirunmap : CaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
			if( calirunsct.second.isExternal ) continue;
			std::string tfilename( calirunsct.second.histMapName );
			//std::string tfilename = calirunmap.first + "_" + std::to_string( calirunsct.first ) + xtalHistMapName;
			//calirunsct.second.histMapName = tfilename;
            std::string sfilename = tfilename + maptypes[2];
            std::string s2filename = tfilename + maptypes[3];
            std::string ofilename = tfilename + maptypes[4];
            std::string mfilename = tfilename + maptypes[0];
            std::string efilename = tfilename + maptypes[1];
            std::string mdfilename = tfilename + "_MeanDist";

            for( auto mapname : maptypes ){
                std::string filename = tfilename + mapname;
				//std::cout << " --- creating TH2F for " << filename << std::endl;
				auto calihist = CaliHists.find(filename);
				if( calihist != CaliHists.end() ){
					//std::cout << " --- erase TH2F for " << filename << std::endl;
					if( CaliHists[filename].h2f != NULL ){ 
						//std::cout << " --- delete TH2F for " << filename << std::endl;
						delete CaliHists[filename].h2f;
						CaliHists.erase(filename);
						//CaliHists[filename].h2f->Scale(0);
						//caliTFile->Delete(sfilename.c_str());
					}//<<>>if( CaliHists[filename].h2f != NULL ) 
					//CaliHists.erase(calihist); 
				}//<<>>if( calihist != CaliHists.end() )
				//} else {
				TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),171,-85,86,360,1,361);
				CaliHists[filename] = { hist, filename, true, false }; // histfile histname isnew isopen lastrun
				//}//<<>>if( calihist != CaliHists.end() )
            }//<<>>for( auto mapname : mtype )

			if( CaliHists.find(mdfilename) != CaliHists.end() ){
                if( CaliHists[mdfilename].h1f != NULL ){
					delete CaliHists[mdfilename].h1f;
					CaliHists.erase(mdfilename);
                }//<<>>if( CaliHists[filename].h1f != NULL ) 
				//TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
				//CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
			}//<<>>if( CaliHists.find(mdfilename) == CaliHists.end() )
			TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2); 
			CaliHists[mdfilename] = { hist1d, mdfilename, false, false };

			//std::cout << " -- filling hists for " << tfilename << std::endl;
			for( auto& entry : calirunsct.second.sumCntMap ){
				DetIDStruct idinfo = DetIDMap[entry.first];
				CaliHists[sfilename].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.sum );
                CaliHists[s2filename].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.sumsqr );
                CaliHists[ofilename].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.cnt );
			}//<<>>for( auto& entry : calirunsct.second.sumCntMap )
            for( auto& entry : calirunsct.second.meanMap ){
                DetIDStruct idinfo = DetIDMap[entry.first];
                CaliHists[mfilename].h2f->Fill( idinfo.i2, idinfo.i1, entry.second );
                CaliHists[mdfilename].h1f->Fill( entry.second );
			}//<<>>for( auto& entry : calirunsct.second.meanMap )
            for( auto& entry : calirunsct.second.errMap ){
                DetIDStruct idinfo = DetIDMap[entry.first];
                CaliHists[efilename].h2f->Fill( idinfo.i2, idinfo.i1, entry.second );
            }//<<>>for( auto& entry : calirunsct.second.errMap )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )


    std::cout << "Making CaliHists from TTCaliRunMapSet" << std::endl;
    for( auto& calirunmap : TTCaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
            std::string tfilename( calirunsct.second.histMapName );
            //std::string tfilename = calirunmap.first + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
            //calirunsct.second.histMapName = tfilename;
            std::string sfilename = tfilename + maptypes[2];
            std::string s2filename = tfilename + maptypes[3];
            std::string ofilename = tfilename + maptypes[4];
            std::string mfilename = tfilename + maptypes[0];
            std::string efilename = tfilename + maptypes[1];
            std::string mdfilename = tfilename + "_MeanDist";

            for( auto mapname : maptypes ){
                std::string filename = tfilename + mapname;
                //std::cout << " --- creating TH2F for " << filename << std::endl;
                auto calihist = CaliHists.find(filename);
                if( calihist != CaliHists.end() ){
                    //std::cout << " --- erase TH2F for " << filename << std::endl;
                    if( CaliHists[filename].h2f != NULL ){
                        //std::cout << " --- delete TH2F for " << filename << std::endl;
                        delete CaliHists[filename].h2f;
                        CaliHists.erase(filename);
                        //CaliHists[filename].h2f->Scale(0);
                        //caliTFile->Delete(sfilename.c_str());
                    }//<<>>if( CaliHists[filename].h2f != NULL ) 
                    //CaliHists.erase(calihist); 
                }//<<>>if( calihist != CaliHists.end() )
                //} else {
                TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),35,-17,18,72,1,73);
                CaliHists[filename] = { hist, filename, true, false }; // histfile histname isnew isopen lastrun
                //}//<<>>if( calihist != CaliHists.end() )
            }//<<>>for( auto mapname : mtype )

            if( CaliHists.find(mdfilename) != CaliHists.end() ){
                if( CaliHists[mdfilename].h1f != NULL ){
                    delete CaliHists[mdfilename].h1f;
                    CaliHists.erase(mdfilename);
                }//<<>>if( CaliHists[filename].h1f != NULL ) 
                //TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
                //CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
            }//<<>>if( CaliHists.find(mdfilename) == CaliHists.end() )
            TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
            CaliHists[mdfilename] = { hist1d, mdfilename, false, false };

            //std::cout << " -- filling hists for " << tfilename << std::endl;
            for( auto& entry : calirunsct.second.sumCntMap ){
				auto idinfo = getTTInfo( entry.first );
                CaliHists[sfilename].h2f->Fill( idinfo.second, idinfo.first, entry.second.sum );
                CaliHists[s2filename].h2f->Fill( idinfo.second, idinfo.first, entry.second.sumsqr );
                CaliHists[ofilename].h2f->Fill( idinfo.second, idinfo.first, entry.second.cnt );
            }//<<>>for( auto& entry : calirunsct.second.sumCntMap )
            for( auto& entry : calirunsct.second.meanMap ){
                auto idinfo = getTTInfo( entry.first );
                CaliHists[mfilename].h2f->Fill( idinfo.second, idinfo.first, entry.second );
                CaliHists[mdfilename].h1f->Fill( entry.second );
            }//<<>>for( auto& entry : calirunsct.second.meanMap )
            for( auto& entry : calirunsct.second.errMap ){
                auto idinfo = getTTInfo( entry.first );
                CaliHists[efilename].h2f->Fill( idinfo.second, idinfo.first, entry.second );
            }//<<>>for( auto& entry : calirunsct.second.errMap )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )


}//<<>>void KUCMSTimeCalibration::makeCaliHists()

inline void KUCMSTimeCalibration::makeCaliMaps(){

    std::cout << " - Making CaliMaps " << std::endl;
	//TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),171,-85.5,85.5,360,0.5,360.5);
    for( auto& calirunmap : CaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
			if( calirunsct.second.isExternal ) continue;
            //std::string tfilename = calirunsct.second.histMapName + "_" + std::to_string( calirunsct.first ) + xtalHistMapName;
			std::string tfilename = calirunsct.second.histMapName;
			//std::cout << " -- opening : " << tfilename << std::endl;
            std::string sfilename = tfilename + "_SumMap";
            std::string s2filename = tfilename + "_Sum2Map";
            std::string ofilename = tfilename + "_OccMap";
            std::string mfilename = tfilename + "_MeanMap";
            std::string efilename = tfilename + "_ErrMap";
            //std::cout << " --- with : " << sfilename << std::endl;
			if( CaliHists.find(sfilename) != CaliHists.end() ){
               //std::cout << " --- found : " << tfilename << std::endl;
				for( int ieta = 1; ieta < 172; ieta++ ){
					for( int iphi = 1; iphi < 361; iphi++ ){
						//std::cout << " -- Processing : " << ieta << " " << iphi;
						int i1 = ieta - 86;
						if( i1 == 0 ) continue;
						uInt detid = InvDetIDMap[iphi][i1][ECAL::EB];
						float sum = CaliHists[sfilename].h2f->GetBinContent(ieta,iphi);
						float sum2 = CaliHists[s2filename].h2f->GetBinContent(ieta,iphi);
                		int occ = CaliHists[ofilename].h2f->GetBinContent(ieta,iphi);
						if( occ == 0 ) continue;
						//std::cout << " " << detid << " " << sum << " " << sum2;
						calirunsct.second.sumCntMap[detid] = { sum, sum2, occ };
                        float mean  = CaliHists[mfilename].h2f->GetBinContent(ieta,iphi);
						float error  = CaliHists[efilename].h2f->GetBinContent(ieta,iphi);
                        //std::cout << " " << mean << " " << error << " " << std::endl;
                        calirunsct.second.meanMap[detid] = mean;
                        calirunsct.second.errMap[detid] = error;
					}//<<>>for( int iphi = 1; iphi < 361; iphi++ )
				}//<<>>for( int ieta = 1; ieta < 172; ieta++ )
			}//<<>>if( CaliHists.find(sfilename) != CaliHists.end() )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

    std::cout << " - Making TTCaliMaps " << std::endl;
	//TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),34,0,34,72,0,72);
    for( auto& calirunmap : TTCaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
            //std::string tfilename = calirunsct.second.histMapName + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
            std::string tfilename = calirunsct.second.histMapName;
            //std::cout << " -- opening : " << tfilename << std::endl;
            std::string sfilename = tfilename + "_SumMap";
            std::string s2filename = tfilename + "_Sum2Map";
            std::string ofilename = tfilename + "_OccMap";
            std::string mfilename = tfilename + "_MeanMap";
            std::string efilename = tfilename + "_ErrMap";
            //std::cout << " --- with : " << sfilename << std::endl;
            if( CaliHists.find(sfilename) != CaliHists.end() ){
				//std::cout << " --- found : " << tfilename << std::endl;
                for( int ieta = 1; ieta < 36; ieta++ ){
                    for( int iphi = 1; iphi < 73; iphi++ ){
                        //std::cout << " -- Processing : " << ieta << " " << iphi;
                        int i1 = ieta - 18;
                        if( i1 == 0 ) continue;
						uInt detid = getInvTTId( iphi, i1 );
                        float sum = CaliHists[sfilename].h2f->GetBinContent(ieta,iphi);
                        float sum2 = CaliHists[s2filename].h2f->GetBinContent(ieta,iphi);
                        int occ = CaliHists[ofilename].h2f->GetBinContent(ieta,iphi);
                        if( occ == 0 ) continue;
                        //std::cout << " " << detid << " " << sum << " " << sum2;
                        calirunsct.second.sumCntMap[detid] = { sum, sum2, occ };
                        float mean  = CaliHists[mfilename].h2f->GetBinContent(ieta,iphi);
                        float error  = CaliHists[efilename].h2f->GetBinContent(ieta,iphi);
                        //std::cout << " " << mean << " " << error << " " << std::endl;
                        calirunsct.second.meanMap[detid] = mean;
                        calirunsct.second.errMap[detid] = error;
                    }//<<>>for( int iphi = 1; iphi < 361; iphi++ )
                }//<<>>for( int ieta = 1; ieta < 172; ieta++ )
            }//<<>>if( CaliHists.find(sfilename) != CaliHists.end() )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

}//<<>>void KUCMSTimeCalibration::makeCaliMaps()

void KUCMSTimeCalibration::LoadExtCali( std::string calihist, std::string mapname, std::string tag, int startr, int endr ){

    std::cout << " - Loading External Cali Map : " << mapname << std::endl;
	//std::cout << " --- found : " << tfilename << std::endl;
	
	TFile* histTFile = TFile::Open( calihist.c_str(), "READ" );

	std::string caliMapName = mapname;
    std::string mfilename = "AveXtalRecTimeEBMap";
    std::string ofilename = "AveXtalOccEBMap";
    std::string efilename = "AveXtalRecTimeErrEBMap";
	TH2F* mhist = (TH2F*)histTFile->Get(mfilename.c_str());
    TH2F* ohist = (TH2F*)histTFile->Get(ofilename.c_str());
    TH2F* ehist = (TH2F*)histTFile->Get(efilename.c_str());

    CaliRunMapSet[tag][startr] = { mapname, startr, endr, endr, 0 };
    CaliRunMapSet[tag][startr].isNew = false;
    CaliRunMapSet[tag][startr].isExternal = true;
    CaliRunMapSet[tag][startr].has2DResMap = false;
    CaliRunMapSet[tag][startr].hasResParams = false;
    CaliRunMapSet[tag][startr].noise = 0;
	CaliRunMapSet[tag][startr].stoch = 0;
    CaliRunMapSet[tag][startr].stant = 0;
	
    for( int ieta = 1; ieta < 172; ieta++ ){
    	for( int iphi = 1; iphi < 361; iphi++ ){
        	//std::cout << " -- Processing : " << ieta << " " << iphi;
            int i1 = ieta - 86;
            if( i1 == 0 ) continue;
            uInt detid = InvDetIDMap[iphi][i1][ECAL::EB];
            int occ = ohist->GetBinContent(ieta,iphi);
            if( occ == 0 ) continue;
            //std::cout << " " << detid << " " << sum << " " << sum2;
            CaliRunMapSet[tag][startr].sumCntMap[detid] = { 0, 0, occ };
            float mean  = mhist->GetBinContent(ieta,iphi);
            float error  = ehist->GetBinContent(ieta,iphi);
            //std::cout << " " << mean << " " << error << " " << std::endl;
            CaliRunMapSet[tag][startr].meanMap[detid] = mean;
            CaliRunMapSet[tag][startr].errMap[detid] = error;
        }//<<>>for( int iphi = 1; iphi < 361; iphi++ )
    }//<<>>for( int ieta = 1; ieta < 172; ieta++ )

}//<<>>void KUCMSTimeCalibration::LoadExtCali( std::string calihist, std::string mapname, std::string tag )

inline uInt KUCMSTimeCalibration::getTTId( uInt detId ){

	int ttphi = 1 + int( DetIDMap[detId].i1 - 1 )/5;
	int tteta = 1 + int( std::abs( DetIDMap[detId].i2 ) - 1 )/5;
	uInt id = ( DetIDMap[detId].i2 < 0 ) ? (ttphi+tteta*100)+2000 : (ttphi+tteta*100);
	//std::cout << " getTTId " << detId << " phi " << DetIDMap[detId].i1 << " eta " << DetIDMap[detId].i2 << std::endl; 
    //std::cout << "  --> p " << ttphi << " e " << tteta << " = " << id <<std::endl;
	return id;

}//<<>>uInt KUCMSTimeCalibration::getTTId( uInt detId )

inline std::pair<int,int> KUCMSTimeCalibration::getTTInfo( uInt ttid ){

	int a = ( ttid > 2000 ) ? ttid - 2000 : ttid;
	int t = a/100;
	int i2 = ( ttid > 2000 ) ? -1*t : t;
	int	i1 = a - t*100;
	return std::make_pair(i1,i2); 

}//<<>>std::pair<int,int> KUCMSTimeCalibration::getTTInfo( uInt ttid )

inline uInt KUCMSTimeCalibration::getInvTTId( int i1, int i2 ){

	return ( i2 < 0 ) ? (i1+std::abs(i2)*100)+2000 : (i1+i2*100);

}//<<>>std::pair<int,int> KUCMSTimeCalibration::getTTInfo( uInt ttid )

inline DetIDStruct& KUCMSTimeCalibration::getDetIdInfo( uInt rhid ){

	return DetIDMap[rhid];

}//<<>>DetIDStruct KUCMSTimeCalibration::getDetIdInfo( uInt rhid )

inline uInt KUCMSTimeCalibration::getDetIdInfo( int i1, int i2, int ecal ){

	return InvDetIDMap[i1][i2][ecal];

}//<<>>DetIDStruct KUCMSTimeCalibration::getDetIdInfo( uInt rhid )

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

inline float KUCMSTimeCalibration::getCalibration( uInt rhid, int run, std::string tag ){

	if( rhid == 0 ) return -999.f;
	//if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
	//if( not isEB ){ std::cout << "XCalibration for EE is not supported." << std::endl; return 0.f; }
    //if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
    if( DetIDMap[rhid].ecal != ECAL::EB ) return 999.f;
    float xtaltime = -1000.f;
	float ttcali = ( externalCali ) ? 0.f : getTTCali( rhid, run, tag );
    for( auto& calirunmap : CaliRunMapSet[tag] ){
		int endrun = calirunmap.second.endRun;
		int startrun = calirunmap.second.startRun;	
		//std::cout << "GEtCali : " << run << " " << tag << " " << rhid << " " << startrun << " " << endrun << std::endl;
		if( ( run >= startrun ) && ( run <= endrun ) ){
            xtaltime = calirunmap.second.meanMap[rhid];
			//std::cout << "Found Cali : " << xtaltime << std::endl;
			break;
		}//<<>>if( run >= calirunmap.second.startRun
	}//<<>>for( auto& calirunmap : CaliRunMapSet )
	if( xtaltime == -1000.f ){ std::cout << "XCalibration period not found for run " << run << std::endl; return 0.f; }
    return xtaltime + ttcali;// proper

}//<<>>float KUCMSTimeCalibration::getCalibration( std::string tag )

inline float KUCMSTimeCalibration::getTTCali( uInt rhid, int run, std::string tag ){

    //if( not isEB ){ std::cout << "TTCalibration for EE is not supported." << std::endl; return 0.f; }
    if( DetIDMap[rhid].ecal != ECAL::EB ) return 999.f;
    float xtaltime = -1000.f;
    uInt ttid = getTTId( rhid );
    for( auto& calirunmap : TTCaliRunMapSet[tag] ){
        int endrun = calirunmap.second.endRun;
        int startrun = calirunmap.second.startRun;
        //std::cout << "TTGetCali : " << run << " " << tag << " " << rhid << " " << ttid << " " << startrun << " " << endrun << std::endl;
        if( ( run >= startrun ) && ( run <= endrun ) ){
            xtaltime = calirunmap.second.meanMap[ttid];
            //std::cout << "Found TTCali : " << xtaltime << std::endl;
			break;
        }//<<>>if( run >= calirunmap.second.startRun
    }//<<>>for( auto& calirunmap : CaliRunMapSet )
    if( xtaltime == -1000.f ){ std::cout << "TTCalibration period not found for run " << run << std::endl; return 0.f; }
    return xtaltime;

}//<<>>float KUCMSTimeCalibration::getCalibration( std::string tag )

//struct smearTagStruct {
//
//    float noise;
//    float stoch;
//    float stant;
//
//};//smearParameters
//std::map<std::string,smearTagStruct> SmearTagSet;
//  do mutiple combos 

inline float KUCMSTimeCalibration::getSmearedTime( float rhtime, float rhamp, std::string stag  ){

    double stnoise = SmearTagSet[stag].noise;
    double ststoch = SmearTagSet[stag].stoch;
    double ststant = SmearTagSet[stag].stant;
    //std::cout << " stag smearing : " << stnoise << " " << ststoch << " " << ststant << " " << srun << std::endl;
	double noise = pow(stnoise,2);
	double stoch = ( ststoch > 0 ) ? pow(ststoch,2) : 0;
    double stant = pow(ststant,2);
    double resolution = std::sqrt( ( noise/pow(rhamp,2) + stoch/rhamp  + 2*stant )/2 );
	//std::cout << " getSmearedTime : " << resolution << " " << noise << " " << stoch << " " << stant << std::endl;
    if( resolution <= 0 ){ std::cout << "No smearing values set for this tag : " << stag << std::endl; return rhtime; }
    float smearedtime = getRandom->Gaus( rhtime, resolution );
    return smearedtime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

inline float KUCMSTimeCalibration::getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int crun, std::string ctag, std::string stag ){

    float crhtime = rhtime - getCalibration( rhid, crun, ctag );
    float smrdCalibTime = getSmearedTime( crhtime, rhamp, stag );
    return smrdCalibTime;

}//<<>>float KUCMSTimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

// assume full runs are present in input files
// input file will give the name of a text file "list" of root files, the path to the root files, and the tag to use for those files 
// check if run range started
// check if run range complete 
//     std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>
//	- if started but not complete - is run needed ( lastrun )
//	- if not started - start new cali map

inline void KUCMSTimeCalibration::makeCaliMapsEGR( std::string inputFileName, bool doTT, bool small, bool doCali ){

	std::string whichstring = ( doTT ) ? "for TT " : "for Xtal ";
	std::cout << "Creating calibration files from EgammaRes Ntuples " << whichstring << std::endl; 

	bool debug = false;
    //bool debug = true;

	// make generic fill with specialized calling functions for input and TT v X
	// assume bool doTT to inicate TT of X Cali set

	const std::string treename("tree/llpgtree");

    // Declaration of leaf types
    UInt_t run;
    std::vector<uInt> *rhID = 0;
    std::vector<float> *rhRtTime = 0;
    std::vector<float> *rhEnergy = 0;
	//std::vector<float> *rhAmp = 0;

    // List of branches
    TBranch *b_run;   //!
    TBranch *b_rhID;   //!
    TBranch *b_rhRtTime;   //!
    TBranch *b_rhEnergy;   //!
    //TBranch *b_rhAmp;

	// make sure iov maps exist for cali maps we are making
	if( doTT && iovMaps.find(curTTIov) == iovMaps.end() ){ std::cout << " TT Iov not valid !!" << std::endl; return; }
    if( not doTT && iovMaps.find(curXIov) == iovMaps.end() ){ std::cout << " X Iov not valid !!" << std::endl; return; }

    std::ifstream infilelist(inputFileName);
    std::string infilestr;
	int finalRun = 0;
    while( std::getline( infilelist, infilestr ) ) {

        std::stringstream ss(infilestr);
        std::string infilename, tag;
        int srun, erun;
		
        ss >> infilename >> srun >> erun >> tag;
		if( infilename[0] == '#' ) continue;
		std::string wichtype = doTT ? "TT" : "X";
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing maps for " << tag << " in " << wichtype << std::endl;

		// insure we have TTCalimaps if we are doing Xtal calibration maps
		if( not doTT && TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() ){ 
			std::cout << " No TT maps for this tag !!" << std::endl; 
			return; 
		}//<<>>if( not doTT && TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        std::ifstream infile(infilename);
        std::string instr;
        auto fInTree = new TChain( treename.c_str() );
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,instr)){
            auto tfilename = eosDir + inDir + instr;
            std::cout << "-";
            fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;

       	run = 0;
        rhID = 0;
        rhRtTime = 0;
        rhEnergy = 0;

		// use special calibration rh set?
        fInTree->SetBranchAddress("run", &run, &b_run);
        fInTree->SetBranchAddress("rhCaliID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhCaliRtTime", &rhRtTime, &b_rhRtTime);
        fInTree->SetBranchAddress("rhCaliEnergy", &rhEnergy, &b_rhEnergy);

        auto treeEntries = fInTree->GetEntries();
		auto nEntries = small ? ( treeEntries < 10000000 ) ? treeEntries : 10000000  : treeEntries;
        std::cout << "Starting entry loops for " << nEntries << " Entries" <<  std::endl;
        if( debug ) nEntries = 10000;
		int modbreak = debug ? 1000 : 1000000;
        for (Long64_t centry = 0; centry < nEntries; centry++){

            if( centry%modbreak == 0 or centry == 0){
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " (" << static_cast<float>((10000*centry)/nEntries)/(100) << "%)" << std::endl;
            }//<<>>if( centry%1000000 == 0 or centry == 0)

            auto entry = fInTree->LoadTree(centry);

            b_run->GetEntry(entry);   //!
            b_rhID->GetEntry(entry);   //!
            b_rhRtTime->GetEntry(entry);   //!
            b_rhEnergy->GetEntry(entry);

			if( debug) std::cout << " processing " << run;
            if( debug) std::cout << " in " << srun << " to " << erun;
            if( run < srun || run > erun ) continue;
            if( finalRun < run ) finalRun = run;
			if( debug) std::cout << " eval w/ final : " << finalRun << std::endl;

			int nRecHits = rhID->size();
            if( debug) std::cout << " - looping over nRecHits: " << nRecHits << std::endl;
			bool tagNotFound( true );
            bool runNotFound( true );
			auto& calirunset = doTT ? TTCaliRunMapSet : CaliRunMapSet;
			if( calirunset.find(tag) != calirunset.end() ){ // tag exists
				if( debug) std::cout << " - tag found : " << tag << std::endl;
				tagNotFound = false;
				auto& runset = calirunset[tag];
				if( debug) std::cout << " - In runset : " << tag << " " << runset[run].startRun << std::endl;
				for( auto& runrange : runset ){
					auto& range = runrange.second;
					if( debug) std::cout << " - run chk : " << run << " " << range.startRun << " " << range.endRun << std::endl;
					if( run >= range.startRun && run <= range.endRun ){ // run range exists 
						runNotFound = false;
						if( debug) std::cout << " - check lastRun : " << run << " last " << range.lastRun;
						if( debug) std::cout  << " nrh " << nRecHits << std::endl; 
						// if end == last : all runs in range completed, if run > last : run in range already filled
						if( run > range.lastRun ){
							for( int idx = 0; idx < nRecHits; idx++ ){
								//if( debug) std::cout << rhEnergy->at(idx) << " " << rhID->at(idx) << " " << rhRtTime->at(idx) << std::endl;
                                //if( rhEnergy->at(idx) < 5.0 || rhEnergy->at(idx) > 160 ) continue;
                                if( not lowEnergy && rhEnergy->at(idx) < 5.0 ) continue;
								if( lowEnergy && rhEnergy->at(idx) > 25.0 ) continue;
								uInt id = rhID->at(idx);
								//if( debug) std::cout << " - EB check --- " << id << " / " << DetIDMap[id].ecal;
								//if( debug) std::cout << " - " << ECAL::EB << std::endl;
								if( DetIDMap[id].ecal != ECAL::EB ) continue;
								float time = rhRtTime->at(idx);
           						uInt crsid = doTT ? getTTId( id ) : id;
           						float crstime = ( not doTT && doCali ) ? time - getTTCali( id, run, tag ) : time ;
								//float crstime = time;
                                //if( crstime == 0 ) std::cout << " -- " << crsid << " = " << id;
                                //if( crstime == 0 ) std::cout << " = e " << DetIDMap[id].i2 << " p " << DetIDMap[id].i1 << std::endl;
								range.fillSumCnt( crsid, crstime ); // lastrun valid?????
							}//<<>>for( int idx = 0; idx < nRecHits; idx++ )
						}//<<>>if( run > range.lastRun )
						break;
					}//<<>>if( run > startRun && run <= endRun )
				}//<<>>for( auto& range : runset )
			}//<<>>if( calirunset.find(tag) != calirunset.end() ) 

			if( tagNotFound || runNotFound ){		
				if( debug) std::cout << " - tag not found : " << tag << std::endl;
				std::string tmpxtalmap = "none"; 
				int tstart = 9999999;
				int tend = 0;
				int last = 0;
				float tlumi = 0;
				auto& iovset = doTT ? iovMaps[curTTIov] : iovMaps[curXIov];
				for( auto& iovmap : iovset ){
					if( debug) std::cout << " checking : " << run << " " << iovmap.first << " " << iovmap.second << std::endl;
					//if( run == iovmap.first ){}// first run in range - start new map  endrun = iovmap.second
					if( run >= iovmap.first && run <= iovmap.second ){
						if( debug) std::cout << " IOV found " << std::endl;
						tstart = iovmap.first;						
						tend = iovmap.second;
						break;
					}// not first run in range no current partial mapping ?
				}//<<>>for( auto& iovmap : iovMaps[tag] )
				if( tstart < 9999999 ) { 
					//std::cout << "Found IOV for this run/tag " << run << " - " << tag << std::endl;
					std::string hfname = tag+"_"+std::to_string(tstart)+"_"+wichtype;
					calirunset[tag][tstart] = { hfname, tstart, tend, last, tlumi };
                    for( int idx = 0; idx < nRecHits; idx++ ){
						//if( debug) std::cout << rhEnergy->at(idx) << " " << rhID->at(idx) << " " << rhRtTime->at(idx) << std::endl;
                        if( not lowEnergy && rhEnergy->at(idx) < 5.0 ) continue;
                        if( lowEnergy && rhEnergy->at(idx) > 25.0 ) continue;
                        uInt id = rhID->at(idx);
						if( DetIDMap[id].ecal != ECAL::EB ) continue;
                        float time = rhRtTime->at(idx);
                        uInt crsid = doTT ? getTTId( id ) : id;
                        float crstime = ( not doTT && doCali ) ? time - getTTCali( id, run, tag ) : time;
						calirunset[tag][tstart].fillSumCnt( crsid, crstime );
                    }//<<>>for( int idx = 0; idx < nRecHits; idx++ )
				} // else { std::cout << "No IOV for this run/tag " << run << " - " << tag << " ? !!!!!!! " << std::endl; }
			}//<<>>if( not tagfound )

		}//<<>>for (Long64_t centry = 0; centry < nEntries; cnentry++)

	}//<<>>while (std::getline(infilelist,infilestr))

	std::cout << "Making Mean Maps from input data : doTT? " << doTT << std::endl;
	// fill map hist -----------------------------------------------
	int count = 0;
	auto& calirunset = doTT ? TTCaliRunMapSet : CaliRunMapSet;
	for( auto& calirun : calirunset ){ 
		for( auto& runrange : calirun.second ){
			auto& range = runrange.second;
			if( range.updated ){
				std::cout << "For: " << std::to_string( runrange.first ) << std::endl; 
				if( finalRun > 100000 ) range.lastRun = ( finalRun < range.endRun ) ? finalRun : range.endRun;
				else range.lastRun = 999999;
				range.makeMeanMap( false );
				range.updated = false;
			}//<<>>if( range.updated )
		}//<<>>for( auto& runrange : calirun.second )
	}//<<>>for( auto& calirun : calirunset )

}//<<>>void KUCMSTimeCalibration::makeTTCaliMap( std::string inputFileName )

inline void KUCMSTimeCalibration::plot2dResbyIovForEGR( std::string inputFileName, bool scale, bool usecali, bool smear, std::string ext ){

    std::cout << "Creating 2D Resolution Hist from EgammaRes Ntuples " << std::endl;
    if( not usecali ) std::cout << " -- NoCali " << std::endl;
    if( smear ) std::cout << " -- Smeared " << std::endl;

    //bool debug = false;
    bool debug = true;
    bool small = false;

    const std::string treename("tree/llpgtree");

    bool useAmp(true);

	float lB = 10; // lower and upper limits of energies for rechits used
	float uB = 120; // lower and upper limits of energies for rechits used
	if( lowEnergy ){ lB = 2; }

	std::cout << " -- use low energy : " << lowEnergy << std::endl;
	std::cout << " -- xbins : " << xBinStr << std::endl;
    std::cout << " -- ybins : " << yBinStr << std::endl;
	std::cout << " -- upper energy bound : " << uB << " lower energy bound : " << lB << std::endl;

    // Declaration of leaf types
    UInt_t          run;
    std::vector<unsigned int> *resRhID;
    std::vector<float>   *resAmp;
    std::vector<float>   *resE;
    std::vector<float>   *resRtTime;
    std::vector<float>   *resTOF;

   	std::vector<unsigned int> *rhID;
   	std::vector<bool>    *rhisGS6;
   	std::vector<bool>    *rhisGS1;

    // List of branches
    TBranch        *b_run;   //!
    TBranch        *b_resRhID;   //!
    TBranch        *b_resAmp;   //!
    TBranch        *b_resE;   //!
    TBranch        *b_resRtTime;   //!
    TBranch        *b_resTOF;   //!

    TBranch        *b_rhID;   //!
  	TBranch        *b_rhisGS6;   //!
   	TBranch        *b_rhisGS1;   //!

	std::string nocalistr("_NoCali");
    std::string smearstr("_Smeared");
    std::string lochist("_SRO_Data_Hist");
    std::string druhist("_DRO_Data_Hist");
    std::string globhist("_ZEE_Data_Hist");
    std::string ehist("_Amp_v_E_Hist");
    std::string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)");
    std::string fXTitle("A_{eff}/#sigma_{n} (EBEB)");
    std::string fYTitle("#Delta(Photon Seed Time) [ns] (EBEB)");
    std::string fZTitle("");
	std::string fEYTitle("E_{eff} [GeV] (EBEB)");

    std::vector<float> fXBins;
    std::vector<float> fYBins;
    setBins(xBinStr,fXBins);
    setBins(yBinStr,fYBins);
    const auto xbins = &fXBins[0];
    const auto ybins = &fYBins[0];
    int nMyBins = fXBins.size()-1;

    double goodlev(0);
    double goodlin(0);
    double goodgev(0);
    double goodgin(0);
    double gevents(0);

    std::cout << "open input files list : " << inputFileName << std::endl;

    std::ifstream infilelist(inputFileName);
    std::string infilestr;
    while( std::getline( infilelist, infilestr ) ) {

        std::stringstream ss(infilestr);
        std::string infilename, tag;
        int srun, erun;

        ss >> infilename >> srun >> erun >> tag;
		if( infilename[0] == '#' ) continue;
        std::cout << "open input file : " << infilename << std::endl;
		std::cout << "with tags : " << tag << " and " << smearTag << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing 2D Resolution Map " << std::endl;

        // insure we have Calimaps if we are doing Xtal calibration maps
        if( CaliRunMapSet.find(tag) == CaliRunMapSet.end() ){
            std::cout << " No Cali maps for this tag !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

		if( smear && ( SmearTagSet.find(smearTag) == SmearTagSet.end() ) ){
            std::cout << " No Smear tag found !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        std::ifstream infile(infilename);
        std::string instr;
        auto fInTree = new TChain( treename.c_str() );
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,instr)){
            auto tfilename = eosDir + inDir + instr;
            std::cout << "-";
            fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;

        run = 0;
        resRhID = 0;
        resAmp = 0;
        resE = 0;
        resRtTime = 0;
        resTOF = 0;

   		rhID = 0;
   		rhisGS6 = 0;
   		rhisGS1 = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "resRhID", &resRhID, &b_resRhID );   //!
        fInTree->SetBranchAddress( "resAmp", &resAmp, &b_resAmp);   //!
        fInTree->SetBranchAddress( "resE", &resE, &b_resE);   //!
        fInTree->SetBranchAddress( "resRtTime", &resRtTime, &b_resRtTime);   //!
        fInTree->SetBranchAddress( "resTOF" , &resTOF, &b_resTOF);   //!

        fInTree->SetBranchAddress("rhID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
        fInTree->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);

        std::cout << " Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
		if( small && nEntries > 1000000 ) nEntries = 1000000;
        if( debug ) nEntries = ( nEntries < 1000 ) ? nEntries : 1000;
        if( debug ) std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
		int report = nEntries/20;
        for (auto centry = 0U; centry < nEntries; centry++){

            if( centry%report == 0 or centry == 0){
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " " << (1000.0*static_cast<float>(centry)/static_cast<float>(nEntries))/10.0 << "%" << std::endl;
            }//<<>>if( centry%10000000 == 0 or centry == 0)

            auto entry = fInTree->LoadTree(centry);

            if(debug) std::cout << "Start loop ----------------------------------- " << std::endl;

            gevents++;

            b_run->GetEntry(entry);   //!
            b_resRhID->GetEntry(entry);   //!
            b_resAmp->GetEntry(entry);   //!
            b_resE->GetEntry(entry);   //!
            b_resRtTime->GetEntry(entry);   //!
            b_resTOF->GetEntry(entry);   //!

			b_rhisGS6->GetEntry(entry);   //!
			b_rhisGS1->GetEntry(entry);   //!
			b_rhID->GetEntry(entry);   //!
			int nRecHits = rhID->size();

            if( ( ((*resRhID)[0]) == 0 ) && ( ((*resRhID)[3]) == 0 ) ) continue;

            //int didx = 0;
            for( int didx = debug?0:4; didx < 4; didx++ ){
                if(debug) std::cout << "Run " << run << " id " << (*resRhID)[didx] << " Amp " << (*resAmp)[didx] << " E " << (*resE)[didx];
                if(debug) std::cout << " Rt " << (*resRtTime)[didx]  << " TOF " << (*resTOF)[didx] << std::endl;
            }//<<>>for( int didx = 0; didx < 4; didx++ ){
            //if( srun != 0 && ( run < srun || run > erun ) ) continue;

            if(debug) std::cout << " - Finshed Get Entry " << std::endl;

            bool tagNotFound( true );
            bool runNotFound( true );
            auto& calirunset = CaliRunMapSet;
            if( calirunset.find(tag) != calirunset.end() ){ // tag exists
                if( debug) std::cout << " - tag found : " << tag << std::endl;
                tagNotFound = false;
                auto& runset = calirunset[tag];
                if( debug) std::cout << " - In runset : " << tag << " " << runset[run].startRun << std::endl;
                for( auto& runrange : runset ){
                    auto& range = runrange.second;
                    if( debug) std::cout << " - run chk : " << run << " " << range.startRun << " " << range.endRun << std::endl;
                    if( run >= range.startRun && run <= range.endRun ){ // run range exists 
                        runNotFound = false;
                        if( debug) std::cout << " - check lastRun : " << run << " last " << range.lastRun << std::endl;
                        // if end == last : all runs in range completed, if run > last : run in range already filled
                        if( range.lastRun >= range.endRun ){ // caliRunRange is ready for 2d res plot
							std::string caliend = usecali ? "" : nocalistr;
                            if( smear ) caliend = caliend + smearstr;
							caliend = caliend + ext;
							std::string lsfname = range.histMapName+lochist+caliend;
                            std::string ldfname = range.histMapName+druhist+caliend;
                            std::string gbfname = range.histMapName+globhist+caliend;
							std::string efname = range.histMapName+ehist+caliend;
							if( CaliHists.find(lsfname) == CaliHists.end() ){ // 2D hist already made? 
    							TH2F* theHistLS = new TH2F(lsfname.c_str(),lsfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    							theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
    							theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
    							theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
                				CaliHists[lsfname] = { theHistLS, lsfname, true, true }; // histfile histname isnew isres
    							TH2F* theHistGB = new TH2F(gbfname.c_str(),gbfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
    							theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
    							theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
    							theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[gbfname] = { theHistGB, gbfname, true, true }; // histfile histname isnew isres
                                TH2F* theHistLD = new TH2F(ldfname.c_str(),ldfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[ldfname] = { theHistLD, ldfname, true, true }; // histfile histname isnew isres
								TH2F* theHistAE = new TH2F(efname.c_str(),efname.c_str(),500,0,1000,100,0,200);
                                theHistAE->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistAE->GetYaxis()->SetTitle(fEYTitle.c_str());
                                theHistAE->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[efname] = { theHistAE, efname, true, true }; // histfile histname isnew isres
								range.has2DResMap = true;
							}//if( CaliHists.find(lsfname) == CaliHists.end() ){ 
							if( CaliHists[lsfname].isNew ){

								// ----  fill res hist

                                auto idinfoL0 = DetIDMap[(*resRhID)[0]];
                                auto idinfoL1 = DetIDMap[(*resRhID)[1]];
                                auto idinfoG0 = DetIDMap[(*resRhID)[2]];
                                auto idinfoG1 = DetIDMap[(*resRhID)[3]];
						
								std::vector<bool> isGainId1 = {true,true,true,true};
								if( useGSwitch ){
									for( int resrhit = 0; resrhit < 4; resrhit++ ){
										for( int rhit = 0; rhit < nRecHits; rhit++ ){
											if( (*resRhID)[resrhit] == (*rhID)[rhit] ){
												isGainId1[resrhit] = not ( (*rhisGS1)[rhit] || (*rhisGS6)[rhit] );
												break;
											}//<>if( (*resRhID)[rhit] == (*resRhID)[resrhit] )
										}//<<>>for( int resrhit = 0; resrhit < nResRecHits; resrhi++ )
									}//<<>>for( int rhit = 0; rhit < 4; rhit++ )
								}//<<>>if( useGSwitch )

                                float seedTimeIC00 = ( usecali ) ? getCalibration( (*resRhID)[0], run, tag ) : 0;
                                float seedTimeIC10 = ( usecali ) ? getCalibration( (*resRhID)[1], run, tag ) : 0;
                                float seedTimeIC01 = ( usecali ) ? getCalibration( (*resRhID)[2], run, tag ) : 0;
                                float seedTimeIC11 = ( usecali ) ? getCalibration( (*resRhID)[3], run, tag ) : 0;
								if(debug) std::cout << " -- IC0l: " << seedTimeIC00 << " IC1l: " << seedTimeIC10;
                                if(debug) std::cout << " IC0g: " << seedTimeIC01 << " IC1g: " << seedTimeIC11 << std::endl;

								//-------------------set for local, repo calcs for global --------------------------------

                                if(debug) std::cout << " - Calc 2D Hist" << std::endl;

                                double leffa0 = useAmp ? (*resAmp)[0] : (*resE)[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double leffa1 = useAmp ? (*resAmp)[1] : (*resE)[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double geffa0 = useAmp ? (*resAmp)[2] : (*resE)[2]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double geffa1 = useAmp ? (*resAmp)[3] : (*resE)[3]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double leffe0 = (*resE)[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double leffe1 = (*resE)[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double geffe0 = (*resE)[2]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double geffe1 = (*resE)[3]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double lxfill = (leffa0*leffa1)/sqrt(pow(leffa0,2)+pow(leffa1,2));
                                double lxfille = (leffe0*leffe1)/sqrt(pow(leffe0,2)+pow(leffe1,2));
                                double gxfill = (geffa0*geffa1)/sqrt(pow(geffa0,2)+pow(geffa1,2));
                                double gxfille = (geffe0*geffe1)/sqrt(pow(geffe0,2)+pow(geffe1,2));

                                double ldTOF = (*resTOF)[0]-(*resTOF)[1]; //phoseedTOF_0-phoseedTOF_1;
                                double gdTOF = (*resTOF)[2]-(*resTOF)[3]; //phoseedTOF_0-phoseedTOF_1;
                                double lyf0 = (*resRtTime)[0]-seedTimeIC00;
                                double lyf1 = (*resRtTime)[1]-seedTimeIC10;
                                double gyf0 = (*resRtTime)[2]-seedTimeIC01;
                                double gyf1 = (*resRtTime)[3]-seedTimeIC11;
								if( smear ){
									std::cout << "Times are smeared !!!!!" << std::endl;
									lyf0 = getSmearedTime( lyf0, (*resAmp)[0] ); 
                                    lyf1 = getSmearedTime( lyf1, (*resAmp)[1] );
                                    gyf0 = getSmearedTime( gyf0, (*resAmp)[2] );
                                    gyf1 = getSmearedTime( gyf1, (*resAmp)[3] );
								}//<<>>if( smear )
                                double lyfill = lyf0-lyf1+ldTOF;
                                double gyfill = gyf0-gyf1+gdTOF;

                                bool le_cut = ((*resE)[0]>=lB)&&((*resE)[0]<=uB)&&((*resE)[1]>=lB)&&((*resE)[1]<=uB);
                                bool ge_cut = ((*resE)[2]>=lB)&&((*resE)[2]<=uB)&&((*resE)[3]>=lB)&&((*resE)[3]<=uB);
								if( useGSwitch ){
									le_cut = isGainId1[0] && isGainId1[1] && (*resE)[0] >= lB && (*resE)[1] >= lB;
									ge_cut = isGainId1[2] && isGainId1[3] && (*resE)[2] >= lB && (*resE)[3] >= lB;
                                	std::cout << "GSL: " << isGainId1[0] << " " << isGainId1[1] << " ";
                                	std::cout << (*resE)[0] << " " << (*resE)[1] << " " << le_cut << std::endl;
                                	std::cout << "GSZ: " << isGainId1[2] << " " << isGainId1[3] << " ";
                                	std::cout << (*resE)[2] << " " << (*resE)[3] << " " << ge_cut << std::endl;
								}//<<>>if( useGSwitch )
                                bool leta_cut = (idinfoL0.ecal == ECAL::EB)&&(idinfoL1.ecal == ECAL::EB);
                                bool geta_cut = (idinfoG0.ecal == ECAL::EB)&&(idinfoG1.ecal == ECAL::EB);
                                //bool goodLocTime = (*resRtTime)[0] != 0 && (*resRtTime)[1] != 0;
                                //bool goodGloTime = (*resRtTime)[2] != 0 && (*resRtTime)[3] != 0;
                                bool goodLocTime = (*resRtTime)[0] != 0 && (*resRtTime)[1] != 0 && (*resRtTime)[0] != (*resRtTime)[1];
                                bool goodGloTime = (*resRtTime)[2] != 0 && (*resRtTime)[3] != 0 && (*resRtTime)[2] != (*resRtTime)[3];
                                bool goodLocRHs = (*resRhID)[0] != 0 && (*resRhID)[1] != 0;
                                bool goodGloRHs = (*resRhID)[2] != 0 && (*resRhID)[3] != 0;
								//if( 

                                bool isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
                                bool levent_good = le_cut && leta_cut && goodLocRHs && goodLocTime;
                                bool gevent_good = ge_cut && geta_cut && goodGloRHs && goodGloTime;

                                if(debug) std::cout << " - Fill 2D Hist" << std::endl;
                                if( levent_good && isd_cut ){ 
									CaliHists[lsfname].h2f->Fill(lxfill,lyfill); 
									CaliHists[efname].h2f->Fill(lxfill,lxfille); 
								}//<<>>if( levent_good && isd_cut )
                                if( gevent_good ){ CaliHists[gbfname].h2f->Fill(gxfill,gyfill); }
                                if( levent_good && not isd_cut ){ CaliHists[ldfname].h2f->Fill(lxfill,lyfill); }
                                if(debug) std::cout << " - Fill hists done" << std::endl;

							}//if( CaliHists[lsfname].isNew )
                        }//<<>>if( run > range.lastRun )
                        //continue;//found the correct run set ( iov range ) no need to look further
                    }//<<>>if( run > startRun && run <= endRun )
                }//<<>>for( auto& range : runset )
            }//<<>>if( calirunset.find(tag) != calirunset.end() )
        } // for (auto entry = 0U; entry < nEntries; entry++)

        if(debug) std::cout << " -------- Next Input file " << std::endl;

    } // while (std::getline(infilelist,infiles))

	// --------------  process new histos ----------------------------------------------

    // ----   scale and set isNew to false

    std::cout << " - Scale & Save 2D Hist" << std::endl;
    cali2DResTFile = TFile::Open( cali2DResPlotsTFileName.c_str(), "UPDATE" );
    cali2DResTFile->cd();
	for( auto& hists : CaliHists ){	

		if( hists.second.isResHist && hists.second.isNew ){

    		if( scale ) scaleHist(hists.second.h2f,false,true,false);// hist, scale up?, varible x bins?, varible y bin?
			hists.second.isNew = false;
			hists.second.h2f->Write( hists.second.h2f->GetName(), TObject::kOverwrite );

		}//<<>>if( hists.second.isResHist && hists.second.isNew )

	}//<<>>for( auto& hists : CaliHists )
	cali2DResTFile->Close();

    std::cout << "Finished making 2D delta t v eff amp plots" << std::endl;

}//<<>> void plot2dResolution( std::string indir, std::string infilelistname, 

/*

inline int KUCMSTimeCalibration::getGainId( unsigned int rhID ){



}
*/

inline void KUCMSTimeCalibration::plot2dResolutionEGR( std::string inputFileName, bool scale, bool usecali, bool smear, std::string ext ){

    std::cout << "Creating 2D Resolution Hist from EgammaRes Ntuples " << std::endl;
    if( not usecali ) std::cout << " -- NoCali " << std::endl;
    if( smear ) std::cout << " -- Smeared " << std::endl;

    bool debug = false;
    //bool debug = true;
    bool small = false;
    //bool small = true;
    //bool doEE = true;
	bool doEE = false;

    const std::string treename("tree/llpgtree");

    bool useAmp = not useEffEnergy;

    float lB = 10; // lower and upper limits of energies for rechits used
    float uB = 120; // lower and upper limits of energies for rechits used
    if( lowEnergy ){ lB = 1; uB = 120; }

    std::cout << " -- use low energy : " << lowEnergy << std::endl;
    std::cout << " -- xbins : " << xBinStr << std::endl;
    std::cout << " -- ybins : " << yBinStr << std::endl;
    std::cout << " -- upper energy bound : " << uB << " lower energy bound : " << lB << std::endl;

    // Declaration of leaf types
    UInt_t          run;
    std::vector<unsigned int> *resRhID;
    std::vector<float>   *resAmp;
    std::vector<float>   *resE;
    std::vector<float>   *resRtTime;
    std::vector<float>   *resTOF;

    std::vector<unsigned int> *rhID;
    std::vector<bool>    *rhisGS6;
    std::vector<bool>    *rhisGS1;

    // List of branches
    TBranch        *b_run;   //!
    TBranch        *b_resRhID;   //!
    TBranch        *b_resAmp;   //!
    TBranch        *b_resE;   //!
    TBranch        *b_resRtTime;   //!
    TBranch        *b_resTOF;   //!

    TBranch        *b_rhID;   //!
    TBranch        *b_rhisGS6;   //!
    TBranch        *b_rhisGS1;   //!

    std::string nocalistr("_NoCali");
    std::string smearstr("_Smeared");
    std::string lochist("_SRO_Data_Hist");
    std::string druhist("_DRO_Data_Hist");
    std::string globhist("_ZEE_Data_Hist");
    std::string ehist("_Amp_v_E_Hist");
    std::string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)");
    std::string fXTitle("A_{eff}/#sigma_{n} (EBEB)");
    std::string fYTitle("#Delta(Photon Seed Time) [ns] (EBEB)");
    std::string fZTitle("");
    std::string fEYTitle("E_{eff} [GeV] (EBEB)");

    std::vector<float> fXBins;
    std::vector<float> fYBins;
    setBins(xBinStr,fXBins);
    setBins(yBinStr,fYBins);
    const auto xbins = &fXBins[0];
    const auto ybins = &fYBins[0];
    int nMyBins = fXBins.size()-1;

    double goodlev(0);
    double goodlin(0);
    double goodgev(0);
    double goodgin(0);
    double gevents(0);

    std::cout << "open input files list : " << inputFileName << std::endl;

    std::ifstream infilelist(inputFileName);
    std::string infilestr;
    while( std::getline( infilelist, infilestr ) ) {

        std::stringstream ss(infilestr);
        std::string infilename, tag;
        int srun, erun;

        ss >> infilename >> srun >> erun >> tag;
        if( infilename[0] == '#' ) continue;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "with tags : " << tag << " and " << smearTag << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing 2D Resolution Map " << std::endl;

        // insure we have Calimaps if we are doing Xtal calibration maps
        if( usecali && CaliRunMapSet.find(tag) == CaliRunMapSet.end() ){
            std::cout << " No Cali maps for this tag !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        if( smear && ( SmearTagSet.find(smearTag) == SmearTagSet.end() ) ){
            std::cout << " No Smear tag found !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        std::ifstream infile(infilename);
        std::string instr;
        auto fInTree = new TChain( treename.c_str() );
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,instr)){
            auto tfilename = eosDir + inDir + instr;
            std::cout << "-";
            fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;

        run = 0;
        resRhID = 0;
        resAmp = 0;
        resE = 0;
        resRtTime = 0;
        resTOF = 0;

        rhID = 0;
        rhisGS6 = 0;
        rhisGS1 = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "resRhID", &resRhID, &b_resRhID );   //!
        fInTree->SetBranchAddress( "resAmp", &resAmp, &b_resAmp);   //!
        fInTree->SetBranchAddress( "resE", &resE, &b_resE);   //!
        fInTree->SetBranchAddress( "resRtTime", &resRtTime, &b_resRtTime);   //!
        fInTree->SetBranchAddress( "resTOF" , &resTOF, &b_resTOF);   //!

        fInTree->SetBranchAddress("rhID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
        fInTree->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);

        std::cout << " Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
        if( small && nEntries > 1000000 ) nEntries = 1000000;
        if( debug ) nEntries = ( nEntries < 10000 ) ? nEntries : 10000;
        if( debug ) std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
        int report = nEntries/20;
        for (auto centry = 0U; centry < nEntries; centry++){

            if( centry%report == 0 or centry == 0){
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " " << (1000.0*static_cast<float>(centry)/static_cast<float>(nEntries))/10.0 << "%" << std::endl;
            }//<<>>if( centry%10000000 == 0 or centry == 0)

            auto entry = fInTree->LoadTree(centry);

            if(debug) std::cout << "Start loop ----------------------------------- " << std::endl;

            gevents++;

            b_run->GetEntry(entry);   //!
            b_resRhID->GetEntry(entry);   //!
            b_resAmp->GetEntry(entry);   //!
            b_resE->GetEntry(entry);   //!
            b_resRtTime->GetEntry(entry);   //!
            b_resTOF->GetEntry(entry);   //!

            b_rhisGS6->GetEntry(entry);   //!
            b_rhisGS1->GetEntry(entry);   //!
            b_rhID->GetEntry(entry);   //!
            int nRecHits = rhID->size();

            if( ( ((*resRhID)[0]) == 0 ) && ( ((*resRhID)[3]) == 0 ) ) continue;

            if(debug) std::cout << " - Finshed Get Entry " << std::endl;

//----------------------------------------------------------------------------------
            //bool tagNotFound( true );
            //bool runNotFound( true );
            //auto& calirunset = CaliRunMapSet;
            //if( calirunset.find(tag) != calirunset.end() ){ // tag exists
                //if( debug) std::cout << " - tag found : " << tag << std::endl;
                //tagNotFound = false;
                //auto& runset = calirunset[tag];
                //if( debug) std::cout << " - In runset : " << tag << std::endl;
                //for( auto& runrange : runset ){
                    //auto& range = runrange.second;
                    //if( debug) std::cout << " - run chk : " << run << " " << range.startRun << " " << range.endRun << std::endl;
                    //if( run >= range.startRun && run <= range.endRun ){ // run range exists 
                    if( run >= srun && run <= erun ){ // run range exists 
                        //runNotFound = false;
                        //if( debug) std::cout << " - check lastRun : " << run << " last " << range.lastRun << std::endl;
                        // if end == last : all runs in range completed, if run > last : run in range already filled
                        //if( range.lastRun >= range.endRun ){ // caliRunRange is ready for 2d res plot
                            std::string caliend = usecali ? "" : nocalistr;
                            if( smear ) caliend = caliend + smearstr;
                            caliend = caliend + ext;
							std::string histMapName = "ResMap_" + std::to_string(srun) + "_" + std::to_string(erun) + "_";
                            std::string lsfname = histMapName+lochist+caliend;
                            std::string ldfname = histMapName+druhist+caliend;
                            std::string gbfname = histMapName+globhist+caliend;
                            std::string efname = histMapName+ehist+caliend;
                            if( CaliHists.find(lsfname) == CaliHists.end() ){ // 2D hist already made? 
                                TH2F* theHistLS = new TH2F(lsfname.c_str(),lsfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[lsfname] = { theHistLS, lsfname, true, true }; // histfile histname isnew isres
                                TH2F* theHistGB = new TH2F(gbfname.c_str(),gbfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[gbfname] = { theHistGB, gbfname, true, true }; // histfile histname isnew isres
                                TH2F* theHistLD = new TH2F(ldfname.c_str(),ldfname.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[ldfname] = { theHistLD, ldfname, true, true }; // histfile histname isnew isres
                                TH2F* theHistAE = new TH2F(efname.c_str(),efname.c_str(),500,0,1000,100,0,200);
                                theHistAE->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistAE->GetYaxis()->SetTitle(fEYTitle.c_str());
                                theHistAE->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[efname] = { theHistAE, efname, true, true }; // histfile histname isnew isres
                                //range.has2DResMap = true;
                            }//if( CaliHists.find(lsfname) == CaliHists.end() ){ 
                            if( CaliHists[lsfname].isNew ){

                                // ----  fill res hist

                                auto idinfoL0 = DetIDMap[(*resRhID)[0]];
                                auto idinfoL1 = DetIDMap[(*resRhID)[1]];
                                auto idinfoG0 = DetIDMap[(*resRhID)[2]];
                                auto idinfoG1 = DetIDMap[(*resRhID)[3]];

                                std::vector<bool> isGainId1 = {true,true,true,true};
                                if( useGSwitch ){
                                    for( int resrhit = 0; resrhit < 4; resrhit++ ){
                                        for( int rhit = 0; rhit < nRecHits; rhit++ ){
                                            if( (*resRhID)[resrhit] == (*rhID)[rhit] ){
                                                isGainId1[resrhit] = not ( (*rhisGS1)[rhit] || (*rhisGS6)[rhit] );
                                                break;
                                            }//<>if( (*resRhID)[rhit] == (*resRhID)[resrhit] )
                                        }//<<>>for( int resrhit = 0; resrhit < nResRecHits; resrhi++ )
                                    }//<<>>for( int rhit = 0; rhit < 4; rhit++ )
                                }//<<>>if( useGSwitch )

                                float seedTimeIC00 = ( usecali ) ? getCalibration( (*resRhID)[0], run, tag ) : 0;
                                float seedTimeIC10 = ( usecali ) ? getCalibration( (*resRhID)[1], run, tag ) : 0;
                                float seedTimeIC01 = ( usecali ) ? getCalibration( (*resRhID)[2], run, tag ) : 0;
                                float seedTimeIC11 = ( usecali ) ? getCalibration( (*resRhID)[3], run, tag ) : 0;
                                //if(debug) std::cout << " -- IC0l: " << seedTimeIC00 << " IC1l: " << seedTimeIC10;
                                //if(debug) std::cout << " IC0g: " << seedTimeIC01 << " IC1g: " << seedTimeIC11 << std::endl;

                                //-------------------set for local, repo calcs for global --------------------------------

                                if(debug) std::cout << " - Calc 2D Hist" << std::endl;

                                double leffa0 = useAmp ? (*resAmp)[0] : (*resE)[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double leffa1 = useAmp ? (*resAmp)[1] : (*resE)[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double geffa0 = useAmp ? (*resAmp)[2] : (*resE)[2]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double geffa1 = useAmp ? (*resAmp)[3] : (*resE)[3]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double leffe0 = (*resE)[0]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double leffe1 = (*resE)[1]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double geffe0 = (*resE)[2]; //(phoseedE_0/phoseedadcToGeV_0)/phoseedpedrms12_0;
                                double geffe1 = (*resE)[3]; //(phoseedE_1/phoseedadcToGeV_1)/phoseedpedrms12_1;
                                double lxfill = (leffa0*leffa1)/sqrt(pow(leffa0,2)+pow(leffa1,2));
                                double lxfille = (leffe0*leffe1)/sqrt(pow(leffe0,2)+pow(leffe1,2));
                                double gxfill = (geffa0*geffa1)/sqrt(pow(geffa0,2)+pow(geffa1,2));
                                double gxfille = (geffe0*geffe1)/sqrt(pow(geffe0,2)+pow(geffe1,2));

                                double ldTOF = (*resTOF)[0]-(*resTOF)[1]; //phoseedTOF_0-phoseedTOF_1;
                                double gdTOF = (*resTOF)[2]-(*resTOF)[3]; //phoseedTOF_0-phoseedTOF_1;
                                double lyf0 = (*resRtTime)[0]-seedTimeIC00;
                                double lyf1 = (*resRtTime)[1]-seedTimeIC10;
                                double gyf0 = (*resRtTime)[2]-seedTimeIC01;
                                double gyf1 = (*resRtTime)[3]-seedTimeIC11;

                                if( smear ){
                                    std::cout << "Times are smeared !!!!!" << std::endl;

                                    lyf0 = getSmearedTime( lyf0, (*resAmp)[0] );
                                    lyf1 = getSmearedTime( lyf1, (*resAmp)[1] );
                                    gyf0 = getSmearedTime( gyf0, (*resAmp)[2] );
                                    gyf1 = getSmearedTime( gyf1, (*resAmp)[3] );

                                }//<<>>if( smear )
                                double lyfill = lyf0-lyf1+ldTOF;
                                double gyfill = gyf0-gyf1+gdTOF;

                                bool le_cut = ((*resE)[0]>=lB)&&((*resE)[0]<=uB)&&((*resE)[1]>=lB)&&((*resE)[1]<=uB);
                                bool ge_cut = ((*resE)[2]>=lB)&&((*resE)[2]<=uB)&&((*resE)[3]>=lB)&&((*resE)[3]<=uB);
                                if( useGSwitch ){
                                    le_cut = isGainId1[0] && isGainId1[1] && (*resE)[0] >= lB && (*resE)[1] >= lB;
                                    ge_cut = isGainId1[2] && isGainId1[3] && (*resE)[2] >= lB && (*resE)[3] >= lB;
                                }//<<>>if( useGSwitch )
                                bool leta_cut = (idinfoL0.ecal == ECAL::EB)&&(idinfoL1.ecal == ECAL::EB);
                                bool geta_cut = (idinfoG0.ecal == ECAL::EB)&&(idinfoG1.ecal == ECAL::EB);
								if( doEE ){
                                	leta_cut = (idinfoL0.ecal != ECAL::EB)&&(idinfoL1.ecal != ECAL::EB);
                                	geta_cut = (idinfoG0.ecal != ECAL::EB)&&(idinfoG1.ecal != ECAL::EB);
								}//<<>>if( doEE )
                                bool goodLocTime = (*resRtTime)[0] != 0 && (*resRtTime)[1] != 0 && (*resRtTime)[0] != (*resRtTime)[1];
                                bool goodGloTime = (*resRtTime)[2] != 0 && (*resRtTime)[3] != 0 && (*resRtTime)[2] != (*resRtTime)[3];
                                bool goodLocRHs = (*resRhID)[0] != 0 && (*resRhID)[1] != 0;
                                bool goodGloRHs = (*resRhID)[2] != 0 && (*resRhID)[3] != 0;

                                bool isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
                                bool levent_good = le_cut && leta_cut && goodLocRHs && goodLocTime;
                                bool gevent_good = ge_cut && geta_cut && goodGloRHs && goodGloTime;

                                if(debug) std::cout << " - Fill 2D Hist" << std::endl;
								
								if( debug && ( levent_good || gevent_good ) ){ for( int didx = 0; didx < 4; didx++ ){
                                //if( ( levent_good || gevent_good ) ){ for( int didx = 0; didx < 4; didx++ ){
                					std::cout << "Run " << run << " id " << (*resRhID)[didx]; 
                                    std::cout << " Amp " << (*resAmp)[didx] << " E " << (*resE)[didx];
                					std::cout << " Rt " << (*resRtTime)[didx]  << " TOF " << (*resTOF)[didx] << std::endl;
                                    std::cout << "GSL: " << isGainId1[0] << " " << isGainId1[1] << " ";
                                    std::cout << (*resE)[0] << " " << (*resE)[1] << " " << le_cut << std::endl;
                                    std::cout << "GSZ: " << isGainId1[2] << " " << isGainId1[3] << " ";
                                    std::cout << (*resE)[2] << " " << (*resE)[3] << " " << ge_cut << std::endl;
            					}}//<<>>for( int didx = 0; didx < 4; didx++ ){

								if( levent_good && debug ){ 
									std::cout << " -- IC0l: " << seedTimeIC00 << " IC1l: " << seedTimeIC10 << std::endl;
								}//<<>>if( levent_good && debug )
                                if( gevent_good && debug ){ 
									std::cout << " -- IC0g: " << seedTimeIC01 << " IC1g: " << seedTimeIC11 << std::endl;
								}//<<>>if( gevent_good && debug )
                                if( levent_good && isd_cut ){
                                    CaliHists[lsfname].h2f->Fill(lxfill,lyfill);
                                    CaliHists[efname].h2f->Fill(lxfill,lxfille);
                                }//<<>>if( levent_good && isd_cut )
                                if( gevent_good ){ CaliHists[gbfname].h2f->Fill(gxfill,gyfill); }
                                if( levent_good && not isd_cut ){ CaliHists[ldfname].h2f->Fill(lxfill,lyfill); }
                                if(debug) std::cout << " - Fill hists done" << std::endl;

                            }//if( CaliHists[lsfname].isNew )
                        //}//<<>>if( run > range.lastRun )
                        //continue;//found the correct run set ( iov range ) no need to look further
                    }//<<>>if( run > srun && run <= erun )
                //}//<<>>for( auto& range : runset )
            //}//<<>>if( calirunset.find(tag) != calirunset.end() )
        } // for (auto entry = 0U; entry < nEntries; entry++)

        if(debug) std::cout << " -------- Next Input file " << std::endl;

    } // while (std::getline(infilelist,infiles))

    // --------------  process new histos ----------------------------------------------

    // ----   scale and set isNew to false

    std::cout << " - Scale & Save 2D Hist" << std::endl;
    cali2DResTFile = TFile::Open( cali2DResPlotsTFileName.c_str(), "UPDATE" );
    cali2DResTFile->cd();
    for( auto& hists : CaliHists ){

        if( hists.second.isResHist && hists.second.isNew ){

            if( scale ) scaleHist(hists.second.h2f,false,true,false);// hist, scale up?, varible x bins?, varible y bin?
            hists.second.isNew = false;
            if( not small ) hists.second.h2f->Write( hists.second.h2f->GetName(), TObject::kOverwrite );

        }//<<>>if( hists.second.isResHist && hists.second.isNew )

    }//<<>>for( auto& hists : CaliHists )
    cali2DResTFile->Close();

    std::cout << "Finished making 2D delta t v eff amp plots" << std::endl;

}//<<>> void plot2dResolution( std::string indir, std::string infilelistname, 

inline SigmaFitResult KUCMSTimeCalibration::runTimeFitter( TH2F* hist2D ){

  	std::string f2DHistName = hist2D->GetName();
	std::cout << "Running time fitter : " << f2DHistName << std::endl;
    //bool doSterm = true;
	bool doSterm = false;
	
  	std::vector<float> fXBins;
	//std::cout << " -- " << xBinStr << std::endl;
 	setBins( xBinStr, fXBins );
  	int fNBinsX = fXBins.size();
	const auto xbins = &fXBins[0];

	std::string outfilename = caliFileDir + f2DHistName + "_resfit.root";
    TFile* resTFile = TFile::Open( outfilename.c_str(), "UPDATE" );
    resTFile->cd();

	int nXbins = fNBinsX - 1;
    std::map<std::string,TH1F*> ResultsMap;
    std::string fTitle = "#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)";
    std::string fXTitle = "A_{eff}/#sigma_{n} (EBEB)";
    std::string fTimeText = "t_{1}-t_{2}";
    std::string fSigmaVarText = "A_{eff}/#sigma_{n}";
    std::string chi2ndfName = f2DHistName+"_chi2ndf";
    std::string chi2ndfTitle = fTitle+" "+"#chi^{2}/NDF."+";"+fXTitle+";"+"#chi^{2}/NDF.";
    ResultsMap["chi2ndf"]  = new TH1F(chi2ndfName.c_str(),chi2ndfName.c_str(),nXbins,xbins);
    ResultsMap["chi2ndf"]->Sumw2();
    std::string chi2probName = f2DHistName+"_chi2prob";
    std::string chi2probTitle = fTitle+" "+"#chi^{2} Prob."+";"+fXTitle+";"+"#chi^{2} Prob.";
    ResultsMap["chi2prob"] = new TH1F(chi2probName.c_str(),chi2probTitle.c_str(),nXbins,xbins);
    ResultsMap["chi2prob"]->Sumw2();
    std::string mutext = "#mu("+fTimeText+") [ns]";
    std::string muName = f2DHistName+"_mu";
    std::string muTitle = fTitle+" "+mutext+";"+fXTitle+";"+mutext;
    ResultsMap["mu"] = new TH1F(muName.c_str(),muTitle.c_str(),nXbins,xbins);
    ResultsMap["mu"]->Sumw2();
    std::string sigtext = "#sigma("+fTimeText+") [ns]";
    std::string sigName = f2DHistName+"_sigma";
    std::string sigTitle = fTitle+" "+sigtext+";"+fXTitle+";"+sigtext;
    ResultsMap["sigma"] = new TH1F(sigName.c_str(),sigTitle.c_str(),nXbins,xbins);
    ResultsMap["sigma"]->Sumw2();
    std::string occtext = "#occ("+fTimeText+")";
    std::string occName = f2DHistName+"_occ";
    std::string occTitle = fTitle+" "+occtext+";"+fXTitle+";"+occtext;
    ResultsMap["occ"] = new TH1F(occName.c_str(),occTitle.c_str(),nXbins,xbins);
    ResultsMap["occ"]->Sumw2();
    std::string rmstext = "#rms("+fTimeText+")";
    std::string rmsName = f2DHistName+"_rms";
    std::string rmsTitle = fTitle+" "+rmstext+";"+fXTitle+";"+rmstext;
    ResultsMap["rms"] = new TH1F(rmsName.c_str(),rmsTitle.c_str(),nXbins,xbins);
    ResultsMap["rms"]->Sumw2();

  	std::map<int,ProfileTimeFit> profileHists;
	//std::vector<ProfileTimeFit> profileHists;
  	// Project out 2D hists into map
  	//-------------TimeFitter::Project2Dto1DHists(FitInfo);
  	std::cout << "Projecting to 1D from 2D plot: " << f2DHistName << std::endl;
  	// get inputs/outputs
  	for (int ibinX = 1; ibinX < fNBinsX; ibinX++){
    	std::string fullHistName = f2DHistName + "_ibin_" + std::to_string(ibinX);
		std::cout << " - projecting bin : " << ibinX << std::endl;
		TH1F* profilehist = (TH1F*)hist2D->ProjectionY(fullHistName.c_str(),ibinX,ibinX);
        profilehist->SetTitle(fullHistName.c_str());
		std::cout << " -- saving to profileHists bin : " << ibinX << std::endl;
		ProfileTimeFit thisprofile( profilehist );
		//auto intgral = thisprofile.profileHist->Integral();
    	profileHists[ibinX] = thisprofile;
		//profileHists.push_back(thisprofile);
        auto intgral = profileHists[ibinX].profileHist->Integral();
    	auto name = profileHists[ibinX].profileHist->GetName();
        std::cout << " -- saved " << name << std::endl;
  	}//<<>>for (int ibinX = 1; ibinX <= fNBinsX; ibinX++)

  	// Fit each 1D hist
  	//--------------TimeFitter::Fit1DHists(FitInfo);
  	std::cout << "Fitting profile hists for: " << f2DHistName << std::endl;
  	// get inputs/outputs
  	//auto mapSize = profileHists.size();
	//auto testhist = profileHists[1].profileHist;
  	//auto intgral = profileHists[1].profileHist->Integral();
	//auto name = profileHists[2].profileHist->GetName();
	//std::cout << " -- Profile Map size : " << name << std::endl;
  	for (int ibinX = 1; ibinX < fNBinsX; ibinX++){
		profileHists[ibinX].DoFit();	
	}//<<>>for (int ibinX = 1; ibinX < fNBinsX; ibinX++)
  	//for( auto& profile : profileHists ){ profile.second.DoFit(); }
    std::cout << "Saving profile hists for: " << f2DHistName << std::endl;

  	// Extract mu and sigma into maps
  	//---------------TimeFitter::ExtractFitResults(FitInfo);
  	std::cout << "Extracting results for: " << f2DHistName << std::endl;
  	// set bin content!
  	for (auto ibinX = 1; ibinX < fNBinsX; ibinX++){

		auto& profile = profileHists[ibinX];
    	//if( profile.isEmpty() ) continue;
    	TimeFitResult results = profile.GetFitResult();

		float herr = std::abs( results.sigmaHigh - results.sigma );
        float lerr = std::abs( results.sigmaLow - results.sigma ); 
		float serr = ( herr > lerr ) ? herr : lerr;
		//float serr = results.esigma;		

    	std::cout << " Bin " << ibinX  <<  " : " << fXBins[ibinX-1] << "-" << fXBins[ibinX];
    	std::cout << " : no guass: " << results.sigma << " err: " << serr << std::endl;
	
		if( results.sigma <= 0 ) continue; 

    	ResultsMap["chi2ndf"]->SetBinContent( ibinX, results.chi2ndf );
    	ResultsMap["chi2prob"]->SetBinContent( ibinX, results.chi2prob );
    	ResultsMap["mu"]->SetBinContent( ibinX, results.mu );
    	ResultsMap["mu"]->SetBinError( ibinX, results.emu );
		if( results.esigma < results.sigma/2 ){
    		ResultsMap["sigma"]->SetBinContent( ibinX, results.sigma );
            ResultsMap["sigma"]->SetBinError( ibinX, serr );
		}//<<>>if( results.esigma < 0.01 )
    	ResultsMap["occ"]->SetBinContent( ibinX, results.occ );
    	ResultsMap["occ"]->SetBinError( ibinX, std::sqrt(results.occ) );
    	ResultsMap["rms"]->SetBinContent( ibinX, results.rms );
    	ResultsMap["rms"]->SetBinError( ibinX, 0 );
 
  	}//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

    for( auto& profile : profileHists ){ 
		profile.second.profileHist->Write( profile.second.profileHist->GetName(), TObject::kOverwrite ); 
		profile.second.deleteHists();
	}//<<>>for( auto& profile : profileHists )

  	// Prep sigma fit
  	//----------------------TimeFitter::PrepSigmaFit(FitInfo);
  	//std::cout << "Prepping sigma fit for: " << f2DHistName << std::endl;
    float nLower = 1; 
	float nVal = 50; 
	float nUpper = 100;
    float cLower = 0.01; 
	float cVal = 1; 
	float cUpper = 5;
  	// get input hist
  	auto& hist = ResultsMap["sigma"];
  	//auto x_low = hist->GetXaxis()->GetBinLowEdge( hist->GetXaxis()->GetFirst() );
  	//auto x_up  = hist->GetXaxis()->GetBinUpEdge( hist->GetXaxis()->GetLast() );
	float x_low = 75;
	float x_up = 1800;
	if( lowEnergy ){ x_low = 5; x_up = 1800; }
    if( lowEnergy && useEffEnergy ){ x_low = 1; x_up = 120; }

  	std::string histname = hist->GetName();
  	std::string formname = histname+"_form";
  	std::string fitname  = histname+"_fit";
    //std::cout << "Prepping sigma fit for: " << histname << std::endl;
	std::string fitformstr;
	if( doSterm ) fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))"; 
	else fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))"; 
  	auto form = new TFormula(formname.c_str(),fitformstr.c_str());
  	auto fit = new TF1(fitname.c_str(),form->GetName(),x_low,x_up);
  	fit->SetParName( 0, "N" ); 
	fit->SetParameter( 0, nVal ); 
	fit->SetParLimits( 0, nLower, nUpper );
  	fit->SetParName( 1, "C" ); 
	fit->SetParameter( 1, cVal ); 
	fit->SetParLimits( 1, cLower, cUpper );
  	if( doSterm ) fit->SetParName( 2, "S" ); 
  	if( doSterm ) fit->SetParameter( 2, 1.0 ); 
  	if( doSterm ) fit->SetParLimits( 2, 0.0, 10.0 );
  	fit->SetLineColor(hist->GetLineColor());
  	hist->Fit(fit->GetName(),"RBQO");

	SigmaFitResult sigmafit;
    sigmafit.noise = fit->GetParameter(0);
    sigmafit.enoise = fit->GetParError (0);
    sigmafit.stant = fit->GetParameter(1);
    sigmafit.estant = fit->GetParError (1);
    sigmafit.stoch = doSterm ? fit->GetParameter(2) : -99;
    sigmafit.estoch = doSterm ? fit->GetParError (2) : -99;
	std::cout << "Sigma fit resutls: " << std::endl;
    std::cout << " noise " << sigmafit.noise << " +/- " << sigmafit.enoise;
    std::cout << " stoch " << sigmafit.stoch << " +/- " << sigmafit.estoch;
    std::cout << " stant " << sigmafit.stant << " +/- " << sigmafit.estant << std::endl;

	for( auto& rmap : ResultsMap ){ rmap.second->Write( rmap.second->GetName(), TObject::kOverwrite ); delete rmap.second;}
	hist2D->Write( hist2D->GetName(), TObject::kOverwrite );
	delete form;
	delete fit;
	resTFile->Close();

	return sigmafit;

}//<<>>SigmaFitResult KUCMSTimeCalibration::runTimeFitter( TH2F* hist2D )

inline void KUCMSTimeCalibration::doResTimeFits( bool doLocal ){

    std::cout << "Making Resolution Time Fits Auto : ";
	if( doLocal ) std::cout << " Local " << std::endl;
	else std::cout << " Global " << std::endl;
	for( auto& runmaptag : CaliRunMapSet ){
		for( auto& runmap : runmaptag.second ){
			bool doResFit = not runmap.second.hasResParams;
			if( runmap.second.has2DResMap && doResFit ){

				std::string histtype = doLocal ? lochist : globhist;
                std::string fname = runmap.second.histMapName+histtype;
				std::cout << " Opening : " << fname << std::endl;
				SigmaFitResult results = runTimeFitter( CaliHists[fname].h2f );
                //std::string gbfname = runmap.second.histMapName+globhist;
				//SigmaFitResult results = runTimeFitter( CaliHists[gbfname].h2f );
				if( not doLocal ){
					runmap.second.hasResParams = true;
					runmap.second.noise = results.noise;
					runmap.second.stoch = results.stoch;
					runmap.second.stant = results.stant;
				}//<<>>if( not doLocal )
			
			}//<<>>if( runmap.second.has2DResMap && not runmap.second.hasResParams )
		}//<<>>for( auto& runmaptag :  CaliRunMapSet )
	}//<<>>for( auto& runmaptag :  CaliRunMapSet )
	if( cali2DResTFile != NULL ) cali2DResTFile->Close();

}//<<>>void KUCMSTimeCalibration::doResTimeFits()

inline void KUCMSTimeCalibration::load2DResHist( std::string histName ){

    std::cout << " --- Load hist from caliTFile: " << histName << std::endl;
	cali2DResTFile = TFile::Open( cali2DResPlotsTFileName.c_str(), "READ" );
	TH2F* lhist = (TH2F*)cali2DResTFile->Get(histName.c_str());
    //std::cout << " --- Found local " << lhist << std::endl;
    if( lhist ) CaliHists[histName] = { lhist, histName, false, true }; // histfile histname isnew isreshist
	else std::cout << " --- No such hist found " << lhist << std::endl;

}//<<>>void KUCMSTimeCalibration::load2DResHist( std::string histName )

inline void KUCMSTimeCalibration::doResTimeFit( std::string histName ){

    std::cout << "Making Resolution Time Fits  : " << histName << std::endl;
	if( CaliHists.find(histName) != CaliHists.end() ){
		if( CaliHists[histName].h2f != NULL ){
    		SigmaFitResult results = runTimeFitter( CaliHists[histName].h2f );
		} else std::cout << " -- No H2F present !! " << std::endl;
	} else std::cout << " -- No Such hist entry !! " << std::endl;
	if( cali2DResTFile != NULL ) cali2DResTFile->Close();

}//<<>>void KUCMSTimeCalibration::doResTimeFit( std::string histName )

inline void KUCMSTimeCalibration::makeSmearTag( std::string sourceName, std::string destName, std::string smearTag ){

    std::cout << "Make Semar Tag : " << smearTag << std::endl;
    bool doSterm = true;

    std::string outfilename = caliFileDir + smearTag + "_resfit.root";
    TFile* smearTFile = TFile::Open( outfilename.c_str(), "UPDATE" );

    std::string sinfilename = caliFileDir + sourceName + "_resfit.root";
    TFile* sourceTFile = TFile::Open( sinfilename.c_str(), "READ" );
    if( not sourceTFile ){ std::cout << "Failed to open : " << sinfilename << std::endl; return; }

    std::string dinfilename = caliFileDir + destName + "_resfit.root";
    TFile* destTFile = TFile::Open( dinfilename.c_str(), "READ" );
	if( not destTFile ){ std::cout << "Failed to open : " << dinfilename << std::endl; return; }

    std::string fTitle = "#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)";
    std::string fXTitle = "A_{eff}/#sigma_{n} (EBEB)";
    std::string fTimeText = "t_{1}-t_{2}";
    std::string fSigmaVarText = "A_{eff}/#sigma_{n}";

    std::vector<float> fXBins;
    setBins( xBinStr, fXBins );
    int fNBinsX = fXBins.size();
    int nXbins = fNBinsX - 1;
    const auto xbins = &fXBins[0];

	smearTFile->cd();
    std::string sigtext = "#sigma("+fTimeText+") [ns]";
    std::string sigName = smearTag+"_sigma";
    std::string sigTitle = fTitle+" "+sigtext+";"+fXTitle+";"+sigtext;
    TH1F* Results = new TH1F(sigName.c_str(),sigTitle.c_str(),nXbins,xbins);
    Results->Sumw2();

	sourceTFile->cd();
    std::string sourceHName = sourceName+"_sigma";
	TH1F* Source = (TH1F*)sourceTFile->Get(sourceHName.c_str());
    if( not Source ){ std::cout << "Failed to open : " << sourceHName << std::endl; return; }

	destTFile->cd();
    std::string destHName = destName+"_sigma";
    TH1F* Dest = (TH1F*)destTFile->Get(destHName.c_str());
    if( not Dest ){ std::cout << "Failed to open : " << destHName << std::endl; return; }

    for (auto ibinX = 1; ibinX < fNBinsX; ibinX++){

		float ssigma = Source->GetBinContent( ibinX );
        float sserr = Source->GetBinError( ibinX );
        float dsigma = Dest->GetBinContent( ibinX );
        float dserr = Dest->GetBinError( ibinX );
		
		if( ssigma <=0 || dsigma <= 0 ) continue;
		float sigma = std::sqrt( pow(dsigma,2) - pow(ssigma,2) );
        float serr = std::sqrt( pow((dsigma*dserr/sigma),2) + pow((ssigma*sserr/sigma),2) );
		std::cout << " - smear sigma: " << sigma << " <- " << dsigma << " qd " << ssigma << " err " << serr << std::endl; 
        if( sigma <= 0 ) continue;

        Results->SetBinContent( ibinX, sigma );
        Results->SetBinError( ibinX, serr );

    }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

    // Prep sigma fit
    //----------------------TimeFitter::PrepSigmaFit(FitInfo);
    //std::cout << "Prepping sigma fit for: " << f2DHistName << std::endl;
    float nLower = 0;
    float nVal = 50;
    float nUpper = 100;
    float cLower = 0;
    float cVal = 0.5;
    float cUpper = 1;
    // get input hist
    float x_low = 25;
    float x_up = 1800;
    std::string histname = Results->GetName();
    std::string formname = histname+"_form";
    std::string fitname  = histname+"_fit";
    //std::cout << "Prepping sigma fit for: " << histname << std::endl;
    std::string fitformstr;
    if( doSterm ) fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))";
    else fitformstr = "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))";
    auto form = new TFormula(formname.c_str(),fitformstr.c_str());
    auto fit = new TF1(fitname.c_str(),form->GetName(),x_low,x_up);
    fit->SetParName( 0, "N" );
    fit->SetParameter( 0, nVal );
    fit->SetParLimits( 0, nLower, nUpper );
    fit->SetParName( 1, "C" );
    fit->SetParameter( 1, cVal );
    fit->SetParLimits( 1, cLower, cUpper );
    if( doSterm ) fit->SetParName( 2, "S" );
    if( doSterm ) fit->SetParameter( 2, 1.0 );
    if( doSterm ) fit->SetParLimits( 2, 0.0, 10.0 );
    fit->SetLineColor(Results->GetLineColor());
    Results->Fit(fit->GetName(),"RBQO");

    float noise = fit->GetParameter(0);
    float enoise = fit->GetParError (0);
    float stant = fit->GetParameter(1);
    float estant = fit->GetParError (1);
    float stoch = doSterm ? fit->GetParameter(2) : -99;
    float estoch = doSterm ? fit->GetParError (2) : -99;
    std::cout << "Sigma fit resutls: " << std::endl;
    std::cout << " noise " << noise << " +/- " << enoise;
    std::cout << " stoch " << stoch << " +/- " << estoch;
    std::cout << " stant " << stant << " +/- " << estant << std::endl;

	SmearTagSet[smearTag] = { noise, stoch, stant };

	smearTFile->cd();
	Results->Write( Results->GetName(), TObject::kOverwrite );

	delete Results;
    delete form;
    delete fit;
    smearTFile->Close();
	sourceTFile->Close();
	destTFile->Close();

}//<<>>void KUCMSTimeCalibration::makeSmearTag( std::string histName )

inline void KUCMSTimeCalibration::plotMeanRunTimeEGR( std::string inputFileName, int srun, int erun, bool usecali ){

    std::cout << "Creating Mean Run Time Hist from EgammaRes Ntuples " << std::endl;

    bool debug = false;
    //bool debug = true;

	std::string calistr = usecali ? "_cali" : "_nocali";
	std::string hname = "meanRunTime_" + std::to_string(srun) + "_to_" + std::to_string(erun)+calistr;
    std::string outfilename = caliFileDir + hname+ "_Hist.root";
    TFile* mrTFile = TFile::Open( outfilename.c_str(), "RECREATE" );
    mrTFile->cd();

    std::cout << " - output filename : " << outfilename << std::endl;

    int rrange = erun - srun;
	std::string title_eb = hname + ";Run;EB Mean Time [ns]";
	std::string hname_eb = hname + "_eb";
	TH1F* hist_eb = new TH1F(hname_eb.c_str(),title_eb.c_str(),rrange,srun,erun);
    std::string title_tt = hname + ";Run;TT 15:25 Mean Time [ns]";//ttid==1525
    std::string hname_tt = hname + "_tt1";
    TH1F* hist_tt = new TH1F(hname_tt.c_str(),title_eb.c_str(),rrange,srun,erun);
    std::string title_xt = hname + ";Run;TT -15:50 Mean Time [ns]";
    std::string hname_xt = hname + "_tt2";
    TH1F* hist_xt = new TH1F(hname_xt.c_str(),title_eb.c_str(),rrange,srun,erun);

	std::map<int,double> sum;
    std::map<int,double> sum2;
    std::map<int,long int> occ;
    std::map<int,double> tsum;
    std::map<int,double> tsum2;
    std::map<int,long int> tocc;
    std::map<int,double> xsum;
    std::map<int,double> xsum2;
    std::map<int,long int> xocc;

    const std::string treename("tree/llpgtree");

    // Declaration of leaf types
    uInt run;
    std::vector<unsigned int> *rhID;
    std::vector<float> *rhTime;
    std::vector<float> *rhEnergy;

    // List of branches
    TBranch        *b_run;   //!
    TBranch        *b_rhID;   //!
    TBranch        *b_rhTime;   //!
    TBranch        *b_rhEnergy;   //!

    std::cout << "open input files list : " << inputFileName << std::endl;

    std::ifstream infilelist(inputFileName);
    std::string infilestr;
    while( std::getline( infilelist, infilestr ) ) {

        std::stringstream ss(infilestr);
        std::string infilename, tag;
        int srun, erun;

        ss >> infilename >> srun >> erun >> tag;
        if( infilename[0] == '#' ) continue;
        std::cout << "open input file : " << infilename << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing Mean Run Time Hist " << std::endl;

        // insure we have Calimaps if we are doing Xtal calibration maps
        if( CaliRunMapSet.find(tag) == CaliRunMapSet.end() ){
            std::cout << " No Cali maps for this tag !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        std::ifstream infile(infilename);
        std::string instr;
        auto fInTree = new TChain( treename.c_str() );
        std::cout << "Adding files to TChain." << std::endl;
        while (std::getline(infile,instr)){
            auto tfilename = eosDir + inDir + instr;
            std::cout << "-";
            fInTree->Add(tfilename.c_str());
        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;

        run = 0;
    	rhID = 0;
        rhTime = 0;
		rhEnergy = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "rhID", &rhID, &b_rhID );   //!
        fInTree->SetBranchAddress( "rhRtTime", &rhTime, &b_rhTime);   //!
        fInTree->SetBranchAddress( "rhEnergy", &rhEnergy, &b_rhEnergy );   //!

        std::cout << " Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
        //if( nEntries > 1000000 ) nEntries = 1000000;
        if( debug ) nEntries = ( nEntries < 1000 ) ? nEntries : 1000;
        if( debug ) std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
        int report = nEntries/20;
        for (auto centry = 0U; centry < nEntries; centry++){

            if( centry%report == 0 or centry == 0){
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " " << (1000.0*static_cast<float>(centry)/static_cast<float>(nEntries))/10.0 << "%" << std::endl;
            }//<<>>if( centry%10000000 == 0 or centry == 0)

            auto entry = fInTree->LoadTree(centry);

            if(debug) std::cout << " - Start loop " << std::endl;

            b_run->GetEntry(entry);   //!
            b_rhID->GetEntry(entry);   //!
            b_rhTime->GetEntry(entry);   //!
			b_rhEnergy->GetEntry(entry);

			if( run > srun && run <= erun ){
				int nRecHits = rhID->size();		
				for( int idx = 0; idx < nRecHits; idx++ ){

					auto idinfo = DetIDMap[rhID->at(idx)];
					double btime = rhTime->at(idx);
					if( idinfo.ecal != ECAL::EB ) continue;
					double calibration = usecali ? getCalibration( rhID->at(idx), run, tag ) : 0;
					double time = btime - calibration;
					double time2 = time*time;
                    //if( rhEnergy->at(idx) < 5.0 || rhEnergy->at(idx) > 150.0 ) continue;
                    //std::cout << " getCalibration : " << rhID->at(idx) << " " << run << " " << tag; 
                    //std::cout << " cali: " << getCalibration( rhID->at(idx), run, tag );
                    //std::cout << " time : " << btime << " caliTime: " << time << std::endl; 

					if( sum.find(run) != sum.end() ){
						sum[run] += time;
                        sum2[run] += time2;
                        occ[run] += 1;
					} else {
                        sum[run] = time;
                        sum2[run] = time2;
                        occ[run] = 1;
					}//<<>>if( sum.find(run) != sum.end() )
					if( getTTId(rhID->at(idx)) == 1525 ){
                    	if( tsum.find(run) != tsum.end() ){
                        	tsum[run] += time;
                        	tsum2[run] += time2;
                        	tocc[run] += 1;
                    	} else {
                        	tsum[run] = time;
                        	tsum2[run] = time2;
                        	tocc[run] = 1;
                    	}//<<>>if( sum.find(run) != sum.end() )
					}//<<>>if( DetIDMap[rhID->at(idx)].i1 == && DetIDMap[rhID->at(idx)].i2 == )
                    if( getTTId(rhID->at(idx)) == 3550 ){
                    	if( xsum.find(run) != xsum.end() ){
                        	xsum[run] += time;
                        	xsum2[run] += time2;
                        	xocc[run] += 1;
                    	} else {
                        	xsum[run] = time;
                        	xsum2[run] = time2;
                        	xocc[run] = 1;
                    	}//<<>>if( sum.find(run) != sum.end() )
                    }//<<>>if( DetIDMap[rhID->at(idx)].i1 == && DetIDMap[rhID->at(idx)].i2 == )

				}//<<>>for( int idx = 0; idx < nRecHits; idx++ )
			}//<>>if( run > srun && run <= erun )
        } // for (auto entry = 0U; entry < nEntries; entry++

        if(debug) std::cout << " -------- Next Input file " << std::endl;

    } // while (std::getline(infilelist,infiles))

	for( auto & runsum : sum ){

		int run = runsum.first;
		int bin = 1 + run - srun;
		float mean = occ[run] > 0 ? sum[run] / occ[run] : -99;
		float err = mean > -99 ? sqrt( (sum2[run]/occ[run] - mean*mean)/occ[run] ) : 0;	
		std::cout << " EB mean by Run : " << run << " sum: " << sum[run] << " occ: " << occ[run]; 
        std::cout << " mean: " << mean << " +/- " << err << std::endl;			
		hist_eb->SetBinContent( bin, mean );
		hist_eb->SetBinError( bin, err );
        mean = tocc[run] > 0 ? tsum[run] / tocc[run] : -99;
        err =  mean > -99 ? sqrt( (tsum2[run]/tocc[run] - mean*mean)/tocc[run] ) : 0;
        hist_tt->SetBinContent( bin, mean );
        hist_tt->SetBinError( bin, err );
        mean = xocc[run] > 0 ? xsum[run] / xocc[run] : -99;
        err = mean > -99  ? sqrt( (xsum2[run]/xocc[run] - mean*mean)/xocc[run] ) : 0;
        hist_xt->SetBinContent( bin, mean );
        hist_xt->SetBinError( bin, err );

	}//<<>>for( auto & runsum : sum )

	hist_eb->Write( hist_eb->GetName(), TObject::kOverwrite );
	delete hist_eb;
    hist_tt->Write( hist_tt->GetName(), TObject::kOverwrite );
    delete hist_tt;
    hist_xt->Write( hist_xt->GetName(), TObject::kOverwrite );
    delete hist_xt;

	mrTFile->Close();

}//<<>> void plotMeanRunTimeEGR( std::string indir, std::string infilelistname, 

inline void KUCMSTimeCalibration::makeTTDiffMaps( int srun, int erun ){

    std::cout << " - Making TTCaliDiffMaps " << std::endl;
    //TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),34,0,34,72,0,72);
    for( auto& calirunmap : TTCaliRunMapSet ){
		int previous = 0;
		std::string prevname = "";
		int first = 0;
        std::string firstname = "";
        for( auto& calirunsct : calirunmap.second ){
			if( previous == 0 ){ 
				previous = calirunsct.first;
				first = calirunsct.first; 
				prevname = calirunsct.second.histMapName; 
				firstname = calirunsct.second.histMapName;
				continue; 
			}
			if( first < srun || first > erun ) continue;
            //std::string tfilename = calirunsct.second.histMapName + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
            std::string tfilename = calirunsct.second.histMapName;
            //std::cout << " -- opening : " << tfilename << std::endl;
            previous = 305112;
            prevname = "EG_EOY_MINI_305112_TT";
            std::string mfilename = tfilename + "_MeanMap";
            std::string pfilename = prevname + "_MeanMap";
			std::cout << " firsts: " << prevname << " current: " << mfilename << std::endl;
            std::string dfilename = tfilename + "_" + std::to_string(previous) + "_DiffMap";
            std::string ddfilename = tfilename + "_" + std::to_string(previous) + "_DiffDist";
            //std::cout << " --- with : " << sfilename << std::endl;
            if( CaliHists.find(dfilename) == CaliHists.end() ){
				TH2F* hist = new TH2F(dfilename.c_str(),dfilename.c_str(),35,-17,18,72,1,73);
				//TH2F* hist = new TH2F(dfilename.c_str(),dfilename.c_str(),34,0,34,72,0,72);
                CaliHists[dfilename] = { hist, dfilename, false, false }; // histfile histname isnew isopen lastrun
				TH1F* dhist = new TH1F(ddfilename.c_str(),ddfilename.c_str(),100,-0.2,0.2);
                CaliHists[ddfilename] = { dhist, ddfilename, false, false }; // histfile histname isnew isopen lastrun
			}//<<>>if( CaliHists.find(dfilename) == CaliHists.end() )
            if( CaliHists.find(mfilename) != CaliHists.end() ){
                //std::cout << " --- found : " << tfilename << std::endl;
                for( int ieta = 1; ieta < 36; ieta++ ){
                    for( int iphi = 1; iphi < 73; iphi++ ){
                        int i1 = ieta - 18;
                        if( i1 == 0 ) continue;
                        //std::cout << " -- Processing : " << ieta << " " << iphi;
                        float mean  = CaliHists[mfilename].h2f->GetBinContent(ieta,iphi);
                        float pmean  = CaliHists[pfilename].h2f->GetBinContent(ieta,iphi);
						float diff = mean - pmean;
						CaliHists[dfilename].h2f->Fill( i1, iphi, diff );
                        CaliHists[ddfilename].h1f->Fill( diff );
                    }//<<>>for( int iphi = 1; iphi < 361; iphi++ )
                }//<<>>for( int ieta = 1; ieta < 172; ieta++ )
            }//<<>>if( CaliHists.find(sfilename) != CaliHists.end() )
			previous = calirunsct.first;
			prevname = calirunsct.second.histMapName;
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

}//<<>>void KUCMSTimeCalibration::makeTTDiffMaps()

inline void KUCMSTimeCalibration::makeTTDriftMaps( std::string tag, int srun, int erun ){

    std::cout << " - Making makeTTDriftMaps " << std::endl;
	if( erun < srun ){ std::cout << " --- End Run less then Start Run : Exiting " << std::endl; return; }
	if( not ( TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() ) ){

		std::cout << " -- Initing TH1F and TGraphErrors " << std::endl;
		auto& calirunmap = TTCaliRunMapSet[tag];
        int first = 0;
        std::string firstname = "";
		std::map<int,TH1F*> histmap;
		//std::map<int,TGraphErrors*> graphmap;
        std::map<int,TGraph*> graphmap;
		for( int ieta = 1; ieta < 37; ieta++ ){
			if( ieta == 18 ) continue; 
			for( int iphi = 1; iphi < 74; iphi++ ){
				int index = ieta*100 + iphi;
				int runrange = erun - srun + 2;
				int erange = erun + 1;
				int srange = srun - 1;
				std::string filename = "TTmeanDrift_" + std::to_string(index);
				histmap[index] = new TH1F(filename.c_str(),filename.c_str(),runrange,srange,erange);
				histmap[index]->Sumw2();
				graphmap[index] = new TGraphErrors(runrange);
			}//<<>>for( int iphi = 1; iphi < 73; iphi++ )
		}//<<>>for( int iphi = 1; iphi < 73; iphi++ )

		std::cout << " -- Filling TH1F and TGraphErrors " << std::endl;
        std::map<int,float> normrun;
		float maxtime = 0;
        for( auto& calirunsct : calirunmap ){

			first = calirunsct.first;
            if( first < srun || first > erun ) continue;
			//std::cout << " - Processing : " << calirunsct.second.histMapName << " in " << first << std::endl;

            for( int ieta = 1; ieta < 37; ieta++ ){
                for( int iphi = 1; iphi < 74; iphi++ ){
                    int i1 = ieta - 18;
                    if( i1 == 0 ) continue;
					uInt detid = getInvTTId( iphi, i1 );
					int index = ieta*100 + iphi;
					int bin = first - srun + 1;
					float mean = calirunsct.second.meanMap[detid];
                    float error = calirunsct.second.errMap[detid];
                    if( bin == 1 ) normrun[index] = mean;
					if( mean == 0 ) continue;
					if( std::abs(mean) > maxtime ) maxtime == std::abs(mean);
                    //std::cout << " -- Filling : " << ieta << " " << iphi << " " << index << " " << first << " ";
					//std::cout << bin << " " << mean << " " << error << " " << normrun[index] << std::endl;
					histmap[index]->SetBinContent( bin, mean );
                    histmap[index]->SetBinError( bin, error );
					graphmap[index]->SetPoint( bin, first, mean - normrun[index] );
                    //graphmap[index]->SetPointError( bin, 0, error);
                }//<<>>for( int iphi = 1; iphi < 361; iphi++ )
            }//<<>>for( int ieta = 1; ieta < 172; ieta++ )

        }//<<>>for( auto& calirunsct : calirunmap )

		std::cout << " -- Filling & Printing TMultiGraph " << std::endl;
		std::string driftHistTCanvasName = "driftHists_" + tag + "_" + std::to_string( srun ) + "_" + std::to_string( erun ) + ".png";
		gROOT->SetBatch(true);
		TCanvas c1 = TCanvas( "c1", "canvas" );
    	c1.SetGridx(1);
    	c1.SetGridy(1);
		TMultiGraph* mg = new TMultiGraph();
        std::map<int,int> colors = {{0,kMagenta+2},{1,kRed+2},{2,kYellow+2},{3,kGreen+2},{4,kCyan+2},{5,kBlue+2},
									{6,kViolet+2},{7,kPink+2},{8,kOrange+2},{9,kSpring+2},{10,kTeal+2},{11,kAzure+2},
									{23,kAzure+7},{12,kMagenta-4},{13,kRed-4},{14,kYellow-4},{15,kGreen-4},{16,kCyan-4},
									{17,kBlue-4},{18,kViolet+7},{19,kPink+7},{20,kOrange+7},{21,kSpring+7},{22,kTeal+7},
                                    {34,kTeal-6},{35,kAzure-6},{24,kMagenta-6},{25,kRed-6},{26,kYellow-6},{27,kGreen-6},
                                    {28,kCyan-6},{29,kBlue-6},{30,kViolet-6},{31,kPink-6},{32,kOrange-6},{33,kSpring-6}
									};
		//std::map<int,int> markers = {{0,20},{1,21},{2,22},{3,23},{4,29},{5,33},{6,34},{7,39},{8,45},{9,47}};
        std::map<int,int> markers = {{0,29},{1,33},{2,34},{3,39},{4,41},{5,43},{6,45},{7,47},{8,48},{9,49},{10,22},{11,23}};	
     	int n = 0, k=0;
		int sEta = 1, eEta = 37; // 35 bins : -17 to 18 ( skip 0 ) -> 1 to 36 ( skip 18 )
        //int sEta = 32, eEta = 33;
        //int sEta = 12, eEta = 13;
        int sPhi = 1, ePhi = 74; // 72 bins : 1 to 73
        //int sPhi = 25, ePhi = 26;
        //int sPhi = 50, ePhi = 51;
        for( int ieta = sEta; ieta < eEta; ieta++ ){
        //for( int ieta = 32; ieta < 33; ieta++ ){
            if( ieta == 18 ) continue;
            for( int iphi = sPhi; iphi < ePhi; iphi++ ){
                int index = ieta*100 + iphi;
                graphmap[index]->SetMarkerStyle(markers[k]);
				graphmap[index]->SetLineColor(colors[n]);
				graphmap[index]->SetMarkerColor(colors[n]);
				graphmap[index]->SetMarkerSize(1.0);
				//graphmap[index] = new TGraphErrors( histmap[index] );
				mg->Add(graphmap[index]);
				n++;
				if( n > 35 ){ n = 0; k++; }
				k++;
				if( k > 11 ){ k = 0; }
            }//<<>>for( int iphi = 1; iphi < 73; iphi++ )
        }//<<>>for( int iphi = 1; iphi < 73; iphi++ )
        mg->GetXaxis()->SetNoExponent();
		mg->GetXaxis()->SetRangeUser(srun,erun);
    	mg->GetXaxis()->CenterTitle(true);
    	mg->GetXaxis()->SetTitle("Run");
        mg->GetYaxis()->SetRangeUser(-1.75,1.75);
    	mg->GetYaxis()->CenterTitle(true);
    	mg->GetYaxis()->SetTitle("normlized TT mean time [ns]");
        //mg->GetXaxis()->SetNoExponent();
        mg->SetTitle(driftHistTCanvasName.c_str());
		c1.Update();
		mg->Draw("AP");
		c1.Print( driftHistTCanvasName.c_str() );
		c1.Close();
		
    	std::cout << " -- Save Drift Hists" << std::endl;
    	std::string driftHistTFileName = "driftHists_" + tag + "_" + std::to_string( srun ) + "_" + std::to_string( erun ) + ".root";
    	TFile* driftHistTFile = TFile::Open( driftHistTFileName.c_str(), "UPDATE" );
    	driftHistTFile->cd();
    	for( auto& hist : histmap ){
        	hist.second->Write( hist.second->GetName(), TObject::kOverwrite );
        	delete hist.second;
    	}//<<>>for( auto& hists : CaliHists )
		//mg->Write("TT_TimeStability", TObject::kOverwrite );
		delete mg;

    	driftHistTFile->Close();

    }//<<>>if( TTCaliRunMapSet.find(tag) != CaliHists.end() )

    std::cout << "Finished making TT Mean Drift Hists" << std::endl;

}//<<>>void KUCMSTimeCalibration::makeTTDriftMaps()

#endif
//-------------------------------------------------------------------------------------------------------------------

