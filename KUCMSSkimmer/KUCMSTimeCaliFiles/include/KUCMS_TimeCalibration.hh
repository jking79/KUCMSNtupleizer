//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

//#include "KUCMS_HelperBaseClass.hh"
#include "KUCMS_RootHelperBaseClass.hh"
#include "KUCMS_EcalDetIDFunctions.hh"
#include "KUCMS_CaliRunClass.hh"
#include "KUCMS_TimeCaliStructs.hh"
#include "KUCMS_ProfileTimeFit.hh"

#include <TRandom.h>
#include "TChain.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include <cmath>
#include <ctime>

#ifndef KUCMS_TimeCalibrationClass_hh
#define KUCMS_TimeCalibrationClass_hh

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

class KUCMS_TimeCalibration : public KUCMS_RootHelperBaseClass {

	public:

    KUCMS_TimeCalibration( bool stayOpen = false, bool makeNew = false  );
	~KUCMS_TimeCalibration();


	private:

	enum ECAL {EB, EM, EP, NONE};

    std::string detIDConfigEB;
    std::string detIDConfigEE;

    std::string caliFileDir;
    std::string caliRunConfig;
    std::string caliTTConfig;
    std::string caliTFileName;
    std::string caliR3TFileName;
    std::string cali2DResPlotsTFileName;
    std::string caliResConfig;

    std::string xtalHistMapName;
    std::string ttHistMapName;
    std::string lochist;
    std::string globhist;

    TFile* caliTFile;
    TFile* caliR3TFile;
    TFile* cali2DResTFile;

	std::string eosDir;
	std::string inDir;

    std::string curTTIov;
    std::string curXIov;
	std::string curTag; // used to change rereco version or campaian calibrations are for
    std::string curLumiTag;

    std::string resTag; // NOTE: resolution used for PDs in set in master file lists, this is resoltion of MC for smearing

	std::string xBinStr;
    std::string yBinStr;

	std::vector<std::string> maptypes;
    std::vector<std::string> ebeemap;

    std::map<UInt_t,kucms_DetIDStruct> DetIDMap; // map of information by detid for each crystal
    std::map<int,std::map<int,std::map<int,uInt>>> InvDetIDMap;

	// do we need TT and X verions ?  yes for set and no for hist maps ?
	std::map<std::string,std::map<int,kucms_CaliRunClass>> CaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map
    std::map<std::string,std::map<int,kucms_CaliRunClass>> TTCaliRunMapSet; // str is tag, int is internal iov/era number  - final cali map

	// interface between TFile and RunMapSet ?  off load to meanMap?
	std::map<std::string,kucms_caliHistStruct> CaliHists; // str is label for specific calibration histogram - stored as ref in calirunmap
    std::map<std::string,std::map<int,kucms_lumiRunStruct>> lumiRunMaps;// str is settag, int is run - ref for maps - give lumi by run
    std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>

    std::map<std::string,kucms_resTagStruct> ResTagSet;

	std::vector<int> ebx_ranges, ebtt_ranges, epmx_ranges, epmtt_ranges;

	bool updated;
	bool lowEnergy;
	bool useEffEnergy;
    bool externalCali;
    bool useGSwitch;
	bool doEE;
	int useGain;
	bool doAllGain;
	bool isCC;
    bool doUnCC;

    TRandom* getRandom;

	public:

    // -------------------------------------------------------------------------
    // Setup and intilization functions --  primarially used internally 
    // --------------------------------------------------------------------------

	void SetupDetIDsEB();
	void SetupDetIDsEE();
	std::map<UInt_t,kucms_DetIDStruct>& getDetIDMapRef(){ return DetIDMap; };

	void SetupIovMaps();
    void SetupIovMap( std::string tag, float lumiMin );

	void ReadCaliRunFile();
    void SaveCaliRunFile();
    void ReadTTRunFile();
    void SaveTTRunFile();
	void ReadResFile();
    void SaveResFile();
	void SetResParamters( std::string tag, float ebnoise, float ebstoch, float ebstant, float eenoise, float eestoch, float eestant );
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


	// -------------------------------------------------------------------------
	// accessors to utilize calibration information  ::  for use in KUSkimmer
    // --------------------------------------------------------------------------

	float getCorrectedTime( float time, float amplitude, unsigned int rechitID, 
			unsigned int Evt_run, std::string dataSetKey, int mctype, int gainID = 1 ); 
	float getCorrectedTime( float time, float amplitude, unsigned int rechitID, 
			unsigned int Evt_run, std::string dataSetKey, int mctype, bool gs6, bool gs1 )
			{ int gid = ( gs6 || gs1 ) ? 2 : 1; return getCorrectedTime( time, amplitude, rechitID, Evt_run, dataSetKey, mctype, gid ); };
    float getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey, int mctype );

	float getCalibration( uInt rhid, int run, std::string tag, int gainID = 1 ); // tag indicates which calibration set to use
    float getCalibration( uInt rhid, int run, int gainID = 1 )
			{ return getCalibration( rhid, run, curTag, gainID ); };
    float getCalibration( uInt rhid, int run, bool gs6, bool gs1 ) 
            { int gid = ( gs6 || gs1 ) ? 2 : 1;  return getCalibration( rhid, run, curTag, gid ); }; 
	float getCalibTime( float rhtime, uInt rhid, int run, int gainID = 1 )
			{ return rhtime - getCalibration( rhid, run, gainID ); };
	float getCalibTime( float rhtime, uInt rhid, int run, bool gs6, bool gs1 )
			{ int gid = ( gs6 || gs1 ) ? 2 : 1; return rhtime - getCalibration( rhid, run, gid ); };

	//tag indicates which smear to use  ---  this is all depreciated - needs to be reworked for current statagy
	float getSmearedTime(  float rhtime, float rhamp, std::string stag );
    float getSmearedTime(  float rhtime, float rhamp )
            { return getSmearedTime( rhtime, rhamp, resTag ); };
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run, std::string ctag, std::string stag );
    float getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int run )
			{ return getSmrdCalibTime( rhtime, rhamp, rhid, run, curTag, resTag ); }; 


	// --------------------------------------------------------------------------
	// accessor function for use in fillTimeCali - Setup and processing functions
	// --------------------------------------------------------------------------

	float getTTCali( uInt rhid, int run, std::string tag, int gainID = 1 );

	kucms_DetIDStruct& getDetIdInfo( uInt rhid );
	std::map<int,std::map<int,std::map<int,uInt>>> getInvDetIDMap();
    std::map<UInt_t,kucms_DetIDStruct> getDetIDMap(); // map of information by detid for each crystal
	uInt getDetIdInfo( int i1, int i2, int ecal );
	uInt getTTId( uInt detId, int gainID = 1 );
	std::pair<int,int> getTTInfo( uInt ttid );
    uInt getInvTTId( int i1, int i2, bool isEB );

	void setTag( std::string tag ){ curTag = tag; }; // used to change rereco version or campaian calibrations are for
    void setXIov( std::string tag ){ curXIov = tag; };
    void setTTIov( std::string tag ){ curTTIov = tag; };
    void setLumiTag( std::string tag ){ curLumiTag = tag; };

    void setMCResTag( std::string tag ){ resTag = tag; }; // NOTE: resolution for PDs set in master file lists,this is resoltion of MC for smearing

    void SetXBinStr( std::string xbins ){ xBinStr = xbins; };
    void SetYBinStr( std::string ybins ){ yBinStr = ybins; };
	void SetEosDir( std::string eosdir ){ eosDir = eosdir; }
    void SetInDir( std::string indir ){ inDir = indir; };

	void Set2DResTFileName( std::string name ){ cali2DResPlotsTFileName = name; };
    void SetCaliTFileName( std::string name ){ caliTFileName = name; };
    void SetCaliFileDir( std::string name ){ caliFileDir = name; };

	// function to create calibration and smear information 
	// make TT maps from egamma res ntuple from lpc/jwk space
	void makeCaliMapsEGR( std::string inputFileName, bool doTT, int GID = 1, bool small = false, bool doCali = true );
	void makeTTCaliMapEGR( std::string inputFileName, int GID = 1, bool small = false, bool doCali = true ) 
			{ makeCaliMapsEGR( inputFileName, true, GID, small, doCali ); }; 
	void makeXCaliMapEGR( std::string inputFileName, int GID = 1, bool small = false, bool doCali = true )
			{ makeCaliMapsEGR( inputFileName, false, GID, small, doCali ); };

	void plot2dResolutionEGR( std::string inputFileName, bool scale = true, bool usecali = true, bool smear = false, std::string ext = "" );
    void plot2dResbyIovForEGR( std::string inputFileName, bool scale = true, bool usecali = true, bool smear = false, std::string ext = "" );
    kucms_SigmaFitResult runTimeFitter( TH2F* hist2D );
	void doResTimeFits( bool doLocal = false );
    void doResTimeFit( std::string histName );
	void load2DResHist( std::string histName );

	TH1F* gprProfileHist( TH1F* hist );

    void makeSmearTag( std::string sourceName, std::string destName, std::string smearTag );

	void plotMeanRunTimeEGR( std::string inputFileName, int srun, int erun, bool usecali = true );
    void makeTTDiffMaps( int srun, int erun );
    void makeTTDriftMaps( std::string tag, int srun, int erun );

	void setUseEffEnergy( bool setting ){ useEffEnergy = setting; };
	void setLowEnergy( bool setting ){ lowEnergy = setting; }; 
	void useExtCali( bool ext ){ externalCali = ext; };
    void useGainSwitch( bool gs, bool all ){ useGSwitch = gs; doAllGain = all; };
    void setDoEE( bool doee ){ doEE = doee; }
    void setUseGainID( int id ){ useGain = id; }
    void setIsCC( bool iscc ){ isCC = iscc; }
	void setDoUnCC( bool douncc ){ doUnCC = douncc; }

};//<<>>class KUCMS_TimeCalibration : KUCMSRootHelperBaseClass

#endif
