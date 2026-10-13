#include "../include/KUCMS_TimeCalibration.hh"

#include <unordered_map>
#include <unordered_set>

namespace {

Long64_t progressEntryForPercent( Long64_t totalEntries, int percent ){
	const Long64_t whole = totalEntries/100;
	const Long64_t remainder = totalEntries%100;
	return whole*percent + ( remainder*percent + 99 )/100;
}

void reportEventLoopProgress( Long64_t processed, Long64_t totalEntries, int percent ){
	std::time_t now = std::time(nullptr);
	std::string curtime = std::ctime(&now);
	if( not curtime.empty() ) curtime.pop_back();
	std::cout << "Processed " << processed << " of " << totalEntries << " entries ("
	          << percent << "%) at " << curtime << std::endl;
}

}//<<>>namespace

#ifndef KUCMS_TimeCalibrationClass_cc
#define KUCMS_TimeCalibrationClass_cc

///////////////////////////////////////////////////////////////////////////////////////////////////
//  Class Object code  --  yes yes this is easier for me, im weird, will divide into hh/cc at end 
//  		--  i lied : leaving it all in one file - kinda of
///////////////////////////////////////////////////////////////////////////////////////////////////

KUCMS_TimeCalibration::KUCMS_TimeCalibration( bool stayOpen, bool makeNew, bool doLocalCali ){

    std::cout << "Initiating KUCMS_TimeCalibrationClass for EE && EB" << std::endl;
	// parameter and setup hardcoded in constructer - some of this needs moved to "run" code

    caliFileDir = "ecal_config/";// move
	detIDConfigEB = caliFileDir + "fullinfo_v2_detids_EB.txt";
    detIDConfigEE = caliFileDir + "fullinfo_v2_detids_EE.txt";
	caliRunConfig = caliFileDir + "caliRunConfig.txt";
    caliTTConfig = caliFileDir + "caliTTConfig.txt";
    caliResConfig = caliFileDir + "caliResConfig.txt";

	isCondorJob = false;

    xtalHistMapName = "_X_";
    ttHistMapName = "_TT_";
    lochist = "_SRO_Data_Hist";
    globhist = "_ZEE_Data_Hist";

	ebx_ranges = {171,-85,86,364,1,365};
	epmx_ranges = {203,-101,102,100,1,101};
	ebtt_ranges = {35,-17,18,74,1,75};
    epmtt_ranges = {23,-11,12,34,1,35};

    //std::cout << " - opening caliHistsTFile " << std::endl;
    //caliTFileName = caliFileDir + "caliHistsTFile.root";// name configurable ?
	//caliR3TFileName = caliFileDir + "caliR3HistsTFile.root";// name configurable ?
	cali2DResPlotsTFileName = caliFileDir + "res2dPlotsTFile.root";
	//caliTFile = NULL;
    caliR3TFile = NULL;
	cali2DResTFile = NULL;

    ///////caliHistFileNames["r2"] = "caliHistsTFile.root";
    ///////caliHistFileNames["r3"] = "caliR3HistsTFile.root";

	caliHistFileNames["r2_ul16"] = "cali_ul16_HistsTFile.root";
    caliHistFileNames["r2_ul17"] = "cali_ul17_HistsTFile.root";
    caliHistFileNames["r2_ul18"] = "cali_ul18_HistsTFile.root";
    ////////caliHistFileNames["r2_ul18_mc"] = "cali_ul18_mc_HistsTFile.root";
    caliHistFileNames["r3_p22"] = "cali_p22_HistsTFile.root";
    caliHistFileNames["r3_p23"] = "cali_p23_HistsTFile.root";
    caliHistFileNames["r3_p24"] = "cali_p24_HistsTFile.root";
    caliHistFileNames["r3_p25"] = "cali_p25_HistsTFile.root";
    caliHistFileNames["r3_p25unc"] = "cali_p25_UnCor_HistsTFile.root";
    caliHistFileNames["r3_p26"] = "cali_p26_HistsTFile.root";

	useEosRootCali =  doLocalCali ? false : true;
	caliMapsValid = true;
	////////useEosRootCali = false;
    eosCaliPath = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/caliconfig/";

	for( auto name : caliHistFileNames ){ caliTFile[name.first] = NULL; }

    curTag = "EG_EOY_MINI";
	// PD: EG_EOY_MINI    MC: RunIIFall17DRPremix 
	curTTIov = "r2ul";
    curXIov = "prompt";
	curLumiTag = "r2ul";

	//resTag = "default";  // NOTE: resolution used for PDs in set in master file lists, this is targ resoltion of MC for smearing
    resTag = "r2_ul18D";

	eosDir = "root://cmseos.fnal.gov//store/user/jaking/";
	inDir = "";

	lowEnergy = true;
	useEffEnergy = false;
	xBinStr = "VARIABLE 25 50 75 100 125 150 175 200 225 250 300 400 600 1000 1800"; 
    yBinStr = "CONSTANT 600 -3 3";

	maptypes = {"_MeanMap","_ErrMap","_SumMap","_Sum2Map","_OccMap"};
	ebeemap = {"_EB","_EE"};

	updated = false;
	externalCali = false;
    useGSwitch = true;

	doEE = false;
	useGain = 1;
	doAllGain = true;
	isCC = false;
	doUnCC = false;

    std::cout << " - setup DetID & IOV Maps " << std::endl;
    SetupDetIDsEB();
    SetupDetIDsEE();

    getRandom = new TRandom();
    getRandom->SetSeed(0);

    std::cout << " - reading in TTrun, Xrun, & smear cali files " << std::endl;
	ReadCaliRunFile();
	ReadTTRunFile();
	ReadResFile();

    //std::cout << " - loading cali hists " << std::endl;
	caliMapsValid = LoadCaliHists( stayOpen, makeNew );

}//<<>>KUCMS_TimeCalibration()   

// move save functions external - iopen close TFile in read/save functions
KUCMS_TimeCalibration::~KUCMS_TimeCalibration(){

    std::cout << "Wrapping KUCMS_TimeCalibrationClass" << std::endl;

    delete getRandom;

    //std::cout << "Finished Wrapping KUCMS_TimeCalibrationClass" << std::endl;

}//<<>>KUCMS_TimeCalibration::~KUCMS_TimeCalibration()

void KUCMS_TimeCalibration::SetupDetIDsEB(){

    std::ifstream infile( detIDConfigEB, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, iphi, ieta, absieta, FED, SM, TT25, iTT, strip5, Xtal, phiSM, etaSM;
    float phi, eta;
    std::string pos;

    while( infile >> cmsswId >> dbID >> hashedId >> iphi >> ieta >> absieta >> pos
                >> FED >> SM >> TT25 >> iTT >> strip5 >> Xtal >> phiSM >> etaSM >> phi >> eta ){

		//std::cout << cmsswId << " " << dbID << " " << hashedId << " " << iphi << " " << ieta;
		//std::cout << " " << absieta << " " << pos << " " << FED; 
		//std::cout << " " << SM << " " << TT25 << " " << iTT << " " << strip5 << " " << Xtal << " " << phiSM;
        //std::cout << " " << etaSM << " " << phi << " " << eta << " " << std::endl;

        DetIDMap[cmsswId] = {iphi,ieta,TT25,0,phi,eta,SM,iTT};
		InvDetIDMap[iphi][ieta][0] = cmsswId;

    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEB( std::map<UInt_t,DetIDStruct> &DetIDMap )

void KUCMS_TimeCalibration::SetupDetIDsEE(){

    std::ifstream infile( detIDConfigEE, std::ios::in);
    unsigned int cmsswId, dbID;
    int hashedId, side, ix, iy, SC, iSC, Fed, TTCCU, strip, Xtal, quadrant;
    float phi, eta;
    std::string EE;

    while( infile >> cmsswId >> dbID >> hashedId >> side >> ix >> iy >> SC
                >> iSC >> Fed >> EE >> TTCCU >> strip >> Xtal >> quadrant >> phi >> eta ){

		//std::cout << cmsswId << " " << dbID << " " << hashedId << " " << side << " " << ix << " " << iy;
        //std::cout << " " << SC << " " << iSC << " " << Fed << " " << EE << " " << TTCCU << " " << strip;
		//std::cout << " " << Xtal << " " << quadrant << " " << phi << " " << eta << std::endl;

        int ec = 1;
        if( side > 0 ) ec = 2;
        DetIDMap[cmsswId] = {ix,iy,TTCCU,ec,phi,eta,quadrant,SC};
        InvDetIDMap[ix][iy][ec] = cmsswId;

		//std::cout << " -- " << DetIDMap[cmsswId].ecal << " " << DetIDMap[cmsswId].i1 << " " << DetIDMap[cmsswId].i2 << std::endl;
	
    }//<<>>while (infile >>

}//<<>>void SetupDetIDsEE( std::map<UInt_t,DetIDStruct> &DetIDMap )

uInt KUCMS_TimeCalibration::getGSRhId( uInt detID, int gianID ){

	if( gianID == 1 ){ std::cout << " !!!!!  In getGSRhId gianID == 1 !!!!" << std::endl; return detID; }
	auto det = DetIDMap.find(detID);
	if( det == DetIDMap.end() ){
		std::cerr << "ERROR: unknown ECAL detID in getGSRhId: " << detID << std::endl;
		return 0;
	}//<<>>if( detID missing )
	const auto& idinfo = det->second;
	if( idinfo.ecal == 0 ){
		if( idinfo.i2 < 0 ){ if( gianID == 2 ) return 206; else return 201; }
		else { if( gianID == 2 ) return 106; else return 101; }
	} else if( idinfo.ecal == 1 ){ 
		if( gianID == 2 ) return 116; else return 111;
	} else if( idinfo.ecal == 2 ){
        if( gianID == 2 ) return 226; else return 221;
	} else { std::cout << " !!!!!  In getGSRhId and idinfo.ecal not = 0, 1, or 2 !!!!" << std::endl; return 0; }   

}//<<>>uInt KUCMS_TimeCalibration::getGSRhId( uInt detID, int gianID )

void KUCMS_TimeCalibration::SetupIovMaps(){

	//Prompt :: These are calibration and era IOVS 
	std::map<int,int> promptIovMap;
	//promptIovMap[1] = 189158;//03/02/2012
	//promptIovMap[189159] = 203827;//26/03/2012
    //promptIovMap[203828] = 208838;//28/09/2012
    //promptIovMap[208839] = 253983;//11/12/2012

    //promptIovMap[253984] = 276811;//11/08/2015
    //promptIovMap[276812] = 293998;//14/07/2016

	promptIovMap[1] = 253983; // 2012
    promptIovMap[253984] = 273157; // 2015

//gldenjson
//2016
//273158
//284044

    promptIovMap[273158] = 275656; 
    promptIovMap[275657] = 283830;
    promptIovMap[283831] = 284044; // 2016


//gldenjson
//2016
//273158
//284044

//2016C  	275657	275836
//2016D  	276315	276775
//2016E  	276831	277096
//2016F  	277992	278239 
//2016G  	278820	278969
//2016H 	281613	283830
// iov  253984 284044
// lumi 273158 284044

    promptIovMap[284045] = 296398; //16-17

//gldenjson
//2017
//297050
//306460

    promptIovMap[296399] = 297049;
    promptIovMap[297050] = 297113;
    promptIovMap[297114] = 297726;
    promptIovMap[297727] = 300201;
    promptIovMap[300202] = 301486;
    promptIovMap[301487] = 304475;
    promptIovMap[304476] = 305081;
    promptIovMap[305082] = 306460;//05/10/2017
  
//gldenjson
//2017
//297050
//306460


//2017B  	297114	297296
//2017C  	299368	301417
//2017D  	302031	302393
//2017E   	303832	304616
//2017F  	305044	305081
// iov  296399 306460
// lumi 297050 306460

    promptIovMap[306461] = 314766;// 17-18

//gldenjson
//2018
//315257
//325172

    promptIovMap[314767] = 315256;
    promptIovMap[315257] = 315343;
    promptIovMap[315344] = 316360;
    promptIovMap[316361] = 316569;
    promptIovMap[316570] = 318744;
    promptIovMap[318745] = 319329;
    promptIovMap[319330] = 321506;
    promptIovMap[321507] = 322356;
    promptIovMap[322357] = 322717;
    promptIovMap[322718] = 323412;
    promptIovMap[323413] = 324305;
    promptIovMap[324306] = 325172;
    promptIovMap[325173] = 327238;

//gldenjson
//2018
//315257
//325172

    promptIovMap[327239] = 355373;//25/11/2018  18-22

//2018A		315257	316993
//2018B		317435	317435
//2018C 	319756	320038
//2018D	    320673	322356
// iov  314767 315343
// lumi 315257 325172

//Run 3

// 2022
// lumi  355374 362760

    promptIovMap[355374] = 355793;
    promptIovMap[355794] = 356513; // Run2022C -v1
    promptIovMap[356514] = 357289;
    promptIovMap[357290] = 357500;
    promptIovMap[357501] = 358883;
    promptIovMap[358884] = 359000;
    promptIovMap[359001] = 359420;
    promptIovMap[359421] = 359560;
    promptIovMap[359561] = 359573;
    promptIovMap[359574] = 359688;
    promptIovMap[359689] = 360089;
    promptIovMap[360090] = 360981;
    promptIovMap[360982] = 361416;
    promptIovMap[361417] = 362180;
    promptIovMap[362181] = 362349;//24/11/2022
    promptIovMap[362350] = 362760;//EOY 2022 Run2022G -v1

// 2022
// lumi  355374 362760

    promptIovMap[362761] = 366364;// 22-23

/// 2023
// lumi 366727 370790 

    promptIovMap[366365] = 366726;
    promptIovMap[366727] = 367079;
    promptIovMap[367080] = 367515;
    promptIovMap[367516] = 367882;
    promptIovMap[367884] = 368825;
    promptIovMap[368826] = 369801;
    promptIovMap[369802] = 369912;//Run2023D -v1
    promptIovMap[369913] = 370496;
    promptIovMap[370497] = 370790;
	promptIovMap[370791] = 373577;
	promptIovMap[373578] = 373861;// EOY 2023  ? Run2023F

/// 2023
// lumi 366727 370790 

    promptIovMap[373862] = 378745;// 23-24

// 2024
// lumi 378985 386951

    promptIovMap[378746] = 378984; 
    promptIovMap[378985] = 381306;
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
    promptIovMap[386663] = 386951;
    promptIovMap[386952] = 387741;// EOY 2024 

// 2024
// lumi 378985 386951 

    promptIovMap[387742] = 390734;// 24-25

// 2025 - (390735) 391531 - 398903

    ////promptIovMap[390735] = 391530;

    promptIovMap[390735] = 390936;
    promptIovMap[390937] = 391260;
    promptIovMap[391261] = 391530;
    promptIovMap[391531] = 392158;
    promptIovMap[392159] = 392405;
    promptIovMap[392406] = 392485;
    promptIovMap[392486] = 393110;
    promptIovMap[393111] = 393146;
    promptIovMap[393147] = 393609;
    promptIovMap[393610] = 394033;
    promptIovMap[394034] = 394217;
    promptIovMap[394218] = 394285;
    promptIovMap[394286] = 394297;
    promptIovMap[394298] = 394458;
    promptIovMap[394459] = 395967;
    promptIovMap[395968] = 396597;
    promptIovMap[396598] = 397304;
    promptIovMap[397305] = 397853;
    promptIovMap[397854] = 398903;

// 2025 - (390735) 391531 - 398903

	////promptIovMap[391531] = 398903;

// END
	promptIovMap[398904] = 999999;
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

	float minttlumigs2 = 12.0;
	float minxlumigs2 = 60.0;

    //float minttlumigs3 = 36.0;
    //float minxlumigs3 = 400.0;

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
    std::string lumi25Config = "brilcalc_lumi_2025.txt";

    ReadLumiFile( caliFileDir+lumiUL16Config, r2ulLumiTag );
    ReadLumiFile( caliFileDir+lumiUL17Config, r2ulLumiTag );
    ReadLumiFile( caliFileDir+lumiUL18Config, r2ulLumiTag );
	setLumiTag(r2ulLumiTag);
    SetupIovMap( r2ulTagTT, minttlumi );
    SetupIovMap( r2ulTagX, minxlumi );
    setXIov("r2ulx");

    //std::string r2ulTagTTGS2( "r2ultt_gs2" );
	//std::string r2ulTagXGS2( "r2ulx_gs2" );

	//SetupIovMap( r2ulTagTTGS2, minttlumigs2 );
	//SetupIovMap( r2ulTagXGS2, minttlumigs2 );

    //std::string r2ulTagTTGS3( "r2ultt_gs3" );
    //std::string r2ulTagXGS3( "r2ulx_gs3" );

    //SetupIovMap( r2ulTagTTGS3, minttlumigs3 );
    //SetupIovMap( r2ulTagXGS3, minttlumigs3 );


    std::string r3TagTT( "r3tt" );
    std::string r3TagX( "r3x" );
    std::string r3LumiTag( "r3lumi" );

    ReadLumiFile( caliFileDir+lumi22Config, r3LumiTag );
    ReadLumiFile( caliFileDir+lumi23Config, r3LumiTag );
    ReadLumiFile( caliFileDir+lumi24Config, r3LumiTag );
    ReadLumiFile( caliFileDir+lumi25Config, r3LumiTag );

    setLumiTag(r3LumiTag);
    SetupIovMap( r3TagTT, minttlumi );
    SetupIovMap( r3TagX, minxlumi );
    //setXIov("r3x");

}//<<>>void KUCMS_TimeCalibration::SetupEraIovMap()

void KUCMS_TimeCalibration::SetupIovMap( std::string tag, float minLumi ){

    std::cout << " - Read SetupIovMap : " << tag << std::endl;
	//bool db = true;
    bool db = false;

	bool newrange( true );
	int start = 0;
	int end = 0;
	int prev = 0;
    int curxiovend = 0;
	float lumisum = 0;
	int currentrun = 0;
	bool inHardJump = false;
	std::map<int,int> theIovMap;
	auto promptIovMap = iovMaps["prompt"];
	// go through runs by lumi and add up lumi > check for end conditions
	for( auto& lumirun : lumiRunMaps[curLumiTag] ){
		currentrun = lumirun.second.run;
		if( newrange ){ // if we are starting a new range - find curent iov period for this run and set end/start of periods
			prev = start;
			start = ( prev == 0 ) ? currentrun : theIovMap[prev] + 1;
			newrange = false;
			for( auto& xiov : promptIovMap ){ 
				if( currentrun >= xiov.first && currentrun <= xiov.second ){ 
                    curxiovend = xiov.second; 
					break; 
				}//<<>>if( start >= xiov.first && start <= xiov.second )
			}//<<>>for( auto& xiov : promptIovMap ) 
			inHardJump = false;
		}//<<>>if( newrange )
        bool yearend = false;
		int yendrun = 0;
		// ::check if current run has past the start of a hard IOV boundary while in a range that started before that boundary
		// the Hard IOV boundary must be the start run of an IOV ( the run listed here must always trigger an IOV jump ) 
		std::vector<int> tthardruns{ 284045, 306461, 325173, 357501, 359001, 359560, 359567, 359574, 359690, 359775, 359870, 359871, 359880, 359895, 359902, 360129, 360131, 360260, 360331, 362761, 373862, 387742, 398904 };
		std::vector<int> xhardruns{ 253984, 284045, 306461, 327239, 362761, 373862, 387742, 398904 };
		if( minLumi < 3.f ) { for( auto hardrun : tthardruns ){ if( currentrun >= hardrun && start < hardrun ){ yearend = true; yendrun = hardrun; } } }
		else { for( auto hardrun : xhardruns ){ if( currentrun >= hardrun && start < hardrun ){ yearend = true; yendrun = hardrun; } } }
		//if( yearend ) std::cout << " --- END OF YEAR : " << currentrun << " start " << start << std::endl; 
		if( ( currentrun > curxiovend ) or yearend ){ // tack small portion to previous and start new run period with this run
			bool hasPrevious = prev > 0;
			bool addToPrevious = (lumisum < ( minLumi*0.4 )) or yearend;
			if( hasPrevious and addToPrevious and not inHardJump ){
				if( db ) std::cout << " IOV JUMP !! Add to Previous  !! Lumi Sum : " <<  lumisum << " min lumi : " << minLumi << std::endl;
				end = yearend ? yendrun - 1  : curxiovend;
				if( end < start ) end = start;// ID-10T check
				theIovMap[prev] = end;
				start = yearend ? yendrun : currentrun;
				lumisum = lumirun.second.lumi;
				inHardJump = yearend;
			}//<<>>if( prev > 0 )
			else if( inHardJump and yearend ){
                if( db ) std::cout << " IOV JUMP !! inHardJump and yearend  !! Lumi Sum : " <<  lumisum << " min lumi : " << minLumi << std::endl;
                end = yendrun - 1;
                if( end < start ) end = start;// ID-10T check
                theIovMap[start] = end;
                start = yendrun;
                lumisum = lumirun.second.lumi;
			}//<<>>else if( inHardJump and yearend ) 
			else {
				if( db ) std::cout << " IOV JUMP !! Continue !! Lumi Sum : " <<  lumisum << " min lumi : " << minLumi << std::endl;
				lumisum += lumirun.second.lumi; 
            }//<<>>else { if( prev > 0 ){ 
            for( auto& xiov : promptIovMap ){ // change IOV range we are in -- if big and not yearend : carry old range forward, if small or yearend : start new range 
                if( currentrun >= xiov.first && currentrun <= xiov.second ){
                    curxiovend = xiov.second;
                    break;
                }//<<>>if( start >= xiov.first && start <= xiov.second )
            }//<<>>for( auto& xiov : promptIovMap )
		} else {
			lumisum += lumirun.second.lumi;
        }//<<>>if( lumirun.second.run > promptIovMap[curxiov] )
		//std::cout << " Lumi Sum : " <<  lumisum << " min lumi : " << minLumi << std::endl;
		if( lumisum > minLumi ){ 
			if( db ) std::cout << " JUMP !!!!! Lumi Sum : " <<  lumisum << " min lumi : " << minLumi << std::endl;
			end = currentrun;
			if( end < start ) end = start;
			theIovMap[start] = end;
			lumisum = 0; 
			newrange = true;
		}//<<>>if( lumisum > minTTLumi )
	}//<<>>for( auto& lumirun : lumiRunMap )
	if( not newrange ) theIovMap[prev] = currentrun;
	iovMaps[tag] = theIovMap;		

    std::cout << " - Read SetupIovMap Finished : " << tag << std::endl;

	//for( auto& iov : theIovMap ){ std::cout << "Iov map " << tag << " : " << iov.first << " " << iov.second << std::endl; }

}//<<>>void KUCMS_TimeCalibration::makeTTIovMap()

void KUCMS_TimeCalibration::ReadLumiFile( std::string lumifile, std::string tag ){

	std::cout << " - Read LumiFile : " << lumifile << std::endl;
    std::ifstream infile( lumifile, std::ios::in);
	if( not infile.is_open() ){ std::cout << " -- " << lumifile << " not opened. "  << std::endl; return; }
    std::string infilestr;
	while( std::getline( infile, infilestr ) ){

		if(!infilestr.empty() && infilestr[0] != '#'){

			std::stringstream ss(infilestr);
    		std::string runfill, date, time;
    		int run, fill, nls, ncms;
    		float delivered, recorded;

			if (!(ss >> runfill >> date >> time >> nls >> ncms >> delivered >> recorded)){
    			std::cout << "Parse error in line: " << infilestr << std::endl;
    			continue;
			}//<<>>if (!(ss >> runfill >> date >> time >> nls >> ncms >> delivered >> recorded))
			//ss >> runfill >> date >> time >> nls >> ncms >> delivered >> recorded;
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
	std::cout << " - Read LumiFile Finished : " << lumifile << std::endl;

}//<<>>void KUCMS_TimeCalibration::upLoadLumiFile( std::string lumifile )

void KUCMS_TimeCalibration::ReadCaliRunFile(){

	std::cout << " - Reading CaliRunFile " << caliRunConfig << std::endl;
    std::ifstream infile( caliRunConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliRunConfig << " not opened. "  << std::endl; return; }
    std::string caliMapName, tag;
	int srun, erun, lrun;
    float lumi, noise, stoch, stant;
	bool resmap, resparms;
    bool hgsc;
    while( infile >> tag >> caliMapName >> srun >> erun >> lrun >> lumi >> resmap >> resparms >> noise >> stoch >> stant >> hgsc ){

		std::cout << " -- " << tag << " " << caliMapName << " " << srun << " " << erun << " " << lrun << " " << lumi;
		std::cout << " " << resmap << " " << resparms << " " << noise << " " << stoch <<  " " << stant << std::endl;

		CaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lrun, lumi };
        CaliRunMapSet[tag][srun].isNew = false;
		CaliRunMapSet[tag][srun].has2DResMap = resmap;
        CaliRunMapSet[tag][srun].hasResParams = resparms;
        CaliRunMapSet[tag][srun].noise = noise;
        CaliRunMapSet[tag][srun].stoch = stoch;
        CaliRunMapSet[tag][srun].stant = stant;
        CaliRunMapSet[tag][srun].hasGID2 = hgsc;

    }//<<>>while (infile >>
	infile.close();

}//<<>>void ReadTimeCaliTagFile()

bool KUCMS_TimeCalibration::SaveCaliRunFile(){

    std::cout << " - Saving CaliRunFile " << std::endl;
    std::ofstream outfile( caliRunConfig, std::ios::out | std::ios::trunc );
	if( not outfile.is_open() ){
		std::cerr << "ERROR: " << caliRunConfig << " not opened for writing." << std::endl;
		return false;
	}//<<>>if( not outfile.is_open() )
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
            bool hgsc = calirunsct.second.hasGID2;

            //std::cout << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
            outfile << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << " " << resmap; 
			outfile << " " << resparms << " " << noise << " " << stoch << " " << stant << " " << hgsc << std::endl;

		}//<<>>for( auto& calirunsct : calirunmap[tag] )

	}//<<>>for( auto& calirunmap : CaliRunMapSet )
	bool success = outfile.good();
	outfile.close();
	if( not success || outfile.fail() ){
		std::cerr << "ERROR: failed while writing " << caliRunConfig << std::endl;
		return false;
	}//<<>>if( write failed )
	return true;

}//<<>>void ReadTimeCaliTagFile()

void KUCMS_TimeCalibration::ReadTTRunFile(){

    std::cout << " - Reading TTCaliRunFile " << caliTTConfig << std::endl;
    std::ifstream infile( caliTTConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliTTConfig << " not opened. "  << std::endl; return; }
    std::string caliMapName, tag;
    int srun, erun, lrun;
    float lumi;
	bool hgsc;
    while( infile >> tag >> caliMapName >> srun >> erun >> lrun >> lumi >> hgsc ){
        //std::cout << " -- " << tag << " " << caliMapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
        TTCaliRunMapSet[tag][srun] = { caliMapName, srun, erun, lrun, lumi };
        TTCaliRunMapSet[tag][srun].isNew = false;
        TTCaliRunMapSet[tag][srun].hasGID2 = hgsc;
    }//<<>>while (infile >>
	infile.close();

}//<<>>void KUCMS_TimeCalibration::ReadTTRunFile()

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

bool KUCMS_TimeCalibration::SaveTTRunFile(){
    
	std::cout << " - Saving TTCaliRunFile " << std::endl;
    std::ofstream outfile( caliTTConfig, std::ios::out | std::ios::trunc );
	if( not outfile.is_open() ){
		std::cerr << "ERROR: " << caliTTConfig << " not opened for writing." << std::endl;
		return false;
	}//<<>>if( not outfile.is_open() )
    for( auto& calirunmap : TTCaliRunMapSet ){
        
        std::string tag = calirunmap.first;
		for( auto& calirunsct : calirunmap.second ){

        	int srun = calirunsct.second.startRun;
        	int erun = calirunsct.second.endRun;
        	std::string mapName = calirunsct.second.histMapName;
        	float lumi = calirunsct.second.lumi;
			int lrun = calirunsct.second.lastRun;
			bool hgsc = calirunsct.second.hasGID2; 
			//std::cout << " -- " << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << std::endl;
        	outfile << tag << " " << mapName << " " << srun << " " << erun << " " << lrun << " " << lumi << " " << hgsc << std::endl;

		}//<<>>for( auto& calirunsct : calirunmap )

    }//<<>>for( auto& calirunmap : CaliRunMapSet )
	bool success = outfile.good();
	outfile.close();
	if( not success || outfile.fail() ){
		std::cerr << "ERROR: failed while writing " << caliTTConfig << std::endl;
		return false;
	}//<<>>if( write failed )
	return true;

}//<<>>void ReadTimeCaliTagFile()

void KUCMS_TimeCalibration::ReadResFile(){

    std::cout << " - Reading ResFile " << caliResConfig << std::endl;
    std::ifstream infile( caliResConfig, std::ios::in);
    if( not infile.is_open() ){ std::cout << " -- " << caliResConfig << " not opened. "  << std::endl; return; }
    std::string tag;
    float ebnoise, ebstoch, ebstant, ebg1he1res, ebg1he2res, ebg1he3res, ebhg1res, ebhg2res, ebhg3res;
    float eenoise, eestoch, eestant, eeg1he1res, eeg1he2res, eeg1he3res, eehg1res, eehg2res, eehg3res;
    while( infile >> tag >> ebnoise >> ebstoch >> ebstant >> ebg1he1res >> ebg1he2res >> ebg1he3res >> ebhg1res >> ebhg2res >> 
				ebhg3res >> eenoise >> eestoch >> eestant >> eeg1he1res >> eeg1he2res >> eeg1he3res >> eehg1res >> eehg2res >> eehg3res ){
        ResTagSet[tag] = { ebnoise, ebstoch, ebstant,ebg1he1res, ebg1he2res, ebg1he3res,ebhg1res, ebhg2res, 
							ebhg3res,eenoise, eestoch, eestant,eeg1he1res, eeg1he2res, eeg1he3res,eehg1res, eehg2res, eehg3res };
		std::cout << " -- restag : " << tag << " " << ebnoise << " " << ebstoch << " " << ebstant << " " << ebg1he1res;
        std::cout << " " << ebg1he2res << " " << ebg1he3res << " " << ebhg1res << " " << ebhg2res << " " << ebhg3res; 
        std::cout << " " << eenoise << " " << eestoch << " " << eestant << " " << eeg1he1res << " " << eeg1he2res; 
        std::cout << " " << eeg1he3res << " " << eehg1res << " " << eehg2res << " " << eehg3res << " " << std::endl;
    }//<<>>while (infile >>
    infile.close();

}//<<>>void KUCMS_TimeCalibration::ReadSmearFile()

void KUCMS_TimeCalibration::SetResParamters( std::string tag, float ebnoise, float ebstoch, float ebstant, float ebg1he1res, 
	float ebg1he2res, float ebg1he3res, float ebhg1res, float ebhg2res, float ebhg3res, float eenoise, float eestoch,
    float eestant, float eeg1he1res, float eeg1he2res, float eeg1he3res, float eehg1res, float eehg2res, float eehg3res ){

	ResTagSet[tag] = { ebnoise, ebstoch, ebstant, ebg1he1res, ebg1he2res, ebg1he3res, ebhg1res, ebhg2res, ebhg3res,
		eenoise, eestoch, eestant, eeg1he1res, eeg1he2res, eeg1he3res, eehg1res, eehg2res, eehg3res };

}//<<>>void KUCMS_TimeCalibration::SetResParamters

void KUCMS_TimeCalibration::SaveResFile(){

    std::cout << " - Saving ResFile " << std::endl;
    std::ofstream outfile( caliResConfig, std::ios::out | std::ios::trunc );
    if( not outfile.is_open() ){ std::cout << "Failed to open file !!! " << std::endl; return; }
    for( auto& restag : ResTagSet ){

        std::string tag = restag.first;
		float ebnoise = restag.second.ebnoise;
        float ebstoch = restag.second.ebstoch;
        float ebstant = restag.second.ebstant;
        float ebg1he1res = restag.second.ebg1he1res;
        float ebg1he2res = restag.second.ebg1he2res;
        float ebg1he3res = restag.second.ebg1he3res;
        float ebhg1res = restag.second.ebhg1res;
        float ebhg2res = restag.second.ebhg2res;
        float ebhg3res = restag.second.ebhg3res;
        float eenoise = restag.second.eenoise;
        float eestoch = restag.second.eestoch;
        float eestant = restag.second.eestant;
        float eeg1he1res = restag.second.eeg1he1res;
        float eeg1he2res = restag.second.eeg1he2res;
        float eeg1he3res = restag.second.eeg1he3res;
        float eehg1res = restag.second.eehg1res;
        float eehg2res = restag.second.eehg2res;
        float eehg3res = restag.second.eehg3res;

        outfile << tag << " " << ebnoise << " " << ebstoch << " " << ebstant << " " << 
					ebg1he1res << " " << ebg1he2res << " " << ebg1he3res << " " <<
					ebhg1res << " " << ebhg2res << " " << ebhg3res << " " <<
					eenoise << " " << eestoch << " " << eestant << " " <<
                    eeg1he1res << " " << eeg1he2res << " " << eeg1he3res << " " <<
                    eehg1res << " " << eehg2res << " " << eehg3res << " " << std::endl;

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
//     ebeemap = {"_EB","_EE"};
// std::string tfilename = calirunmap.first + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
//            //calirunsct.second.histMapName = tfilename;
//            std::string sfilename = tfilename + "SumMap";
//            std::string s2filename = tfilename + "Sum2Map";
//            std::string ofilename = tfilename + "OccMap";
//            std::string mfilename = tfilename + "MeanMap";
//            std::string efilename = tfilename + "ErrMap";

bool KUCMS_TimeCalibration::LoadCaliHists( bool stayOpen, bool makeNew ){

	//bool oldstyle = false;

    std::cout << " - Loading CaliHists " << std::endl;
	//caliHistFileNames  caliFileDir

    //std::cout << " - opening caliHist Files " << std::endl;
    //caliTFileName = caliFileDir + "caliHistsTFile.root";// name configurable ?

	std::string califilepath;
	for( auto filename : caliHistFileNames ){
		if( useEosRootCali && not stayOpen ) califilepath = eosCaliPath + filename.second;
	 	else califilepath = caliFileDir + filename.second; 
		std::cout << " -- Cali file path : " << califilepath << " " << filename.first << std::endl;
        //if( not std::filesystem::exists(califilepath) ) continue; 
		if( makeNew ) caliTFile[filename.first] = TFile::Open( califilepath.c_str(), "RECREATE" );
		if( stayOpen && not makeNew ) caliTFile[filename.first] = TFile::Open( califilepath.c_str(), "UPDATE" );
		if( not stayOpen && not makeNew ) caliTFile[filename.first] = TFile::Open( califilepath.c_str(), "READ" );
		if( not caliTFile[filename.first] || caliTFile[filename.first]->IsZombie() || not caliTFile[filename.first]->IsOpen() ){
			std::cerr << "ERROR: failed to open ROOT file: " << califilepath << std::endl;
		}//<<>>if( caliTFile[filename.first]->IsZombie() or !caliTFile[filename.first] )
	}//<<>>for( auto filename : caliHistFileNames )

    std::cout << " - loading TT hists from caliHist Files " << std::endl;
	for( auto& calirunmap : TTCaliRunMapSet ){
		//std::string prefix = oldstyle ? calirunmap.first.substr(0,2) : getPrefix(calirunmap.first);
        std::string prefix = getPrefix(calirunmap.first);
		for( auto& calirunsct : calirunmap.second ){
			std::string ttfilename( calirunsct.second.histMapName );
			auto file = caliTFile.find(prefix);
			if( file == caliTFile.end() || not file->second || file->second->IsZombie() || not file->second->IsOpen() ){
				std::cerr << "ERROR: cannot load TT calibration histograms for prefix " << prefix << std::endl;
				continue;
			}//<<>>if( input file is unavailable )
			file->second->cd();
			for( auto mapname : maptypes ){
			for( auto ebeehist : ebeemap ){
				std::string filename = ttfilename + ebeehist + mapname;
				//std::cout << " -- loading : " << filename << " ( " << prefix << " ) ";
				if( ttfilename != "none" && CaliHists.find(filename) == CaliHists.end() ){
					TH2F* hist = (TH2F*)file->second->Get(filename.c_str());
					if( hist != NULL) hist->SetDirectory(nullptr);
                    else { std::cout << " Not Found " << filename << std::endl; }
                    //std::cout << " Found " << hist << std::endl;
					if( hist ) CaliHists[filename] = { hist, filename, false, false }; // histfile histname isnew isreshist
				}//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
				else { std::cout << " Duplicate " << filename << std::endl; }
			}//<<>>for( auto ebeehist : ebeemap )
			}//<<>>for( auto mapname : mtype )
		}//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

    std::cout << " - loading X hists from caliHists Files " << std::endl;
    for( auto& calirunmap : CaliRunMapSet ){
		//std::string prefix = oldstyle ? calirunmap.first.substr(0,2) : getPrefix(calirunmap.first);
        std::string prefix = getPrefix(calirunmap.first);
		for( auto& calirunsct : calirunmap.second ){
        	std::string xtfilename( calirunsct.second.histMapName );
			auto file = caliTFile.find(prefix);
			if( file == caliTFile.end() || not file->second || file->second->IsZombie() || not file->second->IsOpen() ){
				std::cerr << "ERROR: cannot load crystal calibration histograms for prefix " << prefix << std::endl;
				continue;
			}//<<>>if( input file is unavailable )
			file->second->cd();
            for( auto mapname : maptypes ){
            for( auto ebeehist : ebeemap ){
                std::string filename = xtfilename + ebeehist + mapname;
                //std::cout << " -- loading : " << filename << " ( " << xtfilename << " ) ";
                if( xtfilename != "none" && CaliHists.find(filename) == CaliHists.end() ){
					TH2F* hist = (TH2F*)file->second->Get(filename.c_str());
                    if( hist != NULL) hist->SetDirectory(nullptr); 
					else { std::cout << " Not Found " << filename << std::endl; }
                    //std::cout << " Found " << hist << std::endl;
                    if( hist ) CaliHists[filename] = { hist, filename, false, false }; // histfile histname isnew isreshist
                }//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
                else { std::cout << " Duplicate " << filename << std::endl; }
			}//<<>>for( auto ebeehist : ebeemap )
            }//<<>>for( auto mapname : mtype )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

	bool success = makeCaliMaps();

	if( not stayOpen ){

    	for( auto& calihist : CaliHists ){
        	if( calihist.second.h2f != NULL ){ delete calihist.second.h2f; }
        	//if( calihist.second.h1f != NULL ) delete calihist.second.h1f;
    	}//<<>>for( auto& calimapsct : TTCaliMaps )	
		CaliHists.clear();
        //caliTFile->Close();
        //caliR3TFile->Close();
		for( auto califile : caliTFile ){ if( califile.second ) califile.second->Close(); }

	}//<<>>if( not stayOpen )
	caliMapsValid = success;
	return success;

}//<<>>void KUCMS_TimeCalibration::SetupCaliHists()

//    TH2F* h2f;
//    std::string histName;
//    bool isNew;

// for updating and/or adding new info to califiles - do last - saving closes TCalifile and clears CaliHists
bool KUCMS_TimeCalibration::SaveCaliHists(){

    std::cout << " - Saving CaliHists " << std::endl;
	bool success = true;
	//caliTFile->cd();
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
			std::string prefix = getPrefix(calihist.second.h2f->GetName()); //.substr(0,2);
            //std::string thename = calihist.second.h2f->GetName();
			//if( thename.rfind("r2_ul18_mc",0) == 0 ){ prefix = "r2_ul18_mc"; std::cout << " found r2_ul18_mc" << std::endl; }
            //std::cout << " hf2 save prefix " << prefix << std::endl;
            //if( prefix == "r2_" ) caliTFile->cd();
            //else if( prefix == "r3_" ) caliR3TFile->cd();
            //else continue;
			auto file = caliTFile.find(prefix);
			if( file == caliTFile.end() || not file->second || file->second->IsZombie() || not file->second->IsOpen() ){
				std::cerr << "ERROR: no writable calibration ROOT file for prefix " << prefix << std::endl;
				success = false;
			} else if( not file->second->cd() ){
				std::cerr << "ERROR: could not select calibration ROOT file for prefix " << prefix << std::endl;
				success = false;
			} else if( notResHist && calihist.second.h2f->Write( calihist.second.h2f->GetName(), TObject::kOverwrite ) <= 0 ){
				std::cerr << "ERROR: failed to write histogram " << calihist.second.h2f->GetName() << std::endl;
				success = false;
			}//<<>>if( output file is valid )
            //if( notResHist ) calihist.second.h2f->Write();
			//std::cout << " -- deleting : " << calihist.second.h2f << std::endl;	
		}//<<>>if( calihist.second.h2f != NULL )
        if( calihist.second.h1f != NULL ){
			//if( notNew ) caliTFile->Delete( command.c_str() );
            std::string prefix = getPrefix(calihist.second.h1f->GetName()); //.substr(0,2);
            //std::string thename = calihist.second.h1f->GetName();
            //if( thename.rfind("r2_ul18_mc",0) == 0 ){ prefix = "r2_ul18_mc"; std::cout << " found r2_ul18_mc" << std::endl; }
            //std::cout << " hf1 save prefix " << prefix << std::endl;
            //if( prefix == "r2_" ) caliTFile->cd();
            //else if( prefix == "r3_" ) caliR3TFile->cd();
            //else continue;
			auto file = caliTFile.find(prefix);
			if( file == caliTFile.end() || not file->second || file->second->IsZombie() || not file->second->IsOpen() ){
				std::cerr << "ERROR: no writable calibration ROOT file for prefix " << prefix << std::endl;
				success = false;
			} else if( not file->second->cd() ){
				std::cerr << "ERROR: could not select calibration ROOT file for prefix " << prefix << std::endl;
				success = false;
			} else if( calihist.second.h1f->Write( calihist.second.h1f->GetName(), TObject::kOverwrite ) <= 0 ){
				std::cerr << "ERROR: failed to write histogram " << calihist.second.h1f->GetName() << std::endl;
				success = false;
			}//<<>>if( output file is valid )
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
    //caliTFile->Close();
	//caliR3TFile->Close();
	for( auto califile : caliTFile ){
		if( califile.second ){
			if( califile.second->TestBit(TFile::kWriteError) ) success = false;
			califile.second->Close();
			if( califile.second->TestBit(TFile::kWriteError) ) success = false;
		}//<<>>if( califile.second )
	}//<<>>for( auto califile : caliTFile )
	return success;

}///<<>>void KUCMS_TimeCalibration::SaveCaliMaps()

//std::string ttfilename = ttHistMapName+name+std::to_string(run);
//      MeanMap + pd/camp/tag name + start run  + TT/X ?
//      ErrMap + pd/camp/tag name + start run
//      SumMap + pd/camp/tag name + start run
//      Sum2Map + pd/camp/tag name + start run
//      OccMap + pd/camp/tag name + start run
// maptypes = {"_MeanMap","_ErrMap","_SumMap","_Sum2Map","_OccMap"};
//     ebeemap = {"_EB","_EE"};
//    std::map<uInt,sumCnt> sumCntMap;
//    std::map<uInt,float> meanMap;
//    std::map<uInt,float> errMap;
// DetIDStruct idinfo = DetIDMap[];


void KUCMS_TimeCalibration::makeCaliHists(){

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

	// do not associate with a TFile at this step
    //std::cout << "Recreate TFile " << std::endl;
    //caliTFile->cd();
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
            std::string sfilename_eb = tfilename + ebeemap[0] + maptypes[2];
            std::string s2filename_eb = tfilename + ebeemap[0] + maptypes[3];
            std::string ofilename_eb = tfilename + ebeemap[0] + maptypes[4];
            std::string mfilename_eb = tfilename + ebeemap[0] + maptypes[0];
            std::string efilename_eb = tfilename + ebeemap[0] + maptypes[1];
            std::string mdfilename_eb = tfilename + ebeemap[0] + "_MeanDist";

            std::string sfilename_ee = tfilename + ebeemap[1] + maptypes[2];
            std::string s2filename_ee = tfilename + ebeemap[1] + maptypes[3];
            std::string ofilename_ee = tfilename + ebeemap[1] + maptypes[4];
            std::string mfilename_ee = tfilename + ebeemap[1] + maptypes[0];
            std::string efilename_ee = tfilename + ebeemap[1] + maptypes[1];
            std::string mdfilename_ee = tfilename + ebeemap[1] + "_MeanDist";

            for( auto mapname : maptypes ){
			for( auto side : ebeemap ){
                std::string filename = tfilename + side + mapname;
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
				TH2F* hist(NULL);
				if( side == ebeemap[0] ) hist = new TH2F(filename.c_str(),filename.c_str(),ebx_ranges[0],ebx_ranges[1],ebx_ranges[2],ebx_ranges[3],ebx_ranges[4],ebx_ranges[5]);
				else hist = new TH2F(filename.c_str(),filename.c_str(),epmx_ranges[0],epmx_ranges[1],epmx_ranges[2],epmx_ranges[3],epmx_ranges[4],epmx_ranges[5]);
				if( hist != NULL) hist->SetDirectory(nullptr);
				CaliHists[filename] = { hist, filename, true, false }; // histfile histname isnew isopen lastrun
				//}//<<>>if( calihist != CaliHists.end() )
			}//<<>>for( auto side : ebeemap )
            }//<<>>for( auto mapname : mtype )

			for( auto& mdfilename : { mdfilename_eb, mdfilename_ee } ){ 
			if( CaliHists.find(mdfilename) != CaliHists.end() ){
                if( CaliHists[mdfilename].h1f != NULL ){
					delete CaliHists[mdfilename].h1f;
					CaliHists.erase(mdfilename);
                }//<<>>if( CaliHists[filename].h1f != NULL ) 
				//TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
				//CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
			}//<<>>if( CaliHists.find(mdfilename) == CaliHists.end() )
			TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2); 
			if( hist1d != NULL) hist1d->SetDirectory(nullptr);
			CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
			}//<<>>for( auto& mdfilename : {std::ref(mdfilename_eb), std::ref(mdfilename_ee) )

			//std::cout << " -- filling hists for " << tfilename << std::endl;
			for( auto& entry : calirunsct.second.sumCntMap ){
				//if( entry.first < 10000 ) std::cout << " -- HG id : " << entry.first << " with " << entry.second.cnt;
				if( DetIDMap.find(entry.first) == DetIDMap.end() ){ 
					std::cout << " -- " << entry.first << " not found in DetIDMap !!!!" << std::endl;
				}//<<>>if( DetIDMap.find(entry.first) == DetIDMap.end() )
				kucms_DetIDStruct idinfo = DetIDMap[entry.first];
                //if( entry.first < 10000 ) std::cout << " ecal " << idinfo.ecal << " i1 " << idinfo.i1 << " i2 " << idinfo.i2 << std::endl;
				if( idinfo.ecal == 0 ){
					CaliHists[sfilename_eb].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.sum );
                	CaliHists[s2filename_eb].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.sumsqr );
                	CaliHists[ofilename_eb].h2f->Fill( idinfo.i2, idinfo.i1, entry.second.cnt );
				} else {
					int mi2 = ( idinfo.ecal == 1 ) ? -1*idinfo.i2 : idinfo.i2;
                	CaliHists[sfilename_ee].h2f->Fill( mi2, idinfo.i1, entry.second.sum );
                	CaliHists[s2filename_ee].h2f->Fill( mi2, idinfo.i1, entry.second.sumsqr );
                	CaliHists[ofilename_ee].h2f->Fill( mi2, idinfo.i1, entry.second.cnt );
				}//<<>>if( idinfo.ecal == 0 )
			}//<<>>for( auto& entry : calirunsct.second.sumCntMap )
            for( auto& entry : calirunsct.second.meanMap ){
                kucms_DetIDStruct idinfo = DetIDMap[entry.first];
                if( idinfo.ecal == 0 ){
                	CaliHists[mfilename_eb].h2f->Fill( idinfo.i2, idinfo.i1, entry.second );
                	CaliHists[mdfilename_eb].h1f->Fill( entry.second );
                } else {
					int mi2 = ( idinfo.ecal == 1 ) ? -1*idinfo.i2 : idinfo.i2;
                	CaliHists[mfilename_ee].h2f->Fill( mi2, idinfo.i1, entry.second );
                	CaliHists[mdfilename_ee].h1f->Fill( entry.second );
                }//<<>>if( idinfo.ecal == 0 )
			}//<<>>for( auto& entry : calirunsct.second.meanMap )
            for( auto& entry : calirunsct.second.errMap ){
                kucms_DetIDStruct idinfo = DetIDMap[entry.first];
				int mi2 = ( idinfo.ecal == 1 ) ? -1*idinfo.i2 : idinfo.i2;
				if( idinfo.ecal == 0 ) CaliHists[efilename_eb].h2f->Fill( idinfo.i2, idinfo.i1, entry.second );
				else CaliHists[efilename_ee].h2f->Fill( mi2, idinfo.i1, entry.second );
            }//<<>>for( auto& entry : calirunsct.second.errMap )

        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )


    std::cout << "Making CaliHists from TTCaliRunMapSet" << std::endl;
    for( auto& calirunmap : TTCaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
            std::string tfilename( calirunsct.second.histMapName );
            //std::string tfilename = calirunmap.first + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
            //calirunsct.second.histMapName = tfilename;

            std::string sfilename_eb = tfilename + ebeemap[0] + maptypes[2];
            std::string s2filename_eb = tfilename + ebeemap[0] + maptypes[3];
            std::string ofilename_eb = tfilename + ebeemap[0] + maptypes[4];
            std::string mfilename_eb = tfilename + ebeemap[0] + maptypes[0];
            std::string efilename_eb = tfilename + ebeemap[0] + maptypes[1];
            std::string mdfilename_eb = tfilename + ebeemap[0] + "_MeanDist";

            std::string sfilename_ee = tfilename + ebeemap[1] + maptypes[2];
            std::string s2filename_ee = tfilename + ebeemap[1] + maptypes[3];
            std::string ofilename_ee = tfilename + ebeemap[1] + maptypes[4];
            std::string mfilename_ee = tfilename + ebeemap[1] + maptypes[0];
            std::string efilename_ee = tfilename + ebeemap[1] + maptypes[1];
            std::string mdfilename_ee = tfilename + ebeemap[1] + "_MeanDist";

            for( auto mapname : maptypes ){
            for( auto side : ebeemap ){
                std::string filename = tfilename + side + mapname;
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
                TH2F* hist(NULL);
                if( side == ebeemap[0] ) hist = new TH2F(filename.c_str(),filename.c_str(),ebtt_ranges[0],ebtt_ranges[1],ebtt_ranges[2],ebtt_ranges[3],ebtt_ranges[4],ebtt_ranges[5]);
                else hist = new TH2F(filename.c_str(),filename.c_str(),epmtt_ranges[0],epmtt_ranges[1],epmtt_ranges[2],epmtt_ranges[3],epmtt_ranges[4],epmtt_ranges[5]);
				if( hist != NULL) hist->SetDirectory(nullptr);
                CaliHists[filename] = { hist, filename, true, false }; // histfile histname isne
                //}//<<>>if( calihist != CaliHists.end() )
            }//<<>>for( auto side : ebeemap )
            }//<<>>for( auto mapname : mtype )

            for( auto& mdfilename : { mdfilename_eb, mdfilename_ee } ){
            if( CaliHists.find(mdfilename) != CaliHists.end() ){
                if( CaliHists[mdfilename].h1f != NULL ){
                    delete CaliHists[mdfilename].h1f;
                    CaliHists.erase(mdfilename);
                }//<<>>if( CaliHists[filename].h1f != NULL ) 
                //TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
                //CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
            }//<<>>if( CaliHists.find(mdfilename) == CaliHists.end() )
            TH1F* hist1d = new TH1F(mdfilename.c_str(),mdfilename.c_str(),400,-2,2);
			if( hist1d != NULL) hist1d->SetDirectory(nullptr);
            CaliHists[mdfilename] = { hist1d, mdfilename, false, false };
            }//<<>>for( auto& mdfilename : {std::ref(mdfilename_eb), std::ref(mdfilename_ee) )

            //std::cout << " -- filling hists for " << tfilename << std::endl;
            for( auto& entry : calirunsct.second.sumCntMap ){
                auto idinfo = getTTInfo( entry.first );
                if( entry.first < 4000 ){
                    CaliHists[sfilename_eb].h2f->Fill( idinfo.second, idinfo.first, entry.second.sum );
                    CaliHists[s2filename_eb].h2f->Fill( idinfo.second, idinfo.first, entry.second.sumsqr );
                    CaliHists[ofilename_eb].h2f->Fill( idinfo.second, idinfo.first, entry.second.cnt );
                } else {
                    CaliHists[sfilename_ee].h2f->Fill( idinfo.second, idinfo.first, entry.second.sum );
                    CaliHists[s2filename_ee].h2f->Fill( idinfo.second, idinfo.first, entry.second.sumsqr );
                    CaliHists[ofilename_ee].h2f->Fill( idinfo.second, idinfo.first, entry.second.cnt );
                }//<<>>if( idinfo.ecal == 0 )
            }//<<>>for( auto& entry : calirunsct.second.sumCntMap )
            for( auto& entry : calirunsct.second.meanMap ){
                auto idinfo = getTTInfo( entry.first );
                if( entry.first < 4000  ){
                    CaliHists[mfilename_eb].h2f->Fill( idinfo.second, idinfo.first, entry.second );
                    CaliHists[mdfilename_eb].h1f->Fill( entry.second );
                } else {
                    CaliHists[mfilename_ee].h2f->Fill( idinfo.second, idinfo.first, entry.second );
                    CaliHists[mdfilename_ee].h1f->Fill( entry.second );
                }//<<>>if( idinfo.ecal == 0 )
            }//<<>>for( auto& entry : calirunsct.second.meanMap )
            for( auto& entry : calirunsct.second.errMap ){
                auto idinfo = getTTInfo( entry.first );
                if( entry.first < 4000 ) CaliHists[efilename_eb].h2f->Fill( idinfo.second, idinfo.first, entry.second );
                else CaliHists[efilename_ee].h2f->Fill( idinfo.second, idinfo.first, entry.second );
            }//<<>>for( auto& entry : calirunsct.second.errMap )

        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )


}//<<>>void KUCMS_TimeCalibration::makeCaliHists()

bool KUCMS_TimeCalibration::makeCaliMaps(){

    std::cout << " - Making CaliMaps " << std::endl;
	bool valid = true;
	auto findConfiguredDetId = [&]( int first, int second, int ecal, uInt& detid ){
		auto firstCoord = InvDetIDMap.find(first);
		if( firstCoord == InvDetIDMap.end() ) return false;
		auto secondCoord = firstCoord->second.find(second);
		if( secondCoord == firstCoord->second.end() ) return false;
		auto ecalCoord = secondCoord->second.find(ecal);
		if( ecalCoord == secondCoord->second.end() ) return false;
		detid = ecalCoord->second;
		return DetIDMap.find(detid) != DetIDMap.end();
	};
	//TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),171,-85.5,85.5,360,0.5,360.5);
    for( auto& calirunmap : CaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
			if( calirunsct.second.isExternal ) continue;
            //std::string tfilename = calirunsct.second.histMapName + "_" + std::to_string( calirunsct.first ) + xtalHistMapName;
			std::string tfilename = calirunsct.second.histMapName;
			for( auto side : ebeemap ){

			//std::cout << std::endl;
			//std::cout << " -- opening : " << tfilename << std::endl;
            std::string sfilename = tfilename + side + "_SumMap";
            std::string s2filename = tfilename + side + "_Sum2Map";
            std::string ofilename = tfilename + side + "_OccMap";
            std::string mfilename = tfilename + side + "_MeanMap";
            std::string efilename = tfilename + side + "_ErrMap";
            //std::cout << " --- with : " << sfilename << std::endl;
            //if( CaliHists.find(sfilename) != CaliHists.end() ) std::cout << " --- Found " << std::endl;
			if( CaliHists.find(sfilename) != CaliHists.end() ){
			   bool isEB( side == ebeemap[0] );
			   int end_range_eta = isEB ? 172 : 204;
               int end_range_phi = isEB ? 365 : 101;
               //std::cout << " --- Found : " << tfilename << std::endl;
				for( int ieta = 1; ieta < end_range_eta; ieta++ ){
					for( int iphi = 1; iphi < end_range_phi; iphi++ ){
						//std::cout << " -- Processing : " << ieta << " " << iphi;
						int tti1 = ieta - 102;
                        int i1 = isEB ? ieta - 86 : std::abs(tti1);
						int ec = isEB ? 0 : ( tti1 < 0 ) ? 1 : 2; 
						if( i1 == 0 ) continue;
						float sum = CaliHists[sfilename].h2f->GetBinContent(ieta,iphi);
						float sum2 = CaliHists[s2filename].h2f->GetBinContent(ieta,iphi);
                		float occ = CaliHists[ofilename].h2f->GetBinContent(ieta,iphi);
						if( not std::isfinite(occ) || not std::isfinite(sum) || not std::isfinite(sum2) ||
						    occ < 0 || sum2 < 0 || ( occ == 0 && ( sum != 0 || sum2 != 0 ) ) ){
							std::cerr << "ERROR: invalid stored crystal moments in " << tfilename << side
							          << " at histogram bin (" << ieta << "," << iphi << "): sum=" << sum
							          << " sum2=" << sum2 << " occ=" << occ << std::endl;
							valid = false;
							continue;
						}//<<>>if( invalid stored moments )
						if( occ == 0 ) continue;
						uInt detid = 0;
						if( not findConfiguredDetId(iphi,i1,ec,detid) ){
							std::cerr << "ERROR: positive-occupancy crystal bin has no configured detector or gain ID in "
							          << tfilename << side << " at histogram bin (" << ieta << "," << iphi
							          << "), occ=" << occ << std::endl;
							valid = false;
							continue;
						}//<<>>if( detector ID missing )
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

			}//<<>>for( auto side : ebeemap )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

    std::cout << " - Making TTCaliMaps " << std::endl;
	//TH2F* hist = new TH2F(filename.c_str(),filename.c_str(),34,0,34,72,0,72);
    for( auto& calirunmap : TTCaliRunMapSet ){
        for( auto& calirunsct : calirunmap.second ){
            //std::string tfilename = calirunsct.second.histMapName + "_" + std::to_string( calirunsct.first ) + ttHistMapName;
            std::string tfilename = calirunsct.second.histMapName;
            for( auto side : ebeemap ){

            //std::cout << std::endl;
            //std::cout << " -- opening : " << tfilename << std::endl;
            std::string sfilename = tfilename + side + "_SumMap";
            std::string s2filename = tfilename + side + "_Sum2Map";
            std::string ofilename = tfilename + side + "_OccMap";
            std::string mfilename = tfilename + side + "_MeanMap";
            std::string efilename = tfilename + side + "_ErrMap";
            //std::cout << " --- with : " << sfilename << std::endl;
            //if( CaliHists.find(sfilename) != CaliHists.end() ) std::cout << " --- Found " << std::endl;
            if( CaliHists.find(sfilename) != CaliHists.end() ){
               	bool isEB( side == ebeemap[0] );
               	int end_range_eta = isEB ? 36 : 24;
               	int end_range_phi = isEB ? 75 : 35;
				//std::cout << " --- found : " << tfilename << std::endl;
                for( int ieta = 1; ieta < end_range_eta; ieta++ ){
                    for( int iphi = 1; iphi < end_range_phi; iphi++ ){
                        //std::cout << " -- Processing : " << ieta << " " << iphi;
                        int i1 = isEB ? ieta - 18 : ieta - 12;
                        if( i1 == 0 ) continue;
						uInt detid = getInvTTId( iphi, i1, isEB );
                        float sum = CaliHists[sfilename].h2f->GetBinContent(ieta,iphi);
                        float sum2 = CaliHists[s2filename].h2f->GetBinContent(ieta,iphi);
                        float occ = CaliHists[ofilename].h2f->GetBinContent(ieta,iphi);
						if( not std::isfinite(occ) || not std::isfinite(sum) || not std::isfinite(sum2) ||
						    occ < 0 || sum2 < 0 || ( occ == 0 && ( sum != 0 || sum2 != 0 ) ) ){
							std::cerr << "ERROR: invalid stored TT moments in " << tfilename << side
							          << " at histogram bin (" << ieta << "," << iphi << "): sum=" << sum
							          << " sum2=" << sum2 << " occ=" << occ << std::endl;
							valid = false;
							continue;
						}//<<>>if( invalid stored moments )
                        if( occ == 0 ) continue;
						if( detid == 0 ){
							std::cerr << "ERROR: positive-occupancy TT bin has no configured tower ID in "
							          << tfilename << side << " at histogram bin (" << ieta << "," << iphi
							          << "), occ=" << occ << std::endl;
							valid = false;
							continue;
						}//<<>>if( TT ID missing )
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

            }//<<>>for( auto side : ebeemap )
        }//<<>>for( auto& calirunsct : calirunmap )
    }//<<>>for( auto& calirunmap : CaliRunMapSet )

	return valid;

}//<<>>bool KUCMS_TimeCalibration::makeCaliMaps()

//  LoadExtCali is EB only right now -- update !!!!
void KUCMS_TimeCalibration::LoadExtCali( std::string calihist, std::string mapname, std::string tag, int startr, int endr ){

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
    	for( int iphi = 1; iphi < 365; iphi++ ){
        	//std::cout << " -- Processing : " << ieta << " " << iphi;
            int i1 = ieta - 86;
            if( i1 == 0 ) continue;
            uInt detid = InvDetIDMap[iphi][i1][ECAL::EB];
         	float occ = ohist->GetBinContent(ieta,iphi);
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

}//<<>>void KUCMS_TimeCalibration::LoadExtCali( std::string calihist, std::string mapname, std::string tag )

uInt KUCMS_TimeCalibration::getTTId( uInt detId, int gainID ){

	uInt id = 9999999;
	if( detId < 10000 ){ std::cout << " -- getTTId with a HG ID number !!!!!!!!!!1" << std::endl; return detId; }
	auto det = DetIDMap.find(detId);
	if( det == DetIDMap.end() ){
		std::cerr << "ERROR: unknown ECAL detID in getTTId: " << detId << std::endl;
		return 0;
	}//<<>>if( detID missing )
	const auto& idinfo = det->second;
	if( idinfo.ecal == 0 ){
		if( gainID == 1 ){
			int ttphi = 1 + int( idinfo.i1 - 1 )/5;
			int tteta = 1 + int( std::abs( idinfo.i2 ) - 1 )/5;
			id = ( idinfo.i2 < 0 ) ? (ttphi+tteta*100)+2000 : (ttphi+tteta*100);
		} else { id = ( idinfo.i2 < 0 ) ? 2174 : 174; }
	} else {
		int side = idinfo.ecal;
		if( gainID == 1 ){
			int mod = idinfo.mod;
			id = 1000*( side + 3 ) + mod; 
		} else { id = 1000*( side + 3 ) + 330; }
	}//<<>>if( DetIDMap[detId].ecal == ECAL::EB ) else
	//std::cout << " getTTId " << detId << " phi " << DetIDMap[detId].i1 << " eta " << DetIDMap[detId].i2 << std::endl; 
    //std::cout << "  --> p " << ttphi << " e " << tteta << " = " << id <<std::endl;
	return id;

}//<<>>uInt KUCMS_TimeCalibration::getTTId( uInt detId )

std::pair<int,int> KUCMS_TimeCalibration::getTTInfo( uInt ttid ){

	int i1 = 0;
	int i2 = 0;
	if( ttid < 4000 ){
		int a = ( ttid > 2000 ) ? ttid - 2000 : ttid;
		int t = a/100;
		i2 = ( ttid > 2000 ) ? -1*t : t;
		i1 = a - t*100;
	} else {
		int s = ttid/1000;
		int a = ttid - (1000*s);
		int ah = a/10;
		int al = a - (10*ah) + 1;
		int pm = ( s == 4 ) ? -1 : 1;
		i1 = ah + 1;
		i2 = al*pm;
		if( ttid > 6000 ){ i1 = 200; i2 = 200; } 
		//std::cout << " -- " << i1 << " " << i2 << " " << ttid << std::endl;
	}//<<>>if( ttid < 4000 ) else
	return std::make_pair(i1,i2); 

}//<<>>std::pair<int,int> KUCMS_TimeCalibration::getTTInfo( uInt ttid )

uInt KUCMS_TimeCalibration::getInvTTId( int i1, int i2, bool isEB ){

	int tt = 0;
	if( isEB ){ 
		if( i2 < 0 ) tt = (i1+std::abs(i2)*100)+2000;
		else tt = i1+(i2*100);
	} else {
		int al = std::abs(i2) - 1;
		int ah = ( i1 - 1 )*10;
		int m = ( i2 < 0 ) ? 4000 : 5000;
        if( al > 9 ) tt = 9000;
		else tt = al + ah + m;
		//std::cout << " -- " << i1 << " " << i2 << " " << tt << " : " << ah << " " << al << " " << m << std::endl;
	}//<<>> if( isEB ) else 
	return tt;

}//<<>>std::pair<int,int> KUCMS_TimeCalibration::getTTInfo( uInt ttid )

kucms_DetIDStruct& KUCMS_TimeCalibration::getDetIdInfo( uInt rhid ){

	return DetIDMap[rhid];

}//<<>>DetIDStruct KUCMS_TimeCalibration::getDetIdInfo( uInt rhid )

uInt KUCMS_TimeCalibration::getDetIdInfo( int i1, int i2, int ecal ){

	return InvDetIDMap[i1][i2][ecal];

}//<<>>DetIDStruct KUCMS_TimeCalibration::getDetIdInfo( uInt rhid )

std::map<int,std::map<int,std::map<int,uInt>>> KUCMS_TimeCalibration::getInvDetIDMap(){

	return InvDetIDMap;

}//<<>>DetIDStruct KUCMS_TimeCalibration::getInvDetIDMap( )

std::map<UInt_t,kucms_DetIDStruct> KUCMS_TimeCalibration::getDetIDMap(){; // map of information by detid for each crystal

	return DetIDMap;

}//<<>>DetIDStruct KUCMS_TimeCalibration::getDetIdInfo( uInt rhid )

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

float KUCMS_TimeCalibration::getCalibration( uInt rhid, int run, const std::string& tag, int gainID ){

	if( rhid == 0 ) return -999.f; // valid rechit id check
	//if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
	//if( not isEB ){ std::cout << "XCalibration for EE is not supported." << std::endl; return 0.f; }
    //if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
    //if( DetIDMap[rhid].ecal != ECAL::EB ) return 0.f;
	if( tag != "None" ) externalCali = false;
	if( externalCali ){ return 0.f; std::cout << " Called getCalibration() while using external cali? " << std::endl; }
    double xtaltime = -1000.f;
	//////float hgcali = ( externalCali ) ? 0.f : ( gainID == 1 ) ? 0.f : getTTCali( rhid, run, tag, 2 );
	double hgcali = 0;
	uInt gsid = ( gainID > 1 ) ? getGSRhId( rhid, gainID ) : 0;
	if( gainID > 1 and ( gsid < 99 or gsid > 300 ) ) std::cout << " Bad GS Id with " << rhid << " gid " << gainID << std::endl;
	double ttcali = getTTCali( rhid, run, tag, 1 );
	auto tagSet = CaliRunMapSet.find(tag);
	if( tagSet != CaliRunMapSet.end() ){
		auto range = tagSet->second.upper_bound(run);
		if( range != tagSet->second.begin() ){
			--range;
			if( run <= range->second.endRun ){
				auto xtalMean = range->second.meanMap.find(rhid);
				xtaltime = xtalMean == range->second.meanMap.end() ? 0.0 : xtalMean->second;
				if( gsid > 0 ){
					auto hgMean = range->second.meanMap.find(gsid);
					if( hgMean != range->second.meanMap.end() ) hgcali = hgMean->second;
				}//<<>>if( high-gain ID )
			}//<<>>if( run in range )
		}//<<>>if( candidate range exists )
	}//<<>>if( tag exists )
	if( xtaltime == -1000.f ){ std::cout << "XCalibration period not found for run " << run << std::endl; return 0.f; }
	//std::cout << " -- For " << rhid << " in " << run << " gainID " << gainID; 
    //std::cout << " Found: " << xtaltime << " + " << ttcali << " + " << hgcali << std::endl;
	return xtaltime + ttcali + hgcali;
	//return ttcali;
	//if( gainID == 1 ) return xtaltime + ttcali + hgcali;
    //else return 0.f;

}//<<>>float KUCMS_TimeCalibration::getCalibration( std::string tag )

float KUCMS_TimeCalibration::getTTCali( uInt rhid, int run, const std::string& tag, int gainID ){

    //if( not isEB ){ std::cout << "TTCalibration for EE is not supported." << std::endl; return 0.f; }
    //if( DetIDMap[rhid].ecal != ECAL::EB ) return 0.f;
    gainID = 1;
    double ttcalitime = -1000.f;
    uInt ttid = getTTId( rhid, gainID );
	auto tagSet = TTCaliRunMapSet.find(tag);
	if( tagSet != TTCaliRunMapSet.end() ){
		auto range = tagSet->second.upper_bound(run);
		if( range != tagSet->second.begin() ){
			--range;
			if( run <= range->second.endRun ){
				auto mean = range->second.meanMap.find(ttid);
				ttcalitime = mean == range->second.meanMap.end() ? 0.0 : mean->second;
			}//<<>>if( run in range )
		}//<<>>if( candidate range exists )
	}//<<>>if( tag exists )
    //std::cout << "TTGetCali : " << run << " " << tag << " " << rhid << " " << ttid << " " << ttcalitime << std::endl;
    if( ttcalitime == -1000.f ){ std::cout << "TTCalibration period not found for run " << run << std::endl; return 0.f; }
    return ttcalitime;

}//<<>>float KUCMS_TimeCalibration::getCalibration( std::string tag )

//struct smearTagStruct {
//
//    float noise;
//    float stoch;
//    float stant;
//
//};//smearParameters
//std::map<std::string,smearTagStruct> SmearTagSet;
//  do mutiple combos 

float KUCMS_TimeCalibration::getSmearedTime( float mtime, float ampres, std::string stag  ){

	if( stag == "debug" ) std::cout << " -- smeairng : " << mtime << " with " << ampres << std::endl;


	return getRandom->Gaus( mtime, ampres );

/*
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
*/
	//return -99.f;

}//<<>>float KUCMS_TimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

float KUCMS_TimeCalibration::getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int crun, std::string ctag, std::string stag ){

/*
    float crhtime = rhtime - getCalibration( rhid, crun, ctag );
    float smrdCalibTime = getSmearedTime( crhtime, rhamp, stag );
    return smrdCalibTime;
*/

	return -99.f;

}//<<>>float KUCMS_TimeCalibration::getSmearedTime( std::string tag , float time, uInt rhid )

//**************  main calibrate &/or smearing function based on tag + info ********************************************************************

///  Returns time Variance -> return sigle crystal variance
//float KUCMS_TimeCalibration::getTimeVariance( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey, int mctype ){
/*
float KUCMS_TimeCalibration::getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey, int mctype ){

	if( amplitude == 0 ) return 100.f;
	bool isEB( DetIDMap[rechitID].ecal == ECAL::EB );
	float res = -1;
	if( mctype < 100 ){ // data

		if( dataSetKey == "r2_ul18" || dataSetKey == "r2_ul18_mc"   ){ 
			res = isEB ? sq2(26.23/amplitude) + sq2(1.72)/amplitude + 2*sq2(0.1391) : sq2(33.36/amplitude) + sq2(1.72)/amplitude + 2*sq2(0.2039); }
        else if( dataSetKey == "r2_ul17"   ){ 
			res = isEB ? sq2(25.98/amplitude ) + sq2(0.67)/amplitude + 2*sq2(0.0947) : sq2(42.17/amplitude ) + sq2(6.29)/amplitude + 2*sq2(0.2986); }
        else if( dataSetKey == "r2_ul16"   ){
            res = isEB ? sq2(25.53/amplitude ) + sq2(1.88)/amplitude + 2*sq2(0.1576) : sq2(31.53/amplitude ) + sq2(1.90)/amplitude + 2*sq2(0.1805); }

        else if( dataSetKey == "r2_eoy17"   ){
            res = isEB ? sq2(31.3/amplitude ) + sq2(2.9)/amplitude + 2*sq2(0.1608) : sq2(31.3/amplitude ) + sq2(2.9)/amplitude + 2*sq2(0.1608); }
		else { res = isEB ? sq2(26.0/amplitude) + sq2(0.6)/amplitude + 2*sq2(0.094) : sq2(42.2/amplitude) + sq2(1.2)/amplitude + 2*sq2(0.298); }

	}//<<>>if( mctype == 0 )
	if( res == -1 ){ std::cout << " -- Resolution for dataSetKey " << dataSetKey << " not found !!!!!!" << std::endl; return 1.f; }
	return res/2;
	//return res;

}//<<>>float getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey )
*/

float KUCMS_TimeCalibration::getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey, int mctype, int gianID ){

    if( amplitude == 0 ) return 100.f;
    bool isEB( DetIDMap[rechitID].ecal == ECAL::EB );
    float var = ( dataSetKey == "None" ) ? -1 : -99;
    if( mctype == 1 ){ // data or MC
		if( dataSetKey == "r2_ul18" ){ if( Evt_run > 320672 ) dataSetKey = "r2_ul18D"; else dataSetKey = "r2_ul18ABC"; }
		if( ResTagSet.find(dataSetKey) != ResTagSet.end() ){
			if( amplitude > 400 ){ 

				float res = 0;
				if( amplitude > 900 ){ 
					if( isEB ) res = ( gianID == 1 ) ? ResTagSet[dataSetKey].ebg1he3res : ResTagSet[dataSetKey].ebhg3res;
					else res = ( gianID == 1 ) ? ResTagSet[dataSetKey].eeg1he3res : ResTagSet[dataSetKey].eehg3res;
				} else if( amplitude > 600 ){ 
					if( isEB ) res = ( gianID == 1 ) ? ResTagSet[dataSetKey].ebg1he2res : ResTagSet[dataSetKey].ebhg2res;
					else res = ( gianID == 1 ) ? ResTagSet[dataSetKey].eeg1he2res : ResTagSet[dataSetKey].eehg2res;
				} else { 
					if( isEB ) res = ( gianID == 1 ) ? ResTagSet[dataSetKey].ebg1he1res : ResTagSet[dataSetKey].ebhg1res;
					else res = ( gianID == 1 ) ? ResTagSet[dataSetKey].eeg1he1res : ResTagSet[dataSetKey].eehg1res;
				}//<<>>if( amplitude > 900 )
				var = sq2( res );

			}//<<>>if( amplitude > 400 )
			if( amplitude <= 400 or var == 0 ){
			
				float noise = isEB ? ResTagSet[dataSetKey].ebnoise : ResTagSet[dataSetKey].eenoise;
            	float stoch = isEB ? ResTagSet[dataSetKey].ebstoch : ResTagSet[dataSetKey].eestoch;
            	float stant = isEB ? ResTagSet[dataSetKey].ebstant : ResTagSet[dataSetKey].eestant;
				var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);

			}//<<>>if( amplitude > 400 )
/* 
 			// test to see if photon wtsig in 24C-E was impacted by incorrect resolutions -  nope 
			if( dataSetKey == "r3_p24" and Evt_run > 379411 and Evt_run < 381943 ){
				std::cout << " -- Resolution swtich to r3_p25 for r3_p24 C-E !!!!!!" << std::endl;
				std::string key = "r3_p25";
				float noise = isEB ? ResTagSet[key].ebnoise : ResTagSet[key].eenoise;
            	float stoch = isEB ? ResTagSet[key].ebstoch : ResTagSet[key].eestoch;
            	float stant = isEB ? ResTagSet[key].ebstant : ResTagSet[key].eestant;
            	var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);
			}//<<>>if( dataSetKey == "run3_24p" and Evt_run > 380253 and Evt_run < 380947 )
*/
		}//<<>>if( ResTagSet.find(tag) != ResTagSet.end() )
    }//<<>>if( mctype == 0 )
    else if( mctype == 0 or mctype == 2 ){ // data or MC   --- resTag is the target resolution for smearing
		std::string themctag = ( resTag == "None" ) ? dataSetKey : resTag;
        if( resTag == "default"){ 
			if( dataSetKey.substr(1,2) == "r2" ){ themctag = "r2_ul18D"; }//<<>>if( dataSetKey.substr(1,2) == "r2" )
        	else if( dataSetKey.substr(1,2) == "r3" ){ themctag = "r2_ul18D"; }//<<>>if( dataSetKey.substr(1,2) == "r3" )
		}//<<>>if( resTag == "default")
        if( ResTagSet.find(themctag) != ResTagSet.end() ){
            float noise = isEB ? ResTagSet[themctag].ebnoise : ResTagSet[themctag].eenoise;
            float stoch = isEB ? ResTagSet[themctag].ebstoch : ResTagSet[themctag].eestoch;
            float stant = isEB ? ResTagSet[themctag].ebstant : ResTagSet[themctag].eestant;
            var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);
			//std::cout << " --- getting MC res for : " << themctag << " " << noise << " " << stoch;
			//std::cout << " " << stant << " " << var << std::endl;
        }//<<>>if( ResTagSet.find(tag) != ResTagSet.end() )
    }//<<>>if( mctype == 0 )
	else var = -9;
    if( var < -10 ) std::cout << " -- Resolution for ResTag : " << dataSetKey << " -- not found !!!!!!" << std::endl;
    if( var < -5 ) std::cout << " -- var < -5 ? MC Type set to : " << mctype << " -- not a valid type ? !!!!!!" << std::endl;
    if( var < 0 ){
		std::cout << " -- Setting resolution to default values !!!!! " << std::endl;
		std::string tag = "default";
        float noise = isEB ? ResTagSet[tag].ebnoise : ResTagSet[tag].eenoise;
        float stoch = isEB ? ResTagSet[tag].ebstoch : ResTagSet[tag].eestoch;
        float stant = isEB ? ResTagSet[tag].ebstant : ResTagSet[tag].eestant;
        var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);
		//var = isEB ? sq2(26.0/amplitude) + sq2(0.6)/amplitude + 2*sq2(0.094) : sq2(42.2/amplitude) + sq2(1.2)/amplitude + 2*sq2(0.298); 
	}//<<>>if( var == -1 )
    return var;

}//<<>>float getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, std::string dataSetKey )

float KUCMS_TimeCalibration::getCorrectedTime( float time, float amplitude, unsigned int rechitID, 
												unsigned int Evt_run, std::string dataSetKey, int mctype, int gainID ){

	//getTimeResoltuion( amplitude, rechitID, Evt_run, dataSetKey, mctype, gianID )
	// NOTE :: dataSetKey is the tctag in this function!!!!    named wring here - please correct
	double rtime = 0;
    //if( Evt_run > 388000 || Evt_run < 253984 ) return time; // within calibrated run range check

    if( mctype == 1 ){ // data
		
		std::string theKey = dataSetKey;
		if( theKey == "r2_ul18D" or theKey == "r2_ul18ABC") theKey = "r2_ul18";
        if( CaliRunMapSet.find(theKey) != CaliRunMapSet.end() ){ 
			//std::cout << " ----- cali : " << getCalibration(rechitID,Evt_run,dataSetKey,gainID) << std::endl;
			rtime = time - getCalibration(rechitID,Evt_run,dataSetKey,gainID);
			return rtime;
		} else { std::cout << "GetCorrectedTime for Data dataSetKey not found : " << dataSetKey << std::endl; return time; }    

    } else { // MC

		double ctime = 0;
		double tarvar = 0;
		double sorvar = 0;

		if( dataSetKey == "r2_ul18_mc" ) dataSetKey = "r2_ul18_sqmc";
		if( CaliRunMapSet.find(dataSetKey) != CaliRunMapSet.end() ) ctime = time - getCalibration(rechitID,Evt_run,dataSetKey,gainID);
        else{ ctime = time; }
		//std::cout << " Pre Check amplitude : " << amplitude << std::endl;
		//if( amplitude < 1 ) return ctime;
        //if( time == 0 ) return ctime;
		if( resTag == "None" ) return ctime; // --- resTag is the target resolution for smearing
		if( dataSetKey == "None" ) return ctime; // dataSetKey is the calibration used && source resolution for smearing
		//std::cout << " Post Check amplitude : " << amplitude << std::endl;

        std::string theKey = dataSetKey;
		if( resTag != "default" && ResTagSet.find(theKey) != ResTagSet.end() ){  // TARGET RESOLUTION

			tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, resTag, mctype, gainID );
			
		}//<<>>if( resTag != "default" )
		else if( resTag == "default" || dataSetKey.substr(1,2) == "r2" ){ // use smear tag?  -- check first for set tag ! 

            std::string tag= "r2_ul18D";
			tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

		}//<<>>if( dataSetKey.substr(1,2) == "r2" )
        else if( dataSetKey.substr(1,2) == "r3" ){ 

            std::string tag= "r2_ul18D";
            tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

		}//<<>>if( dataSetKey.substr(1,2) == "r3" )

		if( ResTagSet.find(dataSetKey) != ResTagSet.end() ){ // SOURCE RESOLTUION

			std::string tag = dataSetKey;
            if( dataSetKey == "r2_ul18_qcdmc" ) tag = "RunIISummer20UL18RECO";
            if( dataSetKey == "r2_ul18_sqmc" ) tag = "RunIISummer20UL18RECO";
            if( dataSetKey == "r3_mc" ) tag = "theRun3MC_campain";
			sorvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

		} else { // 2018 MC Resoltion values

			std::string tag = "default";
			sorvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

		}//<<>>if( ResTagSet.find(tag) != ResTagSet.end() )

		//0.00691415 1.60666 0.133913  smearing for 2017 MC -> Data ( original smearing values )
		double tresv = tarvar/2;
		double sresv = sorvar/2;
        double smvar = ( tresv > sresv ) ? tresv - sresv : 0.0;//      std::max( 0.0, tresv - sresv );
        double resolution = ( smvar > 0.0 ) ? std::sqrt(smvar) : 0.0;
        rtime = getRandom->Gaus( ctime, resolution );
    	//if( tresv - sresv < 0.0 ){ 
        if( false ){
			std::cout << "No smearing values set for this tag : " << dataSetKey << " time " << time << " c " << ctime << std::endl;
			std::cout << " smearing parameters: t " << tresv << " s " << sresv << " amp " << amplitude << std::endl;
        	std::cout << " getSmeared Res: " << resolution << " rtime " << rtime << std::endl;
		}//<<>>if( resolution <= 0 )
		return rtime;

    }//<<>>if( mctype == 0 )

	return rtime;

}//<<>>float KUCMS_TimeCalibration::getCorrectedTime

//*********************************************************************************************************************************************

// assume full runs are present in input files
// input file will give the name of a text file "list" of root files, the path to the root files, and the tag to use for those files 
// check if run range started
// check if run range complete 
//     std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>
//	- if started but not complete - is run needed ( lastrun )
//	- if not started - start new cali map

bool KUCMS_TimeCalibration::makeGain1CaliMapsEGR( std::string inputFileName, bool small ){

	std::cout << "Creating gain-1 TT and crystal calibration maps in one event pass" << std::endl;
	if( not caliMapsValid ){
		std::cerr << "ERROR: cannot process calibration input because stored calibration maps are invalid." << std::endl;
		return false;
	}//<<>>if( stored maps invalid )
	if( iovMaps.find(curTTIov) == iovMaps.end() ){
		std::cerr << "ERROR: TT IOV is not valid: " << curTTIov << std::endl;
		return false;
	}//<<>>if( invalid TT IOV )
	if( iovMaps.find(curXIov) == iovMaps.end() ){
		std::cerr << "ERROR: crystal IOV is not valid: " << curXIov << std::endl;
		return false;
	}//<<>>if( invalid crystal IOV )

	const std::string treename("tree/llpgtree");
	using RawAccumulator = std::unordered_map<uInt,kucms_sumCnt>;
	struct FastDetInfo { uInt ttId; bool isEB; };
	std::unordered_map<uInt,FastDetInfo> fastDetInfo;
	fastDetInfo.reserve(DetIDMap.size());
	for( const auto& det : DetIDMap ){
		const auto& info = det.second;
		uInt ttId = 0;
		if( info.ecal == ECAL::EB ){
			const int ttPhi = 1 + ( info.i1 - 1 )/5;
			const int ttEta = 1 + ( std::abs(info.i2) - 1 )/5;
			ttId = info.i2 < 0 ? (ttPhi + ttEta*100) + 2000 : ttPhi + ttEta*100;
		} else ttId = 1000*(info.ecal + 3) + info.mod;
		fastDetInfo.emplace(det.first,FastDetInfo{ttId,info.ecal == ECAL::EB});
	}//<<>>for( detector IDs )
	// TT and crystal IOV boundaries are independent. Keep raw crystal moments split by
	// both boundaries until the corresponding TT means have been calculated.
	std::map<std::string,std::map<int,std::map<int,RawAccumulator>>> rawCrystalByTTIov;
	std::map<std::string,std::map<int,int>> maxTTRun;
	std::map<std::string,std::map<int,int>> maxCrystalRun;
	std::map<std::string,int> maxRunSeenByTag;
	unsigned long long rejectedSentinelTimes = 0;
	unsigned long long rejectedNonFiniteTimes = 0;
	unsigned int inputSetsProcessed = 0;
	int doLastRun = 0;

	auto findIov = []( const std::map<int,int>& iovs, int run ) -> std::pair<int,int> {
		auto iov = iovs.upper_bound(run);
		if( iov == iovs.begin() ) return {0,0};
		--iov;
		if( run > iov->second ) return {0,0};
		return *iov;
	};
	auto findOrCreateRange = [&]( auto& mapSet, const std::string& tag, const std::map<int,int>& iovs,
	                              int run, const std::string& type ) -> kucms_CaliRunClass* {
		auto bounds = findIov(iovs,run);
		if( bounds.first == 0 ){
			std::cerr << "ERROR: no configured " << type << " IOV contains run " << run
			          << " for tag " << tag << std::endl;
			return nullptr;
		}//<<>>if( IOV not found )
		auto& runSet = mapSet[tag];
		auto range = runSet.find(bounds.first);
		if( range == runSet.end() ){
			std::string histName = tag + "_" + std::to_string(bounds.first) + "_" + type;
			range = runSet.try_emplace(bounds.first, histName, bounds.first, bounds.second, 0, 0).first;
		} else if( range->second.startRun != bounds.first || range->second.endRun != bounds.second ){
			std::cerr << "ERROR: existing " << type << " range disagrees with configured IOV for run "
			          << run << " and tag " << tag << std::endl;
			return nullptr;
		}//<<>>if( range missing or inconsistent )
		return &range->second;
	};

	std::ifstream infilelist(inputFileName);
	if( not infilelist.is_open() ){
		std::cerr << "ERROR: could not open input file list: " << inputFileName << std::endl;
		return false;
	}//<<>>if( input list failed )

	std::string infilestr;
	while( std::getline(infilelist,infilestr) ){
		std::cout << "Processing : " << infilestr;
		if( infilestr.empty() || infilestr[0] == '#' ){
			std::cout << " -> Aborted" << std::endl;
			continue;
		}//<<>>if( blank or comment )
		std::cout << std::endl;

		std::stringstream ss(infilestr);
		std::string subdir, matchstr, tag;
		int srun = 0, erun = 0, dolast = 0;
		if( not ( ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast ) || subdir.empty() ||
		    matchstr.empty() || tag.empty() || srun > erun || ( dolast != 0 && dolast != 1 ) ){
			std::cerr << "ERROR: malformed calibration input row: " << infilestr << std::endl;
			return false;
		}//<<>>if( malformed row )
		doLastRun = dolast;
		if( isCC && doUnCC ) tag += "unc";

		auto fInTree = new TChain(treename.c_str());
		std::string eosHost = "cmseos.fnal.gov";
		std::string scanDir = eosDir + inDir + subdir;
		std::cout << "Finding ROOT files under EOS directory: " << scanDir << std::endl;
		auto files = findEOSRootFiles(eosHost,scanDir,matchstr);
		std::cout << "Found " << files.size() << " candidate ROOT files." << std::endl;
		int nAdded = 0;
		int nBadFiles = 0;
		for( const auto& fileName : files ){
			TFile* testFile = TFile::Open(fileName.c_str(),"READ");
			if( not testFile || testFile->IsZombie() || not testFile->IsOpen() ){
				std::cerr << "ERROR: bad ROOT file, skipping: " << fileName << std::endl;
				if( testFile ){ testFile->Close(); delete testFile; }
				nBadFiles++;
				continue;
			}//<<>>if( bad ROOT file )
			TTree* testTree = nullptr;
			testFile->GetObject(treename.c_str(),testTree);
			if( not testTree ){
				std::cerr << "ERROR: missing tree " << treename << ", skipping: " << fileName << std::endl;
				testFile->Close();
				delete testFile;
				nBadFiles++;
				continue;
			}//<<>>if( tree missing )
			testFile->Close();
			delete testFile;
			if( fInTree->Add(fileName.c_str()) > 0 ) nAdded++;
		}//<<>>for( input files )
		std::cout << "[INFO] filesAdded=" << nAdded << " badFiles=" << nBadFiles << std::endl;
		if( nAdded == 0 ){
			std::cerr << "ERROR: no valid ROOT files were added for " << subdir << " " << matchstr << std::endl;
			delete fInTree;
			return false;
		}//<<>>if( no files added )

		auto treeEntries = fInTree->GetEntries();
		if( treeEntries <= 0 ){
			std::cerr << "ERROR: calibration input chain has no entries." << std::endl;
			delete fInTree;
			return false;
		}//<<>>if( empty chain )

		UInt_t run = 0;
		std::vector<uInt>* rhID = nullptr;
		std::vector<float>* rhTime = nullptr;
		std::vector<float>* rhEnergy = nullptr;
		std::vector<bool>* rhisGS6 = nullptr;
		std::vector<bool>* rhisGS1 = nullptr;
		std::vector<float>* rhTimeError = nullptr;
		TBranch* b_run = nullptr;
		TBranch* b_rhID = nullptr;
		TBranch* b_rhTime = nullptr;
		TBranch* b_rhEnergy = nullptr;
		TBranch* b_rhisGS6 = nullptr;
		TBranch* b_rhisGS1 = nullptr;
		TBranch* b_rhTimeError = nullptr;
		const bool useCCTime = isCC && doUnCC;
		const bool useTimeError = not isCC;
		const char* timeBranchName = useCCTime ? "rhCCTime" : "rhRtTime";

		fInTree->SetBranchStatus("*",0);
		bool branchesValid = true;
		auto requireBranch = [&]( const char* name, auto* address, TBranch** branch ){
			if( not fInTree->GetBranch(name) ){
				std::cerr << "ERROR: required input branch is missing: " << name << std::endl;
				branchesValid = false;
				return;
			}//<<>>if( branch missing )
			fInTree->SetBranchStatus(name,1);
			if( fInTree->SetBranchAddress(name,address,branch) < 0 ){
				std::cerr << "ERROR: failed to bind input branch: " << name << std::endl;
				branchesValid = false;
			}//<<>>if( branch binding failed )
		};
		requireBranch("run",&run,&b_run);
		requireBranch("rhID",&rhID,&b_rhID);
		requireBranch(timeBranchName,&rhTime,&b_rhTime);
		requireBranch("rhEnergy",&rhEnergy,&b_rhEnergy);
		requireBranch("rhisGS6",&rhisGS6,&b_rhisGS6);
		requireBranch("rhisGS1",&rhisGS1,&b_rhisGS1);
		if( useTimeError ) requireBranch("rhTimeError",&rhTimeError,&b_rhTimeError);
		if( not branchesValid ){
			delete fInTree;
			return false;
		}//<<>>if( invalid branch contract )

		fInTree->SetCacheSize(100LL*1024LL*1024LL);
		fInTree->AddBranchToCache("run",true);
		fInTree->AddBranchToCache("rhID",true);
		fInTree->AddBranchToCache(timeBranchName,true);
		fInTree->AddBranchToCache("rhEnergy",true);
		fInTree->AddBranchToCache("rhisGS6",true);
		fInTree->AddBranchToCache("rhisGS1",true);
		if( useTimeError ) fInTree->AddBranchToCache("rhTimeError",true);
		fInTree->StopCacheLearningPhase();

		const auto nEntries = small ? std::min<Long64_t>(treeEntries,10000000) : treeEntries;
		std::cout << "Starting combined entry loop for " << nEntries << " entries" << std::endl;
		int cachedRun = -1;
		kucms_CaliRunClass* ttRange = nullptr;
		kucms_CaliRunClass* crystalRange = nullptr;
		int nextProgressPercent = 0;
		Long64_t nextProgressEntry = 0;
		for( Long64_t centry = 0; centry < nEntries; centry++ ){
			if( centry >= nextProgressEntry ){
				const int completedPercent = static_cast<int>(100.0L*centry/nEntries);
				const int reportPercent = std::min(90,(completedPercent/10)*10);
				reportEventLoopProgress(centry,nEntries,reportPercent);
				nextProgressPercent = reportPercent + 10;
				nextProgressEntry = progressEntryForPercent(nEntries,nextProgressPercent);
			}//<<>>if( progress update )
			auto entry = fInTree->LoadTree(centry);
			if( entry < 0 || b_run->GetEntry(entry) <= 0 ){
				std::cerr << "ERROR: failed to read run at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( run read failed )
			if( run < static_cast<UInt_t>(srun) || run > static_cast<UInt_t>(erun) ) continue;

			if( cachedRun != static_cast<int>(run) ){
				ttRange = findOrCreateRange(TTCaliRunMapSet,tag,iovMaps[curTTIov],run,"TT");
				crystalRange = findOrCreateRange(CaliRunMapSet,tag,iovMaps[curXIov],run,"X");
				if( not ttRange || not crystalRange ){
					delete fInTree;
					return false;
				}//<<>>if( range setup failed )
				cachedRun = run;
			}//<<>>if( run changed )
			maxRunSeenByTag[tag] = std::max(maxRunSeenByTag[tag],static_cast<int>(run));
			const bool fillTT = run > static_cast<UInt_t>(ttRange->lastRun);
			const bool fillCrystal = run > static_cast<UInt_t>(crystalRange->lastRun);
			if( not fillTT && not fillCrystal ) continue;

			bool entryRead = b_rhID->GetEntry(entry) > 0 && b_rhTime->GetEntry(entry) > 0 &&
			                 b_rhEnergy->GetEntry(entry) > 0 && b_rhisGS6->GetEntry(entry) > 0 &&
			                 b_rhisGS1->GetEntry(entry) > 0;
			if( useTimeError ) entryRead = entryRead && b_rhTimeError->GetEntry(entry) > 0;
			if( not entryRead ){
				std::cerr << "ERROR: failed to read rechit branches at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( vector read failed )
			if( not rhID || not rhTime || not rhEnergy || not rhisGS6 || not rhisGS1 ||
			    ( useTimeError && not rhTimeError ) ){
				std::cerr << "ERROR: null rechit vector at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( null vector )
			const auto nRecHits = rhID->size();
			if( rhTime->size() != nRecHits || rhEnergy->size() != nRecHits || rhisGS6->size() != nRecHits ||
			    rhisGS1->size() != nRecHits || ( useTimeError && rhTimeError->size() != nRecHits ) ){
				std::cerr << "ERROR: rechit branch size mismatch at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( vector size mismatch )

			if( fillTT ) maxTTRun[tag][ttRange->startRun] =
				std::max(maxTTRun[tag][ttRange->startRun],static_cast<int>(run));
			RawAccumulator* rawCrystal = nullptr;
			if( fillCrystal ){
				maxCrystalRun[tag][crystalRange->startRun] =
					std::max(maxCrystalRun[tag][crystalRange->startRun],static_cast<int>(run));
				rawCrystal = &rawCrystalByTTIov[tag][crystalRange->startRun][ttRange->startRun];
				if( rawCrystal->empty() ) rawCrystal->reserve(fastDetInfo.size());
			}//<<>>if( crystal fill enabled )
			for( std::size_t idx = 0; idx < nRecHits; idx++ ){
				if( (*rhisGS6)[idx] || (*rhisGS1)[idx] || (*rhEnergy)[idx] < 5.0f ) continue;
				const uInt id = (*rhID)[idx];
				auto det = fastDetInfo.find(id);
				if( det == fastDetInfo.end() ){
					std::cerr << "ERROR: unknown ECAL detID in input: " << id << std::endl;
					delete fInTree;
					return false;
				}//<<>>if( detID missing )
				if( useTimeError ){
					const float errorThreshold = det->second.isEB ? 0.61f : 1.01f;
					const float timeError = (*rhTimeError)[idx];
					if( timeError <= errorThreshold || timeError > 900.0f ) continue;
				}//<<>>if( time error selection enabled )
				const float time = (*rhTime)[idx];
				if( not std::isfinite(time) ){
					rejectedNonFiniteTimes++;
					continue;
				}//<<>>if( non-finite time )
				if( time <= -29.0f ){
					rejectedSentinelTimes++;
					continue;
				}//<<>>if( sentinel time )

				const uInt ttId = det->second.ttId;
				if( fillTT && not ttRange->fillSumCnt(ttId,time) ){
					delete fInTree;
					return false;
				}//<<>>if( TT fill failed )
				if( rawCrystal ){
					auto raw = rawCrystal->try_emplace(id,0.0,0.0,0.0).first;
					raw->second.sum += time;
					raw->second.sumsqr += static_cast<double>(time)*time;
					raw->second.cnt += 1.0;
				}//<<>>if( crystal fill enabled )
			}//<<>>for( rechits )
		}//<<>>for( entries )
		reportEventLoopProgress(nEntries,nEntries,100);

		std::cout << "[INFO] Finished combined batch for " << subdir << " " << matchstr
		          << " entries=" << nEntries << std::endl;
		inputSetsProcessed++;
		delete fInTree;
	}//<<>>while( input rows )

	if( inputSetsProcessed == 0 ){
		std::cerr << "ERROR: input list contains no calibration data sets." << std::endl;
		return false;
	}//<<>>if( no input sets )
	if( doLastRun != 1 ){
		std::cerr << "ERROR: combined gain-1 processing must receive all chunks in one invocation and the final "
		          << "input row must set dolast=1; raw TT-partitioned crystal moments are not persisted." << std::endl;
		return false;
	}//<<>>if( input is incomplete )

	for( const auto& tagEntry : maxRunSeenByTag ){
		auto& runSet = TTCaliRunMapSet[tagEntry.first];
		for( auto& runEntry : runSet ){
			auto& range = runEntry.second;
			bool rangeObserved = false;
			int maximumRangeRun = 0;
			auto observedTag = maxTTRun.find(tagEntry.first);
			if( observedTag != maxTTRun.end() ){
				auto observedRun = observedTag->second.find(runEntry.first);
				if( observedRun != observedTag->second.end() ){
					rangeObserved = true;
					maximumRangeRun = observedRun->second;
				}//<<>>if( range observed )
			}//<<>>if( tag observed )
			if( tagEntry.second >= range.endRun ) range.lastRun = range.endRun;
			else if( rangeObserved ) range.lastRun = std::min(maximumRangeRun,range.endRun);
			if( range.updated && not range.makeMeanMap(false) ) return false;
			range.updated = false;
		}//<<>>for( TT ranges )
	}//<<>>for( TT tags )

	for( auto& tagEntry : rawCrystalByTTIov ){
		for( auto& crystalEntry : tagEntry.second ){
			auto& crystalRange = CaliRunMapSet[tagEntry.first].at(crystalEntry.first);
			for( auto& ttEntry : crystalEntry.second ){
				auto& ttRange = TTCaliRunMapSet[tagEntry.first].at(ttEntry.first);
				for( const auto& rawEntry : ttEntry.second ){
					auto det = fastDetInfo.find(rawEntry.first);
					if( det == fastDetInfo.end() ) return false;
					const uInt ttId = det->second.ttId;
					auto ttMean = ttRange.meanMap.find(ttId);
					if( ttMean == ttRange.meanMap.end() || not std::isfinite(ttMean->second) ){
						std::cerr << "ERROR: no valid TT mean for detID " << rawEntry.first << " in TT IOV "
						          << ttEntry.first << std::endl;
						return false;
					}//<<>>if( TT mean missing )
					const auto& raw = rawEntry.second;
					const double correction = ttMean->second;
					const double correctedSum = raw.sum - raw.cnt*correction;
					const double correctedSumSqr = raw.sumsqr - 2.0*correction*raw.sum +
					                               raw.cnt*correction*correction;
					if( not crystalRange.mergeSumCnt(rawEntry.first,correctedSum,correctedSumSqr,raw.cnt) )
						return false;
				}//<<>>for( crystal accumulators )
			}//<<>>for( TT IOV partitions )
		}//<<>>for( crystal IOV ranges )
	}//<<>>for( tags )

	for( const auto& tagEntry : maxRunSeenByTag ){
		auto& runSet = CaliRunMapSet[tagEntry.first];
		for( auto& runEntry : runSet ){
			auto& range = runEntry.second;
			bool rangeObserved = false;
			int maximumRangeRun = 0;
			auto observedTag = maxCrystalRun.find(tagEntry.first);
			if( observedTag != maxCrystalRun.end() ){
				auto observedRun = observedTag->second.find(runEntry.first);
				if( observedRun != observedTag->second.end() ){
					rangeObserved = true;
					maximumRangeRun = observedRun->second;
				}//<<>>if( range observed )
			}//<<>>if( tag observed )
			if( tagEntry.second >= range.endRun ) range.lastRun = range.endRun;
			else if( rangeObserved ) range.lastRun = std::min(maximumRangeRun,range.endRun);
			if( range.updated && not range.makeMeanMap(false) ) return false;
			range.updated = false;
		}//<<>>for( crystal ranges )
	}//<<>>for( crystal tags )

	std::cout << "Rejected invalid raw times: sentinel=" << rejectedSentinelTimes
	          << " non-finite=" << rejectedNonFiniteTimes << std::endl;
	return true;

}//<<>>bool KUCMS_TimeCalibration::makeGain1CaliMapsEGR( std::string inputFileName, bool small )

bool KUCMS_TimeCalibration::makeCaliMapsEGR( std::string inputFileName, bool doTT, int GID,  bool small, bool doCali ){

	std::string whichstring = ( doTT ) ? "for TT " : "for Xtal ";
	std::cout << "Creating calibration files from EgammaRes Ntuples " << whichstring << std::endl; 
	if( not caliMapsValid ){
		std::cerr << "ERROR: cannot process calibration input because stored calibration maps are invalid." << std::endl;
		return false;
	}//<<>>if( stored maps invalid )

	bool debug = false;
    //bool debug = true;

	//bool doMC = true;
    bool doMC = false;

	// make generic fill with specialized calling functions for input and TT v X
	// assume bool doTT to inicate TT of X Cali set

	const std::string treename("tree/llpgtree");

    // Declaration of leaf types
    UInt_t run;
    std::vector<uInt> *rhID = 0;
    std::vector<float> *rhRtTime = 0;
    std::vector<float> *rhCCTime = 0;
    std::vector<float> *rhEnergy = 0;
    std::vector<bool> *rhisGS6 = 0;
    std::vector<bool> *rhisGS1 = 0;
	std::vector<float> *rhTimeError = 0;

	//std::vector<float> *rhAmp = 0;

    // List of branches
	TBranch *b_run = nullptr;   //!
	TBranch *b_rhID = nullptr;   //!
	TBranch *b_rhRtTime = nullptr;   //!
	TBranch *b_rhCCTime = nullptr;   //!
	TBranch *b_rhEnergy = nullptr;   //!
	TBranch *b_rhisGS6 = nullptr;   //!
	TBranch *b_rhisGS1 = nullptr;   //!
	TBranch *b_rhTimeError = nullptr;
    //TBranch *b_rhAmp; // base for resolutions 

	//Set current IOV tag to corrcet version for this gain
	if( GID == 1 ){ std::cout << " -------- GainID 1 Processing " << std::endl; }
	else if( GID == 2 ){ std::cout << " -------- GainID 2+ Processing " << std::endl; }
	//else if( GID == 2 ){ setTTIov( curTTIov + "_gs2" ); setXIov( curXIov + "_gs2" ); }
    //else if( GID == 3 ){ setTTIov( curTTIov + "_gs3" ); setXIov( curXIov + "_gs3" ); }
    //else if( GID == 4 ){ setTTIov( curTTIov + "_gs4" ); setXIov( curXIov + "_gs4" ); }

	// make sure iov maps exist for cali maps we are making
	if( doTT && iovMaps.find(curTTIov) == iovMaps.end() ){ std::cerr << "ERROR: TT IOV is not valid." << std::endl; return false; }
	if( not doTT && iovMaps.find(curXIov) == iovMaps.end() ){ std::cerr << "ERROR: X IOV is not valid." << std::endl; return false; }

	// basicly a bool (0,1) switch to indicate if we should close processing with this file or more input is coming ( false = not done ) 
	int doLastRun = 0;
	unsigned long long rejectedSentinelTimes = 0;
	unsigned long long rejectedNonFiniteTimes = 0;
	unsigned int inputSetsProcessed = 0;


    std::cout << "open input files list : " << inputFileName << std::endl;
	float fillwgt = 1;
    std::ifstream infilelist(inputFileName);
    if( !infilelist.is_open() ){
        std::cerr << "ERROR: could not open input file list: " << inputFileName << std::endl;
		return false;
    }//<<>>if( !infilelist.is_open() )
    std::string infilestr;
	int finalRun = 0;
    while( std::getline( infilelist, infilestr ) ) {

		std::cout << "Processing : " << infilestr;
		if( infilestr.empty() || infilestr[0] == '#' ){ std::cout << " -> Aborted" << std::endl;  continue; }
		std::cout << std::endl;
        std::stringstream ss(infilestr);
        std::string subdir, matchstr, tag;
        int srun = 0, erun = 0, dolast = 0, wgt = 1;

		bool parsed = doMC ? bool( ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast >> wgt )
		                   : bool( ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast );
		if( not parsed || subdir.empty() || matchstr.empty() || tag.empty() || srun > erun || ( dolast != 0 && dolast != 1 ) ){
			std::cerr << "ERROR: malformed calibration input row: " << infilestr << std::endl;
			return false;
		}//<<>>if( malformed row )
		if( subdir[0] == '#' ){ std::cout << " -- skipping: " << subdir << " " << matchstr << std::endl; continue; }
		std::string wichtype = doTT ? "TT" : "X";

		doLastRun = dolast;
		if( doMC ) fillwgt = wgt;
	
		if( isCC and doUnCC ) tag += "unc";
        std::cout << "open input file : " << subdir << " " << matchstr << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing maps for " << tag << " in " << wichtype << std::endl;

		// insure we have TTCalimaps if we are doing Xtal calibration maps
		if( not doTT && TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() ){ 
			std::cout << " No TT maps for this tag !!" << std::endl; 
			return false;
		}//<<>>if( not doTT && TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

		// insure we have TTCalimaps and XCalimaps if we are doing high gain calibration maps
        if( GID > 1 && doCali ){ 
			if( TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() ){ 
				std::cout << " No TT for this tag in High Gain Calibration!!" << std::endl; 
				return false;
			}//<<>>if( TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )
			if( CaliRunMapSet.find(tag) == CaliRunMapSet.end() ){ 
				std::cout << " No X maps for this tag in High Gain Calibration!!" << std::endl;
					return false;
			}//<<>>if( CaliRunMapSet.find(tag) == CaliRunMapSet.end() )
        }//<<>>if( not doTT && TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

		// ID 10 T error checking
		if( GID > 1 && doTT ){
			std::cerr << "ERROR: high-gain calibration cannot create TT maps." << std::endl;
			return false;
		}//<<>>if( GID > 1 && doTT )

		int nAdded = 0;
        int nLines = 0;
        int nBadFiles = 0;
        auto fInTree = new TChain( treename.c_str() );
        std::cout << "Adding files to TChain." << std::endl;

        std::string eosHost = "cmseos.fnal.gov";
        std::string scanDir = eosDir + inDir + subdir;
        std::cout << "Finding ROOT files under EOS directory: " << scanDir << std::endl;
        std::vector<std::string> files = findEOSRootFiles(eosHost, scanDir, matchstr);
        std::cout << "Found " << files.size() << " candidate ROOT files." << std::endl;

        for (const auto& tfilename : files) {

            nLines++;
            TFile* testFile = TFile::Open(tfilename.c_str(), "READ");

            if( !testFile || testFile->IsZombie() || !testFile->IsOpen() ){
                std::cerr << "\nERROR: bad ROOT file, skipping: " << tfilename << std::endl;
                if( testFile ) testFile->Close();
                nBadFiles++;
                std::cout << "Z";
                continue;
            }//<<>>if( !testFile || testFile->IsZombie() || !testFile->IsOpen() )

            TTree* testTree = nullptr;
            testFile->GetObject(treename.c_str(), testTree);

            if( !testTree ){
                std::cerr << "\nERROR: missing tree " << treename << ", skipping: " << tfilename << std::endl;
                testFile->Close();
                nBadFiles++;
                std::cout << "T";
                continue;
            }//<<>>if( !testTree )

            testFile->Close();
            //std::cout << "opening tfile : " << tfilename << std::endl;
            std::cout << "-";
            int result = fInTree->Add(tfilename.c_str());

            if( result > 0 ){ std::cout << "-"; nAdded++; }
            else std::cout << "0";

        }//<<>>while (std::getline(infile,str))
		std::cout << std::endl;
		std::cout << " -- [INFO] lines=" << nLines << " filesAdded=" << nAdded << " badFiles=" << nBadFiles << std::endl;
		if( nAdded == 0 ){
			std::cerr << "ERROR: no valid ROOT files were added for " << subdir << " " << matchstr << std::endl;
			delete fInTree;
			return false;
		}//<<>>if( nAdded == 0 )

       	run = 0;
        rhID = 0;
        rhRtTime = 0;
        rhCCTime = 0;
        rhEnergy = 0;
		rhisGS6 = 0;
		rhisGS1 = 0;
		rhTimeError = 0;

		// Use one internally aligned rechit collection for IDs, times, energies, flags, and errors.
		const bool useCCTime = isCC && doUnCC;
		const bool useTimeError = not isCC;
		fInTree->SetBranchStatus("*",0);
		bool branchesValid = true;
		auto requireBranch = [&]( const char* name, auto* address, TBranch** branch ){
			if( not fInTree->GetBranch(name) ){
				std::cerr << "ERROR: required input branch is missing: " << name << std::endl;
				branchesValid = false;
				return;
			}//<<>>if( branch missing )
			fInTree->SetBranchStatus(name,1);
			if( fInTree->SetBranchAddress(name, address, branch) < 0 ){
				std::cerr << "ERROR: failed to bind input branch: " << name << std::endl;
				branchesValid = false;
			}//<<>>if( branch binding failed )
		};
		requireBranch("run", &run, &b_run);
		requireBranch("rhID", &rhID, &b_rhID);
		if( useCCTime ) requireBranch("rhCCTime", &rhCCTime, &b_rhCCTime);
		else requireBranch("rhRtTime", &rhRtTime, &b_rhRtTime);
		requireBranch("rhEnergy", &rhEnergy, &b_rhEnergy);
		requireBranch("rhisGS6", &rhisGS6, &b_rhisGS6);
		requireBranch("rhisGS1", &rhisGS1, &b_rhisGS1);
		if( useTimeError ) requireBranch("rhTimeError", &rhTimeError, &b_rhTimeError);
		if( not branchesValid ){
			delete fInTree;
			return false;
		}//<<>>if( not branchesValid )
		fInTree->SetCacheSize(100LL*1024LL*1024LL);
		fInTree->AddBranchToCache("run",true);
		fInTree->AddBranchToCache("rhID",true);
		fInTree->AddBranchToCache(useCCTime ? "rhCCTime" : "rhRtTime",true);
		fInTree->AddBranchToCache("rhEnergy",true);
		fInTree->AddBranchToCache("rhisGS6",true);
		fInTree->AddBranchToCache("rhisGS1",true);
		if( useTimeError ) fInTree->AddBranchToCache("rhTimeError",true);
		fInTree->StopCacheLearningPhase();

		auto treeEntries = fInTree->GetEntries();
		if( treeEntries <= 0 ){
			std::cerr << "ERROR: calibration input chain has no entries." << std::endl;
			delete fInTree;
			return false;
		}//<<>>if( treeEntries <= 0 )
		auto nEntries = small ? ( treeEntries < 10000000 ) ? treeEntries : 10000000  : treeEntries;
        std::cout << "Starting entry loops for " << nEntries << " Entries" <<  std::endl;
        //if( debug ) nEntries = 10000;
		int nextProgressPercent = 0;
		Long64_t nextProgressEntry = 0;
        for (Long64_t centry = 0; centry < nEntries; centry++){

			if( centry >= nextProgressEntry ){
				const int completedPercent = static_cast<int>(100.0L*centry/nEntries);
				const int reportPercent = std::min(90,(completedPercent/10)*10);
				reportEventLoopProgress(centry,nEntries,reportPercent);
				nextProgressPercent = reportPercent + 10;
				nextProgressEntry = progressEntryForPercent(nEntries,nextProgressPercent);
			}//<<>>if( progress update )

			auto entry = fInTree->LoadTree(centry);
			if( entry < 0 ){
				std::cerr << "ERROR: failed to load calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( entry load failed )
			if( b_run->GetEntry(entry) <= 0 ){
				std::cerr << "ERROR: failed to read run at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( run read failed )
			if( run < static_cast<UInt_t>(srun) || run > static_cast<UInt_t>(erun) ) continue;
			TBranch* timeBranch = useCCTime ? b_rhCCTime : b_rhRtTime;
			bool entryRead = b_rhID->GetEntry(entry) > 0 && timeBranch->GetEntry(entry) > 0 &&
			                 b_rhEnergy->GetEntry(entry) > 0 && b_rhisGS6->GetEntry(entry) > 0 &&
			                 b_rhisGS1->GetEntry(entry) > 0;
			if( useTimeError ) entryRead = entryRead && b_rhTimeError->GetEntry(entry) > 0;
			if( not entryRead ){
				std::cerr << "ERROR: failed to read one or more required branches at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( branch read failed )
			std::vector<float>* selectedTime = useCCTime ? rhCCTime : rhRtTime;
			if( not rhID || not selectedTime || not rhEnergy || not rhisGS6 || not rhisGS1 ||
			    ( useTimeError && not rhTimeError ) ){
				std::cerr << "ERROR: null rechit vector at calibration entry " << centry << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( null vector )
			const auto nRecHitsSize = rhID->size();
			bool vectorSizesMatch = selectedTime->size() == nRecHitsSize && rhEnergy->size() == nRecHitsSize &&
			                        rhisGS6->size() == nRecHitsSize && rhisGS1->size() == nRecHitsSize &&
			                        ( not useTimeError || rhTimeError->size() == nRecHitsSize );
			if( not vectorSizesMatch ){
				std::cerr << "ERROR: rechit branch size mismatch at calibration entry " << centry
				          << ": rhID=" << nRecHitsSize << " " << ( useCCTime ? "rhCCTime=" : "rhRtTime=" )
				          << selectedTime->size()
				          << " rhEnergy=" << rhEnergy->size() << " rhisGS6=" << rhisGS6->size()
				          << " rhisGS1=" << rhisGS1->size();
				if( useTimeError ) std::cerr << " rhTimeError=" << rhTimeError->size();
				std::cerr << std::endl;
				delete fInTree;
				return false;
			}//<<>>if( not vectorSizesMatch )
			const int nRecHits = static_cast<int>(nRecHitsSize);

			//if( debug) std::cout << " processing " << run;
            //if( debug) std::cout << " in " << srun << " to " << erun;
			// we are assuming that inputs will run over a continous run range
			// all runs prior ( inclusive ) to the most recent run processed in this input set are !! finished !! and should not be added to.
			// for high gain a one time exception is made to allow adding high gain cali to TT calibration maps
            if( finalRun < run ) finalRun = run; 
			//if( debug) std::cout << " eval w/ final : " << finalRun << std::endl;

			//if( debug) std::cout << " - looping over nRecHits: " << nRecHits << std::endl;
			bool tagNotFound( true );
            bool runNotFound( true );
			auto& calirunset = doTT ? TTCaliRunMapSet : CaliRunMapSet;
			if( calirunset.find(tag) != calirunset.end() ){ // tag exists

				//if( debug) std::cout << " - tag found : " << tag << std::endl;
				tagNotFound = false;
				auto& runset = calirunset[tag];
				//if( debug) std::cout << " - In runset : " << tag << " " << runset[run].startRun << std::endl;
				for( auto& runrange : runset ){
					auto& range = runrange.second;
					//if( debug) std::cout << " - run chk : " << run << " " << range.startRun << " " << range.endRun << std::endl;
					if( run >= range.startRun && run <= range.endRun ){ // run range exists 
						runNotFound = false;
						//if( debug) std::cout << " - check lastRun : " << run << " last " << range.lastRun;
						//if( debug) std::cout  << " nrh " << nRecHits << std::endl; 
						// if end == last : all runs in range completed, if run > last : run in range already filled
						if( GID > 1 ) range.doGID2 = true;
						if( run > range.lastRun || ( not range.hasGID2 && range.doGID2 ) ){

							for( int idx = 0; idx < nRecHits; idx++ ){

								bool gs6 = (*rhisGS6)[idx];
								bool gs1 = (*rhisGS1)[idx];
                        		//if( gs6 || gs1 ) std::cout << " -- : time error " << rhTimeError->at(idx) << std::endl;
                        		//if(debug) std::cout << rhEnergy->at(idx) << " " << rhID->at(idx) << " " << rhRtTime->at(idx) << std::endl;
                        		bool gid1 = not gs6 && not gs1;
                        		bool gid2 = gs6 && not gs1;
                        		bool gid3 = not gs6 && gs1;
                        		bool gid4 = gs6 && gs1;
                        		int whichGS  = gid1 ? 1 : gid2 ? 2 : 3;
                        		if( GID > 1 && gid1 ) continue;
                        		if( GID == 1 && not gid1 ) continue;
                                //if( debug) std::cout << " -- : time error " << rhTimeError->at(idx) << std::endl;
                                //if( debug) std::cout << " -- : GID " << GID << " gs6 " << gs6 << " gs1 " << gs1;
                                //if( debug) std::cout << " e " << rhEnergy->at(idx) << std::endl;
                                
                        		float minRhEnergy = 5.0;
								if( gid2 ){ minRhEnergy = 120.0; }
								else if( gid3 || gid4 ){ minRhEnergy = 240.0; }
								if( (*rhEnergy)[idx] < minRhEnergy ) continue;

	                                bool skiprechit = false;
	                                uInt id = (*rhID)[idx];
								auto det = DetIDMap.find(id);
								if( det == DetIDMap.end() ){
									std::cerr << "ERROR: unknown ECAL detID in input: " << id << std::endl;
									delete fInTree;
									return false;
								}//<<>>if( detID missing )
								if( not isCC ){
									float errorThres = gid1 ? 0.61 : 0.6;
									if( det->second.ecal != ECAL::EB ) errorThres = gid1 ? 1.01 : 1.0;
									if( (*rhTimeError)[idx] <= errorThres ) skiprechit = true;
									if( (*rhTimeError)[idx] > 900 ) skiprechit = true;
								};//<<>>if( isCC ) skiprechit = false; 
								if( skiprechit ) continue;

								// old branch names ( with online prodution post 2024 rhRtTime is cor CC, rhCCTime is uncorr CC )
								float time = (*selectedTime)[idx];
								if( not std::isfinite(time) ){
									rejectedNonFiniteTimes++;
									continue;
								}//<<>>if( not finite time )
								if( time <= -29.0f ){
									rejectedSentinelTimes++;
									continue;
								}//<<>>if( invalid uncorrected CC sentinel )
								uInt crsid = doTT ? getTTId( id, whichGS ) : id;
                        		float bcrstime = ( not doTT && doCali ) ? time - getTTCali( id, run, tag ) : time;

                        		uInt hgcrsid = 0;
                        		float hgcrstime = 0;
                        		bool filter = false;
                        		bool doHgGs = GID > 1;
                        		if( doHgGs ){
                            		hgcrsid = getGSRhId( id, whichGS );
                            		if( hgcrsid < 99 or hgcrsid > 300 ){
                                		std::cout << " doHgGs with " << crsid << " ( " << id << ", " << whichGS << " )" << std::endl;
                            		}//<<>>if( crsid < 0 or )
                            		hgcrstime = time - getCalibration( id, run, tag, 1 );
                            		if( std::abs( hgcrstime ) > 4.f ) filter = true;
                        		}//<<>>if( doHgGs )
                        		if( filter ) continue;

                        		uInt fcrsid = doHgGs ? hgcrsid : crsid;
                        		float fcrstime = doHgGs ? hgcrstime : bcrstime;
								if( not range.fillSumCnt( fcrsid, fcrstime, fillwgt ) ){
									delete fInTree;
									return false;
								}//<<>>if( fill failed )

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
				// ID 10 T check
				if( not doCali && GID > 1 ){
					std::cerr << "ERROR: new IOV range encountered for high-gain calibration with doCali false." << std::endl;
					delete fInTree;
					return false;
				}//<<>>if( not doCali && GID > 1 )

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

	                        bool gs6 = (*rhisGS6)[idx];
	                        bool gs1 = (*rhisGS1)[idx];
                        //if( gs6 || gs1 ) std::cout << " -- : time error " << rhTimeError->at(idx) << std::endl;
                        //if( debug) std::cout << rhEnergy->at(idx) << " " << rhID->at(idx) << " " << rhRtTime->at(idx) << std::endl;
                        bool gid1 = not gs6 && not gs1;
                        bool gid2 = gs6 && not gs1;
                        bool gid3 = not gs6 && gs1;
                        bool gid4 = gs6 && gs1;
                        int whichGS  = gid1 ? 1 : gid2 ? 2 : 3;
						if( GID > 1 && gid1 ) continue;
                        if( GID == 1 && not gid1 ) continue;

                        float minRhEnergy = 5.0;
	                        if( gid2 ){ minRhEnergy = 120.0; }
	                        else if( gid3 || gid4 ){ minRhEnergy = 240.0; }
	                        if( (*rhEnergy)[idx] < minRhEnergy ) continue;

	                        bool skiprechit = false;
						uInt id = (*rhID)[idx];
						auto det = DetIDMap.find(id);
						if( det == DetIDMap.end() ){
							std::cerr << "ERROR: unknown ECAL detID in input: " << id << std::endl;
							delete fInTree;
							return false;
						}//<<>>if( detID missing )
	                        if( not isCC ){
							float errorThres = gid1 ? 0.61 : 0.6;
							if( det->second.ecal != ECAL::EB ) errorThres = gid1 ? 1.01 : 1.0;
							if( (*rhTimeError)[idx] <= errorThres ) skiprechit = true;
							if( (*rhTimeError)[idx] > 900 ) skiprechit = true;
						};//<<>>if( isCC ) skiprechit = false;
                        if( skiprechit ) continue;

						float time = (*selectedTime)[idx];
						if( not std::isfinite(time) ){
							rejectedNonFiniteTimes++;
							continue;
						}//<<>>if( not finite time )
						if( time <= -29.0f ){
							rejectedSentinelTimes++;
							continue;
						}//<<>>if( invalid uncorrected CC sentinel )
						uInt crsid = doTT ? getTTId( id, whichGS ) : id;
                        float bcrstime = ( not doTT && doCali ) ? time - getTTCali( id, run, tag ) : time;

                        uInt hgcrsid = 0;
                        float hgcrstime = 0;
						bool filter = false;
                        bool doHgGs = GID > 1;
                        if( doHgGs ){ 
							hgcrsid = getGSRhId( id, whichGS ); 
                            if( hgcrsid < 99 or hgcrsid > 300 ){
                            	std::cout << " doHgGs with " << crsid << " ( " << id << ", " << whichGS << " )" << std::endl;
                            }//<<>>if( crsid < 0 or )
							hgcrstime = time - getCalibration( id, run, tag, 1 );
							if( std::abs( hgcrstime ) > 4.f ) filter = true;
						}//<<>>if( doHgGs )
						if( filter ) continue;

                        uInt fcrsid = doHgGs ? hgcrsid : crsid;
						float fcrstime = doHgGs ? hgcrstime : bcrstime;
						if( not calirunset[tag][tstart].fillSumCnt( fcrsid, fcrstime, fillwgt ) ){
							delete fInTree;
							return false;
						}//<<>>if( fill failed )

                    }//<<>>for( int idx = 0; idx < nRecHits; idx++ )
					} else {
						std::cerr << "ERROR: no configured IOV contains run " << run << " for tag " << tag << std::endl;
						delete fInTree;
						return false;
					}//<<>>if( configured IOV found )

			}//<<>>if( not tagfound )

		}//<<>>for (Long64_t centry = 0; centry < nEntries; cnentry++)
		reportEventLoopProgress(nEntries,nEntries,100);

        // Finished processing this batch
		std::cout << "[INFO] Finished batch for " << subdir << " " << matchstr << " entries=" << nEntries << std::endl;
		inputSetsProcessed++;
		delete fInTree;
        fInTree = nullptr;

	}//<<>>while (std::getline(infilelist,infilestr))

	if( inputSetsProcessed == 0 ){
		std::cerr << "ERROR: input list contains no calibration data sets." << std::endl;
		return false;
	}//<<>>if( inputSetsProcessed == 0 )
	std::cout << "Rejected invalid raw times: sentinel=" << rejectedSentinelTimes
	          << " non-finite=" << rejectedNonFiniteTimes << std::endl;
	std::cout << "Making Mean Maps from input data : doTT? " << doTT << std::endl;
	// fill map hist -----------------------------------------------
	//int count = 0;
	if( doLastRun == 1 ){
	auto& calirunset = doTT ? TTCaliRunMapSet : CaliRunMapSet;
	for( auto& calirun : calirunset ){ 
		for( auto& runrange : calirun.second ){
			auto& range = runrange.second;
			if( range.updated ){
				std::cout << "For: " << std::to_string( runrange.first ) << std::endl; 
				if( range.lastRun != range.endRun ){
					if( finalRun > 100000 ) range.lastRun = ( finalRun < range.endRun ) ? finalRun : range.endRun;
					else range.lastRun = 999999;
				}//<<>>if( range.lastRun != range.endRun )
				if( not range.makeMeanMap( false ) ) return false;
				range.updated = false;
				if( GID > 1 && ( range.lastRun == range.endRun ) ){ range.hasGID2 = true; range.doGID2 = false; }
			}//<<>>if( range.updated )
		}//<<>>for( auto& runrange : calirun.second )
	}//<<>>for( auto& calirun : calirunset )
	}//<<>>if( dolast == 1 )

	return true;

}//<<>>bool KUCMS_TimeCalibration::makeCaliMapsEGR( std::string inputFileName )

void KUCMS_TimeCalibration::plot2dResolutionEGR( std::string inputFileName, bool scale, bool usecali, bool smear, std::string ext ){

    std::cout << "Creating 2D Resolution Hist from EgammaRes Ntuples " << std::endl;
    if( not usecali ) std::cout << " -- NoCali " << std::endl;
    if( smear ) std::cout << " -- Smeared " << std::endl;

    bool debug = false;
    //bool debug = true;
    bool small = false;
    //bool small = true;

    //bool doMC = false;
    bool doMC = true;
	float fillwgt = 1;

    const std::string treename("tree/llpgtree");

    bool useAmp = not useEffEnergy;

    float lB = 10; // lower and upper limits of energies for rechits used
    float uB = 240; // lower and upper limits of energies for rechits used
    if( lowEnergy ){ lB = 2; uB = 240; }

    std::cout << " -- use low energy : " << lowEnergy << std::endl;
    std::cout << " -- xbins : " << xBinStr << std::endl;
    std::cout << " -- ybins : " << yBinStr << std::endl;
	if( useGSwitch ) std::cout << " -- using GS with high gian included : " << doAllGain << " lower energy bound : " << lB << std::endl;
    else std::cout << " -- upper energy bound : " << uB << " lower energy bound : " << lB << std::endl;

    // Declaration of leaf types
    UInt_t          run;
    std::vector<unsigned int> *resRhID;
    std::vector<float>   *resAmp;
    std::vector<float>   *resE;
    std::vector<float>   *resRtTime;
    std::vector<float>   *resCCTime;
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
    TBranch        *b_resCCTime;   //!
    TBranch        *b_resTOF;   //!

    TBranch        *b_rhID;   //!
    TBranch        *b_rhisGS6;   //!
    TBranch        *b_rhisGS1;   //!

    std::string nocalistr("_NoCali");
    std::string uncorstr("_UnCorr");
    std::string smearstr("_Smeared");
    std::string lochist("_SRO_Data_Hist");
    std::string druhist("_DRO_Data_Hist");
    std::string globhist("_ZEE_Data_Hist");
    std::string ehist("_Amp_v_E_Hist");
    std::string eb_hist("_EB");
    std::string ee_hist("_EE");
    std::string fTitle("#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n}");
    std::string fXTitle("A_{eff}/#sigma_{n}");
    std::string fYTitle("#Delta(Photon Seed Time) [ns]");
    std::string fZTitle("");
    std::string fEYTitle("E_{eff} [GeV]");
    std::string fEBEBTitle(" (EBEB) ");
    std::string fEEEETitle(" (EEEE) ");

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
	if( !infilelist.is_open() ){
    	std::cerr << "ERROR: could not open input file list: " << inputFileName << std::endl;
    	return;
	}//<<>>if( !infilelist.is_open() )
    std::string infilestr;
    while( std::getline( infilelist, infilestr ) ) {

    	// remove leading/trailing whitespace
    	auto first = infilestr.find_first_not_of(" \t\r\n");
    	if( first == std::string::npos ) continue;  // blank or whitespace-only line

    	auto last = infilestr.find_last_not_of(" \t\r\n");
    	infilestr = infilestr.substr(first, last - first + 1);

    	// now comments with leading spaces are also skipped
		if( infilestr.empty() || infilestr[0] == '#' ) continue;

        std::stringstream ss(infilestr);
        std::string subdir, matchstr, tag;
        int srun(0), erun(0), wgt(1), dolast(0);


    	//if( !(ss >> subdir >> matchstr >> srun >> erun >> tag) ){
        //	std::cout << "Bad input-list entry, skipping: [" << infilestr << "]" << std::endl;
        //	continue;
    	//}//<<>>if( !(ss >> subdir >> matchstr >> srun >> erun >> tag) )

        //if( !( ss >> infilename >> srun >> erun >> tag )) std::cout << "[WARN] Parse error in : " << infilestr << std::endl; continue;
        if( doMC ) ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast >> wgt;
        else ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast;
        if( subdir[0] == '#' ){ std::cout << " -- skipping: " << subdir << " " << matchstr << std::endl; continue; }
        if( doMC ) fillwgt = wgt;

        if( subdir == "None" ) subdir = "";
        if( matchstr == "None" ) matchstr = "";

        if( !subdir.empty() && subdir.back() != '/' ){ subdir += "/"; }

        //std::cout << "open input file : " << infilename << std::endl;
        std::cout << "open input entry : "
                  << " subdir=" << subdir
                  << " matchstr=" << matchstr
                  << " srun=" << srun
                  << " erun=" << erun
                  << " tag=" << tag
                  << std::endl;

		if( isCC && doUnCC ){ 
			if( tag == "r3_p25" ) tag = "r3_p25unc";
            else if( tag == "r3_p26" ) tag = "r3_p26unc";
			else { std::cout << "Bad tag found for doCC and use UnCC times " << std::endl; return; }
		}//<<>>if( isCC && doUnCC )

    	if( smear ) std::cout << "with tags : " << tag << " and smearing " << resTag << std::endl;
        else std::cout << "with tags : " << tag << std::endl;
        std::cout << "For Run " << srun << " to Run " << erun << std::endl;
        std::cout << "Producing 2D Resolution Map " << std::endl;

        // insure we have Calimaps if we are doing Xtal calibration maps
        if( usecali && CaliRunMapSet.find(tag) == CaliRunMapSet.end() ){
            std::cout << " No Cali maps for this tag !!" << std::endl;
            return;
        }//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        //if( smear && ( SmearTagSet.find(smearTag) == SmearTagSet.end() ) ){
        //    std::cout << " No Smear tag found !!" << std::endl;
        //    return;
        //}//<<>>if( CaliRunMapSet.find(tag) == TTCaliRunMapSet.end() )

        //int fileskipcnt = 0;
        int nAdded = 0;
        int nLines = 0;
		int nBadFiles = 0;

        //std::ifstream infile(infilename);
		//if( !infile.is_open() ){
    	//	std::cerr << "ERROR: could not open input file list: " << infilename << std::endl;
    	//	continue;
		//}//<<>>if( !f || f->IsZombie() )
        //std::string instr;

        auto fInTree = new TChain( treename.c_str() );

		//std::string indir("KUCMSNtuple/gammares_prmt25/");
		//std::string eosdir("root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/");
		//  = '/store/user/lpcsusylep/jaking/'
    	// eosh example: "cmseos.fnal.gov"
    	// scanDir example: "/store/user/lpcsusylep/jaking/KUCMSNtuple/someDir"
		// std::string scanDir = "/store/user/lpcsusylep/jaking/KUCMSNtuple/some/final/directory";
    	// matchStr example: "somename.root"
    	// std::vector<std::string> files = findEOSRootFiles(eosHost, scanDir, matchString);
		
		std::string eosdir = "/store/user/lpcsusylep/jaking/";
		std::string eosHost = "cmseos.fnal.gov";
        //auto tfilename = eosDir + inDir + instr;
        std::string scanDir = eosdir + inDir + subdir;
		std::cout << "Finding ROOT files under EOS directory: " << scanDir << std::endl;
		std::vector<std::string> files = findEOSRootFiles(eosHost, scanDir, matchstr);
		std::cout << "Found " << files.size() << " candidate ROOT files." << std::endl;

        std::cout << "Adding files to TChain." << std::endl;
        //while (std::getline(infile,instr)){
		for (const auto& tfilename : files) {

			//if( instr.empty() || instr[0] == '#' ) continue;
            nLines++;
            //fileskipcnt++;
            //if( fileskipcnt%10 != 0 ) continue;
            //auto tfilename = eosDir + inDir + instr;

    		TFile* testFile = TFile::Open(tfilename.c_str(), "READ");

    		if( !testFile || testFile->IsZombie() || !testFile->IsOpen() ){
        		std::cerr << "\nERROR: bad ROOT file, skipping: " << tfilename << std::endl;
        		if( testFile ) testFile->Close();
        		nBadFiles++;
        		std::cout << "Z";
        		continue;
    		}//<<>>if( !testFile || testFile->IsZombie() || !testFile->IsOpen() )

    		TTree* testTree = nullptr;
    		testFile->GetObject(treename.c_str(), testTree);

    		if( !testTree ){
        		std::cerr << "\nERROR: missing tree " << treename << ", skipping: " << tfilename << std::endl;
        		testFile->Close();
        		nBadFiles++;
        		std::cout << "T";
        		continue;
    		}//<<>>if( !testTree )

		    testFile->Close();
            //std::cout << "opening tfile : " << tfilename << std::endl;
            std::cout << "-";
            int result = fInTree->Add(tfilename.c_str());


            if( result > 0 ){ std::cout << "-"; nAdded++; }
            else std::cout << "0";


        }//<<>>while (std::getline(infile,str))
        std::cout << std::endl;
		std::cout << " -- [INFO] lines=" << nLines << " filesAdded=" << nAdded << " badFiles=" << nBadFiles << std::endl;

        run = 0;
        resRhID = 0;
        resAmp = 0;
        resE = 0;
        resRtTime = 0;
		resCCTime = 0;
        resTOF = 0;

        rhID = 0;
        rhisGS6 = 0;
        rhisGS1 = 0;

        fInTree->SetBranchAddress( "run", &run, &b_run );   //!
        fInTree->SetBranchAddress( "resRhID", &resRhID, &b_resRhID );   //!
        fInTree->SetBranchAddress( "resAmp", &resAmp, &b_resAmp);   //!
        fInTree->SetBranchAddress( "resE", &resE, &b_resE);   //!
        fInTree->SetBranchAddress( "resRtTime", &resRtTime, &b_resRtTime);   //!
		if( isCC ) fInTree->SetBranchAddress("resCCTime", &resCCTime, &b_resCCTime);
        fInTree->SetBranchAddress( "resTOF" , &resTOF, &b_resTOF);   //!

        fInTree->SetBranchAddress("rhID", &rhID, &b_rhID);
        fInTree->SetBranchAddress("rhisGS6", &rhisGS6, &b_rhisGS6);
        fInTree->SetBranchAddress("rhisGS1", &rhisGS1, &b_rhisGS1);

        std::cout << " Getting calibration values and plotting" << std::endl;

        auto nEntries = fInTree->GetEntries();
        if( small && nEntries > 1000000 ) nEntries = 1000000;
        if( debug ) nEntries = ( nEntries < 10000 ) ? nEntries : 10000;
        //if( debug )
		std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
        int report = nEntries/20;
        for (auto centry = 0U; centry < nEntries; centry++){

            if( centry%report == 0 or centry == 0){
                std::time_t now = std::time(nullptr);
                std::string curtime = std::ctime(&now);
                curtime.pop_back();
                std::cout << "Proccessed " << centry << " of " << nEntries;
                std::cout << " " << curtime;
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
			if( isCC ) b_resCCTime->GetEntry(entry);   //!
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
							if( doUnCC ) caliend = caliend + uncorstr;
                            if( smear ) caliend = caliend + smearstr;
                            caliend = caliend + ext;
							std::string histMapName = "ResMap_" + std::to_string(srun) + "_" + std::to_string(erun);
                            std::string lsfname = histMapName+lochist+caliend;
                            std::string ldfname = histMapName+druhist+caliend;
                            std::string gbfname = histMapName+globhist+caliend;
                            std::string efname = histMapName+ehist+caliend;
                            std::string lsfnameeb = lsfname + eb_hist;
                            std::string ldfnameeb = ldfname + eb_hist;
                            std::string gbfnameeb = gbfname + eb_hist;
                            std::string efnameeb = efname + eb_hist;
                            std::string lsfnameee = lsfname + ee_hist;
                            std::string ldfnameee = ldfname + ee_hist;
                            std::string gbfnameee = gbfname + ee_hist;
                            std::string efnameee = efname + ee_hist;
                            if( CaliHists.find(lsfnameeb) == CaliHists.end() ){ // 2D hist already made? 
                                TH2F* theHistLS = new TH2F(lsfnameeb.c_str(),lsfnameeb.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[lsfnameeb] = { theHistLS, lsfnameeb, true, true }; // histfile histname isnew isres
                                TH2F* theHistGB = new TH2F(gbfnameeb.c_str(),gbfnameeb.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[gbfnameeb] = { theHistGB, gbfnameeb, true, true }; // histfile histname isnew isres
                                TH2F* theHistLD = new TH2F(ldfnameeb.c_str(),ldfnameeb.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[ldfnameeb] = { theHistLD, ldfnameeb, true, true }; // histfile histname isnew isres
                                TH2F* theHistAE = new TH2F(efnameeb.c_str(),efnameeb.c_str(),500,0,1000,100,0,200);
                                theHistAE->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistAE->GetYaxis()->SetTitle(fEYTitle.c_str());
                                theHistAE->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[efnameeb] = { theHistAE, efnameeb, true, true }; // histfile histname isnew isres
                                //range.has2DResMap = true;
                                TH2F* theHistLSEE = new TH2F(lsfnameee.c_str(),lsfnameee.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLS->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLS->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLS->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[lsfnameee] = { theHistLSEE, lsfnameee, true, true }; // histfile histname isnew isres
                                TH2F* theHistGBEE = new TH2F(gbfnameee.c_str(),gbfnameee.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistGB->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistGB->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistGB->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[gbfnameee] = { theHistGBEE, gbfnameee, true, true }; // histfile histname isnew isres
                                TH2F* theHistLDEE = new TH2F(ldfnameee.c_str(),ldfnameee.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
                                theHistLD->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistLD->GetYaxis()->SetTitle(fYTitle.c_str());
                                theHistLD->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[ldfnameee] = { theHistLDEE, ldfnameee, true, true }; // histfile histname isnew isres
                                TH2F* theHistAEEE = new TH2F(efnameee.c_str(),efnameee.c_str(),500,0,1000,100,0,200);
                                theHistAEEE->GetXaxis()->SetTitle(fXTitle.c_str());
                                theHistAEEE->GetYaxis()->SetTitle(fEYTitle.c_str());
                                theHistAEEE->GetZaxis()->SetTitle(fZTitle.c_str());
                                CaliHists[efnameee] = { theHistAEEE, efnameee, true, true }; // histfile histname isnew isres
                            }//if( CaliHists.find(lsfname) == CaliHists.end() ){ 
                            if( CaliHists[lsfnameeb].isNew ){

                                // ----  fill res hist

                                auto idinfoL0 = DetIDMap[(*resRhID)[0]];
                                auto idinfoL1 = DetIDMap[(*resRhID)[1]];
                                auto idinfoG0 = DetIDMap[(*resRhID)[2]];
                                auto idinfoG1 = DetIDMap[(*resRhID)[3]];

								// set gain level in this block ---------------------------------------------------------------
                                std::vector<bool> isGainId1 = {true,true,true,true};
                                std::vector<int> rhGainId = {0,0,0,0};
                                if( useGSwitch ){
                                    for( int resrhit = 0; resrhit < 4; resrhit++ ){
                                        for( int rhit = 0; rhit < nRecHits; rhit++ ){
                                            if( (*resRhID)[resrhit] == (*rhID)[rhit] ){
                                                if( useGain == 1 ) isGainId1[resrhit] = not ( (*rhisGS1)[rhit] || (*rhisGS6)[rhit] );
												else if ( useGain == 2 ) isGainId1[resrhit] = not (*rhisGS1)[rhit] && (*rhisGS6)[rhit];
                                                else if ( useGain == 3 ) isGainId1[resrhit] = (*rhisGS1)[rhit] && not (*rhisGS6)[rhit];
                                                else if ( useGain == 4 ) isGainId1[resrhit] = (*rhisGS1)[rhit] && (*rhisGS6)[rhit];
												else isGainId1[resrhit] = false;
												rhGainId[resrhit] = useGain;
                                                break;
                                            }//<>if( (*resRhID)[rhit] == (*resRhID)[resrhit] )
                                        }//<<>>for( int resrhit = 0; resrhit < nResRecHits; resrhi++ )
                                    }//<<>>for( int rhit = 0; rhit < 4; rhit++ )
                                }//<<>>if( useGSwitch )

								//if( rhGainId[0] == 0 or rhGainId[1] == 0 or rhGainId[2] == 0 or rhGainId[3] == 0 ) continue;
                                float seedTimeIC00 = ( usecali ) ? getCalibration( (*resRhID)[0], run, tag, rhGainId[0] ) : 0;
                                float seedTimeIC10 = ( usecali ) ? getCalibration( (*resRhID)[1], run, tag, rhGainId[1] ) : 0;
                                float seedTimeIC01 = ( usecali ) ? getCalibration( (*resRhID)[2], run, tag, rhGainId[2] ) : 0;
                                float seedTimeIC11 = ( usecali ) ? getCalibration( (*resRhID)[3], run, tag, rhGainId[3] ) : 0;
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

								bool useUnCCTime = isCC && doUnCC;
								double Time0 = ( useUnCCTime ) ? (*resCCTime)[0] : (*resRtTime)[0];
                                double Time1 = ( useUnCCTime ) ? (*resCCTime)[1] : (*resRtTime)[1];
                                double Time2 = ( useUnCCTime ) ? (*resCCTime)[2] : (*resRtTime)[2];
                                double Time3 = ( useUnCCTime ) ? (*resCCTime)[3] : (*resRtTime)[3];

                                double ldTOF = (*resTOF)[0]-(*resTOF)[1]; //phoseedTOF_0-phoseedTOF_1;
                                double gdTOF = (*resTOF)[2]-(*resTOF)[3]; //phoseedTOF_0-phoseedTOF_1;
                                double lyf0 = Time0-seedTimeIC00;
                                double lyf1 = Time1-seedTimeIC10;
                                double gyf0 = Time2-seedTimeIC01;
                                double gyf1 = Time3-seedTimeIC11;
								//std::cout << " -- : " << Time0 << " - " << seedTimeIC00 << " " << usecali << " " << doEE;
								//std::cout << " " << (*resRhID)[0] << std::endl;

                                if( smear ){
                                    //std::cout << "Times are smeared !!!!!" << std::endl;

									float amp0 = ( (*resRhID)[0] > 0 ) ? (*resAmp)[0] : -999;
                                    float amp1 = ( (*resRhID)[1] > 0 ) ? (*resAmp)[1] : -999;
                                    float amp2 = ( (*resRhID)[2] > 0 ) ? (*resAmp)[2] : -999;
                                    float amp3 = ( (*resRhID)[3] > 0 ) ? (*resAmp)[3] : -999;

                                    //getCorrectedTime( time, amplitude, rechitID, Evt_run, dataSetKey, mctype, gainID )
									lyf0 = getCorrectedTime( lyf0, amp0, (*resRhID)[0], run, tag, 0, rhGainId[0] );
                                    lyf1 = getCorrectedTime( lyf1, amp1, (*resRhID)[1], run, tag, 0, rhGainId[1] );
                                    gyf0 = getCorrectedTime( gyf0, amp2, (*resRhID)[2], run, tag, 0, rhGainId[2] );
                                    gyf1 = getCorrectedTime( gyf1, amp3, (*resRhID)[3], run, tag, 0, rhGainId[3] );


                                }//<<>>if( smear )
                                double lyfill = lyf0-lyf1+ldTOF;
                                double gyfill = gyf0-gyf1+gdTOF;

                                bool le_cut = ((*resE)[0]>=lB)&&((*resE)[0]<=uB)&&((*resE)[1]>=lB)&&((*resE)[1]<=uB);
                                bool ge_cut = ((*resE)[2]>=lB)&&((*resE)[2]<=uB)&&((*resE)[3]>=lB)&&((*resE)[3]<=uB);
                                if( useGSwitch ){
									if( doAllGain ){
                                        le_cut = (*resE)[0] >= lB && (*resE)[1] >= lB;
                                        ge_cut = (*resE)[2] >= lB && (*resE)[3] >= lB;
									} else {
                                    	le_cut = isGainId1[0] && isGainId1[1] && (*resE)[0] >= lB && (*resE)[1] >= lB;
                                    	ge_cut = isGainId1[2] && isGainId1[3] && (*resE)[2] >= lB && (*resE)[3] >= lB;
									}//<<>>if( doAllGain )
                                }//<<>>if( useGSwitch )
                                bool leta_cut_eb = (idinfoL0.ecal == ECAL::EB)&&(idinfoL1.ecal == ECAL::EB);
                                bool geta_cut_eb = (idinfoG0.ecal == ECAL::EB)&&(idinfoG1.ecal == ECAL::EB);
                                bool leta_cut_ee = (idinfoL0.ecal != ECAL::EB)&&(idinfoL1.ecal != ECAL::EB);
                                bool geta_cut_ee = (idinfoG0.ecal != ECAL::EB)&&(idinfoG1.ecal != ECAL::EB);
                                bool goodLocTime = Time0 != 0 && Time1 != 0 && Time0 != Time1;
                                bool goodGloTime = Time2 != 0 && Time3 != 0 && Time2 != Time3;
                                bool goodLocRHs = (*resRhID)[0] != 0 && (*resRhID)[1] != 0;
                                bool goodGloRHs = (*resRhID)[2] != 0 && (*resRhID)[3] != 0;

                                bool isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, fasle = different
                                bool levent_good = le_cut && goodLocRHs && goodLocTime;
                                bool gevent_good = ge_cut && goodGloRHs && goodGloTime;

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
                                    if( leta_cut_eb ){ 
										CaliHists[lsfnameeb].h2f->Fill(lxfill,lyfill,fillwgt); 
										CaliHists[efnameeb].h2f->Fill(lxfill,lxfille,fillwgt); }
                                    if( leta_cut_ee ){ 
										CaliHists[lsfnameee].h2f->Fill(lxfill,lyfill,fillwgt); 
										CaliHists[efnameee].h2f->Fill(lxfill,lxfille,fillwgt); }
                                }//<<>>if( levent_good && isd_cut )
                                if( gevent_good ){ 
									if( leta_cut_eb ) CaliHists[gbfnameeb].h2f->Fill(gxfill,gyfill,fillwgt); 
                                    if( leta_cut_ee ) CaliHists[gbfnameee].h2f->Fill(gxfill,gyfill,fillwgt); 
								}//<<>>if( gevent_good )
                                if( levent_good && not isd_cut ){ 
									if( leta_cut_eb ) CaliHists[ldfnameeb].h2f->Fill(lxfill,lyfill,fillwgt); 
                                    if( leta_cut_ee ) CaliHists[ldfnameee].h2f->Fill(lxfill,lyfill,fillwgt); 
								}//<<>>if( levent_good && not isd_cut )
                                if(debug) std::cout << " - Fill hists done" << std::endl;

                            }//if( CaliHists[lsfname].isNew )
                        //}//<<>>if( run > range.lastRun )
                        //continue;//found the correct run set ( iov range ) no need to look further
                    }//<<>>if( run > srun && run <= erun )
                //}//<<>>for( auto& range : runset )
            //}//<<>>if( calirunset.find(tag) != calirunset.end() )
        } // for (auto entry = 0U; entry < nEntries; entry++)

        if(debug) std::cout << " -------- Next Input file " << std::endl;
        std::cout << "[INFO] Finished batch for " << scanDir << " entries=" << nEntries << std::endl;
        delete fInTree;
        fInTree = nullptr;

    } // while (std::getline(infilelist,infiles))

    // --------------  process new histos ----------------------------------------------

    // ----   scale and set isNew to false

    std::cout << " - Scale & Save 2D Hist" << std::endl;
	std::string cali2dfile_condor = "res2dPlotsTFileCondor.root";
	if( isCondorJob ) cali2DResTFile = TFile::Open( cali2dfile_condor.c_str(), "RECREATE" );
    else cali2DResTFile = TFile::Open( cali2DResPlotsTFileName.c_str(), "UPDATE" );
	if( !cali2DResTFile || cali2DResTFile->IsZombie() ){
    	std::cerr << "ERROR: could not open output ROOT file" << std::endl;
    	return;
	}//<<>>if( !cali2DResTFile || cali2DResTFile->IsZombie() )
    cali2DResTFile->cd();
    for( auto& hists : CaliHists ){

        if( hists.second.isResHist && hists.second.isNew ){

            if( scale ) scaleHist(hists.second.h2f,false,true,false);// hist, scale up?, varible x bins?, varible y bin?
            hists.second.isNew = false;
            if( not small && hists.second.h2f ) hists.second.h2f->Write( hists.second.h2f->GetName(), TObject::kOverwrite );

        }//<<>>if( hists.second.isResHist && hists.second.isNew )

    }//<<>>for( auto& hists : CaliHists )
    cali2DResTFile->Close();

    std::cout << "Finished making 2D delta t v eff amp plots" << std::endl;

}//<<>> void plot2dResolution( std::string indir, std::string infilelistname, 

kucms_SigmaFitResult KUCMS_TimeCalibration::runTimeFitter( TH2F* hist2D ){

  	std::string f2DHistName = hist2D->GetName();
	std::cout << "Running time fitter : " << f2DHistName << std::endl;
    bool doSterm = true;
	//bool doSterm = false;
	
  	std::vector<float> fXBins;
	std::cout << " -- " << xBinStr << std::endl;
 	setBins( xBinStr, fXBins );
  	int fNBinsX = fXBins.size();
	const auto xbins = &fXBins[0];

	std::string fittypename = ( doSterm ) ? "_NSC" : "_NC";
	//std::string outfilename = caliFileDir + f2DHistName + fittypename  + "_resfit.root";
    std::string outfilename = f2DHistName + fittypename  + "_resfit.root";
    TFile* resTFile = TFile::Open( outfilename.c_str(), "RECREATE" );
	if( resTFile == NULL ) std::cout << " -- Failed to open : " << f2DHistName << std::endl;
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

  	std::map<int,kucms_ProfileTimeFit> profileHists;
	//std::vector<kucms_ProfileTimeFit> profileHists;
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
		kucms_ProfileTimeFit thisprofile( profilehist );
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
    	kucms_TimeFitResult results = profile.GetFitResult();

		//float herr = std::abs( results.sigmaHigh - results.sigma );
        //float lerr = std::abs( results.sigmaLow - results.sigma ); 
		//float serr = ( herr > lerr ) ? herr : lerr;
		float serr = results.esigma;		

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
		TH1F* gprHist = gprProfileHist( profile.second.profileHist );
		gprHist->Write( gprHist->GetName(), TObject::kOverwrite );
		delete gprHist;
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
	float x_up = 1600;
	if( lowEnergy ){ x_low = 5; x_up = 1600; }
    if( doAllGain ){ x_up = 3200; }
    if( useEffEnergy ){ x_low = 5; x_up = 240; }
    if( lowEnergy && useEffEnergy ){ x_low = 2; }
    if( doAllGain && useEffEnergy ){ x_up = 2000; }

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

	kucms_SigmaFitResult sigmafit;
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

}//<<>>SigmaFitResult KUCMS_TimeCalibration::runTimeFitter( TH2F* hist2D )

TH1F* KUCMS_TimeCalibration::gprProfileHist( TH1F* hist ){


    std::cout << "Make GPR ProfileHist: " << hist->GetName() << std::endl;

    int nBins = hist->GetNbinsX();
	float sigma = setdec( hist->GetStdDev(), 2 );	

	float minbin = 0.04;
	float thres = 100;
	float hill = 10000;
	int firstBin = -1;
	int lastBin = -1;
	bool overHill = false;
    for( int ibinX = 1; ibinX < nBins; ibinX++ ){

		float content = hist->GetBinContent(ibinX);
		if( content > thres && firstBin < 0 ) firstBin = ibinX;
		if( content > hill ) overHill = true; 
		if( overHill && content < thres && lastBin < 0 ) lastBin = ibinX;

    }//<<>>for (auto ibinX = 1; ibinX <= nBins; ibinX++)

	int startbin = ( firstBin > 2 ) ? firstBin-1 : 1;
	float start = setdec( hist->GetBinCenter(startbin), 2 );
    int endbin = ( lastBin < (nBins-2) ) ? lastBin+1 : nBins-1;
	float end = setdec( hist->GetBinCenter(endbin), 2);
	if( std::abs(start) > std::abs(end) ) end = -1*start;
	else start = -1*end;

    int range = int((end-start)/minbin);
	float secorr =  ((float(range+1)*minbin) - (end-start))/2;
	end += secorr;
	start -= secorr;
	float div = sigma/2;
	int divrange = int(div/minbin);
	div = float(divrange+1)*minbin;
	int binrange = int((end-start)/(div));
	float bcorr = ((float(binrange+1)*div) - (end-start))/2;
	end += bcorr;
    start -= bcorr; 
	int nDiv = int((end-start)/(div));
	if( nDiv%2 == 0 ){
		nDiv++;
		float fcorr = div/2;
		end += fcorr;
		start -= fcorr;
	}//<<>>if( nDiv%2 == 0 )
	std::string title = hist->GetName();
	std::string newtitle = title + "_gpr";
	TH1F* newhist = new TH1F( newtitle.c_str(), newtitle.c_str(), nDiv, start, end );
	//newhist->Sumw2();
	std::cout << " -- start : " << start << " end : " << end << " nDiv : " << nDiv << " bin size : " << (end-start)/nDiv << std::endl;
	for( int ibinX = 1; ibinX <= nBins; ibinX++ ){
	
		float content = hist->GetBinContent(ibinX);
		float bin = hist->GetBinCenter(ibinX);
		newhist->Fill(bin,content);

	}//<<>>for( int ibinX = 1; ibinX <= nBins; ibinX++ )

	int nGPRBins = hist->GetNbinsX();
    for( int ibinX = 1; ibinX < nGPRBins; ibinX++ ){

        float content = hist->GetBinContent(ibinX);
		float error = std::sqrt(content);
		if( error < 1 ) error = 1;
        newhist->SetBinError( ibinX, error );

    }//<<>>for( int ibinX = 1; ibinX <= nBins; ibinX++ )

	return newhist;

}//<<>>TH1D* KUCMS_TimeCalibration::gprProfileHist( TH1D* hist )

void KUCMS_TimeCalibration::doResTimeFits( bool doLocal ){

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
				kucms_SigmaFitResult results = runTimeFitter( CaliHists[fname].h2f );
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

}//<<>>void KUCMS_TimeCalibration::doResTimeFits()

void KUCMS_TimeCalibration::doResTimeFitExt( std::string histName, std::string fileName ){

    std::cout << "Making Resolution Time Fits" << std::endl;
	std::cout << " -- With : " << fileName << std::endl;
	auto tFile = TFile::Open( fileName.c_str(), "READ" );
	if (!tFile || tFile->IsZombie()){ std::cout << " -- Failed to Open : " << fileName << std::endl; return; }
	std::cout << " -- For : " << histName << std::endl;
	TH2F* lhist = (TH2F*)tFile->Get(histName.c_str());
    if( lhist == nullptr ){ std::cout << " -- Failed to Open : " << histName << std::endl; return; }
	kucms_SigmaFitResult results = runTimeFitter( lhist );
	tFile->Close();
	delete tFile;

}//<<>>void KUCMS_TimeCalibration::doResTimeFit( std::string histName )

void KUCMS_TimeCalibration::load2DResHist( std::string histName ){

    std::cout << " --- Load hist from caliTFile: " << histName << std::endl;
	cali2DResTFile = TFile::Open( cali2DResPlotsTFileName.c_str(), "READ" );
	TH2F* lhist = (TH2F*)cali2DResTFile->Get(histName.c_str());
    //std::cout << " --- Found local " << lhist << std::endl;
    if( lhist ) CaliHists[histName] = { lhist, histName, false, true }; // histfile histname isnew isreshist
	else std::cout << " --- No such hist found " << lhist << std::endl;

}//<<>>void KUCMS_TimeCalibration::load2DResHist( std::string histName )

void KUCMS_TimeCalibration::doResTimeFit( std::string histName ){

    std::cout << "Making Resolution Time Fits  : " << histName << std::endl;
	if( CaliHists.find(histName) != CaliHists.end() ){
		if( CaliHists[histName].h2f != NULL ){
    		kucms_SigmaFitResult results = runTimeFitter( CaliHists[histName].h2f );
		} else std::cout << " -- No H2F present !! " << std::endl;
	} else std::cout << " -- No Such hist entry !! " << std::endl;
	if( cali2DResTFile != NULL ) cali2DResTFile->Close();

}//<<>>void KUCMS_TimeCalibration::doResTimeFit( std::string histName )

void KUCMS_TimeCalibration::makeSmearTag( std::string sourceName, std::string destName, std::string smearTag ){

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

	//SmearTagSet[smearTag] = { noise, stoch, stant };

	smearTFile->cd();
	Results->Write( Results->GetName(), TObject::kOverwrite );

	delete Results;
    delete form;
    delete fit;
    smearTFile->Close();
	sourceTFile->Close();
	destTFile->Close();

}//<<>>void KUCMS_TimeCalibration::makeSmearTag( std::string histName )


bool KUCMS_TimeCalibration::plotMeanRunTimeEGR( std::string inputFileName, int srun, int erun ){

	std::cout << "Creating calibrated and uncalibrated TT mean-by-run histograms from EgammaRes Ntuples"
	          << std::endl;
	if( erun < srun ){
		std::cerr << "ERROR: invalid requested run range " << srun << " to " << erun << std::endl;
		return false;
	}//<<>>if( invalid requested range )

	struct InputSpec {
		std::string subdir;
		std::string matchstr;
		int startRun;
		int endRun;
		std::string tag;
	};
	std::vector<InputSpec> inputs;
	std::ifstream infilelist(inputFileName);
	if( not infilelist.is_open() ){
		std::cerr << "ERROR: could not open input file list: " << inputFileName << std::endl;
		return false;
	}//<<>>if( input list unavailable )
	std::string infilestr;
	unsigned int lineNumber = 0;
	while( std::getline(infilelist,infilestr) ){
		lineNumber++;
		std::stringstream ss(infilestr);
		std::string subdir;
		if( not (ss >> subdir) ) continue;
		if( subdir[0] == '#' ) continue;
		InputSpec input;
		input.subdir = subdir;
		if( not (ss >> input.matchstr >> input.startRun >> input.endRun >> input.tag) ){
			std::cerr << "ERROR: malformed mean-run input row " << lineNumber << ": " << infilestr << std::endl;
			return false;
		}//<<>>if( malformed input row )
		if( input.endRun < input.startRun ){
			std::cerr << "ERROR: invalid run range on input row " << lineNumber << std::endl;
			return false;
		}//<<>>if( invalid input range )
		if( input.endRun < srun || input.startRun > erun ) continue;
		inputs.push_back(input);
	}//<<>>while( input rows )
	if( inputs.empty() ){
		std::cerr << "ERROR: input list contains no rows overlapping the requested run range." << std::endl;
		return false;
	}//<<>>if( no relevant input rows )

	bool makeCalibrated = caliMapsValid;
	if( not makeCalibrated )
		std::cout << "WARNING: stored calibration maps are invalid; producing uncalibrated histograms only."
		          << std::endl;
	for( const auto& input : inputs ){
		if( CaliRunMapSet.find(input.tag) == CaliRunMapSet.end() ||
		    TTCaliRunMapSet.find(input.tag) == TTCaliRunMapSet.end() ){
			std::cout << "WARNING: no complete crystal/TT calibration for tag " << input.tag
			          << "; producing uncalibrated histograms only." << std::endl;
			makeCalibrated = false;
		}//<<>>if( calibration tag missing )
	}//<<>>for( input rows )

	struct FastDetInfo { uInt ttId; bool isEB; };
	std::unordered_map<uInt,FastDetInfo> fastDetInfo;
	fastDetInfo.reserve(DetIDMap.size());
	for( const auto& det : DetIDMap ){
		const auto& info = det.second;
		uInt ttId = 0;
		if( info.ecal == ECAL::EB ){
			const int ttPhi = 1 + (info.i1 - 1)/5;
			const int ttEta = 1 + (std::abs(info.i2) - 1)/5;
			ttId = info.i2 < 0 ? (ttPhi + ttEta*100) + 2000 : ttPhi + ttEta*100;
		} else ttId = 1000*(info.ecal + 3) + info.mod;
		fastDetInfo.emplace(det.first,FastDetInfo{ttId,info.ecal == ECAL::EB});
	}//<<>>for( detector IDs )

	struct TTMoments {
		double rawSum = 0.0;
		double calibratedSum = 0.0;
		unsigned long long count = 0;
	};
	using TTAccumulator = std::unordered_map<unsigned long long,TTMoments>;
	TTAccumulator ebMoments;
	TTAccumulator eeMoments;
	auto makeTTKey = []( uInt run, uInt ttId ){
		return (static_cast<unsigned long long>(run) << 32) | static_cast<unsigned long long>(ttId);
	};
	auto findRange = []( const auto& mapSet, const std::string& tag, int run )
		-> const kucms_CaliRunClass* {
		auto tagSet = mapSet.find(tag);
		if( tagSet == mapSet.end() ) return nullptr;
		auto range = tagSet->second.upper_bound(run);
		if( range == tagSet->second.begin() ) return nullptr;
		--range;
		return run <= range->second.endRun ? &range->second : nullptr;
	};

	const std::string treename("tree/llpgtree");
	const std::string timeBranchName = isCC && doUnCC ? "rhCCTime" : "rhRtTime";
	const std::string eosHost("cmseos.fnal.gov");
	std::map<std::string,std::vector<std::string>> eosFileCache;
	std::unordered_set<std::string> claimedFiles;
	unsigned long long rejectedSentinelTimes = 0;
	unsigned long long rejectedNonFiniteTimes = 0;
	unsigned long long selectedRechits = 0;
	unsigned int inputSetsProcessed = 0;

	for( const auto& input : inputs ){
		std::cout << "Processing mean-run input: " << input.subdir << " " << input.matchstr
		          << " runs " << input.startRun << " to " << input.endRun << std::endl;
		const std::string scanDir = eosDir + inDir + input.subdir;
		auto cachedFiles = eosFileCache.find(scanDir);
		if( cachedFiles == eosFileCache.end() ){
			std::cout << "Finding ROOT files under EOS directory: " << scanDir << std::endl;
			cachedFiles = eosFileCache.emplace(scanDir,findEOSRootFiles(eosHost,scanDir)).first;
		} else std::cout << "Using cached EOS listing for: " << scanDir << std::endl;
		std::vector<std::string> files;
		for( const auto& path : cachedFiles->second )
			if( input.matchstr.empty() || path.find(input.matchstr) != std::string::npos ) files.push_back(path);
		std::cout << "Found " << files.size() << " candidate ROOT files." << std::endl;

		TChain fInTree(treename.c_str());
		int filesAdded = 0;
		int badFiles = 0;
		for( const auto& tfilename : files ){
			TFile* testFile = TFile::Open(tfilename.c_str(),"READ");
			if( not testFile || testFile->IsZombie() || not testFile->IsOpen() ){
				std::cerr << "ERROR: bad ROOT file, skipping: " << tfilename << std::endl;
				if( testFile ){ testFile->Close(); delete testFile; }
				badFiles++;
				continue;
			}//<<>>if( bad ROOT file )
			TTree* testTree = nullptr;
			testFile->GetObject(treename.c_str(),testTree);
			if( not testTree ){
				std::cerr << "ERROR: missing tree " << treename << ", skipping: " << tfilename << std::endl;
				testFile->Close();
				delete testFile;
				badFiles++;
				continue;
			}//<<>>if( tree missing )
			testFile->Close();
			delete testFile;
			if( not claimedFiles.insert(tfilename).second ){
				std::cerr << "ERROR: input ROOT file matched more than one row: " << tfilename << std::endl;
				return false;
			}//<<>>if( duplicate input file )
			if( fInTree.Add(tfilename.c_str()) > 0 ) filesAdded++;
			else{
				std::cerr << "ERROR: failed to add ROOT file to TChain: " << tfilename << std::endl;
				badFiles++;
			}//<<>>if( chain add failed )
		}//<<>>for( candidate files )
		std::cout << "[INFO] filesAdded=" << filesAdded << " badFiles=" << badFiles << std::endl;
		if( filesAdded == 0 ){
			std::cerr << "ERROR: no valid ROOT files for " << input.subdir << " " << input.matchstr << std::endl;
			return false;
		}//<<>>if( no input files )
		const Long64_t nEntries = fInTree.GetEntries();
		if( nEntries <= 0 ){
			std::cerr << "ERROR: mean-run input chain has no entries." << std::endl;
			return false;
		}//<<>>if( empty input chain )

		const std::vector<std::string> requiredBranches = {
			"run","rhID",timeBranchName,"rhEnergy","rhisGS6","rhisGS1"
		};
		for( const auto& branch : requiredBranches ){
			if( not fInTree.GetBranch(branch.c_str()) ){
				std::cerr << "ERROR: required branch is missing: " << branch << std::endl;
				return false;
			}//<<>>if( branch missing )
		}//<<>>for( required branches )
		if( not isCC && not fInTree.GetBranch("rhTimeError") ){
			std::cerr << "ERROR: required branch is missing: rhTimeError" << std::endl;
			return false;
		}//<<>>if( time error branch missing )

		uInt run = 0;
		std::vector<unsigned int>* rhID = nullptr;
		std::vector<float>* rhTime = nullptr;
		std::vector<float>* rhEnergy = nullptr;
		std::vector<bool>* rhisGS6 = nullptr;
		std::vector<bool>* rhisGS1 = nullptr;
		std::vector<float>* rhTimeError = nullptr;
		fInTree.SetBranchStatus("*",0);
		fInTree.SetBranchStatus("run",1);
		fInTree.SetBranchStatus("rhID",1);
		fInTree.SetBranchStatus(timeBranchName.c_str(),1);
		fInTree.SetBranchStatus("rhEnergy",1);
		fInTree.SetBranchStatus("rhisGS6",1);
		fInTree.SetBranchStatus("rhisGS1",1);
		if( not isCC ) fInTree.SetBranchStatus("rhTimeError",1);
		if( fInTree.SetBranchAddress("run",&run) < 0 ||
		    fInTree.SetBranchAddress("rhID",&rhID) < 0 ||
		    fInTree.SetBranchAddress(timeBranchName.c_str(),&rhTime) < 0 ||
		    fInTree.SetBranchAddress("rhEnergy",&rhEnergy) < 0 ||
		    fInTree.SetBranchAddress("rhisGS6",&rhisGS6) < 0 ||
		    fInTree.SetBranchAddress("rhisGS1",&rhisGS1) < 0 ||
		    (not isCC && fInTree.SetBranchAddress("rhTimeError",&rhTimeError) < 0) ){
			std::cerr << "ERROR: failed to bind one or more mean-run branches." << std::endl;
			return false;
		}//<<>>if( branch binding failed )
		fInTree.SetCacheSize(100LL*1024LL*1024LL);
		for( const auto& branch : requiredBranches ) fInTree.AddBranchToCache(branch.c_str(),true);
		if( not isCC ) fInTree.AddBranchToCache("rhTimeError",true);
		fInTree.StopCacheLearningPhase();

		std::cout << "Starting mean-run entry loop for " << nEntries << " entries" << std::endl;
		int nextProgressPercent = 0;
		Long64_t nextProgressEntry = 0;
		int cachedRun = -1;
		const kucms_CaliRunClass* crystalRange = nullptr;
		const kucms_CaliRunClass* ttRange = nullptr;
		for( Long64_t centry = 0; centry < nEntries; centry++ ){
			if( centry >= nextProgressEntry ){
				const int completedPercent = static_cast<int>(100.0L*centry/nEntries);
				const int reportPercent = std::min(90,(completedPercent/10)*10);
				reportEventLoopProgress(centry,nEntries,reportPercent);
				nextProgressPercent = reportPercent + 10;
				nextProgressEntry = progressEntryForPercent(nEntries,nextProgressPercent);
			}//<<>>if( progress update )
			if( fInTree.GetEntry(centry) <= 0 ){
				std::cerr << "ERROR: failed to read mean-run entry " << centry << std::endl;
				return false;
			}//<<>>if( entry read failed )
			if( run < static_cast<uInt>(srun) || run > static_cast<uInt>(erun) ||
			    run < static_cast<uInt>(input.startRun) || run > static_cast<uInt>(input.endRun) ) continue;
			if( not rhID || not rhTime || not rhEnergy || not rhisGS6 || not rhisGS1 ||
			    (not isCC && not rhTimeError) ){
				std::cerr << "ERROR: null rechit vector at mean-run entry " << centry << std::endl;
				return false;
			}//<<>>if( null branch vector )
			const std::size_t nRecHits = rhID->size();
			if( rhTime->size() != nRecHits || rhEnergy->size() != nRecHits ||
			    rhisGS6->size() != nRecHits || rhisGS1->size() != nRecHits ||
			    (not isCC && rhTimeError->size() != nRecHits) ){
				std::cerr << "ERROR: rechit branch size mismatch at mean-run entry " << centry << std::endl;
				return false;
			}//<<>>if( vector size mismatch )

			if( cachedRun != static_cast<int>(run) ){
				cachedRun = run;
				crystalRange = makeCalibrated ? findRange(CaliRunMapSet,input.tag,run) : nullptr;
				ttRange = makeCalibrated ? findRange(TTCaliRunMapSet,input.tag,run) : nullptr;
				if( makeCalibrated && (not crystalRange || not ttRange || crystalRange->meanMap.empty() ||
				                       ttRange->meanMap.empty()) ){
					std::cout << "WARNING: complete calibration maps are unavailable for tag " << input.tag
					          << " run " << run
					          << "; calibrated histograms will not be produced." << std::endl;
					makeCalibrated = false;
					crystalRange = nullptr;
					ttRange = nullptr;
				}//<<>>if( calibration IOV missing )
			}//<<>>if( run changed )

			for( std::size_t idx = 0; idx < nRecHits; idx++ ){
				if( (*rhisGS6)[idx] || (*rhisGS1)[idx] ) continue;
				const float energy = (*rhEnergy)[idx];
				if( not std::isfinite(energy) || energy < 5.0f ) continue;
				const uInt id = (*rhID)[idx];
				auto det = fastDetInfo.find(id);
				if( det == fastDetInfo.end() ){
					std::cerr << "ERROR: unknown ECAL detID in mean-run input: " << id << std::endl;
					return false;
				}//<<>>if( unknown detector ID )
				if( not isCC ){
					const float timeError = (*rhTimeError)[idx];
					const float errorThreshold = det->second.isEB ? 0.61f : 1.01f;
					if( not std::isfinite(timeError) || timeError <= errorThreshold || timeError > 900.0f ) continue;
				}//<<>>if( time-error selection enabled )
				const float rawTime = (*rhTime)[idx];
				if( not std::isfinite(rawTime) ){
					rejectedNonFiniteTimes++;
					continue;
				}//<<>>if( non-finite raw time )
				if( rawTime <= -29.0f ){
					rejectedSentinelTimes++;
					continue;
				}//<<>>if( sentinel raw time )

				double calibratedTime = rawTime;
				if( makeCalibrated ){
					double correction = 0.0;
					auto crystalMean = crystalRange->meanMap.find(id);
					if( crystalMean != crystalRange->meanMap.end() ) correction += crystalMean->second;
					auto ttMean = ttRange->meanMap.find(det->second.ttId);
					if( ttMean != ttRange->meanMap.end() ) correction += ttMean->second;
					if( not std::isfinite(correction) ){
						std::cerr << "ERROR: non-finite calibration for detID " << id << " run " << run << std::endl;
						return false;
					}//<<>>if( invalid correction )
					calibratedTime -= correction;
				}//<<>>if( calibrated output enabled )
				auto& accumulator = det->second.isEB ? ebMoments : eeMoments;
				auto& moments = accumulator[makeTTKey(run,det->second.ttId)];
				moments.rawSum += rawTime;
				if( makeCalibrated ) moments.calibratedSum += calibratedTime;
				moments.count++;
				selectedRechits++;
			}//<<>>for( rechits )
		}//<<>>for( entries )
		reportEventLoopProgress(nEntries,nEntries,100);
		std::cout << "[INFO] Finished mean-run batch for " << input.subdir << " " << input.matchstr
		          << " entries=" << nEntries << std::endl;
		inputSetsProcessed++;
	}//<<>>for( input rows )

	if( inputSetsProcessed == 0 ){
		std::cerr << "ERROR: no mean-run input sets were processed." << std::endl;
		return false;
	}//<<>>if( no processed inputs )
	std::cout << "Rejected invalid raw times: sentinel=" << rejectedSentinelTimes
	          << " non-finite=" << rejectedNonFiniteTimes << std::endl;
	std::cout << "Selected gain-1 rechits for TT means: " << selectedRechits << std::endl;

	const int nRunBins = erun - srun + 1;
	const double runLowEdge = static_cast<double>(srun) - 0.5;
	const double runHighEdge = static_cast<double>(erun) + 0.5;
	constexpr unsigned long long minimumTTCount = 100;
	auto writeMeanHistograms = [&]( bool calibrated ){
		const std::string suffix = calibrated ? "_cali" : "_nocali";
		const std::string baseName = "meanRunTime_" + std::to_string(srun) + "_to_" +
		                             std::to_string(erun) + suffix;
		const std::string ebName = baseName + "_allTT_eb";
		const std::string eeName = baseName + "_allTT_ee";
		TH2F ebHist(ebName.c_str(),(ebName + ";Run;EB Mean TT Time [ns];TT count").c_str(),
		            nRunBins,runLowEdge,runHighEdge,500,-5.0,5.0);
		TH2F eeHist(eeName.c_str(),(eeName + ";Run;EE Mean TT Time [ns];TT count").c_str(),
		            nRunBins,runLowEdge,runHighEdge,500,-5.0,5.0);
		ebHist.SetDirectory(nullptr);
		eeHist.SetDirectory(nullptr);
		auto fill = [&]( const TTAccumulator& accumulator, TH2F& hist ){
			for( const auto& entry : accumulator ){
				const auto& moments = entry.second;
				// More than 100 selected rechits are required to suppress low-statistics TT means.
				if( moments.count <= minimumTTCount ) continue;
				const uInt run = static_cast<uInt>(entry.first >> 32);
				const double sum = calibrated ? moments.calibratedSum : moments.rawSum;
				const double mean = sum/static_cast<double>(moments.count);
				if( std::isfinite(mean) ) hist.Fill(run,mean);
			}//<<>>for( TT/run accumulators )
		};
		fill(ebMoments,ebHist);
		fill(eeMoments,eeHist);
		const std::string outputName = caliFileDir + baseName + "_Hist.root";
		TFile* output = TFile::Open(outputName.c_str(),"RECREATE");
		if( not output || output->IsZombie() || not output->IsOpen() ){
			std::cerr << "ERROR: could not create mean-run output file: " << outputName << std::endl;
			if( output ){ output->Close(); delete output; }
			return false;
		}//<<>>if( output unavailable )
		output->cd();
		ebHist.Write();
		eeHist.Write();
		output->Close();
		delete output;
		std::cout << "Wrote mean-run histograms: " << outputName << std::endl;
		return true;
	};

	if( not writeMeanHistograms(false) ) return false;
	if( makeCalibrated ){
		if( not writeMeanHistograms(true) ) return false;
	} else std::cout << "Calibrated mean-run histograms were not produced because complete calibration was unavailable."
	                 << std::endl;
	return true;

}//<<>>bool KUCMS_TimeCalibration::plotMeanRunTimeEGR( std::string inputFileName, int srun, int erun )

void KUCMS_TimeCalibration::makeTTDiffMaps( int srun, int erun ){

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

}//<<>>void KUCMS_TimeCalibration::makeTTDiffMaps()

void KUCMS_TimeCalibration::makeTTDriftMaps( std::string tag, int srun, int erun ){

	bool etaMod = true;
	int etaModFac = 6;
	bool phiMod = true;
	int phiModFac = 12;

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
			if( etaMod && ieta%etaModFac != 0 ) continue;
			for( int iphi = 1; iphi < 74; iphi++ ){
				if( phiMod && iphi%phiModFac != 0 ) continue;
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
					if( phiMod && iphi%phiModFac != 0 ) continue;
                    int i1 = ieta - 18;
                    if( i1 == 0 ) continue;
					if( etaMod && ieta%etaModFac != 0 ) continue;
					uInt detid = getInvTTId( iphi, i1, true );
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
									{12,kAzure+7},{13,kMagenta-4},{14,kRed-4},{15,kYellow-4},{16,kGreen-4},{17,kCyan-4},
									{16,kBlue-4},{19,kViolet+7},{20,kPink+7},{21,kOrange+7},{22,kSpring+7},{23,kTeal+7},
                                    {24,kTeal-6},{25,kAzure-6},{26,kMagenta-6},{27,kRed-6},{28,kYellow-6},{29,kGreen-6},
                                    {30,kCyan-6},{31,kBlue-6},{32,kViolet-6},{33,kPink-6},{34,kOrange-6},{35,kSpring-6},
                                    {36,kSpring+4},{37,kTeal+4},{38,kAzure+4},{39,kMagenta+4},{40,kRed+4},{41,kYellow+4},
                                    {42,kGreen+4},{43,kCyan+4},{44,kBlue+4},{45,kViolet+4},{46,kPink+4},{47,kOrange+4}
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
			if( etaMod && ieta%etaModFac != 0 ) continue;
            for( int iphi = sPhi; iphi < ePhi; iphi++ ){
				if( phiMod && iphi%phiModFac != 0 ) continue;
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
        mg->GetYaxis()->SetRangeUser(-0.4,0.1);
        //mg->GetYaxis()->SetRangeUser(-1.75,1.75);
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

}//<<>>void KUCMS_TimeCalibration::makeTTDriftMaps()

#endif
//-------------------------------------------------------------------------------------------------------------------
