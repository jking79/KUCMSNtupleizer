#include "../include/KUCMS_TimeCalibration.hh"

#include <array>
#include <functional>
#include <limits>
#include <memory>
#include <unordered_map>
#include <unordered_set>

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

namespace {

  Long64_t progressEntryForPercent( Long64_t totalEntries, int percent ){
    const Long64_t whole = totalEntries/100;
    const Long64_t remainder = totalEntries%100;
    return whole*percent + ( remainder*percent + 99 )/100;
  }//<<>>Long64_t progressEntryForPercent( Long64_t totalEntries, int percent )

  void reportEventLoopProgress( Long64_t processed, Long64_t totalEntries, int percent ){
    std::time_t now = std::time(nullptr);
    std::string curtime = std::ctime(&now);
    if( not curtime.empty() ) curtime.pop_back();
    std::cout << "Processed " << processed << " of " << totalEntries << " entries ("
    << percent << "%) at " << curtime << std::endl;
  }//<<>>void reportEventLoopProgress( Long64_t processed, Long64_t totalEntries, int percent )

}//<<>>namespace

#ifndef KUCMS_TimeCalibrationClass_cc
#define KUCMS_TimeCalibrationClass_cc

// Calibration workflow overview:
//   1. load detector-ID and luminosity/run-range inputs;
//   2. construct calibration IOVs and accumulate crystal/tower statistics;
//   3. write calibration maps and optional external-map products;
//   4. derive timing-resolution parameters and smearing configurations.

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
  caliHistFileNames["r3_p26unc"] = "cali_p26_UnCor_HistsTFile.root";

  useEosRootCali =  doLocalCali ? false : true;
  caliMapsValid = true;
  ////////useEosRootCali = false;
  eosCaliPath = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/caliconfig/";

  for( auto name : caliHistFileNames ){ caliTFile[name.first] = NULL; }

  curTag = "EG_EOY_MINI";
  // PD: EG_EOY_MINI    MC: RunIIFall17DRPremix
  curTTIov = "r2ul";
  curXIov = "prompt";
  autoSelectGain1Iovs = true;
  curLumiTag = "r2ul";

  //resTag = "default";  // NOTE: resolution used for PDs in set in master file lists, this is targ resoltion of MC for smearing
  resTag = "r2_ul18D";

  eosDir = "root://cmseos.fnal.gov//store/user/jaking/";
  inDir = "";

  lowEnergy = true;
  useEffEnergy = false;
  doGPRHists = false;
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
  doMC = false;

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

  // Define the era/tag boundaries used to translate luminosity and run
  // records into the calibration intervals consumed by later stages.

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

  promptIovMap[398904] = 401629;

  // 2026
  // lumi 401630 403867

  //Commissioning2026, Run2026A, Run2026B  400454 402229
  //Run2026B 402230 402513
  //Run2026C 402514 403422
  //Run2026D 403423 404067

  promptIovMap[401630] = 402513;
  promptIovMap[402514] = 403422;
  promptIovMap[403423] = 403867;
  promptIovMap[403868] = 404067;

  // 2026

  // END
  promptIovMap[404068] = 999999;
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
  std::string lumi26Config = "brilcalc_lumi_2026.txt";

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
  ReadLumiFile( caliFileDir+lumi26Config, r3LumiTag );

  setLumiTag(r3LumiTag);
  SetupIovMap( r3TagTT, minttlumi );
  SetupIovMap( r3TagX, minxlumi );
  //setXIov("r3x");

}//<<>>void KUCMS_TimeCalibration::SetupIovMaps()

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
    std::vector<int> tthardruns{ 284045, 306461, 325173, 357501, 359001, 359560, 359567, 359574, 359690, 359775, 359870, 359871, 359880, 359895, 359902, 360129, 360131, 360260, 360331, 362761, 373862, 387742, 398904, 401630, 402514, 403868 };
    std::vector<int> xhardruns{ 253984, 284045, 306461, 327239, 362761, 373862, 387742, 398904, 401630, 402514, 403868 };
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

}//<<>>void KUCMS_TimeCalibration::SetupIovMap( std::string tag, float minLumi )

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

}//<<>>void KUCMS_TimeCalibration::ReadLumiFile( std::string lumifile, std::string tag )

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

}//<<>>void KUCMS_TimeCalibration::ReadCaliRunFile()

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

}//<<>>bool KUCMS_TimeCalibration::SaveCaliRunFile()

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

}//<<>>bool KUCMS_TimeCalibration::SaveTTRunFile()

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

}//<<>>void KUCMS_TimeCalibration::SetResParamters( std::string tag, ... )

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
  bool success = true;
  std::unordered_set<std::string> requiredPrefixes;
  if( not makeNew ){
    for( const auto& tagEntry : TTCaliRunMapSet ){
      for( const auto& range : tagEntry.second ){
        if( range.second.histMapName != "none" ) requiredPrefixes.insert(getPrefix(tagEntry.first));
      }//<<>>for( TT ranges )
    }//<<>>for( TT tags )
    for( const auto& tagEntry : CaliRunMapSet ){
      for( const auto& range : tagEntry.second ){
        if( not range.second.isExternal && range.second.histMapName != "none" )
        requiredPrefixes.insert(getPrefix(tagEntry.first));
      }//<<>>for( crystal ranges )
    }//<<>>for( crystal tags )
  }//<<>>if( existing maps are expected )
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
      if( makeNew || requiredPrefixes.find(filename.first) != requiredPrefixes.end() ) success = false;
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
        if( not makeNew && ttfilename != "none" ) success = false;
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
            else if( not makeNew ){
              std::cerr << "ERROR: required TT calibration histogram not found: " << filename << std::endl;
              success = false;
            }//<<>>else if( not makeNew )
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
        if( not makeNew && not calirunsct.second.isExternal && xtfilename != "none" ) success = false;
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
            else if( not makeNew ){
              std::cerr << "ERROR: required crystal calibration histogram not found: " << filename << std::endl;
              success = false;
            }//<<>>else if( not makeNew )
            //std::cout << " Found " << hist << std::endl;
            if( hist ) CaliHists[filename] = { hist, filename, false, false }; // histfile histname isnew isreshist
          }//<<>>if( TTMaps.find(calirunmap.second.TTCaliMapName) == TTMaps.end() )
          else { std::cout << " Duplicate " << filename << std::endl; }
        }//<<>>for( auto ebeehist : ebeemap )
      }//<<>>for( auto mapname : mtype )
    }//<<>>for( auto& calirunsct : calirunmap )
  }//<<>>for( auto& calirunmap : CaliRunMapSet )

  if( success ) success = makeCaliMaps();

  if( not stayOpen ){

    for( auto& calihist : CaliHists ){
      if( calihist.second.h2f != NULL ){ delete calihist.second.h2f; }
      //if( calihist.second.h1f != NULL ) delete calihist.second.h1f;
    }//<<>>for( auto& calimapsct : TTCaliMaps )
    CaliHists.clear();
    //caliTFile->Close();
    //caliR3TFile->Close();
    for( auto& califile : caliTFile ){
      if( califile.second ){
        califile.second->Close();
        delete califile.second;
        califile.second = nullptr;
      }//<<>>if( file exists )
    }//<<>>for( calibration files )

  }//<<>>if( not stayOpen )
  caliMapsValid = success;
  return success;

}//<<>>void KUCMS_TimeCalibration::SetupCaliHists()

//    TH2F* h2f;
//    std::string histName;
//    bool isNew;

// for updating and/or adding new info to califiles - do last - saving closes TCalifile and clears CaliHists
bool KUCMS_TimeCalibration::SaveCaliHists(){

  // Persist the accumulated calibration histograms and close files that were
  // opened for update.  The caller uses the return value to detect write errors.

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

}//<<>>bool KUCMS_TimeCalibration::SaveCaliHists()

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

  // Convert accumulated run-range statistics into EB/EE crystal and trigger-
  // tower histograms.  These histograms are the inputs to the map writers.

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
          std::cerr << "ERROR: " << entry.first << " not found in DetIDMap while saving crystal maps."
          << std::endl;
          continue;
        }//<<>>if( DetIDMap.find(entry.first) == DetIDMap.end() )
        const kucms_DetIDStruct& idinfo = DetIDMap.at(entry.first);
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
        auto det = DetIDMap.find(entry.first);
        if( det == DetIDMap.end() ){
          std::cerr << "ERROR: " << entry.first << " not found in DetIDMap while saving crystal means."
          << std::endl;
          continue;
        }//<<>>if( detector ID missing )
        const kucms_DetIDStruct& idinfo = det->second;
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
        auto det = DetIDMap.find(entry.first);
        if( det == DetIDMap.end() ){
          std::cerr << "ERROR: " << entry.first << " not found in DetIDMap while saving crystal errors."
          << std::endl;
          continue;
        }//<<>>if( detector ID missing )
        const kucms_DetIDStruct& idinfo = det->second;
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

  // Materialize calibration maps from the accumulated run-range records,
  // applying the configured IOV boundaries and detector-ID conversions.

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
      uInt detid = getDetIdInfo(iphi,i1,ECAL::EB);
      if( detid == 0 ) continue;
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

}//<<>>void KUCMS_TimeCalibration::LoadExtCali( std::string calihist, std::string mapname, std::string tag, int startr, int endr )

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

}//<<>>uInt KUCMS_TimeCalibration::getTTId( uInt detId, int gainID )

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

}//<<>>uInt KUCMS_TimeCalibration::getInvTTId( int i1, int i2, bool isEB )

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

const kucms_DetIDStruct& KUCMS_TimeCalibration::getDetIdInfo( uInt rhid ){

  auto det = DetIDMap.find(rhid);
  if( det != DetIDMap.end() ) return det->second;
  std::cerr << "ERROR: unknown ECAL detID requested: " << rhid << std::endl;
  static const kucms_DetIDStruct invalidDetID;
  return invalidDetID;

}//<<>>const kucms_DetIDStruct& KUCMS_TimeCalibration::getDetIdInfo( uInt rhid )

uInt KUCMS_TimeCalibration::getDetIdInfo( int i1, int i2, int ecal ){

  auto first = InvDetIDMap.find(i1);
  if( first != InvDetIDMap.end() ){
    auto second = first->second.find(i2);
    if( second != first->second.end() ){
      auto region = second->second.find(ecal);
      if( region != second->second.end() ) return region->second;
    }//<<>>if( second coordinate exists )
  }//<<>>if( first coordinate exists )
  std::cerr << "ERROR: no ECAL detID for coordinates (" << i1 << "," << i2
  << ") in region " << ecal << std::endl;
  return 0;

}//<<>>uInt KUCMS_TimeCalibration::getDetIdInfo( int i1, int i2, int ecal )

std::map<int,std::map<int,std::map<int,uInt>>> KUCMS_TimeCalibration::getInvDetIDMap(){

  return InvDetIDMap;

}//<<>>std::map<int,std::map<int,std::map<int,uInt>>> KUCMS_TimeCalibration::getInvDetIDMap()

std::map<UInt_t,kucms_DetIDStruct> KUCMS_TimeCalibration::getDetIDMap(){; // map of information by detid for each crystal

  return DetIDMap;

}//<<>>std::map<UInt_t,kucms_DetIDStruct> KUCMS_TimeCalibration::getDetIDMap()

//std::map<std::string,std::map<int,CaliRunStruct>> CaliRunMapSet;

float KUCMS_TimeCalibration::getCalibration( uInt rhid, int run, const std::string& tag, int gainID ){

  if( rhid == 0 ) return -999.f; // valid rechit id check
  //if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
  //if( not isEB ){ std::cout << "XCalibration for EE is not supported." << std::endl; return 0.f; }
  //if( not validCurrentTag ){ std::cout << "No current tag set." << std::endl; return 0.f; }
  //if( DetIDMap[rhid].ecal != ECAL::EB ) return 0.f;
  if( tag == "None" ) return 0.f;
  externalCali = false;
  if( externalCali ){ return 0.f; std::cout << " Called getCalibration() while using external cali? " << std::endl; }
  auto tagSet = CaliRunMapSet.find(tag);
  auto ttTagSet = TTCaliRunMapSet.find(tag);
  if( tagSet == CaliRunMapSet.end() or ttTagSet == TTCaliRunMapSet.end() ){
    if( reportedInvalidCalibrationTags.insert(tag).second ){
      std::cerr << "ERROR: calibration tag '" << tag << "' is not available in "
      << ( tagSet == CaliRunMapSet.end() ? "the crystal map" : "" )
      << ( tagSet == CaliRunMapSet.end() and ttTagSet == TTCaliRunMapSet.end() ? " or " : "" )
      << ( ttTagSet == TTCaliRunMapSet.end() ? "the trigger-tower map" : "" )
      << "; no calibration will be applied." << std::endl;
    }//<<>>if( first report for this tag )
    return 0.f;
  }//<<>>if( invalid calibration tag )
  double xtaltime = -1000.f;
  //////float hgcali = ( externalCali ) ? 0.f : ( gainID == 1 ) ? 0.f : getTTCali( rhid, run, tag, 2 );
  double hgcali = 0;
  uInt gsid = ( gainID > 1 ) ? getGSRhId( rhid, gainID ) : 0;
  if( gainID > 1 and ( gsid < 99 or gsid > 300 ) ) std::cout << " Bad GS Id with " << rhid << " gid " << gainID << std::endl;
  double ttcali = getTTCali( rhid, run, tag, 1 );
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
  if( xtaltime == -1000.f ){ std::cout << "XCalibration period not found for run " << run << std::endl; return 0.f; }
  //std::cout << " -- For " << rhid << " in " << run << " gainID " << gainID;
  //std::cout << " Found: " << xtaltime << " + " << ttcali << " + " << hgcali << std::endl;
  return xtaltime + ttcali + hgcali;
  //return ttcali;
  //if( gainID == 1 ) return xtaltime + ttcali + hgcali;
  //else return 0.f;

}//<<>>float KUCMS_TimeCalibration::getCalibration( uInt rhid, int run, const std::string& tag, int gainID )

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

}//<<>>float KUCMS_TimeCalibration::getTTCali( uInt rhid, int run, const std::string& tag, int gainID )

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

}//<<>>float KUCMS_TimeCalibration::getSmearedTime( float mtime, float ampres, std::string stag )

float KUCMS_TimeCalibration::getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int crun, std::string ctag, std::string stag ){

  // This combined calibration-plus-smearing path is currently a placeholder;
  // callers receive the historical sentinel until the workflow is implemented.
  return -99.f;

}//<<>>float KUCMS_TimeCalibration::getSmrdCalibTime( float rhtime, float rhamp, uInt rhid, int crun, std::string ctag, std::string stag )

//**************  main calibrate &/or smearing function based on tag + info ********************************************************************

// Return the configured time variance for one rechit.  The active overload
// resolves the data/MC tag and gain-dependent parameters before applying the
// fallback resolution.
float KUCMS_TimeCalibration::getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, const std::string& dataSetKey, int mctype, int gianID ){

  if( amplitude == 0 ) return 100.f;
  auto det = DetIDMap.find(rechitID);
  if( det == DetIDMap.end() ){
    std::cerr << "ERROR: unknown ECAL detID in getTimeResoltuion: " << rechitID << std::endl;
    return 100.f;
  }//<<>>if( detector ID missing )
  bool isEB( det->second.ecal == ECAL::EB );
  auto reportInvalidTag = [this]( const std::string& tag ){
    if( reportedInvalidResolutionTags.insert(tag).second ){
      std::cerr << "ERROR: resolution tag '" << tag << "' is not configured; "
      << ( tag == "default" ? "returning fallback variance 100."
      : "using the default resolution." ) << std::endl;
    }//<<>>if( first report for this tag )
  };
  float var = ( dataSetKey == "None" ) ? -1 : -99;
  if( mctype == 1 ){ // data or MC
    static const std::string r2UL18ABCTag = "r2_ul18ABC";
    static const std::string r2UL18DTag = "r2_ul18D";
    const std::string* resolvedTagName = &dataSetKey;
    if( dataSetKey == "r2_ul18" ) resolvedTagName = ( Evt_run > 320672 ) ? &r2UL18DTag : &r2UL18ABCTag;
    auto resolvedTag = ResTagSet.find(*resolvedTagName);
    if( resolvedTag != ResTagSet.end() ){
      const auto& parameters = resolvedTag->second;
      if( amplitude > 400 ){

        float res = 0;
        if( amplitude > 900 ){
          if( isEB ) res = ( gianID == 1 ) ? parameters.ebg1he3res : parameters.ebhg3res;
          else res = ( gianID == 1 ) ? parameters.eeg1he3res : parameters.eehg3res;
        } else if( amplitude > 600 ){
          if( isEB ) res = ( gianID == 1 ) ? parameters.ebg1he2res : parameters.ebhg2res;
          else res = ( gianID == 1 ) ? parameters.eeg1he2res : parameters.eehg2res;
        } else {
          if( isEB ) res = ( gianID == 1 ) ? parameters.ebg1he1res : parameters.ebhg1res;
          else res = ( gianID == 1 ) ? parameters.eeg1he1res : parameters.eehg1res;
        }//<<>>if( amplitude > 900 )
        var = sq2( res );

      }//<<>>if( amplitude > 400 )
      if( amplitude <= 400 or var == 0 ){

        float noise = isEB ? parameters.ebnoise : parameters.eenoise;
        float stoch = isEB ? parameters.ebstoch : parameters.eestoch;
        float stant = isEB ? parameters.ebstant : parameters.eestant;
        var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);

      }//<<>>if( amplitude > 400 )
    } else if( *resolvedTagName != "None" ) reportInvalidTag(*resolvedTagName);
  }//<<>>if( mctype == 0 )
  else if( mctype == 0 or mctype == 2 ){ // data or MC   --- resTag is the target resolution for smearing
    static const std::string r2UL18DTag = "r2_ul18D";
    const std::string* themctag = ( resTag == "None" ) ? &dataSetKey : &resTag;
    if( resTag == "default"){
      if( dataSetKey.substr(0,2) == "r2" ){ themctag = &r2UL18DTag; }//<<>>if( Run 2 key )
      else if( dataSetKey.substr(0,2) == "r3" ){ themctag = &r2UL18DTag; }//<<>>if( Run 3 key )
    }//<<>>if( resTag == "default")
    auto resolvedTag = ResTagSet.find(*themctag);
    if( resolvedTag != ResTagSet.end() ){
      const auto& parameters = resolvedTag->second;
      float noise = isEB ? parameters.ebnoise : parameters.eenoise;
      float stoch = isEB ? parameters.ebstoch : parameters.eestoch;
      float stant = isEB ? parameters.ebstant : parameters.eestant;
      var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);
      //std::cout << " --- getting MC res for : " << *themctag << " " << noise << " " << stoch;
      //std::cout << " " << stant << " " << var << std::endl;
    } else if( *themctag != "None" ) reportInvalidTag(*themctag);
  }//<<>>if( mctype == 0 )
  else {
    var = -9;
    if( reportedInvalidMCTypes.insert(mctype).second ){
      std::cerr << "ERROR: MC type " << mctype << " is not valid; using the default resolution." << std::endl;
    }//<<>>if( first report for this type )
  }//<<>>if( invalid MC type )
  if( var < 0 ){
    auto defaultTag = ResTagSet.find("default");
    if( defaultTag == ResTagSet.end() ){
      reportInvalidTag("default");
      return 100.f;
    }//<<>>if( default resolution missing )
    const auto& parameters = defaultTag->second;
    float noise = isEB ? parameters.ebnoise : parameters.eenoise;
    float stoch = isEB ? parameters.ebstoch : parameters.eestoch;
    float stant = isEB ? parameters.ebstant : parameters.eestant;
    var = sq2(noise/amplitude) + sq2(stoch)/amplitude + 2*sq2(stant);
    //var = isEB ? sq2(26.0/amplitude) + sq2(0.6)/amplitude + 2*sq2(0.094) : sq2(42.2/amplitude) + sq2(1.2)/amplitude + 2*sq2(0.298);
  }//<<>>if( var == -1 )
  return var;

}//<<>>float KUCMS_TimeCalibration::getTimeResoltuion( float amplitude, unsigned int rechitID, unsigned int Evt_run, const std::string& dataSetKey, int mctype, int gianID )

float KUCMS_TimeCalibration::getCorrectedTime( float time, float amplitude, unsigned int rechitID,
unsigned int Evt_run, const std::string& dataSetKey, int mctype, int gainID,
bool applyCalibration ){

  //getTimeResoltuion( amplitude, rechitID, Evt_run, dataSetKey, mctype, gainID )
  // NOTE: dataSetKey is the calibration tag in this function.
  double rtime = 0;
  //if( Evt_run > 388000 || Evt_run < 253984 ) return time; // within calibrated run range check

  if( mctype == 1 ){ // data
    if( not applyCalibration ) return time;
    rtime = time - getCalibration(rechitID,Evt_run,dataSetKey,gainID);
    return rtime;

  } else { // MC

    double ctime = 0;
    double tarvar = 0;
    double sorvar = 0;

    static const std::string r2UL18SqmcTag = "r2_ul18_sqmc";
    const std::string& resolvedDataSetKey = ( dataSetKey == "r2_ul18_mc" ) ? r2UL18SqmcTag : dataSetKey;
    if( applyCalibration )
    ctime = time - getCalibration(rechitID,Evt_run,resolvedDataSetKey,gainID);
    else{ ctime = time; }
    //std::cout << " Pre Check amplitude : " << amplitude << std::endl;
    //if( amplitude < 1 ) return ctime;
    //if( time == 0 ) return ctime;
    if( resTag == "None" ) return ctime; // --- resTag is the target resolution for smearing
    if( resolvedDataSetKey == "None" ) return ctime; // dataSetKey is the calibration used && source resolution for smearing
    //std::cout << " Post Check amplitude : " << amplitude << std::endl;

    if( resTag != "default" ){  // TARGET RESOLUTION

      tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, resTag, mctype, gainID );

    }//<<>>if( resTag != "default" )
    else if( resTag == "default" || resolvedDataSetKey.substr(0,2) == "r2" ){ // use smear tag?  -- check first for set tag !

      static const std::string tag = "r2_ul18D";
      tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

    }//<<>>if( Run 2 key or default target )
    else if( resolvedDataSetKey.substr(0,2) == "r3" ){

      static const std::string tag = "r2_ul18D";
      tarvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

    }//<<>>if( Run 3 key )

    if( ResTagSet.find(resolvedDataSetKey) != ResTagSet.end() ){ // SOURCE RESOLTUION

      static const std::string runIIUL18RecoTag = "RunIISummer20UL18RECO";
      static const std::string run3MCTag = "theRun3MC_campain";
      const std::string* sourceTag = &resolvedDataSetKey;
      if( resolvedDataSetKey == "r2_ul18_qcdmc" or resolvedDataSetKey == "r2_ul18_sqmc" ) sourceTag = &runIIUL18RecoTag;
      if( resolvedDataSetKey == "r3_mc" ) sourceTag = &run3MCTag;
      sorvar = getTimeResoltuion( amplitude, rechitID, Evt_run, *sourceTag, mctype, gainID );

    } else { // 2018 MC Resoltion values

      static const std::string tag = "default";
      sorvar = getTimeResoltuion( amplitude, rechitID, Evt_run, tag, mctype, gainID );

    }//<<>>if( ResTagSet.find(tag) != ResTagSet.end() )

    //0.00691415 1.60666 0.133913  smearing for 2017 MC -> Data ( original smearing values )
    double tresv = tarvar/2;
    double sresv = sorvar/2;
    double smvar = ( tresv > sresv ) ? tresv - sresv : 0.0;//      std::max( 0.0, tresv - sresv );
    double resolution = ( smvar > 0.0 ) ? std::sqrt(smvar) : 0.0;
    rtime = getRandom->Gaus( ctime, resolution );
    // Intentionally disabled diagnostic branch for negative smearing variance.
    if( false ){
      std::cout << "No smearing values set for this tag : " << resolvedDataSetKey << " time " << time << " c " << ctime << std::endl;
      std::cout << " smearing parameters: t " << tresv << " s " << sresv << " amp " << amplitude << std::endl;
      std::cout << " getSmeared Res: " << resolution << " rtime " << rtime << std::endl;
    }//<<>>if( false )
    return rtime;

  }//<<>>if( mctype == 0 )

  return rtime;

}//<<>>float KUCMS_TimeCalibration::getCorrectedTime( float time, float amplitude, unsigned int rechitID, unsigned int Evt_run, const std::string& dataSetKey, int mctype, int gainID, bool applyCalibration )

//*********************************************************************************************************************************************

// assume full runs are present in input files
// input file will give the name of a text file "list" of root files, the path to the root files, and the tag to use for those files
// check if run range started
// check if run range complete
//     std::map<std::string,std::map<int,int>> iovMaps; // < tag, < start run, end rn >>
//	- if started but not complete - is run needed ( lastrun )
//	- if not started - start new cali map

bool KUCMS_TimeCalibration::makeGain1CaliMapsEGR( std::string inputFileName, bool small ){

  // Build Gain-1 calibration maps from EgammaRes input, including run-range
  // bookkeeping and the EB/EE occupancy and timing products.

  std::cout << "Creating gain-1 TT and crystal calibration maps in one event pass" << std::endl;
  if( not caliMapsValid ){
    std::cerr << "ERROR: cannot process calibration input because stored calibration maps are invalid." << std::endl;
    return false;
  }//<<>>if( stored maps invalid )
  if( not autoSelectGain1Iovs && iovMaps.find(curTTIov) == iovMaps.end() ){
    std::cerr << "ERROR: TT IOV is not valid: " << curTTIov << std::endl;
    return false;
  }//<<>>if( invalid TT IOV )
  if( not autoSelectGain1Iovs && iovMaps.find(curXIov) == iovMaps.end() ){
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

    std::string ttIovTag = curTTIov;
    std::string xIovTag = curXIov;
    if( autoSelectGain1Iovs ){
      const bool isMCTag = tag == "mc" ||
      ( tag.size() >= 2 && tag.compare(tag.size()-2,2,"mc") == 0 );
      if( isMCTag ){
        ttIovTag = "mc";
        xIovTag = "mc";
      } else if( tag.rfind("r2_",0) == 0 ){
        ttIovTag = "r2ultt";
        xIovTag = "prompt";
      } else if( tag.rfind("r3_",0) == 0 ){
        ttIovTag = "r3tt";
        xIovTag = "r3x";
      } else {
        std::cerr << "ERROR: cannot automatically select gain-1 IOV maps for calibration tag "
        << tag << "; call setAutoSelectGain1Iovs(false), setTTIov(...), and "
        << "setXIov(...) for an R&D override." << std::endl;
        return false;
      }//<<>>if( recognized calibration tag family )
    }//<<>>if( automatic gain-1 IOV selection )
    if( iovMaps.find(ttIovTag) == iovMaps.end() ){
      std::cerr << "ERROR: TT IOV is not valid: " << ttIovTag << " for calibration tag " << tag << std::endl;
      return false;
    }//<<>>if( invalid TT IOV )
    if( iovMaps.find(xIovTag) == iovMaps.end() ){
      std::cerr << "ERROR: crystal IOV is not valid: " << xIovTag << " for calibration tag " << tag << std::endl;
      return false;
    }//<<>>if( invalid crystal IOV )
    const auto& ttIovs = iovMaps.at(ttIovTag);
    const auto& crystalIovs = iovMaps.at(xIovTag);
    std::cout << "Using gain-1 IOV maps TT=" << ttIovTag << " X=" << xIovTag
    << " for calibration tag " << tag << std::endl;

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
        ttRange = findOrCreateRange(TTCaliRunMapSet,tag,ttIovs,run,"TT");
        crystalRange = findOrCreateRange(CaliRunMapSet,tag,crystalIovs,run,"X");
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

  // General EgammaRes map production for crystal or trigger-tower inputs.
  // The flags select the map type, gain, reduced test mode, and calibration use.

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
        if( testFile ){ testFile->Close(); delete testFile; }
        nBadFiles++;
        std::cout << "Z";
        continue;
      }//<<>>if( !testFile || testFile->IsZombie() || !testFile->IsOpen() )

      TTree* testTree = nullptr;
      testFile->GetObject(treename.c_str(), testTree);

      if( !testTree ){
        std::cerr << "\nERROR: missing tree " << treename << ", skipping: " << tfilename << std::endl;
        testFile->Close();
        delete testFile;
        nBadFiles++;
        std::cout << "T";
        continue;
      }//<<>>if( !testTree )

      testFile->Close();
      delete testFile;
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

}//<<>>bool KUCMS_TimeCalibration::makeCaliMapsEGR( std::string inputFileName, bool doTT, int GID, bool small, bool doCali )

bool KUCMS_TimeCalibration::plot2dResolutionEGR( std::string inputFileName, bool scale, bool usecali, bool smear, std::string ext ){

  std::cout << "Creating 2D Resolution Hist from EgammaRes Ntuples " << std::endl;
  if( not usecali ) std::cout << " -- NoCali " << std::endl;
  if( smear ) std::cout << " -- Smeared " << std::endl;

  bool debug = false;
  //bool debug = true;
  bool small = false;
  //bool small = true;

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
  UInt_t          run = 0;
  std::vector<unsigned int> *resRhID = nullptr;
  std::vector<float>   *resAmp = nullptr;
  std::vector<float>   *resE = nullptr;
  std::vector<float>   *resRtTime = nullptr;
  std::vector<float>   *resCCTime = nullptr;
  std::vector<float>   *resTOF = nullptr;

  std::vector<unsigned int> *rhID = nullptr;
  std::vector<bool>    *rhisGS6 = nullptr;
  std::vector<bool>    *rhisGS1 = nullptr;

  // List of branches
  TBranch        *b_run = nullptr;   //!
  TBranch        *b_resRhID = nullptr;   //!
  TBranch        *b_resAmp = nullptr;   //!
  TBranch        *b_resE = nullptr;   //!
  TBranch        *b_resRtTime = nullptr;   //!
  TBranch        *b_resCCTime = nullptr;   //!
  TBranch        *b_resTOF = nullptr;   //!

  TBranch        *b_rhID = nullptr;   //!
  TBranch        *b_rhisGS6 = nullptr;   //!
  TBranch        *b_rhisGS1 = nullptr;   //!

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
  if( fXBins.size() < 2 || fYBins.size() < 2 ){
    std::cerr << "ERROR: invalid resolution histogram bin configuration." << std::endl;
    return false;
  }//<<>>if( invalid bins )
  const auto xbins = &fXBins[0];
  const auto ybins = &fYBins[0];
  std::unordered_map<uInt,const kucms_DetIDStruct*> fastDetInfo;
  fastDetInfo.reserve(DetIDMap.size());
  for( const auto& det : DetIDMap ) fastDetInfo.emplace(det.first,&det.second);
  std::cout << "open input files list : " << inputFileName << std::endl;

  std::ifstream infilelist(inputFileName);
  if( !infilelist.is_open() ){
    std::cerr << "ERROR: could not open input file list: " << inputFileName << std::endl;
    return false;
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
    const bool parsed = doMC ? bool(ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast >> wgt)
    : bool(ss >> subdir >> matchstr >> srun >> erun >> tag >> dolast);
    if( not parsed || subdir.empty() || matchstr.empty() || tag.empty() || srun > erun ||
    (dolast != 0 && dolast != 1) ){
      std::cerr << "ERROR: malformed resolution input row: " << infilestr << std::endl;
      return false;
    }//<<>>if( malformed row )
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
      else { std::cerr << "ERROR: bad tag for uncorrected CC times: " << tag << std::endl; return false; }
    }//<<>>if( isCC && doUnCC )

    if( smear ) std::cout << "with tags : " << tag << " and smearing " << resTag << std::endl;
    else std::cout << "with tags : " << tag << std::endl;
    std::cout << "For Run " << srun << " to Run " << erun << std::endl;
    std::cout << "Producing 2D Resolution Map " << std::endl;

    // insure we have Calimaps if we are doing Xtal calibration maps
    if( usecali && (CaliRunMapSet.find(tag) == CaliRunMapSet.end() ||
    TTCaliRunMapSet.find(tag) == TTCaliRunMapSet.end()) ){
      std::cerr << "ERROR: incomplete calibration maps for tag " << tag << std::endl;
      return false;
    }//<<>>if( calibration tag is incomplete )

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
    for( const auto& tfilename : files ){

      nLines++;
      TFile* testFile = TFile::Open(tfilename.c_str(),"READ");
      if( !testFile || testFile->IsZombie() || !testFile->IsOpen() ){
        std::cerr << "\nERROR: bad ROOT file, skipping: " << tfilename << std::endl;
        if( testFile ){ testFile->Close(); delete testFile; }
        nBadFiles++;
        std::cout << "Z";
        continue;
      }//<<>>if( bad ROOT file )

      TTree* testTree = nullptr;
      testFile->GetObject(treename.c_str(),testTree);
      if( !testTree ){
        std::cerr << "\nERROR: missing tree " << treename << ", skipping: " << tfilename << std::endl;
        testFile->Close();
        delete testFile;
        nBadFiles++;
        std::cout << "T";
        continue;
      }//<<>>if( tree missing )

      testFile->Close();
      delete testFile;
      std::cout << "-";
      const int result = fInTree->Add(tfilename.c_str());
      if( result > 0 ){ std::cout << "-"; nAdded++; }
      else std::cout << "0";

    }//<<>>for( candidate ROOT files )
    std::cout << std::endl;
    std::cout << " -- [INFO] lines=" << nLines << " filesAdded=" << nAdded << " badFiles=" << nBadFiles << std::endl;
    if( nAdded == 0 ){
      std::cerr << "ERROR: no valid ROOT files were added for input row: " << infilestr << std::endl;
      delete fInTree;
      return false;
    }//<<>>if( nAdded == 0 )

    const bool useUnCCTime = isCC && doUnCC;
    const bool needGainInfo = useGSwitch && (not doAllGain || usecali || smear);
    fInTree->SetBranchStatus("*",0);
    bool branchesValid = true;
    auto requireBranch = [&]( const char* name, auto* address, TBranch** branch ){
      if( not fInTree->GetBranch(name) ){
        std::cerr << "ERROR: required resolution input branch is missing: " << name << std::endl;
        branchesValid = false;
        return;
      }//<<>>if( branch missing )
      fInTree->SetBranchStatus(name,1);
      if( fInTree->SetBranchAddress(name,address,branch) < 0 ){
        std::cerr << "ERROR: failed to bind resolution input branch: " << name << std::endl;
        branchesValid = false;
      }//<<>>if( branch binding failed )
    };
    requireBranch("run",&run,&b_run);
    requireBranch("resRhID",&resRhID,&b_resRhID);
    requireBranch("resAmp",&resAmp,&b_resAmp);
    requireBranch("resE",&resE,&b_resE);
    requireBranch("resRtTime",&resRtTime,&b_resRtTime);
    if( useUnCCTime ) requireBranch("resCCTime",&resCCTime,&b_resCCTime);
    requireBranch("resTOF",&resTOF,&b_resTOF);
    if( needGainInfo ){
      requireBranch("rhID",&rhID,&b_rhID);
      requireBranch("rhisGS6",&rhisGS6,&b_rhisGS6);
      requireBranch("rhisGS1",&rhisGS1,&b_rhisGS1);
    }//<<>>if( gain information is required )
    if( not branchesValid ){
      delete fInTree;
      return false;
    }//<<>>if( not branchesValid )
    fInTree->SetCacheSize(100LL*1024LL*1024LL);
    for( const auto& branch : {"run","resRhID","resAmp","resE","resRtTime","resTOF"} )
    fInTree->AddBranchToCache(branch,true);
    if( useUnCCTime ) fInTree->AddBranchToCache("resCCTime",true);
    if( needGainInfo ){
      for( const auto& branch : {"rhID","rhisGS6","rhisGS1"} )
      fInTree->AddBranchToCache(branch,true);
    }//<<>>if( gain information is required )
    fInTree->StopCacheLearningPhase();

    std::cout << " Getting calibration values and plotting" << std::endl;

    auto nEntries = fInTree->GetEntries();
    if( nEntries <= 0 ){
      std::cerr << "ERROR: resolution input chain has no entries." << std::endl;
      delete fInTree;
      return false;
    }//<<>>if( nEntries <= 0 )
    if( small && nEntries > 1000000 ) nEntries = 1000000;
    if( debug ) nEntries = ( nEntries < 10000 ) ? nEntries : 10000;
    //if( debug )

    std::string caliend = usecali ? "" : nocalistr;
    if( doUnCC ) caliend += uncorstr;
    if( smear ) caliend += smearstr;
    caliend += ext;
    const std::string histMapName = "ResMap_" + std::to_string(srun) + "_" + std::to_string(erun);
    const std::string lsfnameeb = histMapName + lochist + caliend + eb_hist;
    const std::string ldfnameeb = histMapName + druhist + caliend + eb_hist;
    const std::string gbfnameeb = histMapName + globhist + caliend + eb_hist;
    const std::string efnameeb = histMapName + ehist + caliend + eb_hist;
    const std::string lsfnameee = histMapName + lochist + caliend + ee_hist;
    const std::string ldfnameee = histMapName + druhist + caliend + ee_hist;
    const std::string gbfnameee = histMapName + globhist + caliend + ee_hist;
    const std::string efnameee = histMapName + ehist + caliend + ee_hist;
    auto getOrCreateResolutionHist = [&]( const std::string& name, bool amplitudeVsEnergy ){
      auto existing = CaliHists.find(name);
      if( existing != CaliHists.end() ) return existing->second.h2f;
      TH2F* hist = amplitudeVsEnergy
      ? new TH2F(name.c_str(),name.c_str(),500,0,1000,100,0,200)
      : new TH2F(name.c_str(),name.c_str(),fXBins.size()-1,xbins,fYBins.size()-1,ybins);
      hist->GetXaxis()->SetTitle(fXTitle.c_str());
      hist->GetYaxis()->SetTitle(( amplitudeVsEnergy ? fEYTitle : fYTitle ).c_str());
      hist->GetZaxis()->SetTitle(fZTitle.c_str());
      CaliHists[name] = { hist, name, true, true };
      return hist;
    };
    TH2F* localSameEB = getOrCreateResolutionHist(lsfnameeb,false);
    TH2F* localDifferentEB = getOrCreateResolutionHist(ldfnameeb,false);
    TH2F* globalEB = getOrCreateResolutionHist(gbfnameeb,false);
    TH2F* amplitudeEnergyEB = getOrCreateResolutionHist(efnameeb,true);
    TH2F* localSameEE = getOrCreateResolutionHist(lsfnameee,false);
    TH2F* localDifferentEE = getOrCreateResolutionHist(ldfnameee,false);
    TH2F* globalEE = getOrCreateResolutionHist(gbfnameee,false);
    TH2F* amplitudeEnergyEE = getOrCreateResolutionHist(efnameee,true);
    const bool fillResolutionHists = CaliHists.find(lsfnameeb)->second.isNew;
    auto findCalibrationRange = []( const auto& mapSet, const std::string& calibrationTag, int eventRun )
    -> const kucms_CaliRunClass* {
      auto tagSet = mapSet.find(calibrationTag);
      if( tagSet == mapSet.end() ) return nullptr;
      auto range = tagSet->second.upper_bound(eventRun);
      if( range == tagSet->second.begin() ) return nullptr;
      --range;
      return eventRun <= range->second.endRun ? &range->second : nullptr;
    };
    int cachedCalibrationRun = -1;
    const kucms_CaliRunClass* cachedCrystalRange = nullptr;
    const kucms_CaliRunClass* cachedTTRange = nullptr;
    auto getCachedCalibration = [&]( uInt rechitID, int gainID, const kucms_DetIDStruct& detInfo ){
      if( not usecali ) return 0.f;
      if( cachedCalibrationRun != static_cast<int>(run) ){
        cachedCalibrationRun = run;
        cachedCrystalRange = findCalibrationRange(CaliRunMapSet,tag,run);
        cachedTTRange = findCalibrationRange(TTCaliRunMapSet,tag,run);
      }//<<>>if( run changed )
      if( not cachedCrystalRange ) return 0.f;
      double correction = 0.0;
      auto crystalMean = cachedCrystalRange->meanMap.find(rechitID);
      if( crystalMean != cachedCrystalRange->meanMap.end() ) correction += crystalMean->second;
      if( cachedTTRange ){
        uInt ttID = 0;
        if( detInfo.ecal == ECAL::EB ){
          const int ttPhi = 1 + (detInfo.i1 - 1)/5;
          const int ttEta = 1 + (std::abs(detInfo.i2) - 1)/5;
          ttID = detInfo.i2 < 0 ? (ttPhi + ttEta*100) + 2000 : ttPhi + ttEta*100;
        } else ttID = 1000*(detInfo.ecal + 3) + detInfo.mod;
        auto ttMean = cachedTTRange->meanMap.find(ttID);
        if( ttMean != cachedTTRange->meanMap.end() ) correction += ttMean->second;
      }//<<>>if( TT calibration range exists )
      if( gainID > 1 ){
        uInt gainSwitchID = 0;
        if( detInfo.ecal == ECAL::EB ) gainSwitchID = detInfo.i2 < 0
        ? ( gainID == 2 ? 206 : 201 ) : ( gainID == 2 ? 106 : 101 );
        else if( detInfo.ecal == ECAL::EM ) gainSwitchID = gainID == 2 ? 116 : 111;
        else if( detInfo.ecal == ECAL::EP ) gainSwitchID = gainID == 2 ? 226 : 221;
        auto gainMean = cachedCrystalRange->meanMap.find(gainSwitchID);
        if( gainMean != cachedCrystalRange->meanMap.end() ) correction += gainMean->second;
      }//<<>>if( gain-switched rechit )
      return static_cast<float>(correction);
    };

    std::cout << "Mf2d Proccessing " << nEntries << " entries : " << std::endl;
    Long64_t malformedResolutionEvents = 0;
    Long64_t recoveredLocalEvents = 0;
    Long64_t recoveredGlobalEvents = 0;
    Long64_t skippedShortResolutionEvents = 0;
    auto areNeighboringCrystals = []( const kucms_DetIDStruct& first, const kucms_DetIDStruct& second ){
      if( first.ecal != second.ecal || first.ecal == ECAL::NONE ) return false;
      int delta1 = std::abs(first.i1-second.i1);
      int delta2 = 0;
      if( first.ecal == ECAL::EB ){
        delta1 = std::min(delta1,360-delta1);
        const int firstEta = first.i2 > 0 ? first.i2-1 : first.i2;
        const int secondEta = second.i2 > 0 ? second.i2-1 : second.i2;
        delta2 = std::abs(firstEta-secondEta);
      } else delta2 = std::abs(first.i2-second.i2);
      return delta1 <= 1 && delta2 <= 1 && (delta1 != 0 || delta2 != 0);
    };
    enum class FirstPairType { Unknown, Local, NotLocal };
    auto classifyFirstPair = [&]( const std::vector<unsigned int>& ids, const std::vector<float>& energies ){
      if( ids.size() < 2 || energies.size() < 2 || ids[0] == 0 || ids[1] == 0 )
      return FirstPairType::Unknown;
      if( not std::isfinite(energies[0]) || not std::isfinite(energies[1]) )
      return FirstPairType::Unknown;
      const auto firstDet = fastDetInfo.find(ids[0]);
      const auto secondDet = fastDetInfo.find(ids[1]);
      if( firstDet == fastDetInfo.end() || secondDet == fastDetInfo.end() )
      return FirstPairType::Unknown;
      const bool energyOrdered = energies[0] > 0.5f && energies[1] > 0.5f &&
      energies[0] > energies[1] && energies[1] > 0.8f*energies[0];
      return energyOrdered && areNeighboringCrystals(*firstDet->second,*secondDet->second)
      ? FirstPairType::Local : FirstPairType::NotLocal;
    };
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

      if( fInTree->GetEntry(centry) <= 0 ){
        std::cerr << "ERROR: failed to read resolution input entry " << centry << std::endl;
        delete fInTree;
        return false;
      }//<<>>if( entry read failed )

      if(debug) std::cout << "Start loop ----------------------------------- " << std::endl;

      if( not resRhID || not resAmp || not resE || not resRtTime || not resTOF ||
      (needGainInfo && (not rhID || not rhisGS6 || not rhisGS1)) ||
      (useUnCCTime && not resCCTime) ){
        std::cerr << "ERROR: null vector in resolution input entry " << centry << std::endl;
        delete fInTree;
        return false;
      }//<<>>if( null input vector )
      const auto nResolutionHits = resRhID->size();
      const bool resolutionVectorsAligned = resAmp->size() == nResolutionHits &&
      resE->size() == nResolutionHits &&
      resRtTime->size() == nResolutionHits &&
      resTOF->size() == nResolutionHits &&
      (not useUnCCTime || resCCTime->size() == nResolutionHits);
      if( not resolutionVectorsAligned ){
        malformedResolutionEvents++;
        if( malformedResolutionEvents <= 5 ){
          std::cerr << "WARNING: skipping entry " << centry << " with unequal resolution-vector sizes:"
          << " ID=" << nResolutionHits << " Amp=" << resAmp->size()
          << " E=" << resE->size() << " Rt=" << resRtTime->size()
          << " TOF=" << resTOF->size();
          if( useUnCCTime ) std::cerr << " CC=" << resCCTime->size();
          std::cerr << std::endl;
        }//<<>>if( first malformed entries )
        continue;
      }//<<>>if( malformed resolution vectors )
      if( needGainInfo && (rhisGS6->size() != rhID->size() || rhisGS1->size() != rhID->size()) ){
        std::cerr << "ERROR: misaligned rechit gain vectors in entry " << centry << std::endl;
        delete fInTree;
        return false;
      }//<<>>if( malformed gain vectors )
      const auto nRecHits = needGainInfo ? rhID->size() : 0;

      std::array<int,4> resolutionSource = {-1,-1,-1,-1};
      if( nResolutionHits >= 4 ) resolutionSource = {0,1,2,3};
      else {
        // GammaResTool stores local0, local1, global0, global1 in that order, but a
        // positive ID missing from miniRhId is omitted and shortens every vector.
        const auto firstPairType = classifyFirstPair(*resRhID,*resE);
        if( nResolutionHits >= 2 && firstPairType == FirstPairType::Local ){
          resolutionSource[0] = 0;
          resolutionSource[1] = 1;
          recoveredLocalEvents++;
        } else if( nResolutionHits == 3 && firstPairType == FirstPairType::NotLocal ){
          // One local slot was omitted, so the last two surviving slots are global0/global1.
          resolutionSource[2] = 1;
          resolutionSource[3] = 2;
          recoveredGlobalEvents++;
        } else {
          skippedShortResolutionEvents++;
          continue;
        }//<<>>if( recoverable shortened resolution record )
      }//<<>>if( complete or recoverable resolution record )

      std::array<unsigned int,4> eventResRhID = {0,0,0,0};
      std::array<float,4> eventResAmp = {-999.f,-999.f,-999.f,-999.f};
      std::array<float,4> eventResE = {-999.f,-999.f,-999.f,-999.f};
      std::array<float,4> eventResRtTime = {-999.f,-999.f,-999.f,-999.f};
      std::array<float,4> eventResCCTime = {-999.f,-999.f,-999.f,-999.f};
      std::array<float,4> eventResTOF = {-999.f,-999.f,-999.f,-999.f};
      for( std::size_t logical = 0; logical < resolutionSource.size(); logical++ ){
        const int source = resolutionSource[logical];
        if( source < 0 ) continue;
        eventResRhID[logical] = (*resRhID)[source];
        eventResAmp[logical] = (*resAmp)[source];
        eventResE[logical] = (*resE)[source];
        eventResRtTime[logical] = (*resRtTime)[source];
        if( useUnCCTime ) eventResCCTime[logical] = (*resCCTime)[source];
        eventResTOF[logical] = (*resTOF)[source];
      }//<<>>for( logical resolution slots )

      if( eventResRhID[0] == 0 && eventResRhID[3] == 0 ) continue;

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
        if( fillResolutionHists ){

          // ----  fill res hist

          kucms_DetIDStruct invalidDetID;
          std::array<const kucms_DetIDStruct*,4> idinfo = {
            &invalidDetID,&invalidDetID,&invalidDetID,&invalidDetID
          };
          for( std::size_t idx = 0; idx < idinfo.size(); idx++ ){
            const auto detID = eventResRhID[idx];
            if( detID == 0 ) continue;
            const auto detIt = fastDetInfo.find(detID);
            if( detIt == fastDetInfo.end() ){
              std::cerr << "ERROR: unknown resolution rechit DetID " << detID
              << " in entry " << centry << std::endl;
              delete fInTree;
              return false;
            }//<<>>if( DetID missing )
            idinfo[idx] = detIt->second;
          }//<<>>for( resolution DetIDs )
          const auto& idinfoL0 = *idinfo[0];
          const auto& idinfoL1 = *idinfo[1];
          const auto& idinfoG0 = *idinfo[2];
          const auto& idinfoG1 = *idinfo[3];

          // set gain level in this block ---------------------------------------------------------------
          std::array<bool,4> isSelectedGain = {true,true,true,true};
          std::array<int,4> rhGainId = {0,0,0,0};
          if( needGainInfo ){
            std::array<bool,4> gainMatchFound = {false,false,false,false};
            std::size_t gainsNeeded = 0;
            for( std::size_t idx = 0; idx < gainMatchFound.size(); idx++ )
            if( eventResRhID[idx] != 0 ) gainsNeeded++;
            std::size_t gainsFound = 0;
            for( std::size_t rhit = 0; rhit < nRecHits && gainsFound < gainsNeeded; rhit++ ){
              const auto rechitID = (*rhID)[rhit];
              for( std::size_t resrhit = 0; resrhit < gainMatchFound.size(); resrhit++ ){
                if( gainMatchFound[resrhit] || eventResRhID[resrhit] == 0 ||
                eventResRhID[resrhit] != rechitID ) continue;
                gainMatchFound[resrhit] = true;
                gainsFound++;
                const bool isGS1 = (*rhisGS1)[rhit];
                const bool isGS6 = (*rhisGS6)[rhit];
                const int actualGainId = not isGS1 ? ( isGS6 ? 2 : 1 ) : ( isGS6 ? 4 : 3 );
                isSelectedGain[resrhit] = actualGainId == useGain;
                rhGainId[resrhit] = actualGainId;
              }//<<>>for( resolution rechits )
            }//<<>>for( full rechits )
            for( std::size_t idx = 0; idx < gainMatchFound.size(); idx++ ){
              if( eventResRhID[idx] != 0 && not gainMatchFound[idx] ){
                std::cerr << "ERROR: resolution rechit " << eventResRhID[idx]
                << " has no matching gain record in entry " << centry << std::endl;
                delete fInTree;
                return false;
              }//<<>>if( gain match missing )
            }//<<>>for( gain matches )
          }//<<>>if( useGSwitch )

          const double Time0 = useUnCCTime ? eventResCCTime[0] : eventResRtTime[0];
          const double Time1 = useUnCCTime ? eventResCCTime[1] : eventResRtTime[1];
          const double Time2 = useUnCCTime ? eventResCCTime[2] : eventResRtTime[2];
          const double Time3 = useUnCCTime ? eventResCCTime[3] : eventResRtTime[3];
          const bool validLocTime = std::isfinite(Time0) && std::isfinite(Time1) &&
          Time0 != -30.0 && Time1 != -30.0;
          const bool validGloTime = std::isfinite(Time2) && std::isfinite(Time3) &&
          Time2 != -30.0 && Time3 != -30.0;
          const bool validLocInputs = std::isfinite(eventResAmp[0]) && std::isfinite(eventResAmp[1]) &&
          std::isfinite(eventResE[0]) && std::isfinite(eventResE[1]) &&
          std::isfinite(eventResTOF[0]) && std::isfinite(eventResTOF[1]);
          const bool validGloInputs = std::isfinite(eventResAmp[2]) && std::isfinite(eventResAmp[3]) &&
          std::isfinite(eventResE[2]) && std::isfinite(eventResE[3]) &&
          std::isfinite(eventResTOF[2]) && std::isfinite(eventResTOF[3]);
          const bool goodLocRHs = eventResRhID[0] != 0 && eventResRhID[1] != 0;
          const bool goodGloRHs = eventResRhID[2] != 0 && eventResRhID[3] != 0;
          const bool goodLocTime = validLocTime && Time0 != 0 && Time1 != 0 && Time0 != Time1;
          const bool goodGloTime = validGloTime && Time2 != 0 && Time3 != 0 && Time2 != Time3;

          bool le_cut = (eventResE[0]>=lB)&&(eventResE[0]<=uB)&&(eventResE[1]>=lB)&&(eventResE[1]<=uB);
          bool ge_cut = (eventResE[2]>=lB)&&(eventResE[2]<=uB)&&(eventResE[3]>=lB)&&(eventResE[3]<=uB);
          if( useGSwitch ){
            if( doAllGain ){
              le_cut = eventResE[0] >= lB && eventResE[1] >= lB;
              ge_cut = eventResE[2] >= lB && eventResE[3] >= lB;
            } else {
              le_cut = isSelectedGain[0] && isSelectedGain[1] && eventResE[0] >= lB && eventResE[1] >= lB;
              ge_cut = isSelectedGain[2] && isSelectedGain[3] && eventResE[2] >= lB && eventResE[3] >= lB;
            }//<<>>if( all gains )
          }//<<>>if( gain selection enabled )
          const bool levent_good = le_cut && goodLocRHs && goodLocTime && validLocInputs;
          const bool gevent_good = ge_cut && goodGloRHs && goodGloTime && validGloInputs;

          const float seedTimeIC00 = usecali && levent_good ? getCachedCalibration(eventResRhID[0],rhGainId[0],idinfoL0) : 0;
          const float seedTimeIC10 = usecali && levent_good ? getCachedCalibration(eventResRhID[1],rhGainId[1],idinfoL1) : 0;
          const float seedTimeIC01 = usecali && gevent_good ? getCachedCalibration(eventResRhID[2],rhGainId[2],idinfoG0) : 0;
          const float seedTimeIC11 = usecali && gevent_good ? getCachedCalibration(eventResRhID[3],rhGainId[3],idinfoG1) : 0;
          double lyf0 = Time0-seedTimeIC00;
          double lyf1 = Time1-seedTimeIC10;
          double gyf0 = Time2-seedTimeIC01;
          double gyf1 = Time3-seedTimeIC11;

          if( smear ){
            //std::cout << "Times are smeared !!!!!" << std::endl;

            float amp0 = eventResRhID[0] > 0 ? eventResAmp[0] : -999;
            float amp1 = eventResRhID[1] > 0 ? eventResAmp[1] : -999;
            float amp2 = eventResRhID[2] > 0 ? eventResAmp[2] : -999;
            float amp3 = eventResRhID[3] > 0 ? eventResAmp[3] : -999;

            //getCorrectedTime( time, amplitude, rechitID, Evt_run, dataSetKey, mctype, gainID )
            if( validLocTime && validLocInputs && eventResRhID[0] != 0 && eventResRhID[1] != 0 ){
              lyf0 = getCorrectedTime( lyf0, amp0, eventResRhID[0], run, tag, 0, rhGainId[0], false );
              lyf1 = getCorrectedTime( lyf1, amp1, eventResRhID[1], run, tag, 0, rhGainId[1], false );
            }//<<>>if( valid local times )
            if( validGloTime && validGloInputs && eventResRhID[2] != 0 && eventResRhID[3] != 0 ){
              gyf0 = getCorrectedTime( gyf0, amp2, eventResRhID[2], run, tag, 0, rhGainId[2], false );
              gyf1 = getCorrectedTime( gyf1, amp3, eventResRhID[3], run, tag, 0, rhGainId[3], false );
            }//<<>>if( valid global times )


          }//<<>>if( smear )
          if( not levent_good && not gevent_good ) continue;
          const double leffa0 = useAmp ? eventResAmp[0] : eventResE[0];
          const double leffa1 = useAmp ? eventResAmp[1] : eventResE[1];
          const double geffa0 = useAmp ? eventResAmp[2] : eventResE[2];
          const double geffa1 = useAmp ? eventResAmp[3] : eventResE[3];
          const double leffe0 = eventResE[0];
          const double leffe1 = eventResE[1];
          const double lxfill = (leffa0*leffa1)/sqrt(leffa0*leffa0+leffa1*leffa1);
          const double lxfille = (leffe0*leffe1)/sqrt(leffe0*leffe0+leffe1*leffe1);
          const double gxfill = (geffa0*geffa1)/sqrt(geffa0*geffa0+geffa1*geffa1);
          const double lyfill = lyf0-lyf1+eventResTOF[0]-eventResTOF[1];
          const double gyfill = gyf0-gyf1+eventResTOF[2]-eventResTOF[3];
          const bool leta_cut_eb = idinfoL0.ecal == ECAL::EB && idinfoL1.ecal == ECAL::EB;
          const bool geta_cut_eb = idinfoG0.ecal == ECAL::EB && idinfoG1.ecal == ECAL::EB;
          const bool leta_cut_ee = idinfoL0.ecal != ECAL::EB && idinfoL1.ecal != ECAL::EB;
          const bool geta_cut_ee = idinfoG0.ecal != ECAL::EB && idinfoG1.ecal != ECAL::EB;
          const bool isd_cut = idinfoL0.TT == idinfoL1.TT; // true = same, false = different

          if(debug) std::cout << " - Fill 2D Hist" << std::endl;

          if( debug && ( levent_good || gevent_good ) ){ for( int didx = 0; didx < 4; didx++ ){
              //if( ( levent_good || gevent_good ) ){ for( int didx = 0; didx < 4; didx++ ){
              std::cout << "Run " << run << " id " << eventResRhID[didx];
              std::cout << " Amp " << eventResAmp[didx] << " E " << eventResE[didx];
              std::cout << " Rt " << eventResRtTime[didx]  << " TOF " << eventResTOF[didx] << std::endl;
              std::cout << "GSL: " << isSelectedGain[0] << " " << isSelectedGain[1] << " ";
              std::cout << eventResE[0] << " " << eventResE[1] << " " << le_cut << std::endl;
              std::cout << "GSZ: " << isSelectedGain[2] << " " << isSelectedGain[3] << " ";
              std::cout << eventResE[2] << " " << eventResE[3] << " " << ge_cut << std::endl;
            }}//<<>>for( int didx = 0; didx < 4; didx++ ){

          if( levent_good && debug ){
            std::cout << " -- IC0l: " << seedTimeIC00 << " IC1l: " << seedTimeIC10 << std::endl;
          }//<<>>if( levent_good && debug )
          if( gevent_good && debug ){
            std::cout << " -- IC0g: " << seedTimeIC01 << " IC1g: " << seedTimeIC11 << std::endl;
          }//<<>>if( gevent_good && debug )
          if( levent_good && isd_cut ){
            if( leta_cut_eb ){
              localSameEB->Fill(lxfill,lyfill,fillwgt);
              amplitudeEnergyEB->Fill(lxfill,lxfille,fillwgt); }
            if( leta_cut_ee ){
              localSameEE->Fill(lxfill,lyfill,fillwgt);
              amplitudeEnergyEE->Fill(lxfill,lxfille,fillwgt); }
          }//<<>>if( levent_good && isd_cut )
          if( gevent_good ){
            if( geta_cut_eb ) globalEB->Fill(gxfill,gyfill,fillwgt);
            if( geta_cut_ee ) globalEE->Fill(gxfill,gyfill,fillwgt);
          }//<<>>if( gevent_good )
          if( levent_good && not isd_cut ){
            if( leta_cut_eb ) localDifferentEB->Fill(lxfill,lyfill,fillwgt);
            if( leta_cut_ee ) localDifferentEE->Fill(lxfill,lyfill,fillwgt);
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
    std::cout << "[INFO] Resolution-record recovery: local=" << recoveredLocalEvents
    << " global=" << recoveredGlobalEvents
    << " skipped-short=" << skippedShortResolutionEvents
    << " skipped-unequal=" << malformedResolutionEvents << std::endl;
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
    return false;
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
  delete cali2DResTFile;
  cali2DResTFile = nullptr;

  std::cout << "Finished making 2D delta t v eff amp plots" << std::endl;
  return true;

}//<<>> bool plot2dResolutionEGR( std::string inputFileName, bool scale, bool usecali, bool smear, std::string ext )

kucms_SigmaFitResult KUCMS_TimeCalibration::runTimeFitter( TH2F* hist2D ){

  kucms_SigmaFitResult sigmafit;
  if( not hist2D ){
    std::cerr << "ERROR: null 2D histogram passed to runTimeFitter." << std::endl;
    return sigmafit;
  }//<<>>if( null input histogram )

  std::string f2DHistName = hist2D->GetName();
  std::cout << "Running time fitter : " << f2DHistName << std::endl;
  bool doSterm = true;
  //bool doSterm = false;

  // Use the input histogram's axis directly so the fitter cannot drift from the
  // binning used when the 2D resolution map was produced.
  const int nXbins = hist2D->GetNbinsX();
  if( nXbins <= 0 ){
    std::cerr << "ERROR: input resolution histogram has no X bins: " << f2DHistName << std::endl;
    return sigmafit;
  }//<<>>if( invalid X axis )
  std::vector<double> fXBins(nXbins+1);
  for( int ibinX = 1; ibinX <= nXbins; ibinX++ ){
    fXBins[ibinX-1] = hist2D->GetXaxis()->GetBinLowEdge(ibinX);
    if( not std::isfinite(fXBins[ibinX-1]) ||
    (ibinX > 1 && fXBins[ibinX-1] <= fXBins[ibinX-2]) ){
      std::cerr << "ERROR: invalid X-axis binning in " << f2DHistName << std::endl;
      return sigmafit;
    }//<<>>if( invalid X-bin edge )
  }//<<>>for( X bins )
  fXBins[nXbins] = hist2D->GetXaxis()->GetBinUpEdge(nXbins);
  if( not std::isfinite(fXBins[nXbins]) || fXBins[nXbins] <= fXBins[nXbins-1] ){
    std::cerr << "ERROR: invalid final X-axis edge in " << f2DHistName << std::endl;
    return sigmafit;
  }//<<>>if( invalid final edge )
  const auto xbins = fXBins.data();
  std::cout << " -- using " << nXbins << " X bins stored in the input histogram" << std::endl;

  std::string fittypename = ( doSterm ) ? "_NSC" : "_NC";
  //std::string outfilename = caliFileDir + f2DHistName + fittypename  + "_resfit.root";
  std::string outfilename = f2DHistName + fittypename  + "_resfit.root";
  TFile* resTFile = TFile::Open( outfilename.c_str(), "RECREATE" );
  if( not resTFile || resTFile->IsZombie() ){
    std::cerr << "ERROR: failed to create time-fit output file: " << outfilename << std::endl;
    if( resTFile ){ resTFile->Close(); delete resTFile; }
    return sigmafit;
  }//<<>>if( output file failed )
  resTFile->cd();

  std::map<std::string,TH1F*> ResultsMap;
  std::string fTitle = "#Delta(Photon Seed Time) [ns] vs. A_{eff}/#sigma_{n} (EBEB)";
  std::string fXTitle = "A_{eff}/#sigma_{n} (EBEB)";
  std::string fTimeText = "t_{1}-t_{2}";
  std::string fSigmaVarText = "A_{eff}/#sigma_{n}";
  std::string chi2ndfName = f2DHistName+"_chi2ndf";
  std::string chi2ndfTitle = fTitle+" "+"#chi^{2}/NDF."+";"+fXTitle+";"+"#chi^{2}/NDF.";
  ResultsMap["chi2ndf"]  = new TH1F(chi2ndfName.c_str(),chi2ndfTitle.c_str(),nXbins,xbins);
  ResultsMap["chi2ndf"]->Sumw2();
  std::string chi2Name = f2DHistName+"_chi2";
  ResultsMap["chi2"] = new TH1F(chi2Name.c_str(),chi2Name.c_str(),nXbins,xbins);
  std::string ndfName = f2DHistName+"_ndf";
  ResultsMap["ndf"] = new TH1F(ndfName.c_str(),ndfName.c_str(),nXbins,xbins);
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
  std::string sigmaRange1Name = f2DHistName+"_sigmaFitRange1";
  ResultsMap["sigmaFitRange1"] = new TH1F(sigmaRange1Name.c_str(),sigmaRange1Name.c_str(),nXbins,xbins);
  std::string sigmaRange4Name = f2DHistName+"_sigmaFitRange4";
  ResultsMap["sigmaFitRange4"] = new TH1F(sigmaRange4Name.c_str(),sigmaRange4Name.c_str(),nXbins,xbins);
  std::string sigmaModelErrName = f2DHistName+"_sigmaModelErr";
  ResultsMap["sigmaModelErr"] = new TH1F(sigmaModelErrName.c_str(),sigmaModelErrName.c_str(),nXbins,xbins);
  std::string sigmaModelVariationsName = f2DHistName+"_sigmaModelVariations";
  ResultsMap["sigmaModelVariations"] = new TH1F(sigmaModelVariationsName.c_str(),sigmaModelVariationsName.c_str(),nXbins,xbins);
  std::string sigmaNSCFitInputName = f2DHistName+"_sigmaNSCFitInput";
  ResultsMap["sigmaNSCFitInput"] = new TH1F(sigmaNSCFitInputName.c_str(),sigTitle.c_str(),nXbins,xbins);
  ResultsMap["sigmaNSCFitInput"]->Sumw2();
  std::string fitStatusName = f2DHistName+"_gaussFitStatus";
  ResultsMap["gaussFitStatus"] = new TH1F(fitStatusName.c_str(),fitStatusName.c_str(),nXbins,xbins);
  std::string covStatusName = f2DHistName+"_gaussCovStatus";
  ResultsMap["gaussCovStatus"] = new TH1F(covStatusName.c_str(),covStatusName.c_str(),nXbins,xbins);
  std::string fitValidName = f2DHistName+"_gaussFitValid";
  ResultsMap["gaussFitValid"] = new TH1F(fitValidName.c_str(),fitValidName.c_str(),nXbins,xbins);
  std::string fitStatusRange1Name = f2DHistName+"_gaussFitStatusRange1";
  ResultsMap["gaussFitStatusRange1"] = new TH1F(fitStatusRange1Name.c_str(),fitStatusRange1Name.c_str(),nXbins,xbins);
  std::string covStatusRange1Name = f2DHistName+"_gaussCovStatusRange1";
  ResultsMap["gaussCovStatusRange1"] = new TH1F(covStatusRange1Name.c_str(),covStatusRange1Name.c_str(),nXbins,xbins);
  std::string fitValidRange1Name = f2DHistName+"_gaussFitValidRange1";
  ResultsMap["gaussFitValidRange1"] = new TH1F(fitValidRange1Name.c_str(),fitValidRange1Name.c_str(),nXbins,xbins);
  std::string fitStatusRange4Name = f2DHistName+"_gaussFitStatusRange4";
  ResultsMap["gaussFitStatusRange4"] = new TH1F(fitStatusRange4Name.c_str(),fitStatusRange4Name.c_str(),nXbins,xbins);
  std::string covStatusRange4Name = f2DHistName+"_gaussCovStatusRange4";
  ResultsMap["gaussCovStatusRange4"] = new TH1F(covStatusRange4Name.c_str(),covStatusRange4Name.c_str(),nXbins,xbins);
  std::string fitValidRange4Name = f2DHistName+"_gaussFitValidRange4";
  ResultsMap["gaussFitValidRange4"] = new TH1F(fitValidRange4Name.c_str(),fitValidRange4Name.c_str(),nXbins,xbins);
  std::string sigmaAtLimitName = f2DHistName+"_gaussSigmaAtLimit";
  ResultsMap["gaussSigmaAtLimit"] = new TH1F(sigmaAtLimitName.c_str(),sigmaAtLimitName.c_str(),nXbins,xbins);
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
  for (int ibinX = 1; ibinX <= nXbins; ibinX++){
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
  for (int ibinX = 1; ibinX <= nXbins; ibinX++){
    profileHists[ibinX].DoFit();
  }//<<>>for (int ibinX = 1; ibinX < fNBinsX; ibinX++)
  //for( auto& profile : profileHists ){ profile.second.DoFit(); }
  std::cout << "Saving profile hists for: " << f2DHistName << std::endl;

  // Extract mu and sigma into maps
  //---------------TimeFitter::ExtractFitResults(FitInfo);
  std::cout << "Extracting results for: " << f2DHistName << std::endl;
  // set bin content!
  for (auto ibinX = 1; ibinX <= nXbins; ibinX++){

    auto& profile = profileHists[ibinX];
    //if( profile.isEmpty() ) continue;
    kucms_TimeFitResult results = profile.GetFitResult();

    //float herr = std::abs( results.sigmaHigh - results.sigma );
    //float lerr = std::abs( results.sigmaLow - results.sigma );
    //float serr = ( herr > lerr ) ? herr : lerr;
    float serr = results.esigma;

    std::cout << " Bin " << ibinX  <<  " : " << fXBins[ibinX-1] << "-" << fXBins[ibinX];
    std::cout << " : no guass: " << results.sigma << " err: " << serr << std::endl;
    ResultsMap["gaussFitStatus"]->SetBinContent(ibinX,results.fitStatus);
    ResultsMap["gaussCovStatus"]->SetBinContent(ibinX,results.covStatus);
    ResultsMap["gaussFitValid"]->SetBinContent(ibinX,results.valid ? 1 : 0);
    ResultsMap["gaussFitStatusRange1"]->SetBinContent(ibinX,results.fitStatusLow);
    ResultsMap["gaussCovStatusRange1"]->SetBinContent(ibinX,results.covStatusLow);
    ResultsMap["gaussFitValidRange1"]->SetBinContent(ibinX,results.validLow ? 1 : 0);
    ResultsMap["gaussFitStatusRange4"]->SetBinContent(ibinX,results.fitStatusHigh);
    ResultsMap["gaussCovStatusRange4"]->SetBinContent(ibinX,results.covStatusHigh);
    ResultsMap["gaussFitValidRange4"]->SetBinContent(ibinX,results.validHigh ? 1 : 0);
    ResultsMap["gaussSigmaAtLimit"]->SetBinContent(ibinX,results.sigmaAtLimit ? 1 : 0);
    ResultsMap["sigmaModelVariations"]->SetBinContent(ibinX,results.modelVariations);

    if( results.sigma <= 0 ) continue;

    ResultsMap["chi2"]->SetBinContent( ibinX, results.chi2 );
    ResultsMap["ndf"]->SetBinContent( ibinX, results.ndf );
    ResultsMap["chi2ndf"]->SetBinContent( ibinX, results.chi2ndf );
    ResultsMap["chi2prob"]->SetBinContent( ibinX, results.chi2prob );
    ResultsMap["mu"]->SetBinContent( ibinX, results.mu );
    ResultsMap["mu"]->SetBinError( ibinX, results.emu );
    if( results.validLow ) ResultsMap["sigmaFitRange1"]->SetBinContent(ibinX,results.sigmaLow);
    if( results.validHigh ) ResultsMap["sigmaFitRange4"]->SetBinContent(ibinX,results.sigmaHigh);
    if( results.modelValid ) ResultsMap["sigmaModelErr"]->SetBinContent(ibinX,results.modelError);
    if( results.esigma < results.sigma/2 ){
      ResultsMap["sigma"]->SetBinContent( ibinX, results.sigma );
      ResultsMap["sigma"]->SetBinError( ibinX, serr );
      if( results.valid && results.modelValid ){
        ResultsMap["sigmaNSCFitInput"]->SetBinContent(ibinX,results.sigma);
        ResultsMap["sigmaNSCFitInput"]->SetBinError(
        ibinX,std::hypot(static_cast<double>(serr),static_cast<double>(results.modelError)));
      }//<<>>if( valid model-aware NSC point )
    }//<<>>if( results.esigma < 0.01 )
    ResultsMap["occ"]->SetBinContent( ibinX, results.occ );
    ResultsMap["occ"]->SetBinError( ibinX, std::sqrt(results.occ) );
    ResultsMap["rms"]->SetBinContent( ibinX, results.rms );
    ResultsMap["rms"]->SetBinError( ibinX, 0 );

  }//<<>>for (auto ibinX = 1; ibinX <= fNBinsX; ibinX++)

  for( auto& profile : profileHists ){
    profile.second.profileHist->Write( profile.second.profileHist->GetName(), TObject::kOverwrite );
    if( doGPRHists ){
      TH1F* gprHist = gprProfileHist( profile.second.profileHist );
      if( gprHist ){
        gprHist->Write( gprHist->GetName(), TObject::kOverwrite );
        delete gprHist;
      }//<<>>if( GPR histogram created )
    }//<<>>if( deprecated GPR output enabled )
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
  auto& hist = ResultsMap["sigmaNSCFitInput"];
  //auto x_low = hist->GetXaxis()->GetBinLowEdge( hist->GetXaxis()->GetFirst() );
  //auto x_up  = hist->GetXaxis()->GetBinUpEdge( hist->GetXaxis()->GetLast() );
  float x_low = 75;
  float x_up = 1600;
  if( lowEnergy ){ x_low = 5; x_up = 1600; }
  if( doAllGain ){ x_up = 3200; }
  if( useEffEnergy ){ x_low = 5; x_up = 240; }
  if( lowEnergy && useEffEnergy ){ x_low = 2; }
  if( doAllGain && useEffEnergy ){ x_up = 2000; }
  for( int ibinX = 1; ibinX <= hist->GetNbinsX(); ibinX++ ){
    const double xvalue = hist->GetBinCenter(ibinX);
    const double content = hist->GetBinContent(ibinX);
    const double error = hist->GetBinError(ibinX);
    if( xvalue >= x_low && xvalue <= x_up && std::isfinite(content) && content > 0 &&
    std::isfinite(error) && error > 0 ) sigmafit.nFitPoints++;
  }//<<>>for( candidate NSC bins )
  const int nParameters = doSterm ? 3 : 2;
  if( sigmafit.nFitPoints <= nParameters ){
    std::cerr << "ERROR: insufficient valid points for NSC fit in " << f2DHistName
    << ": " << sigmafit.nFitPoints << std::endl;
    for( auto& rmap : ResultsMap ){
      rmap.second->Write(rmap.second->GetName(),TObject::kOverwrite);
      delete rmap.second;
    }//<<>>for( auto& rmap : ResultsMap )
    hist2D->Write(hist2D->GetName(),TObject::kOverwrite);
    resTFile->Close();
    delete resTFile;
    return sigmafit;
  }//<<>>if( insufficient NSC points )

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
  TFitResultPtr nscFitResult = hist->Fit(fit,"RBQONS");

  sigmafit.noise = fit->GetParameter(0);
  sigmafit.enoise = fit->GetParError (0);
  sigmafit.stant = fit->GetParameter(1);
  sigmafit.estant = fit->GetParError (1);
  sigmafit.stoch = doSterm ? fit->GetParameter(2) : -99;
  sigmafit.estoch = doSterm ? fit->GetParError (2) : -99;
  sigmafit.chi2 = fit->GetChisquare();
  sigmafit.ndf = fit->GetNDF();
  sigmafit.chi2ndf = sigmafit.ndf > 0 ? sigmafit.chi2/sigmafit.ndf : 0;
  sigmafit.chi2prob = fit->GetProb();
  sigmafit.fitStatus = static_cast<int>(nscFitResult);
  sigmafit.covStatus = nscFitResult.Get() ? nscFitResult->CovMatrixStatus() : -1;
  for( int parameter = 0; parameter < nParameters; parameter++ ){
    double lower = 0;
    double upper = 0;
    fit->GetParLimits(parameter,lower,upper);
    if( upper <= lower ) continue;
    const double value = fit->GetParameter(parameter);
    const double tolerance = 1e-6*std::max(1.0,std::abs(upper-lower));
    if( value <= lower+tolerance || value >= upper-tolerance ) sigmafit.parameterAtLimit = true;
  }//<<>>for( NSC parameters )
  const bool stochasticValid = not doSterm ||
  (std::isfinite(sigmafit.stoch) && std::isfinite(sigmafit.estoch) && sigmafit.estoch >= 0);
  sigmafit.valid = sigmafit.fitStatus == 0 && sigmafit.covStatus >= 2 && sigmafit.ndf > 0 &&
  std::isfinite(sigmafit.noise) && std::isfinite(sigmafit.enoise) && sigmafit.enoise >= 0 &&
  std::isfinite(sigmafit.stant) && std::isfinite(sigmafit.estant) && sigmafit.estant >= 0 &&
  stochasticValid &&
  std::isfinite(sigmafit.chi2) && std::isfinite(sigmafit.chi2prob);

  // Fit the same model-aware sigma points without the stochastic term so the
  // NC and NSC assumptions can be compared directly. The NSC result remains
  // the primary return value used by the calibration class.
  kucms_SigmaFitResult ncfit;
  ncfit.nFitPoints = sigmafit.nFitPoints;
  std::string ncFormName = histname+"_form_NC";
  std::string ncFitName = histname+"_fit_NC";
  auto ncForm = new TFormula(ncFormName.c_str(),"sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))");
  auto ncFit = new TF1(ncFitName.c_str(),ncForm->GetName(),x_low,x_up);
  ncFit->SetParName(0,"N");
  ncFit->SetParameter(0,nVal);
  ncFit->SetParLimits(0,nLower,nUpper);
  ncFit->SetParName(1,"C");
  ncFit->SetParameter(1,cVal);
  ncFit->SetParLimits(1,cLower,cUpper);
  ncFit->SetLineColor(hist->GetLineColor()+1);
  TFitResultPtr ncFitResult = hist->Fit(ncFit,"RBQONS");
  ncfit.noise = ncFit->GetParameter(0);
  ncfit.enoise = ncFit->GetParError(0);
  ncfit.stant = ncFit->GetParameter(1);
  ncfit.estant = ncFit->GetParError(1);
  ncfit.stoch = -99;
  ncfit.estoch = -99;
  ncfit.chi2 = ncFit->GetChisquare();
  ncfit.ndf = ncFit->GetNDF();
  ncfit.chi2ndf = ncfit.ndf > 0 ? ncfit.chi2/ncfit.ndf : 0;
  ncfit.chi2prob = ncFit->GetProb();
  ncfit.fitStatus = static_cast<int>(ncFitResult);
  ncfit.covStatus = ncFitResult.Get() ? ncFitResult->CovMatrixStatus() : -1;
  for( int parameter = 0; parameter < 2; parameter++ ){
    double lower = 0;
    double upper = 0;
    ncFit->GetParLimits(parameter,lower,upper);
    if( upper <= lower ) continue;
    const double value = ncFit->GetParameter(parameter);
    const double tolerance = 1e-6*std::max(1.0,std::abs(upper-lower));
    if( value <= lower+tolerance || value >= upper-tolerance ) ncfit.parameterAtLimit = true;
  }//<<>>for( NC parameters )
  ncfit.valid = ncfit.fitStatus == 0 && ncfit.covStatus >= 2 && ncfit.ndf > 0 &&
  std::isfinite(ncfit.noise) && std::isfinite(ncfit.enoise) && ncfit.enoise >= 0 &&
  std::isfinite(ncfit.stant) && std::isfinite(ncfit.estant) && ncfit.estant >= 0 &&
  std::isfinite(ncfit.chi2) && std::isfinite(ncfit.chi2prob);

  // Repeat both resolution-model fits using only the original statistical
  // uncertainties stored in the unchanged core sigma histogram. These are
  // diagnostics: the model-aware NSC result above remains the primary result.
  auto& statHist = ResultsMap["sigma"];
  kucms_SigmaFitResult statNscFit;
  for( int ibinX = 1; ibinX <= statHist->GetNbinsX(); ibinX++ ){
    const double xvalue = statHist->GetBinCenter(ibinX);
    const double content = statHist->GetBinContent(ibinX);
    const double error = statHist->GetBinError(ibinX);
    if( xvalue >= x_low && xvalue <= x_up && std::isfinite(content) && content > 0 &&
    std::isfinite(error) && error > 0 ) statNscFit.nFitPoints++;
  }//<<>>for( statistical-only candidate bins )
  std::string statNscFormName = f2DHistName+"_sigma_form_StatOnly_NSC";
  std::string statNscFitName = f2DHistName+"_sigma_fit_StatOnly_NSC";
  auto statNscForm = new TFormula(statNscFormName.c_str(),fitformstr.c_str());
  auto statNscFunction = new TF1(statNscFitName.c_str(),statNscForm->GetName(),x_low,x_up);
  statNscFunction->SetParName(0,"N");
  statNscFunction->SetParameter(0,nVal);
  statNscFunction->SetParLimits(0,nLower,nUpper);
  statNscFunction->SetParName(1,"C");
  statNscFunction->SetParameter(1,cVal);
  statNscFunction->SetParLimits(1,cLower,cUpper);
  if( doSterm ){
    statNscFunction->SetParName(2,"S");
    statNscFunction->SetParameter(2,1.0);
    statNscFunction->SetParLimits(2,0.0,10.0);
  }//<<>>if( stochastic term enabled )
  statNscFunction->SetLineColor(statHist->GetLineColor()+2);
  TFitResultPtr statNscFitResult = statHist->Fit(statNscFunction,"RBQONS");
  statNscFit.noise = statNscFunction->GetParameter(0);
  statNscFit.enoise = statNscFunction->GetParError(0);
  statNscFit.stant = statNscFunction->GetParameter(1);
  statNscFit.estant = statNscFunction->GetParError(1);
  statNscFit.stoch = doSterm ? statNscFunction->GetParameter(2) : -99;
  statNscFit.estoch = doSterm ? statNscFunction->GetParError(2) : -99;
  statNscFit.chi2 = statNscFunction->GetChisquare();
  statNscFit.ndf = statNscFunction->GetNDF();
  statNscFit.chi2ndf = statNscFit.ndf > 0 ? statNscFit.chi2/statNscFit.ndf : 0;
  statNscFit.chi2prob = statNscFunction->GetProb();
  statNscFit.fitStatus = static_cast<int>(statNscFitResult);
  statNscFit.covStatus = statNscFitResult.Get() ? statNscFitResult->CovMatrixStatus() : -1;
  for( int parameter = 0; parameter < nParameters; parameter++ ){
    double lower = 0;
    double upper = 0;
    statNscFunction->GetParLimits(parameter,lower,upper);
    if( upper <= lower ) continue;
    const double value = statNscFunction->GetParameter(parameter);
    const double tolerance = 1e-6*std::max(1.0,std::abs(upper-lower));
    if( value <= lower+tolerance || value >= upper-tolerance ) statNscFit.parameterAtLimit = true;
  }//<<>>for( statistical-only NSC parameters )
  const bool statStochasticValid = not doSterm ||
  (std::isfinite(statNscFit.stoch) && std::isfinite(statNscFit.estoch) && statNscFit.estoch >= 0);
  statNscFit.valid = statNscFit.fitStatus == 0 && statNscFit.covStatus >= 2 && statNscFit.ndf > 0 &&
  std::isfinite(statNscFit.noise) && std::isfinite(statNscFit.enoise) && statNscFit.enoise >= 0 &&
  std::isfinite(statNscFit.stant) && std::isfinite(statNscFit.estant) && statNscFit.estant >= 0 &&
  statStochasticValid && std::isfinite(statNscFit.chi2) && std::isfinite(statNscFit.chi2prob);

  kucms_SigmaFitResult statNcFit;
  statNcFit.nFitPoints = statNscFit.nFitPoints;
  std::string statNcFormName = f2DHistName+"_sigma_form_StatOnly_NC";
  std::string statNcFitName = f2DHistName+"_sigma_fit_StatOnly_NC";
  auto statNcForm = new TFormula(statNcFormName.c_str(),"sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))");
  auto statNcFunction = new TF1(statNcFitName.c_str(),statNcForm->GetName(),x_low,x_up);
  statNcFunction->SetParName(0,"N");
  statNcFunction->SetParameter(0,nVal);
  statNcFunction->SetParLimits(0,nLower,nUpper);
  statNcFunction->SetParName(1,"C");
  statNcFunction->SetParameter(1,cVal);
  statNcFunction->SetParLimits(1,cLower,cUpper);
  statNcFunction->SetLineColor(statHist->GetLineColor()+3);
  TFitResultPtr statNcFitResult = statHist->Fit(statNcFunction,"RBQONS");
  statNcFit.noise = statNcFunction->GetParameter(0);
  statNcFit.enoise = statNcFunction->GetParError(0);
  statNcFit.stant = statNcFunction->GetParameter(1);
  statNcFit.estant = statNcFunction->GetParError(1);
  statNcFit.stoch = -99;
  statNcFit.estoch = -99;
  statNcFit.chi2 = statNcFunction->GetChisquare();
  statNcFit.ndf = statNcFunction->GetNDF();
  statNcFit.chi2ndf = statNcFit.ndf > 0 ? statNcFit.chi2/statNcFit.ndf : 0;
  statNcFit.chi2prob = statNcFunction->GetProb();
  statNcFit.fitStatus = static_cast<int>(statNcFitResult);
  statNcFit.covStatus = statNcFitResult.Get() ? statNcFitResult->CovMatrixStatus() : -1;
  for( int parameter = 0; parameter < 2; parameter++ ){
    double lower = 0;
    double upper = 0;
    statNcFunction->GetParLimits(parameter,lower,upper);
    if( upper <= lower ) continue;
    const double value = statNcFunction->GetParameter(parameter);
    const double tolerance = 1e-6*std::max(1.0,std::abs(upper-lower));
    if( value <= lower+tolerance || value >= upper-tolerance ) statNcFit.parameterAtLimit = true;
  }//<<>>for( statistical-only NC parameters )
  statNcFit.valid = statNcFit.fitStatus == 0 && statNcFit.covStatus >= 2 && statNcFit.ndf > 0 &&
  std::isfinite(statNcFit.noise) && std::isfinite(statNcFit.enoise) && statNcFit.enoise >= 0 &&
  std::isfinite(statNcFit.stant) && std::isfinite(statNcFit.estant) && statNcFit.estant >= 0 &&
  std::isfinite(statNcFit.chi2) && std::isfinite(statNcFit.chi2prob);

  // Determine the additional fractional scatter required by the sigma points
  // instead of assuming the fixed 3% used by the final plotting code.  This is
  // a separate diagnostic fit; it does not alter ResultsMap["sigma"], nor does
  // it replace the model-aware NSC result returned by runTimeFitter().
  //
  // For each accepted sigma point y_i at energy x_i, the resolution model is
  //
  //   m_i^2 = N^2/x_i^2 + 2 C^2 + S^2/x_i,
  //
  // and the variance used by this likelihood is
  //
  //   V_i = e_stat,i^2 + (f y_i)^2 .
  //
  // e_stat is the original Gaussian-fit error and f is a fitted fractional
  // scatter.  The objective is the Gaussian -2 log likelihood (up to a common
  // additive constant), sum_i[(y_i-m_i)^2/V_i + log(V_i)].  The log(V_i)
  // term is essential: without it, increasing f would always reduce the
  // residual term and the fit would choose an artificially large uncertainty.
  //
  // The NSC fit has free parameters (N,C,S,f); the NC fit fixes S=0 and has
  // (N,C,f).  The fitted f is reported in percent and companion histograms
  // carry e_eff = hypot(e_stat, f*y).  Since f is estimated from the same
  // points, its chi2 probability is a diagnostic rather than an independent
  // goodness-of-fit test.
  std::vector<std::array<double,3>> scatterPoints;
  for( int ibinX = 1; ibinX <= statHist->GetNbinsX(); ibinX++ ){
    const double xvalue = statHist->GetBinCenter(ibinX);
    const double content = statHist->GetBinContent(ibinX);
    const double error = statHist->GetBinError(ibinX);
    if( xvalue >= x_low && xvalue <= x_up && std::isfinite(content) && content > 0 &&
    std::isfinite(error) && error > 0 ) scatterPoints.push_back({xvalue,content,error});
  }//<<>>for( fractional-scatter candidate bins )

  auto runFractionalScatterFit = [&]( bool includeStochastic, const std::string& functionName ){
    kucms_SigmaFitResult result;
    result.nFitPoints = static_cast<int>(scatterPoints.size());
    const unsigned int dimensions = includeStochastic ? 4 : 3;
    const unsigned int scatterIndex = dimensions-1;
    std::unique_ptr<TF1> resolutionFunction(new TF1(
    functionName.c_str(),
    includeStochastic ? "sqrt((([0]*[0])/(x*x))+(2*[1]*[1])+([2]*[2]/x))"
    : "sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))",
    x_low,x_up));
    if( result.nFitPoints <= static_cast<int>(dimensions) ){
      return std::make_pair(result,std::move(resolutionFunction));
    }//<<>>if( insufficient fractional-scatter points )

    std::function<double(const double*)> objective = [&]( const double* parameters ){
      const double noise = parameters[0];
      const double constant = parameters[1];
      const double stochastic = includeStochastic ? parameters[2] : 0.0;
      const double fractionalScatter = parameters[scatterIndex];
      double value = 0;
      for( const auto& point : scatterPoints ){
        const double xvalue = point[0];
        const double measuredSigma = point[1];
        const double statisticalError = point[2];
        const double modelSquared = (noise*noise)/(xvalue*xvalue) + 2*constant*constant +
        (stochastic*stochastic)/xvalue;
        const double variance = statisticalError*statisticalError +
        fractionalScatter*fractionalScatter*measuredSigma*measuredSigma;
        if( not std::isfinite(modelSquared) || modelSquared <= 0 ||
        not std::isfinite(variance) || variance <= 0 ){
          return std::numeric_limits<double>::max()/1024.0;
        }//<<>>if( invalid likelihood point )
        const double residual = measuredSigma-std::sqrt(modelSquared);
        value += residual*residual/variance + std::log(variance);
      }//<<>>for( sigma points )
      return value;
    };

    std::unique_ptr<ROOT::Math::Minimizer> minimizer(
    ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad"));
    if( not minimizer ) return std::make_pair(result,std::move(resolutionFunction));
    ROOT::Math::Functor objectiveFunctor(objective,dimensions);
    minimizer->SetFunction(objectiveFunctor);
    minimizer->SetErrorDef(1.0);
    minimizer->SetMaxFunctionCalls(100000);
    minimizer->SetMaxIterations(10000);
    minimizer->SetTolerance(1e-4);
    minimizer->SetPrintLevel(0);
    minimizer->SetLimitedVariable(0,"N",nVal,0.1,nLower,nUpper);
    minimizer->SetLimitedVariable(1,"C",0.1,0.001,cLower,cUpper);
    if( includeStochastic ) minimizer->SetLimitedVariable(2,"S",1.0,0.01,0.0,10.0);
    minimizer->SetLimitedVariable(scatterIndex,"fractionalScatter",0.03,0.001,0.0,1.0);
    const bool minimized = minimizer->Minimize();
    if( minimized ) minimizer->Hesse();
    const double* parameters = minimizer->X();
    const double* errors = minimizer->Errors();
    result.noise = parameters[0];
    result.enoise = errors[0];
    result.stant = parameters[1];
    result.estant = errors[1];
    result.stoch = includeStochastic ? parameters[2] : -99;
    result.estoch = includeStochastic ? errors[2] : -99;
    result.fractionalScatter = parameters[scatterIndex];
    result.efractionalScatter = errors[scatterIndex];
    result.objective = minimizer->MinValue();
    result.fitStatus = minimizer->Status();
    result.covStatus = minimizer->CovMatrixStatus();
    result.ndf = result.nFitPoints-static_cast<int>(dimensions);
    for( const auto& point : scatterPoints ){
      const double stochastic = includeStochastic ? result.stoch : 0.0;
      const double model = std::sqrt((result.noise*result.noise)/(point[0]*point[0]) +
      2*result.stant*result.stant + (stochastic*stochastic)/point[0]);
      const double variance = point[2]*point[2] + result.fractionalScatter*result.fractionalScatter*
      point[1]*point[1];
      const double residual = point[1]-model;
      result.chi2 += residual*residual/variance;
    }//<<>>for( fractional-scatter chi2 diagnostic )
    result.chi2ndf = result.ndf > 0 ? result.chi2/result.ndf : 0;
    result.chi2prob = result.ndf > 0 ? TMath::Prob(result.chi2,result.ndf) : 0;
    const std::vector<std::array<double,2>> parameterLimits = includeStochastic ?
    std::vector<std::array<double,2>>{{nLower,nUpper},{cLower,cUpper},{0.0,10.0},{0.0,1.0}} :
    std::vector<std::array<double,2>>{{nLower,nUpper},{cLower,cUpper},{0.0,1.0}};
    for( unsigned int parameter = 0; parameter < dimensions; parameter++ ){
      const double tolerance = 1e-6*std::max(1.0,std::abs(
      parameterLimits[parameter][1]-parameterLimits[parameter][0]));
      if( parameters[parameter] <= parameterLimits[parameter][0]+tolerance ||
      parameters[parameter] >= parameterLimits[parameter][1]-tolerance ){
        result.parameterAtLimit = true;
      }//<<>>if( parameter at limit )
    }//<<>>for( fractional-scatter parameters )
    const bool stochasticValid = not includeStochastic ||
    (std::isfinite(result.stoch) && std::isfinite(result.estoch) && result.estoch >= 0);
    result.valid = minimized && result.fitStatus == 0 && result.covStatus >= 2 && result.ndf > 0 &&
    std::isfinite(result.noise) && std::isfinite(result.enoise) && result.enoise >= 0 &&
    std::isfinite(result.stant) && std::isfinite(result.estant) && result.estant >= 0 &&
    stochasticValid && std::isfinite(result.fractionalScatter) &&
    std::isfinite(result.efractionalScatter) && result.efractionalScatter >= 0 &&
    std::isfinite(result.objective) && std::isfinite(result.chi2);

    resolutionFunction->SetParName(0,"N");
    resolutionFunction->SetParameter(0,result.noise);
    resolutionFunction->SetParError(0,result.enoise);
    resolutionFunction->SetParName(1,"C");
    resolutionFunction->SetParameter(1,result.stant);
    resolutionFunction->SetParError(1,result.estant);
    if( includeStochastic ){
      resolutionFunction->SetParName(2,"S");
      resolutionFunction->SetParameter(2,result.stoch);
      resolutionFunction->SetParError(2,result.estoch);
    }//<<>>if( stochastic term enabled )
    return std::make_pair(result,std::move(resolutionFunction));
  };

  auto scatterNscOutput = runFractionalScatterFit(
  true,f2DHistName+"_sigma_fit_FractionalScatter_NSC");
  auto scatterNcOutput = runFractionalScatterFit(
  false,f2DHistName+"_sigma_fit_FractionalScatter_NC");
  const kucms_SigmaFitResult& scatterNscFit = scatterNscOutput.first;
  const kucms_SigmaFitResult& scatterNcFit = scatterNcOutput.first;
  auto makeScatterInput = [&]( const std::string& name, double fractionalScatter ){
    std::unique_ptr<TH1F> scatterInput(dynamic_cast<TH1F*>(statHist->Clone(name.c_str())));
    if( not scatterInput ) return scatterInput;
    scatterInput->SetDirectory(nullptr);
    for( int ibinX = 1; ibinX <= scatterInput->GetNbinsX(); ibinX++ ){
      const double content = scatterInput->GetBinContent(ibinX);
      const double statisticalError = scatterInput->GetBinError(ibinX);
      if( std::isfinite(content) && content > 0 && std::isfinite(statisticalError) &&
      statisticalError > 0 && std::isfinite(fractionalScatter) && fractionalScatter >= 0 ){
        scatterInput->SetBinError(ibinX,std::hypot(
        statisticalError,fractionalScatter*content));
      }//<<>>if( valid fractional-scatter input bin )
    }//<<>>for( sigma bins )
    return scatterInput;
  };
  auto scatterNscInput = makeScatterInput(
  f2DHistName+"_sigmaFractionalScatterNSCFitInput",scatterNscFit.fractionalScatter);
  auto scatterNcInput = makeScatterInput(
  f2DHistName+"_sigmaFractionalScatterNCFitInput",scatterNcFit.fractionalScatter);

  std::cout << "Sigma fit results:" << std::endl;
  std::cout << " Model-aware NSC: noise " << sigmafit.noise << " +/- " << sigmafit.enoise
  << " stoch " << sigmafit.stoch << " +/- " << sigmafit.estoch
  << " stant " << sigmafit.stant << " +/- " << sigmafit.estant << std::endl;
  std::cout << "      chi2/ndf " << sigmafit.chi2 << "/" << sigmafit.ndf
  << " = " << sigmafit.chi2ndf << " prob " << sigmafit.chi2prob
  << " status " << sigmafit.fitStatus << " cov " << sigmafit.covStatus
  << " atLimit " << sigmafit.parameterAtLimit << " valid " << sigmafit.valid << std::endl;
  std::cout << " Model-aware NC : noise " << ncfit.noise << " +/- " << ncfit.enoise
  << " stant " << ncfit.stant << " +/- " << ncfit.estant << std::endl;
  std::cout << "      chi2/ndf " << ncfit.chi2 << "/" << ncfit.ndf
  << " = " << ncfit.chi2ndf << " prob " << ncfit.chi2prob
  << " status " << ncfit.fitStatus << " cov " << ncfit.covStatus
  << " atLimit " << ncfit.parameterAtLimit << " valid " << ncfit.valid << std::endl;
  std::cout << " Stat-only NSC  : noise " << statNscFit.noise << " +/- " << statNscFit.enoise
  << " stoch " << statNscFit.stoch << " +/- " << statNscFit.estoch
  << " stant " << statNscFit.stant << " +/- " << statNscFit.estant << std::endl;
  std::cout << "                  chi2/ndf " << statNscFit.chi2 << "/" << statNscFit.ndf
  << " = " << statNscFit.chi2ndf << " prob " << statNscFit.chi2prob
  << " status " << statNscFit.fitStatus << " cov " << statNscFit.covStatus
  << " atLimit " << statNscFit.parameterAtLimit << " valid " << statNscFit.valid << std::endl;
  std::cout << " Stat-only NC   : noise " << statNcFit.noise << " +/- " << statNcFit.enoise
  << " stant " << statNcFit.stant << " +/- " << statNcFit.estant << std::endl;
  std::cout << "                  chi2/ndf " << statNcFit.chi2 << "/" << statNcFit.ndf
  << " = " << statNcFit.chi2ndf << " prob " << statNcFit.chi2prob
  << " status " << statNcFit.fitStatus << " cov " << statNcFit.covStatus
  << " atLimit " << statNcFit.parameterAtLimit << " valid " << statNcFit.valid << std::endl;
  std::cout << " Fitted-scatter NSC: noise " << scatterNscFit.noise << " +/- " << scatterNscFit.enoise
  << " stoch " << scatterNscFit.stoch << " +/- " << scatterNscFit.estoch
  << " stant " << scatterNscFit.stant << " +/- " << scatterNscFit.estant << std::endl;
  std::cout << "                     fractional scatter " << scatterNscFit.fractionalScatter
  << " +/- " << scatterNscFit.efractionalScatter << " ("
  << 100*scatterNscFit.fractionalScatter << " +/- "
  << 100*scatterNscFit.efractionalScatter << "%)"
  << " -2logL (constant omitted) " << scatterNscFit.objective << std::endl;
  std::cout << "                     chi2/ndf " << scatterNscFit.chi2 << "/" << scatterNscFit.ndf
  << " = " << scatterNscFit.chi2ndf << " prob " << scatterNscFit.chi2prob
  << " status " << scatterNscFit.fitStatus << " cov " << scatterNscFit.covStatus
  << " atLimit " << scatterNscFit.parameterAtLimit << " valid " << scatterNscFit.valid << std::endl;
  std::cout << " Fitted-scatter NC : noise " << scatterNcFit.noise << " +/- " << scatterNcFit.enoise
  << " stant " << scatterNcFit.stant << " +/- " << scatterNcFit.estant << std::endl;
  std::cout << "                     fractional scatter " << scatterNcFit.fractionalScatter
  << " +/- " << scatterNcFit.efractionalScatter << " ("
  << 100*scatterNcFit.fractionalScatter << " +/- "
  << 100*scatterNcFit.efractionalScatter << "%)"
  << " -2logL (constant omitted) " << scatterNcFit.objective << std::endl;
  std::cout << "                     chi2/ndf " << scatterNcFit.chi2 << "/" << scatterNcFit.ndf
  << " = " << scatterNcFit.chi2ndf << " prob " << scatterNcFit.chi2prob
  << " status " << scatterNcFit.fitStatus << " cov " << scatterNcFit.covStatus
  << " atLimit " << scatterNcFit.parameterAtLimit << " valid " << scatterNcFit.valid << std::endl;

  std::string fitQualityName = f2DHistName+"_NSCFitQuality";
  TH1F fitQuality(fitQualityName.c_str(),fitQualityName.c_str(),9,0.5,9.5);
  fitQuality.SetDirectory(nullptr);
  const std::vector<std::string> qualityLabels = {
    "fitStatus","covStatus","chi2","ndf","chi2ndf","prob","parameterAtLimit","nFitPoints","valid"
  };
  const std::vector<double> qualityValues = {
    static_cast<double>(sigmafit.fitStatus),static_cast<double>(sigmafit.covStatus),
    sigmafit.chi2,static_cast<double>(sigmafit.ndf),sigmafit.chi2ndf,sigmafit.chi2prob,
    sigmafit.parameterAtLimit ? 1.0 : 0.0,static_cast<double>(sigmafit.nFitPoints),
    sigmafit.valid ? 1.0 : 0.0
  };
  for( std::size_t index = 0; index < qualityLabels.size(); index++ ){
    fitQuality.GetXaxis()->SetBinLabel(index+1,qualityLabels[index].c_str());
    fitQuality.SetBinContent(index+1,qualityValues[index]);
  }//<<>>for( NSC quality fields )
  fitQuality.Write(fitQuality.GetName(),TObject::kOverwrite);
  std::string ncFitQualityName = f2DHistName+"_NCFitQuality";
  TH1F ncFitQuality(ncFitQualityName.c_str(),ncFitQualityName.c_str(),9,0.5,9.5);
  ncFitQuality.SetDirectory(nullptr);
  const std::vector<double> ncQualityValues = {
    static_cast<double>(ncfit.fitStatus),static_cast<double>(ncfit.covStatus),
    ncfit.chi2,static_cast<double>(ncfit.ndf),ncfit.chi2ndf,ncfit.chi2prob,
    ncfit.parameterAtLimit ? 1.0 : 0.0,static_cast<double>(ncfit.nFitPoints),
    ncfit.valid ? 1.0 : 0.0
  };
  for( std::size_t index = 0; index < qualityLabels.size(); index++ ){
    ncFitQuality.GetXaxis()->SetBinLabel(index+1,qualityLabels[index].c_str());
    ncFitQuality.SetBinContent(index+1,ncQualityValues[index]);
  }//<<>>for( NC quality fields )
  ncFitQuality.Write(ncFitQuality.GetName(),TObject::kOverwrite);
  auto writeFitQuality = [&]( const std::string& name, const kucms_SigmaFitResult& result ){
    TH1F quality(name.c_str(),name.c_str(),9,0.5,9.5);
    quality.SetDirectory(nullptr);
    const std::vector<double> values = {
      static_cast<double>(result.fitStatus),static_cast<double>(result.covStatus),
      result.chi2,static_cast<double>(result.ndf),result.chi2ndf,result.chi2prob,
      result.parameterAtLimit ? 1.0 : 0.0,static_cast<double>(result.nFitPoints),
      result.valid ? 1.0 : 0.0
    };
    for( std::size_t index = 0; index < qualityLabels.size(); index++ ){
      quality.GetXaxis()->SetBinLabel(index+1,qualityLabels[index].c_str());
      quality.SetBinContent(index+1,values[index]);
    }//<<>>for( statistical-only quality fields )
    quality.Write(quality.GetName(),TObject::kOverwrite);
  };
  writeFitQuality(f2DHistName+"_StatOnlyNSCFitQuality",statNscFit);
  writeFitQuality(f2DHistName+"_StatOnlyNCFitQuality",statNcFit);
  auto writeScatterFitQuality = [&]( const std::string& name, const kucms_SigmaFitResult& result ){
    TH1F quality(name.c_str(),name.c_str(),12,0.5,12.5);
    quality.SetDirectory(nullptr);
    const std::vector<std::string> labels = {
      "fitStatus","covStatus","chi2","ndf","chi2ndf","prob","parameterAtLimit",
      "nFitPoints","valid","fractionalScatter","fractionalScatterError","neg2LogLikelihood"
    };
    const std::vector<double> values = {
      static_cast<double>(result.fitStatus),static_cast<double>(result.covStatus),result.chi2,
      static_cast<double>(result.ndf),result.chi2ndf,result.chi2prob,
      result.parameterAtLimit ? 1.0 : 0.0,static_cast<double>(result.nFitPoints),
      result.valid ? 1.0 : 0.0,result.fractionalScatter,result.efractionalScatter,result.objective
    };
    for( std::size_t index = 0; index < labels.size(); index++ ){
      quality.GetXaxis()->SetBinLabel(index+1,labels[index].c_str());
      quality.SetBinContent(index+1,values[index]);
    }//<<>>for( fractional-scatter quality fields )
    quality.Write(quality.GetName(),TObject::kOverwrite);
  };
  writeScatterFitQuality(f2DHistName+"_FractionalScatterNSCFitQuality",scatterNscFit);
  writeScatterFitQuality(f2DHistName+"_FractionalScatterNCFitQuality",scatterNcFit);
  fit->Write(fit->GetName(),TObject::kOverwrite);
  ncFit->Write(ncFit->GetName(),TObject::kOverwrite);
  statNscFunction->Write(statNscFunction->GetName(),TObject::kOverwrite);
  statNcFunction->Write(statNcFunction->GetName(),TObject::kOverwrite);
  if( scatterNscOutput.second ){
    scatterNscOutput.second->Write(scatterNscOutput.second->GetName(),TObject::kOverwrite);
  }//<<>>if( scatterNscOutput.second )
  if( scatterNcOutput.second ){
    scatterNcOutput.second->Write(scatterNcOutput.second->GetName(),TObject::kOverwrite);
  }//<<>>if( scatterNcOutput.second )
  if( scatterNscInput ) scatterNscInput->Write(scatterNscInput->GetName(),TObject::kOverwrite);
  if( scatterNcInput ) scatterNcInput->Write(scatterNcInput->GetName(),TObject::kOverwrite);
  for( auto& rmap : ResultsMap ){ rmap.second->Write( rmap.second->GetName(), TObject::kOverwrite ); delete rmap.second;}
  hist2D->Write( hist2D->GetName(), TObject::kOverwrite );
  delete form;
  delete fit;
  delete ncForm;
  delete ncFit;
  delete statNscForm;
  delete statNscFunction;
  delete statNcForm;
  delete statNcFunction;
  resTFile->Close();
  delete resTFile;

  return sigmafit;

}//<<>>kucms_SigmaFitResult KUCMS_TimeCalibration::runTimeFitter( TH2F* hist2D )

TH1F* KUCMS_TimeCalibration::gprProfileHist( TH1F* hist ){

  if( not hist ){
    std::cerr << "ERROR: null histogram passed to gprProfileHist." << std::endl;
    return nullptr;
  }//<<>>if( null input histogram )

  std::cout << "Make GPR ProfileHist: " << hist->GetName() << std::endl;
  const int nBins = hist->GetNbinsX();
  const std::string newtitle = std::string(hist->GetName()) + "_gpr";
  if( nBins <= 0 ) return nullptr;

  const double threshold = 100.0;
  const int firstBin = hist->FindFirstBinAbove(threshold);
  const int lastBin = hist->FindLastBinAbove(threshold);
  const double sigma = std::abs(hist->GetStdDev());
  if( firstBin < 1 || lastBin < firstBin || not std::isfinite(sigma) || sigma <= 0 ){
    // Keep a valid diagnostic for empty or sparse slices rather than creating
    // a histogram with reversed or zero-width limits.
    auto fallback = dynamic_cast<TH1F*>(hist->Clone(newtitle.c_str()));
    if( fallback ) fallback->SetDirectory(nullptr);
    return fallback;
  }//<<>>if( adaptive rebinning is not defined )

  const int startBin = std::max(1,firstBin-1);
  const int endBin = std::min(nBins,lastBin+1);
  const double lowEdge = hist->GetXaxis()->GetBinLowEdge(startBin);
  const double highEdge = hist->GetXaxis()->GetBinUpEdge(endBin);
  const double extent = std::max(std::abs(lowEdge),std::abs(highEdge));
  if( not std::isfinite(extent) || extent <= 0 ) return nullptr;

  const double minBinWidth = 0.04;
  const double binWidth = std::max(minBinWidth,std::ceil((sigma/2.0)/minBinWidth)*minBinWidth);
  int nDiv = std::max(1,static_cast<int>(std::ceil((2.0*extent)/binWidth)));
  if( nDiv%2 == 0 ) nDiv++;
  const double symmetricExtent = nDiv*binWidth/2.0;
  auto newhist = new TH1F(newtitle.c_str(),newtitle.c_str(),nDiv,-symmetricExtent,symmetricExtent);
  newhist->SetDirectory(nullptr);
  newhist->Sumw2();
  std::vector<double> error2(nDiv+2,0.0);
  std::cout << " -- start : " << -symmetricExtent << " end : " << symmetricExtent
  << " nDiv : " << nDiv << " bin size : " << binWidth << std::endl;
  for( int ibinX = 1; ibinX <= nBins; ibinX++ ){
    const double content = hist->GetBinContent(ibinX);
    const int targetBin = newhist->FindBin(hist->GetBinCenter(ibinX));
    newhist->AddBinContent(targetBin,content);
    const double error = hist->GetBinError(ibinX);
    error2[targetBin] += error*error;
  }//<<>>for( input bins )
  for( int ibinX = 1; ibinX <= nDiv; ibinX++ )
  newhist->SetBinError(ibinX,std::max(1.0,std::sqrt(error2[ibinX])));

  return newhist;

}//<<>>TH1F* KUCMS_TimeCalibration::gprProfileHist( TH1F* hist )

bool KUCMS_TimeCalibration::doResTimeFits( bool doLocal ){

  // Fit the configured resolution histograms and store the resulting timing
  // parameters for subsequent correction and smearing calls.

  std::cout << "Making Resolution Time Fits Auto : ";
  if( doLocal ) std::cout << " Local " << std::endl;
  else std::cout << " Global " << std::endl;
  bool success = true;
  for( auto& runmaptag : CaliRunMapSet ){
    for( auto& runmap : runmaptag.second ){
      bool doResFit = not runmap.second.hasResParams;
      if( runmap.second.has2DResMap && doResFit ){

        std::string histtype = doLocal ? lochist : globhist;
        std::string fname = runmap.second.histMapName+histtype;
        std::cout << " Opening : " << fname << std::endl;
        auto histEntry = CaliHists.find(fname);
        if( histEntry == CaliHists.end() || not histEntry->second.h2f ){
          std::cerr << "ERROR: no loaded TH2F resolution histogram named " << fname << std::endl;
          success = false;
          continue;
        }//<<>>if( input histogram is unavailable )
        kucms_SigmaFitResult results = runTimeFitter(histEntry->second.h2f);
        //std::string gbfname = runmap.second.histMapName+globhist;
        //SigmaFitResult results = runTimeFitter( CaliHists[gbfname].h2f );
        if( not doLocal && results.valid ){
          runmap.second.hasResParams = true;
          runmap.second.noise = results.noise;
          runmap.second.stoch = results.stoch;
          runmap.second.stant = results.stant;
        }//<<>>if( successful global fit )
        if( not results.valid ) success = false;

      }//<<>>if( runmap.second.has2DResMap && not runmap.second.hasResParams )
    }//<<>>for( auto& runmaptag :  CaliRunMapSet )
  }//<<>>for( auto& runmaptag :  CaliRunMapSet )
  if( cali2DResTFile != NULL ) cali2DResTFile->Close();
  return success;

}//<<>>bool KUCMS_TimeCalibration::doResTimeFits( bool doLocal )

bool KUCMS_TimeCalibration::doResTimeFitExt( std::string histName, std::string fileName ){

  std::cout << "Making Resolution Time Fits" << std::endl;
  std::cout << " -- With : " << fileName << std::endl;
  auto tFile = TFile::Open( fileName.c_str(), "READ" );
  if( not tFile || tFile->IsZombie() ){
    std::cerr << "ERROR: failed to open: " << fileName << std::endl;
    if( tFile ){ tFile->Close(); delete tFile; }
    return false;
  }//<<>>if( input file failed )
  std::cout << " -- For : " << histName << std::endl;
  auto lhist = dynamic_cast<TH2F*>(tFile->Get(histName.c_str()));
  if( lhist == nullptr ){
    std::cout << " -- Failed to Open : " << histName << std::endl;
    tFile->Close();
    delete tFile;
    return false;
  }//<<>>if( histogram missing or wrong type )
  const kucms_SigmaFitResult results = runTimeFitter(lhist);
  tFile->Close();
  delete tFile;
  return results.valid;

}//<<>>bool KUCMS_TimeCalibration::doResTimeFitExt( std::string histName, std::string fileName )

bool KUCMS_TimeCalibration::load2DResHist( std::string histName ){

  std::cout << " --- Load hist from caliTFile: " << histName << std::endl;
  auto inputFile = TFile::Open(cali2DResPlotsTFileName.c_str(),"READ");
  if( not inputFile || inputFile->IsZombie() ){
    std::cerr << "ERROR: failed to open 2D resolution input file: "
    << cali2DResPlotsTFileName << std::endl;
    if( inputFile ){ inputFile->Close(); delete inputFile; }
    return false;
  }//<<>>if( input file failed )
  auto sourceHist = dynamic_cast<TH2F*>(inputFile->Get(histName.c_str()));
  if( not sourceHist ){
    std::cerr << "ERROR: TH2F resolution histogram not found: " << histName << std::endl;
    inputFile->Close();
    delete inputFile;
    return false;
  }//<<>>if( histogram missing or wrong type )
  auto loadedHist = dynamic_cast<TH2F*>(sourceHist->Clone(histName.c_str()));
  if( loadedHist ) loadedHist->SetDirectory(nullptr);
  inputFile->Close();
  delete inputFile;
  if( not loadedHist ){
    std::cerr << "ERROR: failed to clone 2D resolution histogram: " << histName << std::endl;
    return false;
  }//<<>>if( histogram clone failed )
  auto existingHist = CaliHists.find(histName);
  if( existingHist != CaliHists.end() ){
    if( not existingHist->second.isResHist ){
      std::cerr << "ERROR: refusing to replace non-resolution histogram named " << histName << std::endl;
      delete loadedHist;
      return false;
    }//<<>>if( name collides with a calibration histogram )
    delete existingHist->second.h2f;
    CaliHists.erase(existingHist);
  }//<<>>if( replacing a previously loaded resolution histogram )
  CaliHists[histName] = { loadedHist, histName, false, true };
  return true;

}//<<>>bool KUCMS_TimeCalibration::load2DResHist( std::string histName )

bool KUCMS_TimeCalibration::doResTimeFit( std::string histName ){

  std::cout << "Making Resolution Time Fits  : " << histName << std::endl;
  auto histEntry = CaliHists.find(histName);
  if( histEntry == CaliHists.end() ){
    std::cerr << "ERROR: no loaded resolution histogram named " << histName << std::endl;
    return false;
  }//<<>>if( histogram entry missing )
  if( not histEntry->second.h2f ){
    std::cerr << "ERROR: null TH2F for loaded resolution histogram " << histName << std::endl;
    return false;
  }//<<>>if( histogram pointer missing )
  const kucms_SigmaFitResult results = runTimeFitter(histEntry->second.h2f);
  if( not results.valid ){
    std::cerr << "ERROR: resolution time fit failed for " << histName << std::endl;
    return false;
  }//<<>>if( fit failed )
  return true;

}//<<>>bool KUCMS_TimeCalibration::doResTimeFit( std::string histName )

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

}//<<>>void KUCMS_TimeCalibration::makeSmearTag( std::string sourceName, std::string destName, std::string smearTag )


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
      }//<<>>if( previous == 0 )
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

}//<<>>void KUCMS_TimeCalibration::makeTTDiffMaps( int srun, int erun )

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
    //float maxtime = 0;
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
          //if( std::abs(mean) > maxtime ) maxtime == std::abs(mean);
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

}//<<>>void KUCMS_TimeCalibration::makeTTDriftMaps( std::string tag, int srun, int erun )

#endif
//-------------------------------------------------------------------------------------------------------------------
