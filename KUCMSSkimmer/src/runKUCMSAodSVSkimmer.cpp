//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSVSkimmer.hh"
#include <stdexcept>
//#include "KUCMSAodSkimmer_sv_cc.hh"
//#include "KUCMSAodSkimmer_cc_probe.hh"
//#include "KUCMSAodSkimmer_rh_cc.hh"
//#include "KUCMSAodSkimmer_v16_rh_cc.hh"
//
//example calls
//Using a 'master' list of lists:
// ./runKUCMSAodSVSkimmer.obj -i ntuple_master_lists/KUCMS_Ntuple_Master_DataPD_Files_List.txt --evtFirst 0 --evtLast 100 -o test_v41_fullDataList
//
//Using a single list of sample files
// ./runKUCMSAodSVSkimmer.obj -i ntuple_master_lists/MET_R18_SVIPM100_v31_MET_AOD_Run2018A.txt --evtFirst 0 --evtLast 100 -o test_v41
//Notes:
//	- hasGenInfo and mctype is set here based on file name (ie 'AOD' vs 'AODSIM')
//	- the time calibration tag (ttag here) is updated for MC automatically ('_mc' appended)
//
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

  bool hprint = false;
  string in_file;
  int evti = 0;
  int evtj = 0;
  bool useEvtGenWgt = true;
  bool noSVorPho = false; // ( must have a sig pho or a sv ) + base selection - loads photons in rjr
  bool genSigPerfect = false;
  bool hasGenInfo = false;
  string outfilename = "rjrskim_v39"; // included BHC processing and branchs, added triggerpath varibles
  string key = "";
  float xsec = 1;
  float glumass = 0;
  float n2mass = 0;
  int mctype= 0;
  string ttag = "r2_ul18";
  float mcw = 0;
  bool noBHC = false;
  bool noSV = false;
  bool hltpaths = true;
  string branchMaskFile = "";
  string dxySigScaleMode = "off";
  float dxySigScaleDeltaM = 0.f;
  try {
    for( int i = 1; i < argc; i++ ){

      const std::string option(argv[i]);
      auto nextValue = [&]() -> const char* {
        if( i + 1 >= argc ){
          throw std::invalid_argument("missing value for " + option);
        }
        return argv[++i];
      };
      auto parseInt = [&]( const char* value ){
        const std::string text(value);
        std::size_t parsed = 0;
        const int result = std::stoi(text, &parsed);
        if( parsed != text.size() ) throw std::invalid_argument("invalid integer for " + option + ": " + text);
        return result;
      };
      auto parseFloat = [&]( const char* value ){
        const std::string text(value);
        std::size_t parsed = 0;
        const float result = std::stof(text, &parsed);
        if( parsed != text.size() ) throw std::invalid_argument("invalid number for " + option + ": " + text);
        return result;
      };

      if( option == "--help" || option == "-h" ) hprint = true;
      else if( option == "-i" || option == "--input" ) in_file = nextValue();
      else if( option == "-o" || option == "--output" ) outfilename = nextValue();
      else if( option == "--evtFirst" ) evti = parseInt(nextValue());
      else if( option == "--evtLast" ) evtj = parseInt(nextValue());
      else if( option == "--hasGenInfo" ) hasGenInfo = true;
      else if( option == "--genSigPerfect" ) genSigPerfect = true;
      else if( option == "--noSVorPho" ) noSVorPho = true;
      else if( option == "--noBHC" ) noBHC = true;
      else if( option == "--noSV" ) noSV = true;
      else if( option == "--dataSetKey" ) key = nextValue();
      else if( option == "--xsec" ) xsec = parseFloat(nextValue());
      else if( option == "--gluinoMass" ) glumass = parseFloat(nextValue());
      else if( option == "--N2Mass" ) n2mass = parseFloat(nextValue());
      else if( option == "--timeCaliTag" ) ttag = nextValue();
      else if( option == "--MCweight" ) mcw = parseFloat(nextValue());
      else if( option == "--MCtype" ) mctype = parseInt(nextValue());
      else if( option == "--HLTPathsOff" ) hltpaths = false;
      else if( option == "--branchMask" || option == "--branch-mask" ) branchMaskFile = nextValue();
      else if( option == "--dxySigScale" ) dxySigScaleMode = nextValue();
      else if( option == "--svDxyDeltaM" ) dxySigScaleDeltaM = parseFloat(nextValue());
      else throw std::invalid_argument("unknown option: " + option);
    }
  } catch( const std::exception& error ){
    std::cerr << "ERROR: " << error.what() << std::endl;
    std::cerr << "Run '" << argv[0] << " --help' for usage." << std::endl;
    return 1;
  }

  if(hprint){
    cout << "Usage: " << argv[0] << " [options]" << endl;
    cout << "  options:" << endl;
    cout << "   --help(-h)                           print options" << endl;
    cout << "   --input(-i) [file]                   input file" << endl;
    cout << "   --output(-o) [file]                  output file tag" << endl;
    cout << "   --evtFirst [i] --evtLast [j]         skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
    cout << "   --hasGenInfo                         sample has gen info (default = false)" << endl;
    cout << "   --genSigPerfect                      turn on gensigperfect (default = false)" << endl;
    cout << "   --noSVorPho                          selection doesn't have SV or photon (default = false)" << endl;
    cout << "   --noBHC                              will not run BHC objects (default = false - BHC on)" << endl;
    cout << "   --noSV                               will not run SV collection (default = false - SVs on)" << endl;
    cout << "   --HLTPathsOff                        sets whether to do certain HLT path branches (default = true - on)" << endl;
    cout << "   --branchMask [file]                  use output branch mask file" << endl;
    cout << "   --dataSetKey                         set dataset key" << endl;
    cout << "   --xsec                               set cross-section (default = 1)" << endl;
    cout << "   --gluinoMass                         set gluino mass (default = 0)" << endl;
    cout << "   --N2Mass                             set N2 mass (default = 0)" << endl;
    cout << "   --timeCaliTag                        set time calibration tag (default = r2_ul18(_mc))" << endl;
    cout << "   --MCweight                           set MC weight (default = 0)" << endl;
    cout << "   --MCtype [type]			 set MC type (default = 0)" << endl;
    cout << "   --dxySigScale [off|nominal|up|down]  FastSim SV dxySig scale correction (default = off)" << endl;
    cout << "   --svDxyDeltaM [GeV]                  mN2-mN1 for --dxySigScale, required unless off" << endl;
    return 0;
  }

  if( in_file.empty() ){
    std::cerr << "ERROR: no input file was specified. Use --input or -i." << std::endl;
    return 1;
  }

  cout << "outfile tag " << outfilename << endl;
  outfilename = outfilename+".root";
  //const std::string listdir = "ntuple_master_lists/";
  //const std::string infilename = "KUCMS_Ntuple_Master_DataPD_Files_List.txt"; hasGenInfo = false;

  int skipCnt = 0; // used to skip files ( in tchian ) for fast processing - if( nFiles%skipCnt != 0 ) continue; --  disabled in code  --  
  std::string eosdir = "root://cmseos.fnal.gov//store/group/lpcsusylep/jaking/";
    KUCMSAodSkimmer llpgana;
    llpgana.SetEventRange(evti, evtj);
    llpgana.SetGenInfoFlag( hasGenInfo );
    llpgana.SetNoSVorPhoFlag( noSVorPho );
    llpgana.SetUseEvtGenWgtFlag( useEvtGenWgt );
    llpgana.SetGenSigPerfectFlag( genSigPerfect );
    llpgana.SetDoBHC( !noBHC ); //if noBHC == true -> do not run BHC
    llpgana.SetDoSV( !noSV ); //if noSV == true -> do not run SVs
    //from master list
    llpgana.SetDataSetKey(key);
    llpgana.SetCrossSection(xsec);
    llpgana.SetGluinoMass(glumass);
    llpgana.SetN2Mass(n2mass);
    llpgana.SetMCType(mctype);
    llpgana.SetTimeCalibrationTag(ttag);
    llpgana.SetMCWeight(mcw);
    llpgana.SetBranchMaskFile(branchMaskFile);
    llpgana.SetDxySigScale(dxySigScaleMode, dxySigScaleDeltaM);
    bool sigbase = false;
    if(in_file.find("_v33_") != string::npos && in_file.find("SMS") != string::npos)
      sigbase = true;
    llpgana.setNewSigBase(sigbase);
    llpgana.setHTLPathsBase(hltpaths);
    //llpgana.setNewSigBase(true);
    //this method takes in 1 list at a time
    eosdir = eosdir+"KUCMSNtuple/";
    if(in_file.find("://") == string::npos)
      in_file = eosdir+in_file;
    llpgana.kucmsAodSkimmer( in_file, outfilename);
    return 0;


}//<<>>int main ( int argc, char *argv[] )
