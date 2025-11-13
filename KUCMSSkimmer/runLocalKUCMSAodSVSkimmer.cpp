//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSkimmer_Main.hh"

// ------------------------------------------- main function ------------------------------------------------------------

int main ( int argc, char *argv[] ){

    const std::string listdir = "ntuple_master_lists/";
    //const string KUCMSAodSkimmer::eosdir = "root://cmseos.fnal.gov//store/user/jaking/";
    const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";
    //const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/anazario/";
    //const std::string eosdir = "root://cmseos.fnal.gov//store/user/janguian/";
	//const std::string eosdir = "/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/";				

	// set below //  bool noSVorPho = true; // only base selection - will not load phos in rjr
	// set below //  bool noSVorPho = false; // ( flase = must have a sig pho or a sv ) + base selection - loads photons in rjr ( defaults false )
    // set below //  bool genSigPerfect = true; // ( defaults false )
    bool useEvtGenWgt = true; // false for SMS MC only
    bool hasGenInfo = true; // typically false for data only 

    //no//const std::string infilename = "KUCMS_Ntuple_Master_BG_Files_List.txt";
    //no//const std::string infilename = "KUCMS_Ntuple_Master_DEG_Files_List.txt"; hasGenInfo = false;
    //no//const std::string infilename = "KUCMS_Ntuple_Master_GoGoG_Files_List.txt";
	//no//const std::string infilename = "KUCMS_Ntuple_Master_SMS_Sig_Files_List_untracked.txt"; useEvtGenWgt = false;
    //no//const std::string infilename = "KUCMS_Ntuple_Master_GoGoZ_Files_List.txt";
    //no//const std::string infilename = "KUCMS_Ntuple_Master_SqSqG_Files_List.txt";

	//const std::string infilename = "KUCMS_Ntuple_Master_NtupleTest_Files_List.txt"; hasGenInfo = false;

    const std::string infilename = "KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt"; useEvtGenWgt = false;
    //const std::string infilename = "KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt";
    //const std::string infilename = "KUCMS_Ntuple_Master_DataPD_Files_List.txt"; hasGenInfo = false;

    //const std::string outfilename = "_LLPGskim_v21_rjrvars.root"; // skim v21 use new process_name tags 
    //const std::string outfilename = "_LLPGskim_v22_rjrvars.root"; // skim v22 keeps 0 pho events 
    //const std::string outfilename = "_LLPGskim_v23_rjrvars.root"; // fixed 2 sig pho events 
    //const std::string outfilename = "_rjrskim_v26.root"; // jet veto on first 2 Sig Phos only
    //const std::string outfilename = "_rjrskim_v27.root"; // SV matching & SVs & lepton veto
    //const std::string outfilename = "_rjrskim_v28.root"; // MET cleaning flags 
    //const std::string outfilename = "_rjrskim_v29.root"; // Track collection consoldation, Mu - Ele collections, RJR subjets 1a, 2a, ect..
	//const std::string outfilename = "_rjrskim_v31.root"; // alternate RJR setup with invis in comb split rule
    //const std::string outfilename = "_rjrskim_v32.root"; // SV vars + BaysC && TimeCali added
    //const std::string outfilename = "_rjrskim_v36.root"; // expanded hts & hs rjr vars selection
    //const std::string outfilename = "_rjrskim_v37.root"; // added SV info, nIsoTracks, and rjr JXX 4vector information
    //const std::string outfilename = "_rjrskim_v38.root"; // added additional SV information and photon time significance
    //const std::string outfilename = "_rjrskim_v39.root"; // included BHC processing and branchs, added triggerpath varibles
    //const std::string outfilename = "_rjrskim_v40.root"; // added HEM veto, restructured skimmer for condor processing
    //const std::string outfilename = "_rjrskim_v41.root"; // modifided RJR oot photon v met handeling : 0 is now sig photon oot met correction only
    const std::string outfilename = "_rjrskim_v43.root"; // added EE photons to RJR, set looser id reqs for jets to include photons

    //int skipCnt = 0; // used to skip files ( in tchian ) for fast processing - if( nFiles%skipCnt != 0 ) continue; --  disabled in code  --  
	// !!!!!!!!!!!!!!!!!!!  alternate RJR setup with invis in comb split rule
    KUCMSAodSkimmer llpgana;
    //llpgana.SetNoSVorPhoFlag(noSVorPho);
    //llpgana.SetGenSigPerfectFlag(genSigPerfect);
    llpgana.SetDoBHC(false);
	llpgana.SetGenInfoFlag(hasGenInfo); 
	llpgana.SetUseEvtGenWgtFlag(useEvtGenWgt);
    llpgana.kucmsAodSkimmer_local( listdir, eosdir, infilename, outfilename );
    return 1;


}//<<>>int main ( int argc, char *argv[] )
