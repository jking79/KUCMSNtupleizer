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

    std::string infilename = std::string(argv[0]); 
    bool hasGenInfo = ( std::stoi(argv[1]) == 1 ) ? true : false; 
	bool useEvtGenWgt = ( std::stoi(argv[1]) == 1 ) ? true : false;

    const std::string outfilename = "_rjrskim_vTest.root";
    const std::string listdir = "ntuple_master_lists/"; 
    const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/"; 
    bool noSVorPho = false; // ( must have a sig pho or a sv ) + base selection - loads photons in rjr 
    bool genSigPerfect = false; 

    int skipCnt = 0; // used to skip files ( in tchian ) for fast processing - if( nFiles%skipCnt != 0 ) continue; --  disabled in code  --  
	// !!!!!!!!!!!!!!!!!!!  alternate RJR setup with invis in comb split rule
    KUCMSAodSkimmer llpgana;
	llpgana.SetGenInfoFlag(hasGenInfo); 
	llpgana.SetNoSVorPhoFlag(noSVorPho); 
	llpgana.SetUseEvtGenWgtFlag(useEvtGenWgt);
	llpgana.SetGenSigPerfectFlag(genSigPerfect);
    llpgana.kucmsAodSkimmer_local( listdir, eosdir, infilename, outfilename );
    return 1;


}//<<>>int main ( int argc, char *argv[] )
