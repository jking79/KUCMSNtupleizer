//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSkimmer_Main.hh"
//#include "KUCMSAodSkimmer_sv_cc.hh"
//#include "KUCMSAodSkimmer_cc_probe.hh"
//#include "KUCMSAodSkimmer_rh_cc.hh"
//#include "KUCMSAodSkimmer_v16_rh_cc.hh"
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

	for(int i = 0; i < argc; i++){

    	if(strncmp(argv[i],"--help", 6) == 0) 			{ hprint = true; }
        if(strncmp(argv[i],"-h", 2) == 0)     			{ hprint = true; }
        if(strncmp(argv[i],"-i", 2) == 0)     			{ i++; in_file = string(argv[i]); }
        if(strncmp(argv[i],"--input", 7) == 0)			{ i++; in_file = string(argv[i]); }
	    if(strncmp(argv[i],"--output", 8) == 0)			{ i++; outfilename = string(argv[i]); }
        if(strncmp(argv[i],"-o", 2) == 0)				{ i++; outfilename = string(argv[i]); }
		if(strncmp(argv[i],"--evtFirst", 6) == 0)		{ i++; evti = std::atoi(argv[i]); }
        if(strncmp(argv[i],"--evtLast", 6) == 0)		{ i++; evtj = std::atoi(argv[i]); }
		if(strncmp(argv[i],"--hasGenInfo", 12) == 0)	{ hasGenInfo = true; }
		if(strncmp(argv[i],"--genSigPerfect", 15) == 0)	{ genSigPerfect = true; }
		if(strncmp(argv[i],"--noSVorPho", 11) == 0)		{ noSVorPho = true; }
        if(strncmp(argv[i],"--dataSetKey", 12) == 0)	{ i++; key = string(argv[i]); }
        if(strncmp(argv[i],"--xsec", 6) == 0)			{ i++; xsec = std::stof(argv[i]); }
        if(strncmp(argv[i],"--gluinoMass", 12) == 0)	{ i++; glumass = std::stof(argv[i]); }
        if(strncmp(argv[i],"--N2Mass", 8) == 0)			{ i++; n2mass = std::stof(argv[i]); }
        if(strncmp(argv[i],"--timeCaliTag", 13) == 0)	{ i++; ttag = string(argv[i]); }
        if(strncmp(argv[i],"--MCweight", 10) == 0)		{ i++; mcw = std::stof(argv[i]); }
	}

	if(hprint){
    	cout << "Usage: " << argv[0] << " [options]" << endl;
        cout << "  options:" << endl;
        cout << "   --help(-h)                           print options" << endl;
        cout << "   --input(-i) [file]                   input filelist" << endl;
		cout << "   --output(-o) [file]                  output file tag" << endl;
		cout << "   --evtFirst [i] --evtLast [j]         skim from event i to event j (default evtFirst = evtLast = 0 to skim over everything)" << endl;
		cout << "   --hasGenInfo                         sample has gen info (default = false)" << endl;
		cout << "   --genSigPerfect                      turn on gensigperfect (default = false)" << endl;
		cout << "   --noSVorPho                          selection doesn't have SV or photon (default = false)" << endl;
        cout << "   --dataSetKey                         set dataset key" << endl;
        cout << "   --xsec                               set cross-section (default = 1)" << endl;
        cout << "   --gluinoMass                         set gluino mass (default = 0)" << endl;
        cout << "   --N2Mass                             set N2 mass (default = 0)" << endl;
        cout << "   --timeCaliTag                        set time calibration tag (default = r2_ul18(_mc))" << endl;
        cout << "   --MCweight                           set MC weight (default = 0)" << endl;
		return 0;
	}

	cout << "outfile tag " << outfilename << endl;
	outfilename = outfilename+".root";
    //const std::string listdir = "ntuple_master_lists/";
	//const std::string infilename = "KUCMS_Ntuple_Master_DataPD_Files_List.txt"; hasGenInfo = false;

	if(hasGenInfo && in_file.find("Data")){
		cout << "Warning: running over list of data samples " << in_file << " with hasGenInfo == true. Setting hasGenInfo = false." << endl;
		hasGenInfo = false;
	}//<<>>if(hasGenInfo && in_file.find("Data"))

    int skipCnt = 0; // used to skip files ( in tchian ) for fast processing - if( nFiles%skipCnt != 0 ) continue; --  disabled in code  --  
    std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";
    KUCMSAodSkimmer llpgana;
	llpgana.SetEventRange(evti, evtj);
	if(in_file.find("_Master_") != string::npos){
        llpgana.kucmsAodSkimmer_listsOfLists( eosdir, in_file, outfilename, hasGenInfo, genSigPerfect, noSVorPho, skipCnt, useEvtGenWgt );
	} else {
		if(in_file.find("AODSIM") != string::npos){
			mctype = 0;
			ttag += "_mc";
		} else {
			mctype = 1;
		}//<<>>if(in_file.find("AODSIM") != string::npos)
		//from master list
   		llpgana.SetDataSetKey(key);
    	llpgana.SetCrossSection(xsec);
    	llpgana.SetGluinoMass(glumass);
    	llpgana.SetN2Mass(n2mass);
    	llpgana.SetMCType(mctype);
    	llpgana.SetTimeCalibrationTag(ttag);
		llpgana.SetMCWeight(mcw);	
		//this method takes in 1 list at a time
		eosdir = eosdir+"KUCMSNtuple/";
        llpgana.kucmsAodSkimmer( eosdir, in_file, outfilename, hasGenInfo, genSigPerfect, noSVorPho, skipCnt, useEvtGenWgt );
	}//<<>>if(in_file.find("_Master_") != string::npos)
    return 1;


}//<<>>int main ( int argc, char *argv[] )
