//////////////////////////////////////////////////////////////////////
// -*- C++ -*-
//
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//////////////////////////////////////////////////////////////////////

#include "KUCMSAodSkimmer_sv_cc.hh"
//#include "KUCMSAodSkimmer_cc_probe.hh"
//#include "KUCMSAodSkimmer_rh_cc.hh"
//#include "KUCMSAodSkimmer_v16_rh_cc.hh"
// ------------------------------------------- main function ------------------------------------------------------------
int main ( int argc, char *argv[] ){

    //if( argc != 4 ) { std::cout << "Insufficent arguments." << std::endl; }
    //else {
                const std::string listdir = "ntuple_master_lists/";
                //const string KUCMSAodSkimmer::eosdir = "root://cmseos.fnal.gov//store/user/jaking/";
                const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/";
                //const std::string eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/anazario/";
                //const std::string eosdir = "root://cmseos.fnal.gov//store/user/janguian/";
				//const std::string eosdir = "../";				

				bool useEvtGenWgt = true;
                //bool useEvtGenWgt = false;
				//bool doSV = true; // only base selection
				bool doSV = false; // must have a sig pho or a sv
                //bool genSigPerfect = true;
                bool genSigPerfect = false;
                bool hasGenInfo = true;
                //bool hasGenInfo = false;
                ////const std::string infilename = "KUCMS_Ntuple_Master_BG_Files_List.txt"; genSigPerfect = false;
                //const std::string infilename = "KUCMS_Ntuple_Master_IsoPho_GMSB_Files_List.txt"; genSigPerfect = true;
				//const std::string infilename = "KUCMS_Ntuple_Master_DataPD_Files_List.txt"; hasGenInfo = false;
                //const std::string infilename = "KUCMS_Ntuple_Master_BG_Met0_Files_List.txt"; genSigPerfect = false;
                //const std::string infilename = "KUCMS_Ntuple_Master_GMSB_Files_List.txt"; genSigPerfect = true
                //const std::string infilename = "KUCMS_Ntuple_Test_GMSB_Files_List.txt"; genSigPerfect = true;
                //const std::string infilename = "KUCMS_Ntuple_Master_JetHT_Files_List.txt"; hasGenInfo = false;
                //const std::string infilename = "KUCMS_Ntuple_Master_DEG_Files_List.txt"; hasGenInfo = false;
                //const std::string infilename = "KUCMS_Ntuple_Jusin_Signal_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_BG_IPM100_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_GoGoG_Files_List.txt";
                const std::string infilename = "KUCMS_Ntuple_Master_ZBG_SVIPM100_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_GoGoZ_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_SqSqG_Files_List.txt";
                //const std::string infilename = "KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt"; useEvtGenWgt = false;

                //const std::string outfilename = "_LLPGskim_v21_rjrvars.root"; // skim v21 use new process_name tags 
                //const std::string outfilename = "_LLPGskim_v22_rjrvars.root"; // skim v22 keeps 0 pho events 
                //const std::string outfilename = "_LLPGskim_v23_rjrvars.root"; // fixed 2 sig pho events 
                //const std::string outfilename = "_rjrskim_v26.root"; // jet veto on first 2 Sig Phos only
                //const std::string outfilename = "_rjrskim_v27.root"; // SV matching & SVs & lepton veto
                const std::string outfilename = "_rjrskim_v28.root"; // MET cleaning flags 

                int skipCnt = 0;// Normal
				//int skipCnt = 10;// Met0 GMSB
                //int skipCnt = 100;// Met0 QCD
                KUCMSAodSkimmer llpgana;
                llpgana.kucmsAodSkimmer( listdir, eosdir, infilename, outfilename, hasGenInfo, genSigPerfect, doSV, skipCnt, useEvtGenWgt );
    //}
    return 1;


}//<<>>int main ( int argc, char *argv[] )
