// -*- C++ -*-
//
// Package:    KUCMSNtupilizer
// Class:      KUCMSNtupilizer
//
/**\class KUCMSNtupilizer LLPgammaAnalyzer.cc LLPgammaAnalyzer/plugins/KUCMSNtupilizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

CMSSW install for KUCMSNtupilizer:

cmsrel CMSSW_13_0_X
cd CMSSW_13_0_X/src 
cmsenv 

git cms-init

( make sure the analysis package is cloned into its own folder inside scr/  example( timing/timing/…  ) )
( yes ..  its required by CMSSW )
( initial scram b -j8 must be in src/ ?  after maybe scram-ed  from plugins/ or test/ …. )


mkdir KUCMSNtupilizer
cd KUCMSNtupilizer/
git clone https://github.com/jking79/KUCMSNtupilizer.git
cd ../
scram b -j 8

To Run :

kuntuple code in KUNtupleizer/KUNtupleizer/plugins
cd KUNtupleizer/KUNtupleizer/test
edit llpgana_mc_aod.py to change analysis and meta paramters ( number of eventts to proecess, event logging  #, input files .... ) 
to run : 
cmsRun llpgana_mc_aod.py




When make and object keep in mind the order of processing :

    // ---  Collecting objests ( preprocessing object pruning ) ---------------------------------------
    // -------------------------------------------------------------------------------------------------
    // -- Process Event  ---------------------------------------    
    // ** extracted from disphoana : starting point **** not all functios/varibles defined ***************
    // ** for example only -- convert to nano?, use ewkino varibles for output, find rechit information ** 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////

    if( DEBUG ) std::cout << "***********************************************************************************************" << std::endl;

    // clear global event varibles 
    geVar.clear(); // floats

    // -------------------------------------------------------------------------------------------------
    // ---  Collecting objects ( preprocessing object pruning & fill global object vectors )------------
    // -------------------------------------------------------------------------------------------------

    if( DEBUG ) std::cout << "LoadEvent ObjMan" << std::endl;
    ObjMan.LoadEvent( iEvent, iSetup, geVar );

    //------------------------------------------------------------------------------------
    // ----   Object processing ----------------------------------------------------------
    //------------------------------------------------------------------------------------
    // call functions to process collections and fill tree varibles to be saved
    // varibles to be saved to ttree are declared in the header
    // use LoadEvent() for any processing that must be done before crosstalk 
    // use PostProcessEvent() for any processing that must be done after crosstalk <<<  Most work should be done here.

    if( DEBUG ) std::cout << "ProcessEvent ObjMan" << std::endl;
    ObjMan.ProcessEvent( geVar );
    ObjMan.PostProcessEvent( geVar );

    //------------------------------------------------------------------------------------
    //---- Object processing Completed ----------------------------------------------------------
    //------------------------------------------------------------------------------------

    // -- Fill output trees ------------------------------------------

    if( DEBUG ) std::cout << "Select Event and Fill Tree" << std::endl;
    if( selectedEvent() ) outTree->Fill();

    // -- EOFun ------------------------------------------------------
    //     #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    //     ESHandle<SetupData> pSetup;
    //     iSetup.get<SetupRecord>().get(pSetup);
    //     #endif


