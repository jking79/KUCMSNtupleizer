// -*- C++ -*-
//
// Package:    KUCMSNtupleizer
// Class:      KUCMSNtupleizer
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

KUCMSNtupleizer is a plugin module for CMSSW
To utilize first setup a CMSSW version and cone the git repo for KUCMSNtupleizer into CMSSW
Installation instructions using CMSSW version 13.3.0 below :

-------------------------------------------------------------------------------------------------------------------------------------------

CMSSW install for KUCMSNtupleizer:

cmsrel CMSSW_13_3_0
cd CMSSW_13_3_0/src/
cmsenv
git cms-init

mkdir KUCMSNtupleizer
cd KUCMSNtupleizer/
git clone https://github.com/jking79/KUCMSNtupleizer.git
cd ../
scram b -j 8

( make sure the analysis package is cloned into its own folder inside scr/  example( timing/timing/  ) )
( in this case after the git clone you will have src/KUCMSNtupleizer/KUCMSNtupleizer/... )
( yes ..  its required by CMSSW )
( initial scram b -j8 must be in src/ ?  after may be scram-ed  from plugins/ or test/ )

----------------------------------------------------------------------------------------------------------------------------------------------

To Run Locally:

cd KUNtupleizer/KUNtupleizer/test
edit llpgana_mc_aod.py to change analysis and meta paramters ( number of eventts to proecess, event logging  #, input files .... ) 
to run : 
cmsRun llpgana_mc_aod.py

To Run on Crab :

Utilize crab submission script : KUNtupleizer/KUNtupleizer/test/multicrab_kuntuple_cmssw133_aod.py
This is a customized scrpit I have modifed from the code I origanlly started working with many years ago. It was orginally a script intended to run a crab command on each 
crab job in a specifed work directory.  This script defaults to running the "submit" crab command with the configuration set in the multicrab_kuntuple_cmssw133_aod.py script.  
It can still be used to run batch crab commands via the -c command line input, though I use the unmodified script multicrab.py for this.

The script may/will need to be customized for the jobs you want to run.  Following is a quick guide of the important parameters to look at :

On line 76 the working are where the log and jobs files are saved is specified via an input paramter that defalts to myWorkSpace.  

Starting on line 53 are the lumi jason filters for data.  Inusre the correct json is used for the era your processing, these are not used for MC.  
The json is loaded into the crab submission in line 85.  This should be commented out for MC.

On lines 86 - 95 are parmaters for the jobs sizes. This is a bit of an art.  The bigger the jobs are the more resources and time they take to run and you run into limits
on the maximuin time, RAM, ect you are allowed to use.  The smaller the jobs are the more jobs each submission has. We are limited to less then 10K jobs per submission, 
and very small jobs under utilze CPU resources which is frowned upon. There is an 'Automatic' option which is suppose to run tests and set the jobs sizes for you,
but I have had problems with utilizing this feature, so I use 'EventAwareLumiBased' for most data & MC, and 'FileBased' for some select MC.  
If you use these two options we need to set unitsPerJob which can take sime trial and error to find a setting that works well. See script for previously used values. 

In lines 105 to 146 I set request names and output path information based on premutations of the input datasets. 
On lines 107 and 108 choose to use an input file with a list of file names or a list defined later in the script for the DAS file batch to process.   
On lines 114 - 122 I do some string splitting for later use in various names, there is a limit on name length so at times we need to reduce the strign lengths by splitting
In lines 125 - 137 I set a "trial" name that I use in the output directory stucture.
On line 142 we set the output directory to use on LPC. Change this to what you want to use for the output directory 
Starting of line 146 the input paramters for the script running KUCMSNtupleizer are set.  There are four paramters to note :
	example : config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
	globalTag ( CMSSW ) : needs to be set to the approrate tag for the PD/MC your are proccessing, I will try to list this on the spreadsheet
	multicrab ( ours ) : should be set true for all crab script submissions ( used in script to unselect some options used in local running only )
	hasGenInfo ( ours ) : set to true for MC and false for PD 
	eventSkim ( ours ) : select what skim type to use, defaults to "MET100" if not specified, "none" for none, "AL1IsoPho" for at least 1 isolated photon 
Starting on line 207 specify the input CMSDAS datasets to run over in a list.

To run submission script with :

	python3 multicrab_kuntuple_cmssw133_aod.py 
or
	python3 multicrab_kuntuple_cmssw133_aod.py -w yourworkspace

----------------------------------------------------------------------------------------------------

Adding objects to KUCMSNtupleizer :

To modify or extend KUCMSNtupleizer you will need to modify an existiing object class or add a new object class
This code is in KUNtupleizer/KUNtupleizer/plugins
Objct classes are created by inheratance from KUCMSObjectBase.hh - please read instruction in KUCMSObjectBase.hh on how to do this. 
All object classes are substanuated and loaded into the object manager ObjMan in the KUCMSNtupleizer constructor.
 -- add instructions on adding objects to coinstructor --

When adding or modify an object class keep in mind the order of processing :

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


