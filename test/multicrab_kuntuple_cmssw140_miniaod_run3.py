#! /usr/bin/env python

import os
from optparse import OptionParser

from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from httplib import HTTPException

def getOptions():
    """
    Parse and return the arguments provided by the user.
    """
    usage = ("Usage: %prog --crabCmd CMD [--workArea WAD --crabCmdOpts OPTS]"
             "\nThe multicrab command executes 'crab CMD OPTS' for each project directory contained in WAD"
             "\nUse multicrab -h for help")

    parser = OptionParser(usage=usage)

    parser.add_option('-c', '--crabCmd',
                      dest = 'crabCmd',
                      default = 'submit',
                      help = "crab command",
                      metavar = 'CMD')

    parser.add_option('-w', '--workArea',
                      dest = 'workArea',
                      default = 'myWorkSpace',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_option('-o', '--crabCmdOpts',
                      dest = 'crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD",
                      metavar = 'OPTS')

    (options, arguments) = parser.parse_args()

    if arguments:
        parser.error("Found positional argument(s): %s." % (arguments))
    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def docrab( dataset ):

    options = getOptions()

    # The submit command needs special treatment.
    if options.crabCmd == 'submit':

        # External files needed by CRAB
        #inputJSON    = 'certifications/Cert_Collisions2022_355100_362760_Golden.json'	
        #inputJSON    = 'certifications/Cert_Collisions2023_366442_370790_Golden.json'
        #inputJSON    = 'certifications/Cert_Collisions2024_378981_386951_Golden.json'
        inputJSON    = 'certifications/Cert_Collisions2025_391658_398903_Golden.json'
        #inputJSON    = 'None'

        filespjob = int(dataset[1])

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'llpgana_miniaod.py'
        config.JobType.pyCfgParams = None

        config.section_("Site")
        #config.Site.whitelist = ['T1_*','T2_US_*','T2_IT_*','T2_DE_*','T2_ES_*','T2_FR_*','T2_UK_*']
        config.Site.blacklist = ['T2_CH_CSCS','T2_ES_CIEMAT']

        #config.Data.partialDataset = True
        config.Data.inputDataset   = None
        #config.Data.lumiMask       = inputJSON    # Comment out for MC only set for data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.splitting     = 'Automatic' # data
        config.Data.splitting = 'FileBased' # Justin
        config.Data.unitsPerJob  =  filespjob
        #config.Data.unitsPerJob  =  10
        #config.Data.unitsPerJob  =  4 # File based : Wjets 70,100,200
        #config.Data.unitsPerJob  =  3 # File based : QCD 4 - EGamma 3
        #config.Data.unitsPerJob  =  2 # File based : Wjets 400+
        #config.Data.unitsPerJob  =  1
        #config.Data.splitting    = 'EventAwareLumiBased' # MC&Data Set unitsperjob correctly for dataset !!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.unitsPerJob   = 72500 # data  !!!!!! lumimask ?
        #config.Data.splitting    = 'LumiBased'
        #config.Data.unitsPerJob   = 14 # 14 for eg c 18ul ~ 2 files & 19 hours too much disk space -> 1 file worth & break up to keep under 2000 jobs
        #config.Data.unitsPerJob   =  24 # MET C 
        #config.Data.unitsPerJob   =  48 # MET D
        #config.Data.runRange = '369927-369999' #
            
    
        #config.JobType.allowUndistributedCMSSW = True
        config.JobType.maxJobRuntimeMin = 1800
        #config.JobType.maxJobRuntimeMin = 2280
        config.JobType.maxMemoryMB = 2250
        config.Data.publication    = False
        config.Site.storageSite    = 'T3_US_FNALLPC'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [[dataset[0]]]

        for inDO in inputDataAndOpts:

            config.Data.inputDataset     = inDO[0]
            #config.Data.userInputFiles = open("kuntuple_SMS-GlGl_mGl-2000_mN2-1950_mN1-1000_gam_N2ctau-0p1_AOD_v23.txt").readlines()

            print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
            print( 'Input dataset for Crab Job : ' )
            #print( inDO )
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
            primaryDataset = (inDO[0].split('/')[1])
            if '-v2' in inDO[0] :
                primaryDataset = primaryDataset + '_v2'
            #primaryDataset = (inDO[0].split('/')[1]).split('00_Tune')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('_13T')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('-sherpa')[0]
            print( primaryDataset )
            #runEra         = (inDO[0].split('/')[2]) # justin
            runEra         = (inDO[0].split('/')[2]).split('-P')[0]
            #runEra         = (inDO[0].split('/')[2]).split('BPixMiniAOD')[0]
            #runEra         = (inDO[0].split('/')[2]).split('RECO')[0]
            print( runEra )
            dataset        = inDO[0].split('/')[3]
            print( dataset )

            print( 'current KUCMSNtuple version : 34' ) 
            version = '_v34'
            # 25 Adds ECAL TimeError and GS info + Adds SV processing           
            # 26 Adds Muon Object and loose Muon/electron counting
            # 27 Updated SV object
            # 28 Updated SV object + sip3d cut on SV + lepton-SV matching
            # 29 Added IsoTrack obejct, updated gen matching for Muons & Electrons
            # 30 consoladated Track objects in single Tracks object w/ SC matching
            # 31 reduced track pt min to 5 unless SC match, ? fixed SC rechit ID bug ?

            ##trial          = "kucmsntuple_JetMET_Mini_R23_SVHPMet100" + version # 
            ##trial          = "kucmsntuple_JetMET_R23_Met150" + version # 
            ##trial          = "kucmsntuple_JetMET_R23_SVIPM100" + version # 
            ##trial          = "kucmsntuple_JetMET_R24_SVIPM100" + version # 
            ##trial          = "kucmsntuple_JetMET_R23_NoFilter" + version # 
            ##trial          = "kucmsntuple_JetMET_R23_MET100" + version # 
            trial          = "kucmsntuple_JetMET_R25_SVHPM100" + version # 

            #trial          = "kucmsntuple_TTXJets_R23_SVIPM100" + version #
            #trial          = "kucmsntuple_GJets_R23_SVHPM100_MiniAOD" + version
            #trial          = "kucmsntuple_QCD_R23_SVIPM100_MiniAOD" + version
            #trial          = "kucmsntuple_ZJets_R23_SVIPM100_MiniAOD" + version
            #trial          = "kucmsntuple_ZJets_Ext_R23_SVIPM100" + version
            #trial          = "kucmsntuple_WJets_R23_SVIPM100" + version
            

            print( 'processing for : ',trial )

            # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
            config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
            config.General.requestName   = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra

            fSVHPM100 = 'eventFilter=SVHPMet100'
            fSVHPM150 = 'eventFilter=SVHPMet150'
            fSVIPM125 = 'eventFilter=SVIPMet125'
            fSVIPM100 = 'eventFilter=SVIPMet100'
            fM100 = 'eventFilter=MET100'
            fM150 = 'eventFilter=MET150'
            fAL1P = 'eventFilter=AL1IsoPho'
            fIPM100 = 'eventFilter=IsoPhoMet100'
            fNone = 'eventFilter=None'
            fAL1E = 'eventFilter=AL1SelEle'
            fNpSC = 'eventFilter=AL1NpSC'
            fInvMetALP = 'eventFilter=InvMet100IP'
            fInvMetPho30 = 'eventFilter=InvMetPho30'

            geninfo = 'hasGenInfo=True'
            nogeninfo = 'hasGenInfo=False'
            mcrab = 'multicrab=True'

            mNoSV = 'doSV=False'
            mDoDEle = 'doDisEle=true'
            mDoETOnly = 'doECALTrackOnly=true'

            efilter = fSVHPM100
            #efilter = fNpSC
            ##efilter = fM100
            ##efilter = fIPM100
            ##efilter = fAL1E
            ##efilter = fSVIPM100 #!!!!!!!!!!!!!!!!
            #efilter = fNone # for EXO Jet Met
            print( 'using :', efilter )
            print( 'using :', inputJSON )
            #####print( 'using :', geninfo )
            print( 'with files/job :', filespjob )

#  -------  selsect PD/MC dependent paramters
#---DATA--------------------------------------------------------------------------------------------------------------------------

#>>>>>>>>>>>>>>>>>>> RUN 3
#>>>>>>>>>>>>>>>>>>> 
#>>>>>>>>>>>>>>>>>>>        #2022 #globalTag=124X_dataRun3_v15 # prompt
            #gt = 'globalTag=124X_dataRun3_v15'
            #isrun3 = 'runera=Run3_2022' 
            ###config.JobType.pyCfgParams   = [gt,mcrab,efilter,mDoETOnly]
#>>>>>>>>>>>>>>>>>>>        # 2022 & 2023 #globalTag=124X_dataRun3_v17 ## rereco
            #gt = 'globalTag=140X_dataRun3_v17' # rereco
            #isrun3 = 'runera=Run3_2023' 
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter,isrun3]
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter,isrun3,mNoSV]
#>>>>>>>>>>>>>>>>>>>        #2024 #globalTag=140X_dataRun3_v20 
            #gt = 'globalTag=140X_dataRun3_v20' 
            #isrun3 = 'runera=Run3_2024'
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter,isrun3] 
#>>>>>>>>>>>>>>>>>>>        #2025 #globalTag=150X_dataRun3_Prompt_v1
            gt = 'globalTag=150X_dataRun3_Prompt_v1' 
            isrun3 = 'runera=Run3_2025'
            config.JobType.pyCfgParams   = [gt,mcrab,efilter,isrun3] 

#---MC--------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------

#>>>>>>>>>>>>>>>>> Run3 MC 2023 PostBPix
            #gt = 'globalTag=130X_mcRun3_2023_realistic_postBPix_v2'
            #runinfo = 'runera=Run3_2023'
            #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter,runinfo]

#--------------------------------------------------------------------------------------------------------------------------

            # Submit.
            try:
                print( "Submitting for input dataset %s" % primaryDataset + '_' + runEra + '_' + dataset )
                crabCommand(options.crabCmd, config = config, *options.crabCmdOpts.split())
                os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
            #except HTTPException as hte:
            #    print( "Submission for input dataset %s failed: %s" % (inDO[0], hte.headers) )
            except ClientException as cle:
                print( "Submission for input dataset %s failed: %s" % (inDO[0], cle) )

    # All other commands can be simply executed.
    elif options.workArea:

        for dir in os.listdir(options.workArea):
            projDir = os.path.join(options.workArea, dir)
            if not os.path.isdir(projDir):
                continue
            # Execute the crab command.
            msg = "Executing (the equivalent of): crab %s --dir %s %s" % (options.crabCmd, projDir, options.crabCmdOpts)
            print( "-"*len(msg) )
            print( msg )
            print( "-"*len(msg) )
            try:
                crabCommand(options.crabCmd, dir = projDir, *options.crabCmdOpts.split())
            except HTTPException as hte:
                print( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, hte.headers) )
            except ClientException as cle:
                print( "Failed executing command %s for task %s: %s" % (options.crabCmd, projDir, cle) )


##33333333333333333333333333333333333333333333333333333333333

def run_multi():

    TuneCP5MP = 'TuneCP5_13TeV-madgraphMLM-pythia8'
    RunIISummer20UL18RECO = '/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1'
    tunecp5r3 = 'TuneCP5_13p6TeV_madgraphMLM-pythia8'
    run3summer = 'Run3Summer23BPixMiniAODv4-130X_mcRun3_2023_realistic_postBPix'

    runDataset = [

        #["/JetMET/Run2022C-19Dec2023-v1/MINIAOD",1],
        #["/JetMET/Run2022D-19Dec2023-v1/MINIAOD",2],
        #["/JetMET/Run2022E-19Dec2023-v1/MINIAOD",2],
        #["/JetMET/Run2022F-19Dec2023-v2/MINIAOD",2],
        #["/JetMET/Run2022G-19Dec2023-v1/MINIAOD",2],

        #["/JetMET0/Run2023B-19Dec2023-v1/MINIAOD",1],
        #["/JetMET0/Run2023C-19Dec2023-v1/MINIAOD",1],
        #["/JetMET0/Run2023D-19Dec2023-v1/MINIAOD",1],
        #["/JetMET1/Run2023B-19Dec2023-v1/MINIAOD",1],
        #["/JetMET1/Run2023C-19Dec2023-v1/MINIAOD",1],
        #["/JetMET1/Run2023D-19Dec2023-v1/MINIAOD",1],

        #['/JetMET0/Run2023B-19Dec2023-v1/AOD',1],
        #['/JetMET0/Run2023C-19Dec2023-v1/AOD',1],
        #['/JetMET0/Run2023D-19Dec2023-v1/AOD',1],
        #['/JetMET1/Run2023B-19Dec2023-v1/AOD',1],
        #['/JetMET1/Run2023C-19Dec2023-v1/AOD',1],
        #['/JetMET1/Run2023D-19Dec2023-v1/AOD',1],

        ['/JetMET0/Run2025C-PromptReco-v1/MINIAOD',1],
        ['/JetMET0/Run2025C-PromptReco-v2/MINIAOD',1],
        ['/JetMET0/Run2025D-PromptReco-v1/MINIAOD',1],
        #['/JetMET0/Run2025E-PromptReco-v1/MINIAOD',1],
        #['/JetMET0/Run2025F-PromptReco-v1/MINIAOD',1],
        #['/JetMET0/Run2025F-PromptReco-v2/MINIAOD',1],
        #['/JetMET0/Run2025G-PromptReco-v1/MINIAOD',1],

        #['/JetMET1/Run2025C-PromptReco-v1/MINIAOD',1],
        #['/JetMET1/Run2025C-PromptReco-v2/MINIAOD',1],
        #['/JetMET1/Run2025D-PromptReco-v1/MINIAOD',1],
        #['/JetMET1/Run2025E-PromptReco-v1/MINIAOD',1],
        #['/JetMET1/Run2025F-PromptReco-v1/MINIAOD',1],
        #['/JetMET1/Run2025F-PromptReco-v2/MINIAOD',1],
        #['/JetMET1/Run2025G-PromptReco-v1/MINIAOD',1],


        #['/QCD-4Jets_HT-1000to1200_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],
        #['/QCD-4Jets_HT-100to200_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],
        #['/QCD-4Jets_HT-1200to1500_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],
        #['/QCD-4Jets_HT-1500to2000_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/QCD-4Jets_HT-2000_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],
        #['/QCD-4Jets_HT-200to400_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/QCD-4Jets_HT-400to600_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],
        #['/QCD-4Jets_HT-40to70_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/QCD-4Jets_HT-600to800_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/QCD-4Jets_HT-70to100_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/QCD-4Jets_HT-800to1000_'+tunecp5r3+'/'+run3summer+'_v2-v4/MINIAODSIM',1],

        #['/Zto2Nu-4Jets_HT-100to200_'+tunecp5r3+'/'+run3summer+'_v2-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-1500to2500_'+tunecp5r3+'/'+run3summer+'_v2-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-1500to2500_'+tunecp5r3+'/'+run3summer+'_v6_ext1-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-200to400_'+tunecp5r3+'/'+run3summer+'_v2-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-2500_'+tunecp5r3+'/'+run3summer+'_v2-v3/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-2500_'+tunecp5r3+'/'+run3summer+'_v6_ext1-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-400to800_'+tunecp5r3+'/'+run3summer+'_v2-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-400to800_'+tunecp5r3+'/'+run3summer+'_v6_ext1-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-800to1500_'+tunecp5r3+'/'+run3summer+'_v2-v2/MINIAODSIM',1],
        #['/Zto2Nu-4Jets_HT-800to1500_'+tunecp5r3+'/'+run3summer+'_v6_ext1-v2/MINIAODSIM',1],

    ] 

    for dataset in runDataset :
        docrab( dataset )

run_multi()
