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
                      default = 'mySigWorkSpace',
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
        ##inputJSON    = 'golden2016.json'
        ##inputJSON    = 'golden2017.json'
        ##inputJSON    = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
        #inputJSON    = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        ##inputJSON    = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
        ##inputJSON    = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
        #inputJSON    = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
        #inputJSON    = 'Cert_Collisions2022_355100_362760_Golden.json'	

        #--------------------------------------------------------
        # This is the base config:
        #--------------------------------------------------------
        from CRABClient.UserUtilities import config
        config = config()

        config.General.workArea    = options.workArea
        config.General.requestName = None

        config.JobType.pluginName  = 'Analysis'
        config.JobType.psetName    = 'llpgana_mc_aod.py'
        config.JobType.pyCfgParams = None

        config.Data.partialDataset = True
        config.Data.inputDataset   = None
        #config.Data.lumiMask       = inputJSON    # Comment out for MC only set for data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.splitting     = 'Automatic' # data
        config.Data.splitting = 'FileBased' # Justin
        config.Data.unitsPerJob  =  1
        #config.Data.splitting    = 'EventAwareLumiBased' # MC&Data Set unitsperjob correctly for dataset !!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.unitsPerJob   = 45000 # data  !!!!!! lumimask ?
        #config.Data.unitsPerJob  =  1500 # MC GMSB
        #config.Data.unitsPerJob  =  10000 # MC GJet
        #config.Data.unitsPerJob  =  15000 # MC QCD - WJetsToLNu
        #config.Data.unitsPerJob  =  50000 # MC TTJet
        #config.Data.unitsPerJob  =  750 # MC DiPhoBox (DPB)

        #config.section_("Site")
        config.Site.whitelist = ['T1_*','T2_US_*','T2_IT_*','T2_DE_*','T2_FR_*','T2_UK_*']
        config.Site.blacklist = ['T2_CH_CSCS','T2_ES_CIEMAT']

        config.JobType.allowUndistributedCMSSW = True
        config.Data.publication    = False
        config.Site.storageSite    = 'T3_US_FNALLPC'
        #--------------------------------------------------------

        # Will submit one task for each of these input datasets.
        inputDataAndOpts = [dataset]
        print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
        print( inputDataAndOpts )
       
        for inDO in inputDataAndOpts : 

            config.Data.userInputFiles = open(inDO[1]).readlines()

            print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
            print( 'Input dataset for Crab Job : ' )
            #print( inDO )
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
            primaryDataset = (inDO[0].split('/')[1])
            #primaryDataset = (inDO[0].split('/')[1]).split('_13T')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('_Tune')[0]
            print( primaryDataset )
            runEra         = (inDO[0].split('/')[2]) # justin
            #runEra         = (inDO[0].split('/')[2]).split('-')[0]
            #runEra         = (inDO[0].split('/')[2]).split('RECO')[0]
            print( runEra )
            dataset        = inDO[0].split('/')[3]
            print( dataset )

            print( 'with : ',inDO[1] )
            #print( 'current KUCMSNtuple version : 31' )            

            #trial          = "kucmsntuple_gogoG_Sig_IPM100_v26" #
            #trial          = "kucmsntuple_gogoZ_Sig_IPM100_v27" # 
            #trial          = "kucmsntuple_gogoZ_Sig_SVIPM100_v27" # filter on #SVs( lep or had ) - redefinded loose electron
            #trial          = "kucmsntuple_sqsqG_Sig_SVIPM100_v27" #
            #trial          = "kucmsntuple_SMS_Sig_SVIPM100_v31" #
            #trial          = "kucmsntuple_MET_R17_SVIPM100_p9_v31"# 
            trial          = "kucmsntuple_SMS_Sig_SVIPM100_v33" #

            # set trial name - used in output path ?

            print( 'processing for : ',trial )

            # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
            config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
            config.General.requestName   = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra

            fSVIPM100 = 'eventFilter=SVIPMet100'
            fM100 = 'eventFilter=MET100'
            fAL1P = 'eventFilter=AL1IsoPho'
            fIPM100 = 'eventFilter=IsoPhoMet100'
            fNone = 'eventFilter=None'
            fAL1E = 'eventFilter=AL1SelEle'
            fHardPho='eventFilter=Pho30'

            geninfo = 'hasGenInfo=True'
            mcrab = 'multicrab=True'

            runinfo = 'runera=Run3'
            doMetFilters = 'runMETFilters=False'
            doSV = 'doSV=False'

            
            #efilter = fIPM100
            efilter = fSVIPM100 #!!!!!!!!!!!!!!!!!!!!!11
            print( 'using :', efilter )
            print( 'using :', geninfo )


#  -------  selsect PD/MC dependent paramters
#-----------------------------------------------------------------------------------------------------------------------------

#>>>>>>>>>>>>>>>>>>>     #2017UL #globalTag=106X_dataRun2_v20
            #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v20','multicrab=True']
#>>>>>>>>>>>>>>>>>>>     #2018UL #globalTag=106X_dataRun2_v36
            #gt = 'globalTag=106X_dataRun2_v36'
            #config.JobType.pyCfgParams   = [gt,mcrab,skimIPM100]
            #config.JobType.pyCfgParams   = [gt,mcrab,skimM100]
            #config.JobType.pyCfgParams   = [gt,mcrab,skimAL1P]
            #config.JobType.pyCfgParams   = [gt,mcrab,skimNone]
            #config.JobType.pyCfgParams   = [gt,mcrab,skimAL1E]
#>>>>>>>>>>>>>>>>>>>        #2022 #globalTag=124X_dataRun3_v15
            #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_v15','multicrab=True']

#-----------------------------------------------------------------------------------------------------------------------------
#>>>>>      #MC RunIISummer20UL18RECO
            #gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
            #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True']
            #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,skimIPM100]
            #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
            #gt = 'globalTag=106X_dataRun2_v37'
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter]
#>>>>>      #MC GMSB RunIIFall17DRPremix  #globalTag=94X_mc2017_realistic_v14  #  <<< comment/uncomment lumi mask when using/!using MC
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True']
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>>>>>>>>>>>>>>>     #MC RunIIAutumn18DRPremix 102X_upgrade2018_realistic_v15  
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True']
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>      #MC Justin
            #gt = 'globalTag=124X_mcRun3_2022_realistic_v12'
            gt = 'globalTag=124X_mcRun3_2022_realistic_postEE_v1'
            config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter,runinfo,doMetFilters,doSV]
#-----------------------------------------------------------------------------------------------------------------------------

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

    inlist = 'sig_aod_lists/filelist_SMS-'
    runDataset = [

        #['/SMS-GlGl/mGl-2300_mN2-2200_mN1-2100/AODSIM',inlist+'GlGl_mGl-2300_mN2-2200_mN1-2100_GZ_N2ctau-0p5_FASTAOD.txt'],
        ['/SMS-GlGl/mGl-2300_mN2-2250_mN1-2200_ct-0p5/AODSIM',inlist+'GlGl_mGl-2300_mN2-2250_mN1-2200_GZ_N2ctau-0p5_FASTAOD.txt'],
        #['/SMS-GlGl/mGl-1500_mN2-500_mN1-100/AODSIM',inlist+'GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1900_mN1-1000/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1000_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1900_mN1-1500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1500_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1900_mN1-1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1900_mN1-250/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1900_mN1-500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-500_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-1000/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1000_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-1500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1500_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-1900/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1900_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-250/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-250_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1950_mN1-500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-500_gam_N2ctau-0p1_AOD_v23.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1500_mN1-1000/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-1000_gam_N2ctau-0p1_AOD_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1500_mN1-1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-1_gam_N2ctau-0p1_AOD_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1500_mN1-250/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-250_gam_N2ctau-0p1_AOD_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1500_mN1-500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-500_gam_N2ctau-0p1_AOD_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1000_mN1-1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-1_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1000_mN1-250/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-250_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-1000_mN1-500/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-500_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-500_mN1-1/AODSIM',inlist+'GlGl_mGl-2000_mN2-500_mN1-1_v24.txt'],
        #['/SMS-GlGl/mGl-2000_mN2-500_mN1-250/AODSIM',inlist+'GlGl_mGl-2000_mN2-500_mN1-250_v24.txt'],
        #['/SMS-GlGl/mGl-2500_mN2-1500_mN1-1000/AODSIM',inlist+'GlGl_mGl-2500_mN2-1500_mN1-1000_v24.txt'],

        #['/SMS-GlGlZ/mGl-1500_mN2-500_mN1-100-ct0p1/AODSIM',inlist+'GlGl_mGl-1500_mN2-500_mN1-100_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-200-ct0p001/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p001_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-200-ct0p1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-200-ct0p3/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-350-ct0p001/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p001_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-350-ct0p1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1900_mN1-350-ct0p3/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p3_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-1950_mN1-1900-ct0p1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1900_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-200-ct0p001/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p001_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-200-ct0p1/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-200-ct0p3/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p3_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-350-ct0p001/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-350-ct0p1/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p1_AOD_v26.txt'],
        #['/SMS-GlGlZ/mGl-2000_mN2-400_mN1-350-ct0p3/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p3_AOD_v26.txt'],

        #['/SMS-SqSq/mGl-1700_mN2-1500_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-1700_mN2-1500_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-1700_mN2-300_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-1700_mN2-300_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-1850_mN2-1650_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-1850_mN2-1650_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-1850_mN2-300_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-1850_mN2-300_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-2000_mN2-1800_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-2000_mN2-1800_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-2000_mN2-300_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-2000_mN2-300_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-2150_mN2-300_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-2150_mN2-300_mN1-100_v27.txt'],
        #['/SMS-SqSq/mGl-2150_mN2-1950_mN1-100_ct0p1/AODSIM',inlist+'SqSq_mSq-2150_mN2-1950_mN1-100_v27.txt'],

        #['/SMS-GlGlGZ/mGl-2500_mN2-2400_mN1-500_ct0p1/AODSIM','sig_aod_lists/SMS-GlGl_mGl-2500_mN2-2400_mN1-500_GZ_N2ctau-0p1_AOD.txt'],

        #['/MET/Run2017F-09Aug2019_UL2017_rsb-v1/AOD','fileslist_missing_17F.txt'],

        #['/SMS-GlGl-GZ/2300_N2-1300_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-1300_mN1-1000_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-1300_N1-1000_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-1300_mN1-1000_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-1600_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-1600_mN1-1000_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-1600_N1-1000_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-1600_mN1-1000_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-1600_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-1600_mN1-500_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-1600_N1-500_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-1600_mN1-500_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2200_N1-2100_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-2200_mN1-2100_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2200_N1-2100_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-2200_mN1-2100_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2200_N1-2150_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-2200_mN1-2150_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2200_N1-2150_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-2200_mN1-2150_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2250_N1-2150_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-2250_mN1-2150_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2250_N1-2150_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-2250_mN1-2150_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2250_N1-2200_ct1/AODSIM',inlist+'GlGl_mGl-2300_mN2-2250_mN1-2200_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2300_N2-2250_N1-2200_ct5/AODSIM',inlist+'GlGl_mGl-2300_mN2-2250_mN1-2200_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-1200_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-1200_mN1-500_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-1200_N1-500_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-1200_mN1-500_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2000_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2000_mN1-1000_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2000_N1-1000_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2000_mN1-1000_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2000_N1-1500_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2000_mN1-1500_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2000_N1-1500_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2000_mN1-1500_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2400_N1-2300_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2400_mN1-2300_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2400_N1-2300_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2400_mN1-2300_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2400_N1-2350_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2400_mN1-2350_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2400_N1-2350_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2400_mN1-2350_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2450_N1-2350_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2450_mN1-2350_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2450_N1-2350_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2450_mN1-2350_GZ_N2ctau-0p5_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2450_N1-2400_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-2450_mN1-2400_GZ_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GZ/2500_N2-2450_N1-2400_ct5/AODSIM',inlist+'GlGl_mGl-2500_mN2-2450_mN1-2400_GZ_N2ctau-0p5_AOD.txt'],

        #['/SMS-GlGl-GG/1500_N2-1000_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-1500_mN2-1000_mN1-500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/1500_N2-500_N1-100_ct1/AODSIM',inlist+'GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1000_N1-1_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-1_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1000_N1-250_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-250_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1000_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1000_mN1-500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1500_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-1000_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1500_N1-1_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-1_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1500_N1-250_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-250_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1500_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1500_mN1-500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1900_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1000_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1900_N1-1500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1900_N1-1_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-1_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1900_N1-250_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1900_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-1000_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1000_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-1500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-1900_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1900_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-1_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-250_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-250_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-1950_N1-500_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-500_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-500_N1-1_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-500_mN1-1_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2000_N2-500_N1-250_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-500_mN1-250_gam_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-GG/2500_N2-1500_N1-1500_ct1/AODSIM',inlist+'GlGl_mGl-2500_mN2-1500_mN1-1000_gam_N2ctau-0p1_AOD.txt'],

        #['/SMS-GlGl-ZZ/1500_N2-500_N1-100_ct1/AODSIM',inlist+'GlGl_mGl-1500_mN2-500_mN1-100_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/1500_N2-500_N1-100_ll_ct1/AODSIM',inlist+'GlGl_mGl-1500_mN2-500_mN1-100_Zll_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-200_ct001/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p001_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-200_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-200_ct3/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-350_ct001/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p001_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-350_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1900_N1-350_ct3/AODSIM',inlist+'GlGl_mGl-2000_mN2-1900_mN1-350_Zff_N2ctau-0p3_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-1950_N1-1900_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-1950_mN1-1900_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-200_ct001/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p001_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-200_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-200_ct3/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p3_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-350_ct001/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-350_ct1/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p1_AOD.txt'],
        #['/SMS-GlGl-ZZ/2000_N2-400_N1-350_ct3/AODSIM',inlist+'GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p3_AOD.txt'],


    ] 

    for dataset in runDataset :
        docrab( dataset )

run_multi()
