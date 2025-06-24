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
        ##inputJSON    = 'golden2016.json'
        ##inputJSON    = 'golden2017.json'
        ##inputJSON    = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
        #inputJSON    = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        ##inputJSON    = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
        ##inputJSON    = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
        inputJSON    = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
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
        #config.Data.unitsPerJob  =  10
        #config.Data.splitting    = 'EventAwareLumiBased' # MC&Data Set unitsperjob correctly for dataset !!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.unitsPerJob   = 45000 # data  !!!!!! lumimask ?
        #config.Data.unitsPerJob  =  1500 # MC GMSB
        ##config.Data.unitsPerJob  =  20000 # MC GJet
        #config.Data.unitsPerJob  =  15000 # MC QCD - WJetsToLNu
        #config.Data.unitsPerJob  =  25000 # MC QCD - WJetsToLNu large
        #config.Data.unitsPerJob  =  30000 # MC TTJet
        #config.Data.unitsPerJob  =  750 # MC DiPhoBox (DPB)
        #config.Data.unitsPerJob  =  1 # File based ZJets
        config.Data.unitsPerJob  =  4 # File based QCD

        config.JobType.allowUndistributedCMSSW = True
        #config.JobType.maxJobRuntimeMin = 1800
        config.JobType.maxJobRuntimeMin = 2280
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
            #primaryDataset = (inDO[0].split('/')[1])
            #primaryDataset = (inDO[0].split('/')[1]).split('_13T')[0]
            primaryDataset = (inDO[0].split('/')[1]).split('_Tune')[0]
            print( primaryDataset )
            #runEra         = (inDO[0].split('/')[2]) # justin
            #runEra         = (inDO[0].split('/')[2]).split('-')[0]
            runEra         = (inDO[0].split('/')[2]).split('RECO')[0]
            print( runEra )
            dataset        = inDO[0].split('/')[3]
            print( dataset )

            print( 'current KUCMSNtuple version : 27' ) 
            # 25 Adds ECAL TimeError and GS info + Adds SV processing           
            # 26 Adds Muon Object and loose Muon/electron counting
            # 27 Updated SV object


            #trial          = "kucmsntuple_WJetsToLNu_R18_IPM100_v25" #
            #trial          = "kucmsntuple_TTXJets_R18_IPM100_v25" #
            #trial          = "kucmsntuple_GJets_R18_IPM100_v25"

            trial          = "kucmsntuple_QCD_R18_IPM100_v27"
            trial          = "kucmsntuple_QCD_R18_ReTry_IPM100_v27"
            #trial          = "kucmsntuple_ZJets_R18_SVIPM100_v27"

            print( 'processing for : ',trial )

            # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
            config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
            config.General.requestName   = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+dataset+"_"+runEra

            fSVIPM125 = 'eventFilter=SVIPMet125'
            fSVIPM100 = 'eventFilter=SVIPMet100'
            fM100 = 'eventFilter=MET100'
            fAL1P = 'eventFilter=AL1IsoPho'
            fIPM100 = 'eventFilter=IsoPhoMet100'
            fNone = 'eventFilter=None'
            fAL1E = 'eventFilter=AL1SelEle'

            geninfo = 'hasGenInfo=True'
            mcrab = 'multicrab=True'

            efilter = fIPM100
            #efilter = fSVIPM125 #!!!!!!!!!!!!!!!!!!!!!11
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
            gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
            #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True']
            config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter]
            #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>      #MC GMSB RunIIFall17DRPremix  #globalTag=94X_mc2017_realistic_v14  #  <<< comment/uncomment lumi mask when using/!using MC
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True']
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>>>>>>>>>>>>>>>     #MC RunIIAutumn18DRPremix 102X_upgrade2018_realistic_v15  
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True']
            #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>      #MC Justin
            #gt = 'globalTag=94X_mc2017_realistic_v14'
            #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,skimNone]
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

    TuneCP5MP = 'TuneCP5_13TeV-madgraphMLM-pythia8'
    RunIISummer20UL18RECO = '/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1'

    runDataset = [

        #['/MET/Run2018A-15Feb2022_UL2018-v1/AOD'],
        #['/MET/Run2018B-15Feb2022_UL2018-v1/AOD'],
        #['/MET/Run2018C-15Feb2022_UL2018-v1/AOD'],
        #['/MET/Run2018D-15Feb2022_UL2018-v1/AOD'],

        #['/WJetsToLNu_HT-70To100_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM'],
        #['/WJetsToLNu_HT-100To200_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM'],
        #['/WJetsToLNu_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM'],
        #['/WJetsToLNu_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM'],
        #['/WJetsToLNu_HT-600To800_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM'],
        #['/WJetsToLNu_HT-800To1200_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM'],
        #['/WJetsToLNu_HT-1200To2500_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM'],
        #['/WJetsToLNu_HT-2500ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM'],

        #['/ZJetsToNuNu_HT-100To200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-1200To2500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-2500ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-600To800_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-800To1200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],

        #['/TTJets_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/ttWJets_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ttZJets_+TuneCP5MP+/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM'],
        #['/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM'],
        #['/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM'],

        #['/GJets_HT-40To100_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/GJets_HT-100To200_'+TuneCP5MP+'/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM'],
        #['/GJets_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/GJets_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/GJets_HT-600ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],

        #['/QCD_HT2000toInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT1500to2000_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT1000to1500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT50to100_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT100to200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT700to1000_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT200to300_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        ['/QCD_HT300to500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],

    ] 

    for dataset in runDataset :
        docrab( dataset )

run_multi()
