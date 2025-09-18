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
        #inputJSON    = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
        inputJSON    = 'Cert_Collisions2022_355100_362760_Golden.json'	

        filespjob = int(dataset[1])

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

        config.section_("Site")
        #config.Site.whitelist = ['T1_*','T2_US_*','T2_IT_*','T2_DE_*','T2_ES_*','T2_FR_*','T2_UK_*']
        config.Site.blacklist = ['T2_CH_CSCS','T2_ES_CIEMAT']

        config.Data.partialDataset = True
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

        #config.JobType.allowUndistributedCMSSW = True
        #config.JobType.maxJobRuntimeMin = 1800
        #config.JobType.maxJobRuntimeMin = 2280
        #config.JobType.maxMemoryMB = 3000
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
            primaryDataset = (inDO[0].split('/')[1]).split('_13T')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('_Tune')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('-sherpa')[0]
            print( primaryDataset )
            #runEra         = (inDO[0].split('/')[2]) # justin
            runEra         = (inDO[0].split('/')[2]).split('-')[0]
            #runEra         = (inDO[0].split('/')[2]).split('RECO')[0]
            print( runEra )
            dataset        = inDO[0].split('/')[3]
            print( dataset )

            print( 'current KUCMSNtuple version : 31' ) 
            version = '_v31'
            # 25 Adds ECAL TimeError and GS info + Adds SV processing           
            # 26 Adds Muon Object and loose Muon/electron counting
            # 27 Updated SV object
            # 28 Updated SV object + sip3d cut on SV + lepton-SV matching
            # 29 Added IsoTrack obejct, updated gen matching for Muons & Electrons
            # 30 consoladated Track objects in single Tracks object w/ SC matching
            # 31 reduced track pt min to 5 unless SC match, ? fixed SC rechit ID bug ?

            #trial          = "kucmsntuple_MET_R22_AL1NpSC" + version #
            #trial          = "kucmsntuple_MET_R18_AL1NpSC_DEOnly" + version #
            ##trial          = "kucmsntuple_MET_R18_AL1SelPho_DEOnly" + version #
            ##trial          = "kucmsntuple_MET_R18_AL1SelEle_DEOnly" + version #
            ##trial          = "kucmsntuple_EGamma_R18_AL1SelEle_DEOnly" + version #

            #trial          = "kucmsntuple_TTXJets_R18_SVIPM100" + version #
            trial          = "kucmsntuple_GJets_R18_SVIPM100" + version
            #trial          = "kucmsntuple_QCD_R18_SVIPM100" + version
            #trial          = "kucmsntuple_ZJets_R18_SVIPM100" + version
            #trial          = "kucmsntuple_WJets_R18_SVIPM100" + version
            #trial          = "kucmsntuple_DTBoson_R18_SVIPM100"  + version
            #trial          = "kucmsntuple_DiPJBox_R18_SVIPM100" + version

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
            fNpSC = 'eventFilter=AL1NpSC'

            geninfo = 'hasGenInfo=True'
            mcrab = 'multicrab=True'

            mNoSV = 'doSV=False'
            mDoDEle = 'doDisEle=true'
            mDoETOnly = 'doECALTrackOnly=true'

            #efilter = fNpSC
            ##efilter = fM100
            ##efilter = fIPM100
            ##efilter = fAL1E
            efilter = fSVIPM100 #!!!!!!!!!!!!!!!!
            print( 'using :', efilter )
            print( 'using :', geninfo )
            print( 'with files/job :', filespjob )

#  -------  selsect PD/MC dependent paramters
#---DATA--------------------------------------------------------------------------------------------------------------------------

#>>>>>>>>>>>>>>>>>>>     #2017UL #globalTag=106X_dataRun2_v20
            #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v20','multicrab=True']
#>>>>>>>>>>>>>>>>>>>     #2018UL #globalTag=106X_dataRun2_v36
            #gt = 'globalTag=106X_dataRun2_v36'
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter]
            ###config.JobType.pyCfgParams   = [gt,mcrab,efilter,mDoETOnly]
#>>>>>>>>>>>>>>>>>>>        #2022 #globalTag=124X_dataRun3_v15
            #gt = 'globalTag=124X_dataRun3_v15'
            ###config.JobType.pyCfgParams   = [gt,mcrab,efilter,mDoETOnly]

#---MC--------------------------------------------------------------------------------------------------------------------------
#>>>>>      #MC RunIISummer20UL18RECO
            gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
            ##config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True']
            config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter]
            ##config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>      #MC GMSB RunIIFall17DRPremix  #globalTag=94X_mc2017_realistic_v14  #  <<< comment/uncomment lumi mask when using/!using MC
            ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True']
            ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>>>>>>>>>>>>>>>     #MC RunIIAutumn18DRPremix 102X_upgrade2018_realistic_v15  
            ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True']
            ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>      #MC Justin
            #gt = 'globalTag=94X_mc2017_realistic_v14'
            ##config.JobType.pyCfgParams   = [gt,mcrab,geninfo,skimNone]
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

        #['/MET/Run2022A-27Jun2023-v1/AOD',1],
        #['/MET/Run2022B-27Jun2023-v2/AOD',1],
        #['/MET/Run2022C-27Jun2023-v2/AOD',1],

        #['/EGamma/Run2018B-15Feb2022_UL2018-v1/AOD'],#72500 eventsperlumi ( 36.5K per file - 4268 files )
        #['/EGamma/Run2018C-15Feb2022_UL2018-v1/AOD'],#72500 eventsperlumi   ( 36.1K per file - 4051 files )
        #['/EGamma/Run2018D-15Feb2022_UL2018-v1/AOD'],#72500 eventsperlumi

        #['/WJetsToLNu_HT-70To100_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM',4],
        #['/WJetsToLNu_HT-100To200_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM',4],
        #['/WJetsToLNu_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext1-v3/AODSIM',4],
        #['/WJetsToLNu_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM',1],
        #['/WJetsToLNu_HT-600To800_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM',2],
        #['/WJetsToLNu_HT-800To1200_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM',1],
        #['/WJetsToLNu_HT-1200To2500_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM',2],
        #['/WJetsToLNu_HT-2500ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'_ext2-v3/AODSIM',2],

        #['/ZJetsToNuNu_HT-100To200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-600To800_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-800To1200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-1200To2500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],
        #['/ZJetsToNuNu_HT-2500ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v1/AODSIM'],

        #['/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],
        #['/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',1],
        #['/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',1],
        #['/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',2],
        #['/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',1],

        #['/GJets_HT-40To100_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM',2],
        #['/GJets_HT-100To200_'+TuneCP5MP+'/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',2],
        ['/GJets_HT-200To400_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM',4],
        ['/GJets_HT-400To600_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM',2],
        #['/GJets_HT-600ToInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM',2],

        #['/QCD_HT2000toInf_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT1500to2000_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT1000to1500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT50to100_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT100to200_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT700to1000_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT200to300_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],
        #['/QCD_HT300to500_'+TuneCP5MP+RunIISummer20UL18RECO+'-v2/AODSIM'],

        #['/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',2],
        #['/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',1],
        #['/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],

        #['/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],
        #['/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],
        #['/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',1],
        #['/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',1],
        #['/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],
        #['/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],
        #['/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',1],
        #['/WGG_5f_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',1],

    ] 

    for dataset in runDataset :
        docrab( dataset )

run_multi()
