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
        #inputJSON    = '' # MC
        inputJSON    = 'certifications/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt'
        #inputJSON    = 'certifications/Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        #inputJSON    = 'certifications/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
        #inputJSON    = 'certifications/Cert_Collisions2022_355100_362760_Golden.json'	
        #inputJSON    = 'certifications/Cert_Collisions2023_366442_370790_Golden.json'

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
        config.Data.lumiMask       = inputJSON    # !!!!!!!!!!!!  Comment out for MC only set for data !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        #config.Data.splitting     = 'Automatic' # data
        config.Data.splitting = 'FileBased' # Justin
        config.Data.unitsPerJob  =  filespjob
        #config.Data.runRange = '306458-306460' #
            
    
        #config.JobType.allowUndistributedCMSSW = True
        #config.JobType.maxJobRuntimeMin = 2500
        config.JobType.maxJobRuntimeMin = 2280
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
            print( inDO )
            # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
            primaryDataset = (inDO[0].split('/')[1])
            #primaryDataset = (inDO[0].split('/')[1]).split('_13T')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('_Tune')[0]
            #primaryDataset = (inDO[0].split('/')[1]).split('-sherpa')[0]
            print( primaryDataset )
            runEra         = (inDO[0].split('/')[2]) # justin
            runEra         = (inDO[0].split('/')[2]).split('_MiniAOD')[0]
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

            #trial          = "kucmsntuple_MET_R22_AL1NpSC" + version #
            #trial          = "kucmsntuple_MET_R18_AL1NpSC_DEOnly" + version #
            ##trial          = "kucmsntuple_MET_R18_AL1SelPho_DEOnly" + version #
            ##trial          = "kucmsntuple_MET_R18_AL1SelEle_DEOnly" + version #
            ##trial          = "kucmsntuple_EGamma_R18_AL1SelEle_DEOnly" + version #

            #trial          = "kucmsntuple_JetMET_R22_SVHPM100_MiniAOD" + version # 
 
            #trial          = "kucmsntuple_MET_R18_SVHPM100_MiniAOD" + version #  
            #trial          = "kucmsntuple_MET_R17_SVHPM100_MiniAOD" + version #  
            trial          = "kucmsntuple_MET_R16_SVHPM100_MiniAOD" + version # 

            #trial          = "kucmsntuple_MET_R17_SVIPM100_p8" + version #  
            #####trial          = "kucmsntuple_JETMET_R24_SVIPM100" + version # 

            #trial          = "kucmsntuple_test" + version 

            #trial          = "kucmsntuple_EGamma_R18_InvMetPho30_NoSV"
            #trial          = "kucmsntuple_EGamma_R18_InvMet100IP_NoSv" + version # 
            #trial          = "kucmsntuple_DisJet_R18_SVIPM100" + version #
            #trial          = "kucmsntuple_DisJet_R18_SVIPM100_nolumimask" + version #

            #trial          = "kucmsntuple_TTXJets_R18_SVIPM100" + version #
            #trial          = "kucmsntuple_GJets_R18_SVHPM100_MiniAOD" + version
            #trial          = "kucmsntuple_QCD_R18_SVIPM100" + version
            #trial          = "kucmsntuple_ZJets_R18_SVIPM100" + version
            #trial          = "kucmsntuple_WJets_R18_SVIPM100" + version
            #trial          = "kucmsntuple_DTBoson_R18_SVIPM100"  + version
            #trial          = "kucmsntuple_DiPJBox_R18_SVIPM100" + version

            print( 'processing for : ',trial )

            # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
            config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
            config.General.requestName   = trial+"_"+primaryDataset+"_"+runEra+"_request"
            config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+runEra

            fSVHPM100 = 'eventFilter=SVHPMet100'
            fSVHPM150 = 'eventFilter=SVHPMet150'
            fSVIPM125 = 'eventFilter=SVIPMet125'
            fSVIPM100 = 'eventFilter=SVIPMet100'
            fM100 = 'eventFilter=MET100'
            fAL1P = 'eventFilter=AL1IsoPho'
            fIPM100 = 'eventFilter=IsoPhoMet100'
            fNone = 'eventFilter=None'
            fAL1E = 'eventFilter=AL1SelEle'
            fNpSC = 'eventFilter=AL1NpSC'
            fInvMetALP = 'eventFilter=InvMet100IP'
            fInvMetPho30 = 'eventFilter=InvMetPho30'

            geninfo = 'hasGenInfo=True'
            mcrab = 'multicrab=True'

            mNoSV = 'doSV=False'
            mDoDEle = 'doDisEle=true'
            mDoETOnly = 'doECALTrackOnly=true'

            #efilter = fNpSC
            ##efilter = fM100
            ##efilter = fIPM100
            ##efilter = fAL1E
            efilter = fSVHPM100 #!!!!!!!!!!!!!!!!
            print( 'using :', efilter )
            print( 'using :', geninfo )
            print( 'with files/job :', filespjob )

#  -------  selsect PD/MC dependent paramters
#---DATA--------------------------------------------------------------------------------------------------------------------------

            # recomended UL in 13.3
            #2018 UL: 106X_dataRun2_v35
            #2017 UL: 106X_dataRun2_v33
            #2016 UL: 106X_dataRun2_v33 (or 102X_dataRun2_v13 depending on CMSSW version; newer CMSSW 13.x usually uses 106X_dataRun2_v33)
            #Recomended for Run2 ana - Data(ALL UL MiniAODv2): 106X_dataRun2_v37 MC 2018: 106X_upgrade2018_realistic_v16_L1v1
#>>>>>>>>>>>>>>>>>>>     #2016UL 
#>>>>>>>>>>>>>>>>>>>     #2017UL 
#>>>>>>>>>>>>>>>>>>>     #2018UL #globalTag=106X_dataRun2_v37 # https://twiki.cern.ch/twiki/bin/view/CMS/PdmVRun2LegacyAnalysis
            gt = 'globalTag=106X_dataRun2_v37'
            #runinfo = 'runera=Run2_2018'
            #runinfo = 'runera=Run2_2017'
            runinfo = 'runera=Run2_2016'
            config.JobType.pyCfgParams   = [gt,mcrab,efilter,runinfo]
            ###config.JobType.pyCfgParams   = [gt,mcrab,efilter,mDoETOnly]
            ##config.JobType.pyCfgParams   = [gt,mcrab,fInvMetALP,mNoSV]
            ##config.JobType.pyCfgParams   = [gt,mcrab,fInvMetPho30,mNoSV]
#>>>>>>>>>>>>>>>>>>>        #2022 #globalTag=124X_dataRun3_v15
            #gt = 'globalTag=124X_dataRun3_v15'
            #runinfo = 'runera=Run3'
            #runinfo = 'runera=Run3_2022'
            #config.JobType.pyCfgParams   = [gt,mcrab,efilter,runinfo]
            ###config.JobType.pyCfgParams   = [gt,mcrab,efilter,mDoETOnly]

#---MC--------------------------------------------------------------------------------------------------------------------------
#>>>>>      #MC RunIISummer20UL18RECO
            #gt = 'globalTag=106X_upgrade2018_realistic_v15_L1v1'
            #runinfo = 'runera=Run2_2018'
            #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter,runinfo]
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
    RunIIS = 'RunIISummer20UL18MiniAODv2'

    runDataset = [

        #['/MET/Run2018A-15Feb2022_UL2018-v1/MINIAOD',1],
        #['/MET/Run2018B-15Feb2022_UL2018-v1/MINIAOD',1],
        #['/MET/Run2018C-15Feb2022_UL2018-v1/MINIAOD',1],
        #['/MET/Run2018D-15Feb2022_UL2018-v1/MINIAOD',1],

        #['/MET/Run2017B-09Aug2019_UL2017_rsb-v1/MINIAOD',1],
        #['/MET/Run2017C-09Aug2019_UL2017_rsb-v1/MINIAOD',1],
        #['/MET/Run2017D-09Aug2019_UL2017_rsb-v1/MINIAOD',1],
        #['/MET/Run2017E-09Aug2019_UL2017_rsb-v1/MINIAOD',1],
        #['/MET/Run2017F-09Aug2019_UL2017_rsb-v1/MINIAOD',1],

        ##['/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/MINIAOD',1],
        ##['/MET/Run2016F-21Feb2020_UL2016-v1/MINIAOD',1],
        ##['/MET/Run2016G-21Feb2020_UL2016-v1/MINIAOD',1],
        ##['/MET/Run2016H-21Feb2020_UL2016-v2/MINIAOD',1],

        ['/MET/Run2016B-ver1_HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016B-ver2_HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016C-HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016D-HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016E-HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016F-HIPM_UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016F-UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016G-UL2016_MiniAODv2-v2/MINIAOD',1],
        ['/MET/Run2016H-UL2016_MiniAODv2-v2/MINIAOD',1],

        #["/JetMET/Run2022C-19Dec2023-v1/MINIAOD",1],
        #["/JetMET/Run2022D-19Dec2023-v1/MINIAOD",1],
        #["/JetMET/Run2022E-19Dec2023-v1/MINIAOD",1],
        #["/JetMET/Run2022F-19Dec2023-v2/MINIAOD",1],
        #["/JetMET/Run2022G-19Dec2023-v1/MINIAOD",1],

        #['/GJets_HT-100To200_'+TuneCP5MP+'/'+RunIIS+'-4cores5k_106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',1],
        #['/GJets_HT-200To400_'+TuneCP5MP+'/'+RunIIS+'-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',1],
        #['/GJets_HT-400To600_'+TuneCP5MP+'/'+RunIIS+'-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',1],
        #['/GJets_HT-40To100_'+TuneCP5MP+'/'+RunIIS+'-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',1],
        #['/GJets_HT-600ToInf_'+TuneCP5MP+'/'+RunIIS+'-106X_upgrade2018_realistic_v16_L1v1-v2/MINIAODSIM',1],

    ] 

    for dataset in runDataset :
        docrab( dataset )

run_multi()
