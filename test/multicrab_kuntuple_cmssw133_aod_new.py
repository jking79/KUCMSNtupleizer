#! /usr/bin/env python

import os
from optparse import OptionParser
import argparse
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

    options = parser.parse_args()

    if not options.crabCmd:
        parser.error("(-c CMD, --crabCmd=CMD) option not provided.")
    if options.crabCmd != 'submit':
        if not options.workArea:
            parser.error("(-w WAR, --workArea=WAR) option not provided.")
        if not os.path.isdir(options.workArea):
            parser.error("'%s' is not a valid directory." % (options.workArea))

    return options


def docrab( dataset, options ):

    # External files needed by CRAB
    ##inputJSON    = 'golden2016.json'
    ##inputJSON    = 'golden2017.json'
    ##inputJSON    = 'Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
    #inputJSON    = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
    #inputJSON    = 'Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON.txt'
    #inputJSON    = 'Cert_314472-325175_13TeV_PromptReco_Collisions18_JSON.txt'
    inputJSON = ''
    print('dataset',dataset)
    year = ""
    if "17" in dataset:
        year = "2017"
    if "18" in dataset:
        year = "2018"
    elif "Run2022" in dataset:
        year = "2022"
    elif "Run2023" in dataset:
        year = "2023"
    if "SIM" not in dataset:
        if year == "2017":
            inputJSON    = 'Cert_294927-306462_13TeV_UL2017_Collisions17_GoldenJSON.txt'
        if year == "2018":
            inputJSON    = 'Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt'
        elif year == "2022":
            inputJSON    = 'Cert_Collisions2022_355100_362760_Golden.json'
        elif year == "2023":
            inputJSON    = 'Cert_Collisions2023_366442_370790_Golden.json'
        else:
            print("Error: run JSON not provided for year",year)
            exit()
    print("json",inputJSON)

    #--------------------------------------------------------
    # This is the base config:
    #--------------------------------------------------------
    from CRABClient.UserUtilities import config
    config = config()
    config.General.workArea    = options[0]
    config.General.requestName = None
    
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'llpgana_mc_aod.py'
    config.JobType.pyCfgParams = None
    
    config.Data.partialDataset = True
    config.Data.inputDataset   = None
    #config.Data.splitting    = 'LumiBased' # MC&Data Set unitsperjob correctly for dataset !!!!!!!!!!!!!!!!!!!!!!!!!!
    config.Data.splitting    = 'FileBased' # MC&Data Set unitsperjob correctly for dataset !!!!!!!!!!!!!!!!!!!!!!!!!!
    if "SIM" not in dataset and options[3]: #applying after ntuple stage to check for correct # of event processed
        print("Applying lumi mask",inputJSON)
        config.Data.lumiMask       = inputJSON    
    #config.Data.unitsPerJob  =  10
    #config.Data.unitsPerJob   = 50000 # data  !!!!!! lumimask ?
    #config.Data.unitsPerJob  =  1500 # MC GMSB
    #config.Data.unitsPerJob  =  10000 # MC GJet
    #config.Data.unitsPerJob = 10000
    config.Data.unitsPerJob = options[4]
    #if "QCD" in dataset:
    #    config.Data.unitsPerJob  =  15000 # MC QCD - WJetsToLNu
    #config.Data.unitsPerJob  =  50000 # MC TTJet
    #config.Data.unitsPerJob  =  750 # MC DiPhoBox (DPB)
    
    config.JobType.allowUndistributedCMSSW = True
    config.Data.publication    = False
    config.Site.storageSite    = 'T3_US_FNALLPC'
    #--------------------------------------------------------
    
    # Will submit one task for each of these input datasets.
    config.Data.inputDataset     = dataset
    #config.Data.userInputFiles = open("kuntuple_SMS-GlGl_mGl-2000_mN2-1950_mN1-1000_gam_N2ctau-0p1_AOD_v23.txt").readlines()
    
    print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
    print( 'Input dataset for Crab Job : ' )
    print(dataset)
    dataset1 = dataset
    dataset = dataset[1:]
    idx = -1
    if "SIM" in dataset:
        idx = dataset.find("_")
        primaryDataset = dataset[:idx]
        dataset_tier = "AODSIM"
        print("primaryDataset",primaryDataset)
        print("dataset",dataset)
        recotag = dataset.split("/")[1]
        recotag = recotag[:recotag.find("-")]
        request = dataset[:dataset.find("_13T")] + "_" + dataset_tier + "_" + recotag
        print("request",request)
        trial = "kucmsntuple_"+primaryDataset
    else:
        primaryDataset = dataset.replace("/","-")
        idx = primaryDataset.find("/")
        print( primaryDataset )
        test = primaryDataset
        test = (dataset1.split('/')[1]).split('_13T')[0]
        #test = (dataset1.split('/')[1]).split('_Tune')[0]
        #test = (dataset1.split('/')[1]).split('-sherpa')[0]
        print( test )
        #runEra         = (dataset1.split('/')[2]) # justin
        #runEra         = (dataset1.split('/')[2]).split('-')[0]
        runEra         = (dataset1.split('/')[2]).split('RECO')[0]
        print( runEra )
        dataset        = dataset1.split('/')[3]
        print( dataset )
        request = test+"_"+dataset+"_"+runEra
        print( "primaryDataset", primaryDataset )
        #runEra         = dataset[idx:]#(inDO[0].split('/')[2]).split('-')[0]
        #idx = runEra.find("-")
        #runEra = runEra[:idx]
        #print( "runEra", runEra )
        dataset_tier        = 'AOD'
        trial = "kucmsntuple_"+primaryDataset[:primaryDataset.find("-")]
    print( "tier", dataset_tier )
   
    print("opts",options)
    #output tag
    trialtag = "R"+year[-2:]+"_"+options[2]
    print("ntuple tag is",trialtag)
    if trialtag != "" and trialtag[0] != "_":
       trial += "_"+trialtag 
    print( 'processing for : ',trial )
    # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
    #config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
    config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
    print("outLFNdirbase",config.Data.outLFNDirBase)
    config.General.requestName   = trial+"_"+request+"_request"
    print("requestname",config.General.requestName,len(config.General.requestName))
    config.Data.outputDatasetTag = trial+"_"+request
    print("outputDatasetTag",config.Data.outputDatasetTag)
    #print("trial",trial)
    #print("primaryDataset",primaryDataset)
    #print("runEra",runEra)
    #print("dataset",dataset)

    evt_filter = ""
    if "M100" in trialtag and "skimIP" not in trialtag:
        evt_filter = 'MET100'
    if "METl100" in trialtag:
        evt_filter = '100MET'
    if "50M75" in trialtag:
        evt_filter = '50MET75'
    elif "AL1IsoPho" in trialtag:
        evt_filter = 'AL1IsoPho'
    elif "AL1SelEle" in trialtag:
        evt_filter = 'AL1SelEle'
    elif "skimIPM100" in trialtag:
        evt_filter = 'IsoPhoMet100'
    elif "AL1NpSC" in trialtag:
        evt_filter = "AL1NpSC"
    else:
        evt_filter = 'None'
    evt_filter = 'eventFilter='+evt_filter

    mcrab = 'multicrab=True'
    
    gt = ""
    #setting global tag
    #MC
    if "SIM" in dataset1:
        if "RunIISummer20UL18RECO" in dataset1:
            gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
        elif "RunIIAutumn18DRPremix" in dataset1:
            gt = 'globalTag=94X_mc2017_realistic_v11'
        else:
            gt = ''
        geninfo = 'hasGenInfo=True'
    else: #data - TODO: check dataset1 match strings
        if "2017UL" in dataset1:
            gt = 'globalTag=106X_dataRun2_v20'
        elif "2018UL" in dataset1:
            gt = 'globalTag=106X_dataRun2_v36'
        elif "2022" in dataset1:
            #TODO - update to Run3 tag when pCastorRecord issue is fixed (CMSSW)
    #        gt = 'globalTag=140X_dataRun3_v17'
            gt = 'globalTag=106X_dataRun2_v36'
        elif "2023" in dataset1:
            #TODO - update to Run3 tag when pCastorRecord issue is fixed (CMSSW)
    #        gt = 'globalTag=140X_dataRun3_v17'
            gt = 'globalTag=106X_dataRun2_v36'
        else:
            gt = ''
        geninfo = 'hasGenInfo=False'
    print("globalTag",gt)
    params = [gt, mcrab, geninfo,'doSV=False'] #turn off SV collection by default
    params.append(evt_filter)
    print("cfgParams",params)
    config.JobType.pyCfgParams = params
    #--  selsect PD/MC dependent paramters
    #----------------------------------------------------------------------------------------------------------------------
    
   # >>>>>>>>>>>>     #2017UL #globalTag=106X_dataRun2_v20
        #config.JobType.pyCfgParams   = ['globalTag=106X_dataRun2_v20','multicrab=True']
    #>>>>>>>>>>>>     #2018UL #globalTag=106X_dataRun2_v36
        #gt = 'globalTag=106X_dataRun2_v36'
        #config.JobType.pyCfgParams   = [gt,mcrab,skimIPM100]
        #config.JobType.pyCfgParams   = [gt,mcrab,skimM100]
        #config.JobType.pyCfgParams   = [gt,mcrab,skimAL1P]
   # >>>>>>>>>>>>        #2022 #globalTag=124X_dataRun3_v15
        #config.JobType.pyCfgParams   = ['globalTag=124X_dataRun3_v15','multicrab=True']
    
    #----------------------------------------------------------------------------------------------------------------------
    #MC RunIISummer20UL18RECO
    #gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
    #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True']
    #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,evt_filter]
    #config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
    #MC GMSB RunIIFall17DRPremix  #globalTag=94X_mc2017_realistic_v14  #  <<< comment/uncomment lumi mask when using/!using MC
    #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True']
    #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
    #>>>>>>>>>     #MC RunIIAutumn18DRPremix 102X_upgrade2018_realistic_v15  
    #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True']
    #config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
    #MC Justin
    #gt = 'globalTag=94X_mc2017_realistic_v14'
    #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,skimNone]
    #----------------------------------------------------------------------------------------------------------------------
    # Submit.
    print("trial",trial)
    #exit()
    try:
        #print( "Submitting for input dataset %s" % primaryDataset + '_' + runEra + '_' + dataset )
        print( "Submitting for input dataset %s" % trial )
        #ret = crabCommand('submit', config = config, dryrun = dryrunopt)
        ret = crabCommand('submit', config = config)
        print("crabCmd",ret)
        os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
    except ClientException as cle:
        print( "Submission for input dataset %s failed: %s" % (trial, cle) )
    

##33333333333333333333333333333333333333333333333333333333333

def run_multi():
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputSample','-i',help='datasets for PDs/MC',required=True,choices=['DisplacedJet','EGamma','GJets','QCD','MET','JetMET','JetMET0','JetMET1'])

    parser.add_argument('--HT',help='ht bin for MCs',default='100to200')
    parser.add_argument('--era',help='era (run) for PDs, default is set to \'MC\'',default='MC')
    parser.add_argument('--year',help='year',default='2018')
    
    parser.add_argument('-w', '--workArea',
                      dest = 'workArea',
                      default = 'crabSubmissions',
                      help = "work area directory (only if CMD != 'submit')",
                      metavar = 'WAD')

    parser.add_argument('--crabCmdOpts',
                      default = '',
                      help = "options for crab command CMD")
    parser.add_argument('-o','--output',help='output tag for ntuples',default='')
    parser.add_argument('--filter',help='specify event filter',choices=['50M75','METl100','M100','AL1IsoPho','AL1NpSC','AL1SelEle'],default='')
    parser.add_argument('--maskLumi',help='apply lumi mask (default = off)',default=False,action='store_true')
    parser.add_argument('--unitsPerJob',help='number of files to run per job',default=1)
    args = parser.parse_args()

    outfilter = args.output
    if args.output != "":
        if not args.maskLumi:
            outfilter = args.filter+"_nolumimask_"+outfilter
        else:
            outfilter = args.filter+"_"+outfilter
    else:
        if not args.maskLumi:
            outfilter = args.filter+"_nolumimask"
        else:
            outfilter = args.filter
    options = [args.workArea, args.crabCmdOpts, outfilter, args.maskLumi, args.unitsPerJob]


    Tune = 'TuneCP5_13TeV-madgraphMLM-pythia8'

    #dataset
    tune = "TuneCP5_13TeV"
    mc = "madgraphMLM-pythia8"
    reco = ''
    if(args.era == "MC"):
        dataset = args.inputSample+"_HT"
        if(args.inputSample == "GJets"):
            dataset += "-"
        dataset += args.HT+"_"+tune+"-"+mc
        #reco info
        if(args.year == "2018"):
           reco = "RunIISummer20UL18RECO-"
        else:
            reco = "RunII"
        if(args.inputSample == "QCD" or args.inputSample == "GJets"):
            if(args.inputSample == "GJets" and args.HT == "100To200"):
                reco += "4cores5k_"
            reco += "106X_upgrade2018_realistic_v11_L1v1-v2"
            reco += "/AODSIM"
        dataset = "/"+dataset+"/"+reco
    else:
        dataset = args.inputSample+"/Run"+args.year+args.era
        #if args.inputSample == "MET" or args.inputSample == "EGamma":
        if args.year == "2018":
            reco = "15Feb2022"
        elif args.year == "2022":
            reco = "27Jun2023"
        elif args.year == "2023":
            reco = "19Dec2023"
        ver = "v1"
        if args.year == "2022" and args.inputSample == "MET" and args.era != "A":
            ver = "v2"
        if args.inputSample == "JetMET" and args.year == "2022" and args.era == "D":
            ver = "v2"
        if args.inputSample == "JetMET" and args.year == "2022":
            dataset += "-EXODelayedJetMET"
        dataset += "-"+reco
        if(args.year == "2018"):
            dataset += "_UL"+args.year
        dataset += "-"+ver+"/AOD"
        dataset = "/"+dataset
        print("dataset",dataset)
    docrab( dataset, options )

run_multi()
