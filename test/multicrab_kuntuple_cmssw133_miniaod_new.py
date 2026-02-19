#! /usr/bin/env python

import os
from optparse import OptionParser
import argparse
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from sample_lookup import MINIAOD_SAMPLES, MC_PROCS, DATA_NAMES, JSON_DICT, GLOBALTAGS_DICT
#from httplib import HTTPException

def docrab( dataset, options ):
    print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
    # The submit command needs special treatment.
    inputJSON = ''
    print('dataset',dataset)
    year = options['year']
    if "SIM" not in dataset:
        if year not in JSON_DICT.keys():
            print("Error: run JSON not provided for year",year)
            exit()
        inputJSON = JSON_DICT[year]
    
    #--------------------------------------------------------
    # This is the base config:
    #--------------------------------------------------------
    from CRABClient.UserUtilities import config
    config = config()
    
    config.General.workArea    = options['workArea']
    config.General.requestName = None
    
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'llpgana_miniaod.py'
    config.JobType.pyCfgParams = None
    
    config.section_("Site")
    #config.Site.whitelist = ['T1_*','T2_US_*','T2_IT_*','T2_DE_*','T2_ES_*','T2_FR_*','T2_UK_*']
    config.Site.blacklist = ['T2_CH_CSCS','T2_ES_CIEMAT']
    
    config.Data.inputDataset   = dataset
    if "SIM" not in dataset and not options['noMaskLumi']: #applying after ntuple stage to check for correct # of event processed
        print("Applying lumi mask",inputJSON)
        config.Data.lumiMask       = inputJSON
    config.Data.splitting = 'FileBased' # Justin
    config.Data.unitsPerJob = options['unitsPerJob'] 
    #config.JobType.allowUndistributedCMSSW = True
    #config.JobType.maxJobRuntimeMin = 1800
    #config.JobType.maxJobRuntimeMin = 2280
    #config.JobType.maxMemoryMB = 3000
    config.Data.publication    = False
    config.Site.storageSite    = 'T3_US_FNALLPC'
    #--------------------------------------------------------
    
    # Will submit one task for each of these input datasets.
    print( 'Input dataset for Crab Job : ' )
    print(dataset)
    dataset1 = dataset
    dataset = dataset[1:]
    idx = -1
    if "SIM" in dataset:
        idx = dataset.find("_")
        primaryDataset = dataset[:idx]
        dataset_tier = "MINIAODSIM"
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
        #print( primaryDataset )
        test = primaryDataset
        test = (dataset1.split('/')[1]).split('_13T')[0]
        #print( test )
        runEra         = (dataset1.split('/')[2]).split('RECO')[0]
        #print( runEra )
        dataset        = dataset1.split('/')[3]
        #print( dataset )
        request = test+"_"+dataset+"_"+runEra
        #print( "primaryDataset", primaryDataset )
        #print( "runEra", runEra )
        trial = "kucmsntuple_"+primaryDataset[:primaryDataset.find("-")]
    print("opts",options)
    #output tag
    trialtag = "R"+year[-2:]+"_"+options['filter']+"_"+options['output']
    print("ntuple tag is",trialtag)
    if trialtag != "" and trialtag[0] != "_":
       trial += "_"+trialtag
    print( 'processing for : ',trial )
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
    
    
    evt_filter = "eventFilter="+options['filter']
    print("evt_filter",evt_filter)
    gt = ""
    #setting global tag
    #MC
    print("dataset1",dataset1)
    if "SIM" in dataset1:
        if "RunIISummer20UL18RECO" in dataset1:
            gt = 'globalTag=106X_upgrade2018_realistic_v11_L1v1'
        elif "RunIIAutumn18DRPremix" in dataset1:
            gt = 'globalTag=94X_mc2017_realistic_v11'
        else:
            gt = ''
        geninfo = 'hasGenInfo=True'
    else: #data - TODO: check dataset1 match strings
        gt = "globalTag="+GLOBALTAGS_DICT[year]
        geninfo = 'hasGenInfo=False'
    print("globalTag",gt)
    params = [gt, "multicrab=True", geninfo,'doSV='+str(options['doSV']),evt_filter,'runera=Run2'] #turn off SV collection by default
    config.JobType.pyCfgParams   = params

#---MC----------------------------------------------------------------------------------------------------------------------
#>>>>>  #MC RunIISummer20UL18RECO
        #gt = 'globalTag=106X_upgrade2018_realistic_v11_L6v1'
        ##config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True']
        #config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter]
        ##config.JobType.pyCfgParams   = ['globalTag=106X_upgrade2018_realistic_v11_L1v1','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>  #MC GMSB RunIIFall17DRPremix  #globalTag=94X_mc2017_realistic_v14  #  <<< comment/uncomment lumi mask when using/!using MC
        ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True']
        ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v14','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>>>>>>>>>>>     #MC RunIIAutumn18DRPremix 102X_upgrade2018_realistic_v15  
        ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True']
        ##config.JobType.pyCfgParams   = ['globalTag=94X_mc2017_realistic_v11','multicrab=True','hasGenInfo=True','eventSkim=AL1IsoPho']
#>>>>>  #MC Justin
        #gt = 'globalTag=94X_mc2017_realistic_v14'
        ##config.JobType.pyCfgParams   = [gt,mcrab,geninfo,skimNone]
#-----------------------------------------------------------------------------------------------------------------------------

    print("params",params)
    # Submit.
    print("trial",trial)
    if(options['dryRun']):
        exit()
    try:
        print( "Submitting for input dataset %s" % primaryDataset + '_' + runEra + '_' + dataset )
        crabCommand('submit', config = config, *options['crabCmdOpts'].split())
        os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
    #except HTTPException as hte:
    #    print( "Submission for input dataset %s failed: %s" % (inDO[0], hte.headers) )
    except ClientException as cle:
        print( "Submission for input dataset %s failed: %s" % (trial, cle) )


#############################################################

def run_multi():
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputSample','-i',help='datasets for PDs/MC',required=True,choices=['DisplacedJet','EGamma','DoubleEG','GJets','QCD','MET','JetMET','JetMET0','JetMET1','JetHT'])

    parser.add_argument('--HT',help='ht bin for MCs',default='')
    parser.add_argument('--era',help='era (run) for PDs, default is set to \'mc\'',default='mc')
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
    parser.add_argument('--noSV',help='turn off SV collection (default = False)',action='store_true')
    parser.add_argument('--filter',help='specify event filter',choices=['50MET75','InvMET100','Pho30','InvMetPho30','MET100','AL1NpSC'],required=True)
    parser.add_argument('--noMaskLumi',help='do not apply lumi mask (default = off)',default=False,action='store_true')
    parser.add_argument('--unitsPerJob',help='number of files to run per job (default = 1)',default=1)
    parser.add_argument('--dryRun',help='will show info about submissions (name of output path, etc) but won\'t submit',action='store_true',default=False)
    args = parser.parse_args()
    if args.noMaskLumi:
        args.output = "nolumimask_"+args.output
    if args.noSV:
        args.output = "noSV_"+args.output
    options = {"workArea":args.workArea, "crabCmdOpts":args.crabCmdOpts, "filter":args.filter,"noMaskLumi":args.noMaskLumi,"unitsPerJob": args.unitsPerJob,"dryRun":args.dryRun,"doSV": not args.noSV, "year":args.year,"output":args.output}

    if args.inputSample in MC_PROCS:
        sample_type = "mc"
        data = False
    elif args.inputSample in DATA_NAMES:
        sample_type = "data"
        data = True
    sample_key = f"{sample_type}_{args.year}"
    sample_list = MINIAOD_SAMPLES[sample_key]
    samples = []
    for sample in sample_list:
        if args.HT != '' and args.HT not in sample:
            continue
        if args.era != 'mc' and args.era not in sample:
            continue
        if args.inputSample not in sample:
            continue
        samples.append(sample)
    for dataset in samples:
        docrab( dataset, options )

    

run_multi()
