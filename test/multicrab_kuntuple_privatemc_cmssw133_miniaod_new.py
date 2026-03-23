#! /usr/bin/env python

import os
from optparse import OptionParser
import argparse
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
from tools.sample_lookup import SMS_MINI
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


def docrab( dataset, options ):

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
    
    config.General.workArea    = options['workArea']
    config.General.requestName = None
    
    config.JobType.pluginName  = 'Analysis'
    config.JobType.psetName    = 'llpgana_miniaod.py'
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
    
    config.Data.userInputFiles = open(dataset).readlines()
    
    print( '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>' )
    print( 'Input dataset for Crab Job : ' )
    #print( inDO )
    # inDO[0] is of the form /A/B/C. Since A+B is unique for each inDS, use this in the CRAB request name.
    print("dataset og",dataset)
    dataset = dataset.split("/")[1]
    primaryDataset = dataset[dataset.find('filelist_')+9:]
    primaryDataset = primaryDataset[:primaryDataset.find("_")]
    print( "primaryDataset",primaryDataset )
    runEra         = dataset[dataset.find(primaryDataset)+len(primaryDataset)+1:dataset.find("_MINI")] #(inDO[0].split('/')[2]) # justin
    if "N2ctau" in runEra:
        runEra = runEra.replace("N2ctau","ctau")
        #assuming sample name also is in wrong place
        sample = runEra[:runEra.rfind("_")]
        sample = sample[sample.rfind("_")+1:]
        print("sample",sample)
        primaryDataset += "_"+sample
        runEra = runEra.replace("_"+sample,"")
    print( "runEra",runEra )
    sample_type        = "FULLMINI"#inDO[0].split('/')[3]
    print( "dataset",sample_type )
    
    print( 'with : ',dataset )
    
    #trial          = "kucmsntuple_gogoG_Sig_IPM100_v26" #
    #trial          = "kucmsntuple_gogoZ_Sig_IPM100_v27" # 
    #trial          = "kucmsntuple_gogoZ_Sig_SVIPM100_v27" # filter on #SVs( lep or had ) - redefinded loose electron
    #trial          = "kucmsntuple_sqsqG_Sig_SVIPM100_v27" #
    #trial          = "kucmsntuple_SMS_Sig_SVIPM100_v31" #
    #trial          = "kucmsntuple_MET_R17_SVIPM100_p9_v31"# 
    trial          = "kucmsntuple_SMS_"+options["outfilter"] #
    
    # set trial name - used in output path ?
    
    print( 'processing for : ',trial )
    
    # !!!!!!!!!!!!!!!!!    Set output paths !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!1
    config.Data.outLFNDirBase  = "/store/group/lpcsusylep/jaking/KUCMSNtuple/"+trial+"/"
    config.General.requestName   = trial+"_"+primaryDataset+"_"+sample_type+"_"+runEra+"_request"
    config.Data.outputDatasetTag = trial+"_"+primaryDataset+"_"+sample_type+"_"+runEra
    
    print("writing to",config.Data.outputDatasetTag)#,"trial",trial,"primaryDataset",primaryDataset,"dataset",dataset,"runEra",runEra)
    
    #fSVHPM100 = 'eventFilter=SVHPMet100'
    #fSVIPM100 = 'eventFilter=SVIPMet100'
    #fM100 = 'eventFilter=MET100'
    #fAL1P = 'eventFilter=AL1IsoPho'
    #fIPM100 = 'eventFilter=IsoPhoMet100'
    #fNone = 'eventFilter=None'
    #fAL1E = 'eventFilter=AL1SelEle'
   
    nfilter = ""
    if "_" in options['outfilter']:
        nfilter = options["outfilter"].split("_")[0]
    else:
        nfilter = options["outfilter"]


    geninfo = 'hasGenInfo=True'
    mcrab = 'multicrab=True'
    
    runinfo = 'runera=Run3'
    doMetFilters = 'runMETFilters=True'
    
    if options['doSV']:
        mNoSV = 'doSV=True'           
    else:
        mNoSV = 'doSV=False'           
    
    #efilter = fIPM100
    efilter = "eventFilter="+nfilter#fSVHPM100 #!!!!!!!!!!!!!!!!!!!!!11
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
    config.JobType.pyCfgParams   = [gt,mcrab,geninfo,efilter,runinfo,mNoSV]
    print("params",config.JobType.pyCfgParams)
#-----------------------------------------------------------------------------------------------------------------------------

    # Submit.
    try:
        print( "Submitting for input dataset %s" % primaryDataset + '_' + runEra + '_' + sample_type )
        if(options['dryRun']):
            exit()
        crabCommand("submit", config = config, *options.crabCmdOpts.split())
        os.system("rm -rf %s/crab_%s/inputs" % (config.General.workArea, config.General.requestName))
    #except HTTPException as hte:
    #    print( "Submission for input dataset %s failed: %s" % (inDO[0], hte.headers) )
    except ClientException as cle:
        print( "Submission for input dataset %s failed: %s" % (inDO[0], cle) )


##33333333333333333333333333333333333333333333333333333333333

def run_multi():
    parser = argparse.ArgumentParser()
    parser.add_argument('-inputSample','-i',help='SMS process',required=True,choices=['gogoGZ'])
    parser.add_argument('--mGl',help='parent sparticle mass',default='')
    parser.add_argument('--mN2',help='n2 mass',default='')
    parser.add_argument('--mN1',help='lsp mass',default='')
    parser.add_argument('--ct',help='ctau',default='')
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
    parser.add_argument('--filter',help='specify event filter',choices=['SVHPM100'],default='SVHPM100')
    parser.add_argument('--unitsPerJob',help='number of files to run per job (default = 1)',default=1)
    parser.add_argument('--dryRun',help='will show info about submissions (name of output path, etc) but won\'t submit',action='store_true',default=False)
    args = parser.parse_args()
    outfilter = args.filter
    if args.output != "":
        outfilter += "_"+args.output
    options = {'workArea':args.workArea, 'crabCmdOpts':args.crabCmdOpts, 'outfilter':outfilter, 'unitsPerJob':args.unitsPerJob, 'dryRun':args.dryRun,'doSV': not args.noSV}

    
    for dataset in SMS_MINI :
        print("dataset",dataset)
        if args.inputSample not in dataset:
            continue
        if args.mGl != "" and "mGl-"+args.mGl not in dataset:
            continue
        if args.mN2 != "" and "mN2-"+args.mN2 not in dataset:
            continue
        if args.mN1 != "" and "mN1-"+args.mN1 not in dataset:
            continue
        if args.ct != "" and "ctau-"+args.ct not in dataset:
            continue
        print("doing point",dataset)
        docrab( dataset, options )

run_multi()
