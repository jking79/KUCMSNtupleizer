import subprocess
import sys
import os

def bash( bashCommand ):
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
	#process = subprocess.Popen(bashCommand.split())
	output, error = process.communicate()
	return output ,error

def bashout( command ):
	output = subprocess.check_output( command, shell=True)
	return output	

def doCommand( command ):
	output = os.system( command )
	return output

gmsb_list = [
#    'cms:/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-100TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-150TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-200TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-250TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-300TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    'cms:/GMSB_L-350TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-350TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-400TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-0p1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-500TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-0p1cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-10000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-1000cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-10cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-1200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-400cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
#    'cms:/GMSB_L-600TeV_Ctau-800cm_TuneCP5_13TeV-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
]

gjets_list = [
    #'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

    #'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    ##'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_4cores5k_94X_mc2017_realistic_v11-v3/AODSIM',
    'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v2/AODSIM',
    #'cms:/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
    #'cms:/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',

]

qcd_list = [


    ##'cms:/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM', #16.1
    ##'cms:/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #38.5
    ##'cms:/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #25.5
    ##'cms:/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #26.7
    #'#cms:/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #26.8
    ##'cms:/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #23.0
    ##'cms:/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_new_pmx_94X_mc2017_realistic_v11-v1/AODSIM', #8.5
    ##'cms:/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #6.1
    ##'cms:/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v1/AODSIM', #3.2

    #'cms:/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/QCD_HT100to200_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
    #'cms:/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/QCD_HT300to500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/QCD_HT500to700_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
    #'cms:/QCD_HT700to1000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/QCD_HT1000to1500_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_new_pmx_94X_mc2017_realistic_v11-v1/AODSIM',
    #'cms:/QCD_HT1500to2000_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
    #'cms:/QCD_HT2000toInf_TuneCP5_13TeV-madgraph-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',

]

met_list = [

    #'cms:/MET/Run2017B-17Nov2017-v1/AOD', #4.2
    #'cms:/MET/Run2017C-17Nov2017-v1/AOD', #27.9
    #'cms:/MET/Run2017D-17Nov2017-v1/AOD', #5.3
    #'cms:/MET/Run2017E-17Nov2017-v1/AOD', #28.5
    #'cms:/MET/Run2017F-17Nov2017-v1/AOD', #75.3

]

deg_list = [

    'cms:/DoubleEG/Run2017D-17Nov2017-v1/AOD',
    #'cms:/DoubleEG/Run2017E-17Nov2017-v1/AOD', 
    #'cms:/DoubleEG/Run2017F-17Nov2017-v1/AOD',
]

jetht18 = [ 'cms:/JetHT/Run2018C-15Feb2022_UL2018-v1/AOD' ]

dy_list = [ 'cms:/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM' ]

#runlist = dy_list
#runlist = deg_list
#runlist = qcd_list
runlist = gjets_list
#runlist = gmsb_list
#runlist = met_list
#runlist = jetht18
for thefile in runlist :

	#thecommand = 'rucio add-rule --ask-approval --lifetime 5184000 ' + thefile + ' 1 T3_US_FNALLPC'
    thecommand = 'rucio add-rule --ask-approval --lifetime 2592000 ' + thefile + ' 1 T3_US_FNALLPC'
    print( '---------------------------------------------------' )
    print( thecommand )
    doCommand( thecommand )
    print( '---------------------------------------------------' )

#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
