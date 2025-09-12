import subprocess
import sys
import os

# uses python3

def bash( bashCommand ):
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        #process = subprocess.Popen(bashCommand.split())
        output, error = process.communicate()
        return output ,error

def bashout( command ):
        #output = subprocess.check_output( command, shell=True)
        output = subprocess.run( command, shell=True, check=True, stdout=subprocess.PIPE, universal_newlines=True )
        return output.stdout

def doCommand( command ):
        output = os.system( command )
        return output

dl2016mc = [ 

    '/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-4cores5k_106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT50to100_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-2000to4000_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-200to500_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-4000to8000_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-500to1000_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-8000toInf_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80to200_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v2/AODSIM',
    '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v4/AODSIM',
    '/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v4/AODSIM',
    '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext2-v3/AODSIM',
    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v2/AODSIM',
    '/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v2/AODSIM',
    '/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/WGG_5f_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v2/AODSIM',

]

dl2016data = [ 

    '/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD', 
    '/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD', 
    '/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD', 
    '/MET/Run2016G-21Feb2020_UL2016-v1/AOD', 
    '/MET/Run2016H-21Feb2020_UL2016-v2/AOD', 

]




dl2017mc = [

    '/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-4cores5k_106X_mc2017_realistic_v6-v2/AODSIM',
    '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',
    '/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-2000to4000_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-200to500_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-4000to8000_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-500to1000_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-8000toInf_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80to200_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v2/AODSIM',
    '/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v2/AODSIM',
    '/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/WGG_5f_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v2/AODSIM',

]

dl2017data = [

    '/MET/Run2017A-09Aug2019_UL2017_rsb-v1/AOD', 
    '/MET/Run2017B-09Aug2019_UL2017_rsb-v1/AOD', 
    '/MET/Run2017C-09Aug2019_UL2017_rsb-v1/AOD', 
    '/MET/Run2017D-09Aug2019_UL2017_rsb-v1/AOD', 
    '/MET/Run2017E-09Aug2019_UL2017_rsb-v1/AOD',
    '/MET/Run2017F-09Aug2019_UL2017_rsb-v1/AOD', 

]


dl2023data = [ 

    '/JetMET0/Run2023B-19Dec2023-v1/AOD',		
    '/JetMET0/Run2023C-19Dec2023-v1/AOD',		
    '/JetMET0/Run2023D-19Dec2023-v1/AOD',		
    '/JetMET1/Run2023B-19Dec2023-v1/AOD',		
    '/JetMET1/Run2023C-19Dec2023-v1/AOD',		
    '/JetMET1/Run2023D-19Dec2023-v1/AOD',		

]

dl2024data = [

    '/JetMET0/Run2024A-PromptReco-v1/AOD',
    '/JetMET0/Run2024B-PromptReco-v1/AOD',
    '/JetMET0/Run2024C-PromptReco-v1/AOD',
    '/JetMET0/Run2024D-PromptReco-v1/AOD',
    '/JetMET0/Run2024E-PromptReco-v1/AOD',
    '/JetMET0/Run2024E-PromptReco-v2/AOD',
    '/JetMET0/Run2024F-PromptReco-v1/AOD',
    '/JetMET0/Run2024G-PromptReco-v1/AOD',
    '/JetMET0/Run2024H-PromptReco-v1/AOD',
    '/JetMET0/Run2024I-PromptReco-v1/AOD',
    '/JetMET0/Run2024I-PromptReco-v2/AOD',
    '/JetMET0/Run2024J-PromptReco-v1/AOD',
    '/JetMET1/Run2024A-PromptReco-v1/AOD',
    '/JetMET1/Run2024B-PromptReco-v1/AOD',
    '/JetMET1/Run2024C-PromptReco-v1/AOD',
    '/JetMET1/Run2024D-PromptReco-v1/AOD',
    '/JetMET1/Run2024E-PromptReco-v1/AOD',
    '/JetMET1/Run2024E-PromptReco-v2/AOD',
    '/JetMET1/Run2024F-PromptReco-v1/AOD',
    '/JetMET1/Run2024G-PromptReco-v1/AOD',
    '/JetMET1/Run2024H-PromptReco-v1/AOD',
    '/JetMET1/Run2024I-PromptReco-v1/AOD',
    '/JetMET1/Run2024I-PromptReco-v2/AOD',
    '/JetMET1/Run2024J-PromptReco-v1/AOD',

]

dl2022ExoData = [ 

    '/JetMET/Run2022C-EXODelayedJetMET-27Jun2023-v1/AOD',
    '/JetMET/Run2022D-EXODelayedJetMET-27Jun2023-v2/AOD',
    '/JetMET/Run2022E-EXODelayedJetMET-27Jun2023-v1/AOD',

]

dl2023ExoData = [

    '/JetMET/Run2023E-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2023A-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET0/Run2023B-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2023C-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2023C-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET0/Run2023C-EXODelayedJetMET-PromptReco-v3/AOD',
    '/JetMET0/Run2023C-EXODelayedJetMET-PromptReco-v4/AOD',
    '/JetMET0/Run2023D-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2023D-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET1/Run2023A-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET1/Run2023B-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2023C-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2023C-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET1/Run2023C-EXODelayedJetMET-PromptReco-v3/AOD',
    '/JetMET1/Run2023C-EXODelayedJetMET-PromptReco-v4/AOD',
    '/JetMET1/Run2023D-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2023D-EXODelayedJetMET-PromptReco-v2/AOD',

]

dl2024ExoData = [

    '/JetMET0/Run2024A-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024B-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024C-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024D-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024E-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024E-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET0/Run2024F-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024G-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024H-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024I-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET0/Run2024I-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET0/Run2024J-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024A-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024B-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024C-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024D-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024E-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024E-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET1/Run2024F-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024G-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024H-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024I-EXODelayedJetMET-PromptReco-v1/AOD',
    '/JetMET1/Run2024I-EXODelayedJetMET-PromptReco-v2/AOD',
    '/JetMET1/Run2024J-EXODelayedJetMET-PromptReco-v1/AOD',

]

dlRun2SkimsData = [ 

    '/DisplacedJet/Run2018A-15Feb2022_UL2018-v1/AOD',
    '/DisplacedJet/Run2018B-15Feb2022_UL2018-v1/AOD',
    '/DisplacedJet/Run2018C-15Feb2022_UL2018-v1/AOD',
    '/DisplacedJet/Run2018D-15Feb2022_UL2018-v1/AOD',
    '/DisplacedJet/Run2017C-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017D-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017E-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017F-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017C-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017D-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017E-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2017F-09Aug2019_UL2017-v1/AOD',
    '/DisplacedJet/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016F-21Feb2020_UL2016-v1/AOD',
    '/DisplacedJet/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
    '/DisplacedJet/Run2016G-21Feb2020_UL2016-v1/AOD',
    '/DisplacedJet/Run2016H-21Feb2020_UL2016-v1/AOD',
]

dl2018ZJets = [

    '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

]

dl2018GJets = [ 

    '/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 
    '/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 
    '/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 
    '/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 
    '/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 

]

dl2018WJets = [

    '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',

]

dl2018QCD = [

    '/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

dl201DiPho = [

    '/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-2000to4000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-200to500_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-4000to8000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-500to1000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-8000toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80to200_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

dl2018DiPho = [

    '/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-2000to4000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-200to500_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-4000to8000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-500to1000_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-8000toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80to200_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

dlttwz = [

    '/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',
    '/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

]

dlwwwzz = [

    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    '/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    '/WGG_5f_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

dl22 = [ '/MET/Run2022C-27Jun2023-v2/AOD' ]

#query = 'dasgoclient --json -query=\'dataset='
query = 'dasgoclient --json -query=\'site dataset=' 

datalist = dl22

for data in datalist :
    
    command = query + data + '\''
    #output = bashout( command ).splitlines()
    #output = bashout( command ).split('size":')
    #output = bashout( command ).split('num_lumi":')
    #output = bashout( command ).split('files":')
    output = bashout( command ).split('\n')
    print( '--------------------------------------------------------------------------')
    print( data )
    for line in output :
        if 'block_completion":"100.00%' in line :
            print( line )

    #fir = output[1].split('size":')[0]
    #sec = output[1].split('size":')[1] 
    #.split('}],')[0]
    #size = int(output[1].split('}]')[0])/1000000000000
    #size = output[1].split('}]')[0]
    #print( data,  size )
    #print( data, '\n',  fir, '\n', sec )

