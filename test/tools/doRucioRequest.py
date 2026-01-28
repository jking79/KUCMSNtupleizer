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

gjets_list = [
    #'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-4cores5k_106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

    ##'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    ###'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_4cores5k_94X_mc2017_realistic_v11-v3/AODSIM',
    ##'cms:/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-94X_mc2017_realistic_v10-v2/AODSIM',
    ##'cms:/GJets_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',
    ##'cms:/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v2/AODSIM',
    ##'cms:/GJets_HT-600ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM',

]

qcd_list = [

    #'cms:/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    ##'cms:/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

    #'cms:/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

met_list = [

    #'cms:/MET/Run2017B-17Nov2017-v1/AOD', #4.2
    #'cms:/MET/Run2017C-17Nov2017-v1/AOD', #27.9
    #'cms:/MET/Run2017D-17Nov2017-v1/AOD', #5.3
    #'cms:/MET/Run2017E-17Nov2017-v1/AOD', #28.5
    #'cms:/MET/Run2017F-17Nov2017-v1/AOD', #75.3

]

met18UL = [

    'cms:/MET/Run2018A-15Feb2022_UL2018-v1/AOD', #17
    'cms:/MET/Run2018B-15Feb2022_UL2018-v1/AOD', #9
    'cms:/MET/Run2018C-15Feb2022_UL2018-v1/AOD', #10
    'cms:/MET/Run2018D-15Feb2022_UL2018-v1/AOD', #52

]

met17UL = [

    #'cms:/MET/Run2017A-09Aug2019_UL2017_rsb-v1/AOD',
    #'cms:/MET/Run2017B-09Aug2019_UL2017_rsb-v1/AOD',
    ##'cms:/MET/Run2017C-09Aug2019_UL2017_rsb-v1/AOD',
    ##'cms:/MET/Run2017D-09Aug2019_UL2017_rsb-v1/AOD',
    #'cms:/MET/Run2017E-09Aug2019_UL2017_rsb-v1/AOD',
    'cms:/MET/Run2017F-09Aug2019_UL2017_rsb-v1/AOD',

]

deg_list = [

    #'cms:/EGamma/Run2018A-15Feb2022_UL2018-v1/MINIAOD',
    #'cms:/EGamma/Run2018B-15Feb2022_UL2018-v1/MINIAOD',
    #'cms:/EGamma/Run2018C-15Feb2022_UL2018-v1/MINIAOD',
    #'cms:/EGamma/Run2018D-15Feb2022_UL2018-v1/MINIAOD',

]

WjetsToLNU = [

    #'cms:/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',

    'cms:/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',

    #'cms:/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',

]

TandJets = [

    'cms:/TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',
    #'#cms:/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
]

ZjetsToNuNU = [

    #'cms:/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    'cms:/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

]

disjet18 = [

    'cms:/DisplacedJet/Run2018A-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018B-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018C-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018D-15Feb2022_UL2018-v1/AOD',

]

jetmetexo22 = [

    'cms:/JetMET/Run2022C-EXODelayedJetMET-27Jun2023-v1/AOD',
    'cms:/JetMET/Run2022D-EXODelayedJetMET-27Jun2023-v2/AOD',
    'cms:/JetMET/Run2022E-EXODelayedJetMET-27Jun2023-v1/AOD',

]

dipho = [ 
    
    #'cms:/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM', 
    'cms:/DiPhotonJetsBox_MGG-0to40_13TeV-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    'cms:/DiPhotonJetsBox_M40_80-sherpa/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

bosun = [

    #'cms:/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/WWG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/WWW_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    #'cms:/WWZ_4F_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    #'cms:/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/WZG_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/WZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
    #'cms:/WGG_5f_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/ZZZ_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v2/AODSIM',
]

jetht18 = [ 'cms:/JetHT/Run2018C-15Feb2022_UL2018-v1/AOD' ]

dy_list = [ 'cms:/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM' ]

egamma = [ 'cms:/EGamma/Run2018C-15Feb2022_UL2018-v1/AOD' ]

egamma16cali = [

    'cms:/DoubleEG/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/MINIAOD',
    'cms:/DoubleEG/Run2016C-21Feb2020_UL2016_HIPM-v1/MINIAOD',
    'cms:/DoubleEG/Run2016D-21Feb2020_UL2016_HIPM-v1/MINIAOD',
    'cms:/DoubleEG/Run2016E-21Feb2020_UL2016_HIPM-v1/MINIAOD',
    'cms:/DoubleEG/Run2016F-21Feb2020_UL2016-v1/MINIAOD',
    #'cms:/DoubleEG/Run2016G-21Feb2020_UL2016-v1/MINIAOD',
    #'cms:/DoubleEG/Run2016H-21Feb2020_UL2016-v1/MINIAOD',

]

jmet_23 = [

    #'cms:/JetMET0/Run2023B-EXODelayedJetMET-PromptReco-v1/AOD',
    'cms:/JetMET0/Run2023C-EXODelayedJetMET-PromptReco-v4/AOD',
    #'cms:/JetMET0/Run2023D-EXODelayedJetMET-PromptReco-v2/AOD',

    #'cms:/JetMET1/Run2023B-EXODelayedJetMET-PromptReco-v1/AOD',
    'cms:/JetMET1/Run2023C-EXODelayedJetMET-PromptReco-v4/AOD',
    'cms:/JetMET1/Run2023D-EXODelayedJetMET-PromptReco-v2/AOD',

    #'cms:/JetMET0/Run2023B-19Dec2023-v1/AOD',
    #'cms:/JetMET0/Run2023C-19Dec2023-v1/AOD',
    #'cms:/JetMET0/Run2023D-19Dec2023-v1/AOD',
    #'cms:/JetMET1/Run2023B-19Dec2023-v1/AOD',
    #'cms:/JetMET1/Run2023C-19Dec2023-v1/AOD',
    #'cms:/JetMET1/Run2023D-19Dec2023-v1/AOD',

]

#runlist = qcd_list
#runlist = gjets_list
#runlist = met_list
#runlist = met18UL
#runlist = WjetsToLNU
#runlist = TandJets
#runlist = ZjetsToNuNU
#runlist = jetmetexo22
#runlist = dipho
#runlist = bosun
#runlist = egamma16cali
#runlist = met17UL
runlist = jmet_23

for thefile in runlist :

    thecommand = 'rucio add-rule --ask-approval --lifetime 5184000 ' + thefile + ' 1 T3_US_FNALLPC'
    #thecommand = 'rucio rule add --ask-approval --lifetime=5184000 --copies 1 --rses T3_US_FNALLPC --did ' + thefile
    print( '---------------------------------------------------' )
    print( thecommand )
    doCommand( thecommand )
    print( '---------------------------------------------------' )

#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
