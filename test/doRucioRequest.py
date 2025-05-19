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
    'cms:/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
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

    'cms:/QCD_HT1000to1500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    ##'cms:/QCD_HT1500to2000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    'cms:/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT200to300_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT300to500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT500to700_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT50to100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/QCD_HT700to1000_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]

met_list = [

    #'cms:/MET/Run2017B-17Nov2017-v1/AOD', #4.2
    #'cms:/MET/Run2017C-17Nov2017-v1/AOD', #27.9
    #'cms:/MET/Run2017D-17Nov2017-v1/AOD', #5.3
    #'cms:/MET/Run2017E-17Nov2017-v1/AOD', #28.5
    #'cms:/MET/Run2017F-17Nov2017-v1/AOD', #75.3

]

met18UL = [

    #'cms:/MET/Run2018A-15Feb2022_UL2018-v1/AOD', #17
    #'cms:/MET/Run2018B-15Feb2022_UL2018-v1/AOD', #9
    #'cms:/MET/Run2018C-15Feb2022_UL2018-v1/AOD', #10
    #'cms:/MET/Run2018D-15Feb2022_UL2018-v1/AOD', #52

]


deg_list = [

    #'cms:/EGamma/Run2018A-15Feb2022_UL2018-v1/MINIAOD',
    'cms:/EGamma/Run2018B-15Feb2022_UL2018-v1/MINIAOD',
    #'cms:/EGamma/Run2018C-15Feb2022_UL2018-v1/MINIAOD',
    #'cms:/EGamma/Run2018D-15Feb2022_UL2018-v1/MINIAOD',

]

WjetsToLNU = [

    #'cms:/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    #'cms:/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
]

TandJets = [

    #'cms://TTJets_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    #'cms:/TGJets_TuneCP5_13TeV-amcatnlo-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',
    #'#cms:/ttZJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    #'cms:/ttWJets_TuneCP5_13TeV_madgraphMLM_pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
]

jetht18 = [ 'cms:/JetHT/Run2018C-15Feb2022_UL2018-v1/AOD' ]

dy_list = [ 'cms:/DY1JetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17DRPremix-PU2017_94X_mc2017_realistic_v11-v1/AODSIM' ]

#runlist = dy_list
runlist = deg_list
#runlist = qcd_list
#runlist = gjets_list
#runlist = gmsb_list
#runlist = met_list
#runlist = jetht18
#runlist = met18UL
#runlist = WjetsToLNU
#runlist = TandJets

for thefile in runlist :

    #thecommand = 'rucio add-rule --ask-approval --lifetime 5184000 ' + thefile + ' 1 T3_US_FNALLPC'
    thecommand = 'rucio add-rule --ask-approval --lifetime 2592000 ' + thefile + ' 1 T3_US_FNALLPC'
    print( '---------------------------------------------------' )
    print( thecommand )
    doCommand( thecommand )
    print( '---------------------------------------------------' )

#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
