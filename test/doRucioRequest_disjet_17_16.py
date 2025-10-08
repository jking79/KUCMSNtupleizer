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


disjet18 = [

    'cms:/DisplacedJet/Run2018A-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018B-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018C-15Feb2022_UL2018-v1/AOD',
    'cms:/DisplacedJet/Run2018D-15Feb2022_UL2018-v1/AOD',

]

disjet1716 = [

    'cms:/DisplacedJet/Run2017C-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017D-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017E-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017F-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017C-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017D-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017E-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2017F-09Aug2019_UL2017-v1/AOD',
    'cms:/DisplacedJet/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    'cms:/DisplacedJet/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    'cms:/DisplacedJet/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    'cms:/DisplacedJet/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    'cms:/DisplacedJet/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
    'cms:/DisplacedJet/Run2016G-21Feb2020_UL2016-v1/AOD',
    'cms:/DisplacedJet/Run2016H-21Feb2020_UL2016-v1/AOD',

]

runlist = disjet1716

for thefile in runlist :

    #thecommand = 'rucio add-rule --ask-approval --lifetime 5184000 ' + thefile + ' 1 T3_US_FNALLPC'
    thecommand = 'rucio add-rule --ask-approval --lifetime 2592000 ' + thefile + ' 1 T3_US_FNALLPC'
    print( '---------------------------------------------------' )
    print( thecommand )
    doCommand( thecommand )
    print( '---------------------------------------------------' )

#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	
