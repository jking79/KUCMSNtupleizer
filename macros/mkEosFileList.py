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

mspc = '/store/user/jaking/'
mdis = '/store/user/lpcsusylep/jaking/'
eosll = 'eos root://cmseos.fnal.gov ls '
#command = eosll+mspc+'LLPGamma/llpga_GMSB_AOD_v48/'
#command = eosll+mspc+'A/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v58/'
#command = eosll+mdis+'/ecalTiming/tt_KUCCRes_126_Test/EGamma/'
#command = eosll+mdis+'LLPGamma/llpga_GJets_AOD_v57/'
#command = eosll+mdis+'LLPGamma/llpga_GMSB_AOD_v59/'
#command = eosll+mspc+'ecalTiming/'
#command = eosll+mdis+'ecalTiming/EGamma/'
#command = eosll+mspc+'EGamma/'
#command = eosll+mdis+'KUCMSNtuple/GMSB_AOD_v1/'
command = eosll+mdis+'/KUCMSNtuple/kucmsntuple_GJETS_AOD_v5/'
#command = eosll+mdis+'/KUCMSNtuple/WJETS_AOD_v1/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_AOD_v2/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_AOD_v5/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DYTT_AOD_v2/'

version = ''
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'
rootfile = '.root'
#dirselect = 'HTo2LongLivedTo4b'
#dirselect = '_newRtParams4_v26b_'
#dirselect = '_newRtParams3_test_v26_'
#dirselect = 'tt_kurhs_124cc5_cert'
#dirselect = '22eraC_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = '22eraC_CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'
#dirselect = 'CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'

#dirselect = 'GMSB_L-400'
dirselect = 'GJets_HT'
#dirselect = 'WJetsToLNu_HT-800'
#dirselect = 'QCD_HT100to200'
#dirselect = 'GMSB_L-400TeV'
#dirselect = 'DYJetsToLL_M-50'
#dirselect = 'TTJets'

#dirselect = ''

debug = True
#debug = False

#deep = True
deep = False

targdirs = []
subdirlist1 = []
subdirlist2 = []
subdirlist3 = []
filelist = []
theFileList = ''

dirls = bashout( command ).splitlines()
if debug : print( '-------------------------------------------------')
for line in dirls:
	#print( line )
	if dirselect in line : targdirs.append( line )
    #targdirs.append( line )

if debug : print( targdirs )

if deep :
	for mydir in targdirs:
		command1 = command+mydir+'/'
		subdir1 = bashout( command1 ).rstrip().splitlines()
		#print( subdir1 )
		#print( mydir+'/'+subdir1+'/' )
		for line in subdir1 : 
			#print( line )
			if version in line : 
				subdirlist1.append( mydir+'/'+line+'/' )
		#print( subdirlist1 )
else : 
	for mydir in targdirs:
		subdirlist1.append( mydir+'/' )

if debug : print( subdirlist1 )
for thesubdir in subdirlist1 :
	command2 = command+thesubdir+'/'
	subdir2 = bashout( command2 ).rstrip().splitlines()
	#print( thesubdir+subdir2+'/0000/' )
	for subdir in subdir2 : 
		command3 = command+thesubdir+subdir+'/'
		subdir3 = bashout( command3 ).rstrip().splitlines()
		for subsubdir in subdir3 : 
			subdirlist2.append(thesubdir+subdir+'/'+subsubdir)

for thesubdir2 in subdirlist2 :
    command4 = command+thesubdir2+'/'
    subdir4 = bashout( command4 ).rstrip().splitlines()
    #print( thesubdir+subdir2+'/0000/' )
    for subdir in subdir4 :
        subdirlist3.append(thesubdir2+'/'+subdir+'/')
        #command5 = command+thesubdir+subdir+'/'
        #subdir5 = bashout( command5 ).rstrip().splitlines()
        #for subsubdir in subdir5 :
            #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')


if debug : print( subdirlist3 )
for subdir2 in subdirlist3:
	lists = bashout( command+subdir2 ).rstrip().splitlines()
	for line in lists :
		if rootfile in line : filelist.append(subdir2+line)

for thefile in filelist:
	print( thefile )

	#filename = 'tmp_'+subdir2.split('/')[1]+'.root '
	#print( filename )
	#lists = bashout( "eosls "+mspc+"LLPGamma/"+subdir2 ).rstrip()
	#print( subdir2 )
	#haddcommand = "hadd -f "+filename+"`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+" | grep '\.root'`"
	#haddcommand = "`xrdfs root://cmseos.fnal.gov ls -u "+mspc+"LLPGamma/"+subdir2+"`"
	#haddcommand = "hadd "+filename+lists
	#print( mspc+"LLPGamma/"+subdir2 )
	#print( haddcommand )
	#doCommand( haddcommand )
	#print( '---------------------------------------------------' )

	
#print( bashout( 'hadd llpgana_HTo2LongLivedTo4b_t37MC_noele_005_jetht_emf00bc3rh2e_id2pt200nrh5eta15rhe2.root ' + theFileList ) )	