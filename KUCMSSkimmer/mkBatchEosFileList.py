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
justin = '/store/user/janguian/'
andreas = '/store/user/lpcsusylep/anazario/'
eosll = 'eos root://cmseos.fnal.gov ls '

#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_R18_IPM100_v24/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_WJets_R18_SVIPM100_v28/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_SMS_Sig_SVIPM100_v28/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJets_R18_IPM100_v30/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R18_AL1NpSC_DEOnly_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R22_AL1NpSC_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DiPJBox_R18_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_R18_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_TTXJets_R18_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R22_AL1NpSC_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DisJet_R18_SVIPMet100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R22_AL1NpSC_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_SMS_Sig_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJets_R18_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_ZJets_R18_SVIPM100_v31/'
#command = eosll+mdis+'KUCMSNtuple/gammares_ul18/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DisJet_R18_SVIPM100_v31/'
command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R18_SVIPM100_v31/'

#version = ''
#version = '_v11_'
#version = '_noOOTAmp_'
#version = '_wthOOTAmp_'
#version = 'SMS'
#version = 'WJetsToLNu'
#version = 'ZJetsToNuNu'
#version = 'kucmsntuple_gogoG_Sig2_IPM100_v24_SMS-GlGl_AODSIM_mGl-2000_mN2-1500_mN1-1'
#version = 'kucmsntuple_gogoZ'
version = ''

rootfile = '.root'

dirselect = ''
#dirselect = 'HTo2LongLivedTo4b'
#dirselect = '_newRtParams4_v26b_'
#dirselect = '_newRtParams3_test_v26_'
#dirselect = 'tt_kurhs_124cc5_cert'
#dirselect = '22eraC_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = '22eraC_CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357328-357331'
#dirselect = 'noOOTCC_kustc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'
#dirselect = 'CCstc0_EGamma_MINIAOD_Run2022C-PromptReco-v1_357101-357268'

#dirselect = 'CRAB_UserFiles'
#dirselect = 'WJetsToLNu'
#dirselect = 'QCD'
#dirselect = 'TTJets'
#dirselect = 'TTGJets'
#dirselect = 'ttZJets'
#dirselect = 'ttWJets'
#dirselect = 'GJets_HT-200To400'
#dirselect = 'MET'
#dirselect = 'SMS-SqSq'

debug = True
#debug = False

deep = True
#deep = False

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
done = False
for thesubdir2 in subdirlist2 :
    command4 = command+thesubdir2+'/'
    subdir4 = bashout( command4 ).rstrip().splitlines()
    #print( thesubdir+subdir2+'/0000/' )
    for subdir in subdir4 :
        if rootfile in subdir : 
            filelist.append(thesubdir2+'/'+subdir)
            done = True
        else : subdirlist3.append(thesubdir2+'/'+subdir+'/')
        #command5 = command+thesubdir+subdir+'/'
        #subdir5 = bashout( command5 ).rstrip().splitlines()
        #for subsubdir in subdir5 :
            #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')


if debug : 
    if not done : print( subdirlist3 )
    if done : print( filelist )
if not done :
    for subdir2 in subdirlist3:
        lists = bashout( command+subdir2 ).rstrip().splitlines()
        for line in lists :
            if rootfile in line : filelist.append(subdir2+line)


#for thefile in filelist:

        #print( thefile )


#print( filelist )
index = 1
ext = '_v29.txt'
past = ''
select =  filelist[0].split("/")
outfile = select[index] + ext
outf = open( outfile, 'w' )
print( outfile )
for thefile in filelist:
    select =  thefile.split("/")
    if select[index] not in past :
        past = select[index]
        outfile = select[index] + ext
        outf.close()
        outf = open( outfile, 'w' )
        print( outfile )
    outf.write( thefile + '\n' )
outf.close()


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
