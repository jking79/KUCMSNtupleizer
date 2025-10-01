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
justin = '/store/user/janguian/gogoG/'
eosll = 'eos root://cmseos.fnal.gov ls '

#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_JetHT_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJETS_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_QCD_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_DEG_R17_v18/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R17E_MET75_v20/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GMSB_R17_MET100_v23/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_MET_R17_AL1IsoPho_v22/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_SMS_GlGl_v23/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_gogoG_Justin_None_v22/'
#command = eosll+justin+'SMS-GlGl_mGl-2000_mN2-1950_mN1-1900_gam_N2ctau-0p1_AOD/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_gogoG_Sig_IPM100_v24/CRAB_UserFiles/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_WZ_R18_IPM100_v24/WZ_TuneCP5_13TeV-pythia8/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_TTJets_R18_IPM100_v24/'
#command = eosll+mdis+'KUCMSNtuple/kucmsntuple_WJetsToLNu_R18_IPM100_v24/'
command = eosll+mdis+'KUCMSNtuple/kucmsntuple_GJets_R18_IPM100_v24/'

rootfile = '.root'

#dirselect = 'GMSB'
#dirselect = 'GMSB_L-250'
#dirselect = 'AOD'
#dirselect = 'WJetsToLNu_HT-800'
#dirselect = 'QCD'
#dirselect = 'GMSB_L-300TeV'
#dirselect = 'DYJetsToLL_M-50'
#dirselect = 'TTJets'
#dirselect = 'JetHT'
dirselect = 'GJets'
#dirselect = 'GJets_HT-100To200'
#dirselect = 'DiPhotonJetsBox'
#dirselect = 'Run2018C'
#dirselect = 'MET'
#dirselect = 'DoubleEG'
#dirselect = 'MET_R17E_MET75'
#dirselect = 'AOD_Run2017E-17Nov2017-v1'
#dirselect = 'MINIAOD'
#dirselect = 'TT'

#dirselect = ''

debug = True
#debug = False

targdirs = []

def readNextDir( inputdir ) :

    subdirlist = []
    dirls = bashout( inputdir ).splitlines()
    for subdir in dirls :    
        subdirlist.append( inputdir + '/' + subdir )
    return subdirlist


dirls = readNextDir( command )
print( dirls )
for line in dirls:
    if dirselect in line :
        while rootfile not in dirselect :
            nextline 
 
        targdirs.append( line )
print( targdirs )
while  
   
    select =  line2.split("Tune")
    #outfile = select[0] + '_filelist.txt'
    outfile = 'kuntuple_' + select[0] + 'v24.txt'
    #outfile = 'kuntuple_' + select[0] + '_AL1IsoPho_R17_v24.txt'
    #print( outfile )
    outf = open( outfile, 'w' )
    for thefile in filelist:
    	outf.write( thefile + '\n' )
    outf.close()
