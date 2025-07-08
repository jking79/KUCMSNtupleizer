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

#mspc = '/store/user/jaking/'
#justin = '/store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/'
justin = '/store/user/lpcsusylep/jaking/KUCMSNtuple/sqsqG/'
#justin = '/store/user/janguian/gogoG/'
eosll = 'eos root://cmseos.fnal.gov ls '
eosopen = 'root://cmseos.fnal.gov/'

rootdir = eosopen + justin 
command = eosll+justin

rootfile = '.root'

#dirselect = '2500_mN2-1500_mN1-1000'

dirselect = ''

debug = True
#debug = False

targdirs = []

dirls = bashout( command ).splitlines()
print( dirls )
print( '************************************************')
for line in dirls:
	#print( line )
	if dirselect in line : targdirs.append( line )
    #targdirs.append( line )
print( targdirs )

for line2 in targdirs :

    subdirlist1 = []
    subdirlist2 = []
    subdirlist3 = []
    filelist = []
    theFileList = ''

    #for mydir in targdirs:
    print( '-------------------------------------------------')
    print( line2 )
    print( '-------------------------------------------------')
    #subdirlist1.append( line+'/' )
    
    #if debug : print( subdirlist1 )
    #for thesubdir in subdirlist1 :
    thesubdir = line2
    command2 = command+thesubdir+'/'
    if debug : print( command2 )
    subdir2 = bashout( command2 ).rstrip().splitlines()
    for subdir in subdir2 : 
    	command3 = command+thesubdir+'/'+subdir+'/'
    	subdir3 = bashout( command3 ).rstrip().splitlines()
    	for subsubdir in subdir3 : 
    		subdirlist2.append(thesubdir+'/'+subdir+'/'+subsubdir)
   
    #if debug : print( subdirlist2 ) 
    #for thesubdir2 in subdirlist2 :
    #    command4 = command+thesubdir2+'/'
    #    subdir4 = bashout( command4 ).rstrip().splitlines()
    #    #print( thesubdir+subdir2+'/0000/' )
    #    for subdir in subdir4 :
    #        subdirlist3.append(thesubdir2+'/'+subdir+'/')
            #command5 = command+thesubdir+subdir+'/'
            #subdir5 = bashout( command5 ).rstrip().splitlines()
            #for subsubdir in subdir5 :
                #subdirlist3.append(thesubdir+subdir+'/'+subsubdir+'/')
    
    
    if debug : print( subdirlist2 )
    for subdir2 in subdirlist2:
    	lists = bashout( command+subdir2 ).rstrip().splitlines()
    	for lline in lists :
    		if rootfile in lline : filelist.append(rootdir+subdir2+'/'+lline)
   
    select =  line2.split("_gam")
    outfile = 'filelist_' + select[0] + '_v24.txt'
    #outfile = 'kuntuple_' + select[0] + '_AL1IsoPho_R17_v24.txt'
    #print( outfile )
    outf = open( outfile, 'w' )
    for thefile in filelist:
    	outf.write( thefile + '\n' )
    outf.close()



