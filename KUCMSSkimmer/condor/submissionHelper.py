#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import os
import re
import ROOT
import numpy as np

# Create directory if it does not exist.
def makeDir(dir_name):
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)

# Create directories for work area.
def createWorkArea(sampleDir):
	os.makedirs(sampleDir)
	os.mkdir(sampleDir+"/src")
	os.mkdir(sampleDir+"/log")
	os.mkdir(sampleDir+"/out")

def getDataSetName(pathToList):
        tmp = pathToList.split('/')
        tmp = tmp[-1].split('.')
        #print(tmp[0])
        return tmp[0]

# Write the header of the condor submit file
def writeSubmissionBase(subf, dirname, ofilename):
        subf.write("universe = vanilla\n")
        subf.write("executable = execute_script.sh\n")
        subf.write("output = ./"+dirname+"/log/job.$(Process).out\n")
        subf.write("error = ./"+dirname+"/log/job.$(Process).err\n")
        subf.write("log = ./"+dirname+"/log/job.log\n")
        #include tarball with CMSSW environment
        subf.write("transfer_input_files = /uscms/home/mlazarov/nobackup/sandboxes/sandbox-CMSSW_13_3_3.tar.bz2, config.tgz, \n")
        subf.write("should_transfer_files = YES\n")
        #subf.write("request_memory = 4096\n")
        subf.write("when_to_transfer_output = ON_EXIT\n")
        outnames = []
        outnames.append(ofilename+".$(Process).root")
        outname = ""
        for o in outnames:
            outname += o+", "
        #remove last comma and space
        outname = outname[:-2]
        #subf.write("transfer_output_files = "+outname+"\n")
        subf.write("transfer_output_files = "+outname+"\n")
        # need to supply absolute path for remap
        #absCWD = os.path.abspath(".") # these cwd give the wrong abs path, there is something special in the environment
        #absCWD = os.getcwd()
        absCWD = os.popen('pwd').readline().rstrip()
        #print("abs path is "+ absCWD)
        remap = ""
        for o in outnames:
            remap += o+"="+absCWD+"/"+dirname+"/out/"+o+";"
            #print("remap is "+ remap)
	        #print("outname is "+outname)
        #remove last semicolon 
        remap = remap[:-1]
        subf.write("transfer_output_remaps = \""+remap+"\"\n")

#splits by event number
def eventsSplit(eosdir, infile, nChunk, filelist, nevtsmax = -999):
    nevts = -1
    if nevtsmax != -999:
        nevts = nevtsmax
    else:
        if filelist:
            #should split by event number in file
            tree = ROOT.TChain("tree/llpgtree")
            with open(infile,"r") as f:
                lines = f.readlines()
                for line in lines:
                    rfile = ROOT.TFile.Open(eosdir+line[:-1]) #remove new line character
                    tree.AddFile(rfile.GetName())
        else:
            #should split by event number in file
            rfile = ROOT.TFile.Open(infile)
            tree = rfile.Get("tree/llpgtree")
        nevts = tree.GetEntries()
    evts = range(nevts+1)
    if nChunk == 0:
        nChunk += 1
    print("Splitting each file into "+str(nChunk)+" jobs ")
    #return array of pairs of evtFirst and evtLast to pass as args into the exe to run
    #make sure to count first event in every chunk after first
    arr = [[min(i)-1, max(i)] for i in np.array_split(evts,nChunk)]
    #set first entry to 0
    arr[0][0] = 0
    return arr


# Write each job to the condor submit file.
def writeQueueList( subf, inFile, ofilename, evts, flags ):
    if evts == 0 or evts is None:
        print("No events found")
        return
    #.root is set in exe
    outFileArg = ofilename+".$(Process)"
    
    jobCtr=0
    subf.write("\n\n\n")
    subf.write("queue Arguments from (\n")
    for e in evts:
            inFileArg = " -i "+inFile
            #Args = "Arguments ="+inFileArg+" "+flags+" --evtFirst "+str(e[0])+" --evtLast "+str(e[1])+" -o "+outFileArg+"\n"
            Args = inFileArg+" "+flags+" --evtFirst "+str(e[0])+" --evtLast "+str(e[1])+" -o "+outFileArg+"\n"
            #subf.write("\n\n\n")
            subf.write("###### job"+str(jobCtr)+ "######\n")
            subf.write(Args)
            #subf.write("Queue\n")
            jobCtr=jobCtr+1

    subf.write(")")
