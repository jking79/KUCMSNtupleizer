#############################
#Thanks to Justin Anguiano (2023) for the basis 
#for this condor submission script generating code
#############################

import sys
import os
import argparse
import shutil
import submissionHelper as SH



# Create workspace and condor submit files.
def generateSubmission(args):
    # Ensure that the directory includes "/" at the end.
    odir = args.directory
    if odir[-1] != "/":
    	odir += "/"
    
    print("Directory for condor submission: {0}".format(odir))
    print("------------------------------------------------------------")
    # Create output directory for condor results if it does not exist.
    SH.makeDir(odir)
    inputMainList = ""
    ofilename = ""
    dirname = odir
    cmsswpath = os.environ["CMSSW_BASE"]+"/src/KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer/"
    print(cmsswpath)
    if(args.inputList == "data"):
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_DataPD_Files_List.txt"
        dirname += "data/"
    elif(args.inputList == "mcBkg"):
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt"
        dirname += "MCbkg/"
    elif(args.inputList == "mcSig"): #gluino only right now!
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt"
        dirname += "MCsig/"
    else:
    	print("Input list"+args.inputMainList+" not found")
    	exit()
    print("Using file list",inputMainList,"with ofilename",ofilename,"and dirname",dirname)
    listpath = inputMainList.split("/")
    listpath = listpath[-2]+"/"
    #parsing lines in master list for inputs for list-by-list processing
    with open(inputMainList,'r') as f:
        for line in f:
            full_line = line.strip()
            if not full_line:
                continue
            if(line[0] == "#"):
                continue
            data = line.split(" ")
            eospath = data[0]
            inputlist = listpath+data[1]
            key = data[2]
            xsec = data[3]
            gluinomass = 0
            n2mass = 0
            n1mass = 0
            #get signal masses from sample name (if they exist)
            if "SMS" in inputlist:
                gluinomass = data[data.find("mGl-")+4 : data.find("_mN2-")]
                n2mass = data[data.find("mN2-")+4 : data.find("_mN1")]
                n1mass = data[data.find("mN1-")+4 : data.find("_ct")]
            mc_wt = data[6]
            timeCaliTag = data[-1][:-1] #remove new line symbol
            
            #set output name
            samplename = data[1][:-4]
            ofiletag = "rjrskim"
            if args.output is not None:
                ofiletag += "_"+args.output
            ofilename = "condor_"+args.inputList+"_"+samplename+"_"+ofiletag
            

            fulldirname = dirname+samplename+"/"+ofiletag
            print("Preparing sample directory: {0}".format(fulldirname))
            ##### Create a workspace (remove existing directory) #####
            if os.path.exists(fulldirname):
            	print("Removing existing directory: {0}".format(fulldirname))
            	shutil.rmtree(fulldirname)
            # Create directories for work area.
            SH.createWorkArea(fulldirname)
            eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/"+eospath
            eventnums = SH.eventsSplit(eosdir, inputlist, args.split,True, args.maxnevts)
            if eventnums == 0 or eventnums is None:
                print("eventnums",eventnums)
                return
            # grab relevant flags
            flags = " --xsec "+str(xsec)+" --dataSetKey "+key+" --gluinoMass "+str(gluinomass)+" --N2Mass "+str(n2mass)+" --timeCaliTag "+timeCaliTag+" --MCweight "+str(mc_wt)

    


            ##### Create condor submission script in src directory #####
            condorSubmitFile = fulldirname + "/src/submit.sh"
            subf = open(condorSubmitFile, "w")
            print("outputfile name "+ofilename)
            SH.writeSubmissionBase(subf, fulldirname, ofilename)
            #need to remove local lpc path for actual args
            inputlist = inputlist[inputlist.rfind("/",0,inputlist.rfind("/"))+1:]
            print("inputfilelist",inputlist)
            SH.writeQueueList(subf, inputlist, ofilename, eventnums, flags)
            print("------------------------------------------------------------")
            print("Submission ready, to run use:")
            print("condor_submit "+condorSubmitFile+"\n")
                
def main():
    # options
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    parser.add_argument("--inputList",'-i',help="list of sample lists to run over (default is SVIPM100 selection)",required=True,choices=['data','mcBkg','mcSig'])
    parser.add_argument('--output','-o',help='output label',default=None)
    parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--maxnevts',help="maximum number of events to run over",default=-999,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    parser.add_argument('--hasGenInfo',help='set hasGenInfo flag',default=False,action="store_true")
    parser.add_argument('--genSigPerfect',help='set genSigPerfect flag',default=False,action="store_true")
    parser.add_argument('--noSVorPho',help='set noSVorPho flag',default=False,action="store_true")
   
    args = parser.parse_args()
    generateSubmission(args)

if __name__ == "__main__":
    main()
