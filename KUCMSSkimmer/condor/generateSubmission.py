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
    year = ""
    if(args.year == "2018"):
        year = "R18"
    elif(args.year == "2017"):
        year = "R17"
    elif(args.year == "2022"):
        year = "R22"
    else:
        print("Year",args.year,"not found")
        exit()
    reco_date = {}
    reco_date["2017"] = "_17Nov2017"
    reco_date["2017_MET"] = "-09Aug2019_UL2017_rsb-v1"
    reco_date["2017_DEG"] = "-09Aug2019_UL2017-v1"
    reco_date["2018_DEG"] = "-15Feb2022"
    reco_date["2018_MET"] = "-15Feb2022_UL2018-v1"
    reco_date["2022_MET"] = "-27Jun2023-v2"
    reco_date["2018_MC"] = "RunIISummer20UL18"
    reco_date["2017_MC"] = "RunIIFall17DRPremix"
    reco_date["2018_EGamma"] = ""
    reco_date["2018"] = ""
    reco_date[""] = ""
    data = False
    MCbkg = False
    MCsig = False
    inputList = ""
    recotag = ""
    #MC bkgs
    if(args.inputSample == "GJets"):
        MCbkg = True
        inputList = "GJets_"+year+"_SVIPM100_v31_GJets_HT-"+args.slice
    elif(args.inputSample == "QCD"):
        MCbkg = True
        inputList = "QCD_"+year+"_SVIPM100_v31_QCD_HT"+args.slice
    elif(args.inputSample == "DiPJBox"):
        MCbkg = True
        inputList = "DiPJBox_"+year+"_SVIPM100_v31_DiPhotonJetsBox_M"+args.slice
        if args.slice == "40_80":
            inputList += "-sherpa"
    elif(args.inputSample == "TTXJets"):
        MCbkg = True
        inputList = "TTXJets_"+year+"_SVIPM100_v31_"+args.slice
    elif(args.inputSample == "WJets"):
        MCbkg = True
        inputList = "WJets_"+year+"_SVIPM100_v31_WJetsToLNu_HT-"+args.slice
    elif(args.inputSample == "ZJets"):
        MCbkg = True
        inputList = "ZJets_"+year+"_SVIPM100_v31_ZJetsToNuNu_HT-"+args.slice
    #MC signals
    elif(args.inputSample == "gogoG"):
        MCsig = True
        inputList = "gogoG_AODSIM"
        if(args.mGl != "" and args.mN2 != "" and args.mN1 != ""):
            inputList += "_mGl-"+args.mGl+"_mN2-"+args.mN2+"_mN1-"+args.mN1
    elif(args.inputSample == "gogoZ"):
        MCsig = True
        inputList = "gogoZ_AODSIM"
        if(args.mGl != "" and args.mN2 != "" and args.mN1 != ""):
            inputList += "_mGl-"+args.mGl+"_mN2-"+args.mN2+"_mN1-"+args.mN1
            if(args.ctau != ""):
                ctau = args.ctau
                ctau = ctau.replace(".","p")
                inputList += "_ct"+ctau 
    elif(args.inputSample == "sqsqG"):
        MCsig = True
        inputList = "sqsqG_AODSIM"
        if(args.mGl != "" and args.mN2 != "" and args.mN1 != ""):
            inputList += "_mGl-"+args.mGl+"_mN2-"+args.mN2+"_mN1-"+args.mN1
    #data
    elif(args.inputSample == "MET"):
        data = True
        inputList = "MET_"+year+"_SVIPM100_v31_MET_AOD_Run"+args.year+args.era
    elif(args.inputSample == "EGamma"):
        data = True
        inputList = "EGamma_"+year+"_InvMetPho30_NoSV_EGamma_AOD_Run"+args.year+args.era
    else:
        print("Lists for input sample",args.inputSample,"not found")
        exit()
    
    if(data):
        recotag = args.year+"_"+args.inputSample
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_DataPD_Files_List.txt"
    elif(MCbkg):
        recotag = args.year+"_MC"
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt"
    elif(MCsig): #gluino only right now!
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_SMS_Sig_Files_List.txt"
    else:
    	print("Input list"+args.inputMainList+" not found")
    	exit()
    print("Using main file list",inputMainList,"with input list",inputList)
    listpath = inputMainList[:inputMainList.rfind("/")]+"/"
    #parsing master list
    #parsing lines in master list for inputs for list-by-list processing
    condor_subs = []
    with open(inputMainList,'r') as f:
        for line in f:
            full_line = line.strip()
            if not full_line:
                continue
            if(line[0] == "#"):
                continue
            data = line.split(" ")
            eospath = data[0]
            #print("data",data[1],"inputList",inputList)
            #if sublist specified, only look at that one
            if(inputList != "" and inputList not in data[1]):
                continue
            inputlist = listpath+data[1]
            key = data[2]
            xsec = data[3]
            gluinomass = 0
            n2mass = 0
            n1mass = 0
            #get signal masses from sample name (if they exist)
            if "SMS" in inputlist:
                file = data[1]
                gluinomass = file[file.find("mGl-")+4 : file.find("_mN2-")]
                n2mass = file[file.find("mN2-")+4 : file.find("_mN1")]
                n1mass = file[file.find("mN1-")+4 : file.find("_ct")]
            mc_wt = data[6]
            timeCaliTag = data[-1][:-1] #remove new line symbol
           
            #set dirname
            dirname = data[1][:data[1].find("_")]


            #set output name
            #print("recotag",recotag)
            samplename = data[1][:-4]#+reco_date[recotag]
            ofiletag = "rjrskim"
            if args.output is not None:
                ofiletag += "_"+args.output
            ofile_inputList = ""
            #if(inputMainList != ""):
            #    ofile_inputList = inputMainList+"_"
            ofilename = "condor_"+ofile_inputList+samplename+"_"+ofiletag

            fulldirname = odir+dirname+"/"+samplename+"/"+ofiletag
            #print("dirname",dirname,"samplename",samplename)
            print("Preparing sample directory: {0}".format(fulldirname))
            ##### Create a workspace (remove existing directory) #####
            if os.path.exists(fulldirname):
            	print("Removing existing directory: {0}".format(fulldirname))
            	shutil.rmtree(fulldirname)
            # Create directories for work area.
            SH.createWorkArea(fulldirname)
            eosdir = "root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/"+eospath
            #eventnums = SH.eventsSplit(eosdir, inputlist, args.split,True, args.maxnevts)
            filearr = SH.filesSplit(eosdir, inputlist, args.maxnevts, args.maxnfiles)
            # grab relevant flags
            flags = " --xsec "+str(xsec)+" --dataSetKey "+key+" --gluinoMass "+str(gluinomass)+" --N2Mass "+str(n2mass)+" --timeCaliTag "+timeCaliTag+" --MCweight "+str(mc_wt)

            mctype = 1
            if(MCbkg or MCsig):
                flags += " --hasGenInfo"
                mctype = 0
            flags += " --MCtype "+str(mctype)

            if(args.noBHC):
                flags += " --noBHC"
            if(args.noSV):
                flags += " --noSV"

            ##### Create condor submission script in src directory #####
            condorSubmitFile = fulldirname + "/src/submit.sh"
            subf = open(condorSubmitFile, "w")
            #print("outputfile name "+ofilename)
            SH.writeSubmissionBase(subf, fulldirname, ofilename, args.max_mat, args.max_idle, args.request_memory)
            #need to remove local lpc path for actual args
            inputlist = inputlist[inputlist.rfind("/",0,inputlist.rfind("/"))+1:]
            #print("inputfilelist",inputlist)
            SH.writeQueueList(subf, ofilename, filearr, flags)
            condor_subs.append(condorSubmitFile)
            print()
        if len(condor_subs) < 1:
            print("No file list found for",inputList)
            exit()
        elif len(condor_subs) < 2:
            print("------------------------------------------------------------")
            print("Submission ready, to run use:")
            print("condor_submit "+condor_subs[0]+"\n")
        else:
            #write bash script
            #print("fulldirname",fulldirname)
            #print("odir",odir,"dirname",dirname,"samplename",samplename,"ofiletag",ofiletag)
            mult_bash_name = odir+dirname+"/"+ofiletag+"_MultiSub.sh"
            mult_bash = open(mult_bash_name, "w")
            SH.writeMultiSubScript(mult_bash, condor_subs)
            print("------------------------------------------------------------")
            print("Submission ready, to run use:")
            print("source "+mult_bash_name+"\n")
                
def main():
    # options
    parser = argparse.ArgumentParser()
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    parser.add_argument("--max_mat",help='max_materialization condor option (default: off)',default=-1)
    parser.add_argument("--max_idle",help='max_idle condor option (default: off)',default=-1)
    parser.add_argument("--request_memory",help='memory to request from condor scheduler in bits (default = 2048)',default=-1)
    #parser.add_argument("--inputList",help="list of sample lists to run over (default is SVIPM100 selection)",choices=['data','mcBkg','mcSig'])
    parser.add_argument('--inputSample','-i',help='Ntuple sample to create skims from',choices=['DiPJBox', 'DTBoson','GJets','TTXJets','QCD','WJets','ZJets','gogoG','gogoZ','sqsqG','MET','EGamma'])
    parser.add_argument("--year",help="run year",choices = ["2017","2018","2022"],default="2018")
    parser.add_argument('--slice',help='HT slice (ie for GJets or QCD), mass range (ie for DiPJBox), or subprocess (ie for TTXJets)',default='')
    parser.add_argument('--era',help='run era (data only)',default='')
    parser.add_argument('--mGl',help='gluino mass for signal',default='')
    parser.add_argument('--mN2',help='neutralino2 mass for signal',default='')
    parser.add_argument('--mN1',help='neutralino1 mass for signal',default='')
    parser.add_argument('--ctau',help='ctau (only for gogoZ samples)',default='')
    parser.add_argument('--output','-o',help='output label',default=None)
    #parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--maxnevts',help="maximum number of events per job",default=-999,type=int)
    parser.add_argument('--maxnfiles',help="maximum number of files total",default=-999,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    parser.add_argument('--noBHC',help='run without creating BHC objects (default = on)',default=False,action='store_true')
    #parser.add_argument('--hasGenInfo',help='set hasGenInfo flag',default=False,action="store_true")
    parser.add_argument('--genSigPerfect',help='set genSigPerfect flag',default=False,action="store_true")
    parser.add_argument('--noSV',help='do not run SV collection',default=False,action="store_true")
    #parser.add_argument('--noSVorPho',help='set noSVorPho flag',default=False,action="store_true")
   
    args = parser.parse_args()
    generateSubmission(args)

if __name__ == "__main__":
    main()
