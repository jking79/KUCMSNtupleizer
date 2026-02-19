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
    year = ""
    if(args.year == "2016"):
        year = "R16"
    elif(args.year == "2018"):
        year = "R18"
    elif(args.year == "2017"):
        year = "R17"
    elif(args.year == "2022"):
        year = "R22"
    elif(args.year == "2023"):
        year = "R23"
    elif(args.year == "2024"):
        year = "R24"
    elif(args.year == ""):
        year = ""
    else:
        print("Year",args.year,"not found")
        exit()
    data = False
    MCbkg = False
    MCsig = False
    data_samples = ["MET","JetMET","JetHT","EGamma","JetMET0","JetMET1"]
    MCbkg_samples = ["GJets","QCD","DiPJBox","TTXJets","WJets","ZJets"]
    MCsig_samples = ["gogoG","gogoZ","gogoGZ","sqsqG"]
    data = (args.inputSample in data_samples)
    MCbkg = (args.inputSample in MCbkg_samples)
    MCsig = (args.inputSample in MCsig_samples)
    if(data):
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_DataPD_Files_List.txt"
    elif(MCbkg):
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_BG_SVIPM100_Files_List.txt"
    elif(MCsig):
        inputMainList = cmsswpath+"ntuple_master_lists/KUCMS_Ntuple_Master_SMS_Sig_Files_List_v33.txt"
    else:
    	print("Input list for samplei "+args.inputSample+" not found")
    	exit()
    print("Using main file list",inputMainList)
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
            #print("list",data[1])
            #if sublist specified, only look at that one
            #if(inputList != "" and inputList not in data[1]):
            #    continue
            if args.inputSample not in data[1]:
                continue
            if args.mini:
                if "MINI" not in data[1]:
                    continue
            else:
                if "MINI" in data[1]:
                    continue
            if args.fast:
                if "FAST" not in data[1]:
                    continue
            else:
                if "FAST" in data[1]:
                    continue
            if args.filter not in data[1]:
                continue
            if year != "" and year not in data[1]:
                continue
            if args.era != "" and args.era not in data[1]:
                continue
            if args.slice != "" and args.slice != "" not in data[1]:
                continue
            if args.mGl != "" and "mGl-"+args.mGl not in data[1]:
                continue
            if args.mN2 != "" and "mN2-"+args.mN2 not in data[1]:
                continue
            if args.mN1 != "" and "mN1-"+args.mN1 not in data[1]:
                continue
            if args.ctau != "":
                if "p" not in args.ctau and "." in args.ctau:
                    args.ctau = args.ctau.replace("p",".")
                if "ct"+args.ctau not in data[1]:
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
            ofilename = "condor_"+samplename+"_"+ofiletag
            fulldirname = odir+dirname+"/"+samplename+"/"+ofiletag
            #print("dirname",dirname,"samplename",samplename)
            print("Preparing sample directory: {0}".format(fulldirname))
            ##### Create a workspace (remove existing directory) #####
            if os.path.exists(fulldirname):
            	print("Removing existing directory: {0}".format(fulldirname))
            	shutil.rmtree(fulldirname)
            # Create directories for work area.
            SH.createWorkArea(fulldirname)
            eosdir = "root://cmseos.fnal.gov//store/group/lpcsusylep/jaking/"+eospath
            #eventnums = SH.eventsSplit(eosdir, inputlist, args.split,True, args.maxnevts)
            print("inputlist",inputlist)
            filearr = SH.filesSplit(eosdir, inputlist, args.maxnevts, args.maxnfiles)
            # grab relevant flags
            flags = " --xsec "+str(xsec)+" --dataSetKey "+key+" --gluinoMass "+str(gluinomass)+" --N2Mass "+str(n2mass)+" --timeCaliTag "+timeCaliTag+" --MCweight "+str(mc_wt)

            mctype = 1
            if(MCbkg or MCsig):
                flags += " --hasGenInfo"
                mctype = 0
            flags += " --MCtype "+str(mctype)

            if( not args.runPSICHE):
                flags += " --noBHC"
            if(args.noSV or "AL1NpSC" in args.inputSample):
                flags += " --noSV"
            if(args.HLTPathsOff):
                flags += " --HLTPathsOff"

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
            print("No file list found for",args)
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
    parser.add_argument('--inputSample','-i',help='Ntuple sample to create skims from',choices=['DiPJBox', 'DTBoson','GJets','TTXJets','QCD','WJets','ZJets','gogoG','gogoZ','gogoGZ', 'sqsqG','MET','EGamma','JetHT','MET_AL1NpSC','JetMET1','JetMET0'])
    parser.add_argument("--year",help="run year",choices = ["2016","2017","2018","2022","2023"],default="")
    parser.add_argument('--slice',help='HT slice (ie for GJets or QCD), mass range (ie for DiPJBox), or subprocess (ie for TTXJets)',default='')
    parser.add_argument('--era',help='run era (data only)',default='')
    parser.add_argument('--mGl',help='gluino mass for signal',default='')
    parser.add_argument('--mN2',help='neutralino2 mass for signal',default='')
    parser.add_argument('--mN1',help='neutralino1 mass for signal',default='')
    parser.add_argument('--ctau',help='ctau (only for gogoZ and gogoGZ samples)',default='')
    parser.add_argument('--filter',help='filter that was used for ntuple creation',default='SVIPM100')
    parser.add_argument('--output','-o',help='output label',default=None)
    #parser.add_argument('--split','-s',help="condor job split",default=0,type=int)
    parser.add_argument('--maxnevts',help="maximum number of events per job",default=-999,type=int)
    parser.add_argument('--maxnfiles',help="maximum number of files total",default=-999,type=int)
    parser.add_argument('--verbosity','-v',help="verbosity",default=0)
    parser.add_argument('--runPSICHE',help='run with creating PSICHE jets (default = off)',default=False,action='store_true')
    #parser.add_argument('--hasGenInfo',help='set hasGenInfo flag',default=False,action="store_true")
    parser.add_argument('--genSigPerfect',help='set genSigPerfect flag',default=False,action="store_true")
    parser.add_argument('--noSV',help='do not run SV collection',default=False,action="store_true")
    #parser.add_argument('--noSVorPho',help='set noSVorPho flag',default=False,action="store_true")
    parser.add_argument('--fast',help='run over fastsim aod ntuples',default=False,action='store_true')
    parser.add_argument('--mini',help='run over miniaod ntuples',default=False,action='store_true')
    parser.add_argument('--noSmear',help='run master list with null calib/smear tag',default=False,action='store_true')
    parser.add_argument('--HLTPathsOff',help='turn off HLT paths',default=False,action='store_true')
    args = parser.parse_args()
    generateSubmission(args)

if __name__ == "__main__":
    main()
