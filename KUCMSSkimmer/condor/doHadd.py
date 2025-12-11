import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir","-d",help="top directory that has directories of root files to hadd i.e. Output/GMSB/GMSB_AOD_v13_GMSB_L-150TeV_Ctau-200cm_AODSIM_RunIIFall17DRPremix/jets/",required=True)
    parser.add_argument("--force",help='force remake of outfile',action='store_true')
    parser.add_argument("--big",help='run routine for larger rootfiles',action='store_true')
    parser.add_argument("--dryRun",help="dry run, print commands without running",action="store_true")
    args = parser.parse_args()

    scratch_path = "/uscmst1b_scratch/lpc1/3DayLifetime/$USER/"	
    
    cmdHadd = "hadd -k -d "+scratch_path+" -j 4"
    for d in os.scandir(args.dir):
        print(d.path)
        if not os.path.exists(d.path+"/out"):
            continue
        oname = d.name
        proc = d.path[d.path.find("/")+1:]
        proc = proc[proc.find("/")+1:]
        proc = proc[:proc.find("/")]
        proc = proc[:proc.rfind("_AOD")]

        bashfilename = "haddScripts/doHadd_"+proc+"_"+oname+".sh"
        bashfile = open(bashfilename,"w")
        #write cmds to bash script
        if(args.big):
            for i in range(10):
                oname = d.name+"_"+proc+"_"+str(i)+".root"
                oname = "condor_"+oname
                oname = scratch_path+"/"+oname
                cmd = ""
                #check if file exists
                if os.path.exists(oname):
                    if(args.force):
                        cmd = cmdHadd+" -f"
                    else:
                        print(oname+" exists ")
                        continue
                else:
                    cmd = cmdHadd
                #print(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root")
                bashfile.write(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root\n")
                #if not args.dryRun:
                #    os.system(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root")	
                #print("Wrote to "+oname)
            big_oname = d.path
            big_oname = big_oname[:big_oname.find("/"+d.name)]
            big_oname = scratch_path+"/"+big_oname[big_oname.rfind("/")+1:]+"_"+d.name+".root"
            #check if file exists
            if os.path.exists(oname):
            	if(args.force):
            		cmd = cmdHadd+" -f"
            	else:
            		print(oname+" exists ")
            #print(cmd+" "+oname+" "+d.path+"/condor_*.root")	
            #bashfile.write(cmd+" "+oname+" "+scratch_path+"/condor_*.root\n")
            
            bashfile.write(cmd+" "+big_oname+" "+oname[:oname.rfind("_")]+"_*.root\n")
            #if not args.dryRun:
            #    os.system(cmd+" "+oname+" "+d.path+"/condor_*.root")	
            #print("Wrote to "+oname)
            bashfile.close()
        else:
            oname = d.path
            oname = oname[:oname.find("/"+d.name)]
            oname = d.path+"/"+oname[oname.rfind("/")+1:]+"_"+d.name+".root"
            print("oname",oname,"d.name",d.name,"proc",proc,"d.path",d.path)
            cmd = cmdHadd
            #check if file exists
            if os.path.exists(oname):
            	if(args.force):
            		cmd += " -f"
            	else:
            		print(oname+" exists ")
            		continue
            #print(cmd+" "+oname+" "+d.path+"/out/*.root")	
            bashfile.write(cmd+" "+oname+" "+d.path+"/out/*.root\n")	
            #if not args.dryRun:
            #    os.system(cmd+" "+oname+" "+d.path+"/out/*.root")
            #print("Wrote to "+oname)
            bashfile.close()
        print("To hadd run:")
        print("source",bashfilename)

if __name__ == "__main__":
	main()

