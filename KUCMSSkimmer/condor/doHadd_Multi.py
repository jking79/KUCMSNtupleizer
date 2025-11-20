import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir","-d",help="top directory that has directories of root files to had (default = Output/)",default="Output/")
    parser.add_argument("--proc","-p",help='process to hadd',required=True)
    parser.add_argument("--tag","-t",help="root outfile tag to hadd",required=True)
    parser.add_argument("--force",help='force remake of outfile',action='store_true')
    parser.add_argument("--big",help='run routine for larger rootfiles',action='store_true')
    parser.add_argument("--dryRun",help="dry run, print commands without running",action="store_true")
    args = parser.parse_args()

    proc = args.proc
    oname = args.tag
    bashfilename = "haddScripts/doHadd_"+proc+"_"+oname+".sh"
    bashfile = open(bashfilename,"w")
    cmdHadd = "hadd -k -d /uscmst1b_scratch/lpc1/3DayLifetime/mlazarov/ -j 4 -v 0"
    #provide directory Output/proc
    topdir = args.dir + "/" + args.proc
    cmds = []
    for dd in os.scandir(topdir):
        #only look at directories (not bash scripts, csvs, root files, etc)
        if not dd.is_dir():
            continue
        for d in os.scandir(dd):
            #print("path",d.path)
            #only look at output root verions we specify
            if args.tag not in d.path:
                continue
            #only look at directories (not bash scripts, csvs, root files, etc)
            if not d.is_dir():
                continue
            #check that directory has the path to the output root files we need
            if not os.path.exists(d.path+"/out"):
                continue
            oname = d.name
            #print("oname",oname)

            #bashfilename = "haddScripts/doHadd_"+proc+"_"+oname+".sh"
            #bashfile = open(bashfilename,"w")
            #append cmds to write
            if(args.big):
                for i in range(10):
                    oname = d.name+"_"+proc+"_"+str(i)+".root"
                    oname = "condor_"+oname
                    oname = d.path+"/"+oname
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
                    #bashfile.write(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root\n")
                    #write pieces of hadd
                    cmds.append(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root\n")
                    #if not args.dryRun:
                    #    os.system(cmd+" "+oname+" "+d.path+"/out/*."+str(i)+"*.root")	
                    #print("Wrote to "+oname)
                oname = d.path
                oname = oname[:oname.find("/"+d.name)]
                oname = d.path+"/"+oname[oname.rfind("/")+1:]+"_"+d.name+".root"
                #check if file exists
                if(args.force):
                	cmd = cmdHadd+" -f"
                #print(cmd+" "+oname+" "+d.path+"/condor_*.root")	
                #save final hadd
                cmds.append(cmd+" "+oname+" "+d.path+"/condor_*.root\n")
                cmds.append("echo \"Hadded "+oname+"\"\n")
                cmds.append("\n")
                #if not args.dryRun:
                #    os.system(cmd+" "+oname+" "+d.path+"/condor_*.root")	
                #print("Wrote to "+oname)
                #bashfile.close()
            else:
                oname = d.path
                oname = oname[:oname.find("/"+d.name)]
                oname = d.path+"/"+oname[oname.rfind("/")+1:]+"_"+d.name+".root"
                #print("oname",oname,"d.name",d.name,"proc",proc,"d.path",d.path)
                cmd = cmdHadd
                #check if file exists
                if(args.force):
                	cmd += " -f"
                #print(cmd+" "+oname+" "+d.path+"/out/*.root")
                #save hadd cmd
                cmds.append(cmd+" "+oname+" "+d.path+"/out/*.root\n")	
                cmds.append("echo \"Hadded "+oname+"\"\n")
                cmds.append("\n")
                #if not args.dryRun:
                #    os.system(cmd+" "+oname+" "+d.path+"/out/*.root")
                #print("Wrote to "+oname)
    if(len(cmds) < 1):
        print("No files in",args.dir+"/"+args.proc+"/","with root filetag",args.tag)
        exit()
    else:
        for cmd in cmds:
            bashfile.write(cmd)
    bashfile.close()
    print("To hadd run:")
    print("source",bashfilename)

if __name__ == "__main__":
	main()

