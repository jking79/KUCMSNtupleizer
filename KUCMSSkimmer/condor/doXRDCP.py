



#EOSPREFIX=root://cmseos.fnal.gov/
#EOSPATH=/store/user/lpcsusylep/malazaro/KUCMSSkims/
#SKIMDIR=skims_v44


import subprocess
import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dir","-d",help="top directory that has directories of root files to had (default = Output/)",default="Output/")
    parser.add_argument("--proc","-p",help='process to transfer',required=True)
    parser.add_argument("--targetdir", "-td", help="target skim version directory on eos",required=True)
    parser.add_argument("--tag","-t",help="root outfile tag to hadd (skim version num)",default='')
    parser.add_argument("--dryRun",help="dry run, print commands without running",action="store_true")
    args = parser.parse_args()
    proc = args.proc
    oname = args.tag
    hadd_file = ""
    hadd_path = ""
    targetdir = args.targetdir
    destination_path = "root://cmseos.fnal.gov//store/user/lpcsusylep/malazaro/KUCMSSkims/"
    proc = args.proc
    oname = args.tag
    #provide directory Output/proc
    topdir = args.dir + "/" + args.proc
    cmds = []
    for ddd in os.scandir(topdir):
        #only look at directories (not bash scripts, csvs, root files, etc)
        if not ddd.is_dir():
            continue
        for dd in os.scandir(ddd):
            #expect an hadded file in rjrskim
           # print("dd",d.path)
            #only look at output root verions we specify
            if args.tag not in dd.path:
                continue
            #only look at directories (not bash scripts, csvs, root files, etc)
            if not dd.is_dir():
                continue
            #check that directory has the path to the output root files we need
            #if not os.path.exists(d.path+"/rjrskim"):
            #    continue
            for d in os.scandir(dd):
                #we should be in rjrskim here, look for a root file
                if d.is_file() and ('.root' in d.name):        
                    hadd_file = d.name
                    hadd_path = d.path
                    cmd = "xrdcp "
                    source = hadd_path
                    destination = destination_path+targetdir+"/"+hadd_file
                    cmds.append( [cmd, source, destination] )

    if(args.dryRun):
        for cmd in cmds:
            print( cmd[0], cmd[1], cmd[2] )
        exit() 
    
        # This runs the xrdcp command and waits for it to complete
    for cmd in cmds:
        subprocess.run(["xrdcp", cmd[1], cmd[2]], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
if __name__ == "__main__":
    main()

