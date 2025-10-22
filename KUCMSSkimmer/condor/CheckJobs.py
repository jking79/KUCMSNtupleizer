import os
import argparse
import subprocess

# Check condor jobs
# - check number of output root files
# - check sizes of small files 
# - check log files for errors


#scan directories
def scanOutDirs(outputDir, subfolders):
	for f in os.scandir(outputDir):
		if f.is_dir():
			if os.path.isdir(f.path+"/log"):
				subfolders.append(f.path)
				continue
			else:
				scanOutDirs(f, subfolders) 
	return subfolders
# Check condor jobs
def checkJobs(outputDir,match_string):
    print("Running over the directory '{0}'.".format(outputDir))
    print("------------------------------------------------------------")
    # get list of directories
    #subfolders = [f.path for f in os.scandir(outputDir) if f.is_dir()]
    subfolders = []
    subfolders = scanOutDirs(outputDir, subfolders)
    
    # grep Queue for possible jobs submitted
    for folder in subfolders:
        if match_string is not None:
            if match_string not in folder:
                continue
        Strategy = folder.split("/")[-1]
        Object = folder.split("/")[-2]
        DataSetName = folder.split("/")[-3]
        print( "Evaluating Dataset "+DataSetName+" on "+Object+" with strategy "+Strategy)
        bash = "grep -c \"# job\" "+folder+"/src/submit.sh"
        #bash = "grep -c \"Queue\" "+folder+"/src/submit.sh"
        print("N Jobs Queued")
        #os.system(bash)
        output1 = int(subprocess.check_output(['bash','-c', bash]).decode())
        #print("the output"+ str(output))
        #print("N Output Produced")
        bash = "ls "+folder+"/out | wc -l"
        #os.system(bash)
        output2 = int(subprocess.check_output(['bash','-c', bash]).decode())
        print("Number of Jobs Queued: %4d, Number of Output Files: %4d" % (output1,output2) )
        if ( output1 > output2 ):
        	print("--- ERROR: Mismatch!!!!!!!!!! Jobs Failed!!!!!!!!!! ---")
        #print("Checking for files under 15 kB...")
        #bash = "find "+folder+"/out -type f -size -15k"
        #output3 = subprocess.check_output(['bash','-c', bash]).decode()
        #output3 = output3.split("\n")
        #output3.remove('')
        #print("Found "+str(len(output3))+" small files")
        #for smallFile in output3:
        #	print(smallFile)
        
        print("Checking Error Logs...")
        #bash = "grep -v -e \"Warning\" -e \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" "+ folder +"/log/*.err > count.txt"
        #bash = "grep -c \"Warning\" -c \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" "+ folder +"/log/*.err"
        #bash = "grep -v -e \"Warning\" -e \"WARNING\" " + folder +"/log/*.err > count.txt"
        bash = "grep -P \"ERROR | FATAL | Aborted | Break\" " + folder +"/log/*.err > "+folder+"/errcount.txt"
        os.system(bash)
        with open(folder+"/errcount.txt","r") as f:
        	print(len(f.readlines()),"jobs with errors: see "+folder+"/errcount.txt for more")	
        #output4 = subprocess.check_output(['bash','-c',bash]).decode()
        #print(output4)	
        print("------------------------------------------------------------")	
        print("")

def main():
    # options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--directory", "-d", default="Output", help="working directory for condor submission")
    parser.add_argument("--matchString","-s",default="",help="string to match for directory scan")
    options = parser.parse_args()
    directory = options.directory	
    
    checkJobs(directory,options.matchString)

if __name__ == "__main__":
    main()
