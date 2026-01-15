import os
import argparse
import subprocess
import numpy as np

# Check condor jobs
# - check number of output root files
# - check sizes of small files 
# - check log files for errors


#find the job number which needs resubmission
def findMissingJobNums(outputlist, expectedNumFiles):
    expectedJobNums = set(np.arange(expectedNumFiles))
    #print("JOBS EXPECTED", expectedJobNums)
    #print("LS OUTPUT")
    #print(outputlist)
    joblist = outputlist.split()
    jobNumlist = [int(s.split('.')[-2]) for s in joblist]
    jobNumlist = set(jobNumlist)
    missingJobNums = sorted(expectedJobNums - jobNumlist)
    print("Missing Jobs:", missingJobNums)
    return missingJobNums

def locateSubmissionLines(targetFile, targetJobNum):
    f = open(targetFile)
    job_lines =[]
    ###### job0######
    start_string = '###### job'+str(targetJobNum)+'######'
    stop_string = '###### job'+str(targetJobNum+1)+'######'
    #print("searching for:", start_string)
    recording=False
    for line in f:
            # Use .strip() for comparison to ignore leading/trailing whitespaces, including newlines
            if start_string == line.strip() :
                recording = True
                            
            if stop_string == line.strip():
                recording = False
                break  # Stop recording after the first end string

            if recording:
                job_lines.append(line.rstrip('\n')) # Append line, removing only trailing newline

    #replace job lines $Process with jobnumber
    job_lines = [s.replace('$(Process)', str(targetJobNum)) for s in job_lines]
    f.close()
    return job_lines

def getJobHeaderLines(targetFile, targetJobNum):
    f = open(targetFile)
    header_lines=[]
    stop_string = 'queue Arguments from ('
    for line in f:
        header_lines.append(line.rstrip('\n'))
        if line.strip() == stop_string:
        #    print("found header stop")
            break
    header_lines = [s.replace('$(Process)', str(targetJobNum)) for s in header_lines]
    f.close()
    return header_lines

def writeResubFile(targetresub, targetJobNum, header, args):
    resubtarget = targetresub.replace('submit.sh', 'resub_'+str(targetJobNum)+'.sh')

    with open(resubtarget, 'w') as f:
        for line in header:
            f.write(line+'\n')
        for line in args:
            f.write(line+'\n')
        f.write(')\n')
    return resubtarget

def createResubmissionScripts( targetFile, targetMissingList ):
    resubmissionScripts=[]
    #individual files for each job to retain job numbers 
    for targetJob in targetMissingList:
        header = getJobHeaderLines(targetFile, targetJob)
        resubargs = locateSubmissionLines(targetFile, targetJob)
        resubtarget = targetFile
        resubmissionScripts.append( writeResubFile(resubtarget, targetJob, header, resubargs) )
    #print("lines for submission script")
    #print(resubmissionScripts)
    return(resubmissionScripts)

def createMultiSubScript( resubmissionScripts,outputDir, match_string):
    #print("preparing multi resub script")
    resubmultiname = './'+outputDir+'/multi_resub_'+match_string+'.sh'
    with open(resubmultiname, 'w') as f :
        f.write("echo Launching Resubmission\n")
        for line in resubmissionScripts:
            f.write("condor_submit "+line+"\n")
    return resubmultiname

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
    
    multi_resub_scripts = []
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
            bashls = "ls "+folder+"/out"
            outputlist = subprocess.check_output(['bash','-c', bashls]).decode()
            missingJobNums = findMissingJobNums(outputlist, output1)
            targetSubmissionFile= folder+"/src/submit.sh"
            multi_resub_scripts.extend( createResubmissionScripts( targetSubmissionFile, missingJobNums ) )

            #joblines=locateSubmissionLines( targetSubmissionFile, missingJobNums[0])
            #print("JOBLINES")
            #print(joblines)
        #print("Checking for files under 15 kB...")
        #bash = "find "+folder+"/out -type f -size -15k"
        #output3 = subprocess.check_output(['bash','-c', bash]).decode()
        #output3 = output3.split("\n")
        #output3.remove('')
        #print("Found "+str(len(output3))+" small files")
        #for smallFile in output3:
        #	print(smallFile)
     
        print("Bypassing error log check...")
        #print("Checking Error Logs...")
        #bash = "grep -v -e \"Warning\" -e \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" "+ folder +"/log/*.err > count.txt"
        #bash = "grep -c \"Warning\" -c \"WARNING\" -e \"TTree::SetBranchStatus\" -e \"libXrdSecztn.so\" "+ folder +"/log/*.err"
        #bash = "grep -v -e \"Warning\" -e \"WARNING\" " + folder +"/log/*.err > count.txt"
        #bash = "grep -P \"ERROR|FATAL|Error|Aborted|Break\" " + folder +"/log/*.err > "+folder+"/errcount.txt"
        #os.system(bash)
        #with open(folder+"/errcount.txt","r") as f:
        #	print(len(f.readlines()),"jobs with errors: see "+folder+"/errcount.txt for more")	
        #output4 = subprocess.check_output(['bash','-c',bash]).decode()
        #print(output4)	
        print("------------------------------------------------------------")	
        print("")
    multiname = createMultiSubScript( multi_resub_scripts ,outputDir, match_string)
    print("Resubmission script ready, do:")
    print("source "+multiname)


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
