#!/bin/bash
# run_skim.sh
# Usage: ./run_skim.sh <arg1> <arg2> <arg3>

echo ">>> Starting job on $(date)"
echo ">>> Host: $(hostname)"
echo ">>> Args: $1 $2 $3"

# --- Setup CMSSW environment ---
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
SCRAM_ARCH=el9_amd64_gcc12

# Unpack CMSSW tarball
tar -xzf CMSSW_13_3_3.tar.gz
cd CMSSW_13_3_3/src
eval `scramv1 runtime -sh`

# Go to KUCMSSkimmer directory
cd KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer

# --- Run the skimmer ---
echo ">>> Running ./skim.obj $1 $2 $3"
./skim.obj "$1" "$2" "$3"

# --- Move only .root files to top-level directory ---
mv *.root ${_CONDOR_SCRATCH_DIR}/ 2>/dev/null || true

echo ">>> Finished at $(date)"

