#!/bin/bash
# run_skim.sh
# Usage: ./run_skim.sh <arg1> <arg2> <arg3>

echo ">>> Starting job on $(date)"
echo ">>> Host: $(hostname)"
echo ">>> Args: $1 $2 $3"

# --- Setup CMSSW environment ---
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh

SCRAM_ARCH=el8_amd64_gcc12

# Unpack CMSSW tarball
tar -xzf CMSSW_13_3_3.tar.gz
cd CMSSW_13_3_3/src
eval `scramv1 runtime -sh`

# --- Add RestFrames and other local libs ---
echo ">>> Setting up additional library paths..."
# Add RestFrames build lib directory (inside CMSSW tarball)
export LD_LIBRARY_PATH=${PWD}/restframe_build/lib:${LD_LIBRARY_PATH}

# Add BayesianClustering libraries if needed
export LD_LIBRARY_PATH=${PWD}/BayesianClustering/lib:${LD_LIBRARY_PATH}

# (Optional) source RestFrames setup if needed for environment variables
# source ${PWD}/../RestFrames/setup_RestFrames.sh

# Go to KUCMSSkimmer directory
cd KUCMSNtupleizer/KUCMSNtupleizer/KUCMSSkimmer

# --- Run the skimmer ---
echo ">>> Running ./runLocalCKUCMSAodSVSkimmer.obj $1 $2 $3"
./runLocalCKUCMSAodSVSkimmer.obj "$1" "$2" "$3"

# --- Move only .root files to top-level directory ---
mv *.root ${_CONDOR_SCRATCH_DIR}/ 2>/dev/null || true

echo ">>> Finished at $(date)"

