#!/bin/bash
set -euo pipefail

dump_dir () {
    echo "==== DIRECTORY DUMP ($(pwd)) ===="
    ls -alh
    echo "================================"
}

echo "Running on: $(hostname)"
echo "PWD: $(pwd)"
echo "Args: $@"

echo "=== Proxy Debug ==="
echo "X509_USER_PROXY = ${X509_USER_PROXY:-unset}" 
voms-proxy-info -timeleft || true
echo "=== Proxy Debug ==="

echo "Running on: $(hostname)"
echo "PWD: $(pwd)"
echo "Args: $@"

dump_dir


INPUTFILE="${1}"
SCALE="${2}"
USECALI="${3}"
SMEAR="${4}"
NAMEEXT="${5}"

# unpack resolution config
tar -xzf configres.tgz

echo "After extracting configres.tgz:"
dump_dir

# set CMSSW / LPC runtime environment
export SCRAM_ARCH=el9_amd64_gcc11
source /cvmfs/cms.cern.ch/cmsset_default.sh
source ./config/cmssw_setup_connect.sh
cmssw_setup sandbox-CMSSW_13_3_3.tar.bz2
#source ./config/setup_RestFrames_connect.sh

echo "After CMSSW setup:"
echo "PWD: $(pwd)"
dump_dir

# extra runtime libs used by your environment
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMSSW_BASE/src/BayesianClustering/lib

# move config dirs to where executable expects them
mv config/ntuple_master_lists .
mv config/ecal_config .

echo "After moving config dirs:"
dump_dir

chmod +x config/makeKUCMSTimeResolution.obj

echo "Tracking new files created during execution:"
BEFORE=$(mktemp)
AFTER=$(mktemp)

find . -type f | sort > "$BEFORE"

./config/makeKUCMSTimeResolution.obj "${INPUTFILE}" "${SCALE}" "${USECALI}" "${SMEAR}" "${NAMEEXT}"


find . -type f | sort > "$AFTER"

echo "=== New files created ==="
comm -13 "$BEFORE" "$AFTER"
echo "=== End new file list ==="


echo "After executable run:"
dump_dir

OUTFILE="res2dPlotsTFileCondor.root"
RENAMED="res2dPlots${NAMEEXT}.root"

if [[ -f "${OUTFILE}" ]]; then
    mv "${OUTFILE}" "${RENAMED}"
else
    echo "ERROR: expected output file ${OUTFILE} not found"
    echo "Searching for any ROOT files:"
    find . -type f -name "*.root" -ls
    echo "Full directory dump:"
    dump_dir
    exit 1
fi

echo "Finished. Output:"
ls -lh "${RENAMED}"

