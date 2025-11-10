#!/bin/bash
#wget --quiet --no-check-certificate http://stash.osgconnect.net/+zflowers/backup/sandbox-CMSSW_10_6_5-6403d6f.tar.bz2
#includes .x file to execute and cmssw_setup_connect script
tar -xzf ./config.tgz
#
#export SCRAM_ARCH=slc7_amd64_gcc700
export SCRAM_ARCH=el9_amd64_gcc11
##sets CMSSW environment variables
source /cvmfs/cms.cern.ch/cmsset_default.sh
##sets up CMSSW environment from a sandbox
source ./config/cmssw_setup_connect.sh
#may need connect version of this bash script
#source ./config/setup_RestFrames.sh
cmssw_setup sandbox-CMSSW_13_3_3.tar.bz2 
source ./config/setup_RestFrames_connect.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CMSSW_BASE/src/BayesianClustering/lib
mv config/ntuple_master_lists .
mv config/ecal_config .
mv config/fullMask60cm_Data.root .
###
####don't need lhapdf library - but may need to do this for cgal, etc.
####export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/cvmfs/cms.cern.ch/slc7_amd64_gcc700/external/lhapdf/6.2.1-pafccj3/lib/
ls
./config/runKUCMSAodSVSkimmer.obj "$@"
#pwd
#echo "--done run"
#ls
#echo "--config"
#ls config
#echo $@
#./FullClusterSkim.x "$@"

