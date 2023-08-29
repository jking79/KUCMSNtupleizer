// -*- C++ -*-
//
// Package:    KUCMSNtupilizer
// Class:      KUCMSNtupilizer
//
/**\class KUCMSNtupilizer LLPgammaAnalyzer.cc LLPgammaAnalyzer/plugins/KUCMSNtupilizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Jack W King III
//         Created:  Wed, 27 Jan 2021 19:19:35 GMT
//
//

Comments to be added


Ntuples currently stored at: store/user/lpcsusylep/jaking/KUCMSNtuple
- can be accessed via eosls or xrdcp (to copy via xrootd)

Ntuples can be run one of two ways:
- `cmsRun test/llpgana_mc_aod.py` will interactively run the Ntuplizer
	- make sure this is running over the correct global tag, AOD files, and number of events
- `python test/multicrab_mc_llpgana_cmssw130_aod.py` will run the Ntuplizer via condor
	- make sure this is running over the correct file sets and number of events
