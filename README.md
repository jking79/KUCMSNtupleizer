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

CMSSW install for KUCMSNtupilizer:

cmsrel CMSSW_13_0_X
cd CMSSW_13_0_X/src 
cmsenv 
git cms-init

( make sure the analysis package is cloned into its own folder inside scr/  example( timing/timing/…  ) )
( yes ..  its required by CMSSW )
( initial scram b -j8 must be in src/ ?  after maybe scram-ed  from plugins/ or test/ …. )


mkdir KUCMSNtupilizer
cd KUCMSNtupilizer/
git clone https://github.com/jking79/KUCMSNtupilizer.git
cd ../
scram b -j 8
