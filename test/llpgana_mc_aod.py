import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## Flags
options.register('multicrab',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'swtich to use muticrab paramters');
options.register('hasGenInfo',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to get pcalo in mc');
options.register('eventFilter','MET100',VarParsing.multiplicity.singleton,VarParsing.varType.string,'filter to use in event processing');
options.register('doSV',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to run displaced SVs');
options.register('doDisEle',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag ro run displaced electrons');
options.register('doECALTrackOnly',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to run ECAL Tracks only');

## object prep cuts
#options.register('jetpTmin',15.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'jet pT minimum cut');
#options.register('jetEtamax',3.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'jet eta maximum cut');
#options.register('jetIDmin',1,VarParsing.multiplicity.singleton,VarParsing.varType.int,'jet ID minimum cut');
#options.register('rhEmin',1.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'recHit energy minimum cut');
#options.register('phpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'photon pT minimum cut');
#options.register('phIDmin','none',VarParsing.multiplicity.singleton,VarParsing.varType.string,'photon ID minimum cut');

## lepton prep cuts
#options.register('ellowpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'electron low pT min cut');
#options.register('elhighpTmin',50.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'electron high pT min cut');
#options.register('mulowpTmin',20.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'muon low pT minimum cut');
#options.register('muhighpTmin',50.0,VarParsing.multiplicity.singleton,VarParsing.varType.float,'muon high pT minimum cut');

## GT to be used
##------------------ mc gt
#options.register('globalTag','106X_mc2017_realistic_v6',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
options.register('globalTag','106X_upgrade2018_realistic_v11_L1v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','94X_mc2017_realistic_v11',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','124X_mcRun3_2022_realistic_v12',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','124X_mcRun3_2022_realistic_postEE_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','130X_mcRun3_2023_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','133X_mcRun3_2024_realistic_v10',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for GMSB MC');
##options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');

##------------------ data gt  
#options.register('globalTag','106X_dataRun2_v36',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
##options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
##options.register('globalTag','124X_dataRun3_v15',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2022');
##options.register('globalTag','106X_dataRun2_v24',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2018UL');
#112X_mcRun3_2021_realistic_v16

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

#outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_Objectified_Test_v24.root' # ntuplizer test
#outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_Objectified_v23_triggertest3_fullgen_oottrue.root' # ntuplizer test
#outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_JustinAOD_Zll_ct10_fullgen_oottrue.root' # ntuplizer test
#outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_gogol_test_v16.root' # ntuplizer test
outfilename = 'kucmsntuple.root'
if options.multicrab == True : 
    outfilename = 'kucmsntuple.root'

options.register('outputFileName',outfilename,VarParsing.multiplicity.singleton,VarParsing.varType.string,'output file name created by cmsRun');

## parsing command line arguments
options.parseArguments()

## Define the CMSSW process
from Configuration.StandardSequences.Eras import eras
process = cms.Process(options.processName,eras.Run2_2018)

## Load the standard set of configuration modules
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')

#process.load("Geometry.CaloEventSetup.CaloTowerConstituents_cfi")
#process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
#process.load('Configuration.StandardSequences.L1Reco_cff')
#process.load('Configuration.StandardSequences.Reconstruction_Data_cff')

#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load("TrackingTools.TrackAssociator.DetIdAssociatorESProducer_cff")
process.load('Configuration.StandardSequences.EndOfProcess_cff')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 2
process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 100000

from KUCMSNtupleizer.KUCMSNtupleizer.TrackAssociator_cfi import tkAssocParamBlock

process.load('KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi import *

process.load('KUCMSNtupleizer.KUCMSNtupleizer.DisplacedElectrons_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.DisplacedElectrons_cfi import *

process.load('KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi import *

## Define the input source
aodpath_1k_450_100k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_1k_450_10k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath_125_25_15k = '/store/mc/Run3Winter20DRPremixMiniAOD/HTo2LongLivedTo4b_MH-125_MFF-25_CTau-15000mm_TuneCP5_14TeV_pythia8/AODSIM/110X_mcRun3_2021_realistic_v6-v2/'
aodpath21_1k_450_100k = '/store/mc/Run3Summer21DRPremix/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-10000mm_TuneCP5_14TeV-pythia8/AODSIM/120X_mcRun3_2021_realistic_v6-v2/'

lpcpath_350_600 = 'file:/eos/uscms/store/mc/RunIIFall17DRPremix/GMSB_L-350TeV_Ctau-600cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/'
lpcpath_GMSB = 'file:/eos/uscms/store/mc/RunIIFall17DRPremix/'
gmsbaodsim = '_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/'
gmsbaodsim2 = '_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v2/'

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #'/store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-0p001cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/009A3F2E-C7B9-E811-9FCA-441EA171A998.root',

        #'file:Met_UL18B_AOD_973EEF0C-44AB-E94A-8591-04DCD00D8B4B.root',
        #'root://cmseos.fnal.gov//store/user/janguian/gogoZ/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD/250525_004955/0000/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD_4.root',
        #'root://cmseos.fnal.gov//store/user/janguian/gogoZ/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD/250529_191831/0000/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD_6.root',
        #'root://cmseos.fnal.gov//store/user/janguian/gogoZ/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD/250601_053241/0000/SMS-GlGl_mGl-2000_mN2-400_mN1-350_Zff_N2ctau-0p001_AOD_89.root',

        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/SlepSnuCascade_3.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/Gluino_N2N1_LLPN2_TuneCP5_13p6TeV-madgraphMLM_AOD.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_13_0_13/src/PPD-Run3Summer22wmLHEGS-00014_AOD.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/gotest4_AOD.root', 
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/GlGl10k_testing/SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/SMS-GlGl_mGl-1p0_mN2-250_ct5_AOD.root',
        #'file:SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD_1.root',
        #'file:SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD_2.root',
        #'file:SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD_3.root',
        #'file:SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD_4.root',
        #'file:SMS-GlGl_mGl-1p0_mN2-250_ct0_AOD_5.root',
        #'file:/uscms/home/janguian/nobackup/CMSSW_12_4_14_patch3/src/SMS-GlGl_mGl-1p0_mN2-250_Zll_ct10_AOD.root',
        #'root://cmseos.fnal.gov//store/user/janguian/gogo/SMS-GlGl_mGl-1p0_mN2-250_Zll_ct10_AOD/241120_224153/0000/SMS-GlGl_mGl-1p0_mN2-250_Zll_ct10_AOD_80.root',

	#'root://cms-xrd-global.cern.ch//store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-200cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/80000/58DCC006-3FB2-E811-94AE-AC1F6B0DE454.root',	

        #'root://cmsxrootd-site.fnal.gov//store/data/Run2022C/JetHT/AOD/27Jun2023-v2/25310000/161a4ef3-6d8f-4026-bc32-12473677119e.root',
        #'root://cmsxrootd-site.fnal.gov//store/data/Run2018A/JetHT/MINIAOD/12Nov2019_UL2018-v2/100000/07082E78-7805-EB44-BB4D-1D3AB8892FB7.root',

        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/6021A8AE-9B41-7A4F-8E33-470F73F3D4FC.root',

        'file:WJets_72B9C618-FE23-1E41-872E-57314D7CB454.root',

        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/72B9C618-FE23-1E41-872E-57314D7CB454.root',
	#'/store/mc/RunIISummer20UL18RECO/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/06CC60F3-D74E-174E-8AD7-DCDAB14FB26F.root',
	#'/store/data/Run2018C/JetHT/AOD/15Feb2022_UL2018-v1/2530001/B2C0793B-3D7B-104F-B34B-87D010EFDEF4.root',	


        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall17DRPremix/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/110000/1E29BF8B-5F60-E811-AD1D-0663CE00010C.root',

        ## HTo2LongLivedTo4b

        #'/store/mc/Run3Winter21DRMiniAOD/HTo2LongLivedTo4b_MH-1000_MFF-450_CTau-100000mm_TuneCP5_14TeV-pythia8/MINIAODSIM/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/280000/17bd2d71-8a76-46c5-947a-7ea2b1df44b6.root'

	# AOD MC GMSB

        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall17DRPremix/GMSB_L-100TeV_Ctau-0_1cm_TuneCP5_13TeV-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/100000/2C55A98D-E4D7-E811-AC40-002590491B1E.root',


        # AOD Data MET
        #'file:/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/test/Met_UL18B_AOD_973EEF0C-44AB-E94A-8591-04DCD00D8B4B.root'      
        #'file:/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/test/MetPD_003A2484-A2DC-E711-9D0A-02163E019C46.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD/250123_145920/0000/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD_99.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD/250607_191347/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_2.root'

	#lpcpath_350_600+'120000/80762156-99D6-E811-8942-34E6D7E3879B.root',
        #lpcpath_350_600+'120000/322875DC-DDD6-E811-8C5F-001E675A68C4.root',
        #lpcpath_350_600+'120000/603C58F0-A9D6-E811-8907-E0071B74AC00.root',
        #lpcpath_350_600+'120000/66DB15F0-DDD6-E811-A7A3-90E2BACBAD58.root',
        #lpcpath_350_600+'120000/8A9DEBDE-DDD6-E811-838E-D4AE526DF2E1.root',
        #lpcpath_350_600+'120000/A8DE3FA7-99D6-E811-92C5-34E6D7BEAF0E.root'

         # AODSIM DPJB model

         #'/store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/0B4130F1-3785-8043-957B-E0C2F86E1A50.root',

	# AODSIM GMSB model		

         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'10000/162DBEEE-DC29-E911-843B-0CC47A745294.root',

         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'10000/BCE6A6F2-A929-E911-AEA9-24BE05C63651.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/22F81F5B-3D33-E911-906C-0CC47AD24D28.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/3CDA2489-9132-E911-B4E3-0025905D1E08.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10000cm'+gmsbaodsim2+'100000/44349E62-F131-E911-97A6-848F69FD09D7.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/26A7404E-C5DA-E811-92D3-001E675A6AB8.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/66971865-6ADA-E811-A9C2-002590D9D8C0.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/6E716649-6ADA-E811-B9E5-0025901AC0FC.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1000cm'+gmsbaodsim+'270000/A89ECEA5-C5DA-E811-AB4B-D4AE526DDB3F.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/08A3D920-E5B2-E811-85C5-0CC47A4DEEE4.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/0C51F68D-E5B2-E811-8703-002590E39F2E.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-0_1cm'+gmsbaodsim+'100000/2C55A98D-E4D7-E811-AC40-002590491B1E.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-10cm'+gmsbaodsim+'80000/2097D13F-30B3-E811-9FB3-0025905C96EA.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/605B95ED-8AD8-E811-B640-001E67A3F70E.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/68C5DD84-EAD6-E811-9028-44A842CFD60C.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/846C6D28-D4D8-E811-A8BE-00266CFFBEB4.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-1200cm'+gmsbaodsim+'60000/8EC1497F-E9D6-E811-A390-44A842CFC98B.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-200cm'+gmsbaodsim+'120000/44485BE6-E2D7-E811-9916-1866DA89095D.root',# <<<<<
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-200cm'+gmsbaodsim+'120000/AEA9923B-9ED8-E811-87D6-34E6D7BDDECE.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/301550E4-81B5-E811-8362-FA163EF08F5B.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/34D601E8-3CB2-E811-B2BA-24BE05C6B701.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/4ACDB9DD-C0AF-E811-A364-24BE05C6D731.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-400cm'+gmsbaodsim+'80000/7411AE64-81B5-E811-9F51-FA163E301A96.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-600cm'+gmsbaodsim+'60000/E25EAAE9-A5D6-E811-8502-34E6D7E05F0E.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-600cm'+gmsbaodsim+'60000/F65FFA8E-09D7-E811-81E5-A4BF0102A4F5.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'270000/B624F4E1-86D8-E811-BD73-A0369F83633E.root',
#         lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'270000/DE490468-30D8-E811-99AC-0025B3E015D2.root',
         #lpcpath_GMSB+'GMSB_L-100TeV_Ctau-800cm'+gmsbaodsim+'80000/2C3F963F-3BB2-E811-BA81-D8D385AF8902.root',

         #'file:BCB550D6-CAB3-5C4F-8866-77897305A646.root',#aodpath_125_25_15k
         #cmssw12XX only #'file:bc04e7b9-31c7-4bec-a396-258ba40b8bd5.root',#aodpath21_1k_450_100k

		 # HTo2LongLivedTo4b

         #aodpath_1k_450_100k+'10000/07C1D360-FDE3-B04C-8E7A-DDE4275C7F04.root',
         #aodpath_1k_450_100k+'10000/2EBC0785-9656-2B43-A9D2-FFEB5032E66A.root',
         #aodpath_1k_450_100k+'10000/4002AC07-8CA1-0643-9237-F38929578E9D.root',

         #aodpath_1k_450_10k+'10000/06A10966-0ABE-A94D-A5C0-BA43FB64B0AB.root',
         #aodpath_1k_450_10k+'10000/190B92B5-4416-A346-B0A9-B9980F15563C.root',
         #aodpath_1k_450_10k+'50000/4E0AC838-A6FF-B549-8441-A86B04B2BEFA.root',

         #aodpath_125_25_15k+'240000/2543A8DF-4540-FE45-A496-9EF8D7918E35.root',
         #aodpath_125_25_15k+'240000/6299C20E-A359-8B4F-A4C3-858E50BC461E.root',
         #aodpath_125_25_15k+'240000/6EA51191-B27D-B547-920E-66EF94292870.root',

		  ## EGamma

        #'/store/data/Run2017E/DoubleEG/AOD/17Nov2017-v1/20000/0009F8BE-23D3-E711-84AA-02163E011C79.root'
        #'/store/data/Run2017E/DoubleEG/AOD/17Nov2017-v1/20000/02C3735A-7DD3-E711-83C0-0025904C66E4.root'
        #'/store/data/Run2022A/EGamma/AOD/16Jun2023-v1/2820000/11cf75d2-bdf0-4768-bd5e-830c2b03655c.root'
        #'/store/data/Run2023E/EGamma/AOD/PromptReco-v1/000/372/597/00000/c97d22a6-fc8d-416f-92b4-01e27f9776c4.root'
        #'/store/data/Run2017E/MET/AOD/17Nov2017-v1/50000/00864810-19DD-E711-A884-02163E01A63A.root'

        #'/store/data/Run2018B/MET/AOD/15Feb2022_UL2018-v1/2530000/C9433911-6C3B-7740-B758-B01EF63573DA.root'

        ),##<<>>fileNames = cms.untracked.vstring
        #skipEvents=cms.untracked.uint32(300),
)##<<>>process.source = cms.Source("PoolSource",


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#TT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250))#KT
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#KT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#QT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5000))#BT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#SM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#MS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#MD
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#LG
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#FL

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))
		
genInfo = True
#genInfo = False
if options.multicrab == True : genInfo = options.hasGenInfo		   

ecalIsoInputs = 'RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'

filterselect = 'None'
#filterselect = 'AL1DisSV'
#filterselect = 'SVIPMet100'
#filterselect = 'MET100'
#filterselect = 'AL1IsoPho'
#filterselect = 'IsoPhoMet100'
#filterselect = 'AL1SelEle'

dosv = True
#dosv = False
if options.multicrab == True : dosv = options.doSV

#dode = True
dode = False
if options.multicrab == True : dode = options.doDisEle

#doet = True
doet = False
if options.multicrab == True : doet = options.doECALTrackOnly
if doet : 
    dosv = False
    dode = False

#probeout = True
probeout = False
if options.multicrab == True : probeout = False

if options.multicrab == True : filterselect = options.eventFilter


#genMomChase = True
genMomChase = False
if options.multicrab == True : genMomChase = False
#makeTrigList = True
makeTrigList = False
if options.multicrab == True : makeTrigList = False

print( "Using options : mutlicrab = ",options.multicrab," geninfo = ",genInfo," filter = ",filterselect )
print( "Using options : momChase = ",genMomChase," trgiList = ",makeTrigList," probeout = ",probeout ) 
print( "Using options : doDisEle = ",dode," doSVs = ",dosv," doECALTrackOnly = ",doet )
print( "Using options : globalTag = ",options.globalTag )
print( "With output file name : ",options.outputFileName )

# Make the tree 
process.tree = cms.EDAnalyzer("KUCMSNtupilizer",

                              ## default ?  paramter
                              tkAssocParamBlock,

                              ## flags
                              hasGenInfo = cms.bool(genInfo),
                              doGenMotherChase = cms.bool(genMomChase),
                              makeTriggerList =  cms.bool(makeTrigList),
                              doProbeOut = cms.bool(probeout),
                              doSVModule = cms.bool(dosv),
                              doDisEleModule = cms.bool(dode),
                              doECALTrackOnly = cms.bool(doet),

                              ##skim type selectuon 
                              fltrSelection = cms.string(filterselect),
                            
                              #Triggers
                              triggerList = cms.vstring( #"hltPFMET100", "hltMETClean100", "hltHIPhoton20Eta3p1" ),   
                                  "hltL1sSingleEGNonIsoOrWithJetAndTauNoPS",
                                  "hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter",
                                  "hltEG60EtFilter",
                                  "hltEG60HEFilter",
                                  "hltEG60R9Id90CaloIdLIsoLR9IdFilter",
                                  "hltEG60R9Id90CaloIdLIsoLClusterShapeFilter",
                                  "hltEG60R9Id90CaloIdLIsoLEcalPFClusterIsoFilter",
                                  "hltEG60R9Id90CaloIdLIsoLHcalPFClusterIsoFilter",
                                  "hltEG60R9Id90CaloIdLIsoLHollowTrackIsoFilter",
                                  "hltEG60R9Id90CaloIdLIsoLDisplacedIdFilter",
                                  "hltHT175Jet10",
                                  "hltPFHT350Jet15",
                                  "hltL1sSingleAndDoubleEG",
                                  "hltL1sSingleAndDoubleEGNonIsoOr",
                                  "hltL1sSingleAndDoubleEGor",
                                  "hltL1sSingleEG15",
                                  "hltL1sSingleEG18",
                                  "hltL1sSingleEG24",
                                  "hltL1sSingleEG26",
                                  "hltL1sSingleEG34to45",
                                  "hltL1sSingleEG34to50",
                                  "hltL1sSingleEG40to50",
                                  "hltL1sSingleEGor",
                                  "hltL1sTripleEGOrDoubleEGOrSingleEG",
                                  "hltEG20EtFilterLooseHoverE",
                                  "hltEG20EtL1TripleEGFilter",
                                  "hltEG20HEFilterLooseHoverE",
                                  "hltEG20HEL1TripleEGFilter",
                                  "hltEG20L1SingleEGLowETOrEtFilter",
                                  "hltEG20L1SingleEGLowETOrEtaREtFilter",
                                  "hltEG30EBHE10R9Id50b80eHEFilter",
                                  "hltEG30EBL1SingleAndDoubleEGOrEtFilter",
                                  "hltEG30EBR9Id50b80eR9IdFilter",
                                  "hltEG30EIso15HE30EcalIsoLastFilter",
                                  "hltEG30EtFilterLooseHoverE",
                                  "hltEG30EtL1TripleEGFilter",
                                  "hltEG30HE30HEFilter",
                                  "hltEG30HEFilterLooseHoverE",
                                  "hltEG30HEL1TripleEGFilter",
                                  "hltEG30L1IsoEGerJetC34drMin0p3EtFilter",
                                  "hltEG30L1SingleAndDoubleEGOrEtFilter",
                                  "hltEG30L1SingleAndDoubleEGWithTauWithJetEtFilter",
                                  "hltEG30LHE12R9Id50b80eHEFilter",
                                  "hltEG30LR9Id50b80eR9IdFilter",
                                  "hltEG30PVHE10R9Id50b80eHEFilter",
                                  "hltEG30PVR9Id50b80eR9IdFilter",
                                  "hltEG30PVrealANDHE10R9Id50b80eHEFilter",
                                  "hltEG30PVrealANDR9Id50b80eR9IdFilter",
                                  "hltHT130Jet30",
                                  "hltHT200Jet30",
                                  "hltPFHT180Jet30",
                                  "hltPFHT250Jet30",
                                  "hltPFMET50",
                                  "hltPFMET70",
                                  "hltPFMET90",
                                  "hltPFMET100",
                                  "hltPFMET110",
                                  "hltPFMET120",
                                  "hltPFMET130",
                                  "hltPFMET140",
                                  "hltPFMET200",
                                  "hltPFMET250",
                                  "hltPFMET300",
                                  "hltPFMHTNoMuTightID70",
                                  "hltPFMHTNoMuTightID90",
                                  "hltPFMHTNoMuTightID100",
                                  "hltPFMHTNoMuTightID110",
                                  "hltPFMHTTightID120",
                                  "hltPFMHTNoMuTightID120",
                                  "hltPFMHTNoMuTightID130",
                                  "hltPFMHTNoMuTightID140",
                                  "hltPFMHTTightID120",
                                  "hltPFMHTTightID130",
                                  "hltPFMHTTightID140",
                                  "hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilter",
                                  "hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterLowPt",
                                  "hltL4PromptDisplacedDijetFullTracksHLTCaloJetTagFilterMidPt",
                                  "hltPFMETNoMu60",
                                  "hltPFMETNoMu100",
                                  "hltPFMETNoMu110",
                                  "hltPFMETNoMu120",
                                  "hltPFMETNoMu130",
                                  "hltPFMETNoMu140",
                              ),#triggerList = cms.vstring(

                              ## set rechit energy limits
                              minRHEi = cms.double(0.5),
                              minRHEf = cms.double(0.5),

                              ## additional collections
                              ## tracks
                              ogGeneralTracks = cms.InputTag("generalTracks"),
                              ogGsfTracks = cms.InputTag("electronGsfTracks"),
                              ecalTracks = cms.InputTag("ecalTracks", "ecalTracks"),
                              tracks = cms.InputTag("ecalTracks", "ecalGeneralTracks"),
                              gsfTracksSrc = cms.InputTag("ecalTracks", "ecalGsfTracks"),
                              displacedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs"),
                              displacedTracks = cms.InputTag("displacedElectrons", "displacedCandidateTracks"),
                              muonEnhancedTracks = cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
                              combinedMuonTracks = cms.InputTag("muonEnhancedTracks", "combinedMuonTracks"),
                              sip2DMuonEnhancedTracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
                              ## vertices
                              #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              vertices = cms.InputTag("offlinePrimaryVertices"),
                              timedSVs = cms.InputTag("timedSVs", "timedSecondaryVertices"),
                              ## pfcandidates
                              #pfcandidates = cms.InputTag("packedPFCandidates"),
                              pfcandidates = cms.InputTag("particleFlow"),
                              particleflow = cms.InputTag("particleFlow",""),
                              pfcanphomap = cms.InputTag("particleFlow","photons"),
                              pfcanootphomap = cms.InputTag("particleFlow","photons"),
                              pfcanelemap = cms.InputTag("particleFlow","electrons"),
                              ## rho
                              #rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
                              rho = cms.InputTag("fixedGridRhoAll"),
                              ## conversions
                              conversions = cms.InputTag("allConversions"),
                              ## beamSpot
                              beamSpot = cms.InputTag("offlineBeamSpot"),
                              ## trigger
                              triggerResults = cms.InputTag("TriggerResults"),
                              #triggerObjects = cms.InputTag("TriggerResults"),
                              triggerEvent = cms.InputTag("hltTriggerSummaryAOD"),
                              ## METs
                              #mets = cms.InputTag("slimmedMETs"),
                              mets = cms.InputTag("pfMet"),
                              ## jets
                              #Charge Hadron Subtracted : charged particles from non-primary vertices (pileup) are removed before clustering.
                              #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
                              #jets = cms.InputTag("slimmedJets"),
                              #jets = cms.InputTag("ak4PFJets"),
                              jets = cms.InputTag("ak4PFJetsCHS"),
                              #jets = cms.InputTag("ak8PFJetsCHS"),
                              calojets = cms.InputTag("ak4CaloJets",""),
                              ## electrons
                              #electrons = cms.InputTag("slimmedElectrons"),
                              displacedElectrons = cms.InputTag("displacedElectrons", "displacedElectrons"),
                              #signalDisplacedElectrons = cms.InputTag("displacedElectrons", "signalDisplacedElectrons"),
                              electrons = cms.InputTag("gedGsfElectrons"),
                              eleMVAIDLooseMap = cms.InputTag("PhotonIDProdGED", "PhotonCutBasedIDLooseEM"),
                              effAreasConfigFile = cms.FileInPath( ecalIsoInputs ),
                              ## muons
                              #muons = cms.InputTag("slimmedMuons"),
                              muons = cms.InputTag("muons"),
                              ## photons
                              #gedPhotons = cms.InputTag("slimmedPhotons"),
                              gedPhotons = cms.InputTag("gedPhotons"),
                              #gedPhotons = cms.InputTag("photons"),
                              phoCBIDLooseMap = cms.InputTag("PhotonIDProd", "PhotonCutBasedIDLooseEM"),
                              #ootPhotons = cms.InputTag("slimmedOOTPhotons"),
                              ootPhotons = cms.InputTag("ootPhotons"),
                              ## ecal recHits
                              #recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                              #recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                              recHitsEB = cms.InputTag("reducedEcalRecHitsEB"),
                              recHitsEE = cms.InputTag("reducedEcalRecHitsEE"),
                              ## superclusters
                              superClusters = cms.InputTag("ecalTracks", "displacedElectronSCs"),
                              #superClusters = cms.InputTag("particleFlowEGamma"),
                              otherSuperClusters = cms.InputTag("particleFlowEGamma"),
                              #otherSuperClusters = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
                              ##otherSuperClusters = cms.InputTag("hybridSuperClusters", "uncleanOnlyHybridSuperClusters"),#51/1000
                              #otherSuperClusters = cms.InputTag("correctedHybridSuperClusters"),
                              #ootSuperClusters = cms.InputTag("reducedEgamma", "reducedOOTSuperClusters"),
                              ootSuperClusters = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"),
                              ## caloclusters
                              #caloClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
                              caloClusters = cms.InputTag("particleFlowEGamma", "EBEEClusters"),
                              ## gen info
                              genEvt = cms.InputTag("generator", ""),
                              gent0 = cms.InputTag("genParticles", "t0"),
                              genxyz0 = cms.InputTag("genParticles", "xyz0"),
                              pileups = cms.InputTag("addPileupInfo", ""),
                              #Phoronpileups = cms.InputTag("mixData", ""),
                              genParticles = cms.InputTag("genParticles", ""),
                              genjets = cms.InputTag("ak4GenJets","")
                              
)##<<>>process.tree = cms.EDAnalyzer("LLPgammaAnalyzer_aod"

process.load('RecoMET.METFilters.metFilters_cff')
#process.myGlobalSuperTightHalo2016Filter = process.globalSuperTightHalo2016Filter.clone( taggingMode = True ) 
process.Flag_goodVertices = cms.Path( process.goodVertices )
process.Flag_globalSuperTightHalo2016Filter = cms.Path( process.globalSuperTightHalo2016Filter )
process.hbheSequence = cms.Sequence( process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter )
process.Flag_HBHENoiseFilter = cms.Path( process.hbheSequence )
process.Flag_HBHENoiseIsoFilter = cms.Path( process.HBHENoiseIsoFilter )
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path( process.EcalDeadCellTriggerPrimitiveFilter )
process.Flag_BadPFMuonFilter = cms.Path( process.BadPFMuonFilter )
process.Flag_BadPFMuonDzFilter = cms.Path( process.BadPFMuonDzFilter )
process.Flag_hfNoisyHitsFilter = cms.Path( process.hfNoisyHitsFilter )
process.Flag_BadChargedCandidateFilter = cms.Path( process.BadChargedCandidateFilter )
process.Flag_eeBadScFilter = cms.Path( process.eeBadScFilter )
process.Flag_ecalBadCalibFilter = cms.Path( process.ecalBadCalibFilter )

process.setFlags = cms.Sequence( process.goodVertices +
                                 process.globalSuperTightHalo2016Filter +
                                 process.HBHENoiseFilterResultProducer * process.HBHENoiseFilter +
                                 process.HBHENoiseIsoFilter +
                                 process.EcalDeadCellTriggerPrimitiveFilter +
                                 process.BadPFMuonFilter +
                                 process.BadPFMuonDzFilter +
                                 process.hfNoisyHitsFilter +
                                 process.BadChargedCandidateFilter +
                                 process.eeBadScFilter +
                                 process.ecalBadCalibFilter )                               

process.kuEcalTracks = cms.Sequence( ecalTracks )
process.kuSV = cms.Sequence( muonEnhancedTracks )
process.kuDisEle = cms.Sequence( displacedElectrons )

# Set up the path
#process.ecalTracks_path = cms.Path(ecalTracks)
#process.displacedElectrons_path = cms.Path(displacedElectrons)
#process.muonEnhancedTracks_path = cms.Path(muonEnhancedTracks)

process.kuDisplaced_path = cms.Path()
if ( dosv and not dode ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks + process.kuSV )
if ( not dosv and dode ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks + process.kuDisEle )
if ( dosv and dode ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks + process.kuDisEle + process.kuSV )
if ( doet ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks )

process.setFlags_path = cms.Path(process.setFlags)
process.tree_step = cms.EndPath(process.tree)
process.endjob_step = cms.EndPath(process.endOfProcess)

#if dosv :
process.schedule = cms.Schedule( process.kuDisplaced_path, process.setFlags_path, process.tree_step, process.endjob_step )
#else :
#    process.schedule = cms.Schedule( process.setFlags_path, process.tree_step, process.endjob_step )


#process.schedule = cms.Schedule( process.ecalTracks_path,
#                                 process.displacedElectrons_path,
#                                 process.muonEnhancedTracks_path,
#                                 process.Flag_goodVertices,
#                                 process.Flag_globalSuperTightHalo2016Filter,
#                                 process.Flag_HBHENoiseFilter,
#                                 process.Flag_HBHENoiseIsoFilter,
#                                 process.Flag_EcalDeadCellTriggerPrimitiveFilter,
#                                 process.Flag_BadPFMuonFilter,
#                                 process.Flag_BadPFMuonDzFilter,
#                                 process.Flag_hfNoisyHitsFilter,
#                                 process.Flag_BadChargedCandidateFilter,
#                                 process.Flag_eeBadScFilter,
#                                 process.Flag_ecalBadCalibFilter,
#                                 process.tree_step,
#                                 process.endjob_step,
#)#process.schedule

process.options = cms.untracked.PSet()

#do not add changes to your config after this point (unless you know what you are doing)
#from FWCore.ParameterSet.Utilities import convertToUnscheduled
#process=convertToUnscheduled(process)

# customisation of the process.
#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
#process = miniAOD_customizeAllData(process)
# End of customisation functions

process.options = cms.untracked.PSet( 
    #numberOfThreads = cms.untracked.uint32(4), 
    #numberOfStreams = cms.untracked.uint32(4), 
    TryToContinue = cms.untracked.vstring('ProductNotFound'), 
    #wantSummary = cms.untracked.bool(True)
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
)#process.options

