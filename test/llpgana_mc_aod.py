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
options.register('runera','Run3',VarParsing.multiplicity.singleton,VarParsing.varType.string,'filter to use in event processing');

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
#options.register('globalTag','106X_upgrade2018_realistic_v11_L1v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','94X_mc2017_realistic_v11',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
options.register('globalTag','124X_mcRun3_2022_realistic_postEE_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for sig22 MC');
#options.register('globalTag','124X_mcRun3_2022_realistic_postEE_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','130X_mcRun3_2023_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','133X_mcRun3_2024_realistic_v10',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for GMSB MC');
##options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','140X_dataRun3_v17',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for R3 22,23');

##------------------ data gt  
#options.register('globalTag','106X_dataRun2_v36',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
##options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','124X_dataRun3_v15',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2022');
#options.register('globalTag','140X_dataRun3_v17',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
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

process.load('PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 2
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
if options.multicrab == True : process.MessageLogger.cerr.FwkReport.reportEvery = 100000

from KUCMSNtupleizer.KUCMSNtupleizer.TrackAssociator_cfi import tkAssocParamBlock

process.load('KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi import *

#process.load('KUCMSNtupleizer.KUCMSNtupleizer.DisplacedElectrons_cfi')
#from KUCMSNtupleizer.KUCMSNtupleizer.DisplacedElectrons_cfi import *

process.load('KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.MuonEnhancedTracks_cfi import *

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))

## Define the input source

#filelist_path = "met_18D_missing_resubmissions.txt"
#with open(filelist_path) as f:
#    readFiles = cms.untracked.vstring(line.strip() for line in f if line.strip() and not line.startswith("#"))

process.source = cms.Source("PoolSource",
    #fileNames=readFiles,
    fileNames = cms.untracked.vstring(

        #root://cmsxrootd.fnal.gov/
        #'root://xrootd-cms.cern.ch//store/data/Run2022C/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/a346eb6b-170c-4e4c-b403-e57a70cec129.root',
        #'root://xrootd-cms.cern.ch//store/data/Run2022C/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/113cb54d-48b1-43f5-96c1-f4725ee59eca.root',
        #'root://xrootd-cms.infn.it//store/data/Run2022C/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/a346eb6b-170c-4e4c-b403-e57a70cec129.root',

    #'root://cmsxrootd.fnal.gov//store/data/Run2023D/JetMET1/AOD/EXODelayedJetMET-PromptReco-v2/000/370/667/00000/89178b58-f461-478d-a2cf-142842257adc.root',

        #'root://cmsxrootd.fnal.gov//store/data/Run2016G/MET/AOD/21Feb2020_UL2016-v1/20000/090027EF-E362-DA48-9977-6E3CE8BB6410.root'

        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/010FAEF3-D41F-D340-B7E9-E3BBD6498B3E.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/0FB1FACF-237F-1749-9878-CA52BC077610.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/1424BBF2-79D4-3B46-8858-DB666EA0A4A6.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/42079EE6-4686-C141-8FB5-E4590022B1A9.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/4F909A9C-185D-CC49-B9E6-18F7A289ADBD.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/50AAF0F9-800F-CC45-A12F-B8FDB3442335.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/746C3653-29BC-9444-A930-546D498E615D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/80502D38-B7DA-5F41-B466-7178313A127C.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/8CDAFC03-8548-7C45-8632-71CF322749F5.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/9051BAF8-BDFE-FB40-A9BC-53C2D1AE3190.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/93EC3493-21C4-AE4D-8305-B04B6ACAB4B9.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/2430000/F5F89F63-AFFB-9F42-A32A-C76168E02B3E.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/018F6B3E-97D2-524A-853B-D4E6F010623D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/03296F17-45E0-9247-B933-334BEC69B943.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/03C0679E-03B6-6B40-9EA2-4A956BC5D371.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/06104F02-486B-FC47-9BAA-915A44A48F8D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/08626076-8276-2440-8E62-71BF442A3330.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/0986A8F9-C54F-8C43-B5D6-9739B3F9E642.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/099C48A9-0BAA-7F47-B2F6-4F1C8053E679.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/0DF60D20-8D60-2948-BE57-49A5664BF652.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/1259A5A2-D30B-084A-AF88-76E0D622DA4F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/1C62FC6A-1CF4-6D4D-A655-45A966562DB7.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/21C90FF5-14BC-BB45-A0EC-436C7A5A8E6D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/2967E911-0086-A148-B2AE-9BF68ABBD82A.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/2EB077EB-75AC-4C4B-8487-5A738424146F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/2F3C22EE-54B0-B541-A4C5-B9689EC0BE38.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/3059175F-573A-964B-A43E-0D2A10D18E05.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/322411C4-4EF1-AB41-B9E8-881AF25DD27F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/3C89AC3C-EE94-A54A-9A11-F5D1A0AE2403.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/3D87B47C-0CE6-4045-8DDF-028CECDB1EDC.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/45798998-F49C-4141-A371-885E1F8F078E.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/479AA18E-81EA-284D-ACA0-AA16C0F35D45.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/48AFE00C-6512-CC44-B898-D117F3C62C9F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/4B29EC46-B241-474E-B53E-C5EF2DC2FBD3.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/4D7D827C-D1F1-1647-BF99-379AAA4AE824.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/507ADF4C-9788-B64B-A439-E982BEF08F2E.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/56052290-C110-524D-9B5A-6F0A2B5F0B77.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/5C869995-83A5-6C4A-944E-4034AF71E8CD.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/5E398AD2-DB52-7B45-8542-07A9BA285ADF.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/62344AB9-E0FF-C74F-A1B8-2614D1638611.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/62DE1882-7064-0746-AFD9-7C0C150785BF.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/66A05CE4-DFF0-A148-B058-78CE6BDFCBE2.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/67B8DA09-2069-BB4A-A588-7E18CFF739F2.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/691E3B11-BE46-7D4F-926E-7F5C49F3CABB.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/6C6296CA-120A-F64D-A46F-DAC3B0725F99.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/71C3E68D-976F-0D4E-B922-119D33AA34C8.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/72B177C9-BF02-CD4A-A135-8A1030937728.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/792060D9-060B-1342-8773-94A91CE76A2F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/7AB51787-1626-5546-8752-56B2202C2E4D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/7AFB8C74-D6E2-294F-BC40-54A3EBF2680A.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/7E94333D-669E-114F-9ECE-BCDA753CDF26.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/863EAF2E-2B43-D14E-AC88-B857DF48118D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/924754BA-7E72-4746-AE82-17DB6C61EEBA.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/93D16D69-9990-E44A-B36C-1C4BCA6BE408.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/9C3D4D3F-62F4-3248-A057-9457B4A6D8D5.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/A3CB66E8-F548-F448-B519-DB67F4BD16D2.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/AB1FDBC7-EB94-A64C-AE96-6C0DB80EB5AE.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/B16CBED8-3CCF-5C4A-99AE-8388A346A14F.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/B57854FC-00E7-FB49-ABDF-1801AC61D29D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/B9FB0A7E-EB8B-4B47-8423-5B81326FCEA7.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/C3C481CE-2C23-2B4A-8AA9-20D20E8C96E8.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/C462805C-50BE-B24D-948F-763BAB88AB1B.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/D500CC25-14E9-D949-AC1E-D90D22DF4EEE.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/D7832157-555D-404B-B019-22D10EB3A53C.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/D8B15542-4EC6-F945-A882-7D0D4A424909.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/DE318281-DA67-5944-B419-BBF7E69755B9.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/E22EA895-51EF-9A47-AF7B-90C220A71B78.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/E50C2ADC-DF96-074C-A6AD-1BE113A5815B.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/E9289BB4-1970-764D-A256-62776020C585.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/EC6AC7EE-7472-194B-8C77-5B8F46379FC0.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/EE033181-2E78-9B40-9E91-F90E6DB2F75C.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/EF367DBD-267D-8C47-8AB2-1E07ADD26009.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/EF69396A-46A8-7642-8C37-8440868E7F9D.root',
        #'root://cmseos.fnal.gov//store/data/Run2018D/MET/AOD/15Feb2022_UL2018-v1/510000/F2FC21D5-4320-AC44-9028-5FB2A928BDB9.root',

        #'root://cmseos.fnal.gov//store/mc/RunIISummer20UL18RECO/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2520000/0184321F-08CA-3740-9700-89D3C169F771.root',

        #SMS_Sig
        'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/SMS-GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD/250317_234753/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD_41.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD/250123_145920/0000/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD_99.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD/250607_191347/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_2.root',
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD/250607_191347/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_4.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD/250610_030144/0000/SMS-GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD_78.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-1500_mN2-500_mN1-100_Zll_N2ctau-0p1_AOD/250123_153207/0000/SMS-GlGl_mGl-1500_mN2-500_mN1-100_Zll_N2ctau-0p1_AOD_83.root'

        # WJets
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/6021A8AE-9B41-7A4F-8E33-470F73F3D4FC.root',

        #'file:root_files/WJets_72B9C618-FE23-1E41-872E-57314D7CB454.root',

        # DiPhotonJetsBox
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/50000/E1C9ABB4-D708-D14E-8E97-803013F321E6.root',
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2520000/E28F5CDC-017D-2D49-87DA-864981D2E02F.root',

        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/72B9C618-FE23-1E41-872E-57314D7CB454.root',

        #GJets
	#'/store/mc/RunIISummer20UL18RECO/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/06CC60F3-D74E-174E-8AD7-DCDAB14FB26F.root',

        #QCD
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall17DRPremix/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/110000/1E29BF8B-5F60-E811-AD1D-0663CE00010C.root',

        #### AOD Run3   MUST BE IN CMSSW  14 or 15 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ####3'root://cmseos.fnal.gov//store/data/Run2023B/JetMET1/AOD/EXODelayedJetMET-PromptReco-v1/000/366/892/00000/dccdb4dd-2e6e-4a1f-ae81-70af305d9f2e.root',

        # AOD Data MET
        #'file:root_files/Met_UL18B_AOD_973EEF0C-44AB-E94A-8591-04DCD00D8B4B.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018B/MET/AOD/15Feb2022_UL2018-v1/25210000/CD3E4875-8B2E-ED4E-B1C9-43BB906B0555.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018B/MET/AOD/15Feb2022_UL2018-v1/2520000/F672D9DC-50E8-684C-BE01-8E73D82651FF.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018B/MET/AOD/15Feb2022_UL2018-v1/25210000/2325D400-4654-D64A-B3A3-9C96703D9207.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2022C/MET/AOD/27Jun2023-v2/2820000/35c49035-7557-4a77-99c8-065844b3649d.root'#

        # JetMet
        #'root://cmsxrootd-site.fnal.gov//store/data/Run2022F/JetMET/AOD/PromptReco-v1/000/360/390/00000/0401b399-e626-492f-9849-24eb67c7f647.root',# file open error
        #'root://cmsxrootd-site.fnal.gov//store/data/Run2022G/JetMET/AOD/19Dec2023-v1/25410000/2572a8de-86c0-4e95-bc8d-9f2bb8397188.root',
        #'root://cmsxrootd-site.fnal.gov//store/data/Run2022E/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/08fd72c7-1be9-4328-9a38-85979d340331.root',#reco
        #`'root://cms-xrd-global.cern.ch//store/data/Run2022E/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/08fd72c7-1be9-4328-9a38-85979d340331.root',
        #'root://cms-xrd-global.cern.ch//store/data/Run2018C/DisplacedJet/AOD/15Feb2022_UL2018-v1/60000/E11CE7E4-249F-0D42-AD15-344870056EF4.root',
        #'file:JetMet_22G_026f1000-0881-4fc1-a256-385c9cb6cf3d.root',
        #'file:JetMet_23C_0d0f0b4a-8fb0-4777-9aed-25f2416791e5.root',

         # AODSIM DPJB model

         #'/store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/0B4130F1-3785-8043-957B-E0C2F86E1A50.root',

        ),##<<>>fileNames = cms.untracked.vstring
        secondaryFileNames=cms.untracked.vstring()
        #skipEvents=cms.untracked.uint32(300),
)##<<>>process.source = cms.Source("PoolSource",

## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))#ST
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#TT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250))#KT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
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

genInfo = True
#genInfo = False
if options.multicrab == True : genInfo = options.hasGenInfo		   

ecalIsoInputsF17 = 'RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'
ecalruneraIsoInputsW22 = 'RecoEgamma/ElectronIdentification/data/Run3_Winter22/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_122X.txt'

filterselect = 'none'
#filterselect = 'InvMet100IP'
#filterselect = 'AL1NpSC'
#filterselect = 'AL1DisSV'
#filterselect = 'SVIPMet100'
#filterselect = 'MET100'
#filterselect = 'AL1IsoPho'
#filterselect = 'IsoPhoMet100'
#filterselect = 'AL1SelEle'
if options.multicrab == True : filterselect = options.eventFilter

runera = "Run3"  # current siganl model   !!!!!!!!  Run3 must be in CMSSW 14 or 15
#runera = "Run2" # BG models
if options.multicrab == True : runera = options.runera

dosv = True
#dosv = False # with v30 this is the only SV flag that matters, doDisEle and doECALTrackOnly are depreciated
if options.multicrab == True : dosv = options.doSV

############dode = True
dode = False # removed - only false with v30
#if options.multicrab == True : dode = options.doDisEle

doet = True ## with v30 no longer optional when using ( "ecalTracks", "displacedElectronSCs" ) for SCs collection
#####doet = False  ## leave true always for now
#if options.multicrab == True : doet = options.doECALTrackOnly
#if doet : 
#    dosv = False
#    dode = False

#runMETFilters = True
runMETFilters = False

#probeout = True
probeout = False
if options.multicrab == True : probeout = False

triggerSet = "RECO"
ecalIsoInputs = ecalIsoInputsF17
if "Run3" in runera :
    triggerSet = "PAT"
    if genInfo is True : 
        triggerSet = "RECO"
    ecalIsoInputs = ecalruneraIsoInputsW22
    if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos : 
        process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR') 

#genMomChase = True
genMomChase = False
if options.multicrab == True : genMomChase = False
#makeTrigList = True
makeTrigList = False
if options.multicrab == True : makeTrigList = False

print( "Using options : mutlicrab = ",options.multicrab," geninfo = ",genInfo," filter = ",filterselect )
print( "Using options : momChase = ",genMomChase," trgiList = ",makeTrigList," probeout = ",probeout ) 
print( "Using options : doDisEle = (NU) ",dode," doSVs = ",dosv," doECALTrackOnly = (NU) ",doet )
print( "Using options : globalTag = ",options.globalTag )
print( "With output file name : ",options.outputFileName )

#test = cms.vstring( "hltPFMET100", "hltMETClean100", "hltHIPhoton20Eta3p1" )

process.eventCounter = cms.EDAnalyzer("KUCMSCounter")

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
                              RunEra = cms.string(runera),

                              ##skim type selectuon 
                              fltrSelection = cms.string(filterselect),
                            
                              #Triggers
                              triggerList = cms.vstring(),

                              #MetFilters
                              metFilters = cms.vstring(),

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
                              #displacedTracks = cms.InputTag("displacedElectrons", "displacedCandidateTracks"),
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
                              #particleflow = cms.InputTag("particleFlow",""),
                              #pfcanphomap = cms.InputTag("particleFlow","photons"),
                              #pfcanootphomap = cms.InputTag("particleFlow","photons"),
                              #pfcanelemap = cms.InputTag("particleFlow","electrons"),
                              ## rho
                              #rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
                              rho = cms.InputTag("fixedGridRhoAll"),
                              ## conversions
                              conversions = cms.InputTag("allConversions"),
                              ## beamSpot
                              beamSpot = cms.InputTag("offlineBeamSpot"),
                              ## trigger
                              #triggerFlagResults = cms.InputTag("TriggerResults","","RECO"),
                              triggerFlagResults = cms.InputTag("TriggerResults","",triggerSet),
                              #triggerFlagResults = cms.InputTag("TriggerResults"),
                              triggerHLTResults = cms.InputTag("TriggerResults","","HLT"),
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
                              #displacedElectrons = cms.InputTag("displacedElectrons", "displacedElectrons"),
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

## Trigger path list
process.tree.triggerList = cms.vstring(  #"hltPFMET100", "hltMETClean100", "hltHIPhoton20Eta3p1" ),   

    "HLT_PFMET90_PFMHT90_IDTight_v",
    "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight_v",
    "HLT_PFMET100_PFMHT100_IDTight_v",
    "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_v",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_v",
    "HLT_PFMET120_PFMHT120_IDTight_PFHT60_v",
    "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v",
    "HLT_PFMET130_PFMHT130_IDTight_v",
    "HLT_PFMETNoMu130_PFMHTNoMu130_IDTighti_v",
    "HLT_PFMET140_PFMHT140_IDTight_v",
    "HLT_PFMETNoMu140_PFMHTNoMu140_IDTight_v",

)#triglist0 = cms.vstring(

## Metfilters

process.tree.metFilters = cms.vstring( #"bla" )

    "Flag_BadChargedCandidateFilter",
    "Flag_BadChargedCandidateSummer16Filter",
    "Flag_BadPFMuonDzFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonSummer16Filter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_CSCTightHaloFilter",
    "Flag_CSCTightHaloTrkMuUnvetoFilter",
    "Flag_EcalDeadCellBoundaryEnergyFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_HcalStripHaloFilter",
    "Flag_METFilters",
    "Flag_chargedHadronTrackResolutionFilter",
    "Flag_ecalBadCalibFilter",
    "Flag_ecalLaserCorrFilter",
    "Flag_eeBadScFilter",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_globalTightHalo2016Filter",
    "Flag_goodVertices",
    "Flag_hcalLaserEventFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_muonBadTrackFilter",
    "Flag_trkPOGFilters",
    "Flag_trkPOG_logErrorTooManyClusters",
    "Flag_trkPOG_manystripclus53X",
    "Flag_trkPOG_toomanystripclus53X"


)##<<>>process.tree.metFilters = cms.vstring(

process.load('RecoMET.METFilters.metFilters_cff')

process.myGlobalSuperTightHalo2016Filter = process.globalSuperTightHalo2016Filter.clone( taggingMode = True ) 
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

process.pattriger = cms.Path( process.patTrigger )


metFilterPaths = [

    process.Flag_goodVertices,
    process.Flag_globalSuperTightHalo2016Filter,
    process.Flag_HBHENoiseFilter,
    process.Flag_HBHENoiseIsoFilter,
    process.Flag_EcalDeadCellTriggerPrimitiveFilter,
    process.Flag_BadPFMuonFilter,
    process.Flag_BadPFMuonDzFilter,
    process.Flag_hfNoisyHitsFilter,
    process.Flag_BadChargedCandidateFilter,
    process.Flag_eeBadScFilter,
    process.Flag_ecalBadCalibFilter,
    process.pattriger

]

# SVs & ecaltracks aka merged SC collection
process.kuEcalTracks = cms.Sequence( ecalTracks )
process.kuSV = cms.Sequence( muonEnhancedTracks )
#process.kuDisEle = cms.Sequence( displacedElectrons )

#process.kuDisplaced_path = cms.Path()
process.kuDisplaced_path = cms.Path( process.kuEcalTracks )
if ( dosv ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks + process.kuSV )

## Set final paths and schedule

#############3process.setFlags_path = cms.Path(process.setFlags)
process.eventcount_path = cms.Path(process.eventCounter)
process.tree_step = cms.EndPath(process.tree)
process.endjob_step = cms.EndPath(process.endOfProcess)

#done above now
#if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos :
#    process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR')

#########process.schedule = cms.Schedule( process.kuDisplaced_path, process.setFlags_path, process.tree_step, process.endjob_step )
# Flags must be run as a path in the Schedule in order to appear in the trigger results for KUCMSEventInfo 
# the Flag should return 1 if the event passed the filter and 0 if it failed the filter
process.schedule = cms.Schedule(    
    process.eventcount_path,                            
    process.kuDisplaced_path,
    process.tree_step,
    process.endjob_step, 
)#process.schedule

if runMETFilters:
    process.schedule.extend(metFilterPaths)


#process.options = cms.untracked.PSet()
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

