import os, re
import FWCore.ParameterSet.Config as cms
  
### CMSSW command line parameter parser
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('python')

## Flags
options.register('hasGenInfo',True,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to get pcalo in mc');
#options.register('hasGenInfo',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to get pcalo in mc');

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
options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','124X_dataRun3_v15',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2022');
#options.register('globalTag','106X_dataRun2_v24',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2018UL');
#112X_mcRun3_2021_realistic_v16

## processName
options.register('processName','TREE',VarParsing.multiplicity.singleton,VarParsing.varType.string,'process name to be considered');

#outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14.root' # ntuplizer test
outfilename = 'gmsb_AODSIM_KUCMSNtuplizer_Objectified_v14_pfecal_oottrue.root' # ntuplizer test

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

#process.load('RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff')
#from RecoEgamma.EgammaPhotonProducers.conversionTrackSequence_cff import *
#process.load('RecoEgamma.EgammaPhotonProducers.allTimingConversions_cfi')

#process.load('TimingWithSVs.SVsTimingProducer.svproducer_cfg')

## Message Logger settings
process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.destinations = ['cout', 'cerr']
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load('KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi')
from KUCMSNtupleizer.KUCMSNtupleizer.ECALTracks_cfi import *

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
                            #eventsToProcess = cms.untracked.VEventRange("1:55050-1:55200"),
                            fileNames = cms.untracked.vstring(
                                'file:root/GMSB_L-150TeV_Ctau-0p001cm_Fall17_AODSIM-2.root',
                                'file:root/GMSB_L-150TeV_Ctau-0p001cm_Fall17_AODSIM.root',
                            ),##<<>>fileNames = cms.untracked.vstring
)##<<>>process.source = cms.Source("PoolSource",


## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(5))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))#TTi
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#LT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500))#US
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(12500))#VS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(25000))#SM
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100000))#MS
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250000))#MD
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2500000))#LG
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#FL

# Set the global tag depending on the sample type
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = options.globalTag  

## Create output file
## Setup the service to make a ROOT TTree
process.TFileService = cms.Service("TFileService", fileName = cms.string(options.outputFileName))
				   
# Make the tree 
process.tree = cms.EDAnalyzer("KUCMSNtupilizer",
                              tkAssocParamBlock,
                              ## flags
                              hasGenInfo = cms.bool(options.hasGenInfo),
                              #hasGenInfo = cms.bool(False),
                              minEvtMet = cms.double(150.0),
                              #minEvtMet = cms.double(50.0),
                              minRHEi = cms.double(0.0),
                              minRHEf = cms.double(0.2),
                              ## additional collections
                              ## tracks
                              tracks = cms.InputTag("ecalTracks", "ecalGeneralTracks"),
                              gsfTracksSrc = cms.InputTag("ecalTracks", "ecalGsfTracks"),
                              displacedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs"),
                              ## vertices
                              #vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              vertices = cms.InputTag("offlinePrimaryVertices"),
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
                              triggerObjects = cms.InputTag("TriggerResults"),
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
                              electrons = cms.InputTag("gedGsfElectrons"),
                              eleMVAIDLooseMap = cms.InputTag("PhotonIDProdGED", "PhotonCutBasedIDLooseEM"),
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
                              #superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
                              superClusters = cms.InputTag("particleFlowEGamma"),
                              otherSuperClusters = cms.InputTag("particleFlowSuperClusterECAL", "particleFlowSuperClusterECALBarrel"),
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

# Set up the path
process.ecalTracks_path = cms.Path(ecalTracks)
process.content = cms.EDAnalyzer("EventContentAnalyzer")
process.content_step = cms.Path(process.content)
process.tree_step = cms.EndPath(process.tree)

#process.conversion_step = cms.Path(process.allTimingConversions)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.schedule = cms.Schedule(process.ecalTracks_path,
                                process.tree_step)
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
    SkipEvent = cms.untracked.vstring('ProductNotFound'), 
    #wantSummary = cms.untracked.bool(True) 
)

