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
options.register('runera','Run2',VarParsing.multiplicity.singleton,VarParsing.varType.string,'filter to use in event processing');

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
#options.register('globalTag','124X_mcRun3_2022_realistic_postEE_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for sig22 MC');
#options.register('globalTag','124X_mcRun3_2022_realistic_postEE_v1',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','130X_mcRun3_2023_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','133X_mcRun3_2024_realistic_v10',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for QCD MC');
#options.register('globalTag','94X_mc2017_realistic_v14',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gt for GMSB MC');
##options.register('globalTag','112X_mcRun3_2021_realistic_v16',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');

##------------------ data gt  
#options.register('globalTag','106X_dataRun2_v36',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
##options.register('globalTag','106X_dataRun2_v28',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used');
#options.register('globalTag','124X_dataRun3_v15',VarParsing.multiplicity.singleton,VarParsing.varType.string,'gloabl tag to be used 2022');

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
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
#process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 100000

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

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(

        #'root://cmseos.fnal.gov//store/mc/RunIISummer20UL18RECO/QCD_HT2000toInf_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2520000/0184321F-08CA-3740-9700-89D3C169F771.root',

        #SMS_Sig
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/SMS-GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD/250317_234753/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-250_gam_N2ctau-0p1_AOD_41.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoG/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD/250123_145920/0000/SMS-GlGl_mGl-1500_mN2-500_mN1-100_gam_N2ctau-0p1_AOD_99.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD/250607_191347/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_2.root',
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD/250607_191347/0000/SMS-GlGl_mGl-2000_mN2-1900_mN1-200_Zff_N2ctau-0p3_AOD_4.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD/250610_030144/0000/SMS-GlGl_mGl-2000_mN2-400_mN1-200_Zff_N2ctau-0p1_AOD_78.root'
        #'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoZ/SMS-GlGl_mGl-1500_mN2-500_mN1-100_Zll_N2ctau-0p1_AOD/250123_153207/0000/SMS-GlGl_mGl-1500_mN2-500_mN1-100_Zll_N2ctau-0p1_AOD_83.root'

        # WJets
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/6021A8AE-9B41-7A4F-8E33-470F73F3D4FC.root',

        'file:WJets_72B9C618-FE23-1E41-872E-57314D7CB454.root',

        # DiPhotonJetsBox
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/50000/E1C9ABB4-D708-D14E-8E97-803013F321E6.root',
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-80toInf_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2520000/E28F5CDC-017D-2D49-87DA-864981D2E02F.root',

        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18RECO/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1_ext2-v3/40000/72B9C618-FE23-1E41-872E-57314D7CB454.root',

        #GJets
	#'/store/mc/RunIISummer20UL18RECO/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/06CC60F3-D74E-174E-8AD7-DCDAB14FB26F.root',

        #QCD
        #'root://cmsxrootd-site.fnal.gov//store/mc/RunIIFall17DRPremix/QCD_HT200to300_TuneCP5_13TeV-madgraph-pythia8/AODSIM/PU2017_94X_mc2017_realistic_v11-v1/110000/1E29BF8B-5F60-E811-AD1D-0663CE00010C.root',

        # AOD Data MET
        #'file:Met_UL18B_AOD_973EEF0C-44AB-E94A-8591-04DCD00D8B4B.root',
        #'file:/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/test/Met_UL18B_AOD_973EEF0C-44AB-E94A-8591-04DCD00D8B4B.root'      
        #'file:/uscms/home/jaking/nobackup/el8/llpana/CMSSW_13_3_3/src/KUCMSNtupleizer/KUCMSNtupleizer/test/MetPD_003A2484-A2DC-E711-9D0A-02163E019C46.root'

        #'root://cms-xrd-global.cern.ch//store/data/Run2022E/JetMET/AOD/EXODelayedJetMET-27Jun2023-v1/40000/08fd72c7-1be9-4328-9a38-85979d340331.root',

         # AODSIM DPJB model

         #'/store/mc/RunIISummer20UL18RECO/DiPhotonJetsBox_MGG-1000to2000_13TeV-sherpa/AODSIM/106X_upgrade2018_realistic_v11_L1v1-v2/2550000/0B4130F1-3785-8043-957B-E0C2F86E1A50.root',

        ),##<<>>fileNames = cms.untracked.vstring
        #skipEvents=cms.untracked.uint32(300),
)##<<>>process.source = cms.Source("PoolSource",

## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#TT
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
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))#FL

genInfo = True
#genInfo = False
if options.multicrab == True : genInfo = options.hasGenInfo		   

ecalIsoInputsF17 = 'RecoEgamma/ElectronIdentification/data/Fall17/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_94X.txt'
ecalruneraIsoInputsW22 = 'RecoEgamma/ElectronIdentification/data/Run3_Winter22/effAreaElectrons_cone03_pfNeuHadronsAndPhotons_122X.txt'

filterselect = 'None'
#filterselect = 'AL1NpSC'
#filterselect = 'AL1DisSV'
#filterselect = 'SVIPMet100'
#filterselect = 'MET100'
#filterselect = 'AL1IsoPho'
#filterselect = 'IsoPhoMet100'
#filterselect = 'AL1SelEle'
if options.multicrab == True : filterselect = options.eventFilter

#runera = "Run3"  # current siganl model
runera = "Run2" # BG models
if options.multicrab == True : runera = options.runera

dosv = True
#dosv = False
if options.multicrab == True : dosv = options.doSV

############dode = True
dode = False # removed - only false with v30
#if options.multicrab == True : dode = options.doDisEle

doet = True ## with v30 no longer optional when using ( "ecalTracks", "displacedElectronSCs" ) for SCs collection
#doet = False
#if options.multicrab == True : doet = options.doECALTrackOnly
#if doet : 
#    dosv = False
#    dode = False

#probeout = True
probeout = False
if options.multicrab == True : probeout = False

ecalIsoInputs = ecalIsoInputsF17
if "Run3" in runera :  ecalIsoInputs = ecalruneraIsoInputsW22

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

test = cms.vstring( "hltPFMET100", "hltMETClean100", "hltHIPhoton20Eta3p1" )

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

)#triglist0 = cms.vstring(

## Metfilters

process.tree.metFilters = cms.vstring( #"bla" )

    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",    
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_BadChargedCandidateFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"

)##<<>>process.tree.metFilters = cms.vstring(

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

# SVs & ecaltracks aka merged SC collection
process.kuEcalTracks = cms.Sequence( ecalTracks )
process.kuSV = cms.Sequence( muonEnhancedTracks )
#process.kuDisEle = cms.Sequence( displacedElectrons )

#process.kuDisplaced_path = cms.Path()
process.kuDisplaced_path = cms.Path( process.kuEcalTracks )
if ( dosv ) : process.kuDisplaced_path = cms.Path( process.kuEcalTracks + process.kuSV )

## Set final paths and schedule

#############3process.setFlags_path = cms.Path(process.setFlags)
process.tree_step = cms.EndPath(process.tree)
process.endjob_step = cms.EndPath(process.endOfProcess)

#########process.schedule = cms.Schedule( process.kuDisplaced_path, process.setFlags_path, process.tree_step, process.endjob_step )
# Flags must be run as a path in the Schedule in order to appear in the trigger results for KUCMSEventInfo 
# the Flag should return 1 if the event passed the filter and 0 if it failed the filter
process.schedule = cms.Schedule( process.kuDisplaced_path,
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
                                 process.tree_step,
                                 process.endjob_step, )#process.schedule

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

