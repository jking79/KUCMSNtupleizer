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
options.register('runMETFilters',False,VarParsing.multiplicity.singleton,VarParsing.varType.bool,'flag to run Met Filters');

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
        #root://cmseos.fnal.gov/
        #root://cms-xrd-global.cern.ch/
        #file:

        #### AOD Run3   MUST BE IN CMSSW  14 or 15 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        # Data
        #'file:root_files/R23C_JetMet_mini_4a8e9a72-a9af-4448-acbd-bafa6cd684b4.root',
        # MC - Sig
        'root://cmseos.fnal.gov//store/user/lpcsusylep/jaking/KUCMSNtuple/gogoGZ_FullSim_Mini/SMS-GlGl_mGl-2300_mN2-1600_mN1-1000_GZ_N2ctau-0p5_MINI/260203_235322/0000/SMS-GlGl_mGl-2300_mN2-1600_mN1-1000_GZ_N2ctau-0p5_MiniAODv4_99.root',

        ),##<<>>fileNames = cms.untracked.vstring
        secondaryFileNames=cms.untracked.vstring()
        #skipEvents=cms.untracked.uint32(300),
)##<<>>process.source = cms.Source("PoolSource",

## How many events to process
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(2))#ONE
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10))#ST
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(100))#TT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(250))#KT
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(500))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))#KT
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

#filterselect = 'none'
filterselect = 'SVHPMet100'
#filterselect = 'InvMet100IP'
####filterselect = 'AL1NpSC'
####filterselect = 'AL1DisSV'
####3#filterselect = 'SVIPMet100'
#filterselect = 'MET100'
####filterselect = 'AL1IsoPho'
####filterselect = 'IsoPhoMet100'
#####filterselect = 'AL1SelEle'
if options.multicrab == True : filterselect = options.eventFilter

runera = "Run3"  # current siganl model   !!!!!!!!  Run3 must be in CMSSW 14 or 15
#runera = "Run2" # BG models
if options.multicrab == True : runera = options.runera

#dosv = True
dosv = False # with v30 this is the only SV flag that matters, doDisEle and doECALTrackOnly are depreciated
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
if options.multicrab == True : runMETFilters = options.runMETFilters

#probeout = True
probeout = False
if options.multicrab == True : probeout = False

era = "2022"
triggerSet = "PAT"
ecalIsoInputs = ecalIsoInputsF17
if "Run3" in runera :
    if "2022" not in era : triggerSet = "PAT"
    if genInfo : triggerSet = "PAT"
    ecalIsoInputs = ecalruneraIsoInputsW22
    if 'CASTOR' in process.CaloGeometryBuilder.SelectedCalos : 
        process.CaloGeometryBuilder.SelectedCalos.remove('CASTOR') 

if runMETFilters == True : triggerSet = ""  # ?  "SIM" for MC

#genMomChase = True
genMomChase = False
if options.multicrab == True : genMomChase = False
#makeTrigList = True
makeTrigList = False
if options.multicrab == True : makeTrigList = False

print( "Using options : mutlicrab = ",options.multicrab," geninfo = ",genInfo," filter = ",filterselect )
print( "Using options : momChase = ",genMomChase," trgiList = ",makeTrigList," probeout = ",probeout ) 
print( "Using options : doDisEle = (NU) ",dode," doSVs = ",dosv," doECALTrackOnly = (NU) ",doet )
print( "Using options : globalTag = ",options.globalTag," trigger set = ",triggerSet )
print( "With output file name : ",options.outputFileName )

#test = cms.vstring( "hltPFMET100", "hltMETClean100", "hltHIPhoton20Eta3p1" )

process.eventCounter = cms.EDAnalyzer("KUCMSCounter")

# Make the tree 
process.tree = cms.EDAnalyzer("KUCMSNtupilizerMini",

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
                              ##ogGeneralTracks = cms.InputTag("displacedTracks"),
                              ogGsfTracks = cms.InputTag("reducedEgamma"),
                              ecalTracks = cms.InputTag("ecalTracks", "ecalTracks"),
                              tracks = cms.InputTag("ecalTracks", "ecalGeneralTracks"),
                              gsfTracksSrc = cms.InputTag("ecalTracks", "ecalGsfTracks"),
                              displacedSCs = cms.InputTag("ecalTracks", "displacedElectronSCs"),
                              #displacedTracks = cms.InputTag("displacedElectrons", "displacedCandidateTracks"),
                              muonEnhancedTracks = cms.InputTag("muonEnhancedTracks", "muonEnhancedTracks"),
                              combinedMuonTracks = cms.InputTag("muonEnhancedTracks", "combinedMuonTracks"),
                              sip2DMuonEnhancedTracks = cms.InputTag("muonEnhancedTracks", "sip2DMuonEnhancedTracks"),
                              ## vertices
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              timedSVs = cms.InputTag("timedSVs", "timedSecondaryVertices"),
                              ## pfcandidates
                              pfcandidates = cms.InputTag("packedPFCandidates"),
                              #particleflow = cms.InputTag("particleFlow",""),
                              #pfcanphomap = cms.InputTag("particleFlow","photons"),
                              #pfcanootphomap = cms.InputTag("particleFlow","photons"),
                              #pfcanelemap = cms.InputTag("particleFlow","electrons"),
                              ## rho
                              #rho = cms.InputTag("fixedGridRhoFastjetAll"), #fixedGridRhoAll
                              rho = cms.InputTag("fixedGridRhoAll"),
                              ## conversions
                              conversions = cms.InputTag("reducedEgamma","reducedConversions","PAT" ),
                              ## beamSpot
                              beamSpot = cms.InputTag("offlineBeamSpot"),
                              ## trigger
                              #triggerFlagResults = cms.InputTag("TriggerResults","","RECO"),
                              triggerFlagResults = cms.InputTag("TriggerResults","",triggerSet),
                              #triggerFlagResults = cms.InputTag("TriggerResults"),
                              triggerHLTResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent = cms.InputTag(""),
                              ## METs
                              mets = cms.InputTag("slimmedMETs"),
                              ## jets
                              #Charge Hadron Subtracted : charged particles from non-primary vertices (pileup) are removed before clustering.
                              #jets = cms.InputTag("updatedPatJetsUpdatedJEC"),
                              jets = cms.InputTag("slimmedJets"),
                              calojets = cms.InputTag("slimmedCaloJets"),
                              ## electrons
                              electrons = cms.InputTag("slimmedElectrons"),
                              #displacedElectrons = cms.InputTag("displacedElectrons", "displacedElectrons"),
                              #signalDisplacedElectrons = cms.InputTag("displacedElectrons", "signalDisplacedElectrons"),
                              eleMVAIDLooseMap = cms.InputTag("PhotonIDProdGED", "PhotonCutBasedIDLooseEM"),
                              effAreasConfigFile = cms.FileInPath( ecalIsoInputs ),
                              ## muons
                              muons = cms.InputTag("slimmedMuons"),
                              ## photons
                              gedPhotons = cms.InputTag("slimmedPhotons"),
                              phoCBIDLooseMap = cms.InputTag("PhotonIDProd", "PhotonCutBasedIDLooseEM"),
                              ootPhotons = cms.InputTag("slimmedOOTPhotons"),
                              ## ecal recHits
                              recHitsEB = cms.InputTag("reducedEgamma", "reducedEBRecHits"),
                              recHitsEE = cms.InputTag("reducedEgamma", "reducedEERecHits"),
                              ## superclusters
                              #superClusters = cms.InputTag("ecalTracks", "displacedElectronSCs"),
                              superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
                              otherSuperClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
                              ootSuperClusters = cms.InputTag("reducedEgamma", "reducedOOTSuperClusters"),
                              ## caloclusters
                              caloClusters = cms.InputTag("reducedEgamma", "reducedEBEEClusters"),
                              ## gen info
                              genEvt = cms.InputTag("generator", ""),
                              gent0 = cms.InputTag("genParticles", "t0"),
                              genxyz0 = cms.InputTag("genParticles", "xyz0"),
                              pileups = cms.InputTag("slimmedAddPileupInfo", ""),
                              #Phoronpileups = cms.InputTag("mixData", ""),
                              genParticles = cms.InputTag("packedGenParticles", ""),
                              genPartPruned = cms.InputTag("prunedGenParticles",""),
                              genjets = cms.InputTag("slimmedGenJets","")
                              
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

    "Flag_BadChargedCandidateFilter", # 16 - 17 - 18 ?
    "Flag_BadChargedCandidateSummer16Filter",
    "Flag_BadPFMuonDzFilter", # all
    "Flag_BadPFMuonFilter", # all
    "Flag_BadPFMuonSummer16Filter",
    "Flag_CSCTightHalo2015Filter",
    "Flag_CSCTightHaloFilter",
    "Flag_CSCTightHaloTrkMuUnvetoFilter",
    "Flag_EcalDeadCellBoundaryEnergyFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter", # all
    "Flag_HBHENoiseFilter", # 16 - 17 - 18
    "Flag_HBHENoiseIsoFilter", # 16 - 17 - 18
    "Flag_HcalStripHaloFilter",
    "Flag_METFilters",
    "Flag_chargedHadronTrackResolutionFilter",
    "Flag_ecalBadCalibFilter", # all but not good in 22 - 23 - 24
    "Flag_ecalLaserCorrFilter",
    "Flag_eeBadScFilter", # all
    "Flag_globalSuperTightHalo2016Filter", # all
    "Flag_globalTightHalo2016Filter",
    "Flag_goodVertices", # all
    "Flag_hcalLaserEventFilter",
    "Flag_hfNoisyHitsFilter", # all
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
    process.Flag_EcalDeadCellTriggerPrimitiveFilter,
    process.Flag_BadPFMuonFilter,
    process.Flag_BadPFMuonDzFilter,
    process.Flag_hfNoisyHitsFilter,
    process.Flag_eeBadScFilter

]

# Run 2 only except Flag_ecalBadCalibFilter ( Run 2 & 3 ) but is bad in Run 3
metFilterPaths += [

    process.Flag_HBHENoiseFilter,
    process.Flag_HBHENoiseIsoFilter,
    process.Flag_BadChargedCandidateFilter,
    process.Flag_ecalBadCalibFilter

]


# SVs & ecaltracks aka merged SC collection
process.kuEcalTracks = cms.Sequence( ecalTracks )
process.kuSV = cms.Sequence( muonEnhancedTracks )
#process.kuDisEle = cms.Sequence( displacedElectrons )

process.kuDisplaced_path = cms.Path()
#process.kuDisplaced_path = cms.Path( process.kuEcalTracks )
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

