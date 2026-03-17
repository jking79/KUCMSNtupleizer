import FWCore.ParameterSet.Config as cms

miniAODMuonEnhancedTracks = cms.EDProducer("MiniAODMuonEnhancedTracksProducer",
    generalTracks = cms.InputTag("miniAODTrackProducer", "merged"),
    generalTracksAll = cms.InputTag("miniAODTrackProducer", "mergedAll"),
    slimmedDisplacedMuonBestTracks = cms.InputTag("miniAODTrackProducer", "displacedMuonGlobalTracks"),
    dsaMuonTracks = cms.InputTag("displacedStandAloneMuons"),
    dgMuonTracks = cms.InputTag("displacedGlobalMuons"),
    displacedTracks = cms.InputTag("displacedTracks"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
)
