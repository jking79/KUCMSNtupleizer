import FWCore.ParameterSet.Config as cms

muonEnhancedTracks = cms.EDProducer("MuonEnhancedTracksProducer",
                                    generalTracks = cms.InputTag("generalTracks"),
                                    dsaMuonTracks = cms.InputTag("displacedStandAloneMuons"),
                                    dgMuonTracks = cms.InputTag("displacedGlobalMuons"),
                                    displacedTracks = cms.InputTag("displacedTracks"),
                                    displacedMuons = cms.InputTag("displacedMuons"),
                                    pvCollection = cms.InputTag("offlinePrimaryVertices"),
)

