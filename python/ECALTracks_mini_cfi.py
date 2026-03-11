import FWCore.ParameterSet.Config as cms

from KUCMSNtupleizer.KUCMSNtupleizer.TrackAssociator_cfi import tkAssocParamBlock

ecalTracks = cms.EDProducer("ECALTracksProducer",
                            tkAssocParamBlock,
                            generalTracksSrc = cms.InputTag("generalTracks"),
                            gsfElectronTracksSrc = cms.InputTag("electronGsfTracks"),
                            superClusters = cms.InputTag("reducedEgamma", "reducedSuperClusters"),
                            ootSuperClusters = cms.InputTag("reducedEgamma", "reducedOOTSuperClusters"),
                            #superClusters = cms.InputTag("particleFlowEGamma"),
                            #ootSuperClusters = cms.InputTag("particleFlowSuperClusterOOTECAL", "particleFlowSuperClusterOOTECALBarrel"),
)

