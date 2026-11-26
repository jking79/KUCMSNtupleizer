import FWCore.ParameterSet.Config as cms

miniAODTrackProducer = cms.EDProducer("MiniAODTrackProducer",
    # Packed candidate collections
    pfCandidates = cms.InputTag("packedPFCandidates"),
    lostTracks = cms.InputTag("lostTracks"),
    eleLostTracks = cms.InputTag("lostTracks", "eleTracks"),
    # Muon collections (for global track extraction)
    muons = cms.InputTag("slimmedMuons"),
    displacedMuons = cms.InputTag("slimmedDisplacedMuons"),
    # Primary vertex collection
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    # Cut configuration (disabled by default)
    applyCuts = cms.bool(False),
    minPt = cms.double(1.0),              # pT > 1 GeV
    minAbsSip2D = cms.double(4.0),        # |sip2D| >= 4 (select displaced tracks)
    maxNormalizedChi2 = cms.double(5.0),  # chi2/ndof < 5
    # Repair non-positive-definite MiniAOD pseudo-track covariances in mergedMuonSip2D (HYDDRA input):
    # keep variances, raise correlation-matrix eigenvalues below this floor to it.
    repairCovariances = cms.bool(True),
    minCorrelationEigenvalue = cms.double(3e-2),  # ~median of untouched pseudo-tracks
)
