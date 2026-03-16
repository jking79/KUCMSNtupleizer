import FWCore.ParameterSet.Config as cms

leptonicHYDDRA = cms.PSet(
    # Seeding
    seedCosThetaCut = cms.double(0.75),

    # Shared kinematic cuts
    minMass            = cms.double(2.0),
    minPOverE          = cms.double(0.6),
    maxNormChi2        = cms.double(5.0),
    minDxySignificance = cms.double(25.0),

    # Cleaning (track removal for vertices with > 2 tracks)
    maxCompatibility = cms.double(1.5),
    minCleanCosTheta = cms.double(0.5),
    useDiagonalCut   = cms.bool(False),
    cleanCutSlope    = cms.double(0.0),

    # Final filtering (post-disambiguation, 2-track vertices only)
    minTrackCosTheta             = cms.double(0.5),
    maxTrackCosThetaCM_Limit     = cms.double(0.95),
    maxTrackCosThetaCM_Intercept = cms.double(1.8),
    trackCosThetaCM_Slope        = cms.double(-1.0),
    requireChargeNeutrality      = cms.bool(True),
    minVtxCosTheta               = cms.double(-1.0),
    useAbsVtxCosTheta            = cms.bool(False),

    # Stage enable flags (set to False to skip a stage)
    doMerging        = cms.bool(True),
    doCleaning       = cms.bool(True),
    doDisambiguation = cms.bool(True),
    doFiltering      = cms.bool(True),
)

hadronicHYDDRA = cms.PSet(
    # Seeding
    seedCosThetaCut = cms.double(0.0),

    # Shared kinematic cuts
    minMass            = cms.double(2.0),
    minPOverE          = cms.double(0.6),
    maxNormChi2        = cms.double(5.0),
    minDxySignificance = cms.double(40.0),

    # Hadronic vertex cuts
    minSize       = cms.int32(5),
    minCosTheta   = cms.double(0.0),
    maxDecayAngle = cms.double(0.9),

    # Stage enable flags (set to False to skip a stage)
    doMerging        = cms.bool(True),
    doCleaning       = cms.bool(True),
    doDisambiguation = cms.bool(True),
    doFiltering      = cms.bool(True),
)

hyddraSVs = cms.EDProducer("HyddraSVsProducer",
    tracks       = cms.InputTag("miniAODTrackProducer", "mergedMuonSip2D"),
    pvCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
    leptonic     = leptonicHYDDRA,
    hadronic     = hadronicHYDDRA,
)
