import FWCore.ParameterSet.Config as cms

timedSVs = cms.EDProducer('TimedSVsProducer',
                          electronSrc = cms.InputTag("displacedElectrons", "displacedElectrons"),
)

