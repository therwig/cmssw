import FWCore.ParameterSet.Config as cms

L1MetPFHWProducer = cms.EDProducer("L1MetPFHWProducer",
                                 src = cms.InputTag("L1PFProducer","l1pfCandidates"),
                                 maxCands = cms.uint32(128)
                             )
