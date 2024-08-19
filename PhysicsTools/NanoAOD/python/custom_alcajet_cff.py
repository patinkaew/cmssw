import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.alcajet_cff import *

# nanoSequenceMC = cms.Sequence(nanoSequenceCommon + cms.Sequence(cms.Task(genJetTask,puTask)))
HLTPreprocessTask = cms.Task(HLTAK4PFJets, HLTAK8PFJets, HLTAK4PFJetsCorrectedMatchedToCaloJets10, HLTAK8PFJetsCorrectedMatchedToCaloJets10, HLTFixedGridRhoFastjetAll)

HLTJetforJECTask = cms.Task(HLTAK4JetForJEC, HLTAK8JetForJEC)

#nanoTableTaskCommon = cms.Task(HLTAK4PFJetTable, HLTAK4PFJetsCorrectedMatchedToCaloJets10Table, HLTRhoTable, HLTRawAK4PFMET, HLTRawAK4PFMETTable, HLTJetforJECTask, HLTAK4PFJetCorrectedMatchedToCaloJets10ForJECTable, HLTAK4PFJetForJECTable) #HLTAK8PFJetsCorrectedMatchedToCaloJets10Table
nanoTableTaskCommon = cms.Task(HLTPreprocessTask, HLTRhoTable, HLTRawAK4PFMET, HLTRawAK4PFMETTable, HLTJetforJECTask, HLTAK4PFJetCorrectedMatchedToCaloJets10ForJECTable, HLTAK8PFJetCorrectedMatchedToCaloJets10ForJECTable)

nanoSequenceCommon = cms.Sequence(nanoTableTaskCommon)

nanoSequence = cms.Sequence(nanoSequenceCommon)

#nanoSequenceMC = cms.Sequence(nanoSequenceCommon)

def nanoAOD_customizeCommon(process):
    return process
