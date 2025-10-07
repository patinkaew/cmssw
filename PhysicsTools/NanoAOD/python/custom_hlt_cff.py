import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.hlt_cff import *
from PhysicsTools.NanoAOD.globals_cff import *
from PhysicsTools.NanoAOD.genWeightsTable_cfi import *
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from EventFilter.L1TRawToDigi.caloStage2Digis_cfi import caloStage2Digis
from PhysicsTools.NanoAOD.l1trig_cff import *

# L1 objects
l1JetReducedTable = l1JetTable.clone(variables=cms.PSet(l1JetReducedVars))
l1EtSumReducedTable = l1EtSumTable.clone(variables=cms.PSet(l1JetReducedVars))

# HLTNano should not be used with standard NanoAOD
# So, we will replace nanoSequence, nanoSequenceMC, etc in standard NanoAOD

nanoTableTaskCommon = cms.Task(
    # HLT
    patHLTAK4CaloJetCorrFactor,
    patHLTAK4CaloJet,
    patHLTAK4CaloJetTable,

    patHLTAK8CaloJetCorrFactor,
    patHLTAK8CaloJet,
    patHLTAK8CaloJetTable,

    patHLTAK4PFJetCorrFactor,
    patHLTAK4PFJet,
    patHLTAK4PFJetTable,
    
    patHLTAK8PFJetCorrFactor,
    patHLTAK8PFJet,
    patHLTAK8PFJetTable,

    hltRhoTable,
    hltCaloMETTable,
    hltPFMETTable,
    hltMHTTable,
    hltPixelVertexTable,
    hltVertexPFTable,
    hltDeepInclusiveVertexFinderPFTable,
    hltDeepInclusiveMergedVertexPFTable,

    # L1T
    l1bits,
    gtStage2Digis,
    caloStage2Digis,
    l1JetReducedTable,
    l1EtSumReducedTable,
)

nanoSequenceCommon = cms.Sequence(nanoTableTaskCommon)

nanoSequence = cms.Sequence(nanoSequenceCommon)

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetPartons
from PhysicsTools.PatAlgos.slimming.genParticles_cff import prunedGenParticlesWithStatusOne
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import prunedGenParticles
from PhysicsTools.PatAlgos.slimming.packedGenParticles_cfi import packedGenParticles
from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJets, slimmedGenJetsAK8
from PhysicsTools.PatAlgos.slimming.slimmedGenJetsFlavourInfos_cfi import slimmedGenJetsFlavourInfos
from PhysicsTools.JetMCAlgos.AK4GenJetFlavourInfos_cfi import ak4GenJetFlavourInfos
from PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi import selectedHadronsAndPartonsForGenJetsFlavourInfos
from PhysicsTools.NanoAOD.genparticles_cff import *
from PhysicsTools.NanoAOD.jetMC_cff import *

nanoTableTaskFS = cms.Task(
    prunedGenParticlesWithStatusOne,
    prunedGenParticles,
    genParticleTask,
    genParticleTable,
    
    packedGenParticles,
    slimmedGenJets,
    genJetTable,
    patJetPartonsNano,
    genJetFlavourAssociation,
    selectedHadronsAndPartonsForGenJetsFlavourInfos,
    ak4GenJetFlavourInfos,
    slimmedGenJetsFlavourInfos,
    genJetFlavourTable,

    #slimmedGenJetsAK8,
    patJetPartons,

    patHLTAK4PFJetGenJetMatch,
    patHLTAK4PFJetFlavourAssociation,
    patHLTAK4PFJetMCTable,

    # older tables
    #genJetTable,
    #genJetNoNuTable,
    #genJetAK8Table,
    #genJetAK8NoNuTable,
    #genWeightsTableTask,
)

#nanoSequenceFS = cms.Sequence(NanoGenTable+cms.Sequence(nanoTableTaskFS))
nanoSequenceFS = cms.Sequence(nanoTableTaskFS)

#nanoSequenceMC = cms.Sequence(nanoSequenceCommon+cms.Sequence(nanoTableTaskFS))
nanoSequenceMC = cms.Sequence(nanoSequenceFS+nanoSequenceCommon)

def nanoAOD_customizeCommon(process):
    #process.load("PatAlgos.slimming.slimmedGenJets_cfi") # load slimmedGenJets and 
    return process
