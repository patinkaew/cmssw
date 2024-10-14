import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from L1Trigger.Configuration.L1TRawToDigi_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits
from PhysicsTools.NanoAOD.globals_cff import puTable

############################
### Sub Task Definitions ###
############################

# Task contains all dependent tasks
# ExtensionTask must be run on top of another Task

#############################
# Scouting Original Objects #
#############################

# Scouting Muon
scoutingMuonTableTask = cms.Task(scoutingMuonTable)
scoutingMuonDisplacedVertexTableTask = cms.Task(scoutingMuonDisplacedVertexTable)

# from 2024, there are two muon collections
from Configuration.Eras.Modifier_run3_scouting_nanoAOD_post2023_cff import run3_scouting_nanoAOD_post2023
run3_scouting_nanoAOD_post2023.toReplaceWith(scoutingMuonTableTask, cms.Task(scoutingMuonVtxTable, scoutingMuonNoVtxTable))\
    .toReplaceWith(scoutingMuonDisplacedVertexTableTask, cms.Task(scoutingMuonVtxDisplacedVertexTable, scoutingMuonNoVtxDisplacedVertexTable))

# other collections are directly from original Run3Scouting objects, so unnessary to define tasks

############################
# Scouting Derived Objects #
############################

scoutingPFCandidateTask = cms.Task(scoutingPFCandidate, scoutingPFCandidateTable)
scoutingPFJetReclusterTask = cms.Task(scoutingPFCandidate, scoutingPFJetRecluster, scoutingPFJetReclusterTable)
scoutingPFJetReclusterParticleNetTagExtensionTask = cms.Task(scoutingPFJetReclusterParticleNetJetTagInfos, scoutingPFJetReclusterParticleNetJetTags, 
        scoutingPFJetReclusterParticleNetTagExtensionTable)
scoutingPFJetReclusterMatchGenExtensionTask = cms.Task(scoutingPFJetReclusterMatchGen, scoutingPFJetReclusterMatchGenExtensionTable)

scoutingFatPFJetReclusterTask = cms.Task(scoutingPFCandidate, scoutingFatPFJetRecluster, scoutingFatPFJetReclusterTable)
scoutingFatPFJetReclusterParticleNetTagExtensionTask = cms.Task(scoutingFatPFJetReclusterParticleNetJetTagInfos, scoutingFatPFJetReclusterParticleNetJetTags,
        scoutingFatPFJetReclusterParticleNetTagExtensionTable)
scoutingFatPFJetReclusterSoftDropMassExtensionTask = cms.Task(scoutingFatPFJetReclusterSoftDrop, scoutingFatPFJetReclusterSoftDropMass,
        scoutingFatPFJetReclusterSoftDropMassExtensionTable)
scoutingFatPFJetReclusterParticleNetMassExtensionTask = cms.Task(scoutingFatPFJetReclusterParticleNetJetTagInfos, scoutingFatPFJetReclusterParticleNetMassRegressionJetTags,
        scoutingFatPFJetReclusterParticleNetMassExtensionTable)
scoutingFatPFJetReclusterJetSubstructureVariableExtensionTask = cms.Task(scoutingFatPFJetReclusterEcfNbeta1, scoutingFatPFJetReclusterNjettiness,
        scoutingFatPFJetReclusterJetSubstructureVariableExtensionTable)
scoutingFatPFJetReclusterMatchGenExtensionTask = cms.Task(scoutingFatPFJetReclusterMatchGen, scoutingFatPFJetReclusterMatchGenExtensionTable)

############################
# Trigger Bits and Objects #
############################

## L1 decisions
gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting") 

## L1 objects
from PhysicsTools.NanoAOD.l1trig_cff import *
l1MuScoutingTable = l1MuTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Muon"))
l1JetScoutingTable = l1JetTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Jet"))
l1EGScoutingTable = l1EGTable.clone(src=cms.InputTag("gtStage2DigisScouting", "EGamma"))
l1TauScoutingTable = l1TauTable.clone(src=cms.InputTag("gtStage2DigisScouting", "Tau"))
l1EtSumScoutingTable = l1EtSumTable.clone(src=cms.InputTag("gtStage2DigisScouting", "EtSum"))

# reduce the variables to the core variables as only these are available in gtStage2Digis
l1EGScoutingTable.variables = cms.PSet(l1EGReducedVars)
l1MuScoutingTable.variables = cms.PSet(l1MuonReducedVars)
l1JetScoutingTable.variables = cms.PSet(l1JetReducedVars)
l1TauScoutingTable.variables = cms.PSet(l1TauReducedVars)
l1EtSumScoutingTable.variables = cms.PSet(l1EtSumReducedVars)

##############################
### Main Tasks Definitions ###
##############################

# default configuration for ScoutingNano common for both data and MC
def prepareScoutingNanoTaskCommon():
    # Scouting original objects
    # all scouting objects are saved except PF Candidate and Track
    scoutingNanoTaskCommon = cms.Task()
    scoutingNanoTaskCommon.add(scoutingMuonTableTask, scoutingMuonDisplacedVertexTableTask)
    scoutingNanoTaskCommon.add(scoutingElectronTable)
    scoutingNanoTaskCommon.add(scoutingPhotonTable)
    scoutingNanoTaskCommon.add(scoutingPrimaryVertexTable)
    scoutingNanoTaskCommon.add(scoutingPFJetTable)
    scoutingNanoTaskCommon.add(scoutingMETTable, scoutingRhoTable)
    
    # Scouting derived objects
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterParticleNetTagExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterParticleNetTagExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterSoftDropMassExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterParticleNetMassExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatPFJetReclusterJetSubstructureVariableExtensionTask)    

    return scoutingNanoTaskCommon

def prepareScoutingTriggerTask():
    scoutingTriggerTask = cms.Task(gtStage2DigisScouting, l1bitsScouting)
    scoutingTriggerTask.add(cms.Task(l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable))

    return scoutingTriggerTask

def prepareScoutingNanoTaskMC():
    # additional tasks for running on mc
    scoutingNanoTaskMC = cms.Task()
    scoutingNanoTaskMC.add(scoutingPFJetReclusterMatchGenExtensionTask)
    scoutingNanoTaskMC.add(scoutingFatPFJetReclusterMatchGenExtensionTask)

    scoutingNanoTaskMC.add(puTable)
    return scoutingNanoTaskMC

# Common tasks added to main scoutingNanoSequence
scoutingNanoTaskCommon = prepareScoutingNanoTaskCommon()
scoutingNanoSequence = cms.Sequence(scoutingNanoTaskCommon)

# Specific tasks which will be added to sequence during customization
scoutingTriggerTask = prepareScoutingTriggerTask()
scoutingTriggerSequence = cms.Sequence(L1TRawToDigi+cms.Sequence(scoutingTriggerTask))
scoutingNanoTaskMC = prepareScoutingNanoTaskMC()

def customiseScoutingNanoAOD(process):
    # if running with standard NanoAOD, triggerSequence is already added
    # if running standalone, triggerSequence need to be added
    if not ((hasattr(process, "nanoSequence") and process.schedule.contains(process.nanoSequence))
            or hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
        process.trigger_step = cms.Path(process.scoutingTriggerSequence)
        process.schedule.extend([process.trigger_step])

    # specific tasks when running on MC
    runOnMC = hasattr(process,"NANOEDMAODSIMoutput") or hasattr(process,"NANOAODSIMoutput")
    if runOnMC:
        process.scoutingNanoSequence.associate(scoutingNanoTaskMC)
    
    return process

#####################
### Customisation ###
#####################
# these function are designed to be used with --customise flags in cmsDriver.py
# e.g. --customise PhysicsTools/NanoAOD/python/custom_run3scouting_cff.addScoutingPFCandidate

def addScoutingParticle(process):
    # original PF candidate without post-processing
    process.scoutingNanoSequence.associate(scoutingParticleTable)
    return process

def addScoutingPFCandidate(process):
    # PF candidate after translation to reco::PFCandidate
    process.scoutingNanoSequence.associate(scoutingPFCandidateTask)
    return process

def addScoutingTrack(process):
    process.scoutingNanoSequence.associate(scoutingTrackTable)
    return process
