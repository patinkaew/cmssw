import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from PhysicsTools.NanoAOD.triggerObjects_cff import unpackedPatTrigger, triggerObjectTable, l1bits
from L1Trigger.Configuration.L1TRawToDigi_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger

############################
### Sub Task Definitions ###
############################

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

# other collections are directly dumped, so uncessary to define task

############################
# Scouting Derived Objects #
############################

scoutingPFCandidateTask = cms.Task(scoutingPFCandidate, scoutingPFCandidateTable)
scoutingPFJetReclusterTask = cms.Task(scoutingPFCandidate, scoutingPFJetRecluster, scoutingPFJetReclusterTable)
scoutingCHSJetReclusterTask = cms.Task(scoutingPFCHSCandidate, scoutingCHSJetRecluster, scoutingCHSJetReclusterTable)
scoutingFatCHSJetReclusterTask = cms.Task(scoutingPFCHSCandidate, scoutingFatCHSJetRecluster, scoutingFatCHSJetReclusterTable)

scoutingFatCHSJetReclusterMatchGenExtensionTask = cms.Task(scoutingFatCHSJetReclusterMatchGen, scoutingFatCHSJetReclusterMatchGenExtensionTable)

# TODO: add the remaining derived objects

###################
# Trigger Objects #
###################

## L1 decisions
gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting") 
patTriggerScouting = patTrigger.clone(
        l1tAlgBlkInputTag = "gtStage2DigisScouting", 
        l1tExtBlkInputTag = "gtStage2DigisScouting"
        )

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

# Trig objects
selectedPatTriggerScouting = selectedPatTrigger.clone(src="patTriggerScouting")
slimmedPatTriggerScouting = slimmedPatTrigger.clone(src="selectedPatTriggerScouting")
unpackedPatTriggerScouting = unpackedPatTrigger.clone(patTriggerObjectsStandAlone="slimmedPatTriggerScouting")
scoutingTriggerObjectTable = triggerObjectTable.clone(
        src = "unpackedPatTriggerScouting",
        l1EG = cms.InputTag("gtStage2DigisScouting", "EGamma"),
        l1Sum = cms.InputTag("gtStage2DigisScouting", "EtSum"),
        l1Jet = cms.InputTag("gtStage2DigisScouting", "Jet"),
        l1Muon = cms.InputTag("gtStage2DigisScouting", "Muon"),
        l1Tau = cms.InputTag("gtStage2DigisScouting", "Tau"),
        )

##############################
### Main Tasks Definitions ###
##############################

# default configuration for ScoutingNano
def prepareScoutingNanoTaskCommon():
    # Scouting original objects
    scoutingNanoTaskCommon = cms.Task()
    scoutingNanoTaskCommon.add(scoutingMuonTableTask, scoutingMuonDisplacedVertexTableTask)
    scoutingNanoTaskCommon.add(scoutingElectronTable)
    scoutingNanoTaskCommon.add(scoutingPhotonTable)
    scoutingNanoTaskCommon.add(scoutingTrackTable)
    scoutingNanoTaskCommon.add(scoutingPrimaryVertexTable)
    #scoutingNanoTaskCommon.add(scoutingParticleTable)
    scoutingNanoTaskCommon.add(scoutingPFJetTable)
    scoutingNanoTaskCommon.add(scoutingMETTable, scoutingRhoTable)
    
    # Scouting derived objects
    scoutingNanoTaskCommon.add(scoutingPFCandidateTask)
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingCHSJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterTask)

    return scoutingNanoTaskCommon

def prepareScoutingTriggerTask():
    # add necessary tasks for trigger table
    # also add L1 objects and Trig objects
    scoutingTriggerTask = cms.Task(gtStage2DigisScouting)
    scoutingTriggerTask.add(cms.Task(l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable))
    #scoutingTriggerTask.add(cms.Task(l1bitsScouting, patTriggerScouting))
    #scoutingTriggerTask.add(cms.Task(selectedPatTriggerScouting, slimmedPatTriggerScouting, unpackedPatTriggerScouting, scoutingTriggerObjectTable))
    
    return scoutingTriggerTask

def prepareScoutingNanoTaskMC():
    # additional tasks for running on mc
    scoutingNanoTaskMC = cms.Task()
    scoutingNanoTaskMC.add(scoutingFatCHSJetReclusterMatchGenExtensionTask)
    return cms.Task()

# Common tasks are added to sequence
scoutingNanoTaskCommon = prepareScoutingNanoTaskCommon()
scoutingNanoSequence = cms.Sequence(scoutingNanoTaskCommon)

# Specific tasks will be added to sequence during customization
scoutingTriggerTask = prepareScoutingTriggerTask()
scoutingTriggerSequence = cms.Sequence(L1TRawToDigi+patTriggerScouting+selectedPatTriggerScouting+slimmedPatTriggerScouting+cms.Sequence(scoutingTriggerTask))
scoutingNanoTaskMC = prepareScoutingNanoTaskMC()

def customiseScoutingNanoAOD(process):
    # if running with standard, triggerSequence is already added
    # if running standalone, triggerSequence need to be added
    if not ((hasattr(process, "nanoSequence") and process.schedule.contains(process.nanoSequence))
            or hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
        print("adding trigger")
        process.trigger_step = cms.Path(process.scoutingTriggerSequence)
        process.schedule.extend([process.trigger_step])

    # specifics for running on mc or on data
    runOnMC = hasattr(process,"NANOEDMAODSIMoutput") or hasattr(process,"NANOAODSIMoutput")
    if runOnMC:
        #from PhysicsTools.NanoAOD.globals_cff import puTable
        #process.pileupTask = cms.Task(puTable)
        process.scoutingNanoSequence.associate(scoutingNanoTaskMC)
        # if running with standard object, GEN objects are already added
        # if running standalone, GEN objects need to be added
        if not (hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
            pass
    
    return process

#####################
### Customisation ###
#####################
# these function are designed to be used with --customise flags in cmsDriver.py
# e.g. --customise PhysicsTools/NanoAOD/python/custom_run3scouting_cff.addScoutingPFCandidate

def addScoutingPFCandidate(process):
    process.scoutingNanoTask.add(scoutingPFCandidateTask)
    return process

def addScoutingOriginalWithoutParticleAndTrack(process): 
    process.scoutingNanoTask.add(scoutingMuonTableTask, scoutingMuonDisplacedVertexTableTask)
    process.scoutingNanoTask.add(scoutingElectronTable)
    process.scoutingNanoTask.add(scoutingPhotonTable)
    process.scoutingNanoTask.add(scoutingTrackTable)
    process.scoutingNanoTask.add(scoutingPrimaryVertexTable)
    process.scoutingNanoTask.add(scoutingPFJetTable)
    process.scoutingNanoTask.add(scoutingMETTable, scoutingRhoTable)
    return process

def addScoutingOriginalWithoutParticle(process):
    addScoutingOriginalWithoutParticleAndTrack(process)
    process.scoutingNanoTask.add(scoutingTrackTable)
    return process

def addScoutingOriginal(process):
    addScoutingOriginalWithoutParticle(process)
    process.scoutingNanoTask(scoutingParticleTable)
    return process 

def resetScoutingDefault(process):
    # reset scoutingNanoTask to default
    if hasattr(process, "scoutingNanoTask"):
        process.scoutingNanoTask = prepareScoutingNanoTask()
    return process

def resetScoutingTask(process): 
    # reset all tasks related to scoutingNano configuration
    if hasattr(process, "scoutingNanoTask"):
        process.scoutingNanoTask = cms.Task()
    return process

def addOfflinePFCandidate(process):
    process.scoutingNanoSequence.associate(cms.Task(offlinePFCandidateTable))
    return process

def addOfflinePFJet(process):
    process.scoutingNanoSequence.associate(cms.Task(offlinePFJet, offlinePFJetTable))
    return process
