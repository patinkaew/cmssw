import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from L1Trigger.Configuration.L1TRawToDigi_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.NanoAOD.triggerObjects_cff import l1bits
#from PhysicsTools.NanoAOD.triggerObjects_cff import unpackedPatTrigger, triggerObjectTable
#from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
#from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
#from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger

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

# other collections are directly dumped, so unnessary to define tasks

############################
# Scouting Derived Objects #
############################

hltAK4PFCorrectorTask = cms.Task(hltAK4PFFastJetCorrector, hltAK4PFRelativeCorrector, hltAK4PFAbsoluteCorrector, hltAK4PFResidualCorrector, hltAK4PFCorrector)

scoutingPFJetCorrectedTask = cms.Task(scoutingPFJetReco, hltAK4PFCorrectorTask, scoutingPFJetCorrected, scoutingPFJetCorrectedTable)

scoutingPFCandidateTask = cms.Task(scoutingPFCandidate, scoutingPFCandidateTable)
scoutingPFJetReclusterTask = cms.Task(scoutingPFCandidate, scoutingPFJetRecluster, scoutingPFJetReclusterTable)
scoutingPFJetReclusterCorrectionExtensionTask = cms.Task(hltAK4PFCorrectorTask, scoutingPFJetReclusterCorrected, scoutingPFJetReclusterCorrectionExtensionTable)

scoutingCHSJetReclusterTask = cms.Task(scoutingPFCHSCandidate, scoutingCHSJetRecluster, scoutingCHSJetReclusterTable)
scoutingCHSJetReclusterParticleNetTagExtensionTask = cms.Task(scoutingCHSJetReclusterParticleNetJetTagInfos, scoutingCHSJetReclusterParticleNetJetTags, 
        scoutingCHSJetReclusterParticleNetTagExtensionTable)
scoutingCHSJetReclusterMatchGenExtensionTask = cms.Task(scoutingCHSJetReclusterMatchGen, scoutingCHSJetReclusterMatchGenExtensionTable)

scoutingFatCHSJetReclusterTask = cms.Task(scoutingPFCHSCandidate, scoutingFatCHSJetRecluster, scoutingFatCHSJetReclusterTable)
scoutingFatCHSJetReclusterParticleNetTagExtensionTask = cms.Task(scoutingFatCHSJetReclusterParticleNetJetTagInfos, scoutingFatCHSJetReclusterParticleNetJetTags,
        scoutingFatCHSJetReclusterParticleNetTagExtensionTable)
scoutingFatCHSJetReclusterSoftDropMassExtensionTask = cms.Task(scoutingFatCHSJetReclusterSoftDrop, scoutingFatCHSJetReclusterSoftDropMass,
        scoutingFatCHSJetReclusterSoftDropMassExtensionTable)
scoutingFatCHSJetReclusterParticleNetMassExtensionTask = cms.Task(scoutingFatCHSJetReclusterParticleNetJetTagInfos, scoutingFatCHSJetReclusterParticleNetMassRegressionJetTags,
        scoutingFatCHSJetReclusterParticleNetMassExtensionTable)
scoutingFatCHSJetReclusterJetSubstructureVariableExtensionTask = cms.Task(scoutingFatCHSJetReclusterEcfNbeta1, scoutingFatCHSJetReclusterNjettiness,
        scoutingFatCHSJetReclusterJetSubstructureVariableExtensionTable)
scoutingFatCHSJetReclusterMatchGenExtensionTask = cms.Task(scoutingFatCHSJetReclusterMatchGen, scoutingFatCHSJetReclusterMatchGenExtensionTable)


###################
# Trigger Objects #
###################

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

# Trig objects
#patTriggerScouting = patTrigger.clone(
#        l1tAlgBlkInputTag = "gtStage2DigisScouting", 
#        l1tExtBlkInputTag = "gtStage2DigisScouting"
#        )
#selectedPatTriggerScouting = selectedPatTrigger.clone(src="patTriggerScouting")
#slimmedPatTriggerScouting = slimmedPatTrigger.clone(src="selectedPatTriggerScouting")
#unpackedPatTriggerScouting = unpackedPatTrigger.clone(patTriggerObjectsStandAlone="slimmedPatTriggerScouting")
#scoutingTriggerObjectTable = triggerObjectTable.clone(
#        src = "unpackedPatTriggerScouting",
#        l1EG = cms.InputTag("gtStage2DigisScouting", "EGamma"),
#        l1Sum = cms.InputTag("gtStage2DigisScouting", "EtSum"),
#        l1Jet = cms.InputTag("gtStage2DigisScouting", "Jet"),
#        l1Muon = cms.InputTag("gtStage2DigisScouting", "Muon"),
#        l1Tau = cms.InputTag("gtStage2DigisScouting", "Tau"),
#        )

##############################
### Main Tasks Definitions ###
##############################

# default configuration for ScoutingNano
# all Scouting objects are saved except PFCandidate and Track
def prepareScoutingNanoTaskCommon():
    # Scouting original objects
    scoutingNanoTaskCommon = cms.Task()
    scoutingNanoTaskCommon.add(scoutingMuonTableTask, scoutingMuonDisplacedVertexTableTask)
    scoutingNanoTaskCommon.add(scoutingElectronTable)
    scoutingNanoTaskCommon.add(scoutingPhotonTable)
    #scoutingNanoTaskCommon.add(scoutingTrackTable)
    scoutingNanoTaskCommon.add(scoutingPrimaryVertexTable)
    #scoutingNanoTaskCommon.add(scoutingParticleTable)
    #scoutingNanoTaskCommon.add(scoutingPFJetTable)
    scoutingNanoTaskCommon.add(scoutingMETTable, scoutingRhoTable)
    
    # Scouting derived objects
    #scoutingNanoTaskCommon.add(scoutingPFCandidateTask)
    scoutingNanoTaskCommon.add(scoutingPFJetCorrectedTask)
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingPFJetReclusterCorrectionExtensionTask)
    scoutingNanoTaskCommon.add(scoutingCHSJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingCHSJetReclusterParticleNetTagExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterParticleNetTagExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterSoftDropMassExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterParticleNetMassExtensionTask)
    scoutingNanoTaskCommon.add(scoutingFatCHSJetReclusterJetSubstructureVariableExtensionTask)

    return scoutingNanoTaskCommon

def prepareScoutingTriggerTask():
    # add necessary tasks for trigger table
    # also add L1 objects
    scoutingTriggerTask = cms.Task(gtStage2DigisScouting, l1bitsScouting)
    scoutingTriggerTask.add(cms.Task(l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable))
    #scoutingTriggerTask.add(patTriggerScouting)
    #scoutingTriggerTask.add(unpackedPatTriggerScouting)
    #scoutingTriggerTask.add(scoutingTriggerObjectTable)
    #scoutingTriggerTask.add(cms.Task(selectedPatTriggerScouting, slimmedPatTriggerScouting, unpackedPatTriggerScouting))
    #scoutingTriggerTask.add(cms.Task(selectedPatTriggerScouting, slimmedPatTriggerScouting, unpackedPatTriggerScouting, scoutingTriggerObjectTable))
    
    return scoutingTriggerTask

def prepareScoutingNanoTaskMC():
    # additional tasks for running on mc
    scoutingNanoTaskMC = cms.Task()
    scoutingNanoTaskMC.add(scoutingCHSJetReclusterMatchGenExtensionTask)
    scoutingNanoTaskMC.add(scoutingFatCHSJetReclusterMatchGenExtensionTask)
    return scoutingNanoTaskMC

# Common tasks are added to sequence
scoutingNanoTaskCommon = prepareScoutingNanoTaskCommon()
scoutingNanoSequence = cms.Sequence(scoutingNanoTaskCommon)

# Specific tasks will be added to sequence during customization
scoutingTriggerTask = prepareScoutingTriggerTask()
#scoutingTriggerSequence = cms.Sequence(L1TRawToDigi+patTriggerScouting+selectedPatTriggerScouting+slimmedPatTriggerScouting+cms.Sequence(scoutingTriggerTask))
scoutingTriggerSequence = cms.Sequence(L1TRawToDigi+cms.Sequence(scoutingTriggerTask))
scoutingNanoTaskMC = prepareScoutingNanoTaskMC()

def customiseScoutingNanoAOD(process):
    # if running with standard, triggerSequence is already added
    # if running standalone, triggerSequence need to be added
    if not ((hasattr(process, "nanoSequence") and process.schedule.contains(process.nanoSequence))
            or hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
        print("Adding trigger sequence")
        process.trigger_step = cms.Path(process.scoutingTriggerSequence)
        process.schedule.extend([process.trigger_step])

    # specifics for running on mc or on data
    runOnMC = hasattr(process,"NANOEDMAODSIMoutput") or hasattr(process,"NANOAODSIMoutput")
    if runOnMC:
        #from PhysicsTools.NanoAOD.globals_cff import puTable
        process.scoutingNanoSequence.associate(scoutingNanoTaskMC)
        # if running with standard object, GEN objects are already added
        # if running standalone, GEN objects need to be added
        if not (hasattr(process, "nanoSequenceMC") and process.schedule.contains(process.nanoSequenceMC)):
            pass
    
    return process

def customiseScoutingPFNanoAOD(process):
    customiseScoutingNanoAOD(process)
    process.scoutingNanoTask.add(scoutingPFCandidateTask)
    return process

def customiseScoutingFullNanoAOD(process):
    customiseScoutingNanoAOD(process)
    process.scoutingNanoTask.add(scoutingTrackTable)
    process.scoutingNanoTask.add(scoutingPFCandidateTask)
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
