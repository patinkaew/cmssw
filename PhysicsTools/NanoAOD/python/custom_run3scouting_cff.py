import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from PhysicsTools.NanoAOD.globals_cff import puTable
from PhysicsTools.NanoAOD.triggerObjects_cff import unpackedPatTrigger, triggerObjectTable, l1bits
from L1Trigger.Configuration.L1TRawToDigi_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger

################################################
##### Configure collections to be included #####
################################################

# Scouting Original Objects
addScoutingMuon = False
addScoutingMuonVtx = True
addScoutingMuonNoVtx = True
addScoutingDisplacedVertex = True
addScoutingElectron = True
addScoutingPhoton = True
addScoutingTrack = True
addScoutingPrimaryVertex = True
addScoutingParticle = True
addScoutingPFJet = True
addScoutingMET = True
addScoutingRho = True

# Scouting Extra Objects
addL1Objects = True
addTrigObj = False
addScoutingPFJetRecluster = True
addScoutingRawPFMETRecluster = True
addScoutingCHSJetRecluster = True
addScoutingRawCHSMETRecluster = True

##################################################
##### Add tasks according to set flags above #####
##################################################
scoutingNanoTask = cms.Task()

#####################################
##### Scouting Original Objects #####
#####################################

#################
# Scouting Muon #
#################

# before 2024, there is only one collection of scouting muon
if (addScoutingMuon and not (addScoutingMuonVtx or addScoutingMuonNoVtx)):
    scoutingNanoTask.add(scoutingMuonArrayTable, scoutingMuonTable)
    if (addScoutingDisplacedVertex):
        scoutingNanoTask.add(scoutingDisplacedVertexTable)

# from 2024, there are two versions of scouting muon
if (addScoutingMuonVtx):
    scoutingMuonVtxArrayTable = scoutingMuonArrayTable.clone(
            src = cms.InputTag("hltScoutingMuonPackerVtx"),
            vertex_index_collection_name = cms.string("ScoutingMuonVtxVertexIndex"),
            hit_pattern_collection_name = cms.string("ScoutingMuonVtxTrackHitPattern"),
            )
    scoutingMuonVtxTable = scoutingMuonTable.clone(
            src = cms.InputTag("hltScoutingMuonPackerVtx"),
            name = cms.string("ScoutingMuonVtx"),
            doc  = cms.string("Scouting Muon Vtx information"),
            )

    scoutingNanoTask.add(scoutingMuonVtxArrayTable)
    scoutingNanoTask.add(scoutingMuonVtxTable)
    if (addScoutingDisplacedVertex):
        scoutingDisplacedVertexVtxTable = scoutingDisplacedVertexTable.clone(
                src = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx"),
                name = cms.string("ScoutingMuonVtxDisplacedVertex"),
                doc  = cms.string("Scouting Muon Vtx DisplacedVertex information"),
                )
        scoutingNanoTask.add(scoutingDisplacedVertexVtxTable)

if (addScoutingMuonNoVtx):
    scoutingMuonNoVtxArrayTable = scoutingMuonArrayTable.clone(
            src = cms.InputTag("hltScoutingMuonPackerNoVtx"),
            vertex_index_collection_name = cms.string("ScoutingMuonNoVtxVertexIndex"),
            hit_pattern_collection_name = cms.string("ScoutingMuonNoVtxTrackHitPattern"),
            )
    scoutingMuonNoVtxTable = scoutingMuonTable.clone(
            src = cms.InputTag("hltScoutingMuonPackerNoVtx"),
            name = cms.string("ScoutingMuonNoVtx"),
            doc  = cms.string("Scouting Muon NoVtx information"),
            )

    scoutingNanoTask.add(scoutingMuonNoVtxArrayTable)
    scoutingNanoTask.add(scoutingMuonNoVtxTable)

    if (addScoutingDisplacedVertex):
        scoutingDisplacedVertexNoVtxTable = scoutingDisplacedVertexTable.clone(
                src = cms.InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx"),
                    name = cms.string("ScoutingMuonNoVtxDisplacedVertex"),
                    doc  = cms.string("Scouting Muon NoVtx DisplacedVertex information"),
                    )
        scoutingNanoTask.add(scoutingDisplacedVertexNoVtxTable)

###################
# Scouting EGamma #
###################

if (addScoutingElectron):
    scoutingNanoTask.add(scoutingElectronArrayTable)
    scoutingNanoTask.add(scoutingElectronTable)

if (addScoutingPhoton):
    scoutingNanoTask.add(scoutingPhotonTable)

#####################################
# Scouting Track and Primary Vertex #
#####################################

if (addScoutingTrack):
    scoutingNanoTask.add(scoutingTrackTable)

if (addScoutingPrimaryVertex):
    scoutingNanoTask.add(scoutingPrimaryVertexTable)

#####################
# Scouting Particle #
#####################

if (addScoutingParticle):
    scoutingNanoTask.add(scoutingParticleTable)

###################
# Scouting JetMET #
###################

if (addScoutingPFJet):
    scoutingNanoTask.add(scoutingPFJetTable)

if (addScoutingMET):
    scoutingNanoTask.add(scoutingMETTable)

if (addScoutingRho):
    scoutingNanoTask.add(scoutingRhoTable)

##################################
##### Scouting Extra Objects #####
##################################

##################
# L1 Information #
##################

# clone gtStage2Digis for scouting
if (addL1Objects or addTrigObj):
    gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
    scoutingNanoTask.add(gtStage2DigisScouting)

# L1 Objects from L1Nano
from PhysicsTools.NanoAOD.l1trig_cff import *
if (addL1Objects):
    l1MuScoutingTable = l1MuTable.clone(src=cms.InputTag("gtStage2DigisScouting","Muon"))
    l1JetScoutingTable = l1JetTable.clone(src=cms.InputTag("gtStage2DigisScouting","Jet"))
    l1EGScoutingTable = l1EGTable.clone(src=cms.InputTag("gtStage2DigisScouting","EGamma"))
    l1TauScoutingTable = l1TauTable.clone(src=cms.InputTag("gtStage2DigisScouting","Tau"))
    l1EtSumScoutingTable = l1EtSumTable.clone(src=cms.InputTag("gtStage2DigisScouting","EtSum"))

    #reduce the variables to the core variables as only these are available in gtStage2Digis
    l1EGScoutingTable.variables = cms.PSet(l1EGReducedVars)
    l1MuScoutingTable.variables = cms.PSet(l1MuonReducedVars)
    l1JetScoutingTable.variables = cms.PSet(l1JetReducedVars)
    l1TauScoutingTable.variables = cms.PSet(l1TauReducedVars)
    l1EtSumScoutingTable.variables = cms.PSet(l1EtSumReducedVars)

    scoutingNanoTask.add(l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable)

# TrigObjTable similar to standard nano TODO:fix this
if (addTrigObj):
    l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting")
    patTriggerScouting = patTrigger.clone(l1tAlgBlkInputTag="gtStage2DigisScouting",l1tExtBlkInputTag="gtStage2DigisScouting")

    selectedPatTriggerScouting = selectedPatTrigger.clone(src="patTriggerScouting")
    slimmedPatTriggerScouting = slimmedPatTrigger.clone(src="selectedPatTriggerScouting")
    unpackedPatTriggerScouting = unpackedPatTrigger.clone(patTriggerObjectsStandAlone="slimmedPatTriggerScouting")
    triggerObjectTableScouting = triggerObjectTable.clone(src="unpackedPatTriggerScouting")

#################
# Jet Recluster #
#################

if (addScoutingPFJetRecluster or addScoutingPFMETRecluster):
    scoutingNanoTask.add(scoutingPFCandidate)

if (addScoutingPFJetRecluster):
    scoutingNanoTask.add(scoutingPFJetRecluster, scoutingPFJetReclusterTable)

if (addScoutingRawPFMETRecluster):
    scoutingNanoTask.add(scoutingRawPFMETRecluster, scoutingRawPFMETRecluster)

if (addScoutingCHSJetRecluster or addScoutingCHSMETRecluster):
    scoutingNanoTask.add(scoutingPFCHSCandidate)

if (addScoutingCHSJetRecluster):
    scoutingNanoTask.add(scoutingCHSJetRecluster, scoutingCHSJetReclusterTable)

if (addScoutingRawCHSMETRecluster):
    scoutingNanoTask.add(scoutingRawCHSMETRecluster, scoutingRawCHSMETReclusterTable)

###############
# Jet tagging #
###############

# TODO

########################################
##### Specific tasks for Data Only #####
########################################

# TODO

######################################
##### Specific tasks for MC Only #####
######################################

# TODO


nanoSequenceCommon = cms.Sequence(scoutingNanoTask)

nanoSequence = cms.Sequence(nanoSequenceCommon)

nanoSequenceMC = cms.Sequence(nanoSequenceCommon)

def nanoAOD_customizeCommon(process):
    return process
