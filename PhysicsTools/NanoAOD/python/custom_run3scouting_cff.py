import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.run3scouting_cff import *
from PhysicsTools.NanoAOD.globals_cff import puTable
from PhysicsTools.NanoAOD.triggerObjects_cff import unpackedPatTrigger, triggerObjectTable, l1bits
from L1Trigger.Configuration.L1TRawToDigi_cff import *
from EventFilter.L1TRawToDigi.gtStage2Digis_cfi import gtStage2Digis
from PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cfi import patTrigger
from PhysicsTools.PatAlgos.slimming.selectedPatTrigger_cfi import selectedPatTrigger
from PhysicsTools.PatAlgos.slimming.slimmedPatTrigger_cfi import slimmedPatTrigger

#scoutingNanoTask = cms.Task()
nanoTableTaskCommon = cms.Task()
addScoutingMuonTable = True
addScoutingElectronTable = True
addScoutingPhotonTable = True



gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
nanoTableTaskCommon.add(gtStage2DigisScouting)

from PhysicsTools.NanoAOD.l1trig_cff import *
# L1 Objects from L1Nano
if (True or addL1Objects):
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

    nanoTableTaskCommon.add(l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable)


# Scouting Original Objects
def addScoutingMuon(scoutingNanoTask, addDisplacedVertex=True):
    before2024 = False
    if (before2024):
        scoutingNanoTask.add(scoutingMuonArrayTable)
        scoutingNanoTask.add(scoutingMuonTable)
        if (addDisplacedVertex):
            scoutingNanoTask.add(scoutingDisplacedVertexTable)
    else:
       # scoutingMuonVtxArrayTable = scoutingMuonArrayTable.clone(
       #         src = cms.InputTag("hltScoutingMuonPackerVtx"),
       #         vertex_index_collection_name = cms.string("ScoutingMuonVtxVertexIndex"),
       #         hit_pattern_collection_name = cms.string("ScoutingMuonVtxTrackHitPattern"),
       #         )
       # scoutingMuonVtxTable = scoutingMuonTable.clone(
       #         src = cms.InputTag("hltScoutingMuonPackerVtx"),
       #         name = cms.string("ScoutingMuonVtx"),
       #         doc  = cms.string("Scouting Muon Vtx information"),
       #         )
        
        # change parameters to Muon with Vtx
        # scoutingMuonVtxArrayTable = scoutingMuonArrayTable
        scoutingMuonArrayTable.src = cms.InputTag("hltScoutingMuonPackerVtx:")
        scoutingMuonArrayTable.vertex_index_collection_name = cms.string("ScoutingMuonVtxVertexIndex")
        scoutingMuonArrayTable.hit_pattern_collection_name = cms.string("ScoutingMuonVtxTrackHitPattern")

        #scoutingMuonTable = scoutingMuonTable
        scoutingMuonTable.src = cms.InputTag("hltScoutingMuonPackerVtx")
        scoutingMuonTable.name = cms.string("ScoutingMuonVtx")
        scoutingMuonTable.doc = cms.string("Scouting Muon Vtx information")


        scoutingNanoTask.add(scoutingMuonArrayTable)
        scoutingNanoTask.add(scoutingMuonTable)
        if (addDisplacedVertex):
            #scoutingDisplacedVertexVtxTable = scoutingDisplacedVertexTable.clone(
            #        src = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx"),
            #        name = cms.string("ScoutingMuonVtxDisplacedVertex"),
            #        doc  = cms.string("Scouting Muon Vtx DisplacedVertex information"),
            #        )
            #scoutingDisplacedVertexVtxTable = scoutingDisplacedVertexTable
            scoutingDisplacedVertexTable.src = cms.InputTag("hltScoutingMuonPackerVtx", "displacedVtx")
            scoutingDisplacedVertexTable.name = cms.string("ScoutingMuonVtxDisplacedVertex")
            scoutingDisplacedVertexTable.doc = cms.string("Scouting Muon Vtx DisplacedVertex information")
            scoutingNanoTask.add(scoutingDisplacedVertexTable)
        
        # clone and change parameters to Muon with NoVtx
        scoutingMuonNoVtxArrayTable = scoutingMuonArrayTable.clone(
                src = cms.InputTag("hltScoutingMuonPackerPackerNoVtx"),
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

        if (addDisplacedVertex):
            scoutingDisplacedVertexNoVtxTable = scoutingDisplacedVertexTable.clone(
                    src = cms.InputTag("hltScoutingMuonPackerNoVtx", "displacedVtx"),
                    name = cms.string("ScoutingMuonNoVtxDisplacedVertex"),
                    doc  = cms.string("Scouting Muon NoVtx DisplacedVertex information"),
                    )
            scoutingNanoTask.add(scoutingDisplacedVertexNoVtxTable)

    return scoutingNanoTask

def addScoutingElectron(scoutingNanoTask):
    scoutingNanoTask.add(scoutingElectronArrayTable)
    scoutingNanoTask.add(scoutingElectronTable)
    return scoutingNanoTask

def addScoutingPhoton(scoutingNanoTask):
    scoutingNanoTask.add(scoutingPhotonTable)
    return scoutingNanoTask

def addScoutingEGamma(scoutingNanoTask):
    scoutingNanoTask = addScoutingElectron(scoutingNanoTask)
    scoutingNanoTask = addScoutingPhoton(scoutingNanoTask)
    return scoutingNanoTask

def addScoutingParticleTable(scoutingNanoTask):
    scoutingNanoTask.add(scoutingParticleTable)
    return scoutingNanoTask

def addScoutingJetMET(scoutingNanoTask):
    scoutingNanoTask.add(scoutingPFJetTable)
    scoutingNanoTask.add(scoutingMETTable)
    scoutingNanoTask.add(scoutingRhoTable)
    return scoutingNanoTask

def addScoutingTrack(scoutingNanoTask):
    scoutingNanoTask.add(scoutingTrackTable)
    return scoutingNanoTask

def addScoutingPrimaryVertex(scoutingNanoTask):
    scoutingNanoTask.add(scoutingPrimaryVertexTable)
    return scoutingNanoTask

def addOriginalScoutingTables(scoutingNanoTask, addScoutingParticle=True):
    scoutingNanoTask = addScoutingMuon(scoutingNanoTask)
    scoutingNanoTask = addScoutingEGamma(scoutingNanoTask)
    if addScoutingParticle:
        scoutingNanoTask = addScoutingParticleTable(scoutingNanoTask)
    scoutingNanoTask = addScoutingJetMET(scoutingNanoTask)
    scoutingNanoTask = addScoutingTrack(scoutingNanoTask)
    scoutingNanoTask = addScoutingPrimaryVertex(scoutingNanoTask)
    return scoutingNanoTask

####################################
##### Scouting Derived Objects ##### 
####################################

from PhysicsTools.NanoAOD.l1trig_cff import *
def addL1Tables(scoutingNanoTask, addL1Objects=True, addTrigObjTable=False):
    # gtStage2Digis producer
    global gtStage2Digis
    gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")
    scoutingNanoTask.add(gtStage2DigisScouting)

    # L1 Objects from L1Nano
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
    
    # TrigObjTable similar to standard nano
    if (addTrigObjTable):
        l1bitsScouting = l1bits.clone(src="gtStage2DigisScouting")
        patTriggerScouting = patTrigger.clone(l1tAlgBlkInputTag="gtStage2DigisScouting",l1tExtBlkInputTag="gtStage2DigisScouting")

        selectedPatTriggerScouting = selectedPatTrigger.clone(src="patTriggerScouting")
        slimmedPatTriggerScouting = slimmedPatTrigger.clone(src="selectedPatTriggerScouting")
        unpackedPatTriggerScouting = unpackedPatTrigger.clone(patTriggerObjectsStandAlone="slimmedPatTriggerScouting")
        triggerObjectTableScouting = triggerObjectTable.clone(src="unpackedPatTriggerScouting")

    return scoutingNanoTask

def addScoutingJetMETReclusterAK4(scoutingNanoTask):
    return scoutingNanoTask

def addScoutingJetMETReclusterAK8(scoutingNanoTask):
    return scoutingNanoTask

def addMCTables(scoutingNanoTask):
    return scoutingNanoTask


# prepare nano
#scoutingNanoTask = cms.Task()
#scoutingNanoTask = addOriginalScoutingTables(scoutingNanoTask, addScoutingParticle=True)
#scoutingNanoTask = addL1Tables(scoutingNanoTask, addL1Objects=True, addTrigObjTable=False)

#nanoTableTaskCommon = cms.Task(scoutingNanoTask)
#nanoTableTaskCommon = cms.Task()
#nanoTableTaskCommon = addScoutingJetMET(nanoTableTaskCommon)
#nanoTableTaskCommon = addScoutingParticleTable(nanoTableTaskCommon)
#nanoTableTaskCommon = addScoutingMuon(nanoTableTaskCommon)
#nanoTableTaskCommon = addL1Tables(nanoTableTaskCommon)
#nanoTableTaskCommon = addOriginalScoutingTables(nanoTableTaskCommon, addScoutingParticle=True)
nanoSequenceCommon = cms.Sequence(nanoTableTaskCommon)

nanoSequence = cms.Sequence(nanoSequenceCommon)

nanoSequenceMC = cms.Sequence(nanoSequenceCommon)

def nanoAOD_customizeCommon(process):
    #scoutingNanoTask = cms.Task()
    #process = addOriginalScoutingTables(process, addScoutingParticle=True)
    #process = addL1Tables(process, addL1Objects=True, addTrigObjTable=False)
    
    #process.nanoSequence = cms.Sequence(scoutingNanoTask)
    #process.nanoSequenceMC = cms.Sequence(process.nanoSequence + cms.Sequence(scouting))
    
    return process


#particleTask = cms.Task(scoutingPFCands)
#particleTableTask = cms.Task(particleScoutingTable)
#ak4JetTableTask = cms.Task(ak4ScoutingJets,ak4ScoutingJetParticleNetJetTagInfos,ak4ScoutingJetParticleNetJetTags,ak4ScoutingJetTable)
#ak8JetTableTask = cms.Task(ak8ScoutingJets,ak8ScoutingJetsSoftDrop,ak8ScoutingJetsSoftDropMass,ak8ScoutingJetEcfNbeta1,ak8ScoutingJetNjettiness,ak8ScoutingJetParticleNetJetTagInfos,ak8ScoutingJetParticleNetJetTags,ak8ScoutingJetParticleNetMassRegressionJetTags,ak8ScoutingJetTable)

## L1 decisions
#gtStage2DigisScouting = gtStage2Digis.clone(InputLabel="hltFEDSelectorL1")

#triggerTask = cms.Task(
#    gtStage2DigisScouting, l1MuScoutingTable, l1EGScoutingTable, l1TauScoutingTable, l1JetScoutingTable, l1EtSumScoutingTable, 
#    unpackedPatTriggerScouting,triggerObjectTableScouting,l1bitsScouting
#)
#triggerSequence = cms.Sequence(L1TRawToDigi+patTriggerScouting+selectedPatTriggerScouting+slimmedPatTriggerScouting+cms.Sequence(triggerTask))

# MC tasks
#genJetTask = cms.Task(ak4ScoutingJetMatchGen,ak4ScoutingJetExtTable,ak8ScoutingJetMatchGen,ak8ScoutingJetExtTable)
#puTask = cms.Task(puTable)

#nanoTableTaskCommon = cms.Task(photonScoutingTable,muonScoutingTable,electronScoutingTable,trackScoutingTable,primaryvertexScoutingTable,displacedvertexScoutingTable,rhoScoutingTable,metScoutingTable,#particleTask,particleTableTask,ak4JetTableTask,ak8JetTableTask)

#nanoSequenceCommon = cms.Sequence(triggerSequence,nanoTableTaskCommon)

#nanoSequence = cms.Sequence(nanoSequenceCommon)

#nanoSequenceMC = cms.Sequence(nanoSequenceCommon + cms.Sequence(cms.Task(genJetTask,puTask)))
