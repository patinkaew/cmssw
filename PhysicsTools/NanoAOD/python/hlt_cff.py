import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.PatAlgos.recoLayer0.jetCorrFactors_cfi import patJetCorrFactors
from PhysicsTools.PatAlgos.producersLayer1.jetProducer_cfi import _patJets
from PhysicsTools.PatAlgos.producersLayer1.jetUpdater_cfi import updatedPatJets
from PhysicsTools.NanoAOD.simplePATJetFlatTableProducer_cfi import simplePATJetFlatTableProducer

##############
# AK4CaloJet #
##############

# Jet Energy Correction
patHLTAK4CaloJetCorrFactor = patJetCorrFactors.clone(
    src = "hltAK4CaloJets",
    payload = "AK4CaloHLT",
    levels = ["L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"],
    primaryVertices = "hltPixelVertices",
    rho = "hltFixedGridRhoFastjetAllCalo",
)

# convert to PAT
patHLTAK4CaloJet = _patJets.clone(
    jetSource            = "hltAK4CaloJets",
    addJetCorrFactors    = True,
    jetCorrFactorsSource = ["patHLTAK4CaloJetCorrFactor"],
    addBTagInfo          = False,
    addDiscriminators    = False,
    discriminatorSources = [],
    addAssociatedTracks  = False,
    addJetCharge         = False,
    addGenPartonMatch    = False,
    embedGenPartonMatch  = False,
    addGenJetMatch       = False,
    getJetMCFlavour      = False,
    addJetFlavourInfo    = False,
)

# output to table
patHLTAK4CaloJetTable = simplePATJetFlatTableProducer.clone(
    src = cms.InputTag("patHLTAK4CaloJet"),
    name = cms.string("HLTAK4CaloJet"),
    doc = cms.string("from hltAK4CaloJets"),
    cut = cms.string(""),
    variables = cms.PSet(
        P4Vars,
        area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
        emEF = Var("emEnergyFraction()", float, doc="electromagnetic energy fraction", precision=10),
        hEF = Var("energyFractionHadronic()", float, doc="hadronic energy fraction", precision=10),
        rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
    )
)

##############
# AK8CaloJet #
##############

# Jet Energy Correction
patHLTAK8CaloJetCorrFactor = patHLTAK4CaloJetCorrFactor.clone(
    src = "hltAK8CaloJets",
    payload = "AK8CaloHLT",
)

# convert to PAT
patHLTAK8CaloJet = patHLTAK4CaloJet.clone(
    jetSource = "hltAK8CaloJets",
    jetCorrFactorsSource = ["patHLTAK8CaloJetCorrFactor"],
)

# output to table
patHLTAK8CaloJetTable = patHLTAK4CaloJetTable.clone(
    src = cms.InputTag("patHLTAK8CaloJet"),
    name = cms.string("HLTAK8CaloJet"),
    doc = cms.string("from hltAK8CaloJets"),
)

############
# AK4PFJet #
############

# Jet Energy Correction
patHLTAK4PFJetCorrFactor = patJetCorrFactors.clone(
    src = "hltAK4PFJets",
    payload = "AK4PFHLT",
    levels = ["L1FastJet", "L2Relative", "L3Absolute", "L2L3Residual"],
    primaryVertices = "hltPixelVertices",
    rho = "hltFixedGridRhoFastjetAll"
)

# convert to PAT
patHLTAK4PFJet = patHLTAK4CaloJet.clone(
    jetSource = "hltAK4PFJets",
    jetCorrFactorsSource = ["patHLTAK4PFJetCorrFactor"],
    addBTagInfo = True,
    addDiscriminators = True,
    matchJetTagByDeltaR = True,
    maxJetTagDeltaR = cms.double(0.0001),
    discriminatorSources = ["hltParticleNetONNXJetTags:probtauhp",
                            "hltParticleNetONNXJetTags:probtauhm",
                            "hltParticleNetONNXJetTags:probb",
                            "hltParticleNetONNXJetTags:probc",
                            "hltParticleNetONNXJetTags:probuds",
                            "hltParticleNetONNXJetTags:probg"],
    addGenJetMatch = True,
    embedGenJetMatch = True,
    genJetMatch = cms.InputTag("patHLTAK4PFJetGenJetMatch"),
    getJetMCFlavour = True,
    useLegacyJetMCFlavour = False,
    addJetFlavourInfo = True,
    JetFlavourInfoSource = cms.InputTag("patHLTAK4PFJetFlavourAssociation"),
)

# common variables for both AK4PFJet and AK8PFJet
PFJetVariables = cms.PSet(
    P4Vars,
    area = Var("jetArea()", float, doc="jet catchment area, for JECs",precision=10),
    chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision=10),
    neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision=10),
    chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision=10),
    neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision=10),
    hfHEF = Var("HFHadronEnergyFraction()",float,doc="hadronic Energy Fraction in HF",precision=10),
    hfEmEF = Var("HFEMEnergyFraction()",float,doc="electromagnetic Energy Fraction in HF",precision=10),
    muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision=10),
    chHadMultiplicity = Var("chargedHadronMultiplicity()", "int16", doc="number of charged hadrons in the jet"),
    neHadMultiplicity = Var("neutralHadronMultiplicity()", int, doc="number of neutral hadrons in the jet"),
    hfHadMultiplicity = Var("HFHadronMultiplicity()", int, doc="number of HF hadrons in the jet"),
    hfEMMultiplicity = Var("HFEMMultiplicity()", int, doc="number of HF EMs in the jet"),
    muMultiplicity = Var("muonMultiplicity()", int, doc="number of muons in the jet"),
    elMultiplicity = Var("electronMultiplicity()", int, doc="number of electrons in the jet"),
    phMultiplicity = Var("photonMultiplicity()", int, doc="number of photons in the jet"),
    nConstituents = Var("numberOfDaughters()", int, doc="number of particles in the jet"),
)

# output to table
patHLTAK4PFJetTable = simplePATJetFlatTableProducer.clone(
    src = cms.InputTag("patHLTAK4PFJet"),
    name = cms.string("HLTAK4PFJet"),
    doc = cms.string("hltAK4PFJet"),
    cut = cms.string(""),
    variables = cms.PSet(
        PFJetVariables,    
        rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
        hltPNet_probtauhp = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probtauhp'):-1", float, doc="HLT PNet tagger tauhp raw score", precision=12),
        hltPNet_probtauhm = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probtauhm'):-1", float, doc="HLT PNet tagger tauhm raw score", precision=12),
        hltPNet_probb = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probb'):-1", float, doc="HLT PNet tagger b raw score", precision=12),
        hltPNet_probc = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probc'):-1", float, doc="HLT PNet tagger c raw score", precision=12),
        hltPNet_probuds = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probuds'):-1", float, doc="HLT PNet tagger uds raw score", precision=12),
        hltPNet_probg = Var("?(pt>=30)&&(abs(eta)<=2.6)?bDiscriminator('hltParticleNetONNXJetTags:probg'):-1", float, doc="HLT PNet tagger g raw score", precision=12),
    ),
)

############
# AK8PFJet #
############

# Jet Energy Correction
patHLTAK8PFJetCorrFactor = patHLTAK4PFJetCorrFactor.clone(
    src = "hltAK8PFJets",
    payload = "AK8PFHLT",
)

patHLTAK8PFJet = patHLTAK4CaloJet.clone(
    jetSource = "hltAK8PFJets",
    addJetCorrFactors = True,
    jetCorrFactorsSource = ["patHLTAK8PFJetCorrFactor"],
    addBTagInfo = True,
    addDiscriminators = True,
    matchJetTagByDeltaR = True,
    maxJetTagDeltaR = cms.double(0.0001),
    discriminatorSources = [
        "hltParticleNetONNXJetTagsAK8:probHtt",
        "hltParticleNetONNXJetTagsAK8:probHtm",
        "hltParticleNetONNXJetTagsAK8:probHte",
        "hltParticleNetONNXJetTagsAK8:probHbb",
        "hltParticleNetONNXJetTagsAK8:probHcc",
        "hltParticleNetONNXJetTagsAK8:probHqq",
        "hltParticleNetONNXJetTagsAK8:probHgg",
        "hltParticleNetONNXJetTagsAK8:probQCD2hf",
        "hltParticleNetONNXJetTagsAK8:probQCD1hf",
        "hltParticleNetONNXJetTagsAK8:probQCD0hf",
        ],
)

# output to table
patHLTAK8PFJetTable = simplePATJetFlatTableProducer.clone(
    src = cms.InputTag("patHLTAK8PFJet"),
    name = cms.string("HLTAK8PFJet"),
    doc = cms.string("from hltAK8PFJets"),
    cut = cms.string(""),
    variables = cms.PSet(
        PFJetVariables,
        rawFactor = Var("1.-jecFactor('Uncorrected')", float, doc="1 - Factor to get back to raw pT", precision=10),
        hltPNetAK8_probHtt = Var("bDiscriminator('hltParticleNetONNXJetTags:probHtt')", float, doc="HLT PNetAK8 tagger Htt raw score", precision=12),
        hltPNetAK8_probHtm = Var("bDiscriminator('hltParticleNetONNXJetTags:probHtm')", float, doc="HLT PNetAK8 tagger Htm raw score", precision=12),
        hltPNetAK8_probHte = Var("bDiscriminator('hltParticleNetONNXJetTags:probHte')", float, doc="HLT PNetAK8 tagger Hte raw score", precision=12),
        hltPNetAK8_probHbb = Var("bDiscriminator('hltParticleNetONNXJetTags:probHbb')", float, doc="HLT PNetAK8 tagger Hbb raw score", precision=12),
        hltPNetAK8_probHcc = Var("bDiscriminator('hltParticleNetONNXJetTags:probHcc')", float, doc="HLT PNetAK8 tagger Hcc raw score", precision=12),
        hltPNetAK8_probHqq = Var("bDiscriminator('hltParticleNetONNXJetTags:probHqq')", float, doc="HLT PNetAK8 tagger Hqq raw score", precision=12),
        hltPNetAK8_probHgg = Var("bDiscriminator('hltParticleNetONNXJetTags:probHgg')", float, doc="HLT PNetAK8 tagger Hgg raw score", precision=12),
        hltPNetAK8_probQCD2hf = Var("bDiscriminator('hltParticleNetONNXJetTags:probQCD2hf')", float, doc="HLT PNetAK8 tagger QCD2hf raw score", precision=12),
        hltPNetAK8_probQCD1hf = Var("bDiscriminator('hltParticleNetONNXJetTags:probQCD1hf')", float, doc="HLT PNetAK8 tagger QCD1hf raw score", precision=12),
        hltPNetAK8_probQCD0hf = Var("bDiscriminator('hltParticleNetONNXJetTags:probQCD0hf')", float, doc="HLT PNetAK8 tagger QCD0hf raw score", precision=12),
    ),
)

#######
# Rho #
#######

hltRhoTable = cms.EDProducer("GlobalVariablesTableProducer",
    name = cms.string("HLTRho"),
    variables = cms.PSet(
        fixedGridRhoFastjetAll = ExtVar(cms.InputTag("hltFixedGridRhoFastjetAll"), "double", doc = "rho from all PF Candidates, used e.g. for JECs, from hltFixedGridRhoFastjetAll"),
        fixedGridRhoFastjetAllCalo = ExtVar(cms.InputTag("hltFixedGridRhoFastjetAllCalo"), "double", doc = "rho from calo towers, used e.g. for JECs, from hltFixedGridRhoFastjetAllCalo"),
    )
)

###########
# CaloMET #
###########

hltCaloMETTable = cms.EDProducer("SimpleCaloMETFlatTableProducer",
    src = cms.InputTag("hltMet"),
    name = cms.string("HLTCaloMET"),
    doc = cms.string("from hltMet"),
    singleton = cms.bool(True),
    variables = cms.PSet(
        PTVars,
        sumEt = Var("sumEt()", float, doc="scalar sum of Et", precision=10),
    )
)

#########
# PFMET #
#########

hltPFMETTable = cms.EDProducer("SimplePFMETFlatTableProducer",
    src = cms.InputTag("hltPFMETProducer"),
    name = cms.string("HLTPFMET"),
    doc = cms.string("from hltPFMETProducer"),
    singleton = cms.bool(True),
    variables = cms.PSet(
        PTVars,
        sumEt = Var("sumEt()", float, doc="scalar sum of Et", precision=10),
    )
)

#######
# MHT #
#######

hltMHTTable = cms.EDProducer("SimpleMETFlatTableProducer",
    src = cms.InputTag("hltPFHTForMC"),
    name = cms.string("HLTMHT"),
    doc = cms.string("from hltPFHTForMC"),
    singleton = cms.bool(True),
    variables = cms.PSet(
        PTVars,
        sumEt = Var("sumEt()", float, doc="scalar sum of Et", precision=10),
    )
)

##########
# GenJet #
##########

#genJetTable = cms.EDProducer("SimpleGenJetFlatTableProducer",
#    src = cms.InputTag("ak4GenJets"),
#    name = cms.string("GenJet"),
#    doc = cms.string("from ak4GenJets"),
#    cut = cms.string(""),
#    variables = cms.PSet(
#        P4Vars,
#    ),
#)
#
#genJetNoNuTable = genJetTable.clone(
#    src = cms.InputTag("ak4GenJetsNoNu"),
#    name = cms.string("GenJetNoNu"),
#    doc = cms.string("from ak4GenJetsNoNu"),
#)
#
#genJetAK8Table = genJetTable.clone(
#    src = cms.InputTag("ak8GenJets"),
#    name = cms.string("GenJetAK8"),
#    doc = cms.string("from ak8GenJets"),
#)
#
#genJetAK8NoNuTable = genJetTable.clone(
#    src = cms.InputTag("ak8GenJetsNoNu"),
#    name = cms.string("GenJetAK8NoNu"),
#    doc = cms.string("from ak8GenJetsNoNu"),
#)

from PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi import patJetGenJetMatch 
patHLTAK4PFJetGenJetMatch = patJetGenJetMatch.clone(
    src = cms.InputTag("hltAK4PFJets"),
    matched = cms.InputTag("slimmedGenJets"),
    resolveByMatchQuality = cms.bool(True)
)

from PhysicsTools.PatAlgos.mcMatchLayer0.jetFlavourId_cff import patJetFlavourAssociation
patHLTAK4PFJetFlavourAssociation = patJetFlavourAssociation.clone(
    jets = cms.InputTag("hltAK4PFJets")
)

from PhysicsTools.NanoAOD.jetMC_cff import jetMCTable
patHLTAK4PFJetMCTable = jetMCTable.clone(
    src = patHLTAK4PFJetTable.src,
    name = patHLTAK4PFJetTable.name,
)

######
# PV #
######

hltPixelVertexTable = cms.EDProducer("SimpleVertexFlatTableProducer",
    src = cms.InputTag("hltPixelVertices"),
    name = cms.string("HLTPixelVertex"),
    doc = cms.string("from hltPixelVertices"),
    variables = cms.PSet(
        chi2 = Var("chi2()", float, doc="chi2", precision=8),
        ndof = Var("ndof()", float, doc="number of degree of freedom", precision=8),
        x = Var("x()", float, doc="position x coordinate", precision=10),
        y = Var("y()", float, doc="position y coordinate", precision=10),
        z = Var("z()", float, doc="position z coordinate", precision=16),
        cov00 = Var("covariance(0,0)", float, doc="vertex covariance (0,0)", precision=10),
        cov10 = Var("covariance(1,0)", float, doc="vertex covariance (1,0)", precision=10),
        cov11 = Var("covariance(1,1)", float, doc="vertex covariance (1,1)", precision=10),
        cov20 = Var("covariance(2,0)", float, doc="vertex covariance (2,0)", precision=10),
        cov21 = Var("covariance(2,1)", float, doc="vertex covariance (2,1)", precision=10),
        cov22 = Var("covariance(2,2)", float, doc="vertex covariance (2,2)", precision=10),
        isFake = Var("isFake()", bool, doc="is Fake"),
    ),
)

hltVertexPFTable = hltPixelVertexTable.clone(
    src = cms.InputTag("hltVerticesPF"),
    name = cms.string("HLTVertexPF"),
    doc = cms.string("from hltVerticesPF"),
)

######
# SV #
######

hltDeepInclusiveVertexFinderPFTable = cms.EDProducer("SimpleSecondaryVertexFlatTableProducer",
    src = cms.InputTag("hltDeepInclusiveVertexFinderPF"),
    name = cms.string("HLTDeepInclusiveVertexFinderPF"),
    doc = cms.string("from hltDeepInclusiveVertexFinderPF"),
    variables = cms.PSet(
        P4Vars,
        chi2 = Var("vertexChi2()", float, doc="chi2", precision=8),
        ndof = Var("vertexNdof()", float, doc="number of degree of freedom", precision=8),
        x = Var("position().x()", float, doc="position x coordinate", precision=10),
        y = Var("position().y()", float, doc="position y coordinate", precision=10),
        z = Var("position().z()", float, doc="position z coordinate", precision=16),
        cov00 = Var("vertexCovariance(0,0)", float, doc="vertex covariance (0,0)", precision=10),
        cov10 = Var("vertexCovariance(1,0)", float, doc="vertex covariance (1,0)", precision=10),
        cov11 = Var("vertexCovariance(1,1)", float, doc="vertex covariance (1,1)", precision=10),
        cov20 = Var("vertexCovariance(2,0)", float, doc="vertex covariance (2,0)", precision=10),
        cov21 = Var("vertexCovariance(2,1)", float, doc="vertex covariance (2,1)", precision=10),
        cov22 = Var("vertexCovariance(2,2)", float, doc="vertex covariance (2,2)", precision=10),
        ntracks = Var("numberOfDaughters()", "uint8", doc = "number of tracks"),
    ),
)

hltDeepInclusiveMergedVertexPFTable = hltDeepInclusiveVertexFinderPFTable.clone(
    src = cms.InputTag("hltDeepInclusiveMergedVerticesPF"),
    name = cms.string("HLTDeepInclusiveMergedVertexPF"),
    doc = cms.string("from hltDeepInclusiveMergedVerticesPF"),
)
