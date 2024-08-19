import FWCore.ParameterSet.Config as cms
from  PhysicsTools.NanoAOD.common_cff import *
from PhysicsTools.NanoAOD.simpleCandidateFlatTableProducer_cfi import simpleCandidateFlatTableProducer

################
HLTPFJETVARS = cms.PSet(P4Vars,
        area = Var("jetArea()", float, doc="area", precision=10),
        # nMuons = Var("?hasOverlaps('muons')?overlaps('muons').size():0", "uint8", doc="number of muons in the jet"),
        # muonIdx1 = Var("?overlaps('muons').size()>0?overlaps('muons')[0].key():-1", "int16", doc="index of first matching muon"),
        # muonIdx2 = Var("?overlaps('muons').size()>1?overlaps('muons')[1].key():-1", "int16", doc="index of second matching muon"),
        # electronIdx1 = Var("?overlaps('electrons').size()>0?overlaps('electrons')[0].key():-1", "int16", doc="index of first matching electron"),
        # electronIdx2 = Var("?overlaps('electrons').size()>1?overlaps('electrons')[1].key():-1", "int16", doc="index of second matching electron"),
        # nElectrons = Var("?hasOverlaps('electrons')?overlaps('electrons').size():0", "uint8", doc="number of electrons in the jet"),
        # svIdx1 = Var("?overlaps('vertices').size()>0?overlaps('vertices')[0].key():-1", "int16", doc="index of first matching secondary vertex"),
        # svIdx2 = Var("?overlaps('vertices').size()>1?overlaps('vertices')[1].key():-1", "int16", doc="index of second matching secondary vertex"),
        # nSVs = Var("?hasOverlaps('vertices')?overlaps('vertices').size():0", "uint8", doc="number of secondary vertices in the jet"),
        # jetId = Var("userInt('tightId')*2+4*userInt('tightIdLepVeto')", "uint8",doc="Jet ID flag: bit2 is tight, bit3 is tightLepVeto"),
        # hfsigmaEtaEta = Var("userFloat('hfsigmaEtaEta')",float,doc="sigmaEtaEta for HF jets (noise discriminating variable)",precision=10),
        # hfsigmaPhiPhi = Var("userFloat('hfsigmaPhiPhi')",float,doc="sigmaPhiPhi for HF jets (noise discriminating variable)",precision=10),
        # hfcentralEtaStripSize = Var("userInt('hfcentralEtaStripSize')", int, doc="eta size of the central tower strip in HF (noise discriminating variable)"),
        # hfadjacentEtaStripsSize = Var("userInt('hfadjacentEtaStripsSize')", int, doc="eta size of the strips next to the central tower strip in HF (noise discriminating variable)"),
        nConstituents = Var("numberOfDaughters()","uint8",doc="Number of particles in the jet"),
        chMultiplicity = Var("chargedMultiplicity()","uint8",doc="Number of charged particles in the jet"),
        neMultiplicity = Var("neutralMultiplicity()","uint8",doc="Number of neutral particles in the jet"),
        # rawFactor = Var("1.-jecFactor('Uncorrected')",float,doc="1 - Factor to get back to raw pT",precision=6),
        chHEF = Var("chargedHadronEnergyFraction()", float, doc="charged Hadron Energy Fraction", precision= 6),
        neHEF = Var("neutralHadronEnergyFraction()", float, doc="neutral Hadron Energy Fraction", precision= 6),
        chEmEF = Var("chargedEmEnergyFraction()", float, doc="charged Electromagnetic Energy Fraction", precision= 6),
        neEmEF = Var("neutralEmEnergyFraction()", float, doc="neutral Electromagnetic Energy Fraction", precision= 6),
        hfHEF = Var("HFHadronEnergyFraction()",float,doc="hadronic Energy Fraction in HF",precision= 6),
        hfEmEF = Var("HFEMEnergyFraction()",float,doc="electromagnetic Energy Fraction in HF",precision= 6),
        muEF = Var("muonEnergyFraction()", float, doc="muon Energy Fraction", precision= 6),
     
        # variables from JMENano
        chHadMultiplicity = Var("chargedHadronMultiplicity()","uint8",doc="number of charged hadrons in the jet"),
        neHadMultiplicity = Var("neutralHadronMultiplicity()","uint8",doc="number of neutral hadrons in the jet"),
        hfHadMultiplicity = Var("HFHadronMultiplicity()", "uint8",doc="number of HF hadrons in the jet"),
        hfEMMultiplicity  = Var("HFEMMultiplicity()","uint8",doc="number of HF EMs in the jet"),
        muMultiplicity    = Var("muonMultiplicity()","uint8",doc="number of muons in the jet"),
        elMultiplicity    = Var("electronMultiplicity()","uint8",doc="number of electrons in the jet"),
        phoMultiplicity   = Var("photonMultiplicity()","uint8",doc="number of photons in the jet"),
)

# create empty product if not exists
HLTAK4PFJets = cms.EDProducer("PFJetCollectionProductNotFoundEmptyProducer", src=cms.InputTag("hltAK4PFJets"))
HLTAK8PFJets = cms.EDProducer("PFJetCollectionProductNotFoundEmptyProducer", src=cms.InputTag("hltAK8PFJets"))
HLTAK4PFJetsCorrectedMatchedToCaloJets10 = cms.EDProducer("PFJetCollectionProductNotFoundEmptyProducer", src=cms.InputTag("hltPFJetsCorrectedMatchedToCaloJets10"))
HLTAK8PFJetsCorrectedMatchedToCaloJets10  = cms.EDProducer("PFJetCollectionProductNotFoundEmptyProducer", src=cms.InputTag("hltPFJetsCorrectedMatchedToCaloJets10AK8"))
HLTFixedGridRhoFastjetAll = cms.EDProducer("DoubleProductNotFoundEmptyProducer", src=cms.InputTag("hltfixedGridRhoFastjetAll"))

# HLTAK4PFJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("hltAK4PFJets"),
#     cut = cms.string("pt>=-1."),
#     name = cms.string("HLTAK4PFJet"),
#     doc = cms.string("HLT AK4 jets"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS,
# )

# HLTAK8PFJetTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("hltAK8PFJets"),
#     cut = cms.string(""),
#     name = cms.string("HLTAK8PFJet"),
#     doc = cms.string("HLT AK8 jets"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS,
# )

# HLTAK4PFJetsCorrectedMatchedToCaloJets10Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("hltPFJetsCorrectedMatchedToCaloJets10"),
#     cut = cms.string("pt>=-1."),
#     name = cms.string("HLTAK4PFJetCorrectedMatchedToCaloJets10"),
#     doc = cms.string("HLT AK4 PF jets corrected and matched to Calo jets"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS,
# )

# HLTAK8PFJetsCorrectedMatchedToCaloJets10Table = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("hltPFJetsCorrectedMatchedToCaloJets10AK8"),
#     cut = cms.string(""),
#     name = cms.string("HLTAK8PFJetCorrectedMatchedToCaloJets10"),
#     doc = cms.string("HLT AK8 PF jets corrected and matched to Calo jets"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS,
# )

HLTRhoTable = cms.EDProducer("GlobalVariablesTableProducer",
    name = cms.string("HLTRho"),
    variables = cms.PSet(
        fixedGridRhoFastjetAll = ExtVar(cms.InputTag("HLTFixedGridRhoFastjetAll"), "double", doc="rho from all PF Candidates, used e.g. for JECs")
    )
)

HLTRawAK4PFMET = cms.EDProducer("PFMETProducer",
    src = cms.InputTag("HLTAK4PFJets"),
)

HLTRawAK4PFMETTable = cms.EDProducer("SimplePFMETFlatTableProducer",
    src  = cms.InputTag("HLTRawAK4PFMET"),
    name = cms.string("HLTRawAK4PFMET"),
    doc = cms.string("HLT AK4 PF MET from summing HLT AK4 PF jets"),
    singleton = cms.bool(True),
    variables = cms.PSet(PTVars,
        sumEt = Var("sumEt()", float, doc="scalar sum of Et", precision=10),
    )
)

HLTRawAK8PFMET = cms.EDProducer("PFMETProducer",
    src = cms.InputTag("HLTAK8PFJets"),
)

HLTRawAK8PFMETTable = cms.EDProducer("SimplePFMETFlatTableProducer",
    src  = cms.InputTag("HLTRawAK8PFMET"),
    name = cms.string("HLTRawAK8PFMET"),
    doc = cms.string("HLT AK8 PF MET from summing HLT AK8 PF jets"),
    singleton = cms.bool(True),
    variables = cms.PSet(PTVars,
        sumEt = Var("sumEt()", float, doc="scalar sum of Et", precision=10),
    )
)

HLTAK4JetForJEC = cms.EDProducer("HLTJetForJECProducer",
    corrected = cms.InputTag("HLTAK4PFJetsCorrectedMatchedToCaloJets10"),
    raw = cms.InputTag("HLTAK4PFJets"),
    ptMin = cms.double(8.),
    maxDeltaR = cms.double(0.2),
)

# HardestHLTJetForJEC  = cms.EDProducer("LargestPtCandSelector",
#     src = cms.InputTag("HLTJetForJEC", "corrected"),
#     maxNumber = cms.uint32(10)
# )

HLTAK4PFJetCorrectedMatchedToCaloJets10ForJECTable = cms.EDProducer("SimplePFJetFlatTableProducer",
    src = cms.InputTag("HLTAK4JetForJEC", "corrected"),
    cut = cms.string(""),
    name = cms.string("HLTAK4PFJetCorrectedMatchedToCaloJets10ForJEC"),
    doc = cms.string("HLT AK4 PF jets corrected and matched to Calo jets, for JEC"),
    singleton = cms.bool(False),
    extension = cms.bool(False), # this is not an extension
    variables = HLTPFJETVARS,
    externalVariables = cms.PSet(
      rawFactor = ExtVar(cms.InputTag("HLTAK4JetForJEC", "rawFactor"), float, doc="1 - Factor to get back to raw pT", precision=6),
    )
)


# HLTAK4PFJetForJECTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("HLTJetForJEC", "raw"),
#     cut = cms.string("pt"),
#     name = cms.string("HLTAK4PFJetForJEC"),
#     doc = cms.string("HLT AK4 jets, for JEC"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS
# )

HLTAK8JetForJEC = cms.EDProducer("HLTJetForJECProducer",
    corrected = cms.InputTag("HLTAK8PFJetsCorrectedMatchedToCaloJets10"),
    raw = cms.InputTag("HLTAK8PFJets"),
    ptMin = cms.double(15.),
    maxDeltaR = cms.double(0.4),
)

# HardestHLTJetForJEC  = cms.EDProducer("LargestPtCandSelector",
#     src = cms.InputTag("HLTJetForJEC", "corrected"),
#     maxNumber = cms.uint32(10)
# )

HLTAK8PFJetCorrectedMatchedToCaloJets10ForJECTable = cms.EDProducer("SimplePFJetFlatTableProducer",
    src = cms.InputTag("HLTAK8JetForJEC", "corrected"),
    cut = cms.string(""),
    name = cms.string("HLTAK8PFJetCorrectedMatchedToCaloJets10ForJEC"),
    doc = cms.string("HLT AK8 PF jets corrected and matched to Calo jets, for JEC"),
    singleton = cms.bool(False),
    extension = cms.bool(False), # this is not an extension
    variables = HLTPFJETVARS,
    externalVariables = cms.PSet(
      rawFactor = ExtVar(cms.InputTag("HLTAK8JetForJEC", "rawFactor"), float, doc="1 - Factor to get back to raw pT", precision=6),
    )
)


# HLTAK8PFJetForJECTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
#     src = cms.InputTag("HLTAK8JetForJEC", "raw"),
#     cut = cms.string("pt"),
#     name = cms.string("HLTAK8PFJetForJEC"),
#     doc = cms.string("HLT AK8 jets, for JEC"),
#     singleton = cms.bool(False),
#     extension = cms.bool(False), # this is not an extension
#     variables = HLTPFJETVARS
# )
