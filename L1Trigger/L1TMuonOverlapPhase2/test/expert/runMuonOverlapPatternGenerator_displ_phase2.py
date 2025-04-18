# -*- coding: utf-8 -*-
import FWCore.ParameterSet.Config as cms
process = cms.Process("L1TMuonEmulation")
import os
import sys
#import commands
import re
from os import listdir
from os.path import isfile, join
import fnmatch

process.load("FWCore.MessageLogger.MessageLogger_cfi")

verbose = True

filesNameLike = sys.argv[2]

#version = "ExtraplMB1nadMB2DTQualAndEtaFloatP_atan_ValueP1Scale_t20_deltaPhiVsPhiB_SingleMu_iPt_10files"
#version = "noExtrapl_ValueP1Scale_t18_qualConverted_min4_ipT1_deltaPhiVsPhiRef_fixedDTScale"
#version = "ExtraplMB1nadMB2DTQualAndEtaFloatP_atan_ValueP1Scale_t20_deltaPhiVsPhiB_300files"
version = "ExtraplMB1nadMB2DTQualAndEtaFixedP_ValueP1Scale_t20_deltaPhiVsPhiB_500files"

if verbose: 
    process.MessageLogger = cms.Service("MessageLogger",
       #suppressInfo       = cms.untracked.vstring('AfterSource', 'PostModule'),
       destinations   = cms.untracked.vstring(
                                               #'detailedInfo',
                                               #'critical',
                                               #'cout',
                                               #'cerr',
                                               'omtfEventPrint'
                    ),
       categories        = cms.untracked.vstring('l1tOmtfEventPrint', 'OMTFReconstruction'),
       omtfEventPrint = cms.untracked.PSet(    
                         filename  = cms.untracked.string("Patterns_dispalced_test_" + version + "_" + filesNameLike),
                         extension = cms.untracked.string('.txt'),                
                         threshold = cms.untracked.string('INFO'),
                         default = cms.untracked.PSet( limit = cms.untracked.int32(0) ), 
                         #INFO   =  cms.untracked.int32(0),
                         #DEBUG   = cms.untracked.int32(0),
                         l1tOmtfEventPrint = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) ),
                         OMTFReconstruction = cms.untracked.PSet( limit = cms.untracked.int32(1000000000) )
                       ),
       debugModules = cms.untracked.vstring('simOmtfPhase2Digis') 
       #debugModules = cms.untracked.vstring('*')
    )

    #process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(100)
if not verbose:
    process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(-1)
    process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(False), 
                                         #SkipEvent = cms.untracked.vstring('ProductNotFound') 
                                     )
    
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
#process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
#process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.EventContent.EventContent_cff')
#process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtendedRun4D95Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
#process.load('Configuration.StandardSequences.RawToDigi_cff')
#process.load('Configuration.StandardSequences.SimL1Emulator_cff')
#process.load('Configuration.StandardSequences.SimPhase2L1GlobalTriggerEmulator_cff')
#process.load('L1Trigger.Configuration.Phase2GTMenus.SeedDefinitions.prototypeSeeds')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '131X_mcRun4_realistic_v9', '')

chosenFiles = []

cscBx = 8

if filesNameLike == 'displHighPt' : # displaced muon sample
    cscBx = 8
    #path = '/eos/user/c/cericeci/forOMTF/OMTF_PhaseII_FixedTiming/'
    path =  '/eos/user/a/asotorod/Samples/OMTF-L1/OMTF_fixedTiming/'
    
    fileCnt = 20 
    firstFile = 1 #1001            
    for i in range(firstFile, firstFile + fileCnt, 1):
        filePathName = path + "custom_Displaced_" + str(i) + "_numEvent5000.root"
        if isfile(filePathName) :
            #chosenFiles.append('file://' + path + "custom_Displaced_Run3_" + str(i) + "_numEvent1000.root") 
            #chosenFiles.append('file://' + path + "custom_Displaced_Run3_" + str(i) + "_numEvent2000.root") 
            chosenFiles.append('file://' + filePathName)

### N.B. for phase1 samples there are no DT hits
elif filesNameLike == 'mcWaw2023' :
    cscBx = 8
      
    for iPt in [0, 1, 2] :
        for charge in [0, 2] : 
            path = "/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch" + str(charge) + "_iPt" + str(iPt) + "_12_5_2_p1_04_04_2023/12_5_2_p1_04_04_2023/"
            #path = '/eos/user/a/akalinow/Data/SingleMu/12_5_2_p1_04_04_2023/SingleMu_ch0_iPt1_12_5_2_p1_04_04_2023/12_5_2_p1_04_04_2023/230404_084317/0000/'
              
            root_files = []
            for root, dirs, files in os.walk(path):
                for file in fnmatch.filter(files, '*.root'):
                    root_files.append(os.path.join(root, file))

            file_cnt = 10
                   
            file_num = 0    
            for root_file in root_files :
                if isfile(root_file) :
                    chosenFiles.append('file://' + root_file)
                    file_num += 1
                else :
                    print("file not found!!!!!!!: " + root_file)   
                    
                if file_num >= file_cnt :
                    break  
                
elif filesNameLike == 'XTo2LLPTo1Mu' :            
    path = "/eos/user/a/almuhamm/ZMu_Test/ExoticLLP/XTo2LLPTo1Mu_CTau8000_Phase2Exotic/231119_223651/0000/"
    root_files = []
    for root, dirs, files in os.walk(path):
        for file in fnmatch.filter(files, '*.root'):
            root_files.append(os.path.join(root, file))

    file_cnt = 500
           
    file_num = 0    
    for root_file in root_files :
        if isfile(root_file) :
            chosenFiles.append('file://' + root_file)
            file_num += 1
        else :
            print("file not found!!!!!!!: " + root_file)   
            
        if file_num >= file_cnt :
            break  
                
print("chosenFiles")
for chFile in chosenFiles:
    print(chFile)

if len(chosenFiles) == 0 :
    print("no files selected!!!!!!!!!!!!!!! (argumetn should be e.g. 20_p")
    exit 
                            
# input files (up to 255 files accepted)
process.source = cms.Source('PoolSource',
fileNames = cms.untracked.vstring( 
    *(list(chosenFiles)) ),
    skipEvents =  cms.untracked.uint32(0),
    inputCommands=cms.untracked.vstring(
        'keep *',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_CSC_HLT',
        'drop l1tEMTFHit2016Extras_simEmtfDigis_RPC_HLT',
        'drop l1tEMTFHit2016s_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016Extras_simEmtfDigis__HLT',
        'drop l1tEMTFTrack2016s_simEmtfDigis__HLT')
)
                    
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1))


#Calibrate Digis
process.load("L1Trigger.DTTriggerPhase2.CalibratedDigis_cfi")
process.CalibratedDigis.dtDigiTag = "simMuonDTDigis" 
process.CalibratedDigis.scenario = 0

#DTTriggerPhase2
process.load("L1Trigger.DTTriggerPhase2.dtTriggerPhase2PrimitiveDigis_cfi")
process.dtTriggerPhase2PrimitiveDigis.debug = False
process.dtTriggerPhase2PrimitiveDigis.dump = False
process.dtTriggerPhase2PrimitiveDigis.scenario = 0

#process.TFileService = cms.Service("TFileService", fileName = cms.string('omtfAnalysis1_1.root'), closeFileFast = cms.untracked.bool(True) )
                                   
####OMTF Emulator
process.load('L1Trigger.L1TMuonOverlapPhase2.simOmtfPhase2Digis_cfi')

#needed by candidateSimMuonMatcher
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
#process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")

process.simOmtfPhase2Digis.candidateSimMuonMatcher = cms.bool(True)
process.simOmtfPhase2Digis.simTracksTag = cms.InputTag('g4SimHits')
process.simOmtfPhase2Digis.simVertexesTag = cms.InputTag('g4SimHits')
process.simOmtfPhase2Digis.muonMatcherFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/muonMatcherHists_100files_smoothStdDev_withOvf.root")


process.simOmtfPhase2Digis.bxMin = cms.int32(0)
process.simOmtfPhase2Digis.bxMax = cms.int32(0)

process.simOmtfPhase2Digis.dumpResultToXML = cms.bool(False)
process.simOmtfPhase2Digis.eventCaptureDebug = cms.bool(False)

process.simOmtfPhase2Digis.patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuonOverlapPhase1/test/expert/omtf/Patterns_template.xml")
#process.simOmtfPhase2Digis.patternsXMLFiles = cms.VPSet(cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_plus_v1.xml")),
#                                                       cms.PSet(patternsXMLFile = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/GPs_parametrised_minus_v1.xml"))
#)

#process.simOmtfPhase2Digis.patternGenerator = cms.string("patternGen")
process.simOmtfPhase2Digis.patternGenerator = cms.string("2DHists")
#process.simOmtfPhase2Digis.patternGenerator = cms.string("deltaPhiVsPhiRef")

process.simOmtfPhase2Digis.patternType = cms.string("GoldenPatternWithStat")
process.simOmtfPhase2Digis.generatePatterns = cms.bool(True)
process.simOmtfPhase2Digis.optimisedPatsXmlFile = cms.string("Patterns_dispalced_test_" + version + "_" + filesNameLike + ".xml")

process.simOmtfPhase2Digis.rpcMaxClusterSize = cms.int32(3)
process.simOmtfPhase2Digis.rpcMaxClusterCnt = cms.int32(2)
process.simOmtfPhase2Digis.rpcDropAllClustersIfMoreThanMax = cms.bool(True)

process.simOmtfPhase2Digis.minCSCStubRME12 = cms.int32(410) #[cm]
process.simOmtfPhase2Digis.minCSCStubR = cms.int32(500) #[cm]

process.simOmtfPhase2Digis.minDtPhiQuality = cms.int32(2)
process.simOmtfPhase2Digis.minDtPhiBQuality = cms.int32(4)

#process.simOmtfPhase2Digis.dtPhiBUnitsRad = cms.int32(2048) #2048 is the orginal phase2 scale, 512 is the phase1 scale

process.simOmtfPhase2Digis.dtRefHitMinQuality =  cms.int32(4)

process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB1 = cms.bool(True)
process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB2 = cms.bool(True)
#process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB1 = cms.bool(False)
#process.simOmtfPhase2Digis.usePhiBExtrapolationFromMB2 = cms.bool(False)

process.simOmtfPhase2Digis.useStubQualInExtr  = cms.bool(True)
process.simOmtfPhase2Digis.useEndcapStubsRInExtr  = cms.bool(True)
process.simOmtfPhase2Digis.useFloatingPointExtrapolation  = cms.bool(False)
#process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("ExtrapolationFactors_DTQualAndEtaValueP1Scale.xml")
#process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("ExtrapolationFactors_simple.xml")
process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("L1Trigger/L1TMuon/data/omtf_config/ExtrapolationFactors_ExtraplMB1nadMB2DTQual_ValueP1Scale_t20.xml")
#process.simOmtfPhase2Digis.extrapolFactorsFilename = cms.FileInPath("")

process.simOmtfPhase2Digis.stubEtaEncoding = cms.string("valueP1Scale")  
#process.simOmtfPhase2Digis.stubEtaEncoding = cms.string("bits")   

process.simOmtfPhase2Digis.goldenPatternResultFinalizeFunction = cms.int32(3) ## is needed here , becasue it just counts the number of layers with a stub
process.simOmtfPhase2Digis.lctCentralBx = cms.int32(cscBx);#<<<<<<<<<<<<<<<<!!!!!!!!!!!!!!!!!!!!TODO this was changed in CMSSW 10(?) to 8. if the data were generated with the previous CMSSW then you have to use 6



#process.dumpED = cms.EDAnalyzer("EventContentAnalyzer")
#process.dumpES = cms.EDAnalyzer("PrintEventSetupContent")

process.L1TMuonSeq = cms.Sequence( process.simOmtfPhase2Digis 
                                   #+ process.dumpED
                                   #+ process.dumpES
)

process.L1TMuonPath = cms.Path(process.L1TMuonSeq)

#process.out = cms.OutputModule("PoolOutputModule", 
#   fileName = cms.untracked.string("l1tomtf_superprimitives1.root")
#)

#process.output_step = cms.EndPath(process.out)
#process.schedule = cms.Schedule(process.L1TMuonPath)
#process.schedule.extend([process.output_step])
