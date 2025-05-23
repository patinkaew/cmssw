import FWCore.ParameterSet.Config as cms

process = cms.Process("ProcessOne")
process.load("CondCore.CondDB.CondDB_cfi")
process.CondDB.DBParameters.authenticationPath = '/afs/cern.ch/cms/DB/conddb/'
#
# Choose the output database
#
#process.CondDBCommon.connect = 'oracle://cms_orcon_prod/CMS_COND_311X_ECAL_LAS'
process.CondDB.connect = 'sqlite_file:EcalLinPed.db'
#process.CondDBCommon.connect = 'oracle://cms_orcoff_prep/CMS_COND_ECAL'

process.MessageLogger = cms.Service("MessageLogger",
    cerr = cms.untracked.PSet(
        enable = cms.untracked.bool(False)
    ),
    cout = cms.untracked.PSet(
        enable = cms.untracked.bool(True)
    ),
    debugModules = cms.untracked.vstring('*')
)

process.source = cms.Source("EmptyIOVSource",
  firstValue = cms.uint64(1),
  lastValue = cms.uint64(1),
  timetype = cms.string('runnumber'),
  interval = cms.uint64(1)
)

process.PoolDBESSource = cms.ESSource("PoolDBESSource",
  process.CondDB,
  toGet = cms.VPSet(
    cms.PSet(
      record = cms.string('EcalTPGPedestalsRcd'),
      tag = cms.string('EcalTPGPedestals_test')
    )
  )
 )

process.PoolDBOutputService = cms.Service("PoolDBOutputService",
  process.CondDB,
  logconnect = cms.untracked.string('sqlite_file:DBLog.db'),
  timetype = cms.untracked.string('runnumber'),
  toPut = cms.VPSet(
    cms.PSet(
      record = cms.string('EcalTPGPedestalsRcd'),
      tag = cms.string('EcalTPGPedestals_test')
    )
  )
)

process.Test1 = cms.EDAnalyzer("ExTestEcalTPGLinPed_Analyzer",
  record = cms.string('EcalTPGPedestalsRcd'),
  Source = cms.PSet(
    debug = cms.bool(True),
  )
)

process.p = cms.Path(process.Test1)
