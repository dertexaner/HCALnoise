# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO --data -s RAW2DIGI,RECO --scenario cosmics --filein file:5C1B1DE5-9B38-E211-A048-001D09F24FBA.root --fileout DummyOutput.root --conditions GR_R_73_V1A::All --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContentCosmics_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
     wantSummary = cms.untracked.bool(True)
)
process.MessageLogger = cms.Service("MessageLogger",
 cout = cms.untracked.PSet( 
   default = cms.untracked.PSet( ## kill all messages in the log 
   limit = cms.untracked.int32(0) 
  ), 
  FwkJob = cms.untracked.PSet( ## but FwkJob category - those unlimitted 
   limit = cms.untracked.int32(-1) 
  ),
  FwkReport = cms.untracked.PSet(
   reportEvery = cms.untracked.int32(1), ## print event record number
   limit = cms.untracked.int32(-1) 
  ),
  FwkSummary = cms.untracked.PSet(
    optionalPSet = cms.untracked.bool(True),
  #  reportEvery = cms.untracked.int32(100),
  #  limit = cms.untracked.int32(10000000)
  )
 ),
 categories = cms.untracked.vstring('FwkJob','FwkReport','FwkSummary'), 
 destinations = cms.untracked.vstring('cout')
)


# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('/store/data/Commissioning2014/MinimumBias/RAW/v3/000/229/684/00000/BE2D17C9-D469-E411-8A9E-02163E010DF4.root')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('RECO nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition
#process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.RECOSIMEventContent.outputCommands,
#    fileName = cms.untracked.string('DummyOutput.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('')
#    )
#)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'GR_R_73_V1A::All', '')
#Recommended GT:  GR_R_73_V1A::All (both for 7_3_1 and 7_4_0_pre5)


# Hcal noise analyzers
process.HBHENoiseFilterResultProducer = cms.EDProducer(
   'HBHENoiseFilterResultProducer',
   noiselabel = cms.InputTag('hcalnoise'),
   minRatio = cms.double(-999),
   maxRatio = cms.double(999),
   minHPDHits = cms.int32(17),
   minRBXHits = cms.int32(999),
   minHPDNoOtherHits = cms.int32(10),
   minZeros = cms.int32(10),
   minHighEHitTime = cms.double(-9999.0),
   maxHighEHitTime = cms.double(9999.0),
   maxRBXEMF = cms.double(-999.0),
   minNumIsolatedNoiseChannels = cms.int32(10),
   minIsolatedNoiseSumE = cms.double(50.0),
   minIsolatedNoiseSumEt = cms.double(25.0),
   useTS4TS5 = cms.bool(False),
   useRBXRechitR45Loose = cms.bool(False),
   useRBXRechitR45Tight = cms.bool(False),
   IgnoreTS4TS5ifJetInLowBVRegion = cms.bool(True),
   jetlabel = cms.InputTag('ak5PFJets'),
   maxjetindex = cms.int32(0),
   maxNHF = cms.double(0.9)
   )
process.TFileService = cms.Service("TFileService",
   fileName = cms.string("NoiseTree_Commissionig2014_MinBias_v3_229684_R45_Jan23.root")
   )
process.ExportTree = cms.EDAnalyzer("HcalNoiseAnalyzer",
  HBHERecHitCollection = cms.untracked.string('hbhereco'),
  IsCosmic             = cms.untracked.bool(True)
)

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstructionCosmics * process.HBHENoiseFilterResultProducer * process.ExportTree)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,process.endjob_step)


# To prevent errors of this type (see https://hypernews.cern.ch/HyperNews/CMS/get/physTools/3260/1/1/1/1/1/1/1/1/1/1/1.html):
#
#Begin processing the 1703rd record. Run 229684, Event 185885, LumiSection 20 at 27-Jan-2015 13:30:42.210 CET
#Begin processing the 1704th record. Run 229684, Event 186396, LumiSection 20 at 27-Jan-2015 13:30:42.786 CET
#----- Begin Fatal Exception 27-Jan-2015 13:30:42 CET-----------------------
#An exception of category 'LogicError' occurred while
#   [0] Processing run: 229684 lumi: 20 event: 186396
#   [1] Running path 'reconstruction_step'
#   [2] Calling event method for module CSCRecHitDProducer/'csc2DRecHits'
#Exception Message:
#trying to insert duplicate entry
#----- End Fatal Exception -------------------------------------------------
#27-Jan-2015 13:30:42 CET  Closed file root://eoscms//eos/cms/store/data/Commissioning2014/MinimumBias/RAW/v3/000/229/684/00000/BE2D17C9-D469-E411-8A9E-02163E010DF4.root?svcClass=default
#
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)
