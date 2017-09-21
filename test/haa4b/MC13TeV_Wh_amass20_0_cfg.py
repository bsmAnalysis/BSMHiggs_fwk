import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


#load run conditions
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'

process.mainNtuplizer.isMC = cms.bool(True)
process.mainNtuplizer.dtag = cms.string("MC13TeV_Wh_amass20")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(0)
process.mainNtuplizer.verbose = cms.bool(False)


process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_1.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_2.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_3.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_4.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_5.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_6.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_7.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_8.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_9.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_10.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_11.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_12.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_13.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_14.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_15.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_16.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_17.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_18.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_19.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_20.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_21.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_22.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_23.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_24.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_25.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_26.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_27.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_28.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_29.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_30.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_31.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_32.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_33.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_34.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_35.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_36.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_37.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_38.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_39.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_40.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_41.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_42.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_43.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_44.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_45.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_46.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_47.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_48.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_49.root",
                                  "root://eoscms.cern.ch//eos/cms/store/user/georgia/h-aa-Madgraph5/final/Wh_production_20/HIG-RunIISummer16MiniAODv2-01613_50.root"
    )
)


process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("analysis.root")
                                  )

process.p = cms.Path(  process.mainNtuplizer )	

