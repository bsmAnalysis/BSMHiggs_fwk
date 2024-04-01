import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
#from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection 

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#load run conditions
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
#------ Declare the correct global tag ------#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '102X_dataRun2_Sep2018ABC_v2'
#process.GlobalTag.globaltag = '102X_dataRun2_Prompt_v13'
#process.GlobalTag.globaltag = '102X_upgrade2018_realistic_v18'

process.options   = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound')
)

updateJetCollection(                                                                                                                                                                              
   process,                                                                                                                                                                                       
   jetSource = cms.InputTag('slimmedJets'),                                                                                                                                                       
   labelName = 'UpdatedJEC',                                                                                                                                                                      
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual orrections (always set to 1)
)     


#declare producer for ecalBadCalibReducedMINIAODFilter
#https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
baddetEcallist = cms.vuint32(
  [872439604,872422825,872420274,872423218,
   872423215,872416066,872435036,872439336,
   872420273,872436907,872420147,872439731,
   872436657,872420397,872439732,872439339,
   872439603,872422436,872439861,872437051,
   872437052,872420649,872422436,872421950,
   872437185,872422564,872421566,872421695,
   872421955,872421567,872437184,872421951,
   872421694,872437056,872437057,872437313])

process.ecalBadCalibReducedMINIAODFilter = cms.EDFilter(
  "EcalBadCalibFilter",
  EcalRecHitSource = cms.InputTag("reducedEgamma:reducedEERecHits"),
  ecalMinEt        = cms.double(50.),
  baddetEcal       = baddetEcallist,
  taggingMode      = cms.bool(True),
  debug            = cms.bool(False)
  )

#rerun energy correction for electrons
setupEgammaPostRecoSeq(process,
    era='2018-Prompt')

runMetCorAndUncFromMiniAOD(process,
                           #isData=False,
                           isData=True,
                          )

#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

#runMetCorAndUncFromMiniAOD(process,
#    isData=False,
#    isData=True,
#    )

#runMetCorAndUncFromMiniAOD(process,
#    isData = False,
#    isData = True,
#    fixEE2017 = True,
#    fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
#    postfix = "ModifiedMET"
#    )

#process.mainNtuplizer.isMC = cms.bool(True)
process.mainNtuplizer.isMC = cms.bool(False)
#process.mainNtuplizer.dtag = cms.string("MC13TeV_TTToSemiLeptonic_2018")
process.mainNtuplizer.dtag = cms.string("DATA13TeV_SingleElectron_2018B")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(0)
process.mainNtuplizer.verbose = cms.bool(True)
process.mainNtuplizer.jetsTag = cms.InputTag('updatedPatJetsUpdatedJEC')
#process.mainNtuplizer.verbose = cms.bool(True)
#process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
#process.mainNtuplizer.verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.mainNtuplizer.jetsTag = cms.InputTag('slimmedMETTEST')

#process.mainNtuplizer.fatjetsTag = cms.InputTag('selectedPatJetsAK8PFCHS')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked')
process.mainNtuplizer.metFilterBitsTag = cms.InputTag('TriggerResults','','HLT')
#process.mainNtuplizer.metsTag = cms.InputTag("slimmedMETsModifiedMET","","bsmAnalysis")
#process.mainNtuplizer.bits = cms.InputTag('TriggerResultsTest','','HLT')
process.mainNtuplizer.Legacy2016 = cms.bool(False)

process.source = cms.Source("PoolSource",
             fileNames = cms.untracked.vstring("/store/data/Run2018C/EGamma/MINIAOD/17Sep2018-v1/70000/25C41660-D0D7-CF46-A544-D9DE8231BA59.root"),
#             fileNames = cms.untracked.vstring("/store/data/Run2018D/EGamma/MINIAOD/PromptReco-v2/000/325/175/00000/9D0F9360-DD60-314A-BB24-33D62A3CD6BD.root"),
#             fileNames = cms.untracked.vstring("/store/mc/RunIIAutumn18MiniAOD/TTToSemiLeptonic_TuneCP5up_13TeV-powheg-pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/90000/1137389F-C639-E641-B930-AF20537789CA.root"),
                            inputCommands=cms.untracked.vstring( 
        'keep *', 
        'drop *_ctppsLocalTrackLiteProducer_*_RECO'
#CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_RECO'
        )
)

process.TFileService = cms.Service("TFileService",
#			fileName = cms.string("analysis_MC_2J_woFJ.root")
			fileName = cms.string("analysis_data2018.root")
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")  

#jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addSoftDropSubjets=True,addNsub=True, Cut="pt>20") 
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addNsub=True, Cut="pt>20") 
#process.p = cms.Path(process.dump)
process.p = cms.Path( 
    process.ecalBadCalibReducedMINIAODFilter *
    process.patJetCorrFactorsUpdatedJEC *  
    process.updatedPatJetsUpdatedJEC *
#    process.fullPatMetSequenceTEST *
    process.fullPatMetSequence *
#    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.mainNtuplizer 
    ) #* process.dump )	
