import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#load run conditions
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
#process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'
process.GlobalTag.globaltag = '94X_mc2017_realistic_v14'
#process.GlobalTag.globaltag = '94X_dataRun2_v6'

process.options   = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),
   SkipEvent = cms.untracked.vstring('ProductNotFound')
)

#rerun energy correction for electrons
setupEgammaPostRecoSeq(process,
    runVID=False, #saves CPU time by not needlessly re-running VID
    era='2017-Nov17ReReco')
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

# to recluster both jets and MET and then recorrect and get the proper uncertainties
# this example adds a postfix in the slimmedMET collection: postfix="TEST"
#runMetCorAndUncFromMiniAOD(process,
#    isData=False,
#    pfCandColl=cms.InputTag("packedPFCandidates"),
#    recoMetFromPFCs=True,
#    CHS = True, #This is an important step and determines what type of jets to be reclustered
#    reclusterJets = True,
#    postfix="TEST"
#    )

runMetCorAndUncFromMiniAOD(process,
    isData=False,
#    isData=True,
    )

runMetCorAndUncFromMiniAOD(process,
    isData = False,
#    isData = True,
    fixEE2017 = True,
    fixEE2017Params = {'userawPt': True, 'ptThreshold':50.0, 'minEtaThreshold':2.65, 'maxEtaThreshold': 3.139} ,
    postfix = "ModifiedMET"
    )

process.mainNtuplizer.isMC = cms.bool(True)
#process.mainNtuplizer.isMC = cms.bool(False)
process.mainNtuplizer.dtag = cms.string("MC13TeV_SUSYZHToAA_AATo4B_M-50_2016")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(0)
process.mainNtuplizer.verbose = cms.bool(False)
#process.mainNtuplizer.verbose = cms.bool(True)
#process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
#process.mainNtuplizer.verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.mainNtuplizer.jetsTag = cms.InputTag('slimmedMETTEST')

#process.mainNtuplizer.fatjetsTag = cms.InputTag('selectedPatJetsAK8PFCHS')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked')
process.mainNtuplizer.metFilterBitsTag = cms.InputTag('TriggerResults','','HLT')
#process.mainNtuplizer.bits = cms.InputTag('TriggerResultsTest','','HLT')

process.source = cms.Source("PoolSource",
#	       	     fileNames = cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv2/SUSYZHToAA_AATo4B_M-50_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/00000/BEEEFB4E-37C8-E811-A98B-001E6757E05C.root"),
#              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/NMSSM_WHToAATo4b_M-125_M-12_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/00F6F5BE-8742-E811-90B6-B499BAABF212.root"),
             fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/NMSSM_WHToAATo4b_M-125_M-40_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/50FCC250-3143-E811-A2B0-AC1F6B1AF03C.root"),
 #             fileNames = cms.untracked.vstring("/store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/04B05308-0038-E811-99AB-008CFAC94314.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/010000/72A98F53-E53B-E811-8945-0CC47A4C8E2E.root"),
                            inputCommands=cms.untracked.vstring( 
        'keep *', 
        'drop *_ctppsLocalTrackLiteProducer_*_RECO'
#CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_RECO'
        )
)

process.TFileService = cms.Service("TFileService",
#			fileName = cms.string("analysis_MC_2J_woFJ.root")
			fileName = cms.string("analysis_MC.root")
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")  

#jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addSoftDropSubjets=True,addNsub=True, Cut="pt>20") 
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addNsub=True, Cut="pt>20") 
#process.p = cms.Path(process.dump)
process.p = cms.Path( 
#    process.fullPatMetSequenceTEST *
    process.fullPatMetSequence *
    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.mainNtuplizer 
    ) #* process.dump )	
