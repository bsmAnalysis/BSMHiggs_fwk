import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 5

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#load run conditions
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

#------ Declare the correct global tag ------#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
#process.GlobalTag.globaltag = '80X_dataRun2_2016SeptRepro_v7'
process.GlobalTag.globaltag = '80X_mcRun2_asymptotic_2016_TrancheIV_v6'

process.options   = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True)
)
## b-tag discriminators
bTagDiscriminators = [
  'deepFlavourJetTags:probb',
  'deepFlavourJetTags:probbb',
]
from PhysicsTools.PatAlgos.tools.jetTools import *
## Update the slimmedJets in miniAOD: corrections from the chosen Global Tag are applied and the b-tag discriminators are re-evaluated
updateJetCollection(
    process,
    jetSource = cms.InputTag('slimmedJets'),
    jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute']), 'None'),
    btagDiscriminators = bTagDiscriminators
)

runMetCorAndUncFromMiniAOD(process,
                           isData=False,
                           )
process.mainNtuplizer.isMC = cms.bool(True)
process.mainNtuplizer.dtag = cms.string("MC13TeV_SUSYZHToAA_AATo4B_M-50_2016")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(0)
process.mainNtuplizer.verbose = cms.bool(False)
process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')

process.source = cms.Source("PoolSource",
	       	     fileNames = cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv2/SUSYZHToAA_AATo4B_M-50_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/00000/BEEEFB4E-37C8-E811-A98B-001E6757E05C.root"),
                            inputCommands=cms.untracked.vstring( 
        'keep *', 
        'drop *_ctppsLocalTrackLiteProducer_*_RECO'
#CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_RECO'
        )
)

process.TFileService = cms.Service("TFileService",
			fileName = cms.string("analysis.root")
)

jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")

process.dump = cms.EDAnalyzer("EventContentAnalyzer")  
                     		 				             
process.p = cms.Path( process.fullPatMetSequence * process.mainNtuplizer ) #* process.dump )	

