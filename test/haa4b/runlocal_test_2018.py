import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")

#load run conditions
#process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')

updateJetCollection(
   process,
   jetSource = cms.InputTag('slimmedJets'),
   labelName = 'UpdatedJEC',
   jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None')  # Update: Safe to always add 'L2L3Residual' as MC contains dummy L2L3Residual orrections (always set to 1)
)

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from RecoJets.JetProducers.PileupJetID_cfi import _chsalgos_106X_UL18
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId
process.pileupJetIdUpdated = pileupJetId.clone( 
        jets=cms.InputTag('updatedPatJetsUpdatedJEC'), #Your JEC corrected jets here
        inputIsCorrected=True,
        applyJec=False,
        vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
        algos = cms.VPSet(_chsalgos_106X_UL18),
    )

#rerun energy correction for electrons
#Recipe for running scales and smearings using EgammaPostRecoTools
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=False, #saves CPU time by not needlessly re-running VID, if you want the Fall17V2 IDs, set this to True or remove (default is True)
                       era='2018-UL')    
#a sequence egammaPostRecoSeq has now been created and should be added to your path, eg process.p=cms.Path(process.egammaPostRecoSeq)

runMetCorAndUncFromMiniAOD(process,
    isData=True
    )

process.mainNtuplizer.isMC = cms.bool(False)
process.mainNtuplizer.dtag = cms.string("DATA13TeV_SingleElectron_2018B")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(0)
process.mainNtuplizer.verbose = cms.bool(True)
process.mainNtuplizer.metFilterBitsTag = cms.InputTag('TriggerResults','','HLT')
process.mainNtuplizer.Legacy2016 = cms.bool(False)
process.mainNtuplizer.jetsTag = cms.InputTag('updatedPatJetsUpdatedJEC')

#process.mainNtuplizer.verbose = cms.bool(True)
#process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
#process.mainNtuplizer.verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.mainNtuplizer.jetsTag = cms.InputTag('slimmedMETTEST')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('selectedPatJetsAK8PFCHS')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked')
#process.mainNtuplizer.metsTag = cms.InputTag("slimmedMETsModifiedMET","","bsmAnalysis")
#process.mainNtuplizer.bits = cms.InputTag('TriggerResultsTest','','HLT')

###########################################################
## Workflow START
## We try to reproduce slimmedJetsAK8 collection
## like in MiniAOD production but lower the pt cut to as low
## as possible
##
###########################################################
from UserCode.bsmhiggs_fwk.puppiJetMETReclusteringTools import puppiAK8ReclusterFromMiniAOD

from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import _pfDeepBoostedJetTagsAll as pfDeepBoostedJetTagsAll
from RecoBTag.ONNXRuntime.pfHiggsInteractionNet_cff import _pfHiggsInteractionNetTagsProbs as pfHiggsInteractionNetTagsProbs
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetJetTagsAll as pfParticleNetJetTagsAll
from RecoBTag.ONNXRuntime.pfParticleNet_cff import _pfParticleNetMassRegressionOutputs

btagDiscriminatorsAK8 = cms.PSet(names = cms.vstring(
'pfCombinedSecondaryVertexV2BJetTags',
'pfCombinedInclusiveSecondaryVertexV2BJetTags',
# 'pfCombinedMVAV2BJetTags',
'pfDeepCSVJetTags:probb',
'pfDeepCSVJetTags:probc',
'pfDeepCSVJetTags:probudsg',
'pfDeepCSVJetTags:probbb',
'pfBoostedDoubleSecondaryVertexAK8BJetTags',
'pfMassIndependentDeepDoubleBvLV2JetTags:probQCD',
'pfMassIndependentDeepDoubleBvLV2JetTags:probHbb',
'pfMassIndependentDeepDoubleCvLV2JetTags:probQCD',
'pfMassIndependentDeepDoubleCvLV2JetTags:probHcc',
'pfMassIndependentDeepDoubleCvBV2JetTags:probHbb',
'pfMassIndependentDeepDoubleCvBV2JetTags:probHcc',
)
# + pfDeepBoostedJetTagsAll
+ pfParticleNetJetTagsAll
# + pfHiggsInteractionNetTagsProbs
# # + _pfParticleNetMassRegressionOutputs
)

btagDiscriminatorsAK8Subjets = cms.PSet(names = cms.vstring(
   'pfDeepCSVJetTags:probb',
   'pfDeepCSVJetTags:probbb',
 )
)

process = puppiAK8ReclusterFromMiniAOD(process,
   runOnMC=False,
   useExistingWeights=True,
   btagDiscriminatorsAK8=btagDiscriminatorsAK8,
   btagDiscriminatorsAK8Subjets=btagDiscriminatorsAK8Subjets,
   reclusterAK8GenJets=False
)

process.patAlgosToolsTask.add(process.patJetPartons)

################################
## Override some configurations
## for reclustered AK8 Gen jets
################################
#process.ak8GenJetsNoNu.jetPtMin = 10
#process.ak8GenJetsNoNuSoftDrop.jetPtMin = 10
#process.ak8GenJetsNoNuConstituents.cut = "pt > 10"

#################################
## Override some configurations
## for reclustered AK8 Puppi jets
##################################
process.ak8PFJetsPuppi.jetPtMin = 15
process.ak8PFJetsPuppiSoftDrop.jetPtMin = 15
process.ak8PFJetsPuppiConstituents.cut = "pt > 15. && abs(rapidity()) < 2.4"

finalAK8PuppiPt = 30
process.selectedPatJetsAK8Puppi.cut = "pt > {}".format(finalAK8PuppiPt)
process.selectedPatJetsAK8Puppi.cutLoose = ""
process.selectedPatJetsAK8Puppi.nLoose = 0
process.slimmedJetsAK8NoDeepTags.dropDaughters = cms.string("pt < {}".format(finalAK8PuppiPt))
process.slimmedJetsAK8NoDeepTags.dropSpecific = cms.string("pt < {}".format(finalAK8PuppiPt))
process.slimmedJetsAK8NoDeepTags.dropTagInfos = cms.string("pt < {}".format(finalAK8PuppiPt))

##################################
## For reclustered AK8 Puppi jets
##################################
process.mainNtuplizer.fatjetsTag = "slimmedJetsAK8"

#######################################################
##
## Workflow END
##
########################################################

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

#------ Declare the correct global tag ------#
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
process.GlobalTag.globaltag = '102X_dataRun2_Sep2018ABC_v2'

process.options   = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True),
   #SkipEvent = cms.untracked.vstring('ProductNotFound')
)

process.source = cms.Source("PoolSource",
             fileNames = cms.untracked.vstring("/store/data/Run2018C/EGamma/MINIAOD/17Sep2018-v1/70000/25C41660-D0D7-CF46-A544-D9DE8231BA59.root"),
             inputCommands=cms.untracked.vstring( 
    'keep *', 
        'drop *_ctppsLocalTrackLiteProducer_*_RECO'
#CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_RECO'
        )
)

process.TFileService = cms.Service("TFileService",
			fileName = cms.string("analysis_data2018.root")
)

process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)

#jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addSoftDropSubjets=True,addNsub=True, Cut="pt>20") 
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addNsub=True, Cut="pt>20") 
#process.p = cms.Path(process.dump)
process.p = cms.Path( 
    process.patJetCorrFactorsUpdatedJEC *  
    process.updatedPatJetsUpdatedJEC *
    process.pileupJetIdUpdated *
    process.fullPatMetSequence *
    process.egammaPostRecoSeq *
    process.mainNtuplizer 
    )	


