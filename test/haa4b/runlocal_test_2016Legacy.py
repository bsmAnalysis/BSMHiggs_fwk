import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

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
#process.GlobalTag.globaltag = '94X_mcRun2_asymptotic_v3'
process.GlobalTag.globaltag = '94X_dataRun2_v10'

process.options   = cms.untracked.PSet(
   allowUnscheduled = cms.untracked.bool(True)
)

runMetCorAndUncFromMiniAOD(process,  
#			isData=False
			isData=True
			)		

#process.mainNtuplizer.isMC = cms.bool(True)
process.mainNtuplizer.isMC = cms.bool(False)
#process.mainNtuplizer.dtag = cms.string("MC13TeV_TTJets_powheg_2016")
process.mainNtuplizer.dtag = cms.string("Data13TeV_MuEG2016B_ver2")
process.mainNtuplizer.xsec = cms.double(1)
process.mainNtuplizer.mctruthmode = cms.int32(1)
process.mainNtuplizer.verbose = cms.bool(False)
process.mainNtuplizer.metFilterBitsTag = cms.InputTag('TriggerResults','','HLT')
process.mainNtuplizer.Legacy2016 = cms.bool(True)
#process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')

process.source = cms.Source("PoolSource",
#	       	     fileNames = cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv3/TT_TuneCUETP8M2T4_13TeV-powheg-pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/00000/001BCE65-30C3-E811-A357-A4BF0112DD3C.root"),
	       	     fileNames = cms.untracked.vstring("/store/data/Run2016B/MuonEG/MINIAOD/17Jul2018_ver2-v1/20000/FAC231DB-E28B-E811-A5F3-0025905C53F2.root"),
		     inputCommands=cms.untracked.vstring(
			   'keep *', 
 		     	   'drop *_ctppsLocalTrackLiteProducer_*_RECO' 
			   )        
)

#jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")

process.TFileService = cms.Service("TFileService",
			fileName = cms.string("analysis_data.root")
)
                     		 				             
process.p = cms.Path( process.fullPatMetSequence * process.mainNtuplizer )	

