import FWCore.ParameterSet.Config as cms

from UserCode.bsmhiggs_fwk.mainNtuplizer_cfi import *

#from JMEAnalysis.JetToolbox.jetToolbox_cff import jetToolbox
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.load("PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff")
process.load("Configuration.EventContent.EventContent_cff")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

#load run conditions
#process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')
#process.load('Configuration.Geometry.GeometryRecoDB_cff')
#process.load('Configuration.Geometry.GeometryReco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('RecoMET.METFilters.ecalBadCalibFilter_cfi')
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

#run ecalBadCalibReducedMINIAODFilter
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
     baddetEcal	      = baddetEcallist, 
     taggingMode      = cms.bool(True),
     debug	      = cms.bool(False)
     )
#rerun energy correction for electrons
setupEgammaPostRecoSeq(process,
    runVID=True, #saves CPU time by not needlessly re-running VID
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
#process.mainNtuplizer.dtag = cms.string("MC13TeV_Wh_amass40_2017")
process.mainNtuplizer.dtag = cms.string("MC13TeV_TTTo2L2Nu_powheg_2017")
#process.mainNtuplizer.dtag = cms.string("MC13TeV_TTGJets_2017_ext1")
#process.mainNtuplizer.dtag = cms.string("MC13TeV_TTToSemiLeptonic_powheg_2017")
#process.mainNtuplizer.dtag = cms.string("MC13TeV_TTToHadronic_powheg_2017")
#process.mainNtuplizer.dtag = cms.string("DATA13TeV_SingleElectron_2017B")
#process.mainNtuplizer.dtag = cms.string("DATA13TeV_SingleElectron_2017D")
#process.mainNtuplizer.dtag = cms.string("DATA13TeV_SingleElectron_2017F")
#process.mainNtuplizer.dtag = cms.string("Data13TeV_MuEG2017B")
#process.mainNtuplizer.dtag = cms.string("Data13TeV_MuEG2017C")
#process.mainNtuplizer.dtag = cms.string("Data13TeV_MuEG2017D")
#process.mainNtuplizer.dtag = cms.string("Data13TeV_MuEG2017F")
process.mainNtuplizer.xsec = cms.double(1.0)
process.mainNtuplizer.mctruthmode = cms.int32(1)
process.mainNtuplizer.verbose = cms.bool(False)
#process.mainNtuplizer.verbose = cms.bool(True)
#process.mainNtuplizer.jetsTag = cms.InputTag('selectedUpdatedPatJets')
#process.mainNtuplizer.verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices")
#process.mainNtuplizer.jetsTag = cms.InputTag('slimmedMETTEST')

#process.mainNtuplizer.fatjetsTag = cms.InputTag('selectedPatJetsAK8PFCHS')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('packedPatJetsAK8PFCHSSoftDrop')
#process.mainNtuplizer.fatjetsTag = cms.InputTag('slimmedJetsAK8PFPuppiSoftDropPacked')
process.mainNtuplizer.metFilterBitsTag = cms.InputTag('TriggerResults','','HLT')
process.mainNtuplizer.metsTag = cms.InputTag("slimmedMETsModifiedMET","","bsmAnalysis")
process.mainNtuplizer.objects = cms.InputTag("slimmedPatTrigger")
#process.mainNtuplizer.bits = cms.InputTag('TriggerResultsTest','','HLT')

process.source = cms.Source("PoolSource",
#	       	     fileNames = cms.untracked.vstring("/store/mc/RunIISummer16MiniAODv2/SUSYZHToAA_AATo4B_M-50_TuneCUETP8M1_13TeV_madgraph_pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v2/00000/BEEEFB4E-37C8-E811-A98B-001E6757E05C.root"),
#              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/NMSSM_WHToAATo4b_M-125_M-12_madgraph_pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/40000/00F6F5BE-8742-E811-90B6-B499BAABF212.root"),
              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/TTTo2L2Nu_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/FECD59BD-1842-E811-96D7-0242AC130002.root"),
#              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/TTGJets_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/120000/FEED6629-6DD6-E811-8C28-E0071B745DC0.root"),
#              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/30000/ECA4210B-F058-E811-8EBF-509A4C74915C.root"),
#              fileNames = cms.untracked.vstring("/store/mc/RunIIFall17MiniAODv2/TTToHadronic_TuneCP5_PSweights_13TeV-powheg-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/00000/289E3510-5942-E811-9FDF-E0071B6CAD20.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017B/SingleElectron/MINIAOD/31Mar2018-v1/30000/04B05308-0038-E811-99AB-008CFAC94314.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017D/SingleElectron/MINIAOD/31Mar2018-v1/90000/4C22E5DC-1A39-E811-9DD3-24BE05C6E561.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017F/SingleElectron/MINIAOD/31Mar2018-v1/010000/72A98F53-E53B-E811-8945-0CC47A4C8E2E.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017B/MuonEG/MINIAOD/31Mar2018-v1/610000/86B6E7A3-8739-E811-96AF-002590DE6E52.root"),
#               fileNames = cms.untracked.vstring("/store/data/Run2017C/MuonEG/MINIAOD/31Mar2018-v1/80000/FCABD91C-5E37-E811-8379-0CC47A5FBDC5.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017D/MuonEG/MINIAOD/31Mar2018-v1/100000/FE93310F-E837-E811-9CBE-B496910A8618.root"),
#              fileNames = cms.untracked.vstring("/store/data/Run2017F/MuonEG/MINIAOD/31Mar2018-v1/80000/F632D544-F636-E811-A63B-FA163E5D77FE.root"),
                            inputCommands=cms.untracked.vstring( 
        'keep *', 
        'drop *_ctppsLocalTrackLiteProducer_*_RECO'
#CTPPSPixelClusteredmDetSetVector_ctppsPixelClusters_*_RECO'
        )
)

process.TFileService = cms.Service("TFileService",
#			fileName = cms.string("analysis_DATA_RunB.root")
			fileName = cms.string("analysisMC_TTTo2L2Nu.root")
)

process.dump = cms.EDAnalyzer("EventContentAnalyzer")  

#jetToolbox(process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True,addNsubSubjets=True,addSoftDropSubjets=True,addNsub=True,Cut="pt>20")
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addSoftDropSubjets=True,addNsub=True, Cut="pt>20") 
#jetToolbox( process, 'ak8', 'jetSequence', 'out', PUMethod='CHS', addPruning=True, addSoftDrop=True, addNsub=True, Cut="pt>20") 
#process.p = cms.Path(process.dump)
process.p = cms.Path( 
#    process.fullPatMetSequenceTEST *
    process.ecalBadCalibReducedMINIAODFilter *
    process.fullPatMetSequence *
    process.fullPatMetSequenceModifiedMET *
    process.egammaPostRecoSeq *
    process.mainNtuplizer 
    ) # * process.dump )	
