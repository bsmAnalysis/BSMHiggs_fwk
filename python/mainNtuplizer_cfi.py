import FWCore.ParameterSet.Config as cms

from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import correctionType
from EgammaAnalysis.ElectronTools.calibrationTablesRun2 import files

process = cms.Process("bsmAnalysis")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 5000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),
                                        SkipEvent = cms.untracked.vstring('ProductNotFound')
                                        )

process.mainNtuplizer = cms.EDAnalyzer('mainNtuplizer',

    lheInfo = cms.InputTag("externalLHEProducer"),

    beamSpotTag = cms.InputTag("offlineBeamSpot"),
    verticesTag = cms.InputTag("offlineSlimmedPrimaryVertices"),

    rhoAll = cms.InputTag("fixedGridRhoAll"),
    rhoFastjetAll = cms.InputTag("fixedGridRhoFastjetAll"),
    rhoFastjetAllCalo = cms.InputTag("fixedGridRhoFastjetAllCalo"),
    rhoFastjetCentralCalo = cms.InputTag("fixedGridRhoFastjetCentralCalo"),
    rhoFastjetCentralChargedPileUp = cms.InputTag("fixedGridRhoFastjetCentralChargedPileUp"),
    rhoFastjetCentralNeutral = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),

    muonsTag = cms.InputTag("slimmedMuons"),

    electronsTag = cms.InputTag("slimmedElectrons"),

    tausTag = cms.InputTag("slimmedTaus"),

    photonsTag = cms.InputTag("slimmedPhotons"),

    jetsTag = cms.InputTag("slimmedJets"),
    jetsPuppiTag = cms.InputTag("slimmedJetsPuppi"),
    fatjetsTag = cms.InputTag("slimmedJetsAK8"),

    metsTag = cms.InputTag("slimmedMETs","","bsmAnalysis"),
    metsTagData = cms.InputTag("slimmedMETsMuEGClean"),                                        
    metsNoHFTag = cms.InputTag("slimmedMETsNoHF"),
    metsPuppiTag = cms.InputTag("slimmedMETsPuppi"),

    svTag = cms.InputTag("slimmedSecondaryVertices"),                                   

    metFilterBitsTag = cms.InputTag("TriggerResults"),
    packedTag = cms.InputTag("packedGenParticles"),
    prunedTag = cms.InputTag("prunedGenParticles"),
    genJetsTag = cms.InputTag("slimmedGenJets"),

    puInfoTag = cms.InputTag("slimmedAddPileupInfo", "", "PAT"),
    genInfoTag = cms.InputTag("generator", "", "SIM"),

    correctionFile = cms.string(files[correctionType]),                                       
    reducedEcalRecHitsEB = cms.InputTag('reducedEgamma:reducedEBRecHits'),
    reducedEcalRecHitsEE = cms.InputTag('reducedEgamma:reducedEERecHits'),                                       

    ##trigger
    bits = cms.InputTag("TriggerResults","","HLT"),
    prescales = cms.InputTag("patTrigger"),
    objects = cms.InputTag("selectedPatTrigger"),

        DoubleMuTrigs = cms.vstring("HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_v",
					"HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ_v"
				   ),
    	DoubleEleTrigs = cms.vstring("HLT_Ele17_Ele12_CaloIdL_TrackIdL_IsoVL_DZ_v",
					),
    	SingleMuTrigs = cms.vstring("HLT_IsoMu22_v", "HLT_IsoMu24_v",
					"HLT_IsoTkMu22_v", "HLT_IsoTkMu24_v"),
    	SingleEleTrigs = cms.vstring("HLT_Ele27_eta2p1_WPLoose_Gsf_v",
					"HLT_Ele27_WPTight_Gsf_v"
					),
   	MuEGTrigs = cms.vstring("HLT_Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL_v",
					"HLT_Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_v"),

)

