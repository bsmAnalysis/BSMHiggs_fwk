import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.coreTools import runOnData
from PhysicsTools.PatAlgos.tools.jetTools import supportedJetAlgos, addJetCollection, updateJetCollection
from PhysicsTools.PatAlgos.tools.helpers  import getPatAlgosToolsTask, addToProcessAndTask

def puppiAK8ReclusterFromMiniAOD(process, 
    runOnMC, 
    useExistingWeights, 
    btagDiscriminatorsAK8=None, 
    btagDiscriminatorsAK8Subjets=None,
    reclusterAK8GenJets=False
  ):
  task = getPatAlgosToolsTask(process)

  pfLabel = "packedPFCandidates"
  pvLabel = "offlineSlimmedPrimaryVertices"
  svLabel = "slimmedSecondaryVertices"
  muLabel = "slimmedMuons"
  elLabel = "slimmedElectrons"
  gpLabel = "prunedGenParticles"

  genJetsAK8Collection = cms.InputTag("slimmedGenJetsAK8")
  genSubJetsForAK8Collection = cms.InputTag("slimmedGenJetsAK8SoftDropSubJets")

  JETCorrLevels = ["L2Relative", "L3Absolute", "L2L3Residual"]

  if runOnMC and reclusterAK8GenJets:
    addToProcessAndTask("packedGenParticlesForJetsNoNu", cms.EDFilter("CandPtrSelector",
        src = cms.InputTag("packedGenParticles"),
        cut = cms.string("abs(pdgId) != 12 && abs(pdgId) != 14 && abs(pdgId) != 16"),
      ),
      process, task
    )
    from RecoJets.JetProducers.ak8GenJets_cfi import ak8GenJets, ak8GenJetsConstituents, ak8GenJetsSoftDrop
    from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJetsAK8
    addToProcessAndTask('ak8GenJetsNoNu', ak8GenJets.clone(src='packedGenParticlesForJetsNoNu'),process, task)
    addToProcessAndTask('ak8GenJetsNoNuConstituents', ak8GenJetsConstituents.clone(src='ak8GenJetsNoNu'), process, task )
    addToProcessAndTask('ak8GenJetsNoNuSoftDrop',ak8GenJetsSoftDrop.clone(src=cms.InputTag('ak8GenJetsNoNuConstituents', 'constituents')),process,task)
    genJetsAK8Collection=cms.InputTag("ak8GenJetsNoNu")
    genSubJetsForAK8Collection=cms.InputTag("ak8GenJetsNoNuSoftDrop","SubJets")

  #########################
  #
  # Setup puppi weights
  #
  ########################
  process.load('CommonTools.PileupAlgos.Puppi_cff')
  task.add(process.puppi)
  #
  process.puppi.candName = pfLabel
  process.puppi.vertexName = pvLabel
  process.puppi.clonePackedCands = True
  process.puppi.useExistingWeights = useExistingWeights # If true, use the weights in MiniAOD.
  #
  # Just in case we choose to recompute puppi, ensure PuppiProducer is setup for PUPPI tune V15. Taken from here:
  # https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_26/CommonTools/PileupAlgos/python/customizePuppiTune_cff.py#L28-L40
  process.puppi.PtMaxCharged = 20.
  process.puppi.EtaMinUseDeltaZ = 2.4
  process.puppi.PtMaxNeutralsStartSlope = 20.
  process.puppi.NumOfPUVtxsForCharged = 2
  process.puppi.algos[0].etaMin = [-0.01]


  ########################
  #
  # AK8 Puppi jets
  #
  ########################
  #
  # Recluster jets and do soft-drop grooming
  #
  process.load("RecoJets.JetProducers.ak8PFJets_cfi")
  task.add(process.ak8PFJetsPuppi)
  task.add(process.ak8PFJetsPuppiSoftDrop)

  process.ak8PFJetsPuppi.src = "puppi"

  # AK8 jet constituents for softdrop
  addToProcessAndTask("ak8PFJetsPuppiConstituents", cms.EDProducer("MiniAODJetConstituentSelector",
      src = cms.InputTag("ak8PFJetsPuppi"),
      cut = cms.string("pt > 100.0 && abs(rapidity()) < 2.4")
    ),
    process, task
  )

  # Soft-drop grooming
  process.ak8PFJetsPuppiSoftDrop.src = "ak8PFJetsPuppiConstituents:constituents"

  # Soft-drop mass
  process.load("RecoJets.JetProducers.ak8PFJetsPuppi_groomingValueMaps_cfi")
  task.add(process.ak8PFJetsPuppiSoftDropMass)
  process.ak8PFJetsPuppiSoftDropMass.src = "ak8PFJetsPuppi"
  process.ak8PFJetsPuppiSoftDropMass.matched = "ak8PFJetsPuppiSoftDrop"
  #=============================================
  #
  # PATify
  #
  #=============================================
  #
  # AK8 jets
  #
  addJetCollection(
    process,
    labelName          = "AK8Puppi",
    jetSource          = cms.InputTag("ak8PFJetsPuppi"),
    algo               = "ak",
    rParam             = 0.8,
    pfCandidates       = cms.InputTag(pfLabel),
    pvSource           = cms.InputTag(pvLabel),
    svSource           = cms.InputTag(svLabel),
    muSource           = cms.InputTag(muLabel),
    elSource           = cms.InputTag(elLabel),
    genJetCollection   = genJetsAK8Collection,
    genParticles       = cms.InputTag(gpLabel),
    jetCorrections     = ("AK8PFPuppi", cms.vstring(JETCorrLevels), "None"),
    btagDiscriminators = ([
        "pfCombinedSecondaryVertexV2BJetTags",
        "pfCombinedInclusiveSecondaryVertexV2BJetTags",
        "pfCombinedMVAV2BJetTags",
        "pfDeepCSVJetTags:probb",
        "pfDeepCSVJetTags:probc",
        "pfDeepCSVJetTags:probudsg",
        "pfDeepCSVJetTags:probbb",
        "pfBoostedDoubleSecondaryVertexAK8BJetTags"
      ]
    ),
  )
  process.patJetsAK8Puppi.userData.userFloats.src = [] # start with empty list of user floats
  process.patJetsAK8Puppi.userData.userFloats.src += ["ak8PFJetsPuppiSoftDropMass"]
  process.patJetsAK8Puppi.addTagInfos = cms.bool(False)

  process.selectedPatJetsAK8Puppi.cut = cms.string("pt > 100")
  process.selectedPatJetsAK8Puppi.cutLoose = cms.string("pt > 30")
  process.selectedPatJetsAK8Puppi.nLoose = cms.uint32(3)

  #
  # Add AK8 Njetiness
  #
  from RecoJets.JetProducers.nJettinessAdder_cfi import Njettiness
  addToProcessAndTask("NjettinessAK8Puppi", Njettiness.clone(
      src = "ak8PFJetsPuppi",
    ),
    process, task
  )
  process.patJetsAK8Puppi.userData.userFloats.src += [
    "NjettinessAK8Puppi:tau1",
    "NjettinessAK8Puppi:tau2",
    "NjettinessAK8Puppi:tau3",
    "NjettinessAK8Puppi:tau4"
  ]

  #
  # AK8 soft-drop jets
  #
  addJetCollection(
    process,
    labelName = "AK8PFPuppiSoftDrop",
    jetSource = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
    btagDiscriminators = ["None"],
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    genJetCollection = genJetsAK8Collection,
    genParticles = cms.InputTag(gpLabel),
    jetCorrections = ("AK8PFPuppi", cms.vstring(JETCorrLevels), "None"),
    getJetMCFlavour = False # jet flavor disabled regardless if running on MC or data
  )

  #
  # Soft-drop subjets
  #
  addJetCollection(
    process,
    labelName = "AK8PFPuppiSoftDropSubjets",
    jetSource = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
    algo = "ak",  # needed for subjet flavor clustering
    rParam = 0.8, # needed for subjet flavor clustering
    explicitJTA = True,  # needed for subjet b tagging
    svClustering = True, # needed for subjet b tagging
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    genJetCollection = genSubJetsForAK8Collection,
    genParticles = cms.InputTag(gpLabel),
    fatJets = cms.InputTag("ak8PFJetsPuppi"),               # needed for subjet flavor clustering
    groomedFatJets = cms.InputTag("ak8PFJetsPuppiSoftDrop"), # needed for subjet flavor clustering
    jetCorrections = ("AK4PFPuppi", cms.vstring(JETCorrLevels), "None"),
    btagDiscriminators = btagDiscriminatorsAK8Subjets.names.value() if btagDiscriminatorsAK8Subjets is not None else ['None'],
  )


  #=============================================
  #
  #
  #
  #=============================================
  #
  # add groomed ECFs and N-subjettiness to soft dropped pat::Jets for fat jets and subjets
  #
  process.load('RecoJets.JetProducers.ECF_cff')

  addToProcessAndTask('nb1AK8PuppiSoftDrop', process.ecfNbeta1.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
      cuts = cms.vstring('', '', 'pt > 250')
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDrop.userData.userFloats.src += [
    'nb1AK8PuppiSoftDrop:ecfN2',
    'nb1AK8PuppiSoftDrop:ecfN3',
  ]

  addToProcessAndTask('nb2AK8PuppiSoftDrop', process.ecfNbeta2.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop"),
      cuts = cms.vstring('', '', 'pt > 250')
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDrop.userData.userFloats.src += [
    'nb2AK8PuppiSoftDrop:ecfN2',
    'nb2AK8PuppiSoftDrop:ecfN3',
  ]

  #
  # add groomed ECFs and N-subjettiness to soft drop subjets
  #
  addToProcessAndTask("nb1AK8PuppiSoftDropSubjets", process.ecfNbeta1.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
    ),
    process, task
  )

  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    'nb1AK8PuppiSoftDropSubjets:ecfN2',
    'nb1AK8PuppiSoftDropSubjets:ecfN3'
  ]

  addToProcessAndTask("nb2AK8PuppiSoftDropSubjets", process.ecfNbeta2.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
    ),
    process, task
  )

  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    'nb2AK8PuppiSoftDropSubjets:ecfN2',
    'nb2AK8PuppiSoftDropSubjets:ecfN3'
  ]

  addToProcessAndTask("NjettinessAK8Subjets", Njettiness.clone(
      src = cms.InputTag("ak8PFJetsPuppiSoftDrop", "SubJets"),
    ),
    process, task
  )
  process.patJetsAK8PFPuppiSoftDropSubjets.userData.userFloats.src += [
    "NjettinessAK8Subjets:tau1",
    "NjettinessAK8Subjets:tau2",
    "NjettinessAK8Subjets:tau3",
    "NjettinessAK8Subjets:tau4",
  ]

  addToProcessAndTask("slimmedJetsAK8PFPuppiSoftDropSubjets", cms.EDProducer("PATJetSlimmer",
      src = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDropSubjets"),
      packedPFCandidates = cms.InputTag(pfLabel),
      dropJetVars = cms.string("1"),
      dropDaughters = cms.string("0"),
      rekeyDaughters = cms.string("0"),
      dropTrackRefs = cms.string("1"),
      dropSpecific = cms.string("1"),
      dropTagInfos = cms.string("1"),
      modifyJets = cms.bool(True),
      mixedDaughters = cms.bool(False),
      modifierConfig = cms.PSet( modifications = cms.VPSet() )
    ),
    process, task
  )

  ## Establish references between PATified fat jets and subjets using the BoostedJetMerger
  ## Take subjets and put it in the grommed jet
  addToProcessAndTask("slimmedJetsAK8PFPuppiSoftDropPacked", cms.EDProducer("BoostedJetMerger",
      jetSrc    = cms.InputTag("selectedPatJetsAK8PFPuppiSoftDrop"),
      subjetSrc = cms.InputTag("slimmedJetsAK8PFPuppiSoftDropSubjets")
    ),
    process, task
  )

  addToProcessAndTask("packedPatJetsAK8", cms.EDProducer("JetSubstructurePacker",
      jetSrc = cms.InputTag("selectedPatJetsAK8Puppi"),
      distMax = cms.double(0.8),
      algoTags = cms.VInputTag(
        cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked")
      ),
      algoLabels = cms.vstring(
        'SoftDropPuppi'
      ),
      fixDaughters = cms.bool(False),
      packedPFCandidates = cms.InputTag(pfLabel),
    ),
    process, task
  )

  #=============================================
  #
  # Update the selectedPatJet collection.
  # This is where we setup
  # -  JEC
  # -  b-tagging discriminators
  #
  #=============================================
  from PhysicsTools.PatAlgos.slimming.slimmedJets_cfi import slimmedJetsAK8
  addToProcessAndTask("slimmedJetsAK8NoDeepTags", slimmedJetsAK8.clone(rekeyDaughters = "0"), process, task)
  # Reconfigure the slimmedAK8 jet information to keep
  process.slimmedJetsAK8NoDeepTags.dropDaughters = cms.string("pt < 170")
  process.slimmedJetsAK8NoDeepTags.dropSpecific = cms.string("pt < 170")
  process.slimmedJetsAK8NoDeepTags.dropTagInfos = cms.string("pt < 170")


  updateJetCollection(
    process,
    jetSource = cms.InputTag("slimmedJetsAK8NoDeepTags"),
    # updateJetCollection defaults to MiniAOD inputs but
    # here it is made explicit (as in training or MINIAOD redoing)
    pfCandidates = cms.InputTag(pfLabel),
    pvSource = cms.InputTag(pvLabel),
    svSource = cms.InputTag(svLabel),
    muSource = cms.InputTag(muLabel),
    elSource = cms.InputTag(elLabel),
    rParam = 0.8,
    jetCorrections = ('AK8PFPuppi', cms.vstring(JETCorrLevels), 'None'),
    btagDiscriminators = btagDiscriminatorsAK8.names.value() if btagDiscriminatorsAK8 is not None else ['None'],
    postfix = "SlimmedAK8DeepTags",
    printWarning = False
  )

  addToProcessAndTask("slimmedJetsAK8", process.selectedUpdatedPatJetsSlimmedAK8DeepTags.clone(), process, task)
  del process.selectedUpdatedPatJetsSlimmedAK8DeepTags

  ########################
  #
  # Modify JECs when processing real Data
  # Disable any MC-only features.
  #
  ########################
  if not(runOnMC):
    runOnData(process, names=["Jets"], outputModules = [])

  return process

