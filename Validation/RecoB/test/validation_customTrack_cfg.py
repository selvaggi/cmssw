# The following comments couldn't be translated into the new config version:
#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms
process = cms.Process("validation")

import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing ('analysis')

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

options.register ('jets',
                  "ak4PFJetsCHS", # default value, examples : "ak4PFJets", "ak4PFJetsCHS"
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.string,  
                  "jet collection to use")

options.parseArguments()

whichJets  = options.jets 
applyJEC = True
corrLabel = "ak4PFCHS"
from Configuration.AlCa.GlobalTag import GlobalTag
tag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
useTrigger = False
triggerPath = "HLT_PFJet80_v*"
runOnMC    = True
#Flavour plots for MC: "all" = plots for all jets ; "dusg" = plots for d, u, s, dus, g independently ; not mandatory and any combinations are possible                                     
#b, c, light (dusg), non-identified (NI), PU jets plots are always produced
flavPlots = "allbcldusg"

###prints###
print "jet collcetion asked : ", whichJets
print "JEC applied?", applyJEC, ", correction:", corrLabel 
print "trigger will be used ? : ", useTrigger, ", Trigger paths:", triggerPath
print "is it MC ? : ", runOnMC, ", Flavours:", flavPlots
print "Global Tag : ", tag.globaltag
############

process.load("DQMServices.Components.DQMEnvironment_cfi")
process.load("DQMServices.Core.DQM_cfg")

process.load("JetMETCorrections.Configuration.JetCorrectors_cff")
process.load("CommonTools.ParticleFlow.goodOfflinePrimaryVertices_cfi")
process.load("RecoJets.JetAssociationProducers.ak4JTA_cff")
process.load("RecoBTag.Configuration.RecoBTag_cff")
process.load("PhysicsTools.JetMCAlgos.HadronAndPartonSelector_cfi")
process.load("PhysicsTools.JetMCAlgos.AK4PFJetsMCFlavourInfos_cfi")
process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
process.JECseq = cms.Sequence(getattr(process,corrLabel+"L1FastL2L3CorrectorChain"))

newjetID=cms.InputTag(whichJets)
process.ak4JetFlavourInfos.jets               = newjetID
process.ak4JetFlavourInfos.hadronFlavourHasPriority = cms.bool(True)
process.AK4byRef.jets                         = newjetID
if not "ak4PFJetsCHS" in whichJets:
    process.ak4JetTracksAssociatorAtVertexPF.jets = newjetID
    process.pfImpactParameterTagInfos.jets        = newjetID
    process.softPFMuonsTagInfos.jets              = newjetID
    process.softPFElectronsTagInfos.jets          = newjetID
    process.patJetGenJetMatch.src                 = newjetID


from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from DQMOffline.RecoB.bTagTrackIPAnalysis_cff import *
from DQMOffline.RecoB.bTagGenericAnalysis_cff import *

process.customPfImpactParameterTagInfos = pfImpactParameterTagInfos.clone()
process.customPfImpactParameterTagInfos.useMvaSelection = cms.bool(True)
#process.customPfImpactParameterTagInfos.minimumMvaDiscriminant = cms.double(0.04)
#process.customPfImpactParameterTagInfos.weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/TMVAClassification_BDT.weights.xml.gz')

process.customPfCombinedInclusiveSecondaryVertexV2BJetTags = pfCombinedInclusiveSecondaryVertexV2BJetTags.clone()
process.customPfCombinedInclusiveSecondaryVertexV2BJetTags.tagInfos = cms.VInputTag(cms.InputTag("customPfImpactParameterTagInfos"), cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfos"))
process.customCandidateCombinedSecondaryVertexV2Computer = candidateCombinedSecondaryVertexV2Computer.clone()

#process.customCandidateCombinedSecondaryVertexV2Computer.trackSelection = trackSelection
#process.customCandidateCombinedSecondaryVertexV2Computer.trackPseudoSelection = trackPseudoSelection

process.customPfCombinedInclusiveSecondaryVertexV2BJetTags.jetTagComputer = cms.string("customCandidateCombinedSecondaryVertexV2Computer")
process.customTrackBTagging = cms.Sequence(
    (
      # impact parameters and IP-only algorithms
      process.customPfImpactParameterTagInfos *
      ( 
        # SV tag infos depending on IP tag infos, and SV (+IP) based algos
        inclusiveCandidateVertexing *
        pfInclusiveSecondaryVertexFinderTagInfos *
        process.customPfCombinedInclusiveSecondaryVertexV2BJetTags
      )
    )
)

tagConfig = cms.VPSet(
    cms.PSet(
            bTagTrackIPAnalysisBlock,
            type = cms.string('CandIP'),
            label = cms.InputTag("customPfImpactParameterTagInfos"),
            folder = cms.string("IPTagCustom")
        ),
    cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("customPfCombinedInclusiveSecondaryVertexV2BJetTags"),
            folder = cms.string("CSVv2Custom")
        )
    )

process.btagSequence = cms.Sequence(
    process.ak4JetTracksAssociatorAtVertexPF *
    process.customTrackBTagging
    )
process.jetSequences = cms.Sequence(process.goodOfflinePrimaryVertices * process.btagSequence)

###
print "inputTag : ", process.ak4JetTracksAssociatorAtVertexPF.jets
###

if runOnMC:
    process.flavourSeq = cms.Sequence(
        process.selectedHadronsAndPartons *
        process.ak4JetFlavourInfos
        )
    process.load("Validation.RecoB.bTagAnalysis_cfi")
    process.bTagValidation.jetMCSrc = 'ak4JetFlavourInfos'
    if "Calo" in whichJets:
        process.bTagValidation.caloJetMCSrc = 'AK4byValAlgo'
        process.bTagValidation.useOldFlavourTool = True
        process.flavourSeq = cms.Sequence(
            process.myPartons *
            process.AK4Flavour
            )

    process.bTagValidation.tagConfig += tagConfig
    process.bTagHarvestMC.tagConfig += tagConfig

    process.bTagValidation.applyPtHatWeight = False
    process.bTagValidation.doJetID = True
    process.bTagValidation.doJEC = applyJEC
    process.bTagValidation.JECsourceMC = cms.InputTag(corrLabel+"L1FastL2L3Corrector")
    process.bTagValidation.flavPlots = flavPlots
    process.bTagHarvestMC.flavPlots = flavPlots
    #process.bTagValidation.ptRecJetMin = cms.double(20.)
    process.bTagValidation.genJetsMatched = cms.InputTag("patJetGenJetMatch")
    process.bTagValidation.doPUid = cms.bool(True)
    process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
                                             src = cms.InputTag("ak4GenJets"),
                                             cut = cms.string('pt > 8.'),
                                             filter = cms.bool(False)
                                             )
    process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
    process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
    process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)
else:
    process.load("DQMOffline.RecoB.bTagAnalysisData_cfi")
    process.bTagAnalysis.doJEC = applyJEC
    process.bTagAnalysis.JECsourceData = cms.InputTag(corrLabel+"L1FastL2L3ResidualCorrector")
    process.JECseq *= (getattr(process,corrLabel+"ResidualCorrector") * getattr(process,corrLabel+"L1FastL2L3ResidualCorrector"))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring()
)

from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if useTrigger: 
    process.bTagHLT  = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_PFJet40_v*"])
    process.bTagHLT.HLTPaths = [triggerPath]

if runOnMC:
    process.dqmSeq = cms.Sequence(process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.bTagValidation * process.bTagHarvestMC * process.dqmSaver)
else:
    process.dqmSeq = cms.Sequence(process.bTagAnalysis * process.bTagHarvest * process.dqmSaver)

if useTrigger:
    process.plots = cms.Path(process.bTagHLT * process.JECseq * process.jetSequences * process.dqmSeq)
else:
    process.plots = cms.Path(process.JECseq * process.jetSequences * process.dqmSeq)
    
process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd =cms.untracked.bool(True) 
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)
process.PoolSource.fileNames = [
 "file:/afs/cern.ch/work/s/selvaggi/public/00C31A90-2237-E611-9C70-002590D0AFD0.root"
]

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.GlobalTag = tag

open('pydump.py','w').write(process.dumpPython())
