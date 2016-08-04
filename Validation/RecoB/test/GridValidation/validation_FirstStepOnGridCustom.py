# The following comments couldn't be translated into the new config version:
#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms
process = cms.Process("validation")

# load the full reconstraction configuration, to make sure we're getting all needed dependencies
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.Reconstruction_cff")

whichJets  = "ak4PFJetsCHS"
applyJEC = True
corrLabel = 'ak4PFCHS'
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
process.JECseq = cms.Sequence(getattr(process,corrLabel+"L1FastL2L3CorrectorChain"))

newjetID=cms.InputTag(whichJets)
process.ak4JetFlavourInfos.jets = newjetID
process.ak4JetFlavourInfos.hadronFlavourHasPriority = cms.bool(True)
if not "ak4PFJetsCHS" in whichJets:
    process.ak4JetTracksAssociatorAtVertexPF.jets = newjetID
    process.pfImpactParameterTagInfos.jets        = newjetID
    process.softPFMuonsTagInfos.jets              = newjetID
    process.softPFElectronsTagInfos.jets          = newjetID
    process.patJetGenJetMatch.src                 = newjetID

'''process.btagging = cms.Sequence(process.legacyBTagging + process.pfBTagging)
process.btagSequence = cms.Sequence(
    process.ak4JetTracksAssociatorAtVertexPF *
    process.btagging
    )
'''

from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from DQMOffline.RecoB.bTagTrackIPAnalysis_cff import *
from DQMOffline.RecoB.bTagGenericAnalysis_cff import *


process.customCandidateCombinedSecondaryVertexV2Computer = candidateCombinedSecondaryVertexV2Computer.clone()

trackSelection = cms.PSet(
        a_dR = cms.double(-0.001053),
        a_pT = cms.double(0.005263),
        b_dR = cms.double(0.6263),
        b_pT = cms.double(0.3684),
        jetDeltaRMax = cms.double(0.3),
        maxDecayLen = cms.double(5),
#        maxDecayLen = cms.double(9999999999999999),
        maxDistToAxis = cms.double(0.07),
#        maxDistToAxis = cms.double(999999999999999),
        max_pT = cms.double(500),
        max_pT_dRcut = cms.double(0.1),
        max_pT_trackPTcut = cms.double(3),
        min_pT = cms.double(120),
        min_pT_dRcut = cms.double(0.5),
        normChi2Max = cms.double(99999.9),
        pixelHitsMin = cms.uint32(0),
        ptMin = cms.double(0.0),
        qualityClass = cms.string('any'),
        sip2dSigMax = cms.double(99999.9),
        sip2dSigMin = cms.double(-99999.9),
        sip2dValMax = cms.double(99999.9),
        sip2dValMin = cms.double(-99999.9),
        sip3dSigMax = cms.double(99999.9),
        sip3dSigMin = cms.double(-99999.9),
        sip3dValMax = cms.double(99999.9),
        sip3dValMin = cms.double(-99999.9),
        totalHitsMin = cms.uint32(0),
        useVariableJTA = cms.bool(False)
    )

process.customCandidateCombinedSecondaryVertexV2Computer.trackSelection = trackSelection
#process.customCandidateCombinedSecondaryVertexV2Computer.trackPseudoSelection = trackPseudoSelection

process.pfImpactParameterTagInfosSequence = cms.Sequence(pfImpactParameterTagInfos)
process.pfCombinedSecondaryVertexSequence = cms.Sequence(pfCombinedInclusiveSecondaryVertexV2BJetTags)

bTagTrackIPAnalysisBlockCustom = cms.PSet(
    parameters = cms.PSet(
        QualityPlots = cms.bool(False),
        endEffPur = cms.double(1.005),
        nBinEffPur = cms.int32(200),
        startEffPur = cms.double(0.005),
        LowerIPSBound = cms.double(-35.0),
        UpperIPSBound = cms.double(35.0),
        LowerIPBound = cms.double(-0.1),
        UpperIPBound = cms.double(0.1),
        LowerIPEBound = cms.double(0.0),
        UpperIPEBound = cms.double(0.04),
        NBinsIPS = cms.int32(100),
        NBinsIP = cms.int32(100),
        NBinsIPE = cms.int32(100),
        MinDecayLength = cms.double(-9999.0),
        MaxDecayLength = cms.double(999.0),
        MinJetDistance = cms.double(0.0),
        MaxJetDistance = cms.double(9999.0),
    )
)

## 1rst loop here


xmin = 0.0;
xmax = 0.2;
nstep = 20;

tagConfigMain = cms.VPSet(
    cms.PSet(
            bTagTrackIPAnalysisBlock,
            type = cms.string('CandIP'),
            label = cms.InputTag("pfImpactParameterTagInfos"),
            folder = cms.string("IPTag")
        ),
    cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("pfCombinedInclusiveSecondaryVertexV2BJetTags"),
            folder = cms.string("CSVv2")
        )
  )


for i in range(0,nstep+1):

  postfix = str(i)
  cutval = xmin + (xmax - xmin)/nstep*i

  setattr(process, "customPfImpactParameterTagInfos"+postfix, process.pfImpactParameterTagInfos.clone())
  getattr(process, "customPfImpactParameterTagInfos"+postfix).useMvaSelection = cms.bool(True)
  #getattr(process, "customPfImpactParameterTagInfos"+postfix).useMvaSelection = cms.bool(False)
  getattr(process, "customPfImpactParameterTagInfos"+postfix).minimumMvaDiscriminant = cms.double(cutval)
  getattr(process, "customPfImpactParameterTagInfos"+postfix).weightFile = cms.FileInPath('RecoBTag/SecondaryVertex/data/BDT_nosel_8Var.weights.xml.gz')

  setattr(process, "customPfCombinedInclusiveSecondaryVertexV2BJetTags"+postfix, process.pfCombinedInclusiveSecondaryVertexV2BJetTags.clone())
  getattr(process, "customPfCombinedInclusiveSecondaryVertexV2BJetTags"+postfix).tagInfos = cms.VInputTag(cms.InputTag("customPfImpactParameterTagInfos"+postfix), 
                                                                                                        cms.InputTag("pfInclusiveSecondaryVertexFinderTagInfos"))
  getattr(process, "customPfCombinedInclusiveSecondaryVertexV2BJetTags"+postfix).jetTagComputer = cms.string("customCandidateCombinedSecondaryVertexV2Computer")

  process.pfImpactParameterTagInfosSequence *= getattr(process,"customPfImpactParameterTagInfos"+postfix)
  process.pfCombinedSecondaryVertexSequence *= getattr(process,"customPfCombinedInclusiveSecondaryVertexV2BJetTags"+postfix)

  tagConfig = cms.VPSet(
    cms.PSet(
            bTagTrackIPAnalysisBlockCustom,
            type = cms.string('CandIP'),
            label = cms.InputTag("customPfImpactParameterTagInfos"+postfix),
            folder = cms.string("IPTagCustom"+postfix)
        ),
    cms.PSet(
            bTagGenericAnalysisBlock,
            label = cms.InputTag("customPfCombinedInclusiveSecondaryVertexV2BJetTags"+postfix),
            folder = cms.string("CSVv2Custom"+postfix)
        )
  )

  process.load("Validation.RecoB.bTagAnalysis_cfi")
  tagConfigMain += tagConfig

############# end of loop ################


process.bTagValidation.tagConfig = tagConfigMain
process.bTagHarvestMC.tagConfig = tagConfigMain

process.customTrackBTagging = cms.Sequence(
    ( 
      # impact parameters and IP-only algorithms
      process.pfImpactParameterTagInfosSequence *
      ( 
        # SV tag infos depending on IP tag infos, and SV (+IP) based algos
        inclusiveCandidateVertexing *
        pfInclusiveSecondaryVertexFinderTagInfos *
        process.pfCombinedSecondaryVertexSequence
      )
    )
)

process.btagSequence = cms.Sequence(
    process.ak4JetTracksAssociatorAtVertexPF*
    process.customTrackBTagging 
    )

process.jetSequences = cms.Sequence(process.goodOfflinePrimaryVertices * process.btagSequence)

###
print "inputTag : ", process.ak4JetTracksAssociatorAtVertexPF.jets
###

process.load("Validation.RecoB.bTagAnalysis_firststep_cfi")
if runOnMC:
    process.flavourSeq = cms.Sequence(
        process.selectedHadronsAndPartons *
        process.ak4JetFlavourInfos
        )
    process.bTagValidationFirstStep.jetMCSrc = 'ak4JetFlavourInfos'
    process.bTagValidationFirstStep.applyPtHatWeight = False
    process.bTagValidationFirstStep.doJetID = True
    process.bTagValidationFirstStep.doJEC = applyJEC
    process.bTagValidation.JECsourceMC = cms.InputTag(corrLabel+"L1FastL2L3Corrector")
    process.bTagValidationFirstStep.flavPlots = flavPlots
    #process.bTagValidationFirstStep.ptRecJetMin = cms.double(20.)
    process.bTagValidationFirstStep.genJetsMatched = cms.InputTag("patJetGenJetMatch")
    process.bTagValidationFirstStep.doPUid = cms.bool(True)
    process.ak4GenJetsForPUid = cms.EDFilter("GenJetSelector",
                                             src = cms.InputTag("ak4GenJets"),
                                             cut = cms.string('pt > 8.'),
                                             filter = cms.bool(False)
                                             )
    process.load("PhysicsTools.PatAlgos.mcMatchLayer0.jetMatch_cfi")
    process.patJetGenJetMatch.matched = cms.InputTag("ak4GenJetsForPUid")
    process.patJetGenJetMatch.maxDeltaR = cms.double(0.25)
    process.patJetGenJetMatch.resolveAmbiguities = cms.bool(True)
else:
    process.bTagValidationFirstStepData.doJEC = applyJEC
    process.bTagAnalysis.JECsourceData = cms.InputTag(corrLabel+"L1FastL2L3ResidualCorrector")
    process.JECseq *= (getattr(process,corrLabel+"ResidualCorrector") * getattr(process,corrLabel+"L1FastL2L3ResidualCorrector"))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("file:/afs/cern.ch/work/s/selvaggi/public/00C31A90-2237-E611-9C70-002590D0AFD0.root")
)

process.EDM = cms.OutputModule("DQMRootOutputModule",
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      "keep *_MEtoEDMConverter_*_*"),
                               fileName = cms.untracked.string('MEtoEDMConverter.root')
                               )
from HLTrigger.HLTfilters.hltHighLevel_cfi import *
if useTrigger: 
    process.bTagHLT  = hltHighLevel.clone(TriggerResultsTag = "TriggerResults::HLT", HLTPaths = ["HLT_PFJet40_v*"])
    process.bTagHLT.HLTPaths = [triggerPath]

if runOnMC:
    process.dqmSeq = cms.Sequence(process.ak4GenJetsForPUid * process.patJetGenJetMatch * process.flavourSeq * process.bTagValidationFirstStep)
else:
    process.dqmSeq = cms.Sequence(process.bTagValidationFirstStepData)

if useTrigger:
    process.plots = cms.Path(process.bTagHLT * process.JECseq * process.jetSequences * process.dqmSeq)
else:
    process.plots = cms.Path(process.JECseq * process.jetSequences * process.dqmSeq)
    
process.outpath = cms.EndPath(process.EDM)

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

