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

process.load("Validation.RecoB.bTagAnalysis_cfi")

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

#    process.bTagValidation.tagConfig += tagConfig
#    process.bTagHarvestMC.tagConfig += tagConfig

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
    input = cms.untracked.int32(-1)
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
"file:/afs/cern.ch/work/s/selvaggi/public/00C31A90-2237-E611-9C70-002590D0AFD0.root",

'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/00C31A90-2237-E611-9C70-002590D0AFD0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/025349EA-C037-E611-B304-20CF3027A594.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/04BFF711-3237-E611-85C0-20CF3027A5BB.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/061C86C8-4237-E611-97B7-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0636224F-9C37-E611-9C5B-20CF3027A5FD.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/06AF4BD1-4037-E611-9EDA-002590D60026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/08192134-EC37-E611-BD60-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/08E0FE77-1837-E611-97F0-00221982D02D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0C6F6817-4F37-E611-8C1D-00221982B6B6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0E0845DC-9E37-E611-BA79-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0E53BEAF-2537-E611-8309-1418774126FB.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0E6ADA36-3937-E611-A562-F04DA2770C8E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0E9D935D-5337-E611-9F32-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0EA2294B-2F37-E611-A065-E0DB55FC1135.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0ED7CE2D-2737-E611-8803-002590D5FFDA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/0EEA5626-3437-E611-8381-BCAEC509720C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1067A15E-2D37-E611-9692-001C23C0F175.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/10DD9D23-1E37-E611-A3FB-14187741120B.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/12D16CDE-4337-E611-A40D-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1485E240-1037-E611-9D0C-B083FED3EE25.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/14D32FBF-9A38-E611-9C4F-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1607EB2B-9D37-E611-9C55-20CF3019DF00.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/188144AD-1937-E611-B4CC-D4AE526DF801.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1C3D3299-0C37-E611-A76C-001EC94B4F45.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1CDB2E21-FE36-E611-A5B7-002590D60026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1E07C227-A437-E611-9136-002590D60098.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1E694F4F-A038-E611-8BFD-F04DA2770C8E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/1EF5EF29-1337-E611-8EB8-90B11C0BDC70.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/208DC7AE-2B37-E611-8CDF-0019B9CABBF1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/20EB8C19-C438-E611-9672-20CF3019DF00.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/22BEBA7F-3E37-E611-9338-001EC9DE20A3.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/22C8CA6B-1537-E611-8238-141877410E71.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2485D091-2537-E611-9EFF-0026B937D37D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/26126F05-4337-E611-9F88-44A8423C4026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2656AB21-3D37-E611-8858-20CF3027A5D4.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2688B698-9737-E611-AB51-E0DB55FC1055.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2AB0494D-2E37-E611-9437-00221982D02D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2CE1AF19-A537-E611-BF77-002590D60150.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2CEEA030-2C37-E611-8492-90B11C0BD312.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2E198633-8A37-E611-A532-3417EBE885A6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2E2D1138-6937-E611-B506-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/2ED87C40-9438-E611-9BD8-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/30352CAA-9F37-E611-A2C8-20CF3019DF0B.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/306FEEF0-4D37-E611-8953-002590D60004.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/30C499EC-0637-E611-BEF7-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/32E0F945-3437-E611-8A3E-44A84224BE51.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/34AB9453-1937-E611-9C67-549F3525BF58.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/36373CC6-E337-E611-ABC8-485B3989725F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/38B1090E-4937-E611-993B-20CF3027A59A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3A2445F3-1737-E611-9E4B-485B398168B7.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3A521CA7-0937-E611-9F85-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3A5D2190-2E37-E611-B40D-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3A7CAB59-A638-E611-885B-002590D5FFDA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3C484F34-FC36-E611-A83B-002590D600A4.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3C4EA490-1737-E611-82AC-90B11C0BCF43.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3E204A9A-1637-E611-8A4C-1418774121A1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3EE736ED-3737-E611-ABF7-00221982D6CC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/3EFF3213-3537-E611-8208-0019B9CABBF1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4025B92A-4237-E611-A346-44A8423C4026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/402C43BD-A537-E611-88EA-002590D60072.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/408736AC-0E37-E611-B913-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4088ACF7-2337-E611-8AE9-14187741121F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/420A8D03-1D37-E611-9069-001C23C0C7AB.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4256562A-2737-E611-B533-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/42B8E527-0B37-E611-B5B8-44A8423C2A62.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/441E1424-5A37-E611-8346-BCAEC5097210.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/442D5A7D-1A37-E611-BF19-001EC94BE9F4.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/448D2064-9E37-E611-A7BD-00221982B6B6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/44A913D3-2F37-E611-93AD-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/46473F92-1837-E611-ADD2-782BCB206470.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/46BD0F8F-0737-E611-8F56-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4AE25984-0D37-E611-B3DC-1418774108DF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4AF77B0C-4937-E611-B669-002590D60068.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4C65B801-3937-E611-A5A7-44A8423C4026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4CB9B49A-1337-E611-9643-782BCB1CFD1C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4CCD52EC-A037-E611-B54E-20CF305B04D6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/4CE0C2E2-2638-E611-9916-44A8423D7989.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/501F70D2-2937-E611-9505-20CF3027A5FA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/502247BA-9E37-E611-B1DF-F04DA27710A2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/54BCD0DF-3D37-E611-B064-002590D60004.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/54C59813-E837-E611-AA4A-44A8423D7E31.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/56B3D111-CB37-E611-AEBD-002590D0AFD0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/5A7D6971-2437-E611-815E-549F3525AE58.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/5CF8AA6E-2737-E611-B9E8-0019B9CABBF1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/5E5ACAC4-5737-E611-8745-3417EBE5354A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/5E9FFF41-1537-E611-BEA2-BCAEC509720C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/60282836-0E37-E611-9F6D-782BCB538FDF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/603F7CC7-1137-E611-A6AF-002590D6004A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/62CB633B-A537-E611-BAA7-00221982C61D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/640995AF-1537-E611-A1F9-A4BADB1C5E28.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/646A6E9B-1237-E611-8627-20CF3019DF0B.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/66880AE7-5537-E611-91B0-44A8423D7989.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/686C3785-3637-E611-ACC4-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6A36DA86-5737-E611-83C3-485B39897231.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6A4BE09F-3237-E611-A038-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6ADC6FDF-A937-E611-846E-002590D5FFFA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6AEBDDCE-FC36-E611-84B4-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6C266063-0337-E611-93C0-549F3525C0BC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6C8E6273-2E37-E611-8ED4-0019B9CABBF1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6CA0DB22-8537-E611-B3F7-44A8423D7E31.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6CA5B9CB-3A37-E611-AD9A-549F3525A184.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6CB9B1CE-8938-E611-B10F-20CF305B04D6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/6E6849F0-2037-E611-8CC6-B083FED3F4E3.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/70B440B4-1C37-E611-A6FE-20CF3027A5F6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7224A9D8-3137-E611-8D1B-0022197BDD29.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/72259083-3537-E611-A3E2-20CF3027A5FD.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7413AE92-4537-E611-A4EC-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/74413C4A-6337-E611-BD31-90B11C0BD312.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/745E72A9-2437-E611-A53B-782BCB47D0B6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7484BEC9-AF37-E611-995C-44A8423C2A62.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/781145B8-0937-E611-9442-141877410ACD.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7837361D-4737-E611-B0FA-20CF3019DF0B.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7881A918-4337-E611-9642-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/788A5D1C-0F37-E611-8190-842B2B18118F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/78DD8A64-1937-E611-A5C4-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7A3FD606-1837-E611-A85E-842B2B181725.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7A82C181-9B37-E611-AFDB-20CF3027A61E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7A8CC227-1337-E611-AC0D-20CF3027A62D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7C4100A2-2637-E611-89AA-1418774126FB.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7CF9CC61-1E37-E611-A391-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7EE23E1B-A337-E611-B392-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/7EF7A876-A037-E611-8537-3417EBE5354A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/80086627-4337-E611-A0B5-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/804C33CD-1337-E611-A084-782BCB539A14.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/80986B72-B737-E611-9B4D-20CF3027A566.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/80D8BEC1-1137-E611-A1A1-D4AE527EDFE6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/82AABB88-9C37-E611-AD34-3417EBE886B2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8448EA13-4337-E611-8AF7-BCAEC509720E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/844C8A83-3437-E611-A35B-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/84C91E0C-1837-E611-8C69-141877411367.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/84E04081-2337-E611-95BB-549F3525CD78.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/86007447-4437-E611-95A9-002590D0AFD0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8688A87A-1837-E611-9580-F04DA27710A2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/869992C0-3637-E611-B465-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/86B53FCC-0C37-E611-9CF1-BCAEC509720C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/88471360-1B37-E611-ACA6-E0DB55FC1135.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8A135AC9-0637-E611-B121-002590D60004.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8A975B2D-1837-E611-B0C0-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8A9E1F82-3E37-E611-8C36-20CF3019DF00.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8C0DB798-A137-E611-965E-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8C31CA3F-0F37-E611-8549-D4AE527EE0EB.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8CC50873-1337-E611-BFEB-782BCB539226.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/8CDE0945-A437-E611-B8EE-E0DB55FC1135.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/90188C47-2337-E611-AC6E-44A8423C2A62.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/90296960-FE36-E611-B5D6-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/90614D8D-4637-E611-A67D-20CF3027A5D8.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/927CACE0-1C37-E611-9F60-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/92820BED-1937-E611-8583-A4BADB1CF89C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/94622ECE-4C37-E611-8143-002590D60098.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/94956961-1E37-E611-85B9-0026B937D37D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/96466B5B-0137-E611-9B36-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/981C63B2-1B37-E611-9C8D-20CF3027A589.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/9878E977-3737-E611-8F95-F04DA27710A2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/98EBFA0B-4C37-E611-88B8-20CF3027A566.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/9ABFB534-1537-E611-823A-549F3525DFE8.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/9AEDC5D8-FB36-E611-9105-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A0DE503F-3137-E611-8559-002590D60098.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A214ADE6-A937-E611-8D03-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A2243B7F-1E37-E611-B4BF-20CF3027A589.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A265C275-3039-E611-85CC-002590D6004A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A41CAB38-0F37-E611-88D9-782BCB47D0B6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A45AB7AD-9A37-E611-A799-20CF3027A5DE.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A609F78D-A137-E611-BDDF-44A8423D7989.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A68E0539-2237-E611-A219-549F3525DD6C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/A8455C17-F836-E611-902A-90B11C0BDC70.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/ACC7E759-BC37-E611-B344-E0DB55FC1135.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/AE64ED9B-2137-E611-8C69-3417EBE535DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/AE81C99E-3E37-E611-A146-3417EBE539DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/AEB09688-1C37-E611-B7B1-B083FED429D6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B0D3959E-2D37-E611-9C82-549F3525C0BC.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B273FA2A-3B37-E611-8F8B-BCAEC509720E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B2DDBEB4-4137-E611-896A-44A8423D7989.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B2FC82FC-1137-E611-915E-842B2B17F557.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B43D8335-E537-E611-A40A-002590D6004A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B850E392-1137-E611-AD47-B083FED177B1.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B86421E4-FA36-E611-82DA-002590D5FFDA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/B8F2C572-1637-E611-B469-E0DB55FC1135.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BA77E677-5037-E611-ABC5-F04DA27710A2.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BA95682B-1B37-E611-AE19-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BC9A0B85-3037-E611-AE5C-0019B9CADD7A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BCC5A54F-1937-E611-8F26-C81F66DCFE01.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BCFB36BD-9D37-E611-8255-20CF3027A594.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BE9548E3-3837-E611-8535-E0DB55FC100D.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BECF73EA-1737-E611-AD39-782BCB538FDF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/BEF410F6-9E37-E611-A994-20CF3027A5FA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C0362529-4137-E611-9D37-E0DB55FC11A5.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C064A70B-3B37-E611-A59D-E0DB55FC1055.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C0898B15-3737-E611-B0EB-20CF3027A610.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C29500FD-3037-E611-AE9E-20CF3027A566.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C61B4AF6-9A37-E611-B37B-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C660BCF1-1537-E611-A04A-A4BADB1CF89C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C69C8726-9D37-E611-9077-001EC9EAFEBE.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C6BC554A-4037-E611-A1A8-E0DB55FC11A5.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/C8507787-1E37-E611-A60C-549F3525B9A0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/CC1BCC2E-0837-E611-A9E6-C81F66B78FF5.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/CC8F425C-1D37-E611-8DFB-44A8423D7989.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D00F4A94-3F38-E611-9254-3417EBE5361A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D051C035-A037-E611-A0FE-E0DB55FC11A5.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D0B751A4-DA37-E611-8579-20CF3027A5F6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D2668DAF-4237-E611-8103-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D26BE6F1-1A37-E611-A92D-3417EBE5354A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D2C22BB5-9C37-E611-BBD4-20CF3027A59A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D2F54FA3-A237-E611-8AB1-44A8423CE96E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D4119666-2137-E611-BD03-B083FED429D6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D438C3EC-1437-E611-BF70-782BCB538FDF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D46E9F05-F936-E611-B601-485B398168B7.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D676FDC3-2C37-E611-B18C-782BCB539A41.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D68E32D8-F836-E611-B9D3-14187741208F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D6990375-1637-E611-BF78-782BCB47D0B6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D83EEE1A-9D37-E611-9C14-3417EBE5361A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D8727C8D-1D37-E611-B842-001EC94BF6CA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D8C49556-3437-E611-81F1-002590D60026.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/D8DC1049-9A37-E611-93A4-20CF3019DF18.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DA17B5AC-6437-E611-8B2E-3417EBE53662.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DA2CFAFD-2737-E611-BE97-549F3525B9A0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DA7F2F8B-3F37-E611-961C-E0DB55FC11A5.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DE22B94E-A037-E611-9B54-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DE4BC196-2237-E611-8D2F-782BCB539A41.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DE8A2107-1D37-E611-B2C6-141877410E71.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/DEF4B534-B437-E611-9CC4-B083FED00118.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E028A1D0-5637-E611-AA9F-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E07F659B-A237-E611-9C92-485B39897231.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E205572F-9F37-E611-B0B9-3417EBE539DA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E2243794-1837-E611-A9E5-1418774108DF.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E2D24AF7-1137-E611-8630-842B2B17EA37.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E4E2403B-5C37-E611-8CDB-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E6346CDB-3C37-E611-A438-E0DB55FC1055.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E6BCAEFD-9037-E611-8507-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E8580E2C-8D37-E611-9009-20CF3019DF00.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E871BF30-9B37-E611-91CB-E0DB55FC1139.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/E8B24F6A-0437-E611-89C9-3417EBE885A6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/EA010392-1437-E611-A502-D4AE527EDBD4.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/EE1DF1D8-0F37-E611-B62B-44A8423C2A62.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F04ED9BB-ED37-E611-996A-20CF3027A5F6.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F0F4F1D6-BE37-E611-94F7-3417EBE8862E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F2250C81-6238-E611-B681-3417EBE5361A.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F267E7E4-1437-E611-BDD0-14187733AD81.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F45828A7-E137-E611-B914-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F64DF6D5-9D37-E611-AD73-20CF3027A61E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F65FD930-A337-E611-9672-44A8423DE2C0.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F6942F86-A437-E611-B563-002590D60098.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F6DB1BD3-3837-E611-B503-00221982C62E.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F809D730-1D37-E611-8F99-782BCB536C90.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/F8BCE210-3D37-E611-8504-44A8423CF41F.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/FA4EE7A7-1737-E611-8645-A4BADB1CF89C.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/FC234259-4337-E611-B43F-20CF3027A5FD.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/FCA2807D-4337-E611-80DE-20CF3027A5FA.root',
'root://cms-xrd-global.cern.ch//store/mc/RunIISpring16DR80/TT_TuneCUETP8M1_13TeV-powheg-pythia8/GEN-SIM-RECODEBUG/PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/20000/FE40CF67-FC36-E611-94F3-782BCB206470.root'
]

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.GlobalTag = tag

open('pydump.py','w').write(process.dumpPython())
