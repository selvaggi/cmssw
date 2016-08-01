# The following comments couldn't be translated into the new config version:

#! /bin/env cmsRun

import FWCore.ParameterSet.Config as cms




runOnMC = True

process = cms.Process("harvest")
process.load("DQMServices.Components.DQMEnvironment_cfi")

#keep the logging output to a nice level
process.load("FWCore.MessageLogger.MessageLogger_cfi")

process.load("DQMServices.Core.DQM_cfg")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.source = cms.Source("DQMRootSource",
    fileNames = cms.untracked.vstring()
)

process.load("Validation.RecoB.bTagAnalysis_harvesting_cfi")


from RecoBTag.ImpactParameter.impactParameter_cff import *
from RecoBTag.SecondaryVertex.secondaryVertex_cff import *
from RecoVertex.AdaptiveVertexFinder.inclusiveVertexing_cff import *
from DQMOffline.RecoB.bTagTrackIPAnalysis_cff import *
from DQMOffline.RecoB.bTagGenericAnalysis_cff import *

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
nstep = 25;

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

process.bTagValidationHarvest.tagConfig = tagConfigMain

if runOnMC:
    process.dqmSeq = cms.Sequence(process.bTagValidationHarvest * process.dqmSaver)
#    process.dqmSeq = cms.Sequence(process.bTagValidationHarvest)

else:
    process.dqmSeq = cms.Sequence(process.bTagValidationHarvestData * process.dqmSaver)

process.plots = cms.Path(process.dqmSeq)

process.dqmEnv.subSystemFolder = 'BTAG'
process.dqmSaver.producer = 'DQM'
process.dqmSaver.workflow = '/POG/BTAG/BJET'
process.dqmSaver.convention = 'Offline'
process.dqmSaver.saveByRun = cms.untracked.int32(-1)
process.dqmSaver.saveAtJobEnd = cms.untracked.bool(True) 
process.dqmSaver.forceRunNumber = cms.untracked.int32(1)

process.DQMRootSource.fileNames = [
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_10.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_18.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_20.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_22.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_28.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_29.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_31.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_32.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_36.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_39.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_3.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_40.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_41.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_46.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_47.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_56.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_57.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_58.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_60.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_61.root',
'file:/storage/data/cms/store/user/selvaggi/bTag/trackOpti/TT_TuneCUETP8M1_13TeV-powheg-pythia8/crab_Track_optimization_fullTrackSel/160801_092649/0000/MEtoEDMConverter_9.root'
]

open('pydump.py','w').write(process.dumpPython())
