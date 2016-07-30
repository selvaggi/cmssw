from WMCore.Configuration import Configuration

config = Configuration()

config.section_("General")
config.General.requestName   = 'Track_optimization_test_1'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName  = 'ANALYSIS'
# Name of the CMSSW configuration file
config.JobType.psetName    = 'validation_FirstStepOnGridCustom.py'

config.section_("Data")
config.Data.inputDataset = '/TT_TuneCUETP8M1_13TeV-powheg-pythia8/RunIISpring16DR80-PUSpring16_RECODEBUG_80X_mcRun2_asymptotic_2016_v3_ext3-v1/GEN-SIM-RECODEBUG'

config.Data.splitting = 'FileBased'
#config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 3
config.Data.totalUnits = -1
config.Data.outLFNDirBase = "/store/user/selvaggi/bTag/trackOpti"
config.Data.publication = False
# This string is used to construct the output dataset name
#config.Data.publishDataName = ''


config.section_("Site")
# Where the output files will be transmitted to
config.Site.storageSite='T2_BE_UCL'
#config.Site.whitelist=['T2_RU_JINR','T2_IT_Pisa', 'T2_US_Caltech', 'T2_DE_DESY', 'T1_US_FNAL']
