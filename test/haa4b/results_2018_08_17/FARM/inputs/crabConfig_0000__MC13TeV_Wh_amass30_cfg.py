from CRABClient.UserUtilities import config
config = config()
import os

config.General.requestName = "MC13TeV_Wh_amass30_0"
config.General.workArea = "/afs/cern.ch/user/y/yuanc/Analysis/CMSSW_8_0_26_patch1/src/UserCode/bsmhiggs_fwk/test/haa4b/results_2018_08_17/FARM/inputs/"
config.General.transferLogs = False

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = "/afs/cern.ch/user/y/yuanc/Analysis/CMSSW_8_0_26_patch1/src/UserCode/bsmhiggs_fwk/test/haa4b/results_2018_08_17/MC13TeV_Wh_amass30_0_cfg.py"
config.JobType.sendPythonFolder = True
config.JobType.inputFiles = ["/afs/cern.ch/user/y/yuanc/Analysis/CMSSW_8_0_26_patch1/src/UserCode/bsmhiggs_fwk/test/haa4b/results_2018_08_17/FARM/inputs//x509_proxy"]
config.JobType.outputFiles = ["analysis.root"]

config.Data.inputDataset = '/MinBias/sghiasis-crab_Wh_production_m30_Dec4_CMSSW8021-28028af67189b3de7224b79195bd0e1d/USER'
config.Data.lumiMask = ""
config.Data.inputDBS = "phys03"
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
config.Data.totalUnits = -1
config.Data.publication = False
config.Data.ignoreLocality = True
config.Data.outLFNDirBase = '/store/user/yuanc/results_2018_08_17'

config.Site.storageSite = 'T2_CH_CERN'
config.Site.whitelist = ['T2_CH_CERN','T2_US_UCSD','T3_US_FNALLPC']