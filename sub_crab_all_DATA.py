
from CRABAPI.RawCommand import crabCommand
from CRABClient.ClientExceptions import ClientException
#from CRABClient.UserUtilities #import config, getUsernameFromSiteDB
from httplib import HTTPException

from WMCore.Configuration import Configuration


# We want to put all the CRAB project directories from the tasks we submit here into one common directory.
# That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
config = Configuration()
config.section_("General")
config.General.workArea = '/afs/cern.ch/work/m/mthiel/private/ANALYZER/Exclusive_VV_PPS/2017_VV/CMSSW_10_6_2/src/MakeNTuple/MakeNTuple/crab_projects/'
config.General.transferOutputs = True
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'ConfFile_data_cfg.py'
config.JobType.inputFiles = ['MyDataPileupHistogram_2017.root','PileupMC_2017.root']
config.JobType.outputFiles = ['out.root']
config.JobType.allowUndistributedCMSSW = True

   
config.section_("Data")
config.Data.useParent = True
config.Data.ignoreLocality = True
config.Data.splitting = 'LumiBased' #LumiBased'
config.Data.unitsPerJob = 50
config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'  #'Automatic'
#config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/mthiel' #%s/' % (getUsernameFromSiteDB())
config.Data.publication = False
config.Data.outputDatasetTag =  'data_2017'
config.Data.lumiMask = '/afs/cern.ch/work/m/mthiel/private/ANALYZER/Exclusive_VV_PPS/2017_VV/CMSSW_10_6_2/src/MakeNTuple/MakeNTuple/Cert_13TeV_Collisions17_JSON_PPSruns.txt'
#    config.Data.runRange = '271036-284044'

config.section_("Site")
config.Site.whitelist = ['T2_CH_*','T2_DE_*','T2_IT_*','T2_US_*']
config.Site.storageSite = "T2_CH_CERN"



def submit(config):
    try:
        crabCommand('submit', config = config)
    except HTTPException as hte:
        print "Failed submitting task: %s" % (hte.headers)
    except ClientException as cle:
        print "Failed submitting task: %s" % (cle)

#############################################################################################
## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
#############################################################################################

if __name__ == '__main__':


#example
#config.General.requestName = 'FSQJet1Bprompt'
#config.Data.inputDataset = '/FSQJet1/Run2018B-PromptReco-v2/MINIAOD'
#submit(config)


    config.General.requestName = 'SingleMuon_B'
    config.Data.inputDataset = '/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleMuon_C'
    config.Data.inputDataset = '/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleMuon_D'
    config.Data.inputDataset = '/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleMuon_E'
    config.Data.inputDataset = '/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleMuon_F'
    config.Data.inputDataset = '/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleElectron_B'
    config.Data.inputDataset = '/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleElectron_C'
    config.Data.inputDataset = '/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleElectron_D'
    config.Data.inputDataset = '/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleElectron_E'
    config.Data.inputDataset = '/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'SingleElectron_F'
    config.Data.inputDataset = '/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleMuon_B'
    config.Data.inputDataset = '/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleMuon_C'
    config.Data.inputDataset = '/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleMuon_D'
    config.Data.inputDataset = '/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleMuon_E'
    config.Data.inputDataset = '/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleMuon_F'
    config.Data.inputDataset = '/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_B'
    config.Data.inputDataset = '/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_C'
    config.Data.inputDataset = '/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_D'
    config.Data.inputDataset = '/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_E'
    config.Data.inputDataset = '/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD'
    submit(config)

    config.General.requestName = 'DoubleEG_F'
    config.Data.inputDataset = '/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD'
    submit(config)

