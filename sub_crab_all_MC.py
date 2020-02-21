if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from CRABClient.UserUtilities import config, getUsernameFromSiteDB
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config = config()
    config.section_("General")
    config.General.workArea = '/afs/cern.ch/work/m/mthiel/private/ANALYZER/Exclusive_VV_PPS/2017_VV/CMSSW_10_6_2/src/MakeNTuple/MakeNTuple/crab_projects_2017_MC/'
    config.General.transferOutputs = True
    config.General.transferLogs = False

    config.section_("JobType")
    config.JobType.pluginName = 'Analysis'
    config.JobType.psetName = 'ConfFile_MC_BKG_cfg.py'
    config.JobType.inputFiles = ['MyDataPileupHistogram_2017.root','PileupMC_2017.root']
    config.JobType.outputFiles = ['out.root']
    config.JobType.allowUndistributedCMSSW = True
   
    config.section_("Data")
#    config.Data.splitting = 'LumiBased' #LumiBased'
#    config.Data.unitsPerJob = 50
    config.Data.inputDBS = 'global'
    config.Data.splitting = 'Automatic' #'FileBased'
    config.Data.outLFNDirBase = '/store/user/mthiel' #%s/' % (getUsernameFromSiteDB())
    config.Data.publication = False
    config.Data.outputDatasetTag =  'background_2'

    config.section_("Site")
    config.Site.storageSite = "T2_BR_UERJ"

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


#example
#config.General.requestName = 'FSQJet1Bprompt'
#config.Data.inputDataset = '/FSQJet1/Run2018B-PromptReco-v2/MINIAOD'
#submit(config)


config.General.requestName = 'W1JetsToLNu_LHEWpT_50-150'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_50-150_2'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_50-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_400-inf'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_400-inf_2'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_400-inf_3'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_400-inf_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_250-400'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_250-400_2'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_250-400_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_150-250'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_150-250_2'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_150-250_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_100-150'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_100-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_100-150_2'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_100-150_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'W1JetsToLNu_LHEWpT_0-50'
config.Data.inputDataset = '/W1JetsToLNu_LHEWpT_0-50_TuneCP5_13TeV-amcnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'TT_TuneCH3'
config.Data.inputDataset = '/TT_TuneCH3_13TeV-powheg-herwig7/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

config.General.requestName = 'DYJetsToLL_M-50'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)


config.General.requestName = 'DYJetsToLL_M-50_2'
config.Data.inputDataset = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)


config.General.requestName = 'DYJetsToLL_M-10to50'
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v2/MINIAODSIM'
submit(config)


config.General.requestName = 'DYJetsToLL_M-10to50_2'
config.Data.inputDataset = '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)


config.General.requestName = 'QCD_Pt-170to300_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-170to300_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-300to470_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-300to470_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-470to600_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-470to600_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-600to800_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-600to800_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-800to1000_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-800to1000_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-1000toInf_MuEnrichedPt5'
config.Data.inputDataset = '/QCD_Pt-1000toInf_MuEnrichedPt5_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-80to120_EMEnriched'
config.Data.inputDataset = '/QCD_Pt-80to120_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-50to80_EMEnriched'
config.Data.inputDataset = '/QCD_Pt-50to80_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'QCD_Pt-120to170_EMEnriched'
config.Data.inputDataset = '/QCD_Pt-120to170_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)



config.General.requestName = 'QCD_Pt-300toInf_EMEnriched'
config.Data.inputDataset = '/QCD_Pt-300toInf_EMEnriched_TuneCP5_13TeV_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'ST_s-channel_4f'
config.Data.inputDataset = '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

#config.General.requestName = ''
#config.Data.inputDataset = '/ST_s-channel_4f_leptonDecays_TuneCP5_13TeV-amcatnlo-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#submit(config)

config.General.requestName = 'ST_t-channel_antitop_4f'
config.Data.inputDataset = '/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

#config.General.requestName = ''
#config.Data.inputDataset = '/ST_t-channel_antitop_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#submit(config)

config.General.requestName = 'ST_t-channel_top_4f'
config.Data.inputDataset = '/ST_t-channel_top_4f_inclusiveDecays_TuneCP5_13TeV-powhegV2-madspin-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'ST_tW_top_5f'
config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

#config.General.requestName = ''
#config.Data.inputDataset = '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#submit(config)

config.General.requestName = 'ST_tW_antitop_5f'
config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)

#config.General.requestName = ''
#config.Data.inputDataset = '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
#submit(config)

config.General.requestName = 'WWToLNuQQ'
config.Data.inputDataset = '/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'WWToLNuQQ_2'
config.Data.inputDataset = '/WWToLNuQQ_NNPDF31_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'ZZTo2L2Q'
config.Data.inputDataset = '/ZZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'WZTo3LNu'
config.Data.inputDataset = '/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'WZTo2L2Q'
config.Data.inputDataset = '/WZTo2L2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'WZTo1L3Nu'
config.Data.inputDataset = '/WZTo1L3Nu_13TeV_amcatnloFXFX_madspin_pythia8_v2/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
submit(config)

config.General.requestName = 'WZTo1L1Nu2Q'
config.Data.inputDataset = '/WZTo1L1Nu2Q_13TeV_amcatnloFXFX_madspin_pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM'
submit(config)



