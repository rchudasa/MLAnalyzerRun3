from CRABClient.UserUtilities import config
config = config()

Process = 'TauData'

inputDataset_ ={
        'TauData': '/Tau/lpcml-Tau_Run2023C_RAW-AOD_EventAwareLumiv2-8ff03ad286ef16aa6b3d480f18073977/USER',
        'Parking': '/ParkingHH/lpcml-ParkingHH_Run2023C_RAW-AOD_EventAwareLumiv2-8ff03ad286ef16aa6b3d480f18073977/USER',
        }.get(Process, None)


#config.section_('General')
config.General.requestName = '%s_MLAnalyzer_ntuples_v2'%Process
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_data_cfg.py'
config.JobType.maxMemoryMB = 4000
config.JobType.numCores = 4

config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
config.Data.inputDataset =inputDataset_
#config.Data.userInputFiles = open('%s'%inputProcess_).readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10 
#config.Data.outputPrimaryDataset = outputDataset_ 

config.Data.ignoreLocality = True
config.Site.whitelist = [
    'T2_AT_Vienna', 'T2_BE_IIHE', 'T2_BE_UCL', 'T2_BR_SPRACE', 'T2_BR_UERJ',
    'T2_CH_CERN', 'T2_CN_Beijing', 'T2_DE_DESY', 'T2_DE_RWTH',
    'T2_EE_Estonia', 'T2_ES_CIEMAT', 'T2_ES_IFCA', 'T2_FI_HIP',
    'T2_FR_IPHC', 'T2_GR_Ioannina', 'T2_HU_Budapest', 'T2_IN_TIFR',
    'T2_IT_Bari', 'T2_IT_Legnaro', 'T2_IT_Pisa', 'T2_IT_Rome',
    'T2_KR_KISTI', 'T2_PK_NCP', 'T2_PL_Cyfronet',
    'T2_PT_NCG_Lisbon', 'T2_RU_IHEP',
    'T2_TR_METU', 'T2_TW_NCHC', 'T2_UA_KIPT',
    'T2_UK_London_Brunel', 'T2_UK_London_IC', 'T2_UK_SGrid_Bristol',
    'T2_UK_SGrid_RALPP', 'T2_US_Caltech', 'T2_US_Florida',
    'T2_US_MIT', 'T2_US_Nebraska', 'T2_US_Purdue', 'T2_US_UCSD',
    'T2_US_Vanderbilt', 'T2_US_Wisconsin'
]

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/dataRun3'
config.Data.publication = True 
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
