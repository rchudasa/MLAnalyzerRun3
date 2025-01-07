from CRABClient.UserUtilities import config
config = config()

Process = '3p7'
#Process = '10'

inputDataset_ ={
        '3p7': '/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/lpcml-3p7_AODSIM_ignoreLocality-953b1873547799e513f8a43f2c57e3b2/USER',
        '4': '/HToAATo4Tau_hadronic_tauDecay_M4_Run3_2023/lpcml-4_AODSIM_ignoreLocality-953b1873547799e513f8a43f2c57e3b2/USER',
        '6': '/HToAATo4Tau_hadronic_tauDecay_M6_Run3_2023/lpcml-6_AODSIM_ignoreLocality-953b1873547799e513f8a43f2c57e3b2/USER',
        '14':'/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/phys_diffraction-14_AODSIM_multiThreads_8Gb_new-953b1873547799e513f8a43f2c57e3b2/USER',
        'QCD':'/GEN_SIM_QCD_pt15to7000_Run3Summer23GS/lpcml-AOD_QCD_pt15to7000_Run3Summer23GS-953b1873547799e513f8a43f2c57e3b2/USER',
        'WJets': "/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/lpcml-WJets_AODSIM_multiThreads-953b1873547799e513f8a43f2c57e3b2/USER",
        'DYTo2L':'/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/lpcml-DYto2L_AODSIM_check-953b1873547799e513f8a43f2c57e3b2/USER',
        'TTbar': "/TT_TuneCP5_13p6TeV_powheg-pythia8/lpcml-TTbar_AODSIM_oneBlock-953b1873547799e513f8a43f2c57e3b2/USER",
        'HTauTau':'/GluGluHToTauTau_M-125_TuneCP5_13p6TeV_powheg-pythia8/lpcml-HTauTau_AODSIM_oneBlock-953b1873547799e513f8a43f2c57e3b2/USER'
        }.get(Process, None)

outputDataset_ = {
        #'3p7':'HToAATo4Tau_M_3p7_pythia8_2018UL_AOD'
        #'8':'HToAATo4Tau_M_8_pythia8_2018UL_AOD',
        '12':'HToAATo4Tau_M_12_pythia8_2018UL_AOD',
        #'WJets':'WJetsToLNu_TuneCP5_13TeV_madgraphMLM-pythia8',
        #'TTbar':'TTToHadronic_TuneCP5_13TeV_powheg-pythia8'
        }.get(Process, None)

#config.section_('General')
config.General.requestName = '%s_MLAnalyzer_bigProduction'%Process
config.General.workArea = 'crab_bigProduction'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'RecHitAnalyzer/python/ConfFile_cfg.py'
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

config.Data.outLFNDirBase = '/store/group/lpcml/rchudasa/MCGenerationRun3'
config.Data.publication = True 
config.Site.storageSite = 'T3_US_FNALLPC'
config.Data.outputDatasetTag = config.General.requestName
