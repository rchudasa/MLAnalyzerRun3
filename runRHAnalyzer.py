import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='/store/group/lpcml/rchudasa/MCGenerationRun3/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/DY2L_miniAODSIM_RAWAOD-RecHits/250607_220044/0000/MINIAODSIM_RAWAOD_RecHits_M3p7_1.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/14_miniAODSIM_RAWAOD-RHv4/250513_023659/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_1.root'
inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/CMSSW_13_0_14/src/MCProduction/E2E-HToAATo4Tau/MINIAODSIM_RAWAOD_RecHits_ATauTau.root'
#inputFiles_='/store/group/lpcml/rchudasa/MCGenerationRun3/GEN_SIM_ATo2Tau_m3p6To18_pt30To300_v2/ATauTau_miniAODSIM_RAWAOD-RecHits-7Jan2026/260107_155633/0001/MINIAODSIM_RAWAOD_RecHits_ATauTau_1009.root'
#inputFiles_='/store/group/lpcml/bbbam/MCGeneration_run3/HTo2Tau_hadronic/HTo2Tau_hadronic_GEN_SIM/HTo2Tau_hadronic_MiniAOD_v2/251221_190902/0000/MINIAODSIM_RAWAOD_RecHits_HTo2Tau_hadronic_1.root'
#inputFiles_='/store/group/lpcml/rchudasa/MCGenerationRun3/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/DY2L_miniAODSIM_RAWAOD-RecHits/250607_220044/0000/MINIAODSIM_RAWAOD_RecHits_M3p7_4.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/miniAOD_allRecHits/MINIAODSIM_RAWAOD_RecHits_M3p7_1635_QCD.root'
#maxEvents_=1000
#maxEvents_=20
maxEvents_=-1
skipEvents_=0#
outputFile_='aTauTau_Genchecks_final_motherPDGID_checked.root'
#outputFile_='HTauTau_hadronicChecks.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
