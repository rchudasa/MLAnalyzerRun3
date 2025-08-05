import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='/store/group/lpcml/rchudasa/MCGenerationRun3/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/DY2L_miniAODSIM_RAWAOD-RecHits/250607_220044/0000/MINIAODSIM_RAWAOD_RecHits_M3p7_1.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/14_miniAODSIM_RAWAOD-RHv4/250513_023659/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_1.root'
inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/miniAOD_allRecHits/MINIAODSIM_RAWAOD_RecHits_M3p7_1635_QCD.root'
maxEvents_=10
#maxEvents_=-1
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar.root'
#outputFile_='HAA4tau_M14_RAWRH.root'
#outputFile_='HAA4tau_M14_allRH_sameNtuples_reducedEgammaHCALRH_allTaupt30.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'
#outputFile_='Ato2Tau_massreg_sample.root'
#outputFile_='data_Tautrigger.root'
outputFile_='qcd_miniAOD_EBRHOnly.root'
#outputFile_='h2aa4Tau_M14_Tautrigger.root'
#outputFile_='ttbar_Tautrigger.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
