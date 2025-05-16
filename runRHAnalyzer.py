import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/14_miniAODSIM_RAWAOD-RHv4/250513_023659/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_1.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_1.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_2.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_3.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_4.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_5.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_6.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_7.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_8.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_9.root,file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM_RAWAOD-RHv3/250509_193418/0000/MINIAOD_HToAATo4Tau_RAW-AODRH_10.root'
#inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/CMSSW_13_0_14/src/MCProduction/E2E-HToAATo4Tau/MINIAOD_HToAATo4Tau_RAW-AODRH_M3p7.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_miniAODSIM/250318_200049/0000/step4_HAA4Tau_3p7_MiniAOD_allevent_509.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M14_Run3_2023/14_miniAODSIM/250318_201016/0000/step4_HAA4Tau_3p7_MiniAOD_allevent_863.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_AODSIM_hadronic/240930_055332/0000/step3_AODSIM_100.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/GluGluHToTauTau_M-125_TuneCP5_13p6TeV_powheg-pythia8/HTauTau_AODSIM_oneBlock/241115_165719/0000/step3_AODSIM_117.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_AODSIM_newBigProd/250113_144409/0000/step3_AODSIM_2.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/TT_TuneCP5_13p6TeV_powheg-pythia8/TTbar_AODSIM_oneBlock/241115_165608/0000/step3_AODSIM_117.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/dataRun3/Tau/Tau_Run2023C_RAW-AOD_EventAwareLumiv2/241203_111828/0000/RAW2DIGI_L1Reco_RECO_32.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/DYto2L_M-50_TuneCP5_13p6TeV_pythia8/DYto2L_AODSIM_multiThreads/250213_062952/0000/step3_AODSIM_41.root'#DYTo2L

maxEvents_=5
#maxEvents_=500
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar.root'
outputFile_='HAA4tau_M14_RAWRH.root'
#outputFile_='HAA4tau_M3p7_RAWRH.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'
#outputFile_='Ato2Tau_massreg_sample.root'
#outputFile_='data_Tautrigger.root'
#outputFile_='dy2Ldy2L_Tautrigger.root'
#outputFile_='h2aa4Tau_M14_Tautrigger.root'
#outputFile_='ttbar_Tautrigger.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
