import os

#cfg='RecHitAnalyzer/python/ConfFile_data_cfg.py'
cfg='RecHitAnalyzer/python/ConfFile_cfg.py'
#inputFiles_='file:/afs/cern.ch/work/r/rchudasa/private/TauClassification/run3/CMSSW_13_0_14/src/MCProduction/E2E-HToAATo4Tau/step3_AODSIM_HTauTau_checkrechit_collection.root'
#inputFiles_='file:/eos/cms/store/group/phys_diffraction/rchudasa/MCGeneration/HToAATo4Tau_hadronic_tauDecay_M3p7_Run3_2023/3p7_AODSIM_hadronic/240930_055332/0000/step3_AODSIM_100.root'
#inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/GluGluHToTauTau_M-125_TuneCP5_13p6TeV_powheg-pythia8/HTauTau_AODSIM_oneBlock/241115_165719/0000/step3_AODSIM_117.root'
inputFiles_='root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MCGenerationRun3/TT_TuneCP5_13p6TeV_powheg-pythia8/TTbar_AODSIM_oneBlock/241115_165608/0000/step3_AODSIM_117.root'

#maxEvents_=10
maxEvents_=20
#maxEvents_=-1
skipEvents_=0#
#outputFile_='MLAnal_PhaseI_TTbar_13TeVu_trackRefitter.root'
#outputFile_='GJet.root'
#outputFile_='ttbar_secVertex.root'
#outputFile_='DYToTauTau_subJet.root'
#outputFile_='WJets_secVertex.root'
#outputFile_='dyToEE.root'
#outputFile_='acd_EmEnriched.root'
#outputFile_='Ato2Tau_massreg_sample.root'
outputFile_='ttbar_trigger.root'

# cmd="cmsTraceExceptions cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
cmd="cmsRun %s inputFiles=%s maxEvents=%d skipEvents=%d outputFile=%s"%(cfg,inputFiles_,maxEvents_,skipEvents_,outputFile_)
print(f"{cmd}")
os.system(cmd)
