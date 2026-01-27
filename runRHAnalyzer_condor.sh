#!/bin/bash

echo "Starting job on " `date` #Date/time of start of job
echo "Running on: `uname -a`" #Condor job is running on this node
echo "System software: `cat /etc/redhat-release`" #Operating System on that node
# bring in the tarball you created before with caches and large files excluded:
xrdcp -s root://cmseos.fnal.gov//store/group/lpcml/rchudasa/CMSSW_13_0_14.tar.gz .
source /cvmfs/cms.cern.ch/cmsset_default.sh 
tar -xf CMSSW_13_0_14.tgz
rm CMSSW_13_0_14.tgz ###note if you do this locally you remove possibly IMPORTANT FILES always be careful with "rm"
cd CMSSW_13_0_14/src/
scramv1 b MLAnalyzerRun3 # this handles linking the already compiled code - do NOT recompile
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers
echo $CMSSW_BASE "is the CMSSW we have on the local worker node"
export _CONDOR_SCRATCH_DIR=${PWD}
cd ${_CONDOR_SCRATCH_DIR}

inputFile=`sed "${1}q;d" /uscms/home/rchudasa/nobackup/run3Analysis/CMSSW_13_0_14/src/MLAnalyzerRun3/input_list.txt`
outputFile="MLNtuples_ATauTau_${1}.root"

if [ -s ${outputFile} ]
then
  echo "File already exists, skipping"
else
  echo "File doesn't exist or is empty - running"
  cmsRun RecHitAnalyzer/python/ConfFile_cfg.py inputFiles_load=$inputFile outputFile=$outputFile
fi

echo ">> Copying output to EOS..."
xrdcp -f $outputFile root://cmseos.fnal.gov//store/group/lpcml/rchudasa/MLAnalyzer_ntuples/ATauTau_physicalMass

rm -f $outputFile

echo ">> Job ${1} finished successfully."
