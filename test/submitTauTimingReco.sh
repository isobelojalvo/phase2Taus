#!/bin/sh                                                                                                                                                                                                                                                                                                                     
#voms-proxy-init --voms cms --valid 100:00                                                                                                                                                                                                                                                                                   
cat sub-runTauTiming-Reco.py > SUBtiming.py
cat submit.py >> SUBtiming.py


for dir in DY0Jets_0PU  #DY0Jets_200PU QCD_0PU QCD_200PU #TTBar-RelVal-0  TTBar-RelVal-140  TTBar-RelVal-200 ZTT-RelVal-0  ZTT-RelVal-140   ZTT-RelVal-200  QCDMuEnriched_realistic
do 
    echo " "
    echo "====================" $dir "========================="

    rm -r /nfs_scratch/ojalvo/$1-$dir-SUBtiming/
    mkdir /nfs_scratch/ojalvo/$1-$dir-SUBtiming/
    
#make dag dir
    mkdir -p /nfs_scratch/ojalvo/$1-$dir-SUBtiming/dags
    mkdir -p /nfs_scratch/ojalvo/$1-$dir-SUBtiming/dags/daginputs
## outputdir = srm://cmssrm.hep.wisc.edu:8443/srm/$1/server?SFN=/hdfs/store/user/ojalvo/$1-$dir-SUBtiming/
    
    farmoutAnalysisJobs --assume-input-files-exist --input-basenames-not-unique  --input-file-list=samples-reco/$dir.txt \
	--submit-dir=/nfs_scratch/ojalvo/$1-$dir-SUBtiming/submit \
	--output-dag-file=/nfs_scratch/ojalvo/$1-$dir-SUBtiming/dags/dag \
	$1-$dir  \
	$CMSSW_BASE  \
	$CMSSW_BASE/src/RecoTauTag/phase2Taus/test/SUBtiming.py 

done

rm SUBtiming.py