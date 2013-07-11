#!/bin/sh

###############################################################################################
###7TeV Samples
###############################################################################################

###############################################################################################
###WW Analysis Ntuple
###############################################################################################

#do regular merging among different datasets
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/AllNtuple/cern/filefi/023/HWWNtuple/ --OutputPath=/home/sixie/hist/AllNtuple/cern/filefi/023/HWWNtuple/ --FilenameHeader=AllNtuple --DatasetListFile=$CMSSW_BASE/src/EWKAna/HZZ4l/scripts/HZZ4lSampleList.txt

python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/AllNtuple/HZZ4lNtuple/mc/ --OutputPath=/home/sixie/hist/AllNtuple/HZZ4lNtuple/mc/ --FilenameHeader=AllNtuple --DatasetListFile=$CMSSW_BASE/src/EWKAna/HZZ4l/scripts/HZZ4lSampleList.txt


#do normalization first
python $CMSSW_BASE/src/EWKAna/HZZ4l/scripts/NormalizeZZ4lNtuple.py --InputPath=/home/sixie/hist/AllNtuple/HZZ4lNtuple/mc/ --OutputPath=/data/blue/sixie/ntuples/HZZ4l/mc/ --FilenameHeader=AllNtuple --DatasetListFile=$CMSSW_BASE/src/EWKAna/HZZ4l/scripts/HZZ4lSampleList.txt
