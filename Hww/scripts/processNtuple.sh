#!/bin/sh

###############################################################################################
###7TeV Samples
###############################################################################################

###############################################################################################
###WW Analysis Ntuple
###############################################################################################

#do regular merging among different datasets
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/HwwAnalysis/cern/filler/015/ --OutputPath=/home/sixie/hist/HwwAnalysis/merged/ --FilenameHeader=HwwAnalysis --DatasetListFile=$CMSSW_BASE/src/EWKAna/Hww/scripts/WW7TevSampleList.txt
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/HwwAnalysis/cern/filefi/017/ --OutputPath=/home/sixie/hist/HwwAnalysis/merged/ --FilenameHeader=HwwAnalysis --DatasetListFile=$CMSSW_BASE/src/EWKAna/Hww/scripts/WW7TevSampleList.txt
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/HwwAnalysis/cern/filefi/019/ --OutputPath=/home/sixie/hist/HwwAnalysis/merged/ --FilenameHeader=HwwAnalysis --DatasetListFile=$CMSSW_BASE/src/EWKAna/Hww/scripts/WW7TevSampleList.txt
python $CMSSW_BASE/src/MitHiggs/scripts/MergeFilesets.py --InputPath=/home/sixie/hist/HwwAnalysis/cern/filefi/020/ --OutputPath=/home/sixie/hist/HwwAnalysis/merged/ --FilenameHeader=HwwAnalysis --DatasetListFile=$CMSSW_BASE/src/EWKAna/Hww/scripts/WW7TevSampleList.txt

#do normalization first
python $CMSSW_BASE/src/EWKAna/Hww/scripts/NormalizeWWNtuple.py --InputPath=/home/sixie/hist/HwwAnalysis/merged/ --OutputPath=/home/sixie/hist/HwwAnalysis/normalized/ --FilenameHeader=HwwAnalysis --DatasetListFile=$CMSSW_BASE/src/EWKAna/Hww/scripts/WW7TevSampleList.txt
