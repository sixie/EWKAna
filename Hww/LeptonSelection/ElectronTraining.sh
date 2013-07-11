#####################################################################################
#Make Electron Ntuples
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeRealElectronTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeFakeElectronTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeMCElectronTrainingNtuple.C+

#####################################################################################
#Do pt reweighting
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeElectronPtSpectrum.C+
root -l -b -q EWKAna/Hww/LeptonSelection/ReweightElectronPU.C+

#####################################################################################
#Skim Ntuples, Split into bins
#####################################################################################

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0LowPt",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1LowPt",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPt",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0HighPt",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1HighPt",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2HighPt",5)'

#####################################################################################
#MVA Training variables to electron ntuple
#####################################################################################

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V0","V0",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V0","V0",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V0","V0",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V0","V0",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V0","V0",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V0","V0",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V1","V1",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V1","V1",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V1","V1",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V1","V1",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V1","V1",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V1","V1",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V2","V2",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V2","V2",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V2","V2",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V2","V2",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V2","V2",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V2","V2",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V3","V3",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V3","V3",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V3","V3",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V3","V3",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V3","V3",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V3","V3",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V4","V4",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V4","V4",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V4","V4",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V4","V4",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V4","V4",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V4","V4",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V5","V5",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V5","V5",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V5","V5",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V5","V5",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V5","V5",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V5","V5",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V6","V6",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V6","V6",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V6","V6",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V6","V6",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V6","V6",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V6","V6",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V7","V7",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V7","V7",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V7","V7",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V7","V7",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V7","V7",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V7","V7",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V9","V9",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V9","V9",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V9","V9",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V9","V9",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V9","V9",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V9","V9",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V10","V10",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V10","V10",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V10","V10",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V10","V10",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V10","V10",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V10","V10",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V11","V11",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V11","V11",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V11","V11",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V11","V11",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V11","V11",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V11","V11",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V12","V12",kFALSE,0,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V12","V12",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V12","V12",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V12","V12",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V12","V12",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V12","V12",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V13","V13",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V13","V13",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V13","V13",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V13","V13",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V13","V13",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V13","V13",kFALSE,5,"BDTG,Likelihood")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V14","V14",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V14","V14",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V14","V14",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V14","V14",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V14","V14",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V14","V14",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V15","V15",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V15","V15",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V15","V15",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V15","V15",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V15","V15",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V15","V15",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V16","V16",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V16","V16",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V16","V16",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V16","V16",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V16","V16",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V16","V16",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V17","V17",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V17","V17",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V17","V17",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V17","V17",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V17","V17",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V17","V17",kFALSE,5,"BDTG,Likelihood")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","Subdet0LowPt_V18","V18",kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","Subdet1LowPt_V18","V18",kFALSE,1,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","Subdet2LowPt_V18","V18",kFALSE,2,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","Subdet0HighPt_V18","V18",kFALSE,3,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","Subdet1HighPt_V18","V18",kFALSE,4,"BDTG,Likelihood")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","Subdet2HighPt_V18","V18",kFALSE,5,"BDTG,Likelihood")'





#####################################################################################
#MVA Evaluate
#####################################################################################

#Fakes
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0LowPt.root","output/ElectronNtuple.Fake.Subdet0LowPtV0.root","Subdet0LowPt_V0","V0",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV0.root","output/ElectronNtuple.Fake.Subdet0LowPtV1.root","Subdet0LowPt_V1","V1",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV1.root","output/ElectronNtuple.Fake.Subdet0LowPtV2.root","Subdet0LowPt_V2","V2",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV2.root","output/ElectronNtuple.Fake.Subdet0LowPtV3.root","Subdet0LowPt_V3","V3",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV3.root","output/ElectronNtuple.Fake.Subdet0LowPtV4.root","Subdet0LowPt_V4","V4",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV4.root","output/ElectronNtuple.Fake.Subdet0LowPtV5.root","Subdet0LowPt_V5","V5",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV5.root","output/ElectronNtuple.Fake.Subdet0LowPtV6.root","Subdet0LowPt_V6","V6",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV6.root","output/ElectronNtuple.Fake.Subdet0LowPtV7.root","Subdet0LowPt_V7","V7",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV7.root","output/ElectronNtuple.Fake.Subdet0LowPtV8.root","Subdet0LowPt_V8","V8",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV8.root","output/ElectronNtuple.Fake.Subdet0LowPtV9.root","Subdet0LowPt_V9","V9",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV9.root","output/ElectronNtuple.Fake.Subdet0LowPtV10.root","Subdet0LowPt_V10","V10",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV10.root","output/ElectronNtuple.Fake.Subdet0LowPtV11.root","Subdet0LowPt_V11","V11",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV11.root","output/ElectronNtuple.Fake.Subdet0LowPtV12.root","Subdet0LowPt_V12","V12",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV12.root","output/ElectronNtuple.Fake.Subdet0LowPtV13.root","Subdet0LowPt_V13","V13",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV13.root","output/ElectronNtuple.Fake.Subdet0LowPtV14.root","Subdet0LowPt_V14","V14",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV14.root","output/ElectronNtuple.Fake.Subdet0LowPtV15.root","Subdet0LowPt_V15","V15",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV15.root","output/ElectronNtuple.Fake.Subdet0LowPtV16.root","Subdet0LowPt_V16","V16",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV16.root","output/ElectronNtuple.Fake.Subdet0LowPtV17.root","Subdet0LowPt_V17","V17",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtV17.root","output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt_V18","V18",kFALSE,kFALSE,0,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1LowPt.root","output/ElectronNtuple.Fake.Subdet1LowPtV0.root","Subdet1LowPt_V0","V0",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV0.root","output/ElectronNtuple.Fake.Subdet1LowPtV1.root","Subdet1LowPt_V1","V1",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV1.root","output/ElectronNtuple.Fake.Subdet1LowPtV2.root","Subdet1LowPt_V2","V2",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV2.root","output/ElectronNtuple.Fake.Subdet1LowPtV3.root","Subdet1LowPt_V3","V3",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV3.root","output/ElectronNtuple.Fake.Subdet1LowPtV4.root","Subdet1LowPt_V4","V4",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV4.root","output/ElectronNtuple.Fake.Subdet1LowPtV5.root","Subdet1LowPt_V5","V5",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV5.root","output/ElectronNtuple.Fake.Subdet1LowPtV6.root","Subdet1LowPt_V6","V6",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV6.root","output/ElectronNtuple.Fake.Subdet1LowPtV7.root","Subdet1LowPt_V7","V7",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV7.root","output/ElectronNtuple.Fake.Subdet1LowPtV8.root","Subdet1LowPt_V8","V8",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV8.root","output/ElectronNtuple.Fake.Subdet1LowPtV9.root","Subdet1LowPt_V9","V9",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV9.root","output/ElectronNtuple.Fake.Subdet1LowPtV10.root","Subdet1LowPt_V10","V10",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV10.root","output/ElectronNtuple.Fake.Subdet1LowPtV11.root","Subdet1LowPt_V11","V11",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV11.root","output/ElectronNtuple.Fake.Subdet1LowPtV12.root","Subdet1LowPt_V12","V12",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV12.root","output/ElectronNtuple.Fake.Subdet1LowPtV13.root","Subdet1LowPt_V13","V13",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV13.root","output/ElectronNtuple.Fake.Subdet1LowPtV14.root","Subdet1LowPt_V14","V14",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV14.root","output/ElectronNtuple.Fake.Subdet1LowPtV15.root","Subdet1LowPt_V15","V15",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV15.root","output/ElectronNtuple.Fake.Subdet1LowPtV16.root","Subdet1LowPt_V16","V16",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV16.root","output/ElectronNtuple.Fake.Subdet1LowPtV17.root","Subdet1LowPt_V17","V17",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtV17.root","output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt_V18","V18",kFALSE,kFALSE,1,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2LowPt.root","output/ElectronNtuple.Fake.Subdet2LowPtV0.root","Subdet2LowPt_V0","V0",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV0.root","output/ElectronNtuple.Fake.Subdet2LowPtV1.root","Subdet2LowPt_V1","V1",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV1.root","output/ElectronNtuple.Fake.Subdet2LowPtV2.root","Subdet2LowPt_V2","V2",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV2.root","output/ElectronNtuple.Fake.Subdet2LowPtV3.root","Subdet2LowPt_V3","V3",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV3.root","output/ElectronNtuple.Fake.Subdet2LowPtV4.root","Subdet2LowPt_V4","V4",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV4.root","output/ElectronNtuple.Fake.Subdet2LowPtV5.root","Subdet2LowPt_V5","V5",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV5.root","output/ElectronNtuple.Fake.Subdet2LowPtV6.root","Subdet2LowPt_V6","V6",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV6.root","output/ElectronNtuple.Fake.Subdet2LowPtV7.root","Subdet2LowPt_V7","V7",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV7.root","output/ElectronNtuple.Fake.Subdet2LowPtV8.root","Subdet2LowPt_V8","V8",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV8.root","output/ElectronNtuple.Fake.Subdet2LowPtV9.root","Subdet2LowPt_V9","V9",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV9.root","output/ElectronNtuple.Fake.Subdet2LowPtV10.root","Subdet2LowPt_V10","V10",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV10.root","output/ElectronNtuple.Fake.Subdet2LowPtV11.root","Subdet2LowPt_V11","V11",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV11.root","output/ElectronNtuple.Fake.Subdet2LowPtV12.root","Subdet2LowPt_V12","V12",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV12.root","output/ElectronNtuple.Fake.Subdet2LowPtV13.root","Subdet2LowPt_V13","V13",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV13.root","output/ElectronNtuple.Fake.Subdet2LowPtV14.root","Subdet2LowPt_V14","V14",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV14.root","output/ElectronNtuple.Fake.Subdet2LowPtV15.root","Subdet2LowPt_V15","V15",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV15.root","output/ElectronNtuple.Fake.Subdet2LowPtV16.root","Subdet2LowPt_V16","V16",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV16.root","output/ElectronNtuple.Fake.Subdet2LowPtV17.root","Subdet2LowPt_V17","V17",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtV17.root","output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt_V18","V18",kFALSE,kFALSE,2,"BDTG")'




root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet0HighPt.root","output/ElectronNtuple.Fake.Subdet0HighPtV0.root","Subdet0HighPt_V0","V0",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV0.root","output/ElectronNtuple.Fake.Subdet0HighPtV1.root","Subdet0HighPt_V1","V1",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV1.root","output/ElectronNtuple.Fake.Subdet0HighPtV2.root","Subdet0HighPt_V2","V2",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV2.root","output/ElectronNtuple.Fake.Subdet0HighPtV3.root","Subdet0HighPt_V3","V3",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV3.root","output/ElectronNtuple.Fake.Subdet0HighPtV4.root","Subdet0HighPt_V4","V4",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV4.root","output/ElectronNtuple.Fake.Subdet0HighPtV5.root","Subdet0HighPt_V5","V5",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV5.root","output/ElectronNtuple.Fake.Subdet0HighPtV6.root","Subdet0HighPt_V6","V6",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV6.root","output/ElectronNtuple.Fake.Subdet0HighPtV7.root","Subdet0HighPt_V7","V7",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV7.root","output/ElectronNtuple.Fake.Subdet0HighPtV8.root","Subdet0HighPt_V8","V8",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV8.root","output/ElectronNtuple.Fake.Subdet0HighPtV9.root","Subdet0HighPt_V9","V9",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV9.root","output/ElectronNtuple.Fake.Subdet0HighPtV10.root","Subdet0HighPt_V10","V10",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV10.root","output/ElectronNtuple.Fake.Subdet0HighPtV11.root","Subdet0HighPt_V11","V11",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV11.root","output/ElectronNtuple.Fake.Subdet0HighPtV12.root","Subdet0HighPt_V12","V12",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV12.root","output/ElectronNtuple.Fake.Subdet0HighPtV13.root","Subdet0HighPt_V13","V13",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV13.root","output/ElectronNtuple.Fake.Subdet0HighPtV14.root","Subdet0HighPt_V14","V14",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV14.root","output/ElectronNtuple.Fake.Subdet0HighPtV15.root","Subdet0HighPt_V15","V15",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV15.root","output/ElectronNtuple.Fake.Subdet0HighPtV16.root","Subdet0HighPt_V16","V16",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV16.root","output/ElectronNtuple.Fake.Subdet0HighPtV17.root","Subdet0HighPt_V17","V17",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtV17.root","output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt_V18","V18",kFALSE,kFALSE,3,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet1HighPt.root","output/ElectronNtuple.Fake.Subdet1HighPtV0.root","Subdet1HighPt_V0","V0",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV0.root","output/ElectronNtuple.Fake.Subdet1HighPtV1.root","Subdet1HighPt_V1","V1",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV1.root","output/ElectronNtuple.Fake.Subdet1HighPtV2.root","Subdet1HighPt_V2","V2",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV2.root","output/ElectronNtuple.Fake.Subdet1HighPtV3.root","Subdet1HighPt_V3","V3",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV3.root","output/ElectronNtuple.Fake.Subdet1HighPtV4.root","Subdet1HighPt_V4","V4",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV4.root","output/ElectronNtuple.Fake.Subdet1HighPtV5.root","Subdet1HighPt_V5","V5",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV5.root","output/ElectronNtuple.Fake.Subdet1HighPtV6.root","Subdet1HighPt_V6","V6",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV6.root","output/ElectronNtuple.Fake.Subdet1HighPtV7.root","Subdet1HighPt_V7","V7",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV7.root","output/ElectronNtuple.Fake.Subdet1HighPtV8.root","Subdet1HighPt_V8","V8",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV8.root","output/ElectronNtuple.Fake.Subdet1HighPtV9.root","Subdet1HighPt_V9","V9",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV9.root","output/ElectronNtuple.Fake.Subdet1HighPtV10.root","Subdet1HighPt_V10","V10",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV10.root","output/ElectronNtuple.Fake.Subdet1HighPtV11.root","Subdet1HighPt_V11","V11",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV11.root","output/ElectronNtuple.Fake.Subdet1HighPtV12.root","Subdet1HighPt_V12","V12",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV12.root","output/ElectronNtuple.Fake.Subdet1HighPtV13.root","Subdet1HighPt_V13","V13",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV13.root","output/ElectronNtuple.Fake.Subdet1HighPtV14.root","Subdet1HighPt_V14","V14",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV14.root","output/ElectronNtuple.Fake.Subdet1HighPtV15.root","Subdet1HighPt_V15","V15",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV15.root","output/ElectronNtuple.Fake.Subdet1HighPtV16.root","Subdet1HighPt_V16","V16",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV16.root","output/ElectronNtuple.Fake.Subdet1HighPtV17.root","Subdet1HighPt_V17","V17",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtV17.root","output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt_V18","V18",kFALSE,kFALSE,4,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.PtAndPUWeighted.Subdet2HighPt.root","output/ElectronNtuple.Fake.Subdet2HighPtV0.root","Subdet2HighPt_V0","V0",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV0.root","output/ElectronNtuple.Fake.Subdet2HighPtV1.root","Subdet2HighPt_V1","V1",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV1.root","output/ElectronNtuple.Fake.Subdet2HighPtV2.root","Subdet2HighPt_V2","V2",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV2.root","output/ElectronNtuple.Fake.Subdet2HighPtV3.root","Subdet2HighPt_V3","V3",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV3.root","output/ElectronNtuple.Fake.Subdet2HighPtV4.root","Subdet2HighPt_V4","V4",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV4.root","output/ElectronNtuple.Fake.Subdet2HighPtV5.root","Subdet2HighPt_V5","V5",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV5.root","output/ElectronNtuple.Fake.Subdet2HighPtV6.root","Subdet2HighPt_V6","V6",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV6.root","output/ElectronNtuple.Fake.Subdet2HighPtV7.root","Subdet2HighPt_V7","V7",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV7.root","output/ElectronNtuple.Fake.Subdet2HighPtV8.root","Subdet2HighPt_V8","V8",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV8.root","output/ElectronNtuple.Fake.Subdet2HighPtV9.root","Subdet2HighPt_V9","V9",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV9.root","output/ElectronNtuple.Fake.Subdet2HighPtV10.root","Subdet2HighPt_V10","V10",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV10.root","output/ElectronNtuple.Fake.Subdet2HighPtV11.root","Subdet2HighPt_V11","V11",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV11.root","output/ElectronNtuple.Fake.Subdet2HighPtV12.root","Subdet2HighPt_V12","V12",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV12.root","output/ElectronNtuple.Fake.Subdet2HighPtV13.root","Subdet2HighPt_V13","V13",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV13.root","output/ElectronNtuple.Fake.Subdet2HighPtV14.root","Subdet2HighPt_V14","V14",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV14.root","output/ElectronNtuple.Fake.Subdet2HighPtV15.root","Subdet2HighPt_V15","V15",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV15.root","output/ElectronNtuple.Fake.Subdet2HighPtV16.root","Subdet2HighPt_V16","V16",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV16.root","output/ElectronNtuple.Fake.Subdet2HighPtV17.root","Subdet2HighPt_V17","V17",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtV17.root","output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt_V18","V18",kFALSE,kFALSE,5,"BDTG")'



#Real
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","output/ElectronNtuple.Real.Subdet0LowPtV0.root","Subdet0LowPt_V0","V0",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV0.root","output/ElectronNtuple.Real.Subdet0LowPtV1.root","Subdet0LowPt_V1","V1",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV1.root","output/ElectronNtuple.Real.Subdet0LowPtV2.root","Subdet0LowPt_V2","V2",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV2.root","output/ElectronNtuple.Real.Subdet0LowPtV3.root","Subdet0LowPt_V3","V3",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV3.root","output/ElectronNtuple.Real.Subdet0LowPtV4.root","Subdet0LowPt_V4","V4",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV4.root","output/ElectronNtuple.Real.Subdet0LowPtV5.root","Subdet0LowPt_V5","V5",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV5.root","output/ElectronNtuple.Real.Subdet0LowPtV6.root","Subdet0LowPt_V6","V6",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV6.root","output/ElectronNtuple.Real.Subdet0LowPtV7.root","Subdet0LowPt_V7","V7",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV7.root","output/ElectronNtuple.Real.Subdet0LowPtV8.root","Subdet0LowPt_V8","V8",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV8.root","output/ElectronNtuple.Real.Subdet0LowPtV9.root","Subdet0LowPt_V9","V9",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV9.root","output/ElectronNtuple.Real.Subdet0LowPtV10.root","Subdet0LowPt_V10","V10",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV10.root","output/ElectronNtuple.Real.Subdet0LowPtV11.root","Subdet0LowPt_V11","V11",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV11.root","output/ElectronNtuple.Real.Subdet0LowPtV12.root","Subdet0LowPt_V12","V12",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV12.root","output/ElectronNtuple.Real.Subdet0LowPtV13.root","Subdet0LowPt_V13","V13",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV13.root","output/ElectronNtuple.Real.Subdet0LowPtV14.root","Subdet0LowPt_V14","V14",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV14.root","output/ElectronNtuple.Real.Subdet0LowPtV15.root","Subdet0LowPt_V15","V15",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV15.root","output/ElectronNtuple.Real.Subdet0LowPtV16.root","Subdet0LowPt_V16","V16",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV16.root","output/ElectronNtuple.Real.Subdet0LowPtV17.root","Subdet0LowPt_V17","V17",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtV17.root","output/ElectronNtuple.Real.Subdet0LowPtV18.root","Subdet0LowPt_V18","V18",kFALSE,kFALSE,0,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","output/ElectronNtuple.Real.Subdet1LowPtV0.root","Subdet1LowPt_V0","V0",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV0.root","output/ElectronNtuple.Real.Subdet1LowPtV1.root","Subdet1LowPt_V1","V1",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV1.root","output/ElectronNtuple.Real.Subdet1LowPtV2.root","Subdet1LowPt_V2","V2",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV2.root","output/ElectronNtuple.Real.Subdet1LowPtV3.root","Subdet1LowPt_V3","V3",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV3.root","output/ElectronNtuple.Real.Subdet1LowPtV4.root","Subdet1LowPt_V4","V4",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV4.root","output/ElectronNtuple.Real.Subdet1LowPtV5.root","Subdet1LowPt_V5","V5",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV5.root","output/ElectronNtuple.Real.Subdet1LowPtV6.root","Subdet1LowPt_V6","V6",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV6.root","output/ElectronNtuple.Real.Subdet1LowPtV7.root","Subdet1LowPt_V7","V7",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV7.root","output/ElectronNtuple.Real.Subdet1LowPtV8.root","Subdet1LowPt_V8","V8",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV8.root","output/ElectronNtuple.Real.Subdet1LowPtV9.root","Subdet1LowPt_V9","V9",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV9.root","output/ElectronNtuple.Real.Subdet1LowPtV10.root","Subdet1LowPt_V10","V10",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV10.root","output/ElectronNtuple.Real.Subdet1LowPtV11.root","Subdet1LowPt_V11","V11",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV11.root","output/ElectronNtuple.Real.Subdet1LowPtV12.root","Subdet1LowPt_V12","V12",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV12.root","output/ElectronNtuple.Real.Subdet1LowPtV13.root","Subdet1LowPt_V13","V13",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV13.root","output/ElectronNtuple.Real.Subdet1LowPtV14.root","Subdet1LowPt_V14","V14",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV14.root","output/ElectronNtuple.Real.Subdet1LowPtV15.root","Subdet1LowPt_V15","V15",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV15.root","output/ElectronNtuple.Real.Subdet1LowPtV16.root","Subdet1LowPt_V16","V16",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV16.root","output/ElectronNtuple.Real.Subdet1LowPtV17.root","Subdet1LowPt_V17","V17",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtV17.root","output/ElectronNtuple.Real.Subdet1LowPtV18.root","Subdet1LowPt_V18","V18",kFALSE,kFALSE,1,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","output/ElectronNtuple.Real.Subdet2LowPtV0.root","Subdet2LowPt_V0","V0",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV0.root","output/ElectronNtuple.Real.Subdet2LowPtV1.root","Subdet2LowPt_V1","V1",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV1.root","output/ElectronNtuple.Real.Subdet2LowPtV2.root","Subdet2LowPt_V2","V2",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV2.root","output/ElectronNtuple.Real.Subdet2LowPtV3.root","Subdet2LowPt_V3","V3",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV3.root","output/ElectronNtuple.Real.Subdet2LowPtV4.root","Subdet2LowPt_V4","V4",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV4.root","output/ElectronNtuple.Real.Subdet2LowPtV5.root","Subdet2LowPt_V5","V5",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV5.root","output/ElectronNtuple.Real.Subdet2LowPtV6.root","Subdet2LowPt_V6","V6",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV6.root","output/ElectronNtuple.Real.Subdet2LowPtV7.root","Subdet2LowPt_V7","V7",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV7.root","output/ElectronNtuple.Real.Subdet2LowPtV8.root","Subdet2LowPt_V8","V8",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV8.root","output/ElectronNtuple.Real.Subdet2LowPtV9.root","Subdet2LowPt_V9","V9",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV9.root","output/ElectronNtuple.Real.Subdet2LowPtV10.root","Subdet2LowPt_V10","V10",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV10.root","output/ElectronNtuple.Real.Subdet2LowPtV11.root","Subdet2LowPt_V11","V11",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV11.root","output/ElectronNtuple.Real.Subdet2LowPtV12.root","Subdet2LowPt_V12","V12",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV12.root","output/ElectronNtuple.Real.Subdet2LowPtV13.root","Subdet2LowPt_V13","V13",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV13.root","output/ElectronNtuple.Real.Subdet2LowPtV14.root","Subdet2LowPt_V14","V14",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV14.root","output/ElectronNtuple.Real.Subdet2LowPtV15.root","Subdet2LowPt_V15","V15",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV15.root","output/ElectronNtuple.Real.Subdet2LowPtV16.root","Subdet2LowPt_V16","V16",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV16.root","output/ElectronNtuple.Real.Subdet2LowPtV17.root","Subdet2LowPt_V17","V17",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtV17.root","output/ElectronNtuple.Real.Subdet2LowPtV18.root","Subdet2LowPt_V18","V18",kFALSE,kFALSE,2,"BDTG")'




root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","output/ElectronNtuple.Real.Subdet0HighPtV0.root","Subdet0HighPt_V0","V0",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV0.root","output/ElectronNtuple.Real.Subdet0HighPtV1.root","Subdet0HighPt_V1","V1",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV1.root","output/ElectronNtuple.Real.Subdet0HighPtV2.root","Subdet0HighPt_V2","V2",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV2.root","output/ElectronNtuple.Real.Subdet0HighPtV3.root","Subdet0HighPt_V3","V3",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV3.root","output/ElectronNtuple.Real.Subdet0HighPtV4.root","Subdet0HighPt_V4","V4",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV4.root","output/ElectronNtuple.Real.Subdet0HighPtV5.root","Subdet0HighPt_V5","V5",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV5.root","output/ElectronNtuple.Real.Subdet0HighPtV6.root","Subdet0HighPt_V6","V6",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV6.root","output/ElectronNtuple.Real.Subdet0HighPtV7.root","Subdet0HighPt_V7","V7",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV7.root","output/ElectronNtuple.Real.Subdet0HighPtV8.root","Subdet0HighPt_V8","V8",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV8.root","output/ElectronNtuple.Real.Subdet0HighPtV9.root","Subdet0HighPt_V9","V9",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV9.root","output/ElectronNtuple.Real.Subdet0HighPtV10.root","Subdet0HighPt_V10","V10",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV10.root","output/ElectronNtuple.Real.Subdet0HighPtV11.root","Subdet0HighPt_V11","V11",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV11.root","output/ElectronNtuple.Real.Subdet0HighPtV12.root","Subdet0HighPt_V12","V12",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV12.root","output/ElectronNtuple.Real.Subdet0HighPtV13.root","Subdet0HighPt_V13","V13",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV13.root","output/ElectronNtuple.Real.Subdet0HighPtV14.root","Subdet0HighPt_V14","V14",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV14.root","output/ElectronNtuple.Real.Subdet0HighPtV15.root","Subdet0HighPt_V15","V15",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV15.root","output/ElectronNtuple.Real.Subdet0HighPtV16.root","Subdet0HighPt_V16","V16",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV16.root","output/ElectronNtuple.Real.Subdet0HighPtV17.root","Subdet0HighPt_V17","V17",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtV17.root","output/ElectronNtuple.Real.Subdet0HighPtV18.root","Subdet0HighPt_V18","V18",kFALSE,kFALSE,3,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","output/ElectronNtuple.Real.Subdet1HighPtV0.root","Subdet1HighPt_V0","V0",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV0.root","output/ElectronNtuple.Real.Subdet1HighPtV1.root","Subdet1HighPt_V1","V1",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV1.root","output/ElectronNtuple.Real.Subdet1HighPtV2.root","Subdet1HighPt_V2","V2",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV2.root","output/ElectronNtuple.Real.Subdet1HighPtV3.root","Subdet1HighPt_V3","V3",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV3.root","output/ElectronNtuple.Real.Subdet1HighPtV4.root","Subdet1HighPt_V4","V4",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV4.root","output/ElectronNtuple.Real.Subdet1HighPtV5.root","Subdet1HighPt_V5","V5",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV5.root","output/ElectronNtuple.Real.Subdet1HighPtV6.root","Subdet1HighPt_V6","V6",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV6.root","output/ElectronNtuple.Real.Subdet1HighPtV7.root","Subdet1HighPt_V7","V7",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV7.root","output/ElectronNtuple.Real.Subdet1HighPtV8.root","Subdet1HighPt_V8","V8",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV8.root","output/ElectronNtuple.Real.Subdet1HighPtV9.root","Subdet1HighPt_V9","V9",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV9.root","output/ElectronNtuple.Real.Subdet1HighPtV10.root","Subdet1HighPt_V10","V10",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV10.root","output/ElectronNtuple.Real.Subdet1HighPtV11.root","Subdet1HighPt_V11","V11",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV11.root","output/ElectronNtuple.Real.Subdet1HighPtV12.root","Subdet1HighPt_V12","V12",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV12.root","output/ElectronNtuple.Real.Subdet1HighPtV13.root","Subdet1HighPt_V13","V13",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV13.root","output/ElectronNtuple.Real.Subdet1HighPtV14.root","Subdet1HighPt_V14","V14",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV14.root","output/ElectronNtuple.Real.Subdet1HighPtV15.root","Subdet1HighPt_V15","V15",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV15.root","output/ElectronNtuple.Real.Subdet1HighPtV16.root","Subdet1HighPt_V16","V16",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV16.root","output/ElectronNtuple.Real.Subdet1HighPtV17.root","Subdet1HighPt_V17","V17",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtV17.root","output/ElectronNtuple.Real.Subdet1HighPtV18.root","Subdet1HighPt_V18","V18",kFALSE,kFALSE,4,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","output/ElectronNtuple.Real.Subdet2HighPtV0.root","Subdet2HighPt_V0","V0",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV0.root","output/ElectronNtuple.Real.Subdet2HighPtV1.root","Subdet2HighPt_V1","V1",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV1.root","output/ElectronNtuple.Real.Subdet2HighPtV2.root","Subdet2HighPt_V2","V2",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV2.root","output/ElectronNtuple.Real.Subdet2HighPtV3.root","Subdet2HighPt_V3","V3",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV3.root","output/ElectronNtuple.Real.Subdet2HighPtV4.root","Subdet2HighPt_V4","V4",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV4.root","output/ElectronNtuple.Real.Subdet2HighPtV5.root","Subdet2HighPt_V5","V5",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV5.root","output/ElectronNtuple.Real.Subdet2HighPtV6.root","Subdet2HighPt_V6","V6",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV6.root","output/ElectronNtuple.Real.Subdet2HighPtV7.root","Subdet2HighPt_V7","V7",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV7.root","output/ElectronNtuple.Real.Subdet2HighPtV8.root","Subdet2HighPt_V8","V8",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV8.root","output/ElectronNtuple.Real.Subdet2HighPtV9.root","Subdet2HighPt_V9","V9",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV9.root","output/ElectronNtuple.Real.Subdet2HighPtV10.root","Subdet2HighPt_V10","V10",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV10.root","output/ElectronNtuple.Real.Subdet2HighPtV11.root","Subdet2HighPt_V11","V11",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV11.root","output/ElectronNtuple.Real.Subdet2HighPtV12.root","Subdet2HighPt_V12","V12",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV12.root","output/ElectronNtuple.Real.Subdet2HighPtV13.root","Subdet2HighPt_V13","V13",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV13.root","output/ElectronNtuple.Real.Subdet2HighPtV14.root","Subdet2HighPt_V14","V14",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV14.root","output/ElectronNtuple.Real.Subdet2HighPtV15.root","Subdet2HighPt_V15","V15",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV15.root","output/ElectronNtuple.Real.Subdet2HighPtV16.root","Subdet2HighPt_V16","V16",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV16.root","output/ElectronNtuple.Real.Subdet2HighPtV17.root","Subdet2HighPt_V17","V17",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtV17.root","output/ElectronNtuple.Real.Subdet2HighPtV18.root","Subdet2HighPt_V18","V18",kFALSE,kFALSE,5,"BDTG")'






#HWW
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW115.root","output/ElectronNtuple.HWW115_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V0.root","output/ElectronNtuple.HWW115_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V1.root","output/ElectronNtuple.HWW115_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V2.root","output/ElectronNtuple.HWW115_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V3.root","output/ElectronNtuple.HWW115_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V4.root","output/ElectronNtuple.HWW115_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V5.root","output/ElectronNtuple.HWW115_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V6.root","output/ElectronNtuple.HWW115_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V7.root","output/ElectronNtuple.HWW115_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V8.root","output/ElectronNtuple.HWW115_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V9.root","output/ElectronNtuple.HWW115_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V10.root","output/ElectronNtuple.HWW115_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V11.root","output/ElectronNtuple.HWW115_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V12.root","output/ElectronNtuple.HWW115_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V13.root","output/ElectronNtuple.HWW115_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V14.root","output/ElectronNtuple.HWW115_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V15.root","output/ElectronNtuple.HWW115_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V16.root","output/ElectronNtuple.HWW115_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115_V17.root","output/ElectronNtuple.HWW115_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW120.root","output/ElectronNtuple.HWW120_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V0.root","output/ElectronNtuple.HWW120_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V1.root","output/ElectronNtuple.HWW120_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V2.root","output/ElectronNtuple.HWW120_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V3.root","output/ElectronNtuple.HWW120_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V4.root","output/ElectronNtuple.HWW120_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V5.root","output/ElectronNtuple.HWW120_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V6.root","output/ElectronNtuple.HWW120_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V7.root","output/ElectronNtuple.HWW120_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V8.root","output/ElectronNtuple.HWW120_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V9.root","output/ElectronNtuple.HWW120_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V10.root","output/ElectronNtuple.HWW120_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V11.root","output/ElectronNtuple.HWW120_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V12.root","output/ElectronNtuple.HWW120_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V13.root","output/ElectronNtuple.HWW120_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V14.root","output/ElectronNtuple.HWW120_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V15.root","output/ElectronNtuple.HWW116_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V16.root","output/ElectronNtuple.HWW120_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120_V17.root","output/ElectronNtuple.HWW120_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW130.root","output/ElectronNtuple.HWW130_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V0.root","output/ElectronNtuple.HWW130_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V1.root","output/ElectronNtuple.HWW130_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V2.root","output/ElectronNtuple.HWW130_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V3.root","output/ElectronNtuple.HWW130_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V4.root","output/ElectronNtuple.HWW130_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V5.root","output/ElectronNtuple.HWW130_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V6.root","output/ElectronNtuple.HWW130_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V7.root","output/ElectronNtuple.HWW130_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V8.root","output/ElectronNtuple.HWW130_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V9.root","output/ElectronNtuple.HWW130_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V10.root","output/ElectronNtuple.HWW130_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V11.root","output/ElectronNtuple.HWW130_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V12.root","output/ElectronNtuple.HWW130_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V13.root","output/ElectronNtuple.HWW130_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V14.root","output/ElectronNtuple.HWW130_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V15.root","output/ElectronNtuple.HWW116_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V16.root","output/ElectronNtuple.HWW130_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130_V17.root","output/ElectronNtuple.HWW130_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'



#Zee
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fall11Zee.root","output/ElectronNtuple.Fall11Zee_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V0.root","output/ElectronNtuple.Fall11Zee_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V1.root","output/ElectronNtuple.Fall11Zee_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V2.root","output/ElectronNtuple.Fall11Zee_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V3.root","output/ElectronNtuple.Fall11Zee_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V4.root","output/ElectronNtuple.Fall11Zee_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V5.root","output/ElectronNtuple.Fall11Zee_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V6.root","output/ElectronNtuple.Fall11Zee_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V7.root","output/ElectronNtuple.Fall11Zee_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V8.root","output/ElectronNtuple.Fall11Zee_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V9.root","output/ElectronNtuple.Fall11Zee_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V10.root","output/ElectronNtuple.Fall11Zee_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V11.root","output/ElectronNtuple.Fall11Zee_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V12.root","output/ElectronNtuple.Fall11Zee_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V13.root","output/ElectronNtuple.Fall11Zee_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V14.root","output/ElectronNtuple.Fall11Zee_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V15.root","output/ElectronNtuple.Fall11Zee_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V16.root","output/ElectronNtuple.Fall11Zee_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fall11Zee_V17.root","output/ElectronNtuple.Fall11Zee_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'



#WJets
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.WJets_Winter10Skimmed.root","output/ElectronNtuple.WJets_Winter10Skimmed_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V0.root","output/ElectronNtuple.WJets_Winter10Skimmed_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V1.root","output/ElectronNtuple.WJets_Winter10Skimmed_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V2.root","output/ElectronNtuple.WJets_Winter10Skimmed_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V3.root","output/ElectronNtuple.WJets_Winter10Skimmed_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V4.root","output/ElectronNtuple.WJets_Winter10Skimmed_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V5.root","output/ElectronNtuple.WJets_Winter10Skimmed_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V6.root","output/ElectronNtuple.WJets_Winter10Skimmed_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V7.root","output/ElectronNtuple.WJets_Winter10Skimmed_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V8.root","output/ElectronNtuple.WJets_Winter10Skimmed_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V9.root","output/ElectronNtuple.WJets_Winter10Skimmed_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V10.root","output/ElectronNtuple.WJets_Winter10Skimmed_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V11.root","output/ElectronNtuple.WJets_Winter10Skimmed_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V12.root","output/ElectronNtuple.WJets_Winter10Skimmed_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V13.root","output/ElectronNtuple.WJets_Winter10Skimmed_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V14.root","output/ElectronNtuple.WJets_Winter10Skimmed_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V15.root","output/ElectronNtuple.WJets_Winter10Skimmed_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V16.root","output/ElectronNtuple.WJets_Winter10Skimmed_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Winter10Skimmed_V17.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.WJets_Summer11Skimmed.root","output/ElectronNtuple.WJets_Summer11Skimmed_V0.root","V0","V0",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V0.root","output/ElectronNtuple.WJets_Summer11Skimmed_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V1.root","output/ElectronNtuple.WJets_Summer11Skimmed_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V2.root","output/ElectronNtuple.WJets_Summer11Skimmed_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V3.root","output/ElectronNtuple.WJets_Summer11Skimmed_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V4.root","output/ElectronNtuple.WJets_Summer11Skimmed_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V5.root","output/ElectronNtuple.WJets_Summer11Skimmed_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V6.root","output/ElectronNtuple.WJets_Summer11Skimmed_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V7.root","output/ElectronNtuple.WJets_Summer11Skimmed_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V8.root","output/ElectronNtuple.WJets_Summer11Skimmed_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V9.root","output/ElectronNtuple.WJets_Summer11Skimmed_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V10.root","output/ElectronNtuple.WJets_Summer11Skimmed_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V11.root","output/ElectronNtuple.WJets_Summer11Skimmed_V12.root","V12","V12",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V12.root","output/ElectronNtuple.WJets_Summer11Skimmed_V13.root","V13","V13",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V13.root","output/ElectronNtuple.WJets_Summer11Skimmed_V14.root","V14","V14",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V14.root","output/ElectronNtuple.WJets_Summer11Skimmed_V15.root","V15","V15",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V15.root","output/ElectronNtuple.WJets_Summer11Skimmed_V16.root","V16","V16",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V16.root","output/ElectronNtuple.WJets_Summer11Skimmed_V17.root","V17","V17",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.WJets_Summer11Skimmed_V17.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","V18","V18",kFALSE,kTRUE,-1,"BDTG")'


#####################################################################################
#Performance Plots
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet0LowPtV18.root","output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet1LowPtV18.root","output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet2LowPtV18.root","output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet0HighPtV18.root","output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet1HighPtV18.root","output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Real.Subdet2HighPtV18.root","output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt",5)'



#####################################################################################
#Performance Plots with HWW115
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt_HWW115",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt_HWW115",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt_HWW115",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt_HWW115",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt_HWW115",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt_HWW115",5)'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt_HWW130",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt_HWW130",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt_HWW130",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt_HWW130",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt_HWW130",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW130_V18.root","output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt_HWW130",5)'


#####################################################################################
#Performance Plots with Fall11Zee
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt_Fall11Zee",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt_Fall11Zee",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt_Fall11Zee",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt_Fall11Zee",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt_Fall11Zee",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.Fall11Zee_V18.root","output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt_Fall11Zee",5)'


#####################################################################################
#Performance Plots HWW115 Vs WJets MC
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet0LowPt_HWW115Summer11WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet1LowPt_HWW115Summer11WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet2LowPt_HWW115Summer11WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet0HighPt_HWW115Summer11WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet1HighPt_HWW115Summer11WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Summer11Skimmed_V18.root","Subdet2HighPt_HWW115Summer11WJets",5)'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet0LowPt_HWW115Winter10WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet1LowPt_HWW115Winter10WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet2LowPt_HWW115Winter10WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet0HighPt_HWW115Winter10WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet1HighPt_HWW115Winter10WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("output/ElectronNtuple.HWW115_V18.root","output/ElectronNtuple.WJets_Winter10Skimmed_V18.root","Subdet2HighPt_HWW115Winter10WJets",5)'



#####################################################################################
#Comparison Plots
#####################################################################################


root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0LowPt.root","output/ElectronNtuple.HWW115.Subdet0LowPt.root","Data","HWW115","Subdet0LowPt_Real",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1LowPt.root","output/ElectronNtuple.HWW115.Subdet1LowPt.root","Data","HWW115","Subdet1LowPt_Real",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2LowPt.root","output/ElectronNtuple.HWW115.Subdet2LowPt.root","Data","HWW115","Subdet2LowPt_Real",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0HighPt.root","output/ElectronNtuple.HWW115.Subdet0HighPt.root","Data","HWW115","Subdet0HighPt_Real",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1HighPt.root","output/ElectronNtuple.HWW115.Subdet1HighPt.root","Data","HWW115","Subdet1HighPt_Real",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2HighPt.root","output/ElectronNtuple.HWW115.Subdet2HighPt.root","Data","HWW115","Subdet2HighPt_Real",5)'


