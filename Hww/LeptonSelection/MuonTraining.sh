#####################################################################################
#Make Muon Ntuples
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeRealMuonTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeFakeMuonTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeMCMuonTrainingNtuple.C+

#####################################################################################
#Do pt reweighting
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeMuonPtSpectrum.C+
root -l -b -q EWKAna/Hww/LeptonSelection/ReweightMuonPU.C+



#####################################################################################
#Skim Ntuples, Split into bins
#####################################################################################

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("BarrelPtBin0",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("EndcapPtBin0",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("BarrelPtBin1",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("EndcapPtBin1",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("BarrelPtBin2",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimMuonNtuples.C+'("EndcapPtBin2",5)'



#####################################################################################
#MVA Training 
#####################################################################################

#Add Variables
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V0","V0",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V0","V0",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V0","V0",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V0","V0",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V0","V0",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V0","V0",kFALSE,5,"Likelihood,BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V1","V1",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V1","V1",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V1","V1",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V1","V1",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V1","V1",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V1","V1",kFALSE,5,"Likelihood,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V2","V2",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V2","V2",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V2","V2",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V2","V2",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V2","V2",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V2","V2",kFALSE,5,"Likelihood,BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V3","V3",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V3","V3",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V3","V3",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V3","V3",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V3","V3",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V3","V3",kFALSE,5,"Likelihood,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V4","V4",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V4","V4",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V4","V4",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V4","V4",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V4","V4",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V4","V4",kFALSE,5,"Likelihood,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V5","V5",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V5","V5",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V5","V5",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V5","V5",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V5","V5",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V5","V5",kFALSE,5,"Likelihood,BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V6","V6",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V6","V6",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V6","V6",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V6","V6",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V6","V6",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V6","V6",kFALSE,5,"Likelihood,BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V7","V7",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V7","V7",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V7","V7",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V7","V7",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V7","V7",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V7","V7",kFALSE,5,"Likelihood,BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V8","V8",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V8","V8",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V8","V8",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V8","V8",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V8","V8",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V8","V8",kFALSE,5,"Likelihood,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V9","V9",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V9","V9",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V9","V9",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V9","V9",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V9","V9",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V9","V9",kFALSE,5,"Likelihood,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V10","V10",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V10","V10",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V10","V10",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V10","V10",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V10","V10",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V10","V10",kFALSE,5,"Likelihood,BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","BarrelPtBin0_V11","V11",kFALSE,0,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","EndcapPtBin0_V11","V11",kFALSE,1,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","BarrelPtBin1_V11","V11",kFALSE,2,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","EndcapPtBin1_V11","V11",kFALSE,3,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","BarrelPtBin2_V11","V11",kFALSE,4,"Likelihood,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","EndcapPtBin2_V11","V11",kFALSE,5,"Likelihood,BDTG")'




#####################################################################################
#MVA Evaluate
#####################################################################################

#FAKES
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin0.root","output/MuonNtuple.Fake.BarrelPtBin0_V1.root","BarrelPtBin0_V1","V1",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V1.root","output/MuonNtuple.Fake.BarrelPtBin0_V2.root","BarrelPtBin0_V2","V2",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V2.root","output/MuonNtuple.Fake.BarrelPtBin0_V3.root","BarrelPtBin0_V3","V3",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V3.root","output/MuonNtuple.Fake.BarrelPtBin0_V4.root","BarrelPtBin0_V4","V4",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V4.root","output/MuonNtuple.Fake.BarrelPtBin0_V5.root","BarrelPtBin0_V5","V5",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V5.root","output/MuonNtuple.Fake.BarrelPtBin0_V6.root","BarrelPtBin0_V6","V6",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V6.root","output/MuonNtuple.Fake.BarrelPtBin0_V7.root","BarrelPtBin0_V7","V7",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V7.root","output/MuonNtuple.Fake.BarrelPtBin0_V8.root","BarrelPtBin0_V8","V8",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V8.root","output/MuonNtuple.Fake.BarrelPtBin0_V9.root","BarrelPtBin0_V9","V9",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V9.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0_V10","V10",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin0_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V11.root","BarrelPtBin0_V11","V11",kFALSE,kFALSE,0,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin0.root","output/MuonNtuple.Fake.EndcapPtBin0_V1.root","EndcapPtBin0_V1","V1",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V1.root","output/MuonNtuple.Fake.EndcapPtBin0_V2.root","EndcapPtBin0_V2","V2",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V2.root","output/MuonNtuple.Fake.EndcapPtBin0_V3.root","EndcapPtBin0_V3","V3",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V3.root","output/MuonNtuple.Fake.EndcapPtBin0_V4.root","EndcapPtBin0_V4","V4",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V4.root","output/MuonNtuple.Fake.EndcapPtBin0_V5.root","EndcapPtBin0_V5","V5",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V5.root","output/MuonNtuple.Fake.EndcapPtBin0_V6.root","EndcapPtBin0_V6","V6",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V6.root","output/MuonNtuple.Fake.EndcapPtBin0_V7.root","EndcapPtBin0_V7","V7",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V7.root","output/MuonNtuple.Fake.EndcapPtBin0_V8.root","EndcapPtBin0_V8","V8",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V8.root","output/MuonNtuple.Fake.EndcapPtBin0_V9.root","EndcapPtBin0_V9","V9",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin0_V9.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_V10","V10",kFALSE,kFALSE,1,"BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin1.root","output/MuonNtuple.Fake.BarrelPtBin1_V1.root","BarrelPtBin1_V1","V1",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V1.root","output/MuonNtuple.Fake.BarrelPtBin1_V2.root","BarrelPtBin1_V2","V2",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V2.root","output/MuonNtuple.Fake.BarrelPtBin1_V3.root","BarrelPtBin1_V3","V3",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V3.root","output/MuonNtuple.Fake.BarrelPtBin1_V4.root","BarrelPtBin1_V4","V4",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V4.root","output/MuonNtuple.Fake.BarrelPtBin1_V5.root","BarrelPtBin1_V5","V5",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V5.root","output/MuonNtuple.Fake.BarrelPtBin1_V6.root","BarrelPtBin1_V6","V6",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V6.root","output/MuonNtuple.Fake.BarrelPtBin1_V7.root","BarrelPtBin1_V7","V7",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V7.root","output/MuonNtuple.Fake.BarrelPtBin1_V8.root","BarrelPtBin1_V8","V8",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V8.root","output/MuonNtuple.Fake.BarrelPtBin1_V9.root","BarrelPtBin1_V9","V9",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin1_V9.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_V10","V10",kFALSE,kFALSE,2,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin1.root","output/MuonNtuple.Fake.EndcapPtBin1_V1.root","EndcapPtBin1_V1","V1",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V1.root","output/MuonNtuple.Fake.EndcapPtBin1_V2.root","EndcapPtBin1_V2","V2",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V2.root","output/MuonNtuple.Fake.EndcapPtBin1_V3.root","EndcapPtBin1_V3","V3",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V3.root","output/MuonNtuple.Fake.EndcapPtBin1_V4.root","EndcapPtBin1_V4","V4",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V4.root","output/MuonNtuple.Fake.EndcapPtBin1_V5.root","EndcapPtBin1_V5","V5",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V5.root","output/MuonNtuple.Fake.EndcapPtBin1_V6.root","EndcapPtBin1_V6","V6",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V6.root","output/MuonNtuple.Fake.EndcapPtBin1_V7.root","EndcapPtBin1_V7","V7",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V7.root","output/MuonNtuple.Fake.EndcapPtBin1_V8.root","EndcapPtBin1_V8","V8",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V8.root","output/MuonNtuple.Fake.EndcapPtBin1_V9.root","EndcapPtBin1_V9","V9",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin1_V9.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_V10","V10",kFALSE,kFALSE,3,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.BarrelPtBin2.root","output/MuonNtuple.Fake.BarrelPtBin2_V1.root","BarrelPtBin2_V1","V1",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V1.root","output/MuonNtuple.Fake.BarrelPtBin2_V2.root","BarrelPtBin2_V2","V2",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V2.root","output/MuonNtuple.Fake.BarrelPtBin2_V3.root","BarrelPtBin2_V3","V3",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V3.root","output/MuonNtuple.Fake.BarrelPtBin2_V4.root","BarrelPtBin2_V4","V4",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V4.root","output/MuonNtuple.Fake.BarrelPtBin2_V5.root","BarrelPtBin2_V5","V5",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V5.root","output/MuonNtuple.Fake.BarrelPtBin2_V6.root","BarrelPtBin2_V6","V6",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V6.root","output/MuonNtuple.Fake.BarrelPtBin2_V7.root","BarrelPtBin2_V7","V7",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V7.root","output/MuonNtuple.Fake.BarrelPtBin2_V8.root","BarrelPtBin2_V8","V8",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V8.root","output/MuonNtuple.Fake.BarrelPtBin2_V9.root","BarrelPtBin2_V9","V9",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.BarrelPtBin2_V9.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_V10","V10",kFALSE,kFALSE,4,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fake.PtAndPUWeighted.EndcapPtBin2.root","output/MuonNtuple.Fake.EndcapPtBin2_V1.root","EndcapPtBin2_V1","V1",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V1.root","output/MuonNtuple.Fake.EndcapPtBin2_V2.root","EndcapPtBin2_V2","V2",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V2.root","output/MuonNtuple.Fake.EndcapPtBin2_V3.root","EndcapPtBin2_V3","V3",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V3.root","output/MuonNtuple.Fake.EndcapPtBin2_V4.root","EndcapPtBin2_V4","V4",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V4.root","output/MuonNtuple.Fake.EndcapPtBin2_V5.root","EndcapPtBin2_V5","V5",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V5.root","output/MuonNtuple.Fake.EndcapPtBin2_V6.root","EndcapPtBin2_V6","V6",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V6.root","output/MuonNtuple.Fake.EndcapPtBin2_V7.root","EndcapPtBin2_V7","V7",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V7.root","output/MuonNtuple.Fake.EndcapPtBin2_V8.root","EndcapPtBin2_V8","V8",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V8.root","output/MuonNtuple.Fake.EndcapPtBin2_V9.root","EndcapPtBin2_V9","V9",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fake.EndcapPtBin2_V9.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_V10","V10",kFALSE,kFALSE,5,"BDTG")'



#Real
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin0.root","output/MuonNtuple.Real.BarrelPtBin0_V1.root","BarrelPtBin0_V1","V1",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V1.root","output/MuonNtuple.Real.BarrelPtBin0_V2.root","BarrelPtBin0_V2","V2",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V2.root","output/MuonNtuple.Real.BarrelPtBin0_V3.root","BarrelPtBin0_V3","V3",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V3.root","output/MuonNtuple.Real.BarrelPtBin0_V4.root","BarrelPtBin0_V4","V4",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V4.root","output/MuonNtuple.Real.BarrelPtBin0_V5.root","BarrelPtBin0_V5","V5",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V5.root","output/MuonNtuple.Real.BarrelPtBin0_V6.root","BarrelPtBin0_V6","V6",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V6.root","output/MuonNtuple.Real.BarrelPtBin0_V7.root","BarrelPtBin0_V7","V7",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V7.root","output/MuonNtuple.Real.BarrelPtBin0_V8.root","BarrelPtBin0_V8","V8",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V8.root","output/MuonNtuple.Real.BarrelPtBin0_V9.root","BarrelPtBin0_V9","V9",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V9.root","output/MuonNtuple.Real.BarrelPtBin0_V10.root","BarrelPtBin0_V10","V10",kFALSE,kFALSE,0,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin0_V10.root","output/MuonNtuple.Real.BarrelPtBin0_V11.root","BarrelPtBin0_V11","V11",kFALSE,kFALSE,0,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin0.root","output/MuonNtuple.Real.EndcapPtBin0_V1.root","EndcapPtBin0_V1","V1",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V1.root","output/MuonNtuple.Real.EndcapPtBin0_V2.root","EndcapPtBin0_V2","V2",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V2.root","output/MuonNtuple.Real.EndcapPtBin0_V3.root","EndcapPtBin0_V3","V3",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V3.root","output/MuonNtuple.Real.EndcapPtBin0_V4.root","EndcapPtBin0_V4","V4",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V4.root","output/MuonNtuple.Real.EndcapPtBin0_V5.root","EndcapPtBin0_V5","V5",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V5.root","output/MuonNtuple.Real.EndcapPtBin0_V6.root","EndcapPtBin0_V6","V6",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V6.root","output/MuonNtuple.Real.EndcapPtBin0_V7.root","EndcapPtBin0_V7","V7",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V7.root","output/MuonNtuple.Real.EndcapPtBin0_V8.root","EndcapPtBin0_V8","V8",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V8.root","output/MuonNtuple.Real.EndcapPtBin0_V9.root","EndcapPtBin0_V9","V9",kFALSE,kFALSE,1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin0_V9.root","output/MuonNtuple.Real.EndcapPtBin0_V10.root","EndcapPtBin0_V10","V10",kFALSE,kFALSE,1,"BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin1.root","output/MuonNtuple.Real.BarrelPtBin1_V1.root","BarrelPtBin1_V1","V1",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V1.root","output/MuonNtuple.Real.BarrelPtBin1_V2.root","BarrelPtBin1_V2","V2",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V2.root","output/MuonNtuple.Real.BarrelPtBin1_V3.root","BarrelPtBin1_V3","V3",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V3.root","output/MuonNtuple.Real.BarrelPtBin1_V4.root","BarrelPtBin1_V4","V4",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V4.root","output/MuonNtuple.Real.BarrelPtBin1_V5.root","BarrelPtBin1_V5","V5",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V5.root","output/MuonNtuple.Real.BarrelPtBin1_V6.root","BarrelPtBin1_V6","V6",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V6.root","output/MuonNtuple.Real.BarrelPtBin1_V7.root","BarrelPtBin1_V7","V7",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V7.root","output/MuonNtuple.Real.BarrelPtBin1_V8.root","BarrelPtBin1_V8","V8",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V8.root","output/MuonNtuple.Real.BarrelPtBin1_V9.root","BarrelPtBin1_V9","V9",kFALSE,kFALSE,2,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin1_V9.root","output/MuonNtuple.Real.BarrelPtBin1_V10.root","BarrelPtBin1_V10","V10",kFALSE,kFALSE,2,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin1.root","output/MuonNtuple.Real.EndcapPtBin1_V1.root","EndcapPtBin1_V1","V1",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V1.root","output/MuonNtuple.Real.EndcapPtBin1_V2.root","EndcapPtBin1_V2","V2",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V2.root","output/MuonNtuple.Real.EndcapPtBin1_V3.root","EndcapPtBin1_V3","V3",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V3.root","output/MuonNtuple.Real.EndcapPtBin1_V4.root","EndcapPtBin1_V4","V4",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V4.root","output/MuonNtuple.Real.EndcapPtBin1_V5.root","EndcapPtBin1_V5","V5",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V5.root","output/MuonNtuple.Real.EndcapPtBin1_V6.root","EndcapPtBin1_V6","V6",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V6.root","output/MuonNtuple.Real.EndcapPtBin1_V7.root","EndcapPtBin1_V7","V7",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V7.root","output/MuonNtuple.Real.EndcapPtBin1_V8.root","EndcapPtBin1_V8","V8",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V8.root","output/MuonNtuple.Real.EndcapPtBin1_V9.root","EndcapPtBin1_V9","V9",kFALSE,kFALSE,3,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin1_V9.root","output/MuonNtuple.Real.EndcapPtBin1_V10.root","EndcapPtBin1_V10","V10",kFALSE,kFALSE,3,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.BarrelPtBin2.root","output/MuonNtuple.Real.BarrelPtBin2_V1.root","BarrelPtBin2_V1","V1",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V1.root","output/MuonNtuple.Real.BarrelPtBin2_V2.root","BarrelPtBin2_V2","V2",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V2.root","output/MuonNtuple.Real.BarrelPtBin2_V3.root","BarrelPtBin2_V3","V3",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V3.root","output/MuonNtuple.Real.BarrelPtBin2_V4.root","BarrelPtBin2_V4","V4",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V4.root","output/MuonNtuple.Real.BarrelPtBin2_V5.root","BarrelPtBin2_V5","V5",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V5.root","output/MuonNtuple.Real.BarrelPtBin2_V6.root","BarrelPtBin2_V6","V6",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V6.root","output/MuonNtuple.Real.BarrelPtBin2_V7.root","BarrelPtBin2_V7","V7",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V7.root","output/MuonNtuple.Real.BarrelPtBin2_V8.root","BarrelPtBin2_V8","V8",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V8.root","output/MuonNtuple.Real.BarrelPtBin2_V9.root","BarrelPtBin2_V9","V9",kFALSE,kFALSE,4,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.BarrelPtBin2_V9.root","output/MuonNtuple.Real.BarrelPtBin2_V10.root","BarrelPtBin2_V10","V10",kFALSE,kFALSE,4,"BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Real.weighted.EndcapPtBin2.root","output/MuonNtuple.Real.EndcapPtBin2_V1.root","EndcapPtBin2_V1","V1",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V1.root","output/MuonNtuple.Real.EndcapPtBin2_V2.root","EndcapPtBin2_V2","V2",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V2.root","output/MuonNtuple.Real.EndcapPtBin2_V3.root","EndcapPtBin2_V3","V3",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V3.root","output/MuonNtuple.Real.EndcapPtBin2_V4.root","EndcapPtBin2_V4","V4",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V4.root","output/MuonNtuple.Real.EndcapPtBin2_V5.root","EndcapPtBin2_V5","V5",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V5.root","output/MuonNtuple.Real.EndcapPtBin2_V6.root","EndcapPtBin2_V6","V6",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V6.root","output/MuonNtuple.Real.EndcapPtBin2_V7.root","EndcapPtBin2_V7","V7",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V7.root","output/MuonNtuple.Real.EndcapPtBin2_V8.root","EndcapPtBin2_V8","V8",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V8.root","output/MuonNtuple.Real.EndcapPtBin2_V9.root","EndcapPtBin2_V9","V9",kFALSE,kFALSE,5,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Real.EndcapPtBin2_V9.root","output/MuonNtuple.Real.EndcapPtBin2_V10.root","EndcapPtBin2_V10","V10",kFALSE,kFALSE,5,"BDTG")'






#HWW

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.HWW115.root","output/MuonNtuple.HWW115_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V1.root","output/MuonNtuple.HWW115_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V2.root","output/MuonNtuple.HWW115_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V3.root","output/MuonNtuple.HWW115_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V4.root","output/MuonNtuple.HWW115_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V5.root","output/MuonNtuple.HWW115_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V6.root","output/MuonNtuple.HWW115_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V7.root","output/MuonNtuple.HWW115_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V8.root","output/MuonNtuple.HWW115_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V9.root","output/MuonNtuple.HWW115_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.HWW115_V11.root","V11","V11",kFALSE,kTRUE,-1,"BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.HWW120.root","output/MuonNtuple.HWW120_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V1.root","output/MuonNtuple.HWW120_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V2.root","output/MuonNtuple.HWW120_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V3.root","output/MuonNtuple.HWW120_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V4.root","output/MuonNtuple.HWW120_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V5.root","output/MuonNtuple.HWW120_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V6.root","output/MuonNtuple.HWW120_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V7.root","output/MuonNtuple.HWW120_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V8.root","output/MuonNtuple.HWW120_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW120_V9.root","output/MuonNtuple.HWW120_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.HWW130.root","output/MuonNtuple.HWW130_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V1.root","output/MuonNtuple.HWW130_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V2.root","output/MuonNtuple.HWW130_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V3.root","output/MuonNtuple.HWW130_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V4.root","output/MuonNtuple.HWW130_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V5.root","output/MuonNtuple.HWW130_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V6.root","output/MuonNtuple.HWW130_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V7.root","output/MuonNtuple.HWW130_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V8.root","output/MuonNtuple.HWW130_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.HWW130_V9.root","output/MuonNtuple.HWW130_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'



#Zmm
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.Fall11Zmm.root","output/MuonNtuple.Fall11Zmm_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V1.root","output/MuonNtuple.Fall11Zmm_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V2.root","output/MuonNtuple.Fall11Zmm_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V3.root","output/MuonNtuple.Fall11Zmm_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V4.root","output/MuonNtuple.Fall11Zmm_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V5.root","output/MuonNtuple.Fall11Zmm_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V6.root","output/MuonNtuple.Fall11Zmm_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V7.root","output/MuonNtuple.Fall11Zmm_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V8.root","output/MuonNtuple.Fall11Zmm_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.Fall11Zmm_V9.root","output/MuonNtuple.Fall11Zmm_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'


#WJets

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.WJets_Summer11Skimmed.root","output/MuonNtuple.WJets_Summer11Skimmed_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V1.root","output/MuonNtuple.WJets_Summer11Skimmed_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V2.root","output/MuonNtuple.WJets_Summer11Skimmed_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V3.root","output/MuonNtuple.WJets_Summer11Skimmed_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V4.root","output/MuonNtuple.WJets_Summer11Skimmed_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V5.root","output/MuonNtuple.WJets_Summer11Skimmed_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V6.root","output/MuonNtuple.WJets_Summer11Skimmed_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V7.root","output/MuonNtuple.WJets_Summer11Skimmed_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V8.root","output/MuonNtuple.WJets_Summer11Skimmed_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Summer11Skimmed_V9.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("MuonSelectionTraining.WJets_Winter10Skimmed.root","output/MuonNtuple.WJets_Winter10Skimmed_V1.root","V1","V1",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V1.root","output/MuonNtuple.WJets_Winter10Skimmed_V2.root","V2","V2",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V2.root","output/MuonNtuple.WJets_Winter10Skimmed_V3.root","V3","V3",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V3.root","output/MuonNtuple.WJets_Winter10Skimmed_V4.root","V4","V4",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V4.root","output/MuonNtuple.WJets_Winter10Skimmed_V5.root","V5","V5",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V5.root","output/MuonNtuple.WJets_Winter10Skimmed_V6.root","V6","V6",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V6.root","output/MuonNtuple.WJets_Winter10Skimmed_V7.root","V7","V7",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V7.root","output/MuonNtuple.WJets_Winter10Skimmed_V8.root","V8","V8",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V8.root","output/MuonNtuple.WJets_Winter10Skimmed_V9.root","V9","V9",kFALSE,kTRUE,-1,"BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateMuonMVA.C+'("output/MuonNtuple.WJets_Winter10Skimmed_V9.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","V10","V10",kFALSE,kTRUE,-1,"BDTG")'


#####################################################################################
#Performance Plots
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.BarrelPtBin0_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.EndcapPtBin0_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.BarrelPtBin1_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.EndcapPtBin1_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.BarrelPtBin2_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Real.EndcapPtBin2_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2",5)'




#####################################################################################
#Performance Plots with HWW115
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V11.root","output/MuonNtuple.Fake.BarrelPtBin0_V11.root","BarrelPtBin0_HWW115",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_HWW115",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_HWW115",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_HWW115",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_HWW115",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_HWW115",5)'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0_HWW130",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_HWW130",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_HWW130",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_HWW130",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_HWW130",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW130_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_HWW130",5)'


#####################################################################################
#Performance Plots with Fall11Zmm
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0_Fall11Zmm",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0_Fall11Zmm",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1_Fall11Zmm",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1_Fall11Zmm",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2_Fall11Zmm",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.Fall11Zmm_V10.root","output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2_Fall11Zmm",5)'



#####################################################################################
#Performance Plots HWW115 Vs WJets MC
#####################################################################################
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin0_HWW115Summer11WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin0_HWW115Summer11WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin1_HWW115Summer11WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin1_HWW115Summer11WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","BarrelPtBin2_HWW115Summer11WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Summer11Skimmed_V10.root","EndcapPtBin2_HWW115Summer11WJets",5)'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin0_HWW115Winter10WJets",0)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin0_HWW115Winter10WJets",1)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin1_HWW115Winter10WJets",2)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin1_HWW115Winter10WJets",3)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","BarrelPtBin2_HWW115Winter10WJets",4)'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("output/MuonNtuple.HWW115_V10.root","output/MuonNtuple.WJets_Winter10Skimmed_V10.root","EndcapPtBin2_HWW115Winter10WJets",5)'


#####################################################################################
#Comparison Plots
#####################################################################################


root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0LowPt.root","output/ElectronNtuple.HWW115.Subdet0LowPt.root","Data","HWW115","Subdet0LowPt_Real",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1LowPt.root","output/ElectronNtuple.HWW115.Subdet1LowPt.root","Data","HWW115","Subdet1LowPt_Real",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2LowPt.root","output/ElectronNtuple.HWW115.Subdet2LowPt.root","Data","HWW115","Subdet2LowPt_Real",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet0HighPt.root","output/ElectronNtuple.HWW115.Subdet0HighPt.root","Data","HWW115","Subdet0HighPt_Real",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet1HighPt.root","output/ElectronNtuple.HWW115.Subdet1HighPt.root","Data","HWW115","Subdet1HighPt_Real",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtuple.Real.Subdet2HighPt.root","output/ElectronNtuple.HWW115.Subdet2HighPt.root","Data","HWW115","Subdet2HighPt_Real",5)'



root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("output/ElectronNtupleC.Real.Subdet0LowPt.root","output/ElectronNtupleC.HWW115.Subdet0LowPt.root","Data","HWW115","Subdet0LowPt_C_Real",0)'
