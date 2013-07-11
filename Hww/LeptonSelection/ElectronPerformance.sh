#####################################################################################
#Make Electron Ntuples
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeRealElectronTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeFakeElectronTrainingNtuple.C+
root -l -b -q EWKAna/Hww/LeptonSelection/MakeHWWElectronTrainingNtuple.C+

#####################################################################################
#Do pt reweighting
#####################################################################################
root -l -b -q EWKAna/Hww/LeptonSelection/MakeLeptonPtSpectrum.C+



#####################################################################################
#Skim Ntuples, Split into bins
#####################################################################################

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0LowPt",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1LowPt",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPt",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0HighPt",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1HighPt",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2HighPt",5)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0LowPtNBrem0",10)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1LowPtNBrem0",11)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPtNBrem0",12)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0HighPtNBrem0",13)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1HighPtNBrem0",14)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2HighPtNBrem0",15)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0LowPtNBremGreaterThan0",20)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1LowPtNBremGreaterThan0",21)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPtNBremGreaterThan0",22)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0HighPtNBremGreaterThan0",23)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1HighPtNBremGreaterThan0",24)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2HighPtNBremGreaterThan0",25)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPtNBremGreaterThan0Pt10To15",221)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPtNBremGreaterThan0Pt15To20",222)'

#####################################################################################
#MVA Training variables to electron ntuple
#####################################################################################

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPt","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","Subdet1LowPt","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","Subdet2LowPt","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","Subdet0HighPt","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","Subdet1HighPt","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","Subdet2HighPt","","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPtNBrem0.root","Subdet0LowPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPtNBrem0.root","Subdet1LowPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBrem0.root","Subdet2LowPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPtNBrem0.root","Subdet0HighPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPtNBrem0.root","Subdet1HighPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPtNBrem0.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPtNBrem0.root","Subdet2HighPtNBrem0","","Likelihood,LikelihoodD,BDT,BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPtNBremGreaterThan0.root","Subdet0LowPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPtNBremGreaterThan0.root","Subdet1LowPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0.root","Subdet2LowPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPtNBremGreaterThan0.root","Subdet0HighPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPtNBremGreaterThan0.root","Subdet1HighPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPtNBremGreaterThan0.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPtNBremGreaterThan0.root","Subdet2HighPtNBremGreaterThan0","","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0Pt10To15.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0Pt10To15.root","Subdet2LowPtNBremGreaterThan0Pt10To15","","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0Pt15To20.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0Pt15To20.root","Subdet2LowPtNBremGreaterThan0Pt15To20","","Likelihood,LikelihoodD,BDT,BDTG")'


##Train including StandardLikelihood
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","Subdet1LowPtWithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","Subdet2LowPtWithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","Subdet0HighPtWithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","Subdet1HighPtWithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","Subdet2HighPtWithLH","","BDTG")'

##Train including StandardLikelihood, H/E, E/P, d0 (V2)
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","Subdet1LowPtWithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","Subdet2LowPtWithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","Subdet0HighPtWithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","Subdet1HighPtWithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","Subdet2HighPtWithLHV2","V2","BDTG")'

##Train including StandardLikelihood, H/E, E/P, d0 , ... (V3)
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","Subdet1LowPtWithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","Subdet2LowPtWithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","Subdet0HighPtWithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","Subdet1HighPtWithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","Subdet2HighPtWithLHV3","V3","BDTG")'


##Train including StandardLikelihood, H/E, E/P, d0 , ... changing BDTG training parameters to what josh suggested : no bagging, maxdepth=6 (V5)
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLHV4","V4","BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLHV5","V5","BDTG")'



##Train including StandardLikelihood, all variables used (V9)
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","Subdet0LowPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","Subdet1LowPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","Subdet2LowPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","Subdet0HighPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","Subdet1HighPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/TrainElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","Subdet2HighPtWithLHV9","V9","Likelihood,LikelihoodD,BDT,BDTG")'






#####################################################################################
#Fill MVA variables to electron ntuple
#####################################################################################

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","output/ElectronNtuple.Real.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPtNBrem0.root","output/ElectronNtuple.Fake.Subdet0LowPtNBrem0.root","Subdet0LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPtNBrem0.root","output/ElectronNtuple.Real.Subdet0LowPtNBrem0.root","Subdet0LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet0LowPtNBremGreaterThan0.root","Subdet0LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet0LowPtNBremGreaterThan0.root","Subdet0LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","output/ElectronNtuple.Real.Subdet1LowPt.root","Subdet1LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1LowPtNBrem0.root","output/ElectronNtuple.Fake.Subdet1LowPtNBrem0.root","Subdet1LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPtNBrem0.root","output/ElectronNtuple.Real.Subdet1LowPtNBrem0.root","Subdet1LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1LowPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet1LowPtNBremGreaterThan0.root","Subdet1LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet1LowPtNBremGreaterThan0.root","Subdet1LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","output/ElectronNtuple.Real.Subdet2LowPt.root","Subdet2LowPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBrem0.root","output/ElectronNtuple.Fake.Subdet2LowPtNBrem0.root","Subdet2LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBrem0.root","output/ElectronNtuple.Real.Subdet2LowPtNBrem0.root","Subdet2LowPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0.root","Subdet2LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0.root","Subdet2LowPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0Pt10To15.root","output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0Pt10To15.root","Subdet2LowPtNBremGreaterThan0Pt10To15","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0Pt10To15.root","output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0Pt10To15.root","Subdet2LowPtNBremGreaterThan0Pt10To15","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPtNBremGreaterThan0Pt15To20.root","output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0Pt15To20.root","Subdet2LowPtNBremGreaterThan0Pt15To20","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPtNBremGreaterThan0Pt15To20.root","output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0Pt15To20.root","Subdet2LowPtNBremGreaterThan0Pt15To20","Likelihood,LikelihoodD,BDT,BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","output/ElectronNtuple.Real.Subdet0HighPt.root","Subdet0HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0HighPtNBrem0.root","output/ElectronNtuple.Fake.Subdet0HighPtNBrem0.root","Subdet0HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPtNBrem0.root","output/ElectronNtuple.Real.Subdet0HighPtNBrem0.root","Subdet0HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0HighPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet0HighPtNBremGreaterThan0.root","Subdet0HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet0HighPtNBremGreaterThan0.root","Subdet0HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","output/ElectronNtuple.Real.Subdet1HighPt.root","Subdet1HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1HighPtNBrem0.root","output/ElectronNtuple.Fake.Subdet1HighPtNBrem0.root","Subdet1HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPtNBrem0.root","output/ElectronNtuple.Real.Subdet1HighPtNBrem0.root","Subdet1HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1HighPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet1HighPtNBremGreaterThan0.root","Subdet1HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet1HighPtNBremGreaterThan0.root","Subdet1HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","output/ElectronNtuple.Real.Subdet2HighPt.root","Subdet2HighPt","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2HighPtNBrem0.root","output/ElectronNtuple.Fake.Subdet2HighPtNBrem0.root","Subdet2HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPtNBrem0.root","output/ElectronNtuple.Real.Subdet2HighPtNBrem0.root","Subdet2HighPtNBrem0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2HighPtNBremGreaterThan0.root","output/ElectronNtuple.Fake.Subdet2HighPtNBremGreaterThan0.root","Subdet2HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPtNBremGreaterThan0.root","output/ElectronNtuple.Real.Subdet2HighPtNBremGreaterThan0.root","Subdet2HighPtNBremGreaterThan0","Likelihood,LikelihoodD,BDT,BDTG")'





######
#Trained With LH
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW115.root","output/ElectronNtuple.HWW115.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW120.root","output/ElectronNtuple.HWW120.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW130.root","output/ElectronNtuple.HWW130.root","WithLH","","BDTG")'



root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","output/ElectronNtuple.Real.Subdet0LowPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1LowPt.root","output/ElectronNtuple.Fake.Subdet1LowPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1LowPt.root","output/ElectronNtuple.Real.Subdet1LowPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2LowPt.root","output/ElectronNtuple.Fake.Subdet2LowPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2LowPt.root","output/ElectronNtuple.Real.Subdet2LowPtWithLH.root","WithLH","","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0HighPt.root","output/ElectronNtuple.Fake.Subdet0HighPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0HighPt.root","output/ElectronNtuple.Real.Subdet0HighPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet1HighPt.root","output/ElectronNtuple.Fake.Subdet1HighPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet1HighPt.root","output/ElectronNtuple.Real.Subdet1HighPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet2HighPt.root","output/ElectronNtuple.Fake.Subdet2HighPtWithLH.root","WithLH","","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet2HighPt.root","output/ElectronNtuple.Real.Subdet2HighPtWithLH.root","WithLH","","BDTG")'


root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW115.root","output/ElectronNtuple.HWW115.Subdet0LowPtWithLH.root","Subdet0LowPtWithLH","","BDTG")'


######
#Training V2
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115.root","output/ElectronNtuple.HWW115.V2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW120.root","output/ElectronNtuple.HWW120.V2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.HWW130.root","output/ElectronNtuple.HWW130.V2.root","WithLHV2","V2","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLH.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLH.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtWithLH.root","output/ElectronNtuple.Fake.Subdet1LowPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtWithLH.root","output/ElectronNtuple.Real.Subdet1LowPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtWithLH.root","output/ElectronNtuple.Fake.Subdet2LowPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtWithLH.root","output/ElectronNtuple.Real.Subdet2LowPtWithLHV2.root","WithLHV2","V2","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtWithLH.root","output/ElectronNtuple.Fake.Subdet0HighPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtWithLH.root","output/ElectronNtuple.Real.Subdet0HighPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtWithLH.root","output/ElectronNtuple.Fake.Subdet1HighPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtWithLH.root","output/ElectronNtuple.Real.Subdet1HighPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtWithLH.root","output/ElectronNtuple.Fake.Subdet2HighPtWithLHV2.root","WithLHV2","V2","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtWithLH.root","output/ElectronNtuple.Real.Subdet2HighPtWithLHV2.root","WithLHV2","V2","BDTG")'

######
#Training V3
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW115.V2.root","output/ElectronNtuple.HWW115.V3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW120.V2.root","output/ElectronNtuple.HWW120.V3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.HWW130.V2.root","output/ElectronNtuple.HWW130.V3.root","WithLHV3","V3","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV2.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1LowPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet1LowPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1LowPtWithLHV2.root","output/ElectronNtuple.Real.Subdet1LowPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2LowPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet2LowPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2LowPtWithLHV2.root","output/ElectronNtuple.Real.Subdet2LowPtWithLHV3.root","WithLHV3","V3","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0HighPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet0HighPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0HighPtWithLHV2.root","output/ElectronNtuple.Real.Subdet0HighPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet1HighPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet1HighPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet1HighPtWithLHV2.root","output/ElectronNtuple.Real.Subdet1HighPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet2HighPtWithLHV2.root","output/ElectronNtuple.Fake.Subdet2HighPtWithLHV3.root","WithLHV3","V3","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet2HighPtWithLHV2.root","output/ElectronNtuple.Real.Subdet2HighPtWithLHV3.root","WithLHV3","V3","BDTG")'


######
#Training V4

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV3.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV5.root","Subdet0LowPtWithLHV5","V5","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV3.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV5.root","Subdet0LowPtWithLHV5","V5","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV5.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV4.root","Subdet0LowPtWithLHV4","V4","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV5.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV4.root","Subdet0LowPtWithLHV4","V4","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV5.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV6.root","Subdet0LowPtWithLHV6","V6","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV5.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV6.root","Subdet0LowPtWithLHV6","V6","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV6.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV7.root","Subdet0LowPtWithLHV7","V7","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV6.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV7.root","Subdet0LowPtWithLHV7","V7","BDTG")'

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV7.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV8.root","Subdet0LowPtWithLHV8","V8","BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV7.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV8.root","Subdet0LowPtWithLHV8","V8","BDTG")'


######
#Trained With all variables (V9)

root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Fake.Subdet0LowPtWithLHV3.root","output/ElectronNtuple.Fake.Subdet0LowPtWithLHV9.root","Subdet0LowPtWithLHV9","V9","BDT,BDTG")'
root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("output/ElectronNtuple.Real.Subdet0LowPtWithLHV3.root","output/ElectronNtuple.Real.Subdet0LowPtWithLHV9.root","Subdet0LowPtWithLHV9","V9","BDT,BDTG")'






#####################################################################################
#Make MVA Performance Plots 
#####################################################################################

root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPtNBrem0",10)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPtNBremGreaterThan0",20)'


root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt",1)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPtNBrem0",11)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPtNBremGreaterThan0",21)'


root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt",2)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPtNBrem0",12)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPtNBremGreaterThan0",22)'


root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt",3)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPtNBrem0",13)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPtNBremGreaterThan0",23)'


root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt",4)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPtNBrem0",14)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPtNBremGreaterThan0",24)'


root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt",5)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPtNBrem0",15)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPtNBremGreaterThan0",25)'


#####################################################################################
#Make MVA Performance Plots for training with NBrem binning
#####################################################################################

root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPtNBrem0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPtNBrem0.root","Subdet0LowPtNBrem0",10)'
root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPtNBremGreaterThan0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPtNBremGreaterThan0.root","Subdet0LowPtNBremGreaterThan0",20)'

root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPtNBrem0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPtNBrem0.root","Subdet1LowPtNBrem0",11)'
root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPtNBremGreaterThan0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPtNBremGreaterThan0.root","Subdet1LowPtNBremGreaterThan0",21)'

root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtNBrem0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtNBrem0.root","Subdet2LowPtNBrem0",12)'
root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0.root","Subdet2LowPtNBremGreaterThan0",22)'


#####################################################################################
#Make MVA Performance Plots for training with StandardLikelihood 
#####################################################################################
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPtWithLH.root","Subdet0LowPtWithLH",0)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPtWithLH.root","Subdet1LowPtWithLH",1)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtWithLH.root","Subdet2LowPtWithLH",2)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPtWithLH.root","Subdet0HighPtWithLH",3)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPtWithLH.root","Subdet1HighPtWithLH",4)'
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPtWithLH.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPtWithLH.root","Subdet2HighPtWithLH",5)'


#####################################################################################
#Make MVA Performance Plots for all trainings
#####################################################################################
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPtWithLHV3.root","Subdet0LowPtWithLHAllOptions",0)'
root -l -b -q   /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPtWithLHV3.root","Subdet1LowPtWithLHAllOptions",1)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtWithLHV3.root","Subdet2LowPtWithLHAllOptions",2)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPtWithLHV3.root","Subdet0HighPtWithLHAllOptions",3)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPtWithLHV3.root","Subdet1HighPtWithLHAllOptions",4)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPtWithLHV3.root","Subdet2HighPtWithLHAllOptions",5)'

#####################################################################################
#Make MVA Performance Plots for all trainings using MC signal
#####################################################################################
root -l -b -q  /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPtWithLHV3.root","Subdet0LowPtWithLHAllOptions_HWW115Signal",0)'
root -l -b -q   /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPtWithLHV3.root","Subdet1LowPtWithLHAllOptions_HWW115Signal",1)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtWithLHV3.root","Subdet2LowPtWithLHAllOptions_HWW115Signal",2)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPtWithLHV3.root","Subdet0HighPtWithLHAllOptions_HWW115Signal",3)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPtWithLHV3.root","Subdet1HighPtWithLHAllOptions_HWW115Signal",4)'
root -l -b -q    /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPtWithLHV3.root","Subdet2HighPtWithLHAllOptions_HWW115Signal",5)'





#####################################################################################
#Evaluate Likelihood Performance 
#####################################################################################

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt_HWW120MC",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt_DataReweightedToHWW120",0)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt_HWW120MC",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt_DataReweightedToHWW120",1)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt_HWW120MC",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt_DataReweightedToHWW120",2)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt_HWW120MC",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt_DataReweightedToHWW120",3)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt_HWW120MC",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt_DataReweightedToHWW120",4)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW120.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt_HWW120MC",5)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt_DataReweightedToHWW120",5)'





root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt_HWW115MC",0)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt_DataReweightedToHWW115",0)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt_HWW115MC",1)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1LowPt.root","Subdet1LowPt_DataReweightedToHWW115",1)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt_HWW115MC",2)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPt.root","Subdet2LowPt_DataReweightedToHWW115",2)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt_HWW115MC",3)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0HighPt.root","Subdet0HighPt_DataReweightedToHWW115",3)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt_HWW115MC",4)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet1HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet1HighPt.root","Subdet1HighPt_DataReweightedToHWW115",4)'

root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt_HWW115MC",5)'
root -l -b -q /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2HighPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2HighPt.root","Subdet2HighPt_DataReweightedToHWW115",5)'




#####################################################################################
#Make Distribution Plots
#####################################################################################

root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0.root","Subdet2LowPtNBremGreaterThan0",22)'

root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0Pt10To15.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0Pt10To15.root","Subdet2LowPtNBremGreaterThan0Pt10To15",22)'
root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet2LowPtNBremGreaterThan0Pt15To20.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet2LowPtNBremGreaterThan0Pt15To20.root","Subdet2LowPtNBremGreaterThan0Pt15To20",22)'



root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0HighPtWithLHV3.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.HWW115.V3.root","Data","HWW115","Subdet0HighPt",3)'
