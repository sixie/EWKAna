##################################
#2011 Analysis
#################################


#Make Shape Systematics Plots
root -l PlotShapeSystematics.C+


#Make Fake Lepton Bkg Prediction Table at WW Preselection Level
root -l -q SummarizeFakeLeptonBkgPrediction.C+'(0,"/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets_mva/ntuples_130train_0jets_hww_syst_skim3.root","/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets_mva/ntuples_130train_0jets_backgroundC_skim2.root",2)'

#Make Fake Lepton Bkg Prediction Table at WW Preselection Level
root -l -q MakeSameSignControlSampleTables.C+'(0,"/build/sixie/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/data_2l_skim2.root","/build/sixie/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/backgroundC_skim2.root",2)'

#Make DataDriven Vs MC Comparison
root -l MakeDataMCPredictionComparison.C+


##################################
#New MVAIdIsoCombined Analysis
#################################

#Produce Fake Rate Tex Tables
root -l PrintFakeRate.C+

#Produce Fake Rate Plots
root -l plotFakeRate.C+

#Produce Fake Bkg Tables
root -l -q SummarizeFakeLeptonBkgPrediction.C+'(0,"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP//mitf-alljets/backgroundC_skim2.root","/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/backgroundC_skim2.root",13)'

#Produce same sign sample background prediction tables
root -l -q MakeSameSignControlSampleTables.C+'(0,"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP//mitf-alljets/data_2l_skim2.root","/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/backgroundC_skim2.root",13)'

#Make Shape Systematics Plots
root -l PlotShapeSystematics.C+

#Make DataDriven Vs MC Comparison
root -l MakeDataMCPredictionComparison.C+
