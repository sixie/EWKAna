##################
# Met Plots
##################
cd /data/blue/sixie/releases/analysis/CMSSW_4_2_3_patch2/src; cmsenv;
cp /data/blue/sixie/Thesis/EventSelection/ThesisPlots.SignalCharacteristics.root ./
root -l EWKAna/Thesis/Selection/SelectionMETPlots.C+'(1)' 

##################
# Top Tagging Plots
##################
cd /data/blue/sixie/releases/analysis/CMSSW_4_2_3_patch2/src; cmsenv;
cp /data/blue/sixie/Thesis/EventSelection/ThesisPlots.SignalCharacteristics.root ./
root -l EWKAna/Thesis/Selection/SelectionTopTaggingPlots.C+'(1)'




##################
# Make TMVA Plots
##################
cd TMVATest

##################
# 0-Jet Bin Plots
##################
root -l TMVAGui.C+'("/data/blue/sixie/Thesis/HWWMVA/ntuples_42x_160train_0jets.root")'

#click on input variables
cd plots
cp variables_id_c1.eps HWWMVA_PtMax_0Jet.eps
cp variables_id_c2.eps HWWMVA_PtMin_0Jet.eps
cp variables_id_c3.eps HWWMVA_DPhi_0Jet.eps
cp variables_id_c4.eps HWWMVA_DR_0Jet.eps
cp variables_id_c5.eps HWWMVA_DileptonMass_0Jet.eps
cp variables_id_c6.eps HWWMVA_FinalState_0Jet.eps
cp variables_id_c7.eps HWWMVA_MTHiggs_0Jet.eps

#click on output button
mv overtrain_BDTG.eps HWWMVA_Output_0Jet.eps

root -l TMVAGui.C+'("/data/blue/sixie/Thesis/HWWMVA/ntuples_42x_160train_1jets.root")'

#click on input variables
cd plots
cp variables_id_c1.eps HWWMVA_PtMax_1Jet.eps
cp variables_id_c2.eps HWWMVA_PtMin_1Jet.eps
cp variables_id_c3.eps HWWMVA_DPhi_1Jet.eps
cp variables_id_c4.eps HWWMVA_DR_1Jet.eps
cp variables_id_c5.eps HWWMVA_DileptonMass_1Jet.eps
cp variables_id_c6.eps HWWMVA_FinalState_1Jet.eps
cp variables_id_c7.eps HWWMVA_MTHiggs_1Jet.eps
cp variables_id_c8.eps HWWMVA_DPhiDilepMET_1Jet.eps
cp variables_id_c9.eps HWWMVA_DPhiDilepJet_1Jet.eps

#click on output button
mv overtrain_BDTG.eps HWWMVA_Output_1Jet.eps
