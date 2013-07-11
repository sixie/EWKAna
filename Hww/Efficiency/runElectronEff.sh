#--------------------------------------------------------------
#  Electron Efficiencies
#==============================================================


#--------------------------------------------------------------
#  Electron Triggers
#==============================================================

####################################
# 
# Double Electron Triggers
#
# HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL : 150000 -> 170053
# HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL : 150000 -> 170053
# HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL  170054 -> 999999
# 
#
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011\",0\)
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053\",1\)
# root -l -q selectDoubleElectronLeadingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999\",2\)
#
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011\",0\)
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053\",1\)
# root -l -q selectDoubleElectronTrailingLegTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999\",2\)
#
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_170054-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elDoubleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_170054-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011/basic_76_106/","h2_results_electron_double_leadingleg","ElectronMVAIDIsoCombined_DoubleEleLeadingLeg_Full2011")'
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011/basic_76_106/","h2_results_electron_double_trailingleg","ElectronMVAIDIsoCombined_DoubleEleTrailingLeg_Full2011")'
#
#####################################
####################################
# 
# MuEG Triggers : Electron Leg
#
# HLT_Mu17_Ele8_CaloIdL : 150000 -> 173198
# HLT_Mu8_Ele17_CaloIdL : 150000 -> 170053
# HLT_Mu8_Ele17_CaloIdT_CaloIsoVL : 170054 -> 999999
# HLT_Mu17_Ele8_CaloIdT_CaloIsoVL : 173199 -> 999999
# 
#
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011\",0\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053\",1\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199\",2\)
# root -l -q selectMuEGElectronLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999\",3\)
#
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011\",0\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053\",1\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199\",2\)
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999\",3\)
#
# root -l -q selectMuEGElectronTrailingLegTrigEffTP.C+\(\"data_sel_Test.conf\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Test\",3\)
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_170054-173199/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronLeadingLeg_173200-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_170054-173199/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elMuEGEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_MuEGElectronTrailingLeg_173200-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#
#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_DoubleEleSeeded_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_DoubleEleSeeded_Full2011/basic_76_106/","h2_results_electron_double","ElectronMVAIDIsoCombined_DoubleEle_Full2011")'
#
#####################################
#
# Single Electron Triggers
# 
# HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT : 150000 -> 164237
# HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT : 165085 -> 166967
# HLT_Ele52_CaloIdVT_TrkIdT                  : 166968 -> 170053
# HLT_Ele65_CaloIdVT_TrkIdT                  : 170054 -> 178380
# HLT_Ele80_CaloIdVT_TrkIdT                  : 178381 -> 999999
#
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011\",0\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237\",1\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967\",2\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053\",3\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380\",4\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999\",5\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967\",10\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053\",11\)
# root -l -q selectSingleElectronTrigEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380\",12\)
#
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-164237/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_165085-166967/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_166968-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_170054-178380/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_178381-999999/basic_76_106\",\"all\",1,0,\"\"\)

#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-166967/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"elSingleEleTrig.bins\",0,0,0,0,\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380/probes.root\",\"Data_ElectronMVAIDIsoCombined_SingleEleSeeded_150000-178380/basic_76_106\",\"all\",1,0,\"\"\)


#
# root -l -q makeTriggerEfficiency.C+'("Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106/eff.root","Data_ElectronMVAIDIsoCombined_SingleEleSeeded_Full2011/basic_76_106/","h2_results_electron_single","ElectronMVAIDIsoCombined_SingleEle_Full2011")'
#
####################################





#--------------------------------------------------------------
#  BDT WithIPInfo ID + PFIso Efficiency
#==============================================================


####################################
# 
# Monte Carlo
#
# root -l -q selectEleBDTGWithIPInfoEffTP.C+\(\"s11-zee.conf\",\"Summer11_ElectronMVAIDIsoCombined\",1\)
#
####################################




####################################
# 
# Run2011A
#
#root -l -q selectEleBDTGWithIPInfoEffTP.C+\(\"data_el_Run2011A.conf\",\"Data_ElectronMVAIDIsoCombined_Run2011A\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToRun2011A\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Run2011A/probes.root\",\"Data_ElectronMVAIDIsoCombined_Run2011A/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_ElectronMVAIDIsoCombined_Run2011A/basic2_76_106/eff.root","Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToRun2011A/eff.root","Data_ElectronMVAIDIsoCombined_Run2011A/basic2_76_106/","h2_results_electron_selection","ElectronMVAIDIsoCombined_Run2011A")'
#
####################################



####################################
# 
# Run2011B
#
#root -l -q selectEleBDTGWithIPInfoEffTP.C+\(\"data_el_Run2011B.conf\",\"Data_ElectronMVAIDIsoCombined_Run2011B\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToRun2011B\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Run2011B/probes.root\",\"Data_ElectronMVAIDIsoCombined_Run2011B/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_ElectronMVAIDIsoCombined_Run2011B/basic2_76_106/eff.root","Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToRun2011B/eff.root","Data_ElectronMVAIDIsoCombined_Run2011B/basic2_76_106/","h2_results_electron_selection","ElectronMVAIDIsoCombined_Run2011B")'
#
####################################





####################################
# 
# Full2011
#
#root -l -q selectEleBDTGWithIPInfoEffTP.C+\(\"data_el.conf\",\"Data_ElectronMVAIDIsoCombined_Full2011\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elbarrel.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106_Barrel\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elendcap.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106_Endcap\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elLowPt.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106_LowPt\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"elHighPt.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106_HighPt\",\"all\",1,0,\"Summer11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106/eff.root","Summer11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToFull2011/eff.root","Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106/","h2_results_electron_selection", "ElectronMVAIDIsoCombined_Full2011")'
#
####################################







#--------------------------------------------------------------
#  Cutbased ID + PFIso Efficiency
#==============================================================


####################################
# 
# Monte Carlo
#
# root -l -q selectEleWPEffTP.C+\(\"s11-zee.conf\",\"Summer11_ElectronSmurfV6\",1\)
#
####################################




####################################
# 
# Run2011A
#
#root -l -q selectEleWPEffTP.C+\(\"data_el_Run2011A.conf\",\"Data_ElectronSmurfV6_Run2011A\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronSmurfV6/probes.root\",\"Summer11_ElectronSmurfV6/basic_76_106_ReweightedToRun2011A\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronSmurfV6_Run2011A/probes.root\",\"Data_ElectronSmurfV6_Run2011A/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronSmurfV6/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#
# root -l -q makeElectronScaleFactors.C+'("Data_ElectronSmurfV6_Run2011A/basic2_76_106/eff.root","Summer11_ElectronSmurfV6/basic_76_106_ReweightedToRun2011A/eff.root","Data_ElectronSmurfV6_Run2011A/basic2_76_106/","ElectronSmurfV6_Run2011A")'
#
####################################



####################################
# 
# Run2011B
#
#root -l -q selectEleWPEffTP.C+\(\"data_el_Run2011B.conf\",\"Data_ElectronSmurfV6_Run2011B\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronSmurfV6/probes.root\",\"Summer11_ElectronSmurfV6/basic_76_106_ReweightedToRun2011B\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronSmurfV6_Run2011B/probes.root\",\"Data_ElectronSmurfV6_Run2011B/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronSmurfV6/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#
# root -l -q makeElectronScaleFactors.C+'("Data_ElectronSmurfV6_Run2011B/basic2_76_106/eff.root","Summer11_ElectronSmurfV6/basic_76_106_ReweightedToRun2011B/eff.root","Data_ElectronSmurfV6_Run2011B/basic2_76_106/","ElectronSmurfV6_Run2011B")'
#
####################################





####################################
# 
# Full2011
#
#root -l -q selectEleWPEffTP.C+\(\"data_el.conf\",\"Data_ElectronSmurfV6_Full2011\"\)
#
#
#root -l -b -q plotEff.C+\(\"el0.bins\",0,0,0,0,\"Summer11_ElectronSmurfV6/probes.root\",\"Summer11_ElectronSmurfV6/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronSmurfV6_Full2011/probes.root\",\"Data_ElectronSmurfV6_Full2011/basic2_76_106\",\"all\",1,0,\"Summer11_ElectronSmurfV6/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#
# root -l -q makeElectronScaleFactors.C+'("Data_ElectronSmurfV6_Full2011/basic2_76_106/eff.root","Summer11_ElectronSmurfV6/basic_76_106_ReweightedToFull2011/eff.root","Data_ElectronSmurfV6_Full2011/basic2_76_106/","ElectronSmurfV6_Full2011")'
#
####################################




