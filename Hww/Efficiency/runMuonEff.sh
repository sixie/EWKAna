#--------------------------------------------------------------
#  Muon Efficiencies
#==============================================================


#--------------------------------------------------------------
#  Muon Triggers
#==============================================================

####################################
# 
# Double Muon Triggers
#
# HLT_DoubleMu7 : 150000 -> 164237
# HLT_Mu13_Mu8  : 165085 -> 178380
# HLT_Mu17_Mu8  : 178381 -> 999999
# 
#
# root -l -q selectDoubleMuLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011\",0\)
# root -l -q selectDoubleMuLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_150000-164237\",1\)
# root -l -q selectDoubleMuLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_165085-178380\",2\)
# root -l -q selectDoubleMuLeadingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_178381-999999\",3\)
#
# root -l -q selectDoubleMuTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011\",0\)
# root -l -q selectDoubleMuTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_150000-164237\",1\)
# root -l -q selectDoubleMuTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_165085-178380\",2\)
# root -l -q selectDoubleMuTrailingLegTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_178381-999999\",3\)
#
#
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_150000-164237/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_150000-164237/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_165085-178380/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_165085-178380/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_178381-999999/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_178381-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_150000-164237/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_150000-164237/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_165085-178380/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_165085-178380/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muDoubleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_178381-999999/probes.root\",\"Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_178381-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#
# root -l -q makeTriggerEfficiency.C+'("Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011/basic_76_106/eff.root","Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011/basic_76_106/","h2_results_muon_double_leadingleg","Data_MuonMVAIDIsoCombined_DoubleMuLeadingLegTrig_Full2011")'
# root -l -q makeTriggerEfficiency.C+'("Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011/basic_76_106/eff.root","Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011/basic_76_106/","h2_results_muon_double_trailingleg","Data_MuonMVAIDIsoCombined_DoubleMuTrailingLegTrig_Full2011")'
#
####################################
# 
# MuEG Triggers : Muon Leg
#
# HLT_Mu17_Ele8_CaloIdL : 150000 -> 173198
# HLT_Mu8_Ele17_CaloIdL : 150000 -> 170053
# HLT_Mu8_Ele17_CaloIdT_CaloIsoVL : 170054 -> 999999
# HLT_Mu17_Ele8_CaloIdT_CaloIsoVL : 173199 -> 999999
# 
#
# root -l -q selectMuEGMuonLeadingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_Full2011\",0\)
# root -l -q selectMuEGMuonLeadingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_150000-170053\",1\)
# root -l -q selectMuEGMuonLeadingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_170054-173199\",2\)
# root -l -q selectMuEGMuonLeadingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_173199-999999\",3\)
#
# root -l -q selectMuEGMuonTrailingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_Full2011\",0\)
# root -l -q selectMuEGMuonTrailingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_150000-170053\",1\)
# root -l -q selectMuEGMuonTrailingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_170054-173199\",2\)
# root -l -q selectMuEGMuonTrailingLegTrigEffTP.C+\(\"data_sel_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_173199-999999\",3\)
#
#
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_150000-170053/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_170054-173199/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_170054-173199/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_173199-999999/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonLeadingLeg_173199-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_150000-170053/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_150000-170053/basic_76_106\",\"all\",1,0,\"\"\)
#root -l -b -q plotEff.C+\(\"muMuEGTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_170054-999999/probes.root\",\"Data_MuonMVAIDIsoCombined_MuEGMuonTrailingLeg_170054-999999/basic_76_106\",\"all\",1,0,\"\"\)
#
#



# root -l -q makeTriggerEfficiency.C+'("Data_MuonMVAIDIsoCombined_MuEGMuonLeg_Full2011/basic_76_106/eff.root","Data_MuonMVAIDIsoCombined_MuEGMuonLeg_Full2011/basic_76_106/","h2_results_electron_double","MuonMVAIDIsoCombined_DoubleEle_Full2011")'
#
#####################################
#
# Single Muon Triggers
# 
# HLT_Mu15        : 150000 -> 163261
# HLT_Mu24        : 163262 -> 164237
# HLT_Mu30        : 165085 -> 170053
# HLT_Mu40        : 170054 -> 999999
# HLT_IsoMu17XXX  : 163262 -> 170053
# HLT_IsoMu24     : 170054 -> 173198
# HLT_IsoMu30     : 173199 -> 999999
#
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_Full2011\",0\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_150000-164237\",1\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_165085-166967\",2\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_166968-170053\",3\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_170054-178380\",4\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_178381-999999\",5\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_150000-166967\",10\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_150000-170053\",11\)
# root -l -q selectSingleMuTrigEffTP.C+\(\"data_smu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_SingleMu_150000-178380\",12\)
#
#root -l -b -q plotEff.C+\(\"muSingleMuTrig.bins\",0,0,0,0,\"Data_MuonMVAIDIsoCombined_SingleMu_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_SingleMu_Full2011/basic_76_106\",\"all\",1,0,\"\"\)
#
#
# root -l -q makeTriggerEfficiency.C+'("Data_MuonMVAIDIsoCombined_SingleMu_Full2011/basic_76_106/eff.root","Data_MuonMVAIDIsoCombined_SingleMu_Full2011/basic_76_106/","h2_results_muon_single","Data_MuonMVAIDIsoCombined_SingleMuTrig_Full2011")'
#
####################################


#--------------------------------------------------------------
#  Muon Reco Eff
#==============================================================

# root -l -q makeEfficiencyScaleFactors.C+'("Data_MuonReco_Full2011/basic_76_106/eff.root","Summer11_Zmm_MuonRecoEffTP_PU2011B/basic_76_106/eff.root","Data_MuonReco_Full2011/basic_76_106/","h2_results_muon_reco","Data_MuonReco_Full2011")'



#--------------------------------------------------------------
#  Muon ID + PFIso Efficiency
#==============================================================


####################################
# 
# Monte Carlo
#
# root -l -q selectMuonWPEffTP.C+\(\"s11-zmm.conf\",\"Summer11_MuonMVAIDIsoCombined\",1\)
#
####################################




####################################
# 
# Run2011A
#
#root -l -q selectMuonWPEffTP.C+\(\"data_mu_Run2011A.conf\",\"Data_MuonMVAIDIsoCombined_Run2011A\"\)
#
#
#root -l -b -q plotEff.C+\(\"mu0.bins\",0,0,0,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToRun2011A\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#root -l -b -q plotEff.C+\(\"mu0.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Run2011A/probes.root\",\"Data_MuonMVAIDIsoCombined_Run2011A/basic2_76_106\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011A.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_MuonMVAIDIsoCombined_Run2011A/basic2_76_106/eff.root","Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToRun2011A/eff.root","Data_MuonMVAIDIsoCombined_Run2011A/basic2_76_106/","h2_results_muon_selection","MuonMVAIDIsoCombined_Run2011A")'
#
####################################


####################################
# 
# Run2011B
#
#root -l -q selectMuonWPEffTP.C+\(\"data_mu_Run2011B.conf\",\"Data_MuonMVAIDIsoCombined_Run2011B\"\)
#
#
#root -l -b -q plotEff.C+\(\"mu0.bins\",0,0,0,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToRun2011B\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#root -l -b -q plotEff.C+\(\"mu0.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Run2011B/probes.root\",\"Data_MuonMVAIDIsoCombined_Run2011B/basic2_76_106\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Run2011B.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_MuonMVAIDIsoCombined_Run2011B/basic2_76_106/eff.root","Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToRun2011B/eff.root","Data_MuonMVAIDIsoCombined_Run2011B/basic2_76_106/","h2_results_muon_selection","MuonMVAIDIsoCombined_Run2011B")'
#
#
####################################


####################################
# 
# Full2011
#
#root -l -q selectMuonWPEffTP.C+\(\"data_mu_Full2011.conf\",\"Data_MuonMVAIDIsoCombined_Full2011\"\)
#
#
#root -l -b -q plotEff.C+\(\"mu0.bins\",0,0,0,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"mu0.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"mubarrel.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106_Barrel\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"muendcap.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106_Endcap\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"mulowpt.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106_LowPt\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#root -l -b -q plotEff.C+\(\"muhighpt.bins\",2,1,2,1,\"Data_MuonMVAIDIsoCombined_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106_HighPt\",\"all\",1,0,\"Summer11_MuonMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root\"\)
#
# root -l -q makeEfficiencyScaleFactors.C+'("Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106/eff.root","Summer11_MuonMVAIDIsoCombined/basic_76_106_ReweightedToFull2011/eff.root","Data_MuonMVAIDIsoCombined_Full2011/basic2_76_106/","h2_results_muon_selection","MuonMVAIDIsoCombined_Full2011")'
#
####################################


