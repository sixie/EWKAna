
#Produce pT Spectra for training sample
root -l PlotMVATrainingSpectra.C+

########################################################################################################################################################################
#Electrons
########################################################################################################################################################################

#Make Electron MVA ROC Curves
cd /data/blue/sixie/releases/analysis/CMSSW_4_2_3_patch2/src
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet0LowPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet0LowPtV18.root","Subdet0LowPt",0)'
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet1LowPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet1LowPtV18.root","Subdet1LowPt",1)'
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet2LowPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet2LowPtV18.root","Subdet2LowPt",2)'
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet0HighPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet0HighPtV18.root","Subdet0HighPt",3)'
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet1HighPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet1HighPtV18.root","Subdet1HighPt",4)'
root -l -q -b EWKAna/Thesis/EventReconstruction/MakeElectronIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Real.Subdet2HighPtV18.root","/data/blue/sixie/Thesis/ElectronMVA2/output/ElectronNtuple.Fake.Subdet2HighPtV18.root","Subdet2HighPt",5)'

#Make Electron MVA Input Distributions, Correlations, MVA output
cd /data/blue/sixie/Thesis/ElectronMVA2/TMVATest
root -l TMVAGui.C+'("../Subdet2LowPt_V15.root")'
#click on variables (uses my modified variables.C code)

cp variables_id_c1.eps    ElectronMVAInput_SigmaIEtaIEta.eps
cp variables_id_c2.eps	  ElectronMVAInput_DEtaIn.eps
cp variables_id_c3.eps	  ElectronMVAInput_DPhiIn.eps
cp variables_id_c4.eps	  ElectronMVAInput_D0.eps
cp variables_id_c5.eps	  ElectronMVAInput_FBrem.eps
cp variables_id_c6.eps	  ElectronMVAInput_EOverP.eps
cp variables_id_c7.eps	  ElectronMVAInput_ESeedClusterOverPout.eps
cp variables_id_c8.eps	  ElectronMVAInput_SigmaIPhiIPhi.eps
cp variables_id_c9.eps	  ElectronMVAInput_OneOverEMinusOneOverP.eps
cp variables_id_c10.eps   ElectronMVAInput_ESeedClusterOverPIn.eps
cp variables_id_c11.eps   ElectronMVAInput_IP3d.eps
cp variables_id_c12.eps   ElectronMVAInput_IP3dSig.eps
cp variables_id_c13.eps   ElectronMVAInput_GsfTrackChi2OverNdof.eps
cp variables_id_c14.eps   ElectronMVAInput_dEtaCalo.eps
cp variables_id_c15.eps   ElectronMVAInput_dPhiCalo.eps
cp variables_id_c16.eps   ElectronMVAInput_R9.eps
cp variables_id_c17.eps   ElectronMVAInput_SCEtaWidth.eps
cp variables_id_c18.eps   ElectronMVAInput_SCPhiWidth.eps
cp variables_id_c19.eps   ElectronMVAInput_CovIEtaIPhi.eps
cp variables_id_c20.eps   ElectronMVAInput_ChargedIso03.eps
cp variables_id_c21.eps   ElectronMVAInput_NeutralHadronIso03.eps
cp variables_id_c22.eps   ElectronMVAInput_GammaIso03.eps
cp variables_id_c23.eps   ElectronMVAInput_ChargedIso04.eps
cp variables_id_c24.eps   ElectronMVAInput_NeutralHadronIso04.eps
cp variables_id_c25.eps   ElectronMVAInput_GammaIso04.eps


#click on linear correlations (uses my modified correlations.C code)
cp CorrelationMatrixS.eps   ElectronMVA_CorrelationMatrixSignal.eps
cp CorrelationMatrixB.eps   ElectronMVA_CorrelationMatrixBkg.eps

#click on mva output (uses my modified mvas.C code)
cp overtrain_BDTG.eps ElectronBDTOutput.eps


####################################
#
# 
#Make Electron Efficiency Fit Plots
#
#
####################################
#
#
cd  /data/blue/sixie/Thesis/Efficiency
root -l -b -q plotEff_ForThesis_ElectronSelection.C+\(\"el0.bins\",0,0,0,0,\"Fall11_ElectronMVAIDIsoCombined/probes.root\",\"Fall11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
root -l -b -q plotEff_ForThesis_ElectronSelection.C+\(\"el0.bins\",2,1,2,3,\"Data_ElectronMVAIDIsoCombined_Full2011/probes.root\",\"Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106\",\"all\",1,0,\"Fall11_ElectronMVAIDIsoCombined/probes.root\",\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)

# Compute Efficiency Scale Factors
 root -l -q makeEfficiencyScaleFactors.C+'("Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106/eff.root","Fall11_ElectronMVAIDIsoCombined/basic_76_106_ReweightedToFull2011/eff.root","Data_ElectronMVAIDIsoCombined_Full2011/basic2_76_106/","h2_results_electron_selection","ElectronMVAIDIsoCombined_Full2011")'

#
#
####################################
ref




#########################################################################################################################################################################
#Muons
#########################################################################################################################################################################

#Make Muon MVA ROC Curves
cd /data/blue/sixie/releases/analysis/CMSSW_4_2_3_patch2/src
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.BarrelPtBin0_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.BarrelPtBin0_V10.root","BarrelPtBin0",0)'
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.EndcapPtBin0_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.EndcapPtBin0_V10.root","EndcapPtBin0",1)'
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.BarrelPtBin1_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.BarrelPtBin1_V10.root","BarrelPtBin1",2)'
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.EndcapPtBin1_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.EndcapPtBin1_V10.root","EndcapPtBin1",3)'
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.BarrelPtBin2_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.BarrelPtBin2_V10.root","BarrelPtBin2",4)'
root -l  -b -q EWKAna/Thesis/EventReconstruction/MakeMuonIDMVAPerformancePlots.C+'("/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Real.EndcapPtBin2_V10.root","/data/blue/sixie/Thesis/MuonMVA2/output/MuonNtuple.Fake.EndcapPtBin2_V10.root","EndcapPtBin2",5)'


#Make Electron MVA Input Distributions, Correlations, MVA output
cd /data/blue/sixie/Thesis/MuonMVA2/TMVATest
root -l TMVAGui.C+'("../BarrelPtBin2_V10.root")'
#click on variables (uses my modified variables.C code)

cp variables_id_c1.eps    MuonMVAInput_TkNchi2.eps
cp variables_id_c2.eps	  MuonMVAInput_GlobalNchi2.eps
cp variables_id_c3.eps	  MuonMVAInput_NValidHits.eps
cp variables_id_c4.eps	  MuonMVAInput_NTrackerHits.eps
cp variables_id_c5.eps	  MuonMVAInput_NPixelHits.eps
cp variables_id_c6.eps	  MuonMVAInput_NMatches.eps
cp variables_id_c7.eps	  MuonMVAInput_D0.eps
cp variables_id_c8.eps	  MuonMVAInput_IP3d.eps
cp variables_id_c9.eps	  MuonMVAInput_IP3dSig.eps
cp variables_id_c10.eps   MuonMVAInput_TrkKink.eps
cp variables_id_c11.eps   MuonMVAInput_SegmentCompatibility.eps
cp variables_id_c12.eps   MuonMVAInput_CaloCompatibility.eps
cp variables_id_c13.eps   MuonMVAInput_HadEnergyOverPt.eps
cp variables_id_c14.eps   MuonMVAInput_EmEnergyOverPt.eps
cp variables_id_c15.eps   MuonMVAInput_HadS9EnergyOverPt.eps
cp variables_id_c16.eps   MuonMVAInput_EmS9EnergyOverPt.eps
cp variables_id_c17.eps   MuonMVAInput_TrkIso03OverPt.eps
cp variables_id_c18.eps   MuonMVAInput_EMIso03OverPt.eps
cp variables_id_c19.eps   MuonMVAInput_HadIso03OverPt.eps
cp variables_id_c20.eps   MuonMVAInput_TrkIso05OverPt.eps
cp variables_id_c21.eps   MuonMVAInput_EMIso05OverPt.eps
cp variables_id_c22.eps   MuonMVAInput_HadIso05OverPt.eps


#click on linear correlations (uses my modified correlations.C code)
cp CorrelationMatrixS.eps   MuonMVA_CorrelationMatrixSignal.eps
cp CorrelationMatrixB.eps   MuonMVA_CorrelationMatrixBkg.eps

#click on mva output (uses my modified mvas.C code)
cp overtrain_BDTG.eps MuonBDTOutput.eps




####################################
#
# 
#Make Muon Efficiency Fit Plots
#
#
####################################
#
#
cd  /data/blue/sixie/Thesis/Efficiency
root -l -b -q plotEff_ForThesis_MuonSelection.C+\(\"mu0.bins\",0,0,0,0,\"Fall11_MuonMVAIDIsoCombinedDetIso/probes.root\",\"Fall11_MuonMVAIDIsoCombinedDetIso/basic_76_106_ReweightedToFull2011\",\"all\",1,0,\"\",\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)
root -l -b -q plotEff_ForThesis_MuonSelection.C+\(\"mu0.bins\",2,1,2,3,\"Data_MuonMVAIDIsoCombinedDetIso_Full2011/probes.root\",\"Data_MuonMVAIDIsoCombinedDetIso_Full2011/basic2_76_106\",\"all\",1,0,\"Fall11_MuonMVAIDIsoCombinedDetIso/probes.root\",\"/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/PileupReweighting.Fall11DYmm_To_Full2011.root\"\)

# Compute Efficiency Scale Factors
 root -l -q makeEfficiencyScaleFactors.C+'("Data_MuonMVAIDIsoCombinedDetIso_Full2011/basic2_76_106/eff.root","Fall11_MuonMVAIDIsoCombinedDetIso/basic_76_106_ReweightedToFull2011/eff.root","Data_MuonMVAIDIsoCombinedDetIso_Full2011/basic2_76_106/","h2_results_electron_selection","MuonMVAIDIsoCombinedDetIso_Full2011")'

#
#
####################################




#############################################################################################################################
#MET
#############################################################################################################################
root -l MakeMetCorrelationPlot.C+

