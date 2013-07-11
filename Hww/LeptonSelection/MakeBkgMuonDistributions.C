//root -l EWKAna/Hww/LeptonSelection/MakeBkgMuonDistributions.C+\(\)
//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TLegend.h>     
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <MitStyle.h>

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType);
Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void MakeBkgMuon(const string inputFilename, const string label, Int_t SelectionType, Double_t PtMin , Double_t PtMax  );

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

//=== MAIN MACRO =================================================================================================
void MakeBkgMuonDistributions() {

//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10Jet30_r11a-smu-pr-v1_Pt10To20",0, 10, 20);
// //   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10_Pt10To20",1, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Mu10Jet30_Pt10To20",10, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Mu10Jet30_Pt10To20",10, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt10To20",100, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt10To20",100, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt10To20",200, 10, 20);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt10To20",200, 10, 20);



//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10Jet30_r11a-smu-pr-v1_Pt20To35",0, 20, 35);
// //   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10_Pt20To35",1, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Mu10Jet30_Pt20To35",10, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Mu10Jet30_Pt20To35",10, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To35",100, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt20To35",100, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt20To35",200, 20, 35);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt20To35",200, 20, 35);



//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10Jet30_r11a-smu-pr-v1_Pt35To50",0, 35, 50);
// //   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v1_MuonFakeRateTriggerSkim.root","Mu10_Pt35To50",1, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Mu10Jet30_Pt35To50",10, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Mu10Jet30_Pt35To50",10, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt35To50",100, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt35To50",100, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt35To50",200, 35, 50);
//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt35To50",200, 35, 50);

//   MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130",100, 10, 200);






  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h115ww2l-gf-v1g1-pu_noskim_normalized.root","HWW115_Pt10To15", 100, 10, 15);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt10To15", 100, 10, 15);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/old_020_PFIsoStudy/WWAnalysisSkimmed_r11a-smu-pr_MuonFakeRateTriggerSkim.root","Mu10Jet30_Pt10To15", 0, 10, 15);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt10To15", 100, 10, 15);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt10To15", 10, 10, 15);

  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h115ww2l-gf-v1g1-pu_noskim_normalized.root","HWW115_Pt15To20", 100, 15, 20);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt15To20", 100, 15, 20);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/old_020_PFIsoStudy/WWAnalysisSkimmed_r11a-smu-pr_MuonFakeRateTriggerSkim.root","Mu10Jet30_Pt15To20", 0, 15, 20);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt15To20", 100, 15, 20);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt15To20", 10, 15, 20);

  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h115ww2l-gf-v1g1-pu_noskim_normalized.root","HWW115_Pt20To30", 100, 20, 30);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To30", 100, 20, 30);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/2011Data/old_020_PFIsoStudy/WWAnalysisSkimmed_r11a-smu-pr_MuonFakeRateTriggerSkim.root","Mu10Jet30_Pt20To30", 0, 20, 30);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root","Zmm_Pt20To30", 100, 20, 30);
  MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/old_020_PFIsoStudy/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt20To30", 10, 20, 30);




















//    MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_PtSelected",100, 10, 2000);
//    MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt10To20",100, 10, 20);
//    MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To30",100, 20, 30);
//    MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt30To50",100, 30, 50);
//    MakeBkgMuon("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20",100, 20, 2000);




}



void MakeBkgMuon(const string inputFilename,
                 const string Label, Int_t SelectionType, Double_t PtMin, Double_t PtMax) 
{  

  gBenchmark->Start("WWTemplate");

  string label = Label; if (Label != "") label = "_"+Label;
  
  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Double_t jetPtThreshold = 30;

  //*****************************************************************************************
  //Histogram
  //*****************************************************************************************
  vector<Double_t> MuonIsolationEfficiencyNumerator;
  vector<Double_t> MuonL1CorrectedIsolationEfficiencyNumerator;
  vector<Double_t> MuonIsolationEfficiencyDenominator;
  vector<Double_t> MuonL1CorrectedIsolationEfficiencyDenominator;
  for(UInt_t i=0 ; i<21; ++i) {
    MuonIsolationEfficiencyNumerator.push_back(0.0);
    MuonL1CorrectedIsolationEfficiencyNumerator.push_back(0.0);
    MuonIsolationEfficiencyDenominator.push_back(0.0);
    MuonL1CorrectedIsolationEfficiencyDenominator.push_back(0.0);
  }

  vector <TH1F*>  NVertices;
  vector<TH1F*>   Rho;
  vector<TH1F*>   Muon_caloIso_Barrel;
  vector<TH1F*>   Muon_ecalIso_Barrel;
  vector<TH1F*>   Muon_hcalIso_Barrel;
  vector<TH1F*>   Muon_trkIso_Barrel;
  vector<TH1F*>   Muon_caloIso05_Barrel;
  vector<TH1F*>   Muon_ecalIso05_Barrel;
  vector<TH1F*>   Muon_hcalIso05_Barrel;
  vector<TH1F*>   Muon_trkIso05_Barrel;
  vector<TH1F*>   Muon_caloIso_Endcap;
  vector<TH1F*>   Muon_ecalIso_Endcap;
  vector<TH1F*>   Muon_hcalIso_Endcap;
  vector<TH1F*>   Muon_trkIso_Endcap;
  vector<TH1F*>   Muon_caloIso05_Endcap;
  vector<TH1F*>   Muon_ecalIso05_Endcap;
  vector<TH1F*>   Muon_hcalIso05_Endcap;
  vector<TH1F*>   Muon_trkIso05_Endcap;

  vector<TH1F*>   Muon_ChargedIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_NeutralIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_ChargedIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_NeutralIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap;


  vector<TH1F*>   Muon_RelIso_Barrel;
  vector<TH1F*>   Muon_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   Muon_RelIso_Endcap;
  vector<TH1F*>   Muon_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   Muon_RelIso05_Barrel;
  vector<TH1F*>   Muon_RelIso05RhoCorrected_Barrel;
  vector<TH1F*>   Muon_RelIso05_Endcap;
  vector<TH1F*>   Muon_RelIso05RhoCorrected_Endcap;

  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;




  TH1F *Muon_relIso = new TH1F((string("Muon_relIso_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Muon_relIsoL1Corrected = new TH1F((string("Muon_relIsoL1Corrected_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Muon_relIso05 = new TH1F((string("Muon_relIso05_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Muon_relIso05L1Corrected = new TH1F((string("Muon_relIso05L1Corrected_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);


  for (int n=0; n < 20; ++n) {     
    TH1F *tmpRho = new TH1F((string("RhoMuon_") + IntToString(n) + label).c_str(), "; Rho; Number of Events ", 2000, -20, 20);
    TH1F *tmpMuon_caloIso_Barrel = new TH1F((string("Muon_caloIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ecalIso_Barrel = new TH1F((string("Muon_ecalIso_Barrel_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_hcalIso_Barrel = new TH1F((string("Muon_hcalIso_Barrel_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_trkIso_Barrel = new TH1F((string("Muon_trkIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_caloIso_Endcap = new TH1F((string("Muon_caloIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ecalIso_Endcap = new TH1F((string("Muon_ecalIso_Endcap_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_hcalIso_Endcap = new TH1F((string("Muon_hcalIso_Endcap_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_trkIso_Endcap = new TH1F((string("Muon_trkIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_caloIso05_Barrel = new TH1F((string("Muon_caloIso05_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ecalIso05_Barrel = new TH1F((string("Muon_ecalIso05_Barrel_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_hcalIso05_Barrel = new TH1F((string("Muon_hcalIso05_Barrel_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_trkIso05_Barrel = new TH1F((string("Muon_trkIso05_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_caloIso05_Endcap = new TH1F((string("Muon_caloIso05_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ecalIso05_Endcap = new TH1F((string("Muon_ecalIso05_Endcap_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_hcalIso05_Endcap = new TH1F((string("Muon_hcalIso05_Endcap_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_trkIso05_Endcap = new TH1F((string("Muon_trkIso05_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);


    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_01Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_05Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_10Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone03_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFIso_Cone04_15Threshold; Fraction of Events ", 2000, -200, 200.0);




    TH1F *tmpMuon_RelIso_Barrel = new TH1F((string("Muon_RelIso_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIsoRhoCorrected_Barrel = new TH1F((string("Muon_RelIsoRhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIso_Endcap = new TH1F((string("Muon_RelIso_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIsoRhoCorrected_Endcap = new TH1F((string("Muon_RelIsoRhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIso05_Barrel = new TH1F((string("Muon_RelIso05_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso05; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIso05RhoCorrected_Barrel = new TH1F((string("Muon_RelIso05RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso05RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIso05_Endcap = new TH1F((string("Muon_RelIso05_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso05; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_RelIso05RhoCorrected_Endcap = new TH1F((string("Muon_RelIso05RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso05RhoCorrected; Fraction of Events ", 2000, -4, 4.0);



    TH1F *tmpMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);




    tmpRho->StatOverflows(kTRUE);
    tmpMuon_caloIso_Barrel->StatOverflows(kTRUE);
    tmpMuon_ecalIso_Barrel->StatOverflows(kTRUE);
    tmpMuon_hcalIso_Barrel->StatOverflows(kTRUE);
    tmpMuon_trkIso_Barrel->StatOverflows(kTRUE);
    tmpMuon_caloIso_Endcap->StatOverflows(kTRUE);
    tmpMuon_ecalIso_Endcap->StatOverflows(kTRUE);
    tmpMuon_hcalIso_Endcap->StatOverflows(kTRUE);
    tmpMuon_trkIso_Endcap->StatOverflows(kTRUE);
    tmpMuon_caloIso05_Barrel->StatOverflows(kTRUE);
    tmpMuon_ecalIso05_Barrel->StatOverflows(kTRUE);
    tmpMuon_hcalIso05_Barrel->StatOverflows(kTRUE);
    tmpMuon_trkIso05_Barrel->StatOverflows(kTRUE);
    tmpMuon_caloIso05_Endcap->StatOverflows(kTRUE);
    tmpMuon_ecalIso05_Endcap->StatOverflows(kTRUE);
    tmpMuon_hcalIso05_Endcap->StatOverflows(kTRUE);
    tmpMuon_trkIso05_Endcap->StatOverflows(kTRUE);

    tmpMuon_ChargedIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_ChargedIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_NeutralIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);



    tmpMuon_RelIso_Barrel->StatOverflows(kTRUE);
    tmpMuon_RelIsoRhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_RelIso_Endcap->StatOverflows(kTRUE);
    tmpMuon_RelIsoRhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_RelIso05_Barrel->StatOverflows(kTRUE);
    tmpMuon_RelIso05RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_RelIso05_Endcap->StatOverflows(kTRUE);
    tmpMuon_RelIso05RhoCorrected_Endcap->StatOverflows(kTRUE);



    tmpMuon_TotalPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);



    Rho.push_back(tmpRho);
    Muon_caloIso_Barrel.push_back(tmpMuon_caloIso_Barrel);
    Muon_ecalIso_Barrel.push_back(tmpMuon_ecalIso_Barrel);
    Muon_hcalIso_Barrel.push_back(tmpMuon_hcalIso_Barrel);
    Muon_trkIso_Barrel.push_back(tmpMuon_trkIso_Barrel);
    Muon_caloIso_Endcap.push_back(tmpMuon_caloIso_Endcap);
    Muon_ecalIso_Endcap.push_back(tmpMuon_ecalIso_Endcap);
    Muon_hcalIso_Endcap.push_back(tmpMuon_hcalIso_Endcap);
    Muon_trkIso_Endcap.push_back(tmpMuon_trkIso_Endcap);
    Muon_caloIso05_Barrel.push_back(tmpMuon_caloIso05_Barrel);       
    Muon_ecalIso05_Barrel.push_back(tmpMuon_ecalIso05_Barrel);       
    Muon_hcalIso05_Barrel.push_back(tmpMuon_hcalIso05_Barrel);       
    Muon_trkIso05_Barrel.push_back(tmpMuon_trkIso05_Barrel);
    Muon_caloIso05_Endcap.push_back(tmpMuon_caloIso05_Endcap);
    Muon_ecalIso05_Endcap.push_back(tmpMuon_ecalIso05_Endcap);
    Muon_hcalIso05_Endcap.push_back(tmpMuon_hcalIso05_Endcap);
    Muon_trkIso05_Endcap.push_back(tmpMuon_trkIso05_Endcap);


    Muon_ChargedIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    Muon_NeutralIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel);
    Muon_TotalPFIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel);
    Muon_ChargedIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    Muon_NeutralIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap);
    Muon_TotalPFIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap);
    Muon_ChargedIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    Muon_NeutralIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel);
    Muon_TotalPFIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel);
    Muon_ChargedIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    Muon_NeutralIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap);
    Muon_TotalPFIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap);
    Muon_ChargedIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    Muon_NeutralIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel);
    Muon_TotalPFIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel);
    Muon_ChargedIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    Muon_NeutralIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap);
    Muon_TotalPFIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap);
    Muon_ChargedIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    Muon_NeutralIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel);
    Muon_TotalPFIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel);
    Muon_ChargedIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    Muon_NeutralIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap);
    Muon_TotalPFIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap);
    Muon_ChargedIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    Muon_NeutralIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel);
    Muon_TotalPFIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel);
    Muon_ChargedIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    Muon_NeutralIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap);
    Muon_TotalPFIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap);
    Muon_ChargedIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    Muon_NeutralIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel);
    Muon_TotalPFIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel);
    Muon_ChargedIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    Muon_NeutralIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap);
    Muon_TotalPFIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap);
    Muon_ChargedIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    Muon_NeutralIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel);
    Muon_TotalPFIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel);
    Muon_ChargedIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    Muon_NeutralIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap);
    Muon_TotalPFIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap);
    Muon_ChargedIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel);
    Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    Muon_NeutralIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel);
    Muon_TotalPFIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel);
    Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel);
    Muon_ChargedIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap);
    Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    Muon_NeutralIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap);
    Muon_TotalPFIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap);
    Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap);

    Muon_RelIso_Barrel.push_back(tmpMuon_RelIso_Barrel);
    Muon_RelIsoRhoCorrected_Barrel.push_back(tmpMuon_RelIsoRhoCorrected_Barrel);
    Muon_RelIso05_Barrel.push_back(tmpMuon_RelIso05_Barrel);
    Muon_RelIso05RhoCorrected_Barrel.push_back(tmpMuon_RelIso05RhoCorrected_Barrel);
    Muon_RelIso_Endcap.push_back(tmpMuon_RelIso_Endcap);
    Muon_RelIsoRhoCorrected_Endcap.push_back(tmpMuon_RelIsoRhoCorrected_Endcap);
    Muon_RelIso05_Endcap.push_back(tmpMuon_RelIso05_Endcap);
    Muon_RelIso05RhoCorrected_Endcap.push_back(tmpMuon_RelIso05RhoCorrected_Endcap);
    
    Muon_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_01Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_01Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_01Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_01Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_05Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_05Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_05Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_05Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_10Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_10Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_10Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_10Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_15Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_15Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Muon_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_15Threshold_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Muon_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_15Threshold_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);

  }
 


  ofstream eventListFile("eventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt"); 
  hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

  cout << "Total Events: " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear(); 
    jetArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    jetBr->GetEntry(ientry);


    //********************************************************
    // TcMet
    //********************************************************
    TVector3 met;        
    if(info->tcMEx!=0 || info->tcMEy!=0) {       
      met.SetXYZ(info->tcMEx, info->tcMEy, 0);
    }
	
    Double_t PUIsolationEnergy = info->PileupEnergyDensity * 3.14159 * pow(0.3,2) * 1.0 ;
    Double_t PUIsolationEnergy04Cone = info->PileupEnergyDensity * 3.14159 * pow(0.4,2) * 1.0 ;
    Double_t PUIsolationEnergy05Cone = info->PileupEnergyDensity * 3.14159 * pow(0.5,2) * 1.0 ;
//     Double_t PUIsolationEnergy = info->PileupEnergyDensity * 0.0931;
//     Double_t PUIsolationEnergy04Cone = info->PileupEnergyDensity * 0.2;
//     Double_t PUIsolationEnergy05Cone = info->PileupEnergyDensity * 0.28;
    Int_t NVertex = info->nPV0; if (NVertex > 19) NVertex = 19;
//     NVertex = 1;

    //********************************************************
    // Event Selection Cuts
    //********************************************************
    if (SelectionType < 100) {
      if (SelectionType < 10) {
        if (!(
              (info->triggerBits & kHLT_Mu8)
              || (info->triggerBits & kHLT_Mu15)
              )
          ) continue;
      }
      
      if (met.Pt() > 20) continue;
      if (muonArr->GetEntries() > 1) continue;
    }

    Double_t tempLeadingJetPt = 0;
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
      if( jet->pt > tempLeadingJetPt) tempLeadingJetPt = jet->pt;
    }
    
     
    //********************************************************************************************
    // Build 20-10 signal sample
    //********************************************************************************************
    Bool_t passDileptonPreselection = kFALSE;
    vector<Double_t> realLeptonPt;
    vector<Double_t> realLeptonType;
    if (SelectionType == 100) {
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

        Int_t originalvalue = ele->isMCReal;
        Int_t isEle;
        Int_t isMu;
        Int_t isTau;
      
        isEle = originalvalue % 2;
        Int_t tmp0 = floor(double(originalvalue) / 2.0);
        isMu = tmp0 % 2;
        Int_t tmp1 = floor(double(tmp0) / 2.0);
        isTau = tmp1 % 2;
      
        if ((isEle)) {
          realLeptonPt.push_back(ele->pt);
          realLeptonType.push_back(11);
        }
      }
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        Int_t originalvalue = mu->isMCReal;
        Int_t isEle;
        Int_t isMu;
        Int_t isTau;
      
        isEle = originalvalue % 2;
        Int_t tmp0 = floor(double(originalvalue) / 2.0);
        isMu = tmp0 % 2;
        Int_t tmp1 = floor(double(tmp0) / 2.0);
        isTau = tmp1 % 2;
      
        if ((isMu)) {
          realLeptonPt.push_back(mu->pt);
          realLeptonType.push_back(13);
        }
      }
      if (realLeptonPt.size() == 2) {
        Bool_t pass = kTRUE;
        if (!(realLeptonPt[0] > 20 || realLeptonPt[1] > 20)) pass = kFALSE;
        if (!( (realLeptonType[0] == 11 && realLeptonPt[0] > 15)
               ||
               (realLeptonType[0] == 13 && realLeptonPt[0] > 10)
              )
          ) pass = kFALSE;
        if (!( (realLeptonType[1] == 11 && realLeptonPt[1] > 15)
               ||
               (realLeptonType[1] == 13 && realLeptonPt[1] > 10)
              )
          ) pass = kFALSE;
        if (pass) passDileptonPreselection = kTRUE;
      }
    }



     for(Int_t i=0; i<muonArr->GetEntries(); i++) {
       const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
       
       if(!(mu->pt > PtMin && mu->pt <= PtMax )) continue;



       //      //pass HLT selection
//       if (!passHLT(info->triggerBits, info->runNum, SelectionType)) continue;
       
       //pass event selection
       Bool_t passJetSelection = kFALSE;
       for(Int_t j=0; j<jetArr->GetEntries(); j++) {
         const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);
         
         if (jet->pt > jetPtThreshold &&
             mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 0.5) {
           passJetSelection = kTRUE;
           break;
         }
       }




       if (SelectionType < 100) {
         if (SelectionType != 1) {
           if (!passJetSelection) continue;
         }
         if (!(mu->pt > 10 && mu->pt < 30)) continue;
       }
       

       Int_t originalvalue = mu->isMCReal;
       Int_t isEle;
       Int_t isMu;
       Int_t isTau;
       
       isEle = originalvalue % 2;
       Int_t tmp0 = floor(double(originalvalue) / 2.0);
       isMu = tmp0 % 2;
       Int_t tmp1 = floor(double(tmp0) / 2.0);
       isTau = tmp1 % 2;


      if (SelectionType == 100) {
        if (!isMu) continue;
        if (!passDileptonPreselection) continue;
      }

      //WJets MC - pick only fake muons
      if (SelectionType == 200) {
//          if (isEle || isMu || isTau) {
        if (!(mu->isMCReal == 0)) {
          continue;
        }
       }
      
      if (!passMuonDenominatorCuts(mu)) continue;
      
       
       
       Rho[NVertex]->Fill(info->PileupEnergyDensity);
       if (fabs(mu->eta) < 1.5) {
         Muon_caloIso_Barrel[NVertex]->Fill(TMath::Min(double( mu->emIso03 + mu->hadIso03), double(199.99)));
         Muon_ecalIso_Barrel[NVertex]->Fill(TMath::Min(double( mu->emIso03 ), double(199.99)));
         Muon_hcalIso_Barrel[NVertex]->Fill(TMath::Min(double( mu->hadIso03), double(199.99)));
         Muon_trkIso_Barrel[NVertex]->Fill(TMath::Min(double( mu->trkIso03 ), double(199.99)));
         Muon_caloIso05_Barrel[NVertex]->Fill(TMath::Min( double(mu->emIso05 + mu->hadIso05), double(199.99)));
         Muon_ecalIso05_Barrel[NVertex]->Fill(TMath::Min( double(mu->emIso05 ), double(199.99)));
         Muon_hcalIso05_Barrel[NVertex]->Fill(TMath::Min( double(mu->hadIso05), double(199.99)));
         Muon_trkIso05_Barrel[NVertex]->Fill(TMath::Min(double( mu->trkIso05 ), double(199.99)));


         Muon_ChargedIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_05Threshold+mu->ChargedIso03FromOtherVertices_05Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_05Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_05Threshold+mu->ChargedIso04FromOtherVertices_05Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_05Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_10Threshold+mu->ChargedIso03FromOtherVertices_10Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_10Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_10Threshold+mu->ChargedIso04FromOtherVertices_10Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_10Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_15Threshold+mu->ChargedIso03FromOtherVertices_15Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_15Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_15Threshold+mu->ChargedIso04FromOtherVertices_15Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_15Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_15Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold) , double(199.99)));


    
         Muon_RelIso_Barrel[NVertex]->Fill(TMath::Min( double((mu->emIso03 + mu->hadIso03 + mu->trkIso03)/mu->pt) , double(3.99)));
         Muon_RelIsoRhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->emIso03 + mu->hadIso03 + mu->trkIso03 - info->PileupEnergyDensity*0.0884)/mu->pt) , double(3.99)));
         Muon_RelIso05_Barrel[NVertex]->Fill(TMath::Min( double((mu->emIso05 + mu->hadIso05 + mu->trkIso05)/mu->pt) , double(3.99)));
         Muon_RelIso05RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->emIso05 + mu->hadIso05 + mu->trkIso05 - info->PileupEnergyDensity*0.2302)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



       } else {
         Muon_caloIso_Endcap[NVertex]->Fill(TMath::Min(double( mu->emIso03 + mu->hadIso03), double(199.99)));
         Muon_ecalIso_Endcap[NVertex]->Fill(TMath::Min(double( mu->emIso03 ), double(199.99)));
         Muon_hcalIso_Endcap[NVertex]->Fill(TMath::Min(double( mu->hadIso03), double(199.99)));
         Muon_trkIso_Endcap[NVertex]->Fill(TMath::Min(double( mu->trkIso03 ), double(199.99)));
         Muon_caloIso05_Endcap[NVertex]->Fill(TMath::Min( double(mu->emIso05 + mu->hadIso05), double(199.99)));
         Muon_ecalIso05_Endcap[NVertex]->Fill(TMath::Min( double(mu->emIso05 ), double(199.99)));
         Muon_hcalIso05_Endcap[NVertex]->Fill(TMath::Min( double( mu->hadIso05), double(199.99)));
         Muon_trkIso05_Endcap[NVertex]->Fill(TMath::Min( double(mu->trkIso05 ) , double(199.99)));


         Muon_ChargedIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_05Threshold+mu->ChargedIso03FromOtherVertices_05Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_05Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_05Threshold+mu->ChargedIso04FromOtherVertices_05Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_05Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_10Threshold+mu->ChargedIso03FromOtherVertices_10Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_10Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_10Threshold+mu->ChargedIso04FromOtherVertices_10Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_10Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold) , double(199.99)));
         Muon_ChargedIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_15Threshold+mu->ChargedIso03FromOtherVertices_15Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_15Threshold) , double(199.99)));
         Muon_NeutralIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold) , double(199.99)));
         Muon_ChargedIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_15Threshold+mu->ChargedIso04FromOtherVertices_15Threshold) , double(199.99)));
         Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_15Threshold) , double(199.99)));
         Muon_NeutralIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->NeutralIso04_15Threshold) , double(199.99)));
         Muon_TotalPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold) , double(199.99)));
         Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold) , double(199.99)));



         Muon_RelIso_Endcap[NVertex]->Fill(TMath::Min( double((mu->emIso03 + mu->hadIso03 + mu->trkIso03)/mu->pt) , double(3.99)));
         Muon_RelIsoRhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->emIso03 + mu->hadIso03 + mu->trkIso03 - info->PileupEnergyDensity*0.0633)/mu->pt) , double(3.99)));
         Muon_RelIso05_Endcap[NVertex]->Fill(TMath::Min( double((mu->emIso05 + mu->hadIso05 + mu->trkIso05)/mu->pt) , double(3.99)));
         Muon_RelIso05RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->emIso05 + mu->hadIso05 + mu->trkIso05 - info->PileupEnergyDensity*0.1564)/mu->pt) , double(3.99)));


         Muon_TotalPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_01Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_01Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_01Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_01Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_05Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_05Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_05Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_05Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_10Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_10Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_10Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_10Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));



         Muon_TotalPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->ChargedIso03FromOtherVertices_01Threshold+mu->NeutralIso03_15Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso03_01Threshold+mu->NeutralIso03_15Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold)/mu->pt) , double(3.99)));
         Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->ChargedIso04FromOtherVertices_01Threshold+mu->NeutralIso04_15Threshold - info->PileupEnergyDensity*0.2054)/mu->pt) , double(3.99)));
         Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((mu->ChargedIso04_01Threshold+mu->NeutralIso04_15Threshold - info->PileupEnergyDensity*0.0455)/mu->pt) , double(3.99)));
       }
       
       Muon_relIso->Fill(TMath::Min(double( (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt), double(1.99)));
       Muon_relIsoL1Corrected->Fill(TMath::Min(double( (mu->trkIso03 + TMath::Max(mu->emIso03 + mu->hadIso03 - PUIsolationEnergy,0.0)) / mu->pt), double(1.99)));
       Muon_relIso05->Fill(TMath::Min( double((mu->trkIso05 + mu->emIso05 + mu->hadIso05) / mu->pt), double(1.99)));
       Muon_relIso05L1Corrected->Fill(TMath::Min(double((mu->trkIso05 + TMath::Max(mu->emIso05 + mu->hadIso05 - PUIsolationEnergy05Cone,0.0)) / mu->pt), double(1.99)));
       
       
       MuonIsolationEfficiencyDenominator[NVertex]++; 
       MuonIsolationEfficiencyDenominator[20]++; 
       if ((mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15) {
         MuonIsolationEfficiencyNumerator[NVertex]++;
         MuonIsolationEfficiencyNumerator[20]++;
       }         

       MuonL1CorrectedIsolationEfficiencyDenominator[NVertex]++;
       MuonL1CorrectedIsolationEfficiencyDenominator[20]++;
       if ((mu->trkIso03 + TMath::Max(mu->emIso03 + mu->hadIso03 - PUIsolationEnergy,0.0)) / mu->pt < 0.15) {
         MuonL1CorrectedIsolationEfficiencyNumerator[NVertex]++;
         MuonL1CorrectedIsolationEfficiencyNumerator[20]++;
       }
       

    }

  } //end loop over data     


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  
  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  const int nPoints = 20;
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.0;     
  }

  double MuonIsolationEff[nPoints];
  double MuonIsolationEffErrLow[nPoints];
  double MuonIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonIsolationEff[i] = ratio;
    MuonIsolationEffErrLow[i] = errLow;
    MuonIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonIsolationEff " << i << " : " << MuonIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonIsolationEff, NPileupError, NPileupError, MuonIsolationEffErrLow, MuonIsolationEffErrHigh);
  MuonIsolationEffVsNPileup->SetMarkerColor(kRed);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);


  Double_t r = 0;
  Double_t eLow = 0;
  Double_t eHigh = 0;     
  Double_t n1 = TMath::Nint(MuonIsolationEfficiencyNumerator[20]);
  Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[20]);
  mithep::MathUtils::CalcRatio(n1 , n2, r, eLow, eHigh, 2);
  cout << "MuonIsolationEff Averaged : " << r << " + " << eHigh << " - " << eLow << endl;
  


  double MuonL1CorrectedIsolationEff[nPoints];
  double MuonL1CorrectedIsolationEffErrLow[nPoints];
  double MuonL1CorrectedIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonL1CorrectedIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonL1CorrectedIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonL1CorrectedIsolationEff[i] = ratio;
    MuonL1CorrectedIsolationEffErrLow[i] = errLow;
    MuonL1CorrectedIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonL1CorrectedIsolationEff " << i << " : " << MuonL1CorrectedIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonL1CorrectedIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonL1CorrectedIsolationEff, NPileupError, NPileupError,MuonL1CorrectedIsolationEffErrLow,MuonL1CorrectedIsolationEffErrHigh  );
  MuonL1CorrectedIsolationEffVsNPileup->SetMarkerColor(kBlue);
  MuonL1CorrectedIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonL1CorrectedIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  Double_t ymin = 0.0;
  Double_t ymax = 1.1;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  legend = new TLegend(0.20,0.55,0.43,0.70);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.0;
  legend->AddEntry(MuonIsolationEffVsNPileup, "NoCorrection", "LP");
  legend->AddEntry(MuonL1CorrectedIsolationEffVsNPileup, "FastJetCorrected", "LP");
  MuonIsolationEffVsNPileup->SetTitle("");
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  MuonIsolationEffVsNPileup->Draw("AP");
  MuonL1CorrectedIsolationEffVsNPileup->Draw("Psame");
  MuonIsolationEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);

  legend->Draw();
  cv->SaveAs(("MuonIsolationEfficiency_BkgData_vs_NVertices"+label+".gif").c_str());



  //*****************************************************************************************
  //Save Efficiency Plots
  //*****************************************************************************************
  TFile *file = new TFile("HwwSelectionPlots_LeptonEfficiency.muons.root", "UPDATE");
  file->cd();
  file->WriteTObject(MuonIsolationEffVsNPileup, ("MuonIsolationEffVsNVertices" +label).c_str() , "WriteDelete");
  file->WriteTObject(MuonL1CorrectedIsolationEffVsNPileup, ("MuonFastjetCorrectedIsolationEffVsNVertices" +label).c_str() , "WriteDelete");

  file->WriteTObject(Muon_relIso, Muon_relIso->GetName(), "WriteDelete");
  file->WriteTObject(Muon_relIsoL1Corrected, Muon_relIsoL1Corrected->GetName(), "WriteDelete");
  file->WriteTObject(Muon_relIso05, Muon_relIso05->GetName(), "WriteDelete");
  file->WriteTObject(Muon_relIso05L1Corrected, Muon_relIso05L1Corrected->GetName(), "WriteDelete");

   for (int n=0; n < 20 ; ++n) {
      file->WriteTObject(Rho[n],Rho[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso_Barrel[n],Muon_caloIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso_Barrel[n],Muon_ecalIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso_Barrel[n],Muon_hcalIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso_Barrel[n],Muon_trkIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso05_Barrel[n],Muon_caloIso05_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso05_Barrel[n],Muon_ecalIso05_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso05_Barrel[n],Muon_hcalIso05_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso05_Barrel[n],Muon_trkIso05_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso_Endcap[n],Muon_caloIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso_Endcap[n],Muon_ecalIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso_Endcap[n],Muon_hcalIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso_Endcap[n],Muon_trkIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso05_Endcap[n],Muon_caloIso05_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso05_Endcap[n],Muon_ecalIso05_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso05_Endcap[n],Muon_hcalIso05_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso05_Endcap[n],Muon_trkIso05_Endcap[n]->GetName(), "WriteDelete") ;


      file->WriteTObject(Muon_ChargedIso_Cone03_01Threshold_Barrel[n],Muon_ChargedIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_01Threshold_Barrel[n],Muon_NeutralIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_01Threshold_Barrel[n],Muon_TotalPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_01Threshold_Endcap[n],Muon_ChargedIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_01Threshold_Endcap[n],Muon_NeutralIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_01Threshold_Endcap[n],Muon_TotalPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_01Threshold_Barrel[n],Muon_ChargedIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_01Threshold_Barrel[n],Muon_NeutralIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_01Threshold_Barrel[n],Muon_TotalPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_01Threshold_Endcap[n],Muon_ChargedIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_01Threshold_Endcap[n],Muon_NeutralIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_01Threshold_Endcap[n],Muon_TotalPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_05Threshold_Barrel[n],Muon_ChargedIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_05Threshold_Barrel[n],Muon_NeutralIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_05Threshold_Barrel[n],Muon_TotalPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_05Threshold_Endcap[n],Muon_ChargedIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_05Threshold_Endcap[n],Muon_NeutralIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_05Threshold_Endcap[n],Muon_TotalPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_05Threshold_Barrel[n],Muon_ChargedIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_05Threshold_Barrel[n],Muon_NeutralIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_05Threshold_Barrel[n],Muon_TotalPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_05Threshold_Endcap[n],Muon_ChargedIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_05Threshold_Endcap[n],Muon_NeutralIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_05Threshold_Endcap[n],Muon_TotalPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_10Threshold_Barrel[n],Muon_ChargedIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_10Threshold_Barrel[n],Muon_NeutralIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_10Threshold_Barrel[n],Muon_TotalPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_10Threshold_Endcap[n],Muon_ChargedIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_10Threshold_Endcap[n],Muon_NeutralIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_10Threshold_Endcap[n],Muon_TotalPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_10Threshold_Barrel[n],Muon_ChargedIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_10Threshold_Barrel[n],Muon_NeutralIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_10Threshold_Barrel[n],Muon_TotalPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_10Threshold_Endcap[n],Muon_ChargedIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_10Threshold_Endcap[n],Muon_NeutralIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_10Threshold_Endcap[n],Muon_TotalPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_15Threshold_Barrel[n],Muon_ChargedIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_15Threshold_Barrel[n],Muon_NeutralIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_15Threshold_Barrel[n],Muon_TotalPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone03_15Threshold_Endcap[n],Muon_ChargedIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone03_15Threshold_Endcap[n],Muon_NeutralIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone03_15Threshold_Endcap[n],Muon_TotalPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_15Threshold_Barrel[n],Muon_ChargedIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel[n],Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_15Threshold_Barrel[n],Muon_NeutralIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_15Threshold_Barrel[n],Muon_TotalPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel[n],Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Cone04_15Threshold_Endcap[n],Muon_ChargedIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap[n],Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Cone04_15Threshold_Endcap[n],Muon_NeutralIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Cone04_15Threshold_Endcap[n],Muon_TotalPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap[n],Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;


      file->WriteTObject(Muon_RelIso_Barrel[n],Muon_RelIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIsoRhoCorrected_Barrel[n],Muon_RelIsoRhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso_Endcap[n],Muon_RelIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIsoRhoCorrected_Endcap[n],Muon_RelIsoRhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso05_Barrel[n],Muon_RelIso05_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso05RhoCorrected_Barrel[n],Muon_RelIso05RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso05_Endcap[n],Muon_RelIso05_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso05RhoCorrected_Endcap[n],Muon_RelIso05RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;



      file->WriteTObject(Muon_TotalPFRelIso_Cone03_01Threshold_Barrel[n],Muon_TotalPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_01Threshold_Endcap[n],Muon_TotalPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_01Threshold_Barrel[n],Muon_TotalPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_01Threshold_Endcap[n],Muon_TotalPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_05Threshold_Barrel[n],Muon_TotalPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_05Threshold_Endcap[n],Muon_TotalPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_05Threshold_Barrel[n],Muon_TotalPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_05Threshold_Endcap[n],Muon_TotalPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_10Threshold_Barrel[n],Muon_TotalPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_10Threshold_Endcap[n],Muon_TotalPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_10Threshold_Barrel[n],Muon_TotalPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_10Threshold_Endcap[n],Muon_TotalPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_15Threshold_Barrel[n],Muon_TotalPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_15Threshold_Endcap[n],Muon_TotalPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_15Threshold_Barrel[n],Muon_TotalPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_15Threshold_Endcap[n],Muon_TotalPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;

    }


  file->Close();
  delete file;

    
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType) {


  Bool_t pass = kFALSE;


  return pass;

}


Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

        && (mu->nSeg > 1 || mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->pterr / mu->pt < 0.1)
        && (fabs(mu->dz) < 0.2)
        )
    ) pass = kFALSE;

  return pass;
}




Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
//         && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 1.0
        && (mu->pterr / mu->pt < 0.1)
        && (fabs(mu->dz) < 0.2)
        )
    ) pass = kFALSE;

  return pass;
}






//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs <<   runNum << " " ;
  ofs <<  lumiSec << " ";
  ofs << evtNum<< " ";
  ofs << mass<< " ";

//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
