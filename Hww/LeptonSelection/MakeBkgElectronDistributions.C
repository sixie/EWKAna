//root -l EWKAna/Hww/LeptonSelection/MakeBkgElectronDistributions.C+\(\)
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
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void MakeBkgElectron(const string inputFilename,
                     const string label, Int_t SelectionType, Double_t PtMin , Double_t PtMax );

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
void MakeBkgElectronDistributions() {

//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt10To20", 0, 10, 20);
// //   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/HwwAnalysis_r11a-del-pr-v1_EleFakeRateTriggerSkim.root","Ele10_Pt10To20", 1, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Ele10Jet30_Pt10To20", 10, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt10To20", 10, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt10To20", 100, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt10To20", 100, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt10To20", 200, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt10To20", 250, 10, 20);
  
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt20To35", 0, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Ele10Jet30_Pt20To35", 10, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt20To35", 10, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To35", 100, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt20To35", 100, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt20To35", 200, 20, 35);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt20To35", 250, 20, 35);

//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt35To50", 0, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd30_50-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD30To50MC_Ele10Jet30_Pt35To50", 10, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt35To50", 10, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt35To50", 100, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt35To50", 100, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root","WJetsMC_Pt35To50", 200, 35, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_WJetsMCCombined.root","WJetsMCCombined_Pt35To50", 250, 35, 50);
  

//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt10To20", 100, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt10To20", 0, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt10To20", 100, 10, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt10To20", 10, 10, 20);

  MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt15To20", 100, 15, 20);
  MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt15To20", 0, 15, 20);
  MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt15To20", 100, 15, 20);
  MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt15To20", 10, 15, 20);

//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To30", 100, 20, 30);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","Ele10Jet30_Pt20To30", 0, 20, 30);  
//    MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root","Zee_Pt20To30", 100, 20, 30);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-qcd15_3000-v1g1-pu-skimmed_RecoLeptonPlusJetSkim_normalized.root","QCD15To3000MC_Ele10Jet30_Pt20To30", 10, 20, 30);





//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_PtSelected", 100, 15, 2000);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt15To20", 100, 15, 20);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20To30", 100, 20, 30);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt30To50", 100, 30, 50);
//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root","HWW130_Pt20", 100, 20, 2000);
  
}



void MakeBkgElectron(const string inputFilename,
                     const string Label, Int_t SelectionType, Double_t PtMin, Double_t PtMax) 
{  

  string label = Label; 
  
  if (Label != "") label = "_"+Label;
  
  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Double_t jetPtThreshold = 30;

  //*****************************************************************************************
  //Histogram
  //*****************************************************************************************
  vector<Double_t> ElectronIsolationEfficiencyNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyNumerator;
  vector<Double_t> ElectronIsolationEfficiencyDenominator;
  for(UInt_t i=0 ; i<21; ++i) {
    ElectronIsolationEfficiencyNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyNumerator.push_back(0.0);
    ElectronIsolationEfficiencyDenominator.push_back(0.0);
  }

//   TH1F *Electron_relIso_Barrel = new TH1F((string("Electron_relIso_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIsoL1Corrected_Barrel = new TH1F((string("Electron_relIsoL1Corrected_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIso04_Barrel = new TH1F((string("Electron_relIso04_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIso04L1Corrected_Barrel = new TH1F((string("Electron_relIso04L1Corrected_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIso_Endcap = new TH1F((string("Electron_relIso_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIsoL1Corrected_Endcap = new TH1F((string("Electron_relIsoL1Corrected_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIso04_Endcap = new TH1F((string("Electron_relIso04_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//   TH1F *Electron_relIso04L1Corrected_Endcap = new TH1F((string("Electron_relIso04L1Corrected_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);


  vector<TH1F*>   Rho;
  vector<TH1F*>   Electron_caloIso_Barrel;
  vector<TH1F*>   Electron_caloIso_Endcap;
  vector<TH1F*>   Electron_ecalIso_Barrel;
  vector<TH1F*>   Electron_ecalIso_Endcap;
  vector<TH1F*>   Electron_hcalIso_Barrel;
  vector<TH1F*>   Electron_hcalIso_Endcap;
  vector<TH1F*>   Electron_trkIso_Barrel;
  vector<TH1F*>   Electron_trkIso_Endcap;
  vector<TH1F*>   Electron_caloIso04_Barrel;
  vector<TH1F*>   Electron_caloIso04_Endcap;
  vector<TH1F*>   Electron_ecalIso04_Barrel;
  vector<TH1F*>   Electron_ecalIso04_Endcap;
  vector<TH1F*>   Electron_hcalIso04_Barrel;
  vector<TH1F*>   Electron_hcalIso04_Endcap;
  vector<TH1F*>   Electron_trkIso04_Barrel;
  vector<TH1F*>   Electron_trkIso04_Endcap;

  vector<TH1F*>   Electron_ChargedIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_GammaIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_ChargedIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralHadronIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_GammaIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap;




  vector<TH1F*>   Electron_RelIso_Barrel;
  vector<TH1F*>   Electron_RelIso_Endcap;
  vector<TH1F*>   Electron_RelIso_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_RelIso_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_RelIso04_Barrel;
  vector<TH1F*>   Electron_RelIso04_Endcap;
  vector<TH1F*>   Electron_RelIso04_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_RelIso04_RhoCorrected_Endcap;


  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;





  for (int n=0; n < 20; ++n) {     
    TH1F *tmpRho = new TH1F((string("RhoElectron_") + IntToString(n) + label).c_str(), "; Rho; Number of Events ", 2000, -20, 20);
    TH1F *tmpElectron_caloIso_Barrel = new TH1F((string("Electron_caloIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso_Endcap = new TH1F((string("Electron_caloIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ecalIso_Barrel = new TH1F((string("Electron_ecalIso_Barrel_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ecalIso_Endcap = new TH1F((string("Electron_ecalIso_Endcap_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_hcalIso_Barrel = new TH1F((string("Electron_hcalIso_Barrel_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_hcalIso_Endcap = new TH1F((string("Electron_hcalIso_Endcap_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso_Barrel = new TH1F((string("Electron_trkIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso_Endcap = new TH1F((string("Electron_trkIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso04_Barrel = new TH1F((string("Electron_caloIso04_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso04_Endcap = new TH1F((string("Electron_caloIso04_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ecalIso04_Barrel = new TH1F((string("Electron_ecalIso04_Barrel_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ecalIso04_Endcap = new TH1F((string("Electron_ecalIso04_Endcap_")+ IntToString(n) +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_hcalIso04_Barrel = new TH1F((string("Electron_hcalIso04_Barrel_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_hcalIso04_Endcap = new TH1F((string("Electron_hcalIso04_Endcap_")+ IntToString(n) +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso04_Barrel = new TH1F((string("Electron_trkIso04_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso04_Endcap = new TH1F((string("Electron_trkIso04_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);

    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; ChargedNoPUIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedTotalPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFIso; Fraction of Events ", 2000, -200, 200.0);









    TH1F *tmpElectron_RelIso_Barrel = new TH1F((string("Electron_RelIso_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso_RhoCorrected_Barrel = new TH1F((string("Electron_RelIso_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso_Endcap = new TH1F((string("Electron_RelIso_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso_RhoCorrected_Endcap = new TH1F((string("Electron_RelIso_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso04_Barrel = new TH1F((string("Electron_RelIso04_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso04_RhoCorrected_Barrel = new TH1F((string("Electron_RelIso04_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; RelIso04_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso04_Endcap = new TH1F((string("Electron_RelIso04_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso_Cone04_01Threshold; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_RelIso04_RhoCorrected_Endcap = new TH1F((string("Electron_RelIso04_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; RelIso04_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);



    TH1F *tmpElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; TotalPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_")+ IntToString(n) +label).c_str(), "; VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected; Fraction of Events ", 2000, -4, 4.0);






    tmpRho->StatOverflows(kTRUE);
    tmpElectron_caloIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_caloIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_ecalIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_ecalIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_hcalIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_hcalIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_trkIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_trkIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_caloIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_caloIso04_Endcap->StatOverflows(kTRUE);
    tmpElectron_ecalIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_ecalIso04_Endcap->StatOverflows(kTRUE);
    tmpElectron_hcalIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_hcalIso04_Endcap->StatOverflows(kTRUE);
    tmpElectron_trkIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_trkIso04_Endcap->StatOverflows(kTRUE);

    tmpElectron_ChargedIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_ChargedIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_GammaIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);





    tmpElectron_RelIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_RelIso_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_RelIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_RelIso_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_RelIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_RelIso04_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_RelIso04_Endcap->StatOverflows(kTRUE);
     tmpElectron_RelIso04_RhoCorrected_Endcap->StatOverflows(kTRUE);




    tmpElectron_TotalPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_01Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_05Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_10Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone03_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap->StatOverflows(kTRUE);
    tmpElectron_TotalPFRelIso_Cone04_15Threshold_Endcap->StatOverflows(kTRUE);




    Rho.push_back(tmpRho);
    Electron_caloIso_Barrel.push_back(tmpElectron_caloIso_Barrel);
    Electron_ecalIso_Barrel.push_back(tmpElectron_ecalIso_Barrel);
    Electron_hcalIso_Barrel.push_back(tmpElectron_hcalIso_Barrel);
    Electron_trkIso_Barrel.push_back(tmpElectron_trkIso_Barrel);
    Electron_caloIso04_Barrel.push_back(tmpElectron_caloIso04_Barrel);       
    Electron_ecalIso04_Barrel.push_back(tmpElectron_ecalIso04_Barrel);       
    Electron_hcalIso04_Barrel.push_back(tmpElectron_hcalIso04_Barrel);       
    Electron_trkIso04_Barrel.push_back(tmpElectron_trkIso04_Barrel);       
    Electron_caloIso_Endcap.push_back(tmpElectron_caloIso_Endcap);
    Electron_ecalIso_Endcap.push_back(tmpElectron_ecalIso_Endcap);
    Electron_hcalIso_Endcap.push_back(tmpElectron_hcalIso_Endcap);
    Electron_trkIso_Endcap.push_back(tmpElectron_trkIso_Endcap);
    Electron_caloIso04_Endcap.push_back(tmpElectron_caloIso04_Endcap);       
    Electron_ecalIso04_Endcap.push_back(tmpElectron_ecalIso04_Endcap);       
    Electron_hcalIso04_Endcap.push_back(tmpElectron_hcalIso04_Endcap);       
    Electron_trkIso04_Endcap.push_back(tmpElectron_trkIso04_Endcap);

    Electron_ChargedIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    Electron_NeutralHadronIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel);
    Electron_GammaIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Barrel);
    Electron_TotalPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel);
    Electron_ChargedIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    Electron_NeutralHadronIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap);
    Electron_GammaIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Endcap);
    Electron_TotalPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap);
    Electron_ChargedIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    Electron_NeutralHadronIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel);
    Electron_GammaIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Barrel);
    Electron_TotalPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel);
    Electron_ChargedIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    Electron_NeutralHadronIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap);
    Electron_GammaIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Endcap);
    Electron_TotalPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap);
    Electron_ChargedIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    Electron_NeutralHadronIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel);
    Electron_GammaIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Barrel);
    Electron_TotalPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel);
    Electron_ChargedIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    Electron_NeutralHadronIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap);
    Electron_GammaIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Endcap);
    Electron_TotalPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap);
    Electron_ChargedIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    Electron_NeutralHadronIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel);
    Electron_GammaIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Barrel);
    Electron_TotalPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel);
    Electron_ChargedIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    Electron_NeutralHadronIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap);
    Electron_GammaIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Endcap);
    Electron_TotalPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap);
    Electron_ChargedIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    Electron_NeutralHadronIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel);
    Electron_GammaIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Barrel);
    Electron_TotalPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel);
    Electron_ChargedIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    Electron_NeutralHadronIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap);
    Electron_GammaIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Endcap);
    Electron_TotalPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap);
    Electron_ChargedIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    Electron_NeutralHadronIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel);
    Electron_GammaIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Barrel);
    Electron_TotalPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel);
    Electron_ChargedIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    Electron_NeutralHadronIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap);
    Electron_GammaIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Endcap);
    Electron_TotalPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap);
    Electron_ChargedIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    Electron_NeutralHadronIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel);
    Electron_GammaIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Barrel);
    Electron_TotalPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel);
    Electron_ChargedIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    Electron_NeutralHadronIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap);
    Electron_GammaIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Endcap);
    Electron_TotalPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap);
    Electron_ChargedIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel);
    Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    Electron_NeutralHadronIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel);
    Electron_GammaIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Barrel);
    Electron_TotalPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel);
    Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel);
    Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel);
    Electron_ChargedIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap);
    Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    Electron_NeutralHadronIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap);
    Electron_GammaIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Endcap);
    Electron_TotalPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap);
    Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap);
    Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap);








    Electron_RelIso_Barrel.push_back(tmpElectron_RelIso_Barrel);
    Electron_RelIso_RhoCorrected_Barrel.push_back(tmpElectron_RelIso_RhoCorrected_Barrel);
    Electron_RelIso_Endcap.push_back(tmpElectron_RelIso_Endcap);
    Electron_RelIso_RhoCorrected_Endcap.push_back(tmpElectron_RelIso_RhoCorrected_Endcap);
    Electron_RelIso04_Barrel.push_back(tmpElectron_RelIso04_Barrel);
    Electron_RelIso04_RhoCorrected_Barrel.push_back(tmpElectron_RelIso04_RhoCorrected_Barrel);
    Electron_RelIso04_Endcap.push_back(tmpElectron_RelIso04_Endcap);
    Electron_RelIso04_RhoCorrected_Endcap.push_back(tmpElectron_RelIso04_RhoCorrected_Endcap);


    Electron_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_01Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel);
    Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_01Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap);
    Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_01Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel);
    Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_01Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap);
    Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_05Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel);
    Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_05Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap);
    Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_05Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel);
    Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_05Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap);
    Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_10Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel);
    Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_10Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap);
    Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_10Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel);
    Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_10Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap);
    Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_15Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel);
    Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_15Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap);
    Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    Electron_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_15Threshold_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel);
    Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    Electron_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_15Threshold_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap);
    Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);




  }

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
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
   rlrm.AddJSONFile("Cert_160404-163757_7TeV_PromptReco_Collisions11_JSON.txt"); 
//    hasJSON = kFALSE;

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
    if(SelectionType < 10 && hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
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
	
    Int_t NVertex = info->nPV0; if (NVertex > 19) NVertex = 19;

    //********************************************************
    // Event Selection Cuts
    //********************************************************
    if (SelectionType < 100) {
      if (SelectionType < 10) {
        if (!(
              (info->triggerBits & kHLT_Ele8)
              || 
              (info->triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
              || 
              (info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
//            || 
//           (info->triggerBits & kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL)
              )
          ) continue;
      }

      if (met.Pt() > 20) continue;
      if (electronArr->GetEntries() > 1) continue;
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


    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

      //pt cut
      if(!(ele->pt > PtMin && ele->pt <= PtMax )) continue;

 //      //pass HLT selection
//       if (!passHLT(info->triggerBits, info->runNum, SelectionType)) continue;
      
      //pass event selection
      Bool_t passJetSelection = kFALSE;
      for(Int_t j=0; j<jetArr->GetEntries(); j++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);
        
        if (jet->pt > jetPtThreshold &&
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 0.5) {
          passJetSelection = kTRUE;
          break;
        }
      }

      if (SelectionType < 100) {
        if (SelectionType != 1) {
          if (!passJetSelection) continue;
        }
      }

      //WJets MC - pick only fake muons
      Int_t originalvalue = ele->isMCReal;
      Int_t isEle;
      Int_t isMu;
      Int_t isTau;
      
      isEle = originalvalue % 2;
      Int_t tmp0 = floor(double(originalvalue) / 2.0);
      isMu = tmp0 % 2;
      Int_t tmp1 = floor(double(tmp0) / 2.0);
      isTau = tmp1 % 2;
      
      if (SelectionType == 100) {
        if (!isEle) continue;
//         if (!passDileptonPreselection) continue;
      }
      
      if (SelectionType == 200) {
//          if (isEle || isMu || isTau) {
        if (!(ele->isMCReal == 0)) {
          continue;
        }
      }

      if (!passElectronDenominatorCuts(ele)) continue;


//       Double_t relIso = (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt;
//       if (fabs(ele->eta) > 1.5) relIso = (ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt;
//       Double_t relIso_Cone04_01Threshold = (ele->trkIso04 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt;
//       if (fabs(ele->eta) > 1.5) relIso_Cone04_01Threshold = (ele->trkIso04 + ele->emIso04 + ele->hadIso04) / ele->pt;
//       Double_t relIsoL1Corrected = (ele->trkIso03 + TMath::Max(TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy, 0.0)) / ele->pt;
//       if (fabs(ele->eta) > 1.5) relIso = (ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt;
//       Double_t relIso_Cone04_01ThresholdL1Corrected = (ele->trkIso04 + TMath::Max(TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 - PUIsolationEnergy04Cone,0.0)) / ele->pt;
//       if (fabs(ele->eta) > 1.5) relIso_Cone04_01ThresholdL1Corrected = (ele->trkIso04 + ele->emIso04 + ele->hadIso04) / ele->pt;

      Rho[NVertex]->Fill(info->PileupEnergyDensity);
      if (fabs(ele->eta) < 1.5) {
        Electron_caloIso_Barrel[NVertex]->Fill(TMath::Min(double( TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03), double(199.99)));
        Electron_ecalIso_Barrel[NVertex]->Fill(TMath::Min(double( TMath::Max(ele->emIso03 - 1.0, 0.0)), double(199.99)));
        Electron_hcalIso_Barrel[NVertex]->Fill(TMath::Min(double( ele->hadIso03 ), double(199.99)));
        Electron_trkIso_Barrel[NVertex]->Fill(TMath::Min(double( ele->trkIso03 ), double(199.99)));
        Electron_caloIso04_Barrel[NVertex]->Fill(TMath::Min( double(TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04), double(199.99)));
        Electron_ecalIso04_Barrel[NVertex]->Fill(TMath::Min( double(TMath::Max(ele->emIso04 - 1.0, 0.0) ), double(199.99)));
        Electron_hcalIso04_Barrel[NVertex]->Fill(TMath::Min( double(ele->hadIso04), double(199.99)));
        Electron_trkIso04_Barrel[NVertex]->Fill(TMath::Min( double(ele->trkIso04 ) , double(199.99)));

        Electron_ChargedIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_01Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_01Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_05Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_05Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_10Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_10Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_15Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_15Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));





  
        Electron_RelIso_Barrel[NVertex]->Fill(TMath::Min( double((TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 + ele->trkIso03)/ele->pt) , double(3.99)));
        Electron_RelIso_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 + ele->trkIso03 - info->PileupEnergyDensity*0.0698)/ele->pt) , double(3.99)));
        Electron_RelIso04_Barrel[NVertex]->Fill(TMath::Min( double((TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 + ele->trkIso04)/ele->pt) , double(3.99)));
        Electron_RelIso04_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 + ele->trkIso04 - info->PileupEnergyDensity*0.1917)/ele->pt) , double(3.99)));



        Electron_TotalPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold - info->PileupEnergyDensity*0.2759)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.2653)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.2647)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold - info->PileupEnergyDensity*0.0659)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.0558)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.0551)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold - info->PileupEnergyDensity*0.5058)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.4909)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.4902)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold - info->PileupEnergyDensity*0.1319)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.1165)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.1164)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold - info->PileupEnergyDensity*0.2574)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.2475)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_05Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.2468)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold - info->PileupEnergyDensity*0.0475)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0382)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0373)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold - info->PileupEnergyDensity*0.4703)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.4567)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.4561)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_05Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold - info->PileupEnergyDensity*0.0964)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0825)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0822)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold - info->PileupEnergyDensity*0.2310)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.2220)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->ChargedIso03FromOtherVertices_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.2213)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold - info->PileupEnergyDensity*0.0211)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0126)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_10Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0118)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold - info->PileupEnergyDensity*0.4201)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.4083)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.4077)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold - info->PileupEnergyDensity*0.0462)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0340)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0338)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold - info->PileupEnergyDensity*0.2186)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.2109)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->ChargedIso03FromOtherVertices_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.2103)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold - info->PileupEnergyDensity*0.0087)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0012)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_15Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0007)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold - info->PileupEnergyDensity*0.3974)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.3874)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.3868)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold - info->PileupEnergyDensity*0.0235)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0128)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0128)/ele->pt) , double(3.99)));









//         Electron_relIso_Barrel->Fill( TMath::Min(double(relIso), double(1.99)) );
//         Electron_relIsoL1Corrected_Barrel->Fill(TMath::Min(double(relIsoL1Corrected), double(1.99)));
//         Electron_relIso04_Barrel->Fill(TMath::Min(double(relIso04), double(1.99)) );
//         Electron_relIso04L1Corrected_Barrel->Fill(TMath::Min(double(relIso04L1Corrected), double(1.99)));



      } else {
        Electron_caloIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->emIso03  + ele->hadIso03), double(199.99)));
        Electron_ecalIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->emIso03 ), double(199.99)));
        Electron_hcalIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->hadIso03 ), double(199.99)));
        Electron_trkIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->trkIso03 ), double(199.99)));
        Electron_caloIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->emIso04  + ele->hadIso04), double(199.99)));
        Electron_ecalIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->emIso04 ), double(199.99)));
        Electron_hcalIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->hadIso04 ), double(199.99)));
        Electron_trkIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->trkIso04 ) , double(199.99)));




        Electron_ChargedIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_01Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_01Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_05Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_05Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_10Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_10Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold) , double(199.99)));
        Electron_ChargedIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso03_15Threshold) , double(199.99)));
        Electron_GammaIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_ChargedIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold ) , double(199.99)));
        Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold ) , double(199.99)));
        Electron_NeutralHadronIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->NeutralHadronIso04_15Threshold) , double(199.99)));
        Electron_GammaIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_TotalPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold) , double(199.99)));
        Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));
        Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double(ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold) , double(199.99)));









        Electron_RelIso_Endcap[NVertex]->Fill(TMath::Min( double((ele->emIso03 + ele->hadIso03 + ele->trkIso03)/ele->pt) , double(3.99)));
        Electron_RelIso_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->emIso03 + ele->hadIso03 + ele->trkIso03 - info->PileupEnergyDensity*0.1001)/ele->pt) , double(3.99)));
        Electron_RelIso04_Endcap[NVertex]->Fill(TMath::Min( double((ele->emIso04 + ele->hadIso04 + ele->trkIso04)/ele->pt) , double(3.99)));
        Electron_RelIso04_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->emIso04 + ele->hadIso04 + ele->trkIso04 - info->PileupEnergyDensity*0.1901)/ele->pt) , double(3.99)));





        Electron_TotalPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold - info->PileupEnergyDensity*0.3562)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.2513)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.2487)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold - info->PileupEnergyDensity*0.1637)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.0573)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_01Threshold+ele->GammaIso03_01Threshold-ele->GammaIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.0573)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold - info->PileupEnergyDensity*0.5447)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.4390)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.4367)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold - info->PileupEnergyDensity*0.2184)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.1104)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_01Threshold+ele->GammaIso04_01Threshold-ele->GammaIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_01Threshold - info->PileupEnergyDensity*0.1108)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold - info->PileupEnergyDensity*0.3493)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.2447)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.2421)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold - info->PileupEnergyDensity*0.1568)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0510)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_05Threshold+ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0509)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold - info->PileupEnergyDensity*0.5297)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.4247)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.4224)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold - info->PileupEnergyDensity*0.2036)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0967)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_05Threshold+ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold-ele->NeutralHadronIso007_05Threshold - info->PileupEnergyDensity*0.0970)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold - info->PileupEnergyDensity*0.3311)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.2266)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.2241)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold - info->PileupEnergyDensity*0.1386)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0325)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_10Threshold+ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0329)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold - info->PileupEnergyDensity*0.4925)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.3879)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.3856)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold - info->PileupEnergyDensity*0.1665)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0596)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold - info->PileupEnergyDensity*0.0604)/ele->pt) , double(3.99)));




        Electron_TotalPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold - info->PileupEnergyDensity*0.3158)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.2147)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->ChargedIso03FromOtherVertices_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.2121)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold - info->PileupEnergyDensity*0.1233)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->ChargedEMIsoVetoEtaStrip03_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0203)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso03_01Threshold+ele->NeutralHadronIso03_15Threshold+ele->GammaIso03_15Threshold-ele->GammaIsoVetoEtaStrip03_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0208)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold)/ele->pt) , double(3.99)));
        Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold - info->PileupEnergyDensity*0.4663)/ele->pt) , double(3.99)));
        Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.3652)/ele->pt) , double(3.99)));
        Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->ChargedIso04FromOtherVertices_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.3629)/ele->pt) , double(3.99)));
        Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold - info->PileupEnergyDensity*0.1402)/ele->pt) , double(3.99)));
        Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->ChargedEMIsoVetoEtaStrip04_01Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0365)/ele->pt) , double(3.99)));
        Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[NVertex]->Fill(TMath::Min( double((ele->ChargedIso04_01Threshold+ele->NeutralHadronIso04_15Threshold+ele->GammaIso04_15Threshold-ele->GammaIsoVetoEtaStrip04_15Threshold-ele->NeutralHadronIso007_15Threshold - info->PileupEnergyDensity*0.0374)/ele->pt) , double(3.99)));
 
 



//         Electron_relIso_Endcap->Fill( TMath::Min(double(relIso), double(1.99)) );
//         Electron_relIsoL1Corrected_Endcap->Fill(TMath::Min(double(relIsoL1Corrected), double(1.99)));
//         Electron_relIso04_Endcap->Fill(TMath::Min(double(relIso04), double(1.99)) );
//         Electron_relIso04L1Corrected_Endcap->Fill(TMath::Min(double(relIso04L1Corrected), double(1.99)));
        
      }


      ElectronIsolationEfficiencyDenominator[NVertex]++;  
      ElectronIsolationEfficiencyDenominator[20]++;  
      
      if (fabs(ele->eta) < 1.5) {
        if ((ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.10) {
          ElectronIsolationEfficiencyNumerator[NVertex]++;
          ElectronIsolationEfficiencyNumerator[20]++;
        }
//         if ((ele->trkIso03 + TMath::Max( TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt < 0.10) {
//           ElectronL1CorrectedIsolationEfficiencyNumerator[NVertex]++;
//         }
      } else {
        if ((ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt < 0.10) {
          ElectronIsolationEfficiencyNumerator[NVertex]++;
          ElectronIsolationEfficiencyNumerator[20]++;
        }
//         if ( (ele->trkIso03 + TMath::Max(ele->emIso03 + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt < 0.10 ) {
//           ElectronL1CorrectedIsolationEfficiencyNumerator[NVertex]++;
//         }
      }
      
    }

  } //end loop over data     


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  
  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  const int nPoints = 20;
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.0;     
  }

  double ElectronIsolationEff[nPoints];
  double ElectronIsolationEffErrLow[nPoints];
  double ElectronIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationEff[i] = ratio;
    ElectronIsolationEffErrLow[i] = errLow;
    ElectronIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronIsolationEff " << i << " : " << ElectronIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationEff, NPileupError, NPileupError, ElectronIsolationEffErrLow, ElectronIsolationEffErrHigh);
  ElectronIsolationEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  Double_t r = 0;
  Double_t eLow = 0;
  Double_t eHigh = 0;     
  Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyNumerator[20]);
  Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyDenominator[20]);
  mithep::MathUtils::CalcRatio(n1 , n2, r, eLow, eHigh, 2);
  cout << "ElectronIsolationEff Averaged : " << r << " + " << eHigh << " - " << eLow << endl;
  



//   double ElectronL1CorrectedIsolationEff[nPoints];
//   double ElectronL1CorrectedIsolationEffErrLow[nPoints];
//   double ElectronL1CorrectedIsolationEffErrHigh[nPoints];
//   for (UInt_t i=0; i<nPoints; ++i) {
    
//     Double_t ratio;
//     Double_t errLow;
//     Double_t errHigh;     
    
//     Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyNumerator[i]);
//     Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyDenominator[i]);
//     mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
//     ElectronL1CorrectedIsolationEff[i] = ratio;
//     ElectronL1CorrectedIsolationEffErrLow[i] = errLow;
//     ElectronL1CorrectedIsolationEffErrHigh[i] = errHigh;
//     NPileup[i] = i;
//     cout << "ElectronL1CorrectedIsolationEff " << i << " : " << ElectronL1CorrectedIsolationEff[i] << endl;
//   }
//   TGraphAsymmErrors *ElectronL1CorrectedIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationEff, NPileupError, NPileupError,ElectronL1CorrectedIsolationEffErrLow,ElectronL1CorrectedIsolationEffErrHigh  );
//   ElectronL1CorrectedIsolationEffVsNPileup->SetMarkerColor(kBlue);
//   ElectronL1CorrectedIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
//   ElectronL1CorrectedIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  Double_t ymin = 0.0;
  Double_t ymax = 1.1;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.0;
  legend->AddEntry(ElectronIsolationEffVsNPileup, "NoCorrection", "LP");
//   legend->AddEntry(ElectronL1CorrectedIsolationEffVsNPileup, "FastJetCorrected", "LP");
  ElectronIsolationEffVsNPileup->SetTitle("");
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronIsolationEffVsNPileup->Draw("AP");
//   ElectronL1CorrectedIsolationEffVsNPileup->Draw("Psame");
  ElectronIsolationEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);

  legend->Draw();
  cv->SaveAs(("ElectronIsolationEfficiency_BkgData_vs_NVertices"+label+".gif").c_str());




  //*****************************************************************************************
  //Save Efficiency Plots
  //*****************************************************************************************
  TFile *file = new TFile("HwwSelectionPlots_LeptonEfficiency.electrons.root", "UPDATE");
  file->cd();
  file->WriteTObject(ElectronIsolationEffVsNPileup, ("ElectronIsolationEffVsNVertices" +label).c_str() , "WriteDelete");
//   file->WriteTObject(ElectronL1CorrectedIsolationEffVsNPileup, ("ElectronFastjetCorrectedIsolationEffVsNVertices" +label).c_str() , "WriteDelete");

//   file->WriteTObject(Electron_relIso_Cone03_01Threshold_Barrel, Electron_relIso_Cone03_01Threshold_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIsoL1Corrected_Barrel, Electron_relIsoL1Corrected_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIso_Cone04_01Threshold_Barrel, Electron_relIso_Cone04_01Threshold_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIso_Cone04_01ThresholdL1Corrected_Barrel, Electron_relIso_Cone04_01ThresholdL1Corrected_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIso_Cone03_01Threshold_Endcap, Electron_relIso_Cone03_01Threshold_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIsoL1Corrected_Endcap, Electron_relIsoL1Corrected_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIso_Cone04_01Threshold_Endcap, Electron_relIso_Cone04_01Threshold_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(Electron_relIso_Cone04_01ThresholdL1Corrected_Endcap, Electron_relIso_Cone04_01ThresholdL1Corrected_Endcap->GetName(), "WriteDelete");

   for (int n=0; n < 20 ; ++n) {
      file->WriteTObject(Rho[n],Rho[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Barrel[n],Electron_caloIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso_Barrel[n],Electron_ecalIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso_Barrel[n],Electron_hcalIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Barrel[n],Electron_trkIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Barrel[n],Electron_caloIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso04_Barrel[n],Electron_ecalIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso04_Barrel[n],Electron_hcalIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Barrel[n],Electron_trkIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Endcap[n],Electron_caloIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso_Endcap[n],Electron_ecalIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso_Endcap[n],Electron_hcalIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Endcap[n],Electron_trkIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Endcap[n],Electron_caloIso04_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso04_Endcap[n],Electron_ecalIso04_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso04_Endcap[n],Electron_hcalIso04_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Endcap[n],Electron_trkIso04_Endcap[n]->GetName(), "WriteDelete") ;

      file->WriteTObject(Electron_ChargedIso_Cone03_01Threshold_Barrel[n],Electron_ChargedIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_01Threshold_Barrel[n],Electron_NeutralHadronIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_01Threshold_Barrel[n],Electron_GammaIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_01Threshold_Barrel[n],Electron_TotalPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone03_01Threshold_Endcap[n],Electron_ChargedIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_01Threshold_Endcap[n],Electron_NeutralHadronIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_01Threshold_Endcap[n],Electron_GammaIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_01Threshold_Endcap[n],Electron_TotalPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_01Threshold_Barrel[n],Electron_ChargedIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_01Threshold_Barrel[n],Electron_NeutralHadronIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_01Threshold_Barrel[n],Electron_GammaIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_01Threshold_Barrel[n],Electron_TotalPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_01Threshold_Endcap[n],Electron_ChargedIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_01Threshold_Endcap[n],Electron_NeutralHadronIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_01Threshold_Endcap[n],Electron_GammaIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_01Threshold_Endcap[n],Electron_TotalPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete");


      file->WriteTObject(Electron_ChargedIso_Cone03_05Threshold_Barrel[n],Electron_ChargedIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_05Threshold_Barrel[n],Electron_NeutralHadronIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_05Threshold_Barrel[n],Electron_GammaIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_05Threshold_Barrel[n],Electron_TotalPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone03_05Threshold_Endcap[n],Electron_ChargedIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_05Threshold_Endcap[n],Electron_NeutralHadronIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_05Threshold_Endcap[n],Electron_GammaIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_05Threshold_Endcap[n],Electron_TotalPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_05Threshold_Barrel[n],Electron_ChargedIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_05Threshold_Barrel[n],Electron_NeutralHadronIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_05Threshold_Barrel[n],Electron_GammaIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_05Threshold_Barrel[n],Electron_TotalPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_05Threshold_Endcap[n],Electron_ChargedIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_05Threshold_Endcap[n],Electron_NeutralHadronIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_05Threshold_Endcap[n],Electron_GammaIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_05Threshold_Endcap[n],Electron_TotalPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete");



      file->WriteTObject(Electron_ChargedIso_Cone03_10Threshold_Barrel[n],Electron_ChargedIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_10Threshold_Barrel[n],Electron_NeutralHadronIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_10Threshold_Barrel[n],Electron_GammaIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_10Threshold_Barrel[n],Electron_TotalPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone03_10Threshold_Endcap[n],Electron_ChargedIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_10Threshold_Endcap[n],Electron_NeutralHadronIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_10Threshold_Endcap[n],Electron_GammaIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_10Threshold_Endcap[n],Electron_TotalPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_10Threshold_Barrel[n],Electron_ChargedIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_10Threshold_Barrel[n],Electron_NeutralHadronIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_10Threshold_Barrel[n],Electron_GammaIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_10Threshold_Barrel[n],Electron_TotalPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_10Threshold_Endcap[n],Electron_ChargedIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_10Threshold_Endcap[n],Electron_NeutralHadronIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_10Threshold_Endcap[n],Electron_GammaIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_10Threshold_Endcap[n],Electron_TotalPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete");



      file->WriteTObject(Electron_ChargedIso_Cone03_15Threshold_Barrel[n],Electron_ChargedIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_15Threshold_Barrel[n],Electron_NeutralHadronIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_15Threshold_Barrel[n],Electron_GammaIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_15Threshold_Barrel[n],Electron_TotalPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone03_15Threshold_Endcap[n],Electron_ChargedIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone03_15Threshold_Endcap[n],Electron_NeutralHadronIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone03_15Threshold_Endcap[n],Electron_GammaIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone03_15Threshold_Endcap[n],Electron_TotalPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_15Threshold_Barrel[n],Electron_ChargedIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel[n],Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_15Threshold_Barrel[n],Electron_NeutralHadronIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_15Threshold_Barrel[n],Electron_GammaIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_15Threshold_Barrel[n],Electron_TotalPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel[n],Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[n],Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedIso_Cone04_15Threshold_Endcap[n],Electron_ChargedIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap[n],Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralHadronIso_Cone04_15Threshold_Endcap[n],Electron_NeutralHadronIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_GammaIso_Cone04_15Threshold_Endcap[n],Electron_GammaIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_TotalPFIso_Cone04_15Threshold_Endcap[n],Electron_TotalPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap[n],Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[n],Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete");




      file->WriteTObject(Electron_RelIso_Barrel[n],Electron_RelIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso_RhoCorrected_Barrel[n],Electron_RelIso_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso_Endcap[n],Electron_RelIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso_RhoCorrected_Endcap[n],Electron_RelIso_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso04_Barrel[n],Electron_RelIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso04_RhoCorrected_Barrel[n],Electron_RelIso04_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso04_Endcap[n],Electron_RelIso04_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso04_RhoCorrected_Endcap[n],Electron_RelIso04_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;


      file->WriteTObject(Electron_TotalPFRelIso_Cone03_01Threshold_Barrel[n],Electron_TotalPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_01Threshold_Endcap[n],Electron_TotalPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_01Threshold_Barrel[n],Electron_TotalPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_01Threshold_Endcap[n],Electron_TotalPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_05Threshold_Barrel[n],Electron_TotalPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_05Threshold_Endcap[n],Electron_TotalPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_05Threshold_Barrel[n],Electron_TotalPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_05Threshold_Endcap[n],Electron_TotalPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_10Threshold_Barrel[n],Electron_TotalPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_10Threshold_Endcap[n],Electron_TotalPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_10Threshold_Barrel[n],Electron_TotalPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_10Threshold_Endcap[n],Electron_TotalPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_15Threshold_Barrel[n],Electron_TotalPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_15Threshold_Endcap[n],Electron_TotalPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_15Threshold_Barrel[n],Electron_TotalPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_15Threshold_Endcap[n],Electron_TotalPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_FPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_NeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n],Electron_VertexSelectedNeutralFPRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[n]->GetName(), "WriteDelete") ;




    }


  file->Close();
  delete file;

    
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType) {


  Bool_t pass = kFALSE;


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
Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }

  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.025
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }
  } else {
    pass = kFALSE;
    return pass;
  }

  if (ele->pt < 20) {
    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 )) pass = kFALSE;
      }
    }
  }

  return pass;
}


Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }

  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
//             && ele->HoverE < 0.04
//             && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
//              && ele->HoverE < 0.025
//              && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }
  } else {
    pass = kFALSE;
    return pass;
  }

  
  return pass;
}
