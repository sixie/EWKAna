//root -l -b -q EWKAna/Hww/Selection/WWFakePrediction_MC.C+\(\"\"\)
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
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "THStack.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"

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
#include "MitAna/Utils/interface/SimpleTable.h"
#include "MitPhysics/FakeMods/interface/FakeRate.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"

#endif

#define LUMINOSITY 35.0

//=== FUNCTION DECLARATIONS ======================================================================================

Double_t getNormalizationWeight(string filename, string datasetName);
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DrawPrediction(TH1F* histPrediction,TH1F* histSimulation, string histname,
                    Double_t minY , Double_t maxY,
                    Bool_t useLogY, TFile *saveFile);
string DoubleToString(double i);

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
      cout << "Cannot get Directory ZeeAnalysisMod from file " << infname << endl;
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





//=== MAIN MACRO =================================================================================================

void WWFakePrediction_MC(const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  vector<vector<string> > processInputFiles;
   processInputFiles.push_back(vector<string>());
   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w0jets-0-100-v26_noskim_normalized.root");
//    processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w1jets-0-100-v26_noskim_normalized.root");
// //     processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w2jets-0-100-v26_noskim_normalized.root");
// //   processInputFiles.push_back(vector<string>());
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w3jets-0-100-v26_noskim_normalized.root");
// //     processInputFiles.push_back(vector<string>());
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w4jets-0-100-v26_noskim_normalized.root");
// //     processInputFiles.push_back(vector<string>());
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-w5jets-0-100-v26_noskim_normalized.root");
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-wc0jets-v26_noskim_normalized.root");
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-wcc0jets-v26_noskim_normalized.root");
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-wbb0jets-v26_noskim_normalized.root");

//     processInputFiles.push_back(vector<string>());
//       processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-wjets-mg-v26_noskim_normalized.root");

//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ww2l-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ttbar2l-v26_noskim.root");
//     processInputFiles.push_back(vector<string>());
//     processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-wjets-mg-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-wtop-mg-v26_noskim.root");
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-stop-mg-v26_noskim.root");
//   //processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ttop-mg-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-wz3l-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-zz4l-v26_noskim.root");
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-zz2l-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ztt-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-zeem20-c66-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-zmmm20-c66-v26_noskim.root");


  vector<vector<string> > processNames;
   processNames.push_back(vector<string>());
    processNames.back().push_back("p10-w0jets-0-100-v26");
//     processNames.back().push_back("p10-w1jets-0-100-v26");
// //    processNames.push_back(vector<string>());
//     processNames.back().push_back("p10-w2jets-0-100-v26");
// //   processNames.push_back(vector<string>());
//     processNames.back().push_back("p10-w3jets-0-100-v26");
// //   processNames.push_back(vector<string>());
//     processNames.back().push_back("p10-w4jets-0-100-v26");
// //   processNames.push_back(vector<string>());
//     processNames.back().push_back("p10-w5jets-0-100-v26");
//    processNames.back().push_back("p10-wc0jets-v26");
//    processNames.back().push_back("p10-wcc0jets-v26");
//    processNames.back().push_back("p10-wbb0jets-v26");

//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-ww2l-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-ttbar2l-v26");
//     processNames.push_back(vector<string>());
//     processNames.back().push_back("p10-wjets-mg-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-wtop-mg-v26");
//   processNames.back().push_back("p10-stop-mg-v26");
//   //processNames.back().push_back("p10-ttop-mg-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-wz3l-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-zz4l-v26");
//   processNames.back().push_back("p10-zz2l-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-ztt-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-zeem20-c66-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-zmmm20-c66-v26");

  vector<string> processLabels;
//   processLabels.push_back("WW");
//   processLabels.push_back("ttbar");
   processLabels.push_back("W1Jets");
//   processLabels.push_back("W2Jets");
//   processLabels.push_back("W3Jets");
//   processLabels.push_back("W4Jets");
//   processLabels.push_back("W5Jets");
//   processLabels.push_back("SingleTop");
//   processLabels.push_back("WZ");
//   processLabels.push_back("ZZ");
//   processLabels.push_back("Ztt");
//   processLabels.push_back("Zee");  
//   processLabels.push_back("Zmm");
  vector<Int_t> processColors;  
//   processColors.push_back(kBlue);
//   processColors.push_back(kRed);
    processColors.push_back(kGreen);
//    processColors.push_back(kGreen+1);
//    processColors.push_back(kGreen+2);
//    processColors.push_back(kGreen+3);
//    processColors.push_back(kGreen+4);
//    processColors.push_back(kCyan);
//   processColors.push_back(kMagenta);
//   processColors.push_back(kYellow);
//   processColors.push_back(kGray);
//   processColors.push_back(kAzure+4);
//   processColors.push_back(kRed+4);

  assert(processNames.size() == processInputFiles.size());
  assert(processNames.size() == processLabels.size());
  assert(processNames.size() == processColors.size());

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
  Int_t ChargeSelection = 0;
  //Int_t ChargeSelection = 1;

  //--------------------------------------------------------------------------------------------------------------
  // Set up Fake Rate
  //==============================================================================================================
  Bool_t use2DFakeRate = kTRUE;
  Bool_t useFitFunction = kFALSE;
//   mithep::FakeRate *fFakeRate = new mithep::FakeRate("MitPhysics/data/ElectronFakeRates.20101116.root",
//                            "MitPhysics/data/FakeRates.20100323.root",
//                            "", "",
//                            "RecoElectronEfficiency_RecoDenominator_WWTightIdIsoNumerator_WithSystematics_MetCut_data_Jet30_PtEta",
//                            "TrackerMuonFakeRate_PtEta_Pythia_Jet30",
//                            use2DFakeRate, useFitFunction );

  mithep::FakeRate *fMCFakeRate = new mithep::FakeRate("MitPhysics/data/ElectronFakeRates.20101116.root",
                                                       "MitPhysics/data/FakeRates.20100323.root",
                                                       "", "",
//                                                         "RecoElectronEfficiency_RecoDenominator_WWTightIdIsoNumerator_WMCMadgraph_PtEta",
//                                                         "RecoElectronEfficiency_RecoDenominator_WWTightIdIsoNumerator_QCD20MCv26a_Jet50_PtEta",
//                                                        "RecoElectronEfficiency_IsoDenominator_WWTightIdIsoNumerator_WMCMadgraph_PtEta",
                                                         "RecoElectronEfficiency_IsoDenominator_WWTightIdIsoNumerator_QCDMC_Jet30_PtEta",
                                                       "TrackerMuonFakeRate_PtEta_Pythia_Jet30",
                                                       use2DFakeRate, useFitFunction );
  

  //--------------------------------------------------------------------------------------------------------------
  //Prediction Histograms
  //==============================================================================================================  
  THStack *fHWWSelectionStack = new THStack("HWWSelectionStack", "HWWSelectionStack");
  THStack *fHWWToEESelectionStack = new THStack("HWWToEESelectionStack", "HWWToEESelectionStack");
  THStack *fHWWToMuMuSelectionStack = new THStack("HWWToMuMuSelectionStack", "HWWToMuMuSelectionStack");
  THStack *fHWWToEMuSelectionStack = new THStack("HWWToEMuSelectionStack", "HWWToEMuSelectionStack");

  THStack *fLeptonEtaStack = new THStack("LeptonEtaStack", "LeptonEtaStack");
  THStack *fLeptonPtMaxStack = new THStack("LeptonPtMaxStack", "LeptonPtMaxStack");
  THStack *fLeptonPtMinStack = new THStack("LeptonPtMinStack", "LeptonPtMinStack");
  THStack *fMetPtHistStack = new THStack("MetPtHistStack", "MetPtHistStack");
  THStack *fMetPhiHistStack = new THStack("MetPhiHistStack", "MetPhiHistStack");
  THStack *fDeltaPhiLeptonsStack = new THStack("DeltaPhiLeptonsStack", "DeltaPhiLeptonsStack");
  THStack *fDeltaEtaLeptonsStack = new THStack("DeltaEtaLeptonsStack", "DeltaEtaLeptonsStack");
  THStack *fDileptonMassStack = new THStack("DileptonMassStack", "DileptonMassStack");
  THStack *fDileptonMass_eeStack = new THStack("DileptonMass_eeStack", "DileptonMass_eeStack");
  THStack *fDileptonMass_mumuStack = new THStack("DileptonMass_mumuStack", "DileptonMass_mumuStack");
  THStack *fDileptonMass_emuStack = new THStack("DileptonMass_emuStack", "DileptonMass_emuStack");

  THStack *fMinDeltaPhiLeptonMet_afterCutsStack = new THStack("MinDeltaPhiLeptonMet_afterCutsStack", "MinDeltaPhiLeptonMet_afterCutsStack");
  THStack *fMtLepton1_afterCutsStack = new THStack("MtLepton1_afterCutsStack", "MtLepton1_afterCutsStack");
  THStack *fMtLepton2_afterCutsStack = new THStack("MtLepton2_afterCutsStack", "MtLepton2_afterCutsStack");
  THStack *fMtHiggs_afterCutsStack = new THStack("MtHiggs_afterCutsStack", "MtHiggs_afterCutsStack");
  THStack *fLeptonPtPlusMet_afterCutsStack = new THStack("LeptonPtPlusMet_afterCutsStack", "LeptonPtPlusMet_afterCutsStack");

  vector<TH1F*> fHWWSelectionDistributions;
  vector<TH1F*> fHWWToEESelectionDistributions;
  vector<TH1F*> fHWWToMuMuSelectionDistributions;
  vector<TH1F*> fHWWToEMuSelectionDistributions;

  vector<TH1F*> fLeptonEtaDistributions;
  vector<TH1F*> fLeptonPtMaxDistributions;
  vector<TH1F*> fLeptonPtMinDistributions;
  vector<TH1F*> fMetPtHistDistributions;
  vector<TH1F*> fMetPhiHistDistributions;
  vector<TH1F*> fDeltaPhiLeptonsDistributions;
  vector<TH1F*> fDeltaEtaLeptonsDistributions;
  vector<TH1F*> fDileptonMassDistributions;
  vector<TH1F*> fDileptonMass_eeDistributions;
  vector<TH1F*> fDileptonMass_mumuDistributions;
  vector<TH1F*> fDileptonMass_emuDistributions;

  vector<TH1F*> fMinDeltaPhiLeptonMet_afterCutsDistributions;
  vector<TH1F*> fMtLepton1_afterCutsDistributions;
  vector<TH1F*> fMtLepton2_afterCutsDistributions;
  vector<TH1F*> fMtHiggs_afterCutsDistributions;
  vector<TH1F*> fLeptonPtPlusMet_afterCutsDistributions;

  vector<TH1F*> fHWWSelectionDistributions_systematicErr;
  vector<TH1F*> fHWWToEESelectionDistributions_systematicErr;
  vector<TH1F*> fHWWToMuMuSelectionDistributions_systematicErr;
  vector<TH1F*> fHWWToEMuSelectionDistributions_systematicErr;

  vector<TH1F*> fLeptonEtaDistributions_systematicErr;
  vector<TH1F*> fLeptonPtMaxDistributions_systematicErr;
  vector<TH1F*> fLeptonPtMinDistributions_systematicErr;
  vector<TH1F*> fMetPtHistDistributions_systematicErr;
  vector<TH1F*> fMetPhiHistDistributions_systematicErr;
  vector<TH1F*> fDeltaPhiLeptonsDistributions_systematicErr;
  vector<TH1F*> fDeltaEtaLeptonsDistributions_systematicErr;
  vector<TH1F*> fDileptonMassDistributions_systematicErr;
  vector<TH1F*> fDileptonMass_eeDistributions_systematicErr;
  vector<TH1F*> fDileptonMass_mumuDistributions_systematicErr;
  vector<TH1F*> fDileptonMass_emuDistributions_systematicErr;

  vector<TH1F*> fMinDeltaPhiLeptonMet_afterCutsDistributions_systematicErr;
  vector<TH1F*> fMtLepton1_afterCutsDistributions_systematicErr;
  vector<TH1F*> fMtLepton2_afterCutsDistributions_systematicErr;
  vector<TH1F*> fMtHiggs_afterCutsDistributions_systematicErr;
  vector<TH1F*> fLeptonPtPlusMet_afterCutsDistributions_systematicErr;


  cout << "here2\n";
  for (int i=0; i<processLabels.size(); ++i) {
    TH1F *tempHWWSelection = new TH1F(("HWWSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
   cout << "here211\n";
   tempHWWSelection->SetFillStyle(1001);
    tempHWWSelection->SetFillColor(processColors[i]);
   cout << "here212\n";
    tempHWWSelection->SetLineWidth(1);
    fHWWSelectionDistributions.push_back(tempHWWSelection);
   cout << "here213\n";
    fHWWSelectionStack->Add(tempHWWSelection);
   cout << "here214\n";
    TH1F *tempHWWToEESelection = new TH1F(("HWWToEESelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToEESelection->SetFillStyle(1001);
    tempHWWToEESelection->SetFillColor(processColors[i]);
    tempHWWToEESelection->SetLineWidth(1);
    fHWWToEESelectionDistributions.push_back(tempHWWToEESelection);
    fHWWToEESelectionStack->Add(tempHWWToEESelection);
    tempHWWToEESelection->Sumw2();
    TH1F *tempHWWToMuMuSelection = new TH1F(("HWWToMuMuSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToMuMuSelection->SetFillStyle(1001);
    tempHWWToMuMuSelection->SetFillColor(processColors[i]);
    tempHWWToMuMuSelection->SetLineWidth(1);
    tempHWWToMuMuSelection->Sumw2();
    fHWWToMuMuSelectionDistributions.push_back(tempHWWToMuMuSelection);
    fHWWToMuMuSelectionStack->Add(tempHWWToMuMuSelection);
    TH1F *tempHWWToEMuSelection = new TH1F(("HWWToEMuSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToEMuSelection->SetFillStyle(1001);
    tempHWWToEMuSelection->SetFillColor(processColors[i]);
    tempHWWToEMuSelection->SetLineWidth(1);
    tempHWWToEMuSelection->Sumw2();
    fHWWToEMuSelectionDistributions.push_back(tempHWWToEMuSelection);
    fHWWToEMuSelectionStack->Add(tempHWWToEMuSelection);

    TH1F *tempLeptonEta = new TH1F(("LeptonEtaStack_"+processLabels[i]).c_str(),";LeptonEta;Number of Events", 100, -5, 5);    
    tempLeptonEta->SetFillStyle(1001);
    tempLeptonEta->SetFillColor(processColors[i]);
    tempLeptonEta->SetLineWidth(1);
    tempLeptonEta->Sumw2();
    fLeptonEtaDistributions.push_back(tempLeptonEta);
    fLeptonEtaStack->Add(tempLeptonEta);
    TH1F *tempLeptonPtMax = new TH1F(("LeptonPtMaxStack_"+processLabels[i]).c_str(),";Lepton P_t Max;Number of Events", 20, 0, 100);    
    tempLeptonPtMax->SetFillStyle(1001);
    tempLeptonPtMax->SetFillColor(processColors[i]);
    tempLeptonPtMax->SetLineWidth(1);
    tempLeptonPtMax->Sumw2();
    fLeptonPtMaxDistributions.push_back(tempLeptonPtMax);
    fLeptonPtMaxStack->Add(tempLeptonPtMax);
    TH1F *tempLeptonPtMin = new TH1F(("LeptonPtMinStack_"+processLabels[i]).c_str(),";Lepton P_t Min;Number of Events", 20, 0, 100);    
    tempLeptonPtMin->SetFillStyle(1001);
    tempLeptonPtMin->SetFillColor(processColors[i]);
    tempLeptonPtMin->SetLineWidth(1);
    tempLeptonPtMin->Sumw2();
    fLeptonPtMinDistributions.push_back(tempLeptonPtMin);
    fLeptonPtMinStack->Add(tempLeptonPtMin);
    TH1F *tempMetPtHist = new TH1F(("MetPtHistStack_"+processLabels[i]).c_str(),";Met;Number of Events", 50, 0, 100);    
    tempMetPtHist->SetFillStyle(1001);
    tempMetPtHist->SetFillColor(processColors[i]);
    tempMetPtHist->SetLineWidth(1);
    tempMetPtHist->Sumw2();
    fMetPtHistDistributions.push_back(tempMetPtHist);
    fMetPtHistStack->Add(tempMetPtHist);
    TH1F *tempMetPhiHist = new TH1F(("MetPhiHistStack_"+processLabels[i]).c_str(),";#phi;Number of Events", 50, -3.5, 3.5);    
    tempMetPhiHist->SetFillStyle(1001);
    tempMetPhiHist->SetFillColor(processColors[i]);
    tempMetPhiHist->SetLineWidth(1);
    tempMetPhiHist->Sumw2();
    fMetPhiHistDistributions.push_back(tempMetPhiHist);
    fMetPhiHistStack->Add(tempMetPhiHist);
    TH1F *tempDeltaPhiLeptons = new TH1F(("DeltaPhiLeptonsStack_"+processLabels[i]).c_str(),";#Delta#phi_{ll};Number of Events", 90, 0, 180);    
    tempDeltaPhiLeptons->SetFillStyle(1001);
    tempDeltaPhiLeptons->SetFillColor(processColors[i]);
    tempDeltaPhiLeptons->SetLineWidth(1);
    tempDeltaPhiLeptons->Sumw2();
    fDeltaPhiLeptonsDistributions.push_back(tempDeltaPhiLeptons);
    fDeltaPhiLeptonsStack->Add(tempDeltaPhiLeptons);
    TH1F *tempDileptonMass = new TH1F(("DileptonMassStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass->SetFillStyle(1001);
    tempDileptonMass->SetFillColor(processColors[i]);
    tempDileptonMass->SetLineWidth(1);
    tempDileptonMass->Sumw2();
    fDileptonMassDistributions.push_back(tempDileptonMass);
    fDileptonMassStack->Add(tempDileptonMass);
    TH1F *tempDileptonMass_ee = new TH1F(("DileptonMass_eeStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_ee->SetFillStyle(1001);
    tempDileptonMass_ee->SetFillColor(processColors[i]);
    tempDileptonMass_ee->SetLineWidth(1);
    tempDileptonMass_ee->Sumw2();
    fDileptonMass_eeDistributions.push_back(tempDileptonMass_ee);
    fDileptonMass_eeStack->Add(tempDileptonMass_ee);
    TH1F *tempDileptonMass_mumu = new TH1F(("DileptonMass_mumuStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_mumu->SetFillStyle(1001);
    tempDileptonMass_mumu->SetFillColor(processColors[i]);
    tempDileptonMass_mumu->SetLineWidth(1);
    tempDileptonMass_mumu->Sumw2();
    fDileptonMass_mumuDistributions.push_back(tempDileptonMass_mumu);
    fDileptonMass_mumuStack->Add(tempDileptonMass_mumu);
    TH1F *tempDileptonMass_emu = new TH1F(("DileptonMass_emuStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_emu->SetFillStyle(1001);
    tempDileptonMass_emu->SetFillColor(processColors[i]);
    tempDileptonMass_emu->SetLineWidth(1);
    tempDileptonMass_emu->Sumw2();
    fDileptonMass_emuDistributions.push_back(tempDileptonMass_emu);
    fDileptonMass_emuStack->Add(tempDileptonMass_emu);

    TH1F *tempMinDeltaPhiLeptonMet_afterCuts = new TH1F(("MinDeltaPhiLeptonMet_afterCutsStack_"+processLabels[i]).c_str(),";Min #Delta#phi_{l,Met};Number of Events", 90, 0, 180);    
    tempMinDeltaPhiLeptonMet_afterCuts->SetFillStyle(1001);
    tempMinDeltaPhiLeptonMet_afterCuts->SetFillColor(processColors[i]);
    tempMinDeltaPhiLeptonMet_afterCuts->SetLineWidth(1);
    tempMinDeltaPhiLeptonMet_afterCuts->Sumw2();
    fMinDeltaPhiLeptonMet_afterCutsDistributions.push_back(tempMinDeltaPhiLeptonMet_afterCuts);
    fMinDeltaPhiLeptonMet_afterCutsStack->Add(tempMinDeltaPhiLeptonMet_afterCuts);
    TH1F *tempMtLepton1_afterCuts = new TH1F(("MtLepton1_afterCutsStack_"+processLabels[i]).c_str(),";M_t (Lepton1,Met);Number of Events", 50,0.,200);    
    tempMtLepton1_afterCuts->SetFillStyle(1001);
    tempMtLepton1_afterCuts->SetFillColor(processColors[i]);
    tempMtLepton1_afterCuts->SetLineWidth(1);
    tempMtLepton1_afterCuts->Sumw2();
    fMtLepton1_afterCutsDistributions.push_back(tempMtLepton1_afterCuts);
    fMtLepton1_afterCutsStack->Add(tempMtLepton1_afterCuts);
    TH1F *tempMtLepton2_afterCuts = new TH1F(("MtLepton2_afterCutsStack_"+processLabels[i]).c_str(),";M_t (Lepton2,Met);Number of Events", 50, 0, 200);    
    tempMtLepton2_afterCuts->SetFillStyle(1001);
    tempMtLepton2_afterCuts->SetFillColor(processColors[i]);
    tempMtLepton2_afterCuts->SetLineWidth(1);
    tempMtLepton2_afterCuts->Sumw2();
    fMtLepton2_afterCutsDistributions.push_back(tempMtLepton2_afterCuts);
    fMtLepton2_afterCutsStack->Add(tempMtLepton2_afterCuts);
    TH1F *tempMtHiggs_afterCuts = new TH1F(("MtHiggs_afterCutsStack_"+processLabels[i]).c_str(),";M_t (l1+l2+Met);Number of Events", 50, 0, 300);    
    tempMtHiggs_afterCuts->SetFillStyle(1001);
    tempMtHiggs_afterCuts->SetFillColor(processColors[i]);
    tempMtHiggs_afterCuts->SetLineWidth(1);
    tempMtHiggs_afterCuts->Sumw2();
    fMtHiggs_afterCutsDistributions.push_back(tempMtHiggs_afterCuts);
    fMtHiggs_afterCutsStack->Add(tempMtHiggs_afterCuts);
    TH1F *tempLeptonPtPlusMet_afterCuts = new TH1F(("LeptonPtPlusMet_afterCutsStack_"+processLabels[i]).c_str(),";LeptonPtPlusMet;Number of Events", 50, 0, 300);    
    tempLeptonPtPlusMet_afterCuts->SetFillStyle(1001);
    tempLeptonPtPlusMet_afterCuts->SetFillColor(processColors[i]);
    tempLeptonPtPlusMet_afterCuts->SetLineWidth(1);
    tempLeptonPtPlusMet_afterCuts->Sumw2();
    fLeptonPtPlusMet_afterCutsDistributions.push_back(tempLeptonPtPlusMet_afterCuts);
    fLeptonPtPlusMet_afterCutsStack->Add(tempLeptonPtPlusMet_afterCuts);




    TH1F *tempHWWSelection_systematicErr = new TH1F(("HWWSelection_systematicErrStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWSelection_systematicErr->SetFillStyle(1001);
    tempHWWSelection_systematicErr->SetFillColor(processColors[i]);
    tempHWWSelection_systematicErr->SetLineWidth(1);
    fHWWSelectionDistributions_systematicErr.push_back(tempHWWSelection_systematicErr);
    TH1F *tempHWWToEESelection_systematicErr = new TH1F(("HWWToEESelection_systematicErrStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToEESelection_systematicErr->SetFillStyle(1001);
    tempHWWToEESelection_systematicErr->SetFillColor(processColors[i]);
    tempHWWToEESelection_systematicErr->SetLineWidth(1);
    fHWWToEESelectionDistributions_systematicErr.push_back(tempHWWToEESelection_systematicErr);
    TH1F *tempHWWToMuMuSelection_systematicErr = new TH1F(("HWWToMuMuSelection_systematicErrStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToMuMuSelection_systematicErr->SetFillStyle(1001);
    tempHWWToMuMuSelection_systematicErr->SetFillColor(processColors[i]);
    tempHWWToMuMuSelection_systematicErr->SetLineWidth(1);
    fHWWToMuMuSelectionDistributions_systematicErr.push_back(tempHWWToMuMuSelection_systematicErr);
    TH1F *tempHWWToEMuSelection_systematicErr = new TH1F(("HWWToEMuSelection_systematicErrStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToEMuSelection_systematicErr->SetFillStyle(1001);
    tempHWWToEMuSelection_systematicErr->SetFillColor(processColors[i]);
    tempHWWToEMuSelection_systematicErr->SetLineWidth(1);
    fHWWToEMuSelectionDistributions_systematicErr.push_back(tempHWWToEMuSelection_systematicErr);

    TH1F *tempLeptonEta_systematicErr = new TH1F(("LeptonEta_systematicErrStack_"+processLabels[i]).c_str(),";LeptonEta_systematicErr;Number of Events", 100, -5, 5);    
    tempLeptonEta_systematicErr->SetFillStyle(1001);
    tempLeptonEta_systematicErr->SetFillColor(processColors[i]);
    tempLeptonEta_systematicErr->SetLineWidth(1);
    fLeptonEtaDistributions_systematicErr.push_back(tempLeptonEta_systematicErr);
    TH1F *tempLeptonPtMax_systematicErr = new TH1F(("LeptonPtMax_systematicErrStack_"+processLabels[i]).c_str(),";Lepton P_t Max;Number of Events", 20, 0, 100);    
    tempLeptonPtMax_systematicErr->SetFillStyle(1001);
    tempLeptonPtMax_systematicErr->SetFillColor(processColors[i]);
    tempLeptonPtMax_systematicErr->SetLineWidth(1);
    fLeptonPtMaxDistributions_systematicErr.push_back(tempLeptonPtMax_systematicErr);
    TH1F *tempLeptonPtMin_systematicErr = new TH1F(("LeptonPtMin_systematicErrStack_"+processLabels[i]).c_str(),";Lepton P_t Min;Number of Events", 20, 0, 100);    
    tempLeptonPtMin_systematicErr->SetFillStyle(1001);
    tempLeptonPtMin_systematicErr->SetFillColor(processColors[i]);
    tempLeptonPtMin_systematicErr->SetLineWidth(1);
    fLeptonPtMinDistributions_systematicErr.push_back(tempLeptonPtMin_systematicErr);
    TH1F *tempMetPtHist_systematicErr = new TH1F(("MetPtHist_systematicErrStack_"+processLabels[i]).c_str(),";Met;Number of Events", 50, 0, 100);    
    tempMetPtHist_systematicErr->SetFillStyle(1001);
    tempMetPtHist_systematicErr->SetFillColor(processColors[i]);
    tempMetPtHist_systematicErr->SetLineWidth(1);
    fMetPtHistDistributions_systematicErr.push_back(tempMetPtHist_systematicErr);
    TH1F *tempMetPhiHist_systematicErr = new TH1F(("MetPhiHist_systematicErrStack_"+processLabels[i]).c_str(),";#phi;Number of Events", 50, -3.5, 3.5);    
    tempMetPhiHist_systematicErr->SetFillStyle(1001);
    tempMetPhiHist_systematicErr->SetFillColor(processColors[i]);
    tempMetPhiHist_systematicErr->SetLineWidth(1);
    fMetPhiHistDistributions_systematicErr.push_back(tempMetPhiHist_systematicErr);
    TH1F *tempDeltaPhiLeptons_systematicErr = new TH1F(("DeltaPhiLeptons_systematicErrStack_"+processLabels[i]).c_str(),";#Delta#phi_{ll};Number of Events", 90, 0, 180);    
    tempDeltaPhiLeptons_systematicErr->SetFillStyle(1001);
    tempDeltaPhiLeptons_systematicErr->SetFillColor(processColors[i]);
    tempDeltaPhiLeptons_systematicErr->SetLineWidth(1);
    fDeltaPhiLeptonsDistributions_systematicErr.push_back(tempDeltaPhiLeptons_systematicErr);
    TH1F *tempDileptonMass_systematicErr = new TH1F(("DileptonMass_systematicErrStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_systematicErr->SetFillStyle(1001);
    tempDileptonMass_systematicErr->SetFillColor(processColors[i]);
    tempDileptonMass_systematicErr->SetLineWidth(1);
    fDileptonMassDistributions_systematicErr.push_back(tempDileptonMass_systematicErr);
    TH1F *tempDileptonMass_ee_systematicErr = new TH1F(("DileptonMass_ee_systematicErrStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_ee_systematicErr->SetFillStyle(1001);
    tempDileptonMass_ee_systematicErr->SetFillColor(processColors[i]);
    tempDileptonMass_ee_systematicErr->SetLineWidth(1);
    fDileptonMass_eeDistributions_systematicErr.push_back(tempDileptonMass_ee_systematicErr);
   TH1F *tempDileptonMass_mumu_systematicErr = new TH1F(("DileptonMass_mumu_systematicErrStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_mumu_systematicErr->SetFillStyle(1001);
    tempDileptonMass_mumu_systematicErr->SetFillColor(processColors[i]);
    tempDileptonMass_mumu_systematicErr->SetLineWidth(1);
    fDileptonMass_mumuDistributions_systematicErr.push_back(tempDileptonMass_mumu_systematicErr);
    TH1F *tempDileptonMass_emu_systematicErr = new TH1F(("DileptonMass_emu_systematicErrStack_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_emu_systematicErr->SetFillStyle(1001);
    tempDileptonMass_emu_systematicErr->SetFillColor(processColors[i]);
    tempDileptonMass_emu_systematicErr->SetLineWidth(1);
    fDileptonMass_emuDistributions_systematicErr.push_back(tempDileptonMass_emu_systematicErr);

    TH1F *tempMinDeltaPhiLeptonMet_afterCuts_systematicErr = new TH1F(("MinDeltaPhiLeptonMet_afterCuts_systematicErrStack_"+processLabels[i]).c_str(),";Min #Delta#phi_{l,Met};Number of Events", 90, 0, 180);    
    tempMinDeltaPhiLeptonMet_afterCuts_systematicErr->SetFillStyle(1001);
    tempMinDeltaPhiLeptonMet_afterCuts_systematicErr->SetFillColor(processColors[i]);
    tempMinDeltaPhiLeptonMet_afterCuts_systematicErr->SetLineWidth(1);
    fMinDeltaPhiLeptonMet_afterCutsDistributions_systematicErr.push_back(tempMinDeltaPhiLeptonMet_afterCuts_systematicErr);
    TH1F *tempMtLepton1_afterCuts_systematicErr = new TH1F(("MtLepton1_afterCuts_systematicErrStack_"+processLabels[i]).c_str(),";M_t (Lepton1,Met);Number of Events", 50,0.,200);    
    tempMtLepton1_afterCuts_systematicErr->SetFillStyle(1001);
    tempMtLepton1_afterCuts_systematicErr->SetFillColor(processColors[i]);
    tempMtLepton1_afterCuts_systematicErr->SetLineWidth(1);
    fMtLepton1_afterCutsDistributions_systematicErr.push_back(tempMtLepton1_afterCuts_systematicErr);
    TH1F *tempMtLepton2_afterCuts_systematicErr = new TH1F(("MtLepton2_afterCuts_systematicErrStack_"+processLabels[i]).c_str(),";M_t (Lepton2,Met);Number of Events", 50, 0, 200);    
    tempMtLepton2_afterCuts_systematicErr->SetFillStyle(1001);
    tempMtLepton2_afterCuts_systematicErr->SetFillColor(processColors[i]);
    tempMtLepton2_afterCuts_systematicErr->SetLineWidth(1);
    fMtLepton2_afterCutsDistributions_systematicErr.push_back(tempMtLepton2_afterCuts_systematicErr);
    TH1F *tempMtHiggs_afterCuts_systematicErr = new TH1F(("MtHiggs_afterCuts_systematicErrStack_"+processLabels[i]).c_str(),";M_t (l1+l2+Met);Number of Events", 50, 0, 300);    
    tempMtHiggs_afterCuts_systematicErr->SetFillStyle(1001);
    tempMtHiggs_afterCuts_systematicErr->SetFillColor(processColors[i]);
    tempMtHiggs_afterCuts_systematicErr->SetLineWidth(1);
    fMtHiggs_afterCutsDistributions_systematicErr.push_back(tempMtHiggs_afterCuts_systematicErr);
    TH1F *tempLeptonPtPlusMet_afterCuts_systematicErr = new TH1F(("LeptonPtPlusMet_afterCuts_systematicErrStack_"+processLabels[i]).c_str(),";LeptonPtPlusMet;Number of Events", 50, 0, 300);    
    tempLeptonPtPlusMet_afterCuts_systematicErr->SetFillStyle(1001);
    tempLeptonPtPlusMet_afterCuts_systematicErr->SetFillColor(processColors[i]);
    tempLeptonPtPlusMet_afterCuts_systematicErr->SetLineWidth(1);
    fLeptonPtPlusMet_afterCutsDistributions_systematicErr.push_back(tempLeptonPtPlusMet_afterCuts_systematicErr);
 
  }
  cout << "here3\n";
 



  //--------------------------------------------------------------------------------------------------------------
  //Simulation Histograms
  //==============================================================================================================  
  THStack *fHWWSelectionStack_Simulation = new THStack("HWWSelectionStack_Simulation", "HWWSelectionStack_Simulation");
  THStack *fHWWToEESelectionStack_Simulation = new THStack("HWWToEESelectionStack_Simulation", "HWWToEESelectionStack_Simulation");
  THStack *fHWWToMuMuSelectionStack_Simulation = new THStack("HWWToMuMuSelectionStack_Simulation", "HWWToMuMuSelectionStack_Simulation");
  THStack *fHWWToEMuSelectionStack_Simulation = new THStack("HWWToEMuSelectionStack_Simulation", "HWWToEMuSelectionStack_Simulation");

  THStack *fLeptonEtaStack_Simulation = new THStack("LeptonEtaStack_Simulation", "LeptonEtaStack_Simulation");
  THStack *fLeptonPtMaxStack_Simulation = new THStack("LeptonPtMaxStack_Simulation", "LeptonPtMaxStack_Simulation");
  THStack *fLeptonPtMinStack_Simulation = new THStack("LeptonPtMinStack_Simulation", "LeptonPtMinStack_Simulation");
  THStack *fMetPtHistStack_Simulation = new THStack("MetPtHistStack_Simulation", "MetPtHistStack_Simulation");
  THStack *fMetPhiHistStack_Simulation = new THStack("MetPhiHistStack_Simulation", "MetPhiHistStack_Simulation");
  THStack *fDeltaPhiLeptonsStack_Simulation = new THStack("DeltaPhiLeptonsStack_Simulation", "DeltaPhiLeptonsStack_Simulation");
  THStack *fDeltaEtaLeptonsStack_Simulation = new THStack("DeltaEtaLeptonsStack_Simulation", "DeltaEtaLeptonsStack_Simulation");
  THStack *fDileptonMassStack_Simulation = new THStack("DileptonMassStack_Simulation", "DileptonMassStack_Simulation");
  THStack *fDileptonMass_eeStack_Simulation = new THStack("DileptonMass_eeStack_Simulation", "DileptonMass_eeStack_Simulation");
  THStack *fDileptonMass_mumuStack_Simulation = new THStack("DileptonMass_mumuStack_Simulation", "DileptonMass_mumuStack_Simulation");
  THStack *fDileptonMass_emuStack_Simulation = new THStack("DileptonMass_emuStack_Simulation", "DileptonMass_emuStack_Simulation");

  THStack *fMinDeltaPhiLeptonMet_afterCutsStack_Simulation = new THStack("MinDeltaPhiLeptonMet_afterCutsStack_Simulation", "MinDeltaPhiLeptonMet_afterCutsStack_Simulation");
  THStack *fMtLepton1_afterCutsStack_Simulation = new THStack("MtLepton1_afterCutsStack_Simulation", "MtLepton1_afterCutsStack_Simulation");
  THStack *fMtLepton2_afterCutsStack_Simulation = new THStack("MtLepton2_afterCutsStack_Simulation", "MtLepton2_afterCutsStack_Simulation");
  THStack *fMtHiggs_afterCutsStack_Simulation = new THStack("MtHiggs_afterCutsStack_Simulation", "MtHiggs_afterCutsStack_Simulation");
  THStack *fLeptonPtPlusMet_afterCutsStack_Simulation = new THStack("LeptonPtPlusMet_afterCutsStack_Simulation", "LeptonPtPlusMet_afterCutsStack_Simulation");

  vector<TH1F*> fHWWSelectionDistributions_Simulation;
  vector<TH1F*> fHWWToEESelectionDistributions_Simulation;
  vector<TH1F*> fHWWToMuMuSelectionDistributions_Simulation;
  vector<TH1F*> fHWWToEMuSelectionDistributions_Simulation;

  vector<TH1F*> fLeptonEtaDistributions_Simulation;
  vector<TH1F*> fLeptonPtMaxDistributions_Simulation;
  vector<TH1F*> fLeptonPtMinDistributions_Simulation;
  vector<TH1F*> fMetPtHistDistributions_Simulation;
  vector<TH1F*> fMetPhiHistDistributions_Simulation;
  vector<TH1F*> fDeltaPhiLeptonsDistributions_Simulation;
  vector<TH1F*> fDeltaEtaLeptonsDistributions_Simulation;
  vector<TH1F*> fDileptonMassDistributions_Simulation;
  vector<TH1F*> fDileptonMass_eeDistributions_Simulation;
  vector<TH1F*> fDileptonMass_mumuDistributions_Simulation;
  vector<TH1F*> fDileptonMass_emuDistributions_Simulation;

  vector<TH1F*> fMinDeltaPhiLeptonMet_afterCutsDistributions_Simulation;
  vector<TH1F*> fMtLepton1_afterCutsDistributions_Simulation;
  vector<TH1F*> fMtLepton2_afterCutsDistributions_Simulation;
  vector<TH1F*> fMtHiggs_afterCutsDistributions_Simulation;
  vector<TH1F*> fLeptonPtPlusMet_afterCutsDistributions_Simulation;

  for (int i=0; i<processLabels.size(); ++i) {
    TH1F *tempHWWSelection_Simulation = new TH1F(("HWWSelectionStack_Simulation_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWSelection_Simulation->SetFillStyle(1001);
    tempHWWSelection_Simulation->SetFillColor(processColors[i]);
    tempHWWSelection_Simulation->SetLineWidth(1);
    fHWWSelectionDistributions_Simulation.push_back(tempHWWSelection_Simulation);
    fHWWSelectionStack_Simulation->Add(tempHWWSelection_Simulation);
    TH1F *tempHWWToEESelection_Simulation = new TH1F(("HWWToEESelectionStack_Simulation_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToEESelection_Simulation->SetFillStyle(1001);
    tempHWWToEESelection_Simulation->SetFillColor(processColors[i]);
    tempHWWToEESelection_Simulation->SetLineWidth(1);
    fHWWToEESelectionDistributions_Simulation.push_back(tempHWWToEESelection_Simulation);
    fHWWToEESelectionStack_Simulation->Add(tempHWWToEESelection_Simulation);
    TH1F *tempHWWToMuMuSelection_Simulation = new TH1F(("HWWToMuMuSelectionStack_Simulation_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToMuMuSelection_Simulation->SetFillStyle(1001);
    tempHWWToMuMuSelection_Simulation->SetFillColor(processColors[i]);
    tempHWWToMuMuSelection_Simulation->SetLineWidth(1);
    fHWWToMuMuSelectionDistributions_Simulation.push_back(tempHWWToMuMuSelection_Simulation);
    fHWWToMuMuSelectionStack_Simulation->Add(tempHWWToMuMuSelection_Simulation);
    TH1F *tempHWWToEMuSelection_Simulation = new TH1F(("HWWToEMuSelectionStack_Simulation_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 11, -1.5, 9.5);    
    tempHWWToEMuSelection_Simulation->SetFillStyle(1001);
    tempHWWToEMuSelection_Simulation->SetFillColor(processColors[i]);
    tempHWWToEMuSelection_Simulation->SetLineWidth(1);
    fHWWToEMuSelectionDistributions_Simulation.push_back(tempHWWToEMuSelection_Simulation);
    fHWWToEMuSelectionStack_Simulation->Add(tempHWWToEMuSelection_Simulation);
  
    TH1F *tempLeptonEta_Simulation = new TH1F(("LeptonEtaStack_Simulation_"+processLabels[i]).c_str(),";LeptonEta;Number of Events", 100, -5, 5);    
    tempLeptonEta_Simulation->SetFillStyle(1001);
    tempLeptonEta_Simulation->SetFillColor(processColors[i]);
    tempLeptonEta_Simulation->SetLineWidth(1);
    fLeptonEtaDistributions_Simulation.push_back(tempLeptonEta_Simulation);
    fLeptonEtaStack_Simulation->Add(tempLeptonEta_Simulation);
    TH1F *tempLeptonPtMax_Simulation = new TH1F(("LeptonPtMaxStack_Simulation_"+processLabels[i]).c_str(),";Lepton P_t Max;Number of Events", 20, 0, 100);    
    tempLeptonPtMax_Simulation->SetFillStyle(1001);
    tempLeptonPtMax_Simulation->SetFillColor(processColors[i]);
    tempLeptonPtMax_Simulation->SetLineWidth(1);
    fLeptonPtMaxDistributions_Simulation.push_back(tempLeptonPtMax_Simulation);
    fLeptonPtMaxStack_Simulation->Add(tempLeptonPtMax_Simulation);
    TH1F *tempLeptonPtMin_Simulation = new TH1F(("LeptonPtMinStack_Simulation_"+processLabels[i]).c_str(),";Lepton P_t Min;Number of Events", 20, 0, 100);    
    tempLeptonPtMin_Simulation->SetFillStyle(1001);
    tempLeptonPtMin_Simulation->SetFillColor(processColors[i]);
    tempLeptonPtMin_Simulation->SetLineWidth(1);
    fLeptonPtMinDistributions_Simulation.push_back(tempLeptonPtMin_Simulation);
    fLeptonPtMinStack_Simulation->Add(tempLeptonPtMin_Simulation);
    TH1F *tempMetPtHist_Simulation = new TH1F(("MetPtHistStack_Simulation_"+processLabels[i]).c_str(),";Met;Number of Events", 50, 0, 100);    
    tempMetPtHist_Simulation->SetFillStyle(1001);
    tempMetPtHist_Simulation->SetFillColor(processColors[i]);
    tempMetPtHist_Simulation->SetLineWidth(1);
    fMetPtHistDistributions_Simulation.push_back(tempMetPtHist_Simulation);
    fMetPtHistStack_Simulation->Add(tempMetPtHist_Simulation);
    TH1F *tempMetPhiHist_Simulation = new TH1F(("MetPhiHistStack_Simulation_"+processLabels[i]).c_str(),";#phi;Number of Events", 50, -3.5, 3.5);    
    tempMetPhiHist_Simulation->SetFillStyle(1001);
    tempMetPhiHist_Simulation->SetFillColor(processColors[i]);
    tempMetPhiHist_Simulation->SetLineWidth(1);
    fMetPhiHistDistributions_Simulation.push_back(tempMetPhiHist_Simulation);
    fMetPhiHistStack_Simulation->Add(tempMetPhiHist_Simulation);
    TH1F *tempDeltaPhiLeptons_Simulation = new TH1F(("DeltaPhiLeptonsStack_Simulation_"+processLabels[i]).c_str(), ";#Delta#phi_{ll};Number of Events", 90, 0, 180); 
    tempDeltaPhiLeptons_Simulation->SetFillStyle(1001);
    tempDeltaPhiLeptons_Simulation->SetFillColor(processColors[i]);
    tempDeltaPhiLeptons_Simulation->SetLineWidth(1);
    fDeltaPhiLeptonsDistributions_Simulation.push_back(tempDeltaPhiLeptons_Simulation);
    fDeltaPhiLeptonsStack_Simulation->Add(tempDeltaPhiLeptons_Simulation);
    TH1F *tempDileptonMass_Simulation = new TH1F(("DileptonMassStack_Simulation_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_Simulation->SetFillStyle(1001);
    tempDileptonMass_Simulation->SetFillColor(processColors[i]);
    tempDileptonMass_Simulation->SetLineWidth(1);
    fDileptonMassDistributions_Simulation.push_back(tempDileptonMass_Simulation);
    fDileptonMassStack_Simulation->Add(tempDileptonMass_Simulation);
    TH1F *tempDileptonMass_ee_Simulation = new TH1F(("DileptonMass_eeStack_Simulation_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_ee_Simulation->SetFillStyle(1001);
    tempDileptonMass_ee_Simulation->SetFillColor(processColors[i]);
    tempDileptonMass_ee_Simulation->SetLineWidth(1);
    fDileptonMass_eeDistributions_Simulation.push_back(tempDileptonMass_ee_Simulation);
    fDileptonMass_eeStack_Simulation->Add(tempDileptonMass_ee_Simulation);
    TH1F *tempDileptonMass_mumu_Simulation = new TH1F(("DileptonMass_mumuStack_Simulation_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_mumu_Simulation->SetFillStyle(1001);
    tempDileptonMass_mumu_Simulation->SetFillColor(processColors[i]);
    tempDileptonMass_mumu_Simulation->SetLineWidth(1);
    fDileptonMass_mumuDistributions_Simulation.push_back(tempDileptonMass_mumu_Simulation);
    fDileptonMass_mumuStack_Simulation->Add(tempDileptonMass_mumu_Simulation);
    TH1F *tempDileptonMass_emu_Simulation = new TH1F(("DileptonMass_emuStack_Simulation_"+processLabels[i]).c_str(),";Mass [GeV/c^{2}];Number of Events", 50, 0, 200);    
    tempDileptonMass_emu_Simulation->SetFillStyle(1001);
    tempDileptonMass_emu_Simulation->SetFillColor(processColors[i]);
    tempDileptonMass_emu_Simulation->SetLineWidth(1);
    fDileptonMass_emuDistributions_Simulation.push_back(tempDileptonMass_emu_Simulation);
    fDileptonMass_emuStack_Simulation->Add(tempDileptonMass_emu_Simulation);

  
    TH1F *tempMinDeltaPhiLeptonMet_afterCuts_Simulation = new TH1F(("MinDeltaPhiLeptonMet_afterCutsStack_Simulation_"+processLabels[i]).c_str(),";Min #Delta#phi_{l,Met};Number of Events", 90, 0, 180);    
    tempMinDeltaPhiLeptonMet_afterCuts_Simulation->SetFillStyle(1001);
    tempMinDeltaPhiLeptonMet_afterCuts_Simulation->SetFillColor(processColors[i]);
    tempMinDeltaPhiLeptonMet_afterCuts_Simulation->SetLineWidth(1);
    fMinDeltaPhiLeptonMet_afterCutsDistributions_Simulation.push_back(tempMinDeltaPhiLeptonMet_afterCuts_Simulation);
    fMinDeltaPhiLeptonMet_afterCutsStack_Simulation->Add(tempMinDeltaPhiLeptonMet_afterCuts_Simulation);
    TH1F *tempMtLepton1_afterCuts_Simulation = new TH1F(("MtLepton1_afterCutsStack_Simulation_"+processLabels[i]).c_str(),";M_t (Lepton1,Met);Number of Events", 50,0.,200);    
    tempMtLepton1_afterCuts_Simulation->SetFillStyle(1001);
    tempMtLepton1_afterCuts_Simulation->SetFillColor(processColors[i]);
    tempMtLepton1_afterCuts_Simulation->SetLineWidth(1);
    fMtLepton1_afterCutsDistributions_Simulation.push_back(tempMtLepton1_afterCuts_Simulation);
    fMtLepton1_afterCutsStack_Simulation->Add(tempMtLepton1_afterCuts_Simulation);
    TH1F *tempMtLepton2_afterCuts_Simulation = new TH1F(("MtLepton2_afterCutsStack_Simulation_"+processLabels[i]).c_str(),";M_t (Lepton2,Met);Number of Events", 50, 0, 200);    
    tempMtLepton2_afterCuts_Simulation->SetFillStyle(1001);
    tempMtLepton2_afterCuts_Simulation->SetFillColor(processColors[i]);
    tempMtLepton2_afterCuts_Simulation->SetLineWidth(1);
    fMtLepton2_afterCutsDistributions_Simulation.push_back(tempMtLepton2_afterCuts_Simulation);
    fMtLepton2_afterCutsStack_Simulation->Add(tempMtLepton2_afterCuts_Simulation);
    TH1F *tempMtHiggs_afterCuts_Simulation = new TH1F(("MtHiggs_afterCutsStack_Simulation_"+processLabels[i]).c_str(),";M_t (l1+l2+Met);Number of Events", 50, 0, 300);    
    tempMtHiggs_afterCuts_Simulation->SetFillStyle(1001);
    tempMtHiggs_afterCuts_Simulation->SetFillColor(processColors[i]);
    tempMtHiggs_afterCuts_Simulation->SetLineWidth(1);
    fMtHiggs_afterCutsDistributions_Simulation.push_back(tempMtHiggs_afterCuts_Simulation);
    fMtHiggs_afterCutsStack_Simulation->Add(tempMtHiggs_afterCuts_Simulation);
    TH1F *tempLeptonPtPlusMet_afterCuts_Simulation = new TH1F(("LeptonPtPlusMet_afterCutsStack_Simulation_"+processLabels[i]).c_str(),";LeptonPtPlusMet;Number of Events", 50, 0, 300);    
    tempLeptonPtPlusMet_afterCuts_Simulation->SetFillStyle(1001);
    tempLeptonPtPlusMet_afterCuts_Simulation->SetFillColor(processColors[i]);
    tempLeptonPtPlusMet_afterCuts_Simulation->SetLineWidth(1);
    fLeptonPtPlusMet_afterCutsDistributions_Simulation.push_back(tempLeptonPtPlusMet_afterCuts_Simulation);
    fLeptonPtPlusMet_afterCutsStack_Simulation->Add(tempLeptonPtPlusMet_afterCuts_Simulation);
  
  }



  TH1F *denominator_Pt = new TH1F("denominator_Pt" , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
  TH1F *numerator_Pt = new TH1F("numerator_Pt" , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
  TGraphAsymmErrors *efficiency_pt = 0;

  ofstream eventListFile("mcEventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TFile *inputFile=0;
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;

  cout << "here\n";

  for (int p=0; p<processInputFiles.size() ; ++p) {
    for (int f=0; f<processInputFiles[p].size() ; ++f) {

      cout << "here " << p << " " << f << endl;
//       inputFile = new TFile(processInputFiles[p][f].c_str());
//       assert(inputFile);


      //********************************************************
      // Get Tree
      //********************************************************
      
//       eventTree = (TTree*)inputFile->Get("Events"); assert(eventTree);
      eventTree = getTreeFromFile(processInputFiles[p][f].c_str(),"Events");

      //*****************************************************************************************
      //Loop over Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
  
      //********************************************************
      // Calculate event weight
      //********************************************************
      assert(eventTree->GetEntries() > 0);
      infoBr->GetEntry(0);
//       cout << processNames[p][f] << endl;
//       Double_t eventWeight = getNormalizationWeight(("/home/sixie/hist/WWAnalysis/WWAnalysis_"+processNames[p][f]+"_noskim.root").c_str(), processNames[p][f]) * LUMINOSITY;
      Double_t eventWeight = info->eventweight*LUMINOSITY;
      //eventWeight = 1;
      cout << eventWeight << endl;


      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
		
//       if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;
        if (ientry % 100000 == 0) {
          cout << "Event " << ientry << endl;
        }
  
        //veto VGamma events
        if (info->VGammaEvent) continue;

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
	
        //********************************************************
        // TcMet
        //********************************************************

        Int_t NSoftMuons = 0;
        Int_t NLeptons = 0;
        vector<Int_t> leptonType;
        vector<Int_t> leptonIndex;
        vector<Double_t> leptonPt;
        vector<Double_t> leptonEta;
        vector<Double_t> leptonPhi;
        vector<Int_t> leptonCharge;

        Int_t NJets = 0;
        const mithep::TJet *leadingJet = 0;
    
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
          
          if ( (0==0)
               &&
               passMuonCuts(mu)
               &&
               fabs(mu->eta) < 2.4
               && 
               mu->pt > 20.0
            ) {
            leptonPt.push_back(mu->pt);
            leptonEta.push_back(mu->eta);
            leptonPhi.push_back(mu->phi);
            leptonType.push_back(13);
            leptonIndex.push_back(i);  
            leptonCharge.push_back(mu->q);
          }
        }

        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
          Bool_t isMuonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isMuonOverlap = kTRUE; 
              break;
            }        
          }

          if ( (0==0)
               && 
               passElectronCuts(ele)
               &&
               fabs(ele->eta) < 2.5
               && 
               ele->pt > 20.0
               &&
               !isMuonOverlap
            ) {
            leptonPt.push_back(ele->pt);
            leptonEta.push_back(ele->eta);
            leptonPhi.push_back(ele->phi);
            leptonType.push_back(11);
            leptonIndex.push_back(i);
            leptonCharge.push_back(ele->q);
          }
        }

 //        if (leptonPt.size() >= 2) {
//           cout << "Dilepton Event:" << info->runNum << " " << info->evtNum << endl;
//           for (int i=0; i<leptonPt.size() ; ++i) {
//             cout << "Lepton " << i << " " << leptonPt[i] << " " << leptonEta[i] << " " << leptonPhi[i] << " " << leptonType[i] << " " << leptonCharge[i] << endl;
//           }
//         }

        //sort leptons
        Int_t tempType;
        Int_t tempIndex;
        Double_t tempPt;
        Double_t tempEta;
        Double_t tempPhi;
        Int_t tempCharge;
        for (int l=0; l<leptonIndex.size(); l++) {
          for (int k=0; k < leptonIndex.size() - 1; k++) {
            if (leptonPt[k+1] > leptonPt[k]) {
              tempType = leptonType[k];
              tempIndex = leptonIndex[k];
              tempPt = leptonPt[k];
              tempEta = leptonEta[k];
              tempPhi = leptonPhi[k];
              tempCharge = leptonCharge[k];
          
              leptonType[k] = leptonType[k+1];
              leptonIndex[k] = leptonIndex[k+1];
              leptonPt[k] = leptonPt[k+1];
              leptonEta[k] = leptonEta[k+1];
              leptonPhi[k] = leptonPhi[k+1];
              leptonCharge[k] = leptonCharge[k+1];

              leptonType[k+1] = tempType;
              leptonIndex[k+1] = tempIndex;
              leptonPt[k+1] = tempPt;
              leptonEta[k+1] = tempEta;
              leptonPhi[k+1] = tempPhi;
              leptonCharge[k+1] = tempCharge;
          
            }
          }
        }    


        //******************************************************************************
        //Debug Test fake rates
        //******************************************************************************
        
        for(int i = 0; i < leptonPt.size(); ++i) {
          if (leptonType[i] == 13) {

            for(Int_t j=0; j<electronArr->GetEntries(); j++) {
              const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[j]);
              
              if (!(ele->pt > 20.0 && fabs(ele->eta) < 2.5)) continue;
              if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[i],leptonEta[i]) < 0.1 )
                continue;
              
              //select denominators that fail final selection
              if (!passElectronDenominatorCuts(ele)) continue;

              denominator_Pt->Fill(ele->pt);
              if (passElectronCuts(ele)) {
                numerator_Pt->Fill(ele->pt);
              }                
            }
          }
        }

      


        


        //******************************************************************************
        //Fake Rate Method Selection : select events with 1lepton only
        //******************************************************************************
        if (leptonPt.size() >= 1 && leptonPt[0] > 20.0 
            && (leptonPt.size() < 2 || leptonPt[1] < 10.0)
          ) {
          
          //loop over electron denominators
          for(Int_t i=0; i<electronArr->GetEntries(); i++) {
            const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
            
            if (!(ele->pt > 20.0 && fabs(ele->eta) < 2.5)) continue;
            if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[0],leptonEta[0]) < 0.1 )
              continue;
            
            //select denominators that fail final selection
            if (!passElectronDenominatorCuts(ele)) continue;
            if (passElectronCuts(ele)) continue;
            
            //preselection
            if ((ChargeSelection == 0 && leptonCharge[0] == ele->q) || (ChargeSelection == 1 && leptonCharge[0] != ele->q)) continue;
            
            
            Int_t finalState = -1;
            if (leptonType[0] == 11) {
              finalState = 0;
            } else if (leptonType[0] == 13) {
              finalState = 2;
            } 
           
            //******************************************************************************
            //Calculate Fake Rate
            //******************************************************************************            
            double fakeRate = 0.0;
            double fakeRateErrorLow = 0.0;
            double fakeRateErrorHigh = 0.0;

            if (fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi) < 1.0) {
              fakeRate = fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi) / (1-fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi));
              fakeRateErrorLow = fMCFakeRate->ElectronFakeRateErrorLow(ele->pt, ele->eta, ele->phi) / pow((1- fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi)),2);
              fakeRateErrorHigh = fMCFakeRate->ElectronFakeRateErrorHigh(ele->pt, ele->eta, ele->phi) / pow((1- fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi)),2);
            } else {
              fakeRate = fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi);
              fakeRateErrorLow = fMCFakeRate->ElectronFakeRateErrorLow(ele->pt, ele->eta, ele->phi);
              fakeRateErrorHigh = fMCFakeRate->ElectronFakeRateErrorHigh(ele->pt, ele->eta, ele->phi);
            }

            Double_t fakeEventWeight = eventWeight*fakeRate; 
            //eventWeight = 1.0; 
            Double_t fakeEventWeightError = eventWeight * (  fakeRateErrorLow +   fakeRateErrorHigh) / 2;
          
    //         cout << "Fake EVENT:\n";
//             cout << "Lepton " <<  leptonPt[0] << " " << leptonEta[0] << " " << leptonPhi[0] << endl;
//             cout << "Fake Electron " <<  ele->pt << " " << ele->eta << " " << ele->phi << endl;
//             cout << "fake rate" << fMCFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi) << " " << fMCFakeRate->ElectronFakeRateErrorLow(ele->pt, ele->eta, ele->phi) << " " << fMCFakeRate->ElectronFakeRateErrorHigh(ele->pt, ele->eta, ele->phi) << " , " << fakeRate << " " << fakeRateErrorLow << " " << fakeRateErrorHigh << endl;

            //******************************************************************************
            //Count Jets
            //******************************************************************************
            double maxBtag = -99999;
            NJets = 0;
            for(Int_t i=0; i<jetArr->GetEntries(); i++) {
              const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
        
              Bool_t leptonOverlap = kFALSE;
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[0],leptonEta[0]) < 0.3 ||
                  mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi,ele->eta) < 0.3 ) {
                leptonOverlap = kTRUE;
              }
        
              if (!leptonOverlap) {
                if (jet->et > 25 && fabs(jet->eta) < 5.0 ) {
                  if (!leadingJet || jet->et > leadingJet->et) {
                    leadingJet = jet;
                  }
                  NJets++;
                } else {
                  if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
                }
              }
            }
      
      
            //******************************************************************************
            //soft muons
            //******************************************************************************
            NSoftMuons = 0;
            for(Int_t j=0; j<muonArr->GetEntries(); j++) {
              const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[j]);
              Bool_t isCleanMuon = kFALSE;        
              if ( leptonType[0] == 13 
                   && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[0],leptonEta[0]) < 0.1
                ) {
                isCleanMuon = kTRUE;         
              }
              Bool_t overlapWithElectronDenominator = kFALSE;
              if ( leptonType[0] == 13 
                   && mithep::MathUtils::DeltaR(mu->phi, mu->eta, ele->phi, ele->eta) < 0.3
                ) {
                overlapWithElectronDenominator = kTRUE;         
              }
              if ( fabs(mu->eta) < 2.4                   
                   && mu->pt > 3.0
                   && (mu->qualityBits & kTMLastStationAngTight)
                   && mu->nTkHits > 10
                   && fabs(mu->d0) < 0.2
                   && (mu->typeBits & kTracker)
                   && !isCleanMuon
                   && !overlapWithElectronDenominator
                ) {
                NSoftMuons++;
              }
            }
      
            //***********************************************************************************************
            //|Z_vert-Z_l| maximum
            //***********************************************************************************************
            double zDiffMax = 0.0;
            
            double dzLepton;
            if (leptonType[0] == 11) {
              dzLepton = ((mithep::TElectron*)((*electronArr)[leptonIndex[0]]))->dz;
            } else {
              dzLepton = ((mithep::TMuon*)((*muonArr)[leptonIndex[0]]))->dz;
            }
            
            double dz = ele->dz;
//       if (dz > zDiffMax) zDiffMax = dz;
            zDiffMax = fabs(dz - dzLepton);
            


            //******************************************************************************
            //construct event variables
            //******************************************************************************
            mithep::FourVectorM lepton1;
            mithep::FourVectorM lepton2;
            if (leptonPt[0] > ele->pt) {
              if (leptonType[0] == 11) {
                lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
              } else {
                lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
              }
              lepton2.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
            } else {
              if (leptonType[0] == 11) {
                lepton2.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
              } else {
                lepton2.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
              }
              lepton1.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
            }

            mithep::FourVectorM dilepton = lepton1+lepton2;

            double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(lepton1.Phi(), 
                                                                 lepton2.Phi())* 180.0 / TMath::Pi();    
            double deltaPhiDileptonMet = mithep::MathUtils::DeltaPhi(met.Phi(), 
                                                                     dilepton.Phi())*180.0 / TMath::Pi();    
            double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt() * met.Phi()*
                                         (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

            //angle between MET and closest lepton
            double deltaPhiMetLepton[2] = {mithep::MathUtils::DeltaPhi(met.Phi(), lepton1.Phi()),
                                           mithep::MathUtils::DeltaPhi(met.Phi(), lepton2.Phi())};
  
            double mTW[2] = {TMath::Sqrt(2.0*lepton1.Pt()*met.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[0]))),
                             TMath::Sqrt(2.0*lepton2.Pt()*met.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[1])))};

            double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
              deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

            double METdeltaPhilEt = met.Pt();
            if(minDeltaPhiMetLepton < TMath::Pi()/2.)
              METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton);

    
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 10;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false};
  
            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 20.0) passCut[0] = true;
    
            if(zDiffMax < 1.0)                    passCut[1] = true;
            
            if(met.Pt()    > 20.0)               passCut[2] = true;
            
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
              if(METdeltaPhilEt > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(METdeltaPhilEt > 20) passCut[5] = true;
            }

            if(NJets     < 1)              passCut[6] = true;
        
            if (NSoftMuons == 0 )      passCut[7] = true;

            passCut[8] = true;

            if(maxBtag < 2.1)          passCut[9] = true;

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
            //Cut Selection Histograms
            fHWWSelectionDistributions[p]->Fill(-1, fakeEventWeight);
            fHWWSelectionDistributions_systematicErr[p]->Fill(-1, fakeEventWeightError);
            if (finalState == 1 ) {
              fHWWToMuMuSelectionDistributions[p]->Fill(-1, fakeEventWeight);
              fHWWToMuMuSelectionDistributions_systematicErr[p]->Fill(-1, fakeEventWeightError);
            } else if(finalState == 0 ) {
              fHWWToEESelectionDistributions[p]->Fill(-1, fakeEventWeight);
              fHWWToEESelectionDistributions_systematicErr[p]->Fill(-1, fakeEventWeightError);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelectionDistributions[p]->Fill(-1, fakeEventWeight);
              fHWWToEMuSelectionDistributions_systematicErr[p]->Fill(-1, fakeEventWeightError);
            }
    
            for (int k=0;k<nCuts;k++) {
              bool pass = true;
              bool passPreviousCut = true;
              for (int q=0;q<=k;q++) {
                pass = (pass && passCut[q]);
                if (q<k)
                  passPreviousCut = (passPreviousCut&& passCut[q]);
              }
      
              if (pass) {
                fHWWSelectionDistributions[p]->Fill(k, fakeEventWeight);
                fHWWSelectionDistributions_systematicErr[p]->Fill(k, fakeEventWeightError);
               if (finalState == 1 ) {
                  fHWWToMuMuSelectionDistributions[p]->Fill(k, fakeEventWeight);
                  fHWWToMuMuSelectionDistributions_systematicErr[p]->Fill(k, fakeEventWeightError);
                } else if(finalState == 0 ) {
                  fHWWToEESelectionDistributions[p]->Fill(k, fakeEventWeight);
                  fHWWToEESelectionDistributions_systematicErr[p]->Fill(k, fakeEventWeightError);
                } else if(finalState == 2 || finalState == 3 ) {
                  fHWWToEMuSelectionDistributions[p]->Fill(k, fakeEventWeight);
                 fHWWToEMuSelectionDistributions_systematicErr[p]->Fill(k, fakeEventWeightError);
                }
              }
            }

            //*****************************************************************************************
            //Make Preselection Histograms  
            //*****************************************************************************************
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] ) {
              fLeptonEtaDistributions[p]->Fill(lepton1.Eta(), fakeEventWeight); 
              fLeptonEtaDistributions[p]->Fill(lepton2.Eta(), fakeEventWeight);
              fLeptonPtMaxDistributions[p]->Fill(lepton1.Pt(), fakeEventWeight);
              fLeptonPtMinDistributions[p]->Fill(lepton2.Pt(), fakeEventWeight);
              fMetPtHistDistributions[p]->Fill(met.Pt(), fakeEventWeight);                             
              fMetPhiHistDistributions[p]->Fill(met.Phi(), fakeEventWeight);                            
              fDeltaPhiLeptonsDistributions[p]->Fill(deltaPhiLeptons, fakeEventWeight);
              fLeptonEtaDistributions_systematicErr[p]->Fill(lepton1.Eta(), fakeEventWeightError); 
              fLeptonEtaDistributions_systematicErr[p]->Fill(lepton2.Eta(), fakeEventWeightError);
              fLeptonPtMaxDistributions_systematicErr[p]->Fill(lepton1.Pt(), fakeEventWeightError);
              fLeptonPtMinDistributions_systematicErr[p]->Fill(lepton2.Pt(), fakeEventWeightError);
              fMetPtHistDistributions_systematicErr[p]->Fill(met.Pt(), fakeEventWeightError);                             
              fMetPhiHistDistributions_systematicErr[p]->Fill(met.Phi(), fakeEventWeightError);                            
              fDeltaPhiLeptonsDistributions_systematicErr[p]->Fill(deltaPhiLeptons, fakeEventWeightError);
              
              fDileptonMassDistributions[p]->Fill(dilepton.M(), fakeEventWeight);
              fDileptonMassDistributions_systematicErr[p]->Fill(dilepton.M(), fakeEventWeightError);
              if (finalState == 0) {
                fDileptonMass_eeDistributions[p]->Fill(dilepton.M(), fakeEventWeight);
                fDileptonMass_eeDistributions_systematicErr[p]->Fill(dilepton.M(), fakeEventWeightError);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              } else if (finalState == 1) {
                fDileptonMass_mumuDistributions[p]->Fill(dilepton.M(), fakeEventWeight);
                fDileptonMass_mumuDistributions_systematicErr[p]->Fill(dilepton.M(), fakeEventWeightError);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              } else {
                fDileptonMass_emuDistributions[p]->Fill(dilepton.M(), fakeEventWeight);
                fDileptonMass_emuDistributions_systematicErr[p]->Fill(dilepton.M(), fakeEventWeightError);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              }
            }

            //*********************************************************************************************
            //Plots after all Cuts
            //*********************************************************************************************
            if (passAllCuts) {
              fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->Fill(minDeltaPhiMetLepton, fakeEventWeight);
              fMtLepton1_afterCutsDistributions[p]->Fill(mTW[0], fakeEventWeight);
              fMtLepton2_afterCutsDistributions[p]->Fill(mTW[1], fakeEventWeight);
              fMtHiggs_afterCutsDistributions[p]->Fill(mtHiggs, fakeEventWeight);
              fLeptonPtPlusMet_afterCutsDistributions[p]->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt(), fakeEventWeight);
              fMinDeltaPhiLeptonMet_afterCutsDistributions_systematicErr[p]->Fill(minDeltaPhiMetLepton, fakeEventWeightError);
              fMtLepton1_afterCutsDistributions_systematicErr[p]->Fill(mTW[0], fakeEventWeightError);
              fMtLepton2_afterCutsDistributions_systematicErr[p]->Fill(mTW[1], fakeEventWeightError);
              fMtHiggs_afterCutsDistributions_systematicErr[p]->Fill(mtHiggs, fakeEventWeightError);
              fLeptonPtPlusMet_afterCutsDistributions_systematicErr[p]->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt(), fakeEventWeightError);
    
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);

            }
          } //loop over electron denominators
        } //end if fake rate selection

        
        //******************************************************************************
        //Standard WW Selection for Simulation based prediction
        //******************************************************************************

        for(int i = 0; i < leptonPt.size(); ++i) {
          for(int j = i+1; j < leptonPt.size(); ++j) {
              
            //impose pt cuts
            if (!(leptonPt[i] > 20.0 && fabs(leptonEta[i]) < 2.5
                  && leptonPt[j] > 20.0 && fabs(leptonEta[j]) < 2.5)) continue;
   
            //require opposite sign
            if (fabs(leptonCharge[i] + leptonCharge[j]) != ChargeSelection) continue;
                         
            Int_t finalState = -1;
            if (leptonType[i] == 11 && leptonType[j] == 11) {
              finalState = 0;
            } else if (leptonType[i] == 13 && leptonType[j] == 13) {
              finalState = 1;
              continue;
            } else if (leptonType[i] == 11 && leptonType[j] == 13) {
              finalState = 2;
            } else if (leptonType[i] == 13 && leptonType[j] == 11) {
              finalState = 3;
            }
     
            //******************************************************************************
            //Count Jets
            //******************************************************************************
            double maxBtag = -99999;
            for(Int_t i=0; i<jetArr->GetEntries(); i++) {
              const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

              Bool_t leptonOverlap = kFALSE;
              for (int k=0; k<leptonPt.size(); ++k) {
                if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
                  leptonOverlap = kTRUE;
                }
              }

              if (!leptonOverlap) {
                if (jet->et > 25 && fabs(jet->eta) < 5.0 ) {
                  if (!leadingJet || jet->et > leadingJet->et) {
                    leadingJet = jet;
                  }
                  NJets++;
                } else {
                  if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
                }
              }
            }

            //******************************************************************************
            //soft muons
            //******************************************************************************
            NSoftMuons = 0;
            //soft muons
            for(Int_t l=0; l<muonArr->GetEntries(); l++) {
              const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[l]);
              Bool_t isCleanMuon = kFALSE;
              for (int k=0; k<leptonPt.size(); ++k) {
                if ( leptonType[k] == 13 
                     && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[k],leptonEta[k]) < 0.1
                  ) {
                  isCleanMuon = kTRUE; 
                  break;
                }
              }
              if ( fabs(mu->eta) < 2.4                   
                   && mu->pt > 3.0
                   && (mu->qualityBits & kTMLastStationAngTight)
                   && mu->nTkHits > 10
                   && fabs(mu->d0) < 0.2
                   && (mu->typeBits & kTracker)
                   && !isCleanMuon
                ) {
                NSoftMuons++;
              }
            }
            
            //***********************************************************************************************
            //|Z_vert-Z_l| maximum
            //***********************************************************************************************
            double zDiffMax = 0.0;
       
            double dz_i = 0;
            if (leptonType[0] == 11) {
              dz_i = ((mithep::TElectron*)((*electronArr)[leptonIndex[i]]))->dz;
            } else {
              dz_i = ((mithep::TMuon*)((*muonArr)[leptonIndex[i]]))->dz;
            }
            if (dz_i > zDiffMax) zDiffMax = dz_i;
    
            double dz_j;
            if (leptonType[j] == 11) {
              dz_j = ((mithep::TElectron*)((*electronArr)[leptonIndex[j]]))->dz;
            } else {
              dz_j = ((mithep::TMuon*)((*muonArr)[leptonIndex[j]]))->dz;
            }
            if (dz_j > zDiffMax) zDiffMax = dz_j;
            zDiffMax = fabs(dz_i - dz_j);

            //******************************************************************************
            //construct event variables
            //******************************************************************************
            mithep::FourVectorM lepton1;
            mithep::FourVectorM lepton2;
            if (leptonType[i] == 11) {
              lepton1.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], 0.51099892e-3 );
            } else {
              lepton1.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], 105.658369e-3 );
            }
            if (leptonType[j] == 11) {
              lepton2.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], 0.51099892e-3 );
            } else {
              lepton2.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], 105.658369e-3 );
            }
            mithep::FourVectorM dilepton = lepton1+lepton2;

            double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(lepton1.Phi(), 
                                                                 lepton2.Phi())* 180.0 / TMath::Pi();    
            double deltaPhiDileptonMet = mithep::MathUtils::DeltaPhi(met.Phi(), 
                                                                     dilepton.Phi())*180.0 / TMath::Pi();    
            double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt() * met.Phi()*
                                         (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

            //angle between MET and closest lepton
            double deltaPhiMetLepton[2] = {mithep::MathUtils::DeltaPhi(met.Phi(), lepton1.Phi()),
                                           mithep::MathUtils::DeltaPhi(met.Phi(), lepton2.Phi())};
  
            double mTW[2] = {TMath::Sqrt(2.0*lepton1.Pt()*met.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[0]))),
                             TMath::Sqrt(2.0*lepton2.Pt()*met.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[1])))};

            double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
              deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

            double METdeltaPhilEt = met.Pt();
            if(minDeltaPhiMetLepton < TMath::Pi()/2.)
              METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton);

    
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 10;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false};
  
            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 20.0) passCut[0] = true;
            
            if(zDiffMax < 1.0)                    passCut[1] = true;
  
            if(met.Pt()    > 20.0)               passCut[2] = true;
  
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
              if(METdeltaPhilEt > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(METdeltaPhilEt > 20) passCut[5] = true;
            }

            if(NJets     < 1)              passCut[6] = true;
        
            if (NSoftMuons == 0 )      passCut[7] = true;

            if (!(leptonPt.size() >= 3 && leptonPt[2] > 20.0)) passCut[8] = true;

            if(maxBtag < 2.1)                     passCut[9] = true;

//               if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] ) {
//               if (passCut[0] && dilepton.M() < 40.0) {
//                 cout << "Sim EVENT:\n";
//                 cout << "Lepton " <<  leptonPt[i] << " " << leptonEta[i] << " " << leptonPhi[i] << " " << leptonType[i] << " " << leptonCharge[i] << endl;
//                 cout << "Lepton " <<  leptonPt[j] << " " << leptonEta[j] << " " << leptonPhi[j] << " " << leptonType[j] << " " << leptonCharge[j] << endl;
//                 cout << fabs(leptonCharge[i] + leptonCharge[j]) << " " << ChargeSelection << endl;
//                 for(Int_t p=0; p<electronArr->GetEntries(); p++) {
//                   const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[p]);
//                   cout << "Electron " << p << " " << ele->pt << " " << ele->eta << " " << ele->phi 
//                        << endl;                                   
//                 }
//                 for(Int_t p=0; p<muonArr->GetEntries(); p++) {
//                   const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[p]);
//                   cout << "Muon " << p << " " << mu->pt << " " << mu->eta << " " << mu->phi << endl;
//                 }
//                 for(Int_t p=0; p<jetArr->GetEntries(); p++) {
//                   const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[p]);
//                   cout << "Jet " << p << " " << jet->et << " " << jet->eta << " " << jet->phi << " " << endl;
//                 }

//               }

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
            //Cut Selection Histograms
            fHWWSelectionDistributions_Simulation[p]->Fill(-1, eventWeight);
            if (finalState == 1 ) {
              fHWWToMuMuSelectionDistributions_Simulation[p]->Fill(-1, eventWeight);
            } else if(finalState == 0 ) {
              fHWWToEESelectionDistributions_Simulation[p]->Fill(-1, eventWeight);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelectionDistributions_Simulation[p]->Fill(-1, eventWeight);
            }
    
//             cout << nCuts << endl;
            for (int k=0;k<nCuts;k++) {
  //             cout << k << " : " << passCut[k] << endl;
              bool pass = true;
              bool passPreviousCut = true;
              for (int q=0;q<=k;q++) {
                pass = (pass && passCut[q]);
                if (q<k)
                  passPreviousCut = (passPreviousCut&& passCut[q]);
              }
      
 //              cout << pass << endl;
//              if (!passPreviousCut && pass) { cout << "STOP!!!\n"; assert(false);}
//               if (k==7) cout << "7: " << pass << " " << fHWWSelectionDistributions_Simulation[0]->GetBinContent(9) << endl;
//               if (k==8) cout << "8: " << pass << " " << fHWWSelectionDistributions_Simulation[0]->GetBinContent(10) << endl;
              if (pass) {
                fHWWSelectionDistributions_Simulation[p]->Fill(k, eventWeight);
                if (finalState == 1 ) {
                  fHWWToMuMuSelectionDistributions_Simulation[p]->Fill(k, eventWeight);
                } else if(finalState == 0 ) {
                  fHWWToEESelectionDistributions_Simulation[p]->Fill(k, eventWeight);
                } else if(finalState == 2 || finalState == 3 ) {
                  fHWWToEMuSelectionDistributions_Simulation[p]->Fill(k, eventWeight);
                }
              }
 //               cout << endl;
              
            }

//             cout << "check: " << fHWWToEMuSelectionDistributions_Simulation[0]->GetBinContent(9) << " " << fHWWToEMuSelectionDistributions_Simulation[0]->GetBinContent(10) << endl;


            //*****************************************************************************************
            //Make Preselection Histograms  
            //*****************************************************************************************
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] ) {
              fLeptonEtaDistributions_Simulation[p]->Fill(lepton1.Eta(), eventWeight); 
              fLeptonEtaDistributions_Simulation[p]->Fill(lepton2.Eta(), eventWeight);
              fLeptonPtMaxDistributions_Simulation[p]->Fill(lepton1.Pt(), eventWeight);
              fLeptonPtMinDistributions_Simulation[p]->Fill(lepton2.Pt(), eventWeight);
              fMetPtHistDistributions_Simulation[p]->Fill(met.Pt(), eventWeight);                             
              fMetPhiHistDistributions_Simulation[p]->Fill(met.Phi(), eventWeight);                            
              fDeltaPhiLeptonsDistributions_Simulation[p]->Fill(deltaPhiLeptons, eventWeight);


              fDileptonMassDistributions_Simulation[p]->Fill(dilepton.M(), eventWeight);
              if (finalState == 0) {
                fDileptonMass_eeDistributions_Simulation[p]->Fill(dilepton.M(), eventWeight);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              } else if (finalState == 1) {
                fDileptonMass_mumuDistributions_Simulation[p]->Fill(dilepton.M(), eventWeight);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              } else {
                fDileptonMass_emuDistributions_Simulation[p]->Fill(dilepton.M(), eventWeight);
                if (dilepton.M() > 60 && dilepton.M() < 120) {
                  if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                  } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                  } else {
                  }
                }
              }
            }

            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4]) {
              eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
                        lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);
            }



            //*********************************************************************************************
            //Plots after all Cuts
            //*********************************************************************************************
            if (passAllCuts) {
              fMinDeltaPhiLeptonMet_afterCutsDistributions_Simulation[p]->Fill(minDeltaPhiMetLepton, eventWeight);
              fMtLepton1_afterCutsDistributions_Simulation[p]->Fill(mTW[0], eventWeight);
              fMtLepton2_afterCutsDistributions_Simulation[p]->Fill(mTW[1], eventWeight);
              fMtHiggs_afterCutsDistributions_Simulation[p]->Fill(mtHiggs, eventWeight);
              fLeptonPtPlusMet_afterCutsDistributions_Simulation[p]->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt(), eventWeight);
    
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);

            }
          }
        }
      
      } //end loop over data           

      vector<double> ptbins;
      ptbins.push_back(10);  
      ptbins.push_back(12.5);  
      ptbins.push_back(15);  
      ptbins.push_back(20);  
      ptbins.push_back(25);  
      ptbins.push_back(30);  
      ptbins.push_back(50);  
      ptbins.push_back(80);  
      Int_t ErrorType = 2; //Clopper Pearson errors
      efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(numerator_Pt, denominator_Pt, "Efficiency_Pt", ptbins, ErrorType, -99, -99, 0, 1);
      
   
    }
  }
  

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;



//   //--------------------------------------------------------------------------------------------------------------
//   // Add systematic errors
//   //==============================================================================================================
//   for (int p=0; p < processNames.size(); ++p) {
      
//     for (int bin = 0; bin < fHWWSelectionDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fHWWSelectionDistributions[p]->GetBinError(bin),2) +
//                                          pow(fHWWSelectionDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fHWWSelectionDistributions[p]->GetBinError(bin) << ", " << fHWWSelectionDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fHWWSelectionDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fHWWSelectionDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fHWWToEESelectionDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fHWWToEESelectionDistributions[p]->GetBinError(bin),2) +
//                                          pow(fHWWToEESelectionDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fHWWToEESelectionDistributions[p]->GetBinError(bin) << ", " << fHWWToEESelectionDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fHWWToEESelectionDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fHWWToEESelectionDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fHWWToMuMuSelectionDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fHWWToMuMuSelectionDistributions[p]->GetBinError(bin),2) +
//                                          pow(fHWWToMuMuSelectionDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fHWWToMuMuSelectionDistributions[p]->GetBinError(bin) << ", " << fHWWToMuMuSelectionDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fHWWToMuMuSelectionDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fHWWToMuMuSelectionDistributions[p]->SetBinError(bin, totalError);
//     }
      
//     for (int bin = 0; bin < fHWWToEMuSelectionDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fHWWToEMuSelectionDistributions[p]->GetBinError(bin),2) +
//                                          pow(fHWWToEMuSelectionDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fHWWToEMuSelectionDistributions[p]->GetBinError(bin) << ", " << fHWWToEMuSelectionDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fHWWToEMuSelectionDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fHWWToEMuSelectionDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fLeptonEtaDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fLeptonEtaDistributions[p]->GetBinError(bin),2) +
//                                          pow(fLeptonEtaDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fLeptonEtaDistributions[p]->GetBinError(bin) << ", " << fLeptonEtaDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fLeptonEtaDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fLeptonEtaDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fLeptonPtMaxDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fLeptonPtMaxDistributions[p]->GetBinError(bin),2) +
//                                          pow(fLeptonPtMaxDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fLeptonPtMaxDistributions[p]->GetBinError(bin) << ", " << fLeptonPtMaxDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fLeptonPtMaxDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fLeptonPtMaxDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fLeptonPtMinDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fLeptonPtMinDistributions[p]->GetBinError(bin),2) +
//                                          pow(fLeptonPtMinDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fLeptonPtMinDistributions[p]->GetBinError(bin) << ", " << fLeptonPtMinDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fLeptonPtMinDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fLeptonPtMinDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fMetPtHistDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMetPtHistDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMetPtHistDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMetPtHistDistributions[p]->GetBinError(bin) << ", " << fMetPtHistDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMetPtHistDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMetPtHistDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fMetPhiHistDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMetPhiHistDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMetPhiHistDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMetPhiHistDistributions[p]->GetBinError(bin) << ", " << fMetPhiHistDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMetPhiHistDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMetPhiHistDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fDeltaPhiLeptonsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fDeltaPhiLeptonsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fDeltaPhiLeptonsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fDeltaPhiLeptonsDistributions[p]->GetBinError(bin) << ", " << fDeltaPhiLeptonsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fDeltaPhiLeptonsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fDeltaPhiLeptonsDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fDileptonMassDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fDileptonMassDistributions[p]->GetBinError(bin),2) +
//                                          pow(fDileptonMassDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fDileptonMassDistributions[p]->GetBinError(bin) << ", " << fDileptonMassDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fDileptonMassDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fDileptonMassDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fDileptonMass_eeDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fDileptonMass_eeDistributions[p]->GetBinError(bin),2) +
//                                          pow(fDileptonMass_eeDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fDileptonMass_eeDistributions[p]->GetBinError(bin) << ", " << fDileptonMass_eeDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fDileptonMass_eeDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fDileptonMass_eeDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fDileptonMass_mumuDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fDileptonMass_mumuDistributions[p]->GetBinError(bin),2) +
//                                          pow(fDileptonMass_mumuDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fDileptonMass_mumuDistributions[p]->GetBinError(bin) << ", " << fDileptonMass_mumuDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fDileptonMass_mumuDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fDileptonMass_mumuDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fDileptonMass_emuDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fDileptonMass_emuDistributions[p]->GetBinError(bin),2) +
//                                          pow(fDileptonMass_emuDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fDileptonMass_emuDistributions[p]->GetBinError(bin) << ", " << fDileptonMass_emuDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fDileptonMass_emuDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fDileptonMass_emuDistributions[p]->SetBinError(bin, totalError);
//     }

//     for (int bin = 0; bin < fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMinDeltaPhiLeptonMet_afterCutsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->GetBinError(bin) << ", " << fMinDeltaPhiLeptonMet_afterCutsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->SetBinError(bin, totalError);
//     }


//     for (int bin = 0; bin < fMtLepton1_afterCutsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMtLepton1_afterCutsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMtLepton1_afterCutsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMtLepton1_afterCutsDistributions[p]->GetBinError(bin) << ", " << fMtLepton1_afterCutsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMtLepton1_afterCutsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMtLepton1_afterCutsDistributions[p]->SetBinError(bin, totalError);
//     }


//     for (int bin = 0; bin < fMtLepton2_afterCutsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMtLepton2_afterCutsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMtLepton2_afterCutsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMtLepton2_afterCutsDistributions[p]->GetBinError(bin) << ", " << fMtLepton2_afterCutsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMtLepton2_afterCutsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMtLepton2_afterCutsDistributions[p]->SetBinError(bin, totalError);
//     }


//    for (int bin = 0; bin < fMtHiggs_afterCutsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fMtHiggs_afterCutsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fMtHiggs_afterCutsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fMtHiggs_afterCutsDistributions[p]->GetBinError(bin) << ", " << fMtHiggs_afterCutsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fMtHiggs_afterCutsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fMtHiggs_afterCutsDistributions[p]->SetBinError(bin, totalError);
//     }


//    for (int bin = 0; bin < fLeptonPtPlusMet_afterCutsDistributions[p]->GetNbinsX()+1; ++bin) {
//       Double_t totalError = TMath::Sqrt( pow(fLeptonPtPlusMet_afterCutsDistributions[p]->GetBinError(bin),2) +
//                                          pow(fLeptonPtPlusMet_afterCutsDistributions_systematicErr[p]->GetBinContent(bin), 2) );
// //       cout << "bin " << bin << " : " << totalError << ", " << fLeptonPtPlusMet_afterCutsDistributions[p]->GetBinError(bin) << ", " << fLeptonPtPlusMet_afterCutsDistributions_systematicErr[p]->GetBinContent(bin) << endl;
//       if (TMath::IsNaN(fLeptonPtPlusMet_afterCutsDistributions[p]->GetBinContent(bin))) cout << "NAN!!\n";
//       fLeptonPtPlusMet_afterCutsDistributions[p]->SetBinError(bin, totalError);
//     }

 
//   }


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TFile *file = new TFile("WWMCSelectionPlots.root", "RECREATE");

  for (UInt_t p=0;p<processNames.size();p++) {          

    DrawPrediction(fHWWSelectionDistributions[p], fHWWSelectionDistributions_Simulation[p], "HWWSelection_"+processLabels[p] , 0 , 0.3*LUMINOSITY, kFALSE, file );
    DrawPrediction(fHWWToEESelectionDistributions[p],fHWWToEESelectionDistributions_Simulation[p], "HWWToEESelection_"+processLabels[p] , 0 , 0.1*LUMINOSITY, kFALSE, file );
    DrawPrediction(fHWWToMuMuSelectionDistributions[p],fHWWToMuMuSelectionDistributions_Simulation[p],  "HWWToMuMuSelection_"+processLabels[p] , 0 , 0.1*LUMINOSITY, kFALSE, file );
    DrawPrediction(fHWWToEMuSelectionDistributions[p],fHWWToEMuSelectionDistributions_Simulation[p],  "HWWToEMuSelection_"+processLabels[p] , 0 , 0.1*LUMINOSITY, kFALSE, file );
    DrawPrediction(fLeptonEtaDistributions[p],fLeptonEtaDistributions_Simulation[p], "LeptonEta_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fLeptonPtMaxDistributions[p],fLeptonPtMaxDistributions_Simulation[p],  "LeptonPtMax_"+processLabels[p] , -999 , 0.05*LUMINOSITY, kFALSE, file );
    DrawPrediction(fLeptonPtMinDistributions[p],fLeptonPtMinDistributions_Simulation[p],  "LeptonPtMin_"+processLabels[p] , -999 , 0.05*LUMINOSITY, kFALSE, file );
    DrawPrediction(fMetPtHistDistributions[p],fMetPtHistDistributions_Simulation[p],  "MetPtHist_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fMetPhiHistDistributions[p],fMetPhiHistDistributions_Simulation[p],  "MetPhiHist_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fDeltaPhiLeptonsDistributions[p],fDeltaPhiLeptonsDistributions_Simulation[p],  "DeltaPhiLeptons_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fDileptonMassDistributions[p],fDileptonMassDistributions_Simulation[p],  "DileptonMass_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fDileptonMass_eeDistributions[p],fDileptonMass_eeDistributions_Simulation[p],  "DileptonMass_ee_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fDileptonMass_mumuDistributions[p],fDileptonMass_mumuDistributions_Simulation[p],  "DileptonMass_mumu_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fDileptonMass_emuDistributions[p],fDileptonMass_emuDistributions_Simulation[p],  "DileptonMass_emu_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fMinDeltaPhiLeptonMet_afterCutsDistributions[p],fMinDeltaPhiLeptonMet_afterCutsDistributions_Simulation[p],  "MinDeltaPhiLeptonMet_afterCuts_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fMtLepton1_afterCutsDistributions[p],fMtLepton1_afterCutsDistributions_Simulation[p],  "MtLepton1_afterCuts_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fMtLepton2_afterCutsDistributions[p],fMtLepton2_afterCutsDistributions_Simulation[p],  "MtLepton2_afterCuts_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fMtHiggs_afterCutsDistributions[p],fMtHiggs_afterCutsDistributions_Simulation[p],  "MtHiggs_afterCuts_"+processLabels[p] , -999 , -999, kFALSE, file );
    DrawPrediction(fLeptonPtPlusMet_afterCutsDistributions[p],fLeptonPtPlusMet_afterCutsDistributions_Simulation[p],  "LeptonPtPlusMet_afterCuts_"+processLabels[p] , -999 , -999, kFALSE, file );
    
  }


  cout << "check: " << fHWWToEMuSelectionDistributions_Simulation[0]->GetBinContent(9) << " " << fHWWToEMuSelectionDistributions_Simulation[0]->GetBinContent(11) << endl;


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  cout << "WW -> ee\n";

  for (UInt_t p=0;p<processNames.size();p++) {          
    cout.width(10);
    cout << left << processLabels[p];
  }
  cout << endl;

  for (int i=1; i < fHWWToEESelectionDistributions[0]->GetXaxis()->GetNbins()+1; ++i) {
    for (UInt_t p=0;p<processNames.size();p++) {          
      cout.width(10);
      cout << left << fHWWToEESelectionDistributions[p]->GetBinContent(i);      
      cout << " : ";
      cout << fHWWToEESelectionDistributions_Simulation[p]->GetBinContent(i) << " " ;      
    }
    cout << endl;
  }


  cout << "WW -> mm\n";

  for (UInt_t p=0;p<processNames.size();p++) {          
    cout.width(10);
    cout << left << processLabels[p];
  }
  cout << endl;

  for (int i=1; i < fHWWToMuMuSelectionDistributions[0]->GetXaxis()->GetNbins()+1; ++i) {
    for (UInt_t p=0;p<processNames.size();p++) {          
      cout.width(10);
      cout << left << fHWWToMuMuSelectionDistributions[p]->GetBinContent(i);      
      cout << " : ";
      cout << fHWWToMuMuSelectionDistributions_Simulation[p]->GetBinContent(i);      
    }
    cout << endl;
  }




  cout << "WW -> em\n";

  for (UInt_t p=0;p<processNames.size();p++) {          
    cout.width(10);
    cout << left << processLabels[p];
  }
  cout << endl;

  for (int i=1; i < fHWWToEMuSelectionDistributions[0]->GetXaxis()->GetNbins()+1; ++i) {
    for (UInt_t p=0;p<processNames.size();p++) {          
      cout.width(10);
      cout << left << fHWWToEMuSelectionDistributions[p]->GetBinContent(i);      
      cout << " : ";
      cout << fHWWToEMuSelectionDistributions_Simulation[p]->GetBinContent(i);      
    }
    cout << endl;
  }



  cout << "WW -> ll\n";

  for (UInt_t p=0;p<processNames.size();p++) {          
    cout.width(10);
    cout << left << processLabels[p];
  }
  cout << endl;

  for (int i=1; i < fHWWSelectionDistributions[0]->GetXaxis()->GetNbins()+1; ++i) {
    for (UInt_t p=0;p<processNames.size();p++) {          
      cout.width(10);
      cout << left << fHWWSelectionDistributions[p]->GetBinContent(i);      
      cout << " : ";
      cout << fHWWSelectionDistributions_Simulation[p]->GetBinContent(i);      
    }
    cout << endl;
  }

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  efficiency_pt->Draw("AP");
  cv->SaveAs("fakerate.gif");


        
  gBenchmark->Show("WWTemplate");       
} 


string DoubleToString(double i) {
  char temp[100];
  sprintf(temp, "%.4f", i);
  string str = temp;
  return str;
}


//*************************************************************************************************
//Draws the signal and background histograms together and makes gif file
//*************************************************************************************************
void DrawPrediction(TH1F* histPrediction,TH1F* histSimulation, string histname,
                    Double_t minY , Double_t maxY,
                    Bool_t useLogY, TFile *saveFile) {

  TCanvas *cv = new TCanvas(histname.c_str(), histname.c_str(), 0,0,800,600);
  if (useLogY) cv->SetLogy();

  TLegend *legend = new TLegend(0.73,0.55,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(histPrediction, "Prediction", "LP"); 
  legend->AddEntry(histSimulation, "Simulation", "F"); 

  Double_t tmpMaxY = histSimulation->GetMaximum();
  if (histPrediction->GetMaximum() > tmpMaxY) tmpMaxY = histPrediction->GetMaximum();
  tmpMaxY = tmpMaxY * 1.20;
  if (maxY != -999) tmpMaxY = maxY;  
  histSimulation->SetMaximum(tmpMaxY);  
  histPrediction->SetMaximum(tmpMaxY);
  
  Double_t tmpMinY = 0.0;
  if (minY != -999) tmpMinY = minY;
  histSimulation->SetMinimum(tmpMinY);  
  histPrediction->SetMinimum(tmpMinY);
  
  histSimulation->Draw("hist");
  histPrediction->Draw("E1,same");
  legend->Draw();

  //CMS Preliminary label
  TPaveText *prelimLabel = new TPaveText(0.21,0.85,0.41,0.90,"NDC");
  prelimLabel->SetTextColor(kBlack);
  prelimLabel->SetFillColor(kWhite);
  prelimLabel->SetBorderSize(0);
  prelimLabel->SetTextAlign(12);
  prelimLabel->SetTextSize(0.03);
  prelimLabel->AddText("CMS Preliminary 2010 #sqrt{s} = 7 TeV");
  prelimLabel->Draw();

  //Luminosity label
  TPaveText *tb = new TPaveText(0.21,0.77,0.41,0.82,"NDC");
  tb->SetTextColor(kBlack);
  tb->SetFillColor(kWhite);
  tb->SetBorderSize(0);
  tb->SetTextAlign(12);
  string lumi = DoubleToString(LUMINOSITY);
  tb->AddText((string("#int#font[12]{L}dt = ") + lumi + string(" pb^{ -1}")).c_str());
  tb->Draw();
  
  cv->SaveAs((histname+".gif").c_str());
  saveFile->WriteTObject(cv ,cv->GetName(), "WriteDelete");

}


//--------------------------------------------------------------------------------------------------
// Get Total Number of Events in the sample
//--------------------------------------------------------------------------------------------------
Double_t getNormalizationWeight(string filename, string datasetName) {
  // Get Normalization Weight Factor

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  if (!dir) {
    cout << "Could not find directory AnaFwkMod"
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }  
  Double_t NEvents = hist->Integral();

  //Get CrossSection
  Double_t Weight = 1.0;
  if (datasetName.find("data") == string::npos) {
    mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
    Double_t CrossSection = xstab.Get(datasetName.c_str());
    Weight = CrossSection / NEvents;
    cerr << datasetName << " : " << CrossSection << " " << NEvents << " " << Weight << endl;
  }

  delete dir;
  delete file;
  return Weight;


}




Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {

  Bool_t pass = kFALSE;
  if (isMuonData) {
    if ((runNum >= 136033) && (runNum <= 147116)) {
      if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
    } 
    if ((runNum >= 136033) && (runNum <= 139980)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140058) && (runNum <= 141882)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 141956) && (runNum <= 144114)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 146428) && (runNum <= 147116)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147196) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15) ) pass = kTRUE;
    }
    if ((runNum >= 147196) && (runNum <= 148058)) {
      if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 148819) && (runNum <= 149442)) {
      if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
    }
  } else {
    //it's electron data
    if ((runNum >= 136033) && (runNum <= 139980)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140058) && (runNum <= 141882)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 141956) && (runNum <= 144114)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 146428) && (runNum <= 147116)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147196) && (runNum <= 148058)) {
      if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 148819) && (runNum <= 149442)) {
      if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
    }
  }


  return pass;

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
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.10
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


Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
    return pass;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->HoverE < 0.15          
            && ele->sigiEtaiEta < 0.014          
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
             && ele->HoverE < 0.15          
             && ele->sigiEtaiEta < 0.034
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
          )
      ) {
      pass = kFALSE;
    }
  }

  return pass;
}

Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

       && (mu->nSeg > 1)
        && (mu->nPixHits > 0)
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
