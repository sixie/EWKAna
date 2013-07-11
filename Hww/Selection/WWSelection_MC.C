//root -l EWKAna/Hww/Selection/WWSelection_MC.C+\(\"\"\)
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

#endif

#define LUMINOSITY 1.0

//=== FUNCTION DECLARATIONS ======================================================================================

Double_t getNormalizationWeight(string filename, string datasetName);
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DrawDataSigBkgHistogram(THStack* stack, vector<string> processLabels, string histname, 
                             Double_t minY , Double_t maxY,
                             Bool_t useLogY, TFile *saveFile);

string DoubleToString(double i) {
  char temp[100];
  sprintf(temp, "%.4f", i);
  string str = temp;
  return str;
}

//=== MAIN MACRO =================================================================================================

void WWSelection_MC(const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  vector<vector<string> > processInputFiles;
    processInputFiles.push_back(vector<string>());
    processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-w1jets-0-100-v26_noskim.old.root");

//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ww2l-v26_noskim.root");
//   processInputFiles.push_back(vector<string>());
//   processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-ttbar2l-v26_noskim.root");
//    processInputFiles.push_back(vector<string>());
//    processInputFiles.back().push_back("/home/sixie/hist/WWAnalysis/WWAnalysis_p10-wjets-mg-v26_noskim.root");
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
   processNames.back().push_back("p10-w1jets-0-100-v26");

//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-ww2l-v26");
//   processNames.push_back(vector<string>());
//   processNames.back().push_back("p10-ttbar2l-v26");
//    processNames.push_back(vector<string>());
//    processNames.back().push_back("p10-wjets-mg-v26");
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
  processLabels.push_back("WJets");
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
    
  //--------------------------------------------------------------------------------------------------------------
  // Histograms
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
  
  cout << "here2\n";
  for (int i=0; i<processLabels.size(); ++i) {
    TH1F *tempHWWSelection = new TH1F(("HWWSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
   cout << "here211\n";
   tempHWWSelection->SetFillStyle(1001);
    tempHWWSelection->SetFillColor(processColors[i]);
   cout << "here212\n";
    tempHWWSelection->SetLineWidth(1);
    fHWWSelectionDistributions.push_back(tempHWWSelection);
   cout << "here213\n";
    fHWWSelectionStack->Add(tempHWWSelection);
   cout << "here214\n";
    TH1F *tempHWWToEESelection = new TH1F(("HWWToEESelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToEESelection->SetFillStyle(1001);
    tempHWWToEESelection->SetFillColor(processColors[i]);
    tempHWWToEESelection->SetLineWidth(1);
    fHWWToEESelectionDistributions.push_back(tempHWWToEESelection);
    fHWWToEESelectionStack->Add(tempHWWToEESelection);
    TH1F *tempHWWToMuMuSelection = new TH1F(("HWWToMuMuSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToMuMuSelection->SetFillStyle(1001);
    tempHWWToMuMuSelection->SetFillColor(processColors[i]);
    tempHWWToMuMuSelection->SetLineWidth(1);
    fHWWToMuMuSelectionDistributions.push_back(tempHWWToMuMuSelection);
    fHWWToMuMuSelectionStack->Add(tempHWWToMuMuSelection);
    TH1F *tempHWWToEMuSelection = new TH1F(("HWWToEMuSelectionStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempHWWToEMuSelection->SetFillStyle(1001);
    tempHWWToEMuSelection->SetFillColor(processColors[i]);
    tempHWWToEMuSelection->SetLineWidth(1);
    fHWWToEMuSelectionDistributions.push_back(tempHWWToEMuSelection);
    fHWWToEMuSelectionStack->Add(tempHWWToEMuSelection);
    
    TH1F *tempLeptonEta = new TH1F(("LeptonEtaStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempLeptonEta->SetFillStyle(1001);
    tempLeptonEta->SetFillColor(processColors[i]);
    tempLeptonEta->SetLineWidth(1);
    fLeptonEtaDistributions.push_back(tempLeptonEta);
    fLeptonEtaStack->Add(tempLeptonEta);
    TH1F *tempLeptonPtMax = new TH1F(("LeptonPtMaxStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempLeptonPtMax->SetFillStyle(1001);
    tempLeptonPtMax->SetFillColor(processColors[i]);
    tempLeptonPtMax->SetLineWidth(1);
    fLeptonPtMaxDistributions.push_back(tempLeptonPtMax);
    fLeptonPtMaxStack->Add(tempLeptonPtMax);
    TH1F *tempLeptonPtMin = new TH1F(("LeptonPtMinStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempLeptonPtMin->SetFillStyle(1001);
    tempLeptonPtMin->SetFillColor(processColors[i]);
    tempLeptonPtMin->SetLineWidth(1);
    fLeptonPtMinDistributions.push_back(tempLeptonPtMin);
    fLeptonPtMinStack->Add(tempLeptonPtMin);
    TH1F *tempMetPtHist = new TH1F(("MetPtHistStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMetPtHist->SetFillStyle(1001);
    tempMetPtHist->SetFillColor(processColors[i]);
    tempMetPtHist->SetLineWidth(1);
    fMetPtHistDistributions.push_back(tempMetPtHist);
    fMetPtHistStack->Add(tempMetPtHist);
    TH1F *tempMetPhiHist = new TH1F(("MetPhiHistStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMetPhiHist->SetFillStyle(1001);
    tempMetPhiHist->SetFillColor(processColors[i]);
    tempMetPhiHist->SetLineWidth(1);
    fMetPhiHistDistributions.push_back(tempMetPhiHist);
    fMetPhiHistStack->Add(tempMetPhiHist);
    TH1F *tempDeltaPhiLeptons = new TH1F(("DeltaPhiLeptonsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempDeltaPhiLeptons->SetFillStyle(1001);
    tempDeltaPhiLeptons->SetFillColor(processColors[i]);
    tempDeltaPhiLeptons->SetLineWidth(1);
    fDeltaPhiLeptonsDistributions.push_back(tempDeltaPhiLeptons);
    fDeltaPhiLeptonsStack->Add(tempDeltaPhiLeptons);
    TH1F *tempDileptonMass = new TH1F(("DileptonMassStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempDileptonMass->SetFillStyle(1001);
    tempDileptonMass->SetFillColor(processColors[i]);
    tempDileptonMass->SetLineWidth(1);
    fDileptonMassDistributions.push_back(tempDileptonMass);
    fDileptonMassStack->Add(tempDileptonMass);
    TH1F *tempDileptonMass_ee = new TH1F(("DileptonMass_eeStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempDileptonMass_ee->SetFillStyle(1001);
    tempDileptonMass_ee->SetFillColor(processColors[i]);
    tempDileptonMass_ee->SetLineWidth(1);
    fDileptonMass_eeDistributions.push_back(tempDileptonMass_ee);
    fDileptonMass_eeStack->Add(tempDileptonMass_ee);
    TH1F *tempDileptonMass_mumu = new TH1F(("DileptonMass_mumuStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempDileptonMass_mumu->SetFillStyle(1001);
    tempDileptonMass_mumu->SetFillColor(processColors[i]);
    tempDileptonMass_mumu->SetLineWidth(1);
    fDileptonMass_mumuDistributions.push_back(tempDileptonMass_mumu);
    fDileptonMass_mumuStack->Add(tempDileptonMass_mumu);
    TH1F *tempDileptonMass_emu = new TH1F(("DileptonMass_emuStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempDileptonMass_emu->SetFillStyle(1001);
    tempDileptonMass_emu->SetFillColor(processColors[i]);
    tempDileptonMass_emu->SetLineWidth(1);
    fDileptonMass_emuDistributions.push_back(tempDileptonMass_emu);
    fDileptonMass_emuStack->Add(tempDileptonMass_emu);


    TH1F *tempMinDeltaPhiLeptonMet_afterCuts = new TH1F(("MinDeltaPhiLeptonMet_afterCutsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMinDeltaPhiLeptonMet_afterCuts->SetFillStyle(1001);
    tempMinDeltaPhiLeptonMet_afterCuts->SetFillColor(processColors[i]);
    tempMinDeltaPhiLeptonMet_afterCuts->SetLineWidth(1);
    fMinDeltaPhiLeptonMet_afterCutsDistributions.push_back(tempMinDeltaPhiLeptonMet_afterCuts);
    fMinDeltaPhiLeptonMet_afterCutsStack->Add(tempMinDeltaPhiLeptonMet_afterCuts);
    TH1F *tempMtLepton1_afterCuts = new TH1F(("MtLepton1_afterCutsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMtLepton1_afterCuts->SetFillStyle(1001);
    tempMtLepton1_afterCuts->SetFillColor(processColors[i]);
    tempMtLepton1_afterCuts->SetLineWidth(1);
    fMtLepton1_afterCutsDistributions.push_back(tempMtLepton1_afterCuts);
    fMtLepton1_afterCutsStack->Add(tempMtLepton1_afterCuts);
    TH1F *tempMtLepton2_afterCuts = new TH1F(("MtLepton2_afterCutsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMtLepton2_afterCuts->SetFillStyle(1001);
    tempMtLepton2_afterCuts->SetFillColor(processColors[i]);
    tempMtLepton2_afterCuts->SetLineWidth(1);
    fMtLepton2_afterCutsDistributions.push_back(tempMtLepton2_afterCuts);
    fMtLepton2_afterCutsStack->Add(tempMtLepton2_afterCuts);
    TH1F *tempMtHiggs_afterCuts = new TH1F(("MtHiggs_afterCutsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempMtHiggs_afterCuts->SetFillStyle(1001);
    tempMtHiggs_afterCuts->SetFillColor(processColors[i]);
    tempMtHiggs_afterCuts->SetLineWidth(1);
    fMtHiggs_afterCutsDistributions.push_back(tempMtHiggs_afterCuts);
    fMtHiggs_afterCutsStack->Add(tempMtHiggs_afterCuts);
    TH1F *tempLeptonPtPlusMet_afterCuts = new TH1F(("LeptonPtPlusMet_afterCutsStack_"+processLabels[i]).c_str(),";Cut Number;Number of Events", 9, -1.5, 7.5);    
    tempLeptonPtPlusMet_afterCuts->SetFillStyle(1001);
    tempLeptonPtPlusMet_afterCuts->SetFillColor(processColors[i]);
    tempLeptonPtPlusMet_afterCuts->SetLineWidth(1);
    fLeptonPtPlusMet_afterCutsDistributions.push_back(tempLeptonPtPlusMet_afterCuts);
    fLeptonPtPlusMet_afterCutsStack->Add(tempLeptonPtPlusMet_afterCuts);

  }
  cout << "here3\n";
 

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
      inputFile = new TFile(processInputFiles[p][f].c_str());
      assert(inputFile);

      //********************************************************
      // Calculate event weight
      //********************************************************
      cout << processNames[p][f] << endl;
      Double_t eventWeight = getNormalizationWeight(processInputFiles[p][f], processNames[p][f]) * LUMINOSITY;
      cout << eventWeight << endl;


      //********************************************************
      // Get Tree
      //********************************************************
      eventTree = (TTree*)inputFile->Get("Events"); assert(eventTree);

      //*****************************************************************************************
      //Loop over Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
		
//       if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;
        if (ientry % 100000 == 0) {
          cout << "Event " << ientry << endl;
        }
  
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
        //soft muons
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
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

    
        for(Int_t i=0; i<jetArr->GetEntries(); i++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
              leptonOverlap = kTRUE;
            }
          }

          if (!(jet->et > 25 && fabs(jet->eta) < 5.0 && !leptonOverlap)) continue;
          if (!leadingJet || jet->et > leadingJet->et) {
            leadingJet = jet;
          }
          NJets++;
        }


        //******************************************************************************
        //dilepton preselection
        //******************************************************************************
        if (leptonPt.size() < 2) continue;
        if (!(leptonPt[0] > 20.0 && leptonPt[1] > 20.0)) continue;

        for(int i = 0; i < leptonPt.size(); ++i) {
          for(int j = i+1; j < leptonPt.size(); ++j) {

            //require opposite sign
            if (leptonCharge[i] == leptonCharge[j]) continue;


            Int_t finalState = -1;
            if (leptonType[i] == 11 && leptonType[j] == 11) {
              finalState = 0;
            } else if (leptonType[i] == 13 && leptonType[j] == 13) {
              finalState = 1;
            } else if (leptonType[i] == 11 && leptonType[j] == 13) {
              finalState = 2;
            } else if (leptonType[i] == 13 && leptonType[j] == 11) {
              finalState = 3;
            }




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
            const int nCuts = 8;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false};
  
            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 20.0) passCut[0] = true;
    
            if(met.Pt()    > 20.0)               passCut[1] = true;
  
            if(dilepton.M() > 12.0)            passCut[2] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[3] = true;
              if(METdeltaPhilEt > 35) passCut[4] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[3] = true;
              if(METdeltaPhilEt > 20) passCut[4] = true;
            }

            if(NJets     < 1)              passCut[5] = true;
        
            if (NSoftMuons == 0 )      passCut[6] = true;

            if (!(leptonPt.size() >= 3 && leptonPt[2] > 20.0)) passCut[7] = true;


            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
            //Cut Selection Histograms
            fHWWSelectionDistributions[p]->Fill(-1);
            if (finalState == 1 ) {
              fHWWToMuMuSelectionDistributions[p]->Fill(-1);
            } else if(finalState == 0 ) {
              fHWWToEESelectionDistributions[p]->Fill(-1);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelectionDistributions[p]->Fill(-1);
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
                fHWWSelectionDistributions[p]->Fill(k);
                if (finalState == 1 ) {
                  fHWWToMuMuSelectionDistributions[p]->Fill(k);
                } else if(finalState == 0 ) {
                  fHWWToEESelectionDistributions[p]->Fill(k);
                } else if(finalState == 2 || finalState == 3 ) {
                  fHWWToEMuSelectionDistributions[p]->Fill(k);
                }
              }
            }

            //*****************************************************************************************
            //Make Preselection Histograms  
            //*****************************************************************************************
            fLeptonEtaDistributions[p]->Fill(lepton1.Eta()); 
            fLeptonEtaDistributions[p]->Fill(lepton2.Eta());
            fLeptonPtMaxDistributions[p]->Fill(lepton1.Pt());
            fLeptonPtMinDistributions[p]->Fill(lepton2.Pt());
            fMetPtHistDistributions[p]->Fill(met.Pt());                             
            fMetPhiHistDistributions[p]->Fill(met.Phi());                            
            fDeltaPhiLeptonsDistributions[p]->Fill(deltaPhiLeptons);


            fDileptonMassDistributions[p]->Fill(dilepton.M());
            if (finalState == 0) {
              fDileptonMass_eeDistributions[p]->Fill(dilepton.M());
              if (dilepton.M() > 60 && dilepton.M() < 120) {
                if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                } else {
                }
              }
            } else if (finalState == 1) {
              fDileptonMass_mumuDistributions[p]->Fill(dilepton.M());
              if (dilepton.M() > 60 && dilepton.M() < 120) {
                if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                } else {
                }
              }
            } else {
              fDileptonMass_emuDistributions[p]->Fill(dilepton.M());
              if (dilepton.M() > 60 && dilepton.M() < 120) {
                if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                } else {
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
              fMinDeltaPhiLeptonMet_afterCutsDistributions[p]->Fill(minDeltaPhiMetLepton);
              fMtLepton1_afterCutsDistributions[p]->Fill(mTW[0]);
              fMtLepton2_afterCutsDistributions[p]->Fill(mTW[1]);
              fMtHiggs_afterCutsDistributions[p]->Fill(mtHiggs);
              fLeptonPtPlusMet_afterCutsDistributions[p]->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt());
    
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);

            }
          }
        }
      } //end loop over data     
    
    }
  }

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TFile *file = new TFile("WWMCSelectionPlots.root", "RECREATE");

  DrawDataSigBkgHistogram(fHWWSelectionStack, processLabels, "HWWSelection" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fHWWToEESelectionStack,processLabels,  "HWWToEESelection" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fHWWToMuMuSelectionStack,processLabels,  "HWWToMuMuSelection" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fHWWToEMuSelectionStack,processLabels,  "HWWToEMuSelection" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fLeptonEtaStack,processLabels,  "LeptonEta" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fLeptonPtMaxStack,processLabels,  "LeptonPtMax" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fLeptonPtMinStack, processLabels, "LeptonPtMin" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMetPtHistStack,processLabels,  "MetPtHist" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMetPhiHistStack, processLabels, "MetPhiHist" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fDeltaPhiLeptonsStack,processLabels,  "DeltaPhiLeptons" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fDileptonMassStack,processLabels,  "DileptonMass" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fDileptonMass_eeStack,processLabels,  "DileptonMass_ee" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fDileptonMass_mumuStack,processLabels,  "DileptonMass_mumu" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fDileptonMass_emuStack,processLabels,  "DileptonMass_emu" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMinDeltaPhiLeptonMet_afterCutsStack,processLabels,  "MinDeltaPhiLeptonMet_afterCuts" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMtLepton1_afterCutsStack, processLabels, "MtLepton1_afterCuts" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMtLepton2_afterCutsStack,processLabels,  "MtLepton2_afterCuts" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fMtHiggs_afterCutsStack,processLabels,  "MtHiggs_afterCuts" , -999 , -999, kFALSE, file );
  DrawDataSigBkgHistogram(fLeptonPtPlusMet_afterCutsStack,processLabels,  "LeptonPtPlusMet_afterCuts" , -999 , -999, kFALSE, file );



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
    }
    cout << endl;
  }





        
  gBenchmark->Show("WWTemplate");       
} 


//*************************************************************************************************
//Draws the signal and background histograms together and makes gif file
//*************************************************************************************************
void DrawDataSigBkgHistogram(THStack* stack, vector<string> processLabels, string histname,
                             Double_t minY , Double_t maxY,
                             Bool_t useLogY, TFile *saveFile) {

  TCanvas *cv = new TCanvas(histname.c_str(), histname.c_str(), 0,0,800,600);
  if (useLogY) cv->SetLogy();

  TLegend *legend = new TLegend(0.73,0.55,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  for (int j=0;j<processLabels.size();++j) {
    legend->AddEntry(stack->GetHists()->At(j), processLabels[j].c_str(), "F"); 
  }

  stack->Draw();
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
  tb->AddText((string("#int#font[12]{L}dt = ") + lumi + string(" nb^{ -1}")).c_str());
  tb->Draw();
  
  cv->SaveAs((histname+".gif").c_str());
  saveFile->WriteTObject(cv ,cv->GetName(), "WriteDelete");

}



//--------------------------------------------------------------------------------------------------
// Get Total Number of Events in the sample
//--------------------------------------------------------------------------------------------------
Double_t getNormalizationWeight(string filename, string datasetName) {
//   // Get Normalization Weight Factor

//   //Get Number of Events in the Sample
//   TFile *file = new TFile(filename.c_str(),"READ");
//   if (!file) {
//     cout << "Could not open file " << filename << endl;
//     return 0;
//   }

//   TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
//   if (!dir) {
//     cout << "Could not find directory AnaFwkMod"
//          << " in file " << filename << endl;
//     delete file;
//     return 0;
//   }

//   TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
//   if (!hist) {
//     cout << "Could not find histogram hDEvents in directory AnaFwkMod"
//          << " in file " << filename << endl;
//     delete dir;
//     delete file;
//     return 0;
//   }  
//   Double_t NEvents = hist->Integral();

//   //Get CrossSection
//   Double_t Weight = 1.0;
//   if (datasetName.find("data") == string::npos) {
//     mithep::SimpleTable xstab("$CMSSW_BASE/src/MitPhysics/data/xs.dat");
//     Double_t CrossSection = xstab.Get(datasetName.c_str());
//     Weight = CrossSection / NEvents;
//     cerr << datasetName << " : " << CrossSection << " " << NEvents << " " << Weight << endl;
//   }

//   delete dir;
//   delete file;
//   return Weight;

  return 1.0;
}




Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;

  if (isMuonData) {
    if ((runNum >= 132440) && (runNum <= 147119)) {
      if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
    } 
    if ((runNum >= 132440) && (runNum <= 135058)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 135059) && (runNum <= 140401)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140042) && (runNum <= 141900)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 141901) && (runNum <= 146427)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 146428) && (runNum <= 147119)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15_v1) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
    }
  } else {
    //it's electron data
    if ((runNum >= 132440) && (runNum <= 135058)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 135059) && (runNum <= 140041)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140042) && (runNum <= 141900)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 141901) && (runNum <= 146427)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 146428) && (runNum <= 147119)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( !(triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
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
            //&& ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
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
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.06
             //&& ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
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
            && ele->sigiEtaiEta < 0.014          
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.5
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.034
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.5
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
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

//         && (mu->nSeg >= 2)
//         && (mu->nValidHits >= 1)
//         && (mu->nPixHits >= 1)
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
