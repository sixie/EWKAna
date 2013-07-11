//root -l EWKAna/Hww/LeptonSelection/PlotIsolationPileupStudies.C+
 
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 1D histograms
#include <TF1.h>                   // 1D histograms
#include <TPaveLabel.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <TGraphAsymmErrors.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLegend.h"
#include "TStyle.h"
#endif




//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string DoubleToString(double i) {
  char temp[100];
  sprintf(temp, "%.4f", i);
  string str = temp;
  return str;
}


void ComputeEffectiveAreaElectrons(string Label = "") {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
//   TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.bkg.root", "READ");
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.electrons.root", "READ");
  
  TLegend *legend = 0;
  
  vector<Double_t> RhoDensityElectrons;
  vector<Double_t> caloIsoAverageElectronBarrel;
  vector<Double_t> ecalIsoAverageElectronBarrel;
  vector<Double_t> hcalIsoAverageElectronBarrel;
  vector<Double_t> trkIsoAverageElectronBarrel;
  vector<Double_t> caloIsoAverageElectronEndcap;
  vector<Double_t> ecalIsoAverageElectronEndcap;
  vector<Double_t> hcalIsoAverageElectronEndcap;
  vector<Double_t> trkIsoAverageElectronEndcap;
  vector<Double_t> caloIso04AverageElectronBarrel;
  vector<Double_t> ecalIso04AverageElectronBarrel;
  vector<Double_t> hcalIso04AverageElectronBarrel;
  vector<Double_t> trkIso04AverageElectronBarrel;
  vector<Double_t> caloIso04AverageElectronEndcap;
  vector<Double_t> ecalIso04AverageElectronEndcap;
  vector<Double_t> hcalIso04AverageElectronEndcap;
  vector<Double_t> trkIso04AverageElectronEndcap;


 
  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;



 



  vector<Double_t> RhoDensityElectrons_Error;
  vector<Double_t> caloIsoAverageElectronBarrel_Error;
  vector<Double_t> ecalIsoAverageElectronBarrel_Error;
  vector<Double_t> hcalIsoAverageElectronBarrel_Error;
  vector<Double_t> trkIsoAverageElectronBarrel_Error;
  vector<Double_t> caloIsoAverageElectronEndcap_Error;
  vector<Double_t> ecalIsoAverageElectronEndcap_Error;
  vector<Double_t> hcalIsoAverageElectronEndcap_Error;
  vector<Double_t> trkIsoAverageElectronEndcap_Error;
  vector<Double_t> caloIso04AverageElectronBarrel_Error;
  vector<Double_t> ecalIso04AverageElectronBarrel_Error;
  vector<Double_t> hcalIso04AverageElectronBarrel_Error;
  vector<Double_t> trkIso04AverageElectronBarrel_Error;
  vector<Double_t> caloIso04AverageElectronEndcap_Error;
  vector<Double_t> ecalIso04AverageElectronEndcap_Error;
  vector<Double_t> hcalIso04AverageElectronEndcap_Error;
  vector<Double_t> trkIso04AverageElectronEndcap_Error;


  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;







  for (int n=0; n < 20; ++n) {     
    TCanvas *cv = new TCanvas("cv","cv", 800,600);

    TH1F *tmpRhoElectron = 0 ;
    TH1F *tmpElectron_caloIso_Barrel = 0;
    TH1F *tmpElectron_caloIso_Endcap = 0;
    TH1F *tmpElectron_ecalIso_Barrel = 0;
    TH1F *tmpElectron_ecalIso_Endcap = 0;
    TH1F *tmpElectron_hcalIso_Barrel = 0;
    TH1F *tmpElectron_hcalIso_Endcap = 0;
    TH1F *tmpElectron_trkIso_Barrel = 0;
    TH1F *tmpElectron_trkIso_Endcap = 0;
    TH1F *tmpElectron_caloIso04_Barrel = 0;
    TH1F *tmpElectron_caloIso04_Endcap = 0;
    TH1F *tmpElectron_ecalIso04_Barrel = 0;
    TH1F *tmpElectron_ecalIso04_Endcap = 0;
    TH1F *tmpElectron_hcalIso04_Barrel = 0;
    TH1F *tmpElectron_hcalIso04_Endcap = 0;
    TH1F *tmpElectron_trkIso04_Barrel = 0;
    TH1F *tmpElectron_trkIso04_Endcap = 0;

    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = 0;






 


    
    cout << "here " << n << endl;

//     if (Label == "JetData") {
//       cout << label << " " << n << endl;
//       tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());


//       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());


  

 


//      } else if (Label == "WJetsMC") {
//       cout << label << " " << n << endl;


//       tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
 

//       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());


  

 

//     } else if (Label == "QCDMC") {
//       cout << label << " " << n << endl;

//       tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//        tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());


//       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());


 

//     } else 
        if (Label == "ZMC") {
      cout << label << " " << n << endl;

      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
 
      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee_Pt10To20").c_str());








    } else if (Label == "HWW130_Pt10To20") {
      cout << label << " " << n << endl;



      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
 
       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
     tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());
      tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130_Pt10To20").c_str());





    } 

//         if (Label == "DataTagAndProbe") {
//       cout << label << " " << n << endl;



//       tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());

//       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());


 

//     } else if (Label == "ZMCTagAndProbe") {
//       cout << label << " " << n << endl;


//       tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());


//       tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
//       tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());


  

//     }


    cout << "load " << n << endl;
    assert( tmpRhoElectron );
    assert( tmpElectron_caloIso_Barrel); 
    assert( tmpElectron_caloIso_Endcap); 
    assert( tmpElectron_ecalIso_Barrel); 
    assert( tmpElectron_ecalIso_Endcap); 
    assert( tmpElectron_hcalIso_Barrel); 
    assert( tmpElectron_hcalIso_Endcap); 
    assert( tmpElectron_trkIso_Barrel);
    assert( tmpElectron_trkIso_Endcap); 
    assert( tmpElectron_caloIso04_Barrel); 
    assert( tmpElectron_caloIso04_Endcap); 
    assert( tmpElectron_ecalIso04_Barrel); 
    assert( tmpElectron_ecalIso04_Endcap); 
    assert( tmpElectron_hcalIso04_Barrel); 
    assert( tmpElectron_hcalIso04_Endcap); 
    assert( tmpElectron_trkIso04_Barrel);
    assert( tmpElectron_trkIso04_Endcap); 
   

    assert(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap);


 



    RhoDensityElectrons.push_back(tmpRhoElectron->GetMean());
    caloIsoAverageElectronBarrel.push_back(tmpElectron_caloIso_Barrel->GetMean());
    ecalIsoAverageElectronBarrel.push_back(tmpElectron_ecalIso_Barrel->GetMean());
    hcalIsoAverageElectronBarrel.push_back(tmpElectron_hcalIso_Barrel->GetMean());
    trkIsoAverageElectronBarrel.push_back(tmpElectron_trkIso_Barrel->GetMean());
    caloIsoAverageElectronEndcap.push_back(tmpElectron_caloIso_Endcap->GetMean());
    ecalIsoAverageElectronEndcap.push_back(tmpElectron_ecalIso_Endcap->GetMean());
    hcalIsoAverageElectronEndcap.push_back(tmpElectron_hcalIso_Endcap->GetMean());
    trkIsoAverageElectronEndcap.push_back(tmpElectron_trkIso_Endcap->GetMean());
    caloIso04AverageElectronBarrel.push_back(tmpElectron_caloIso04_Barrel->GetMean());
    ecalIso04AverageElectronBarrel.push_back(tmpElectron_ecalIso04_Barrel->GetMean());
    hcalIso04AverageElectronBarrel.push_back(tmpElectron_hcalIso04_Barrel->GetMean());
    trkIso04AverageElectronBarrel.push_back(tmpElectron_trkIso04_Barrel->GetMean());
    caloIso04AverageElectronEndcap.push_back(tmpElectron_caloIso04_Endcap->GetMean());
    ecalIso04AverageElectronEndcap.push_back(tmpElectron_ecalIso04_Endcap->GetMean());
    hcalIso04AverageElectronEndcap.push_back(tmpElectron_hcalIso04_Endcap->GetMean());
    trkIso04AverageElectronEndcap.push_back(tmpElectron_trkIso04_Endcap->GetMean());

    ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel->GetMean());
    GammaIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel->GetMean());
    FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap->GetMean());
    GammaIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap->GetMean());
    FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());
    VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap->GetMean());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());



 


    RhoDensityElectrons_Error.push_back(tmpRhoElectron->GetMeanError());
    caloIsoAverageElectronBarrel_Error.push_back(tmpElectron_caloIso_Barrel->GetMeanError());
    ecalIsoAverageElectronBarrel_Error.push_back(tmpElectron_ecalIso_Barrel->GetMeanError());
    hcalIsoAverageElectronBarrel_Error.push_back(tmpElectron_hcalIso_Barrel->GetMeanError());
    trkIsoAverageElectronBarrel_Error.push_back(tmpElectron_trkIso_Barrel->GetMeanError());
    caloIsoAverageElectronEndcap_Error.push_back(tmpElectron_caloIso_Endcap->GetMeanError());
    ecalIsoAverageElectronEndcap_Error.push_back(tmpElectron_ecalIso_Endcap->GetMeanError());
    hcalIsoAverageElectronEndcap_Error.push_back(tmpElectron_hcalIso_Endcap->GetMeanError());
    trkIsoAverageElectronEndcap_Error.push_back(tmpElectron_trkIso_Endcap->GetMeanError());
    caloIso04AverageElectronBarrel_Error.push_back(tmpElectron_caloIso04_Barrel->GetMeanError());
    ecalIso04AverageElectronBarrel_Error.push_back(tmpElectron_ecalIso04_Barrel->GetMeanError());
    hcalIso04AverageElectronBarrel_Error.push_back(tmpElectron_hcalIso04_Barrel->GetMeanError());
    trkIso04AverageElectronBarrel_Error.push_back(tmpElectron_trkIso04_Barrel->GetMeanError());
    caloIso04AverageElectronEndcap_Error.push_back(tmpElectron_caloIso04_Endcap->GetMeanError());
    ecalIso04AverageElectronEndcap_Error.push_back(tmpElectron_ecalIso04_Endcap->GetMeanError());
    hcalIso04AverageElectronEndcap_Error.push_back(tmpElectron_hcalIso04_Endcap->GetMeanError());
    trkIso04AverageElectronEndcap_Error.push_back(tmpElectron_trkIso04_Endcap->GetMeanError());


    ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel->GetMeanError());
    GammaIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap->GetMeanError());
    GammaIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedTotalPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error.push_back(tmpElectron_VertexSelectedNeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());








  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  ofstream outputfile("ElectronEffectiveArea.txt");
  TCanvas *cv1 = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
  TPaveLabel *LabelText;

  const int nPoints = 20;
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.25;     
  }

  double y[nPoints];
  double yErrLow[nPoints];
  double yErrHigh[nPoints];


  TF1 *f1= 0;

  //*******************************************************************************************

  Double_t ElectronRho = 0;

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {
    
    y[i] = RhoDensityElectrons[i];
    yErrLow[i] = RhoDensityElectrons_Error[i];
    yErrHigh[i] = RhoDensityElectrons_Error[i];
    NPileup[i] = i;
  }
  TGraphAsymmErrors *GraphRhoDensityElectrons = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphRhoDensityElectrons->SetName(("GraphRhoDensityElectrons"+label).c_str());
  GraphRhoDensityElectrons->SetTitle("");
  GraphRhoDensityElectrons->SetMarkerColor(kRed);
  GraphRhoDensityElectrons->GetXaxis()->SetTitleOffset(1.02);
  GraphRhoDensityElectrons->GetYaxis()->SetTitleOffset(1.05);
  GraphRhoDensityElectrons->Draw("AP");
  
  f1 = new TF1(("RhoDensityElectronsFit"+label).c_str(), "pol1", 0.5, 9);
  GraphRhoDensityElectrons->Fit(("RhoDensityElectronsFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  ElectronRho = f1->GetParameter(1);
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphRhoDensityElectrons" + label + ".gif").c_str());
  outputfile << "ElectronRho: " << ElectronRho << endl;

  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = caloIsoAverageElectronBarrel[i];
    yErrLow[i] = caloIsoAverageElectronBarrel_Error[i];
    yErrHigh[i] = caloIsoAverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIsoAverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIsoAverageElectronBarrel->SetName(("GraphCaloIsoAverageElectronBarrel"+label).c_str());
  GraphCaloIsoAverageElectronBarrel->SetTitle("");
  GraphCaloIsoAverageElectronBarrel->SetMarkerColor(kRed);
  GraphCaloIsoAverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIsoAverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIsoAverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("CaloIsoAverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIsoAverageElectronBarrel->Fit(("CaloIsoAverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIsoAverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphCaloIsoAverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************
  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ecalIsoAverageElectronBarrel[i];
    yErrLow[i] = ecalIsoAverageElectronBarrel_Error[i];
    yErrHigh[i] = ecalIsoAverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIsoAverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIsoAverageElectronBarrel->SetName(("GraphEcalIsoAverageElectronBarrel"+label).c_str());
  GraphEcalIsoAverageElectronBarrel->SetTitle("");
  GraphEcalIsoAverageElectronBarrel->SetMarkerColor(kRed);
  GraphEcalIsoAverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIsoAverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIsoAverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("EcalIsoAverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIsoAverageElectronBarrel->Fit(("EcalIsoAverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIsoAverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphEcalIsoAverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


  //*******************************************************************************************
  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = hcalIsoAverageElectronBarrel[i];
    yErrLow[i] = hcalIsoAverageElectronBarrel_Error[i];
    yErrHigh[i] = hcalIsoAverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIsoAverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIsoAverageElectronBarrel->SetName(("GraphHcalIsoAverageElectronBarrel"+label).c_str());
  GraphHcalIsoAverageElectronBarrel->SetTitle("");
  GraphHcalIsoAverageElectronBarrel->SetMarkerColor(kRed);
  GraphHcalIsoAverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIsoAverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIsoAverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("HcalIsoAverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIsoAverageElectronBarrel->Fit(("HcalIsoAverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIsoAverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphHcalIsoAverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = trkIsoAverageElectronBarrel[i];
    yErrLow[i] = trkIsoAverageElectronBarrel_Error[i];
    yErrHigh[i] = trkIsoAverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIsoAverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIsoAverageElectronBarrel->SetName(("GraphTrkIsoAverageElectronBarrel"+label).c_str());
  GraphTrkIsoAverageElectronBarrel->SetTitle("");
  GraphTrkIsoAverageElectronBarrel->SetMarkerColor(kRed);
  GraphTrkIsoAverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIsoAverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIsoAverageElectronBarrel->Draw("AP");

  f1 = new TF1(("TrkIsoAverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIsoAverageElectronBarrel->Fit(("TrkIsoAverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIsoAverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTrkIsoAverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = caloIsoAverageElectronEndcap[i];
    yErrLow[i] = caloIsoAverageElectronEndcap_Error[i];
    yErrHigh[i] = caloIsoAverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIsoAverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIsoAverageElectronEndcap->SetName(("GraphCaloIsoAverageElectronEndcap"+label).c_str());
  GraphCaloIsoAverageElectronEndcap->SetTitle("");
  GraphCaloIsoAverageElectronEndcap->SetMarkerColor(kRed);
  GraphCaloIsoAverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIsoAverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIsoAverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("CaloIsoAverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIsoAverageElectronEndcap->Fit(("CaloIsoAverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIsoAverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphCaloIsoAverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ecalIsoAverageElectronEndcap[i];
    yErrLow[i] = ecalIsoAverageElectronEndcap_Error[i];
    yErrHigh[i] = ecalIsoAverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIsoAverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIsoAverageElectronEndcap->SetName(("GraphEcalIsoAverageElectronEndcap"+label).c_str());
  GraphEcalIsoAverageElectronEndcap->SetTitle("");
  GraphEcalIsoAverageElectronEndcap->SetMarkerColor(kRed);
  GraphEcalIsoAverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIsoAverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIsoAverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("EcalIsoAverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIsoAverageElectronEndcap->Fit(("EcalIsoAverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIsoAverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphEcalIsoAverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = hcalIsoAverageElectronEndcap[i];
    yErrLow[i] = hcalIsoAverageElectronEndcap_Error[i];
    yErrHigh[i] = hcalIsoAverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIsoAverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIsoAverageElectronEndcap->SetName(("GraphHcalIsoAverageElectronEndcap"+label).c_str());
  GraphHcalIsoAverageElectronEndcap->SetTitle("");
  GraphHcalIsoAverageElectronEndcap->SetMarkerColor(kRed);
  GraphHcalIsoAverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIsoAverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIsoAverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("HcalIsoAverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIsoAverageElectronEndcap->Fit(("HcalIsoAverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIsoAverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphHcalIsoAverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = trkIsoAverageElectronEndcap[i];
    yErrLow[i] = trkIsoAverageElectronEndcap_Error[i];
    yErrHigh[i] = trkIsoAverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIsoAverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIsoAverageElectronEndcap->SetName(("GraphTrkIsoAverageElectronEndcap"+label).c_str());
  GraphTrkIsoAverageElectronEndcap->SetTitle("");
 GraphTrkIsoAverageElectronEndcap->SetMarkerColor(kRed);
  GraphTrkIsoAverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIsoAverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIsoAverageElectronEndcap->Draw("AP");

  f1 = new TF1(("TrkIsoAverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIsoAverageElectronEndcap->Fit(("TrkIsoAverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIsoAverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTrkIsoAverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************



  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = caloIso04AverageElectronBarrel[i];
    yErrLow[i] = caloIso04AverageElectronBarrel_Error[i];
    yErrHigh[i] = caloIso04AverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIso04AverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIso04AverageElectronBarrel->SetName(("GraphCaloIso04AverageElectronBarrel"+label).c_str());
  GraphCaloIso04AverageElectronBarrel->SetTitle("");
  GraphCaloIso04AverageElectronBarrel->SetMarkerColor(kRed);
  GraphCaloIso04AverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIso04AverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIso04AverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("CaloIso04AverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIso04AverageElectronBarrel->Fit(("CaloIso04AverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIso04AverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphCaloIso04AverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************



  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ecalIso04AverageElectronBarrel[i];
    yErrLow[i] = ecalIso04AverageElectronBarrel_Error[i];
    yErrHigh[i] = ecalIso04AverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIso04AverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIso04AverageElectronBarrel->SetName(("GraphEcalIso04AverageElectronBarrel"+label).c_str());
  GraphEcalIso04AverageElectronBarrel->SetTitle("");
  GraphEcalIso04AverageElectronBarrel->SetMarkerColor(kRed);
  GraphEcalIso04AverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIso04AverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIso04AverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("EcalIso04AverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIso04AverageElectronBarrel->Fit(("EcalIso04AverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIso04AverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphEcalIso04AverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************



  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = hcalIso04AverageElectronBarrel[i];
    yErrLow[i] = hcalIso04AverageElectronBarrel_Error[i];
    yErrHigh[i] = hcalIso04AverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIso04AverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIso04AverageElectronBarrel->SetName(("GraphHcalIso04AverageElectronBarrel"+label).c_str());
  GraphHcalIso04AverageElectronBarrel->SetTitle("");
  GraphHcalIso04AverageElectronBarrel->SetMarkerColor(kRed);
  GraphHcalIso04AverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIso04AverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIso04AverageElectronBarrel->Draw("AP");
  
  f1 = new TF1(("HcalIso04AverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIso04AverageElectronBarrel->Fit(("HcalIso04AverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIso04AverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphHcalIso04AverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = trkIso04AverageElectronBarrel[i];
    yErrLow[i] = trkIso04AverageElectronBarrel_Error[i];
    yErrHigh[i] = trkIso04AverageElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIso04AverageElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIso04AverageElectronBarrel->SetName(("GraphTrkIso04AverageElectronBarrel"+label).c_str());
  GraphTrkIso04AverageElectronBarrel->SetTitle("");
  GraphTrkIso04AverageElectronBarrel->SetMarkerColor(kRed);
  GraphTrkIso04AverageElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIso04AverageElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIso04AverageElectronBarrel->Draw("AP");

  f1 = new TF1(("TrkIso04AverageElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIso04AverageElectronBarrel->Fit(("TrkIso04AverageElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIso04AverageElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTrkIso04AverageElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = caloIso04AverageElectronEndcap[i];
    yErrLow[i] = caloIso04AverageElectronEndcap_Error[i];
    yErrHigh[i] = caloIso04AverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIso04AverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIso04AverageElectronEndcap->SetName(("GraphCaloIso04AverageElectronEndcap"+label).c_str());
  GraphCaloIso04AverageElectronEndcap->SetTitle("");
  GraphCaloIso04AverageElectronEndcap->SetMarkerColor(kRed);
  GraphCaloIso04AverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIso04AverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIso04AverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("CaloIso04AverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIso04AverageElectronEndcap->Fit(("CaloIso04AverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIso04AverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphCaloIso04AverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ecalIso04AverageElectronEndcap[i];
    yErrLow[i] = ecalIso04AverageElectronEndcap_Error[i];
    yErrHigh[i] = ecalIso04AverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIso04AverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIso04AverageElectronEndcap->SetName(("GraphEcalIso04AverageElectronEndcap"+label).c_str());
  GraphEcalIso04AverageElectronEndcap->SetTitle("");
  GraphEcalIso04AverageElectronEndcap->SetMarkerColor(kRed);
  GraphEcalIso04AverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIso04AverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIso04AverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("EcalIso04AverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIso04AverageElectronEndcap->Fit(("EcalIso04AverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIso04AverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphEcalIso04AverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = hcalIso04AverageElectronEndcap[i];
    yErrLow[i] = hcalIso04AverageElectronEndcap_Error[i];
    yErrHigh[i] = hcalIso04AverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIso04AverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIso04AverageElectronEndcap->SetName(("GraphHcalIso04AverageElectronEndcap"+label).c_str());
  GraphHcalIso04AverageElectronEndcap->SetTitle("");
  GraphHcalIso04AverageElectronEndcap->SetMarkerColor(kRed);
  GraphHcalIso04AverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIso04AverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIso04AverageElectronEndcap->Draw("AP");
  
  f1 = new TF1(("HcalIso04AverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIso04AverageElectronEndcap->Fit(("HcalIso04AverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIso04AverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphHcalIso04AverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = trkIso04AverageElectronEndcap[i];
    yErrLow[i] = trkIso04AverageElectronEndcap_Error[i];
    yErrHigh[i] = trkIso04AverageElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIso04AverageElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIso04AverageElectronEndcap->SetName(("GraphTrkIso04AverageElectronEndcap"+label).c_str());
  GraphTrkIso04AverageElectronEndcap->SetTitle("");
  GraphTrkIso04AverageElectronEndcap->SetMarkerColor(kRed);
  GraphTrkIso04AverageElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIso04AverageElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIso04AverageElectronEndcap->Draw("AP");

  f1 = new TF1(("TrkIso04AverageElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIso04AverageElectronEndcap->Fit(("TrkIso04AverageElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIso04AverageElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTrkIso04AverageElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;










  //*******************************************************************************************
  // 01Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


















  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;



























  //*******************************************************************************************
  // 05Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


















  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;
















  //*******************************************************************************************
  // 10Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


















  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;























  //*******************************************************************************************
  // 15Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


















  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("ChargedIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("GammaIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;





  //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("ChargedIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;

 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("GammaIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("GammaIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;




 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("VertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("VertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


 //*******************************************************************************************

  for (UInt_t i=0; i<UInt_t(nPoints); ++i) {    
    y[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    yErrHigh[i] = VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetName(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap"+label).c_str());
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetTitle("");
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->Fit(("VertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap" + label << " : " << f1->GetParameter(1) / ElectronRho << endl;


































































































  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************
  //*******************************************************************************************



















  TFile *file = new TFile("RhoPlots.root", "UPDATE");
  file->cd();

  file->WriteTObject(GraphRhoDensityElectrons,GraphRhoDensityElectrons->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIsoAverageElectronBarrel,GraphCaloIsoAverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIsoAverageElectronBarrel,GraphEcalIsoAverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIsoAverageElectronBarrel,GraphHcalIsoAverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIsoAverageElectronBarrel,GraphTrkIsoAverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIsoAverageElectronEndcap,GraphCaloIsoAverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIsoAverageElectronEndcap,GraphEcalIsoAverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIsoAverageElectronEndcap,GraphHcalIsoAverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIsoAverageElectronEndcap,GraphTrkIsoAverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIso04AverageElectronBarrel,GraphCaloIso04AverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIso04AverageElectronBarrel,GraphEcalIso04AverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIso04AverageElectronBarrel,GraphHcalIso04AverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIso04AverageElectronBarrel,GraphTrkIso04AverageElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIso04AverageElectronEndcap,GraphCaloIso04AverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIso04AverageElectronEndcap,GraphEcalIso04AverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIso04AverageElectronEndcap,GraphHcalIso04AverageElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIso04AverageElectronEndcap,GraphTrkIso04AverageElectronEndcap->GetName(), "WriteDelete");



  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphVertexSelectedTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphVertexSelectedFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphVertexSelectedNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");







  file->Close();


}




void ComputeElectronIsolationEffectiveArea() {


//   ComputeEffectiveAreaElectrons("HWW130_Pt10To20");
  ComputeEffectiveAreaElectrons("ZMC");


//    ComputeEffectiveArea("JetData");    
//    ComputeEffectiveArea("QCDMC");    
//    ComputeEffectiveArea("DataTagAndProbe");
//    ComputeEffectiveArea("ZMCTagAndProbe");
//    ComputeEffectiveArea("WJetsMC");
    


}
