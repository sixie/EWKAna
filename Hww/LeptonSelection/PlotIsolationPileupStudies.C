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


  vector<Double_t> RhoDensityMuons;
  vector<Double_t> caloIsoAverageMuonBarrel;
  vector<Double_t> ecalIsoAverageMuonBarrel;
  vector<Double_t> hcalIsoAverageMuonBarrel;
  vector<Double_t> trkIsoAverageMuonBarrel;
  vector<Double_t> caloIsoAverageMuonEndcap;
  vector<Double_t> ecalIsoAverageMuonEndcap;
  vector<Double_t> hcalIsoAverageMuonEndcap;
  vector<Double_t> trkIsoAverageMuonEndcap;
  vector<Double_t> caloIso05AverageMuonBarrel;
  vector<Double_t> ecalIso05AverageMuonBarrel;
  vector<Double_t> hcalIso05AverageMuonBarrel;
  vector<Double_t> trkIso05AverageMuonBarrel;
  vector<Double_t> caloIso05AverageMuonEndcap;
  vector<Double_t> ecalIso05AverageMuonEndcap;
  vector<Double_t> hcalIso05AverageMuonEndcap;
  vector<Double_t> trkIso05AverageMuonEndcap;


  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap;

 
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonEndcap;





  vector<Double_t> RhoDensityElectrons_Error;
  vector<Double_t> RhoDensityMuons_Error;
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
  vector<Double_t> caloIsoAverageMuonBarrel_Error;
  vector<Double_t> ecalIsoAverageMuonBarrel_Error;
  vector<Double_t> hcalIsoAverageMuonBarrel_Error;
  vector<Double_t> trkIsoAverageMuonBarrel_Error;
  vector<Double_t> caloIsoAverageMuonEndcap_Error;
  vector<Double_t> ecalIsoAverageMuonEndcap_Error;
  vector<Double_t> hcalIsoAverageMuonEndcap_Error;
  vector<Double_t> trkIsoAverageMuonEndcap_Error;
  vector<Double_t> caloIso05AverageMuonBarrel_Error;
  vector<Double_t> ecalIso05AverageMuonBarrel_Error;
  vector<Double_t> hcalIso05AverageMuonBarrel_Error;
  vector<Double_t> trkIso05AverageMuonBarrel_Error;
  vector<Double_t> caloIso05AverageMuonEndcap_Error;
  vector<Double_t> ecalIso05AverageMuonEndcap_Error;
  vector<Double_t> hcalIso05AverageMuonEndcap_Error;
  vector<Double_t> trkIso05AverageMuonEndcap_Error;


  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel_Error;
  vector<Double_t> ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> GammaIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;
  vector<Double_t> NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap_Error;

 
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonEndcap_Error;



  for (int n=0; n < 20; ++n) {     
    TCanvas *cv = new TCanvas("cv","cv", 800,600);

    TH1F *tmpRhoElectron = 0 ;
    TH1F *tmpRhoMuon = 0;
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
    TH1F *tmpMuon_caloIso_Barrel = 0;
    TH1F *tmpMuon_ecalIso_Barrel = 0;
    TH1F *tmpMuon_hcalIso_Barrel = 0;
    TH1F *tmpMuon_trkIso_Barrel = 0;
    TH1F *tmpMuon_caloIso05_Barrel = 0;
    TH1F *tmpMuon_ecalIso05_Barrel = 0;
    TH1F *tmpMuon_hcalIso05_Barrel = 0;
    TH1F *tmpMuon_trkIso05_Barrel = 0;
    TH1F *tmpMuon_caloIso_Endcap = 0;
    TH1F *tmpMuon_ecalIso_Endcap = 0;
    TH1F *tmpMuon_hcalIso_Endcap = 0;
    TH1F *tmpMuon_trkIso_Endcap = 0;
    TH1F *tmpMuon_caloIso05_Endcap = 0;
    TH1F *tmpMuon_ecalIso05_Endcap = 0;
    TH1F *tmpMuon_hcalIso05_Endcap = 0;
    TH1F *tmpMuon_trkIso05_Endcap = 0;


    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_GammaIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = 0;
 

 

    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = 0;


    
    cout << "here " << n << endl;

    if (Label == "JetData") {
      cout << label << " " << n << endl;
      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_Ele10Jet30").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_Mu10Jet30").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Ele10Jet30").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Mu10Jet30").c_str());
 



     } else if (Label == "WJetsMC") {
      cout << label << " " << n << endl;


      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_WJetsMC").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_WJetsMC").c_str());
 


    } else if (Label == "QCDMC") {
      cout << label << " " << n << endl;

      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Ele10Jet30").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_QCDMC_Mu10Jet30").c_str());
 

    } else if (Label == "ZMC") {
      cout << label << " " << n << endl;

      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_Zee").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_Zmm").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_Zmm").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zee").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm").c_str());
 




    } else if (Label == "HWW130") {
      cout << label << " " << n << endl;



      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_HWW130").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_HWW130").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_HWW130").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_HWW130").c_str());
 


    } else if (Label == "DataTagAndProbe") {
      cout << label << " " << n << endl;



      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_DataTagAndProbe").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_DataTagAndProbe").c_str());
 


    } else if (Label == "ZMCTagAndProbe") {
      cout << label << " " << n << endl;


      tmpRhoElectron = (TH1F*)f->Get((string("RhoElectron_") + IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_caloIso_Barrel = (TH1F*)f->Get((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_caloIso_Endcap = (TH1F*)f->Get((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ecalIso_Barrel = (TH1F*)f->Get((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ecalIso_Endcap = (TH1F*)f->Get((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_hcalIso_Barrel = (TH1F*)f->Get((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_hcalIso_Endcap = (TH1F*)f->Get((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_trkIso_Barrel = (TH1F*)f->Get((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_trkIso_Endcap = (TH1F*)f->Get((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_caloIso04_Barrel = (TH1F*)f->Get((string("Electron_caloIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_caloIso04_Endcap = (TH1F*)f->Get((string("Electron_caloIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ecalIso04_Barrel = (TH1F*)f->Get((string("Electron_ecalIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ecalIso04_Endcap = (TH1F*)f->Get((string("Electron_ecalIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_hcalIso04_Barrel = (TH1F*)f->Get((string("Electron_hcalIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_hcalIso04_Endcap = (TH1F*)f->Get((string("Electron_hcalIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_trkIso04_Barrel = (TH1F*)f->Get((string("Electron_trkIso04_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_trkIso04_Endcap = (TH1F*)f->Get((string("Electron_trkIso04_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());



      tmpElectron_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralHadronIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_GammaIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_GammaIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());


  

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str()); 
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_ZMCTagAndProbe").c_str());
 


    }


    cout << "load " << n << endl;
    assert( tmpRhoElectron );
    assert( tmpRhoMuon); 
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
    assert( tmpMuon_caloIso_Barrel);
    assert( tmpMuon_ecalIso_Barrel);
    assert( tmpMuon_hcalIso_Barrel);
    assert( tmpMuon_trkIso_Barrel);
    assert( tmpMuon_caloIso05_Barrel);
    assert( tmpMuon_ecalIso05_Barrel);
    assert( tmpMuon_hcalIso05_Barrel);
    assert( tmpMuon_trkIso05_Barrel);
    assert( tmpMuon_caloIso_Endcap);
    assert( tmpMuon_ecalIso_Endcap);
    assert( tmpMuon_hcalIso_Endcap);
    assert( tmpMuon_trkIso_Endcap);
    assert( tmpMuon_caloIso05_Endcap);
    assert( tmpMuon_ecalIso05_Endcap);
    assert( tmpMuon_hcalIso05_Endcap);
    assert( tmpMuon_trkIso05_Endcap);

   

    assert(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_GammaIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_GammaIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap);
    assert(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap);


  

    assert(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap);






    RhoDensityElectrons.push_back(tmpRhoElectron->GetMean());
    RhoDensityMuons.push_back(tmpRhoMuon->GetMean());
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
    caloIsoAverageMuonBarrel.push_back(tmpMuon_caloIso_Barrel->GetMean());
    ecalIsoAverageMuonBarrel.push_back(tmpMuon_ecalIso_Barrel->GetMean());
    hcalIsoAverageMuonBarrel.push_back(tmpMuon_hcalIso_Barrel->GetMean());
    trkIsoAverageMuonBarrel.push_back(tmpMuon_trkIso_Barrel->GetMean());
    caloIso05AverageMuonBarrel.push_back(tmpMuon_caloIso05_Barrel->GetMean());
    ecalIso05AverageMuonBarrel.push_back(tmpMuon_ecalIso05_Barrel->GetMean());
    hcalIso05AverageMuonBarrel.push_back(tmpMuon_hcalIso05_Barrel->GetMean());
    trkIso05AverageMuonBarrel.push_back(tmpMuon_trkIso05_Barrel->GetMean());
    caloIsoAverageMuonEndcap.push_back(tmpMuon_caloIso_Endcap->GetMean());
    ecalIsoAverageMuonEndcap.push_back(tmpMuon_ecalIso_Endcap->GetMean());
    hcalIsoAverageMuonEndcap.push_back(tmpMuon_hcalIso_Endcap->GetMean());
    trkIsoAverageMuonEndcap.push_back(tmpMuon_trkIso_Endcap->GetMean());
    caloIso05AverageMuonEndcap.push_back(tmpMuon_caloIso05_Endcap->GetMean());
    ecalIso05AverageMuonEndcap.push_back(tmpMuon_ecalIso05_Endcap->GetMean());
    hcalIso05AverageMuonEndcap.push_back(tmpMuon_hcalIso05_Endcap->GetMean());
    trkIso05AverageMuonEndcap.push_back(tmpMuon_trkIso05_Endcap->GetMean());



    ChargedIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel->GetMean());
    GammaIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone03_01ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap->GetMean());
    GammaIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone03_01ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel->GetMean());
    GammaIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone04_01ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap->GetMean());
    GammaIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone04_01ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel->GetMean());
    GammaIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone03_05ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap->GetMean());
    GammaIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone03_05ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel->GetMean());
    GammaIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone04_05ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap->GetMean());
    GammaIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone04_05ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel->GetMean());
    GammaIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone03_10ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap->GetMean());
    GammaIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone03_10ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel->GetMean());
    GammaIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone04_10ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap->GetMean());
    GammaIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone04_10ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel->GetMean());
    GammaIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone03_15ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap->GetMean());
    GammaIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone03_15ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralHadronIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel->GetMean());
    GammaIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel->GetMean());
    FPRemovedPFIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralFPRemovedPFIso_Cone04_15ThresholdAverageElectronBarrel.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralHadronIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap->GetMean());
    GammaIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap->GetMean());
    FPRemovedPFIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralFPRemovedPFIso_Cone04_15ThresholdAverageElectronEndcap.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMean());


    ChargedIso_Cone03_01ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_01ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralIso_Cone03_01ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_01ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedIso_Cone03_01ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_01ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralIso_Cone03_01ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_01ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedIso_Cone04_01ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_01ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralIso_Cone04_01ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_01ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedIso_Cone04_01ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_01ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralIso_Cone04_01ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_01ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedIso_Cone03_05ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_05ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralIso_Cone03_05ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_05ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedIso_Cone03_05ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_05ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralIso_Cone03_05ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_05ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedIso_Cone04_05ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_05ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralIso_Cone04_05ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_05ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedIso_Cone04_05ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_05ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralIso_Cone04_05ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_05ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedIso_Cone03_10ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_10ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralIso_Cone03_10ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_10ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedIso_Cone03_10ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_10ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralIso_Cone03_10ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_10ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedIso_Cone04_10ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_10ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralIso_Cone04_10ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_10ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedIso_Cone04_10ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_10ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralIso_Cone04_10ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_10ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedIso_Cone03_15ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone03_15ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralIso_Cone03_15ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel->GetMean());
    TotalPFIso_Cone03_15ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedIso_Cone03_15ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone03_15ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralIso_Cone03_15ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap->GetMean());
    TotalPFIso_Cone03_15ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedIso_Cone04_15ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedNoPUIso_Cone04_15ThresholdAverageMuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralIso_Cone04_15ThresholdAverageMuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel->GetMean());
    TotalPFIso_Cone04_15ThresholdAverageMuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedIso_Cone04_15ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap->GetMean());
    ChargedNoPUIso_Cone04_15ThresholdAverageMuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralIso_Cone04_15ThresholdAverageMuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap->GetMean());
    TotalPFIso_Cone04_15ThresholdAverageMuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap->GetMean());

 


    RhoDensityElectrons_Error.push_back(tmpRhoElectron->GetMeanError());
    RhoDensityMuons_Error.push_back(tmpRhoMuon->GetMeanError());
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
    caloIsoAverageMuonBarrel_Error.push_back(tmpMuon_caloIso_Barrel->GetMeanError());
    ecalIsoAverageMuonBarrel_Error.push_back(tmpMuon_ecalIso_Barrel->GetMeanError());
    hcalIsoAverageMuonBarrel_Error.push_back(tmpMuon_hcalIso_Barrel->GetMeanError());
    trkIsoAverageMuonBarrel_Error.push_back(tmpMuon_trkIso_Barrel->GetMeanError());    
    caloIso05AverageMuonBarrel_Error.push_back(tmpMuon_caloIso05_Barrel->GetMeanError());
    ecalIso05AverageMuonBarrel_Error.push_back(tmpMuon_ecalIso05_Barrel->GetMeanError());
    hcalIso05AverageMuonBarrel_Error.push_back(tmpMuon_hcalIso05_Barrel->GetMeanError());
    trkIso05AverageMuonBarrel_Error.push_back(tmpMuon_trkIso05_Barrel->GetMeanError()); 
    caloIsoAverageMuonEndcap_Error.push_back(tmpMuon_caloIso_Endcap->GetMeanError());
    ecalIsoAverageMuonEndcap_Error.push_back(tmpMuon_ecalIso_Endcap->GetMeanError());
    hcalIsoAverageMuonEndcap_Error.push_back(tmpMuon_hcalIso_Endcap->GetMeanError());
    trkIsoAverageMuonEndcap_Error.push_back(tmpMuon_trkIso_Endcap->GetMeanError());    
    caloIso05AverageMuonEndcap_Error.push_back(tmpMuon_caloIso05_Endcap->GetMeanError());
    ecalIso05AverageMuonEndcap_Error.push_back(tmpMuon_ecalIso05_Endcap->GetMeanError());
    hcalIso05AverageMuonEndcap_Error.push_back(tmpMuon_hcalIso05_Endcap->GetMeanError());
    trkIso05AverageMuonEndcap_Error.push_back(tmpMuon_trkIso05_Endcap->GetMeanError());    



    ChargedIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Barrel->GetMeanError());
    GammaIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_01Threshold_Endcap->GetMeanError());
    GammaIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_01Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Barrel->GetMeanError());
    GammaIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_01ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_01Threshold_Endcap->GetMeanError());
    GammaIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_01Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_01ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Barrel->GetMeanError());
    GammaIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_05Threshold_Endcap->GetMeanError());
    GammaIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_05Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Barrel->GetMeanError());
    GammaIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_05ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_05Threshold_Endcap->GetMeanError());
    GammaIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_05Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_05ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Barrel->GetMeanError());
    GammaIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_10Threshold_Endcap->GetMeanError());
    GammaIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_10Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Barrel->GetMeanError());
    GammaIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_10ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_10Threshold_Endcap->GetMeanError());
    GammaIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_10Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_10ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Barrel->GetMeanError());
    GammaIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone03_15Threshold_Endcap->GetMeanError());
    GammaIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone03_15Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone03_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralHadronIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Barrel->GetMeanError());
    GammaIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    FPRemovedPFIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_15ThresholdAverageElectronBarrel_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedIso_Cone04_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralHadronIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralHadronIso_Cone04_15Threshold_Endcap->GetMeanError());
    GammaIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_GammaIso_Cone04_15Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_TotalPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    FPRemovedPFIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_FPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralFPRemovedPFIso_Cone04_15ThresholdAverageElectronEndcap_Error.push_back(tmpElectron_NeutralFPRemovedPFIso_Cone04_15Threshold_Endcap->GetMeanError());






    ChargedIso_Cone03_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralIso_Cone03_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralIso_Cone03_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralIso_Cone04_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_01ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralIso_Cone04_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_01ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralIso_Cone03_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralIso_Cone03_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralIso_Cone04_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_05ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralIso_Cone04_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_05ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralIso_Cone03_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralIso_Cone03_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralIso_Cone04_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_10ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralIso_Cone04_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_10ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedIso_Cone03_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone03_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralIso_Cone03_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone03_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedIso_Cone03_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone03_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralIso_Cone03_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone03_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedIso_Cone04_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIso_Cone04_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralIso_Cone04_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel->GetMeanError());
    TotalPFIso_Cone04_15ThresholdAverageMuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedIso_Cone04_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIso_Cone04_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralIso_Cone04_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap->GetMeanError());
    TotalPFIso_Cone04_15ThresholdAverageMuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap->GetMeanError());


//     tmpMuon_caloIso->Draw("hist");
//     cv->SaveAs((string(tmpMuon_caloIso->GetName()) + ".gif").c_str());

  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  TCanvas *cv1 = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
  TPaveLabel *LabelText;

  const int nPoints = 20;
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.25;     
  }

  double y[nPoints];
  double yErrLow[nPoints];
  double yErrHigh[nPoints];


  TF1 *f1= 0;

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {
    
    y[i] = RhoDensityElectrons[i];
    yErrLow[i] = RhoDensityElectronsError[i];
    yErrHigh[i] = RhoDensityElectronsError[i];
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
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphRhoDensityElectrons" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIsoAverageElectronBarrel[i];
    yErrLow[i] = caloIsoAverageElectronBarrelError[i];
    yErrHigh[i] = caloIsoAverageElectronBarrelError[i];
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


  //*******************************************************************************************
  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIsoAverageElectronBarrel[i];
    yErrLow[i] = ecalIsoAverageElectronBarrelError[i];
    yErrHigh[i] = ecalIsoAverageElectronBarrelError[i];
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


  //*******************************************************************************************
  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIsoAverageElectronBarrel[i];
    yErrLow[i] = hcalIsoAverageElectronBarrelError[i];
    yErrHigh[i] = hcalIsoAverageElectronBarrelError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIsoAverageElectronBarrel[i];
    yErrLow[i] = trkIsoAverageElectronBarrelError[i];
    yErrHigh[i] = trkIsoAverageElectronBarrelError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIsoAverageElectronEndcap[i];
    yErrLow[i] = caloIsoAverageElectronEndcapError[i];
    yErrHigh[i] = caloIsoAverageElectronEndcapError[i];
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


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIsoAverageElectronEndcap[i];
    yErrLow[i] = ecalIsoAverageElectronEndcapError[i];
    yErrHigh[i] = ecalIsoAverageElectronEndcapError[i];
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

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIsoAverageElectronEndcap[i];
    yErrLow[i] = hcalIsoAverageElectronEndcapError[i];
    yErrHigh[i] = hcalIsoAverageElectronEndcapError[i];
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


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIsoAverageElectronEndcap[i];
    yErrLow[i] = trkIsoAverageElectronEndcapError[i];
    yErrHigh[i] = trkIsoAverageElectronEndcapError[i];
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


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIso04AverageElectronBarrel[i];
    yErrLow[i] = caloIso04AverageElectronBarrelError[i];
    yErrHigh[i] = caloIso04AverageElectronBarrelError[i];
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

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIso04AverageElectronBarrel[i];
    yErrLow[i] = ecalIso04AverageElectronBarrelError[i];
    yErrHigh[i] = ecalIso04AverageElectronBarrelError[i];
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

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIso04AverageElectronBarrel[i];
    yErrLow[i] = hcalIso04AverageElectronBarrelError[i];
    yErrHigh[i] = hcalIso04AverageElectronBarrelError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIso04AverageElectronBarrel[i];
    yErrLow[i] = trkIso04AverageElectronBarrelError[i];
    yErrHigh[i] = trkIso04AverageElectronBarrelError[i];
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


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIso04AverageElectronEndcap[i];
    yErrLow[i] = caloIso04AverageElectronEndcapError[i];
    yErrHigh[i] = caloIso04AverageElectronEndcapError[i];
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

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIso04AverageElectronEndcap[i];
    yErrLow[i] = ecalIso04AverageElectronEndcapError[i];
    yErrHigh[i] = ecalIso04AverageElectronEndcapError[i];
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

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIso04AverageElectronEndcap[i];
    yErrLow[i] = hcalIso04AverageElectronEndcapError[i];
    yErrHigh[i] = hcalIso04AverageElectronEndcapError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIso04AverageElectronEndcap[i];
    yErrLow[i] = trkIso04AverageElectronEndcapError[i];
    yErrHigh[i] = trkIso04AverageElectronEndcapError[i];
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










  //*******************************************************************************************
  // 01Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcapError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcapError[i];
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










  //*******************************************************************************************
  // 05Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcapError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcapError[i];
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








  //*******************************************************************************************
  // 10Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcapError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcapError[i];
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








  //*******************************************************************************************
  // 15Threshold
  //*******************************************************************************************



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcapError[i];
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



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrelError[i];
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








  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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




 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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

 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = GammaIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = FPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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


 //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap[i];
    yErrLow[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
    yErrHigh[i] = NeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcapError[i];
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

















  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = RhoDensityMuons[i];
    yErrLow[i] = RhoDensityMuonsError[i];
    yErrHigh[i] = RhoDensityMuonsError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphRhoDensityMuons = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphRhoDensityMuons->SetName(("GraphRhoDensityMuons"+label).c_str());
  GraphRhoDensityMuons->SetTitle("");
  GraphRhoDensityMuons->SetMarkerColor(kRed);
  GraphRhoDensityMuons->GetXaxis()->SetTitleOffset(1.02);
  GraphRhoDensityMuons->GetXaxis()->SetTitle("Number of Primary Vertices");
  GraphRhoDensityMuons->GetYaxis()->SetTitleOffset(1.05);
  GraphRhoDensityMuons->GetYaxis()->SetTitle("Average Rho [GeV]");
  GraphRhoDensityMuons->Draw("AP");

  f1 = new TF1(("RhoDensityMuonsFit"+label).c_str(), "pol1", 0.5, 9);
  GraphRhoDensityMuons->Fit(("RhoDensityMuonsFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphRhoDensityMuons" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIsoAverageMuonBarrel[i];
    yErrLow[i] = caloIsoAverageMuonBarrelError[i];
    yErrHigh[i] = caloIsoAverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIsoAverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIsoAverageMuonBarrel->SetName(("GraphCaloIsoAverageMuonBarrel"+label).c_str());
  GraphCaloIsoAverageMuonBarrel->SetTitle("");
  GraphCaloIsoAverageMuonBarrel->SetMarkerColor(kRed);
  GraphCaloIsoAverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIsoAverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIsoAverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("CaloIsoAverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIsoAverageMuonBarrel->Fit(("CaloIsoAverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIsoAverageMuonBarrel" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIsoAverageMuonBarrel[i];
    yErrLow[i] = ecalIsoAverageMuonBarrelError[i];
    yErrHigh[i] = ecalIsoAverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIsoAverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIsoAverageMuonBarrel->SetName(("GraphEcalIsoAverageMuonBarrel"+label).c_str());
  GraphEcalIsoAverageMuonBarrel->SetTitle("");
  GraphEcalIsoAverageMuonBarrel->SetMarkerColor(kRed);
  GraphEcalIsoAverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIsoAverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIsoAverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("EcalIsoAverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIsoAverageMuonBarrel->Fit(("EcalIsoAverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIsoAverageMuonBarrel" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIsoAverageMuonBarrel[i];
    yErrLow[i] = hcalIsoAverageMuonBarrelError[i];
    yErrHigh[i] = hcalIsoAverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIsoAverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIsoAverageMuonBarrel->SetName(("GraphHcalIsoAverageMuonBarrel"+label).c_str());
  GraphHcalIsoAverageMuonBarrel->SetTitle("");
  GraphHcalIsoAverageMuonBarrel->SetMarkerColor(kRed);
  GraphHcalIsoAverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIsoAverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIsoAverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("HcalIsoAverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIsoAverageMuonBarrel->Fit(("HcalIsoAverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIsoAverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIsoAverageMuonBarrel[i];
    yErrLow[i] = trkIsoAverageMuonBarrelError[i];
    yErrHigh[i] = trkIsoAverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIsoAverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIsoAverageMuonBarrel->SetName(("GraphTrkIsoAverageMuonBarrel"+label).c_str());
  GraphTrkIsoAverageMuonBarrel->SetTitle("");
  GraphTrkIsoAverageMuonBarrel->SetMarkerColor(kRed);
  GraphTrkIsoAverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIsoAverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIsoAverageMuonBarrel->Draw("AP");

  f1 = new TF1(("TrkIsoAverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIsoAverageMuonBarrel->Fit(("TrkIsoAverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIsoAverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIsoAverageMuonEndcap[i];
    yErrLow[i] = caloIsoAverageMuonEndcapError[i];
    yErrHigh[i] = caloIsoAverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIsoAverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIsoAverageMuonEndcap->SetName(("GraphCaloIsoAverageMuonEndcap"+label).c_str());
  GraphCaloIsoAverageMuonEndcap->SetTitle("");
  GraphCaloIsoAverageMuonEndcap->SetMarkerColor(kRed);
  GraphCaloIsoAverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIsoAverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIsoAverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("CaloIsoAverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIsoAverageMuonEndcap->Fit(("CaloIsoAverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIsoAverageMuonEndcap" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIsoAverageMuonEndcap[i];
    yErrLow[i] = ecalIsoAverageMuonEndcapError[i];
    yErrHigh[i] = ecalIsoAverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIsoAverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIsoAverageMuonEndcap->SetName(("GraphEcalIsoAverageMuonEndcap"+label).c_str());
  GraphEcalIsoAverageMuonEndcap->SetTitle("");
  GraphEcalIsoAverageMuonEndcap->SetMarkerColor(kRed);
  GraphEcalIsoAverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIsoAverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIsoAverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("EcalIsoAverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIsoAverageMuonEndcap->Fit(("EcalIsoAverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIsoAverageMuonEndcap" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIsoAverageMuonEndcap[i];
    yErrLow[i] = hcalIsoAverageMuonEndcapError[i];
    yErrHigh[i] = hcalIsoAverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIsoAverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIsoAverageMuonEndcap->SetName(("GraphHcalIsoAverageMuonEndcap"+label).c_str());
  GraphHcalIsoAverageMuonEndcap->SetTitle("");
  GraphHcalIsoAverageMuonEndcap->SetMarkerColor(kRed);
  GraphHcalIsoAverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIsoAverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIsoAverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("HcalIsoAverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIsoAverageMuonEndcap->Fit(("HcalIsoAverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIsoAverageMuonEndcap" + label + ".gif").c_str());




  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIsoAverageMuonEndcap[i];
    yErrLow[i] = trkIsoAverageMuonEndcapError[i];
    yErrHigh[i] = trkIsoAverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIsoAverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIsoAverageMuonEndcap->SetName(("GraphTrkIsoAverageMuonEndcap"+label).c_str());
  GraphTrkIsoAverageMuonEndcap->SetTitle("");
  GraphTrkIsoAverageMuonEndcap->SetMarkerColor(kRed);
  GraphTrkIsoAverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIsoAverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIsoAverageMuonEndcap->Draw("AP");

  f1 = new TF1(("TrkIsoAverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIsoAverageMuonEndcap->Fit(("TrkIsoAverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIsoAverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIso05AverageMuonBarrel[i];
    yErrLow[i] = caloIso05AverageMuonBarrelError[i];
    yErrHigh[i] = caloIso05AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIso05AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIso05AverageMuonBarrel->SetName(("GraphCaloIso05AverageMuonBarrel"+label).c_str());
  GraphCaloIso05AverageMuonBarrel->SetTitle("");
  GraphCaloIso05AverageMuonBarrel->SetMarkerColor(kRed);
  GraphCaloIso05AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIso05AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIso05AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("CaloIso05AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIso05AverageMuonBarrel->Fit(("CaloIso05AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIso05AverageMuonBarrel" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIso05AverageMuonBarrel[i];
    yErrLow[i] = ecalIso05AverageMuonBarrelError[i];
    yErrHigh[i] = ecalIso05AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIso05AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIso05AverageMuonBarrel->SetName(("GraphEcalIso05AverageMuonBarrel"+label).c_str());
  GraphEcalIso05AverageMuonBarrel->SetTitle("");
  GraphEcalIso05AverageMuonBarrel->SetMarkerColor(kRed);
  GraphEcalIso05AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIso05AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIso05AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("EcalIso05AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIso05AverageMuonBarrel->Fit(("EcalIso05AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIso05AverageMuonBarrel" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIso05AverageMuonBarrel[i];
    yErrLow[i] = hcalIso05AverageMuonBarrelError[i];
    yErrHigh[i] = hcalIso05AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIso05AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIso05AverageMuonBarrel->SetName(("GraphHcalIso05AverageMuonBarrel"+label).c_str());
  GraphHcalIso05AverageMuonBarrel->SetTitle("");
  GraphHcalIso05AverageMuonBarrel->SetMarkerColor(kRed);
  GraphHcalIso05AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIso05AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIso05AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("HcalIso05AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIso05AverageMuonBarrel->Fit(("HcalIso05AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIso05AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = caloIso05AverageMuonEndcap[i];
    yErrLow[i] = caloIso05AverageMuonEndcapError[i];
    yErrHigh[i] = caloIso05AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphCaloIso05AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphCaloIso05AverageMuonEndcap->SetName(("GraphCaloIso05AverageMuonEndcap"+label).c_str());
  GraphCaloIso05AverageMuonEndcap->SetTitle("");
  GraphCaloIso05AverageMuonEndcap->SetMarkerColor(kRed);
  GraphCaloIso05AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphCaloIso05AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphCaloIso05AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("CaloIso05AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphCaloIso05AverageMuonEndcap->Fit(("CaloIso05AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphCaloIso05AverageMuonEndcap" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ecalIso05AverageMuonEndcap[i];
    yErrLow[i] = ecalIso05AverageMuonEndcapError[i];
    yErrHigh[i] = ecalIso05AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphEcalIso05AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphEcalIso05AverageMuonEndcap->SetName(("GraphEcalIso05AverageMuonEndcap"+label).c_str());
  GraphEcalIso05AverageMuonEndcap->SetTitle("");
  GraphEcalIso05AverageMuonEndcap->SetMarkerColor(kRed);
  GraphEcalIso05AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphEcalIso05AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphEcalIso05AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("EcalIso05AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphEcalIso05AverageMuonEndcap->Fit(("EcalIso05AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphEcalIso05AverageMuonEndcap" + label + ".gif").c_str());

  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = hcalIso05AverageMuonEndcap[i];
    yErrLow[i] = hcalIso05AverageMuonEndcapError[i];
    yErrHigh[i] = hcalIso05AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphHcalIso05AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphHcalIso05AverageMuonEndcap->SetName(("GraphHcalIso05AverageMuonEndcap"+label).c_str());
  GraphHcalIso05AverageMuonEndcap->SetTitle("");
  GraphHcalIso05AverageMuonEndcap->SetMarkerColor(kRed);
  GraphHcalIso05AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphHcalIso05AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphHcalIso05AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("HcalIso05AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphHcalIso05AverageMuonEndcap->Fit(("HcalIso05AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphHcalIso05AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIso05AverageMuonBarrel[i];
    yErrLow[i] = trkIso05AverageMuonBarrelError[i];
    yErrHigh[i] = trkIso05AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIso05AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIso05AverageMuonBarrel->SetName(("GraphTrkIso05AverageMuonBarrel"+label).c_str());
  GraphTrkIso05AverageMuonBarrel->SetTitle("");
  GraphTrkIso05AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTrkIso05AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIso05AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIso05AverageMuonBarrel->Draw("AP");

  f1 = new TF1(("TrkIso05AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIso05AverageMuonBarrel->Fit(("TrkIso05AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIso05AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = trkIso05AverageMuonEndcap[i];
    yErrLow[i] = trkIso05AverageMuonEndcapError[i];
    yErrHigh[i] = trkIso05AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTrkIso05AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTrkIso05AverageMuonEndcap->SetName(("GraphTrkIso05AverageMuonEndcap"+label).c_str());
  GraphTrkIso05AverageMuonEndcap->SetTitle("");
  GraphTrkIso05AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTrkIso05AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTrkIso05AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTrkIso05AverageMuonEndcap->Draw("AP");

  f1 = new TF1(("TrkIso05AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTrkIso05AverageMuonEndcap->Fit(("TrkIso05AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTrkIso05AverageMuonEndcap" + label + ".gif").c_str());







 //*******************************************************************************************
  // 01Threshold
  //*******************************************************************************************




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_AverageMuonEndcap" + label + ".gif").c_str());



 //*******************************************************************************************
  // 05Threshold
  //*******************************************************************************************




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_AverageMuonEndcap" + label + ".gif").c_str());



 //*******************************************************************************************
  // 10Threshold
  //*******************************************************************************************




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_AverageMuonEndcap" + label + ".gif").c_str());



 //*******************************************************************************************
  // 15Threshold
  //*******************************************************************************************




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());




  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Fit(("ChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Fit(("NeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel->Fit(("TotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonBarrel" + label + ".gif").c_str());


  //*******************************************************************************************


  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Fit(("ChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Fit(("NeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap->Fit(("TotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_AverageMuonEndcap" + label + ".gif").c_str());





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

  file->WriteTObject(GraphRhoDensityMuons,GraphRhoDensityMuons->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIsoAverageMuonBarrel,GraphCaloIsoAverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIsoAverageMuonBarrel,GraphEcalIsoAverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIsoAverageMuonBarrel,GraphHcalIsoAverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIsoAverageMuonBarrel,GraphTrkIsoAverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIsoAverageMuonEndcap,GraphCaloIsoAverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIsoAverageMuonEndcap,GraphEcalIsoAverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIsoAverageMuonEndcap,GraphHcalIsoAverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIsoAverageMuonEndcap,GraphTrkIsoAverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIso05AverageMuonBarrel,GraphCaloIso05AverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIso05AverageMuonBarrel,GraphEcalIso05AverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIso05AverageMuonBarrel,GraphHcalIso05AverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIso05AverageMuonBarrel,GraphTrkIso05AverageMuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphCaloIso05AverageMuonEndcap,GraphCaloIso05AverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphEcalIso05AverageMuonEndcap,GraphEcalIso05AverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphHcalIso05AverageMuonEndcap,GraphHcalIso05AverageMuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTrkIso05AverageMuonEndcap,GraphTrkIso05AverageMuonEndcap->GetName(), "WriteDelete");



  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_01Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_05Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_10Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone03_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphChargedIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphGammaIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel,GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphChargedIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphChargedNoPUIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphNeutralHadronIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphGammaIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphTotalPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap,GraphNeutralFPRemovedPFIsoAverage_Cone04_15Threshold_ElectronEndcap->GetName(), "WriteDelete");



  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");

  file->Close();


}




void PlotIsolationPileupStudies() {


//    ComputeEffectiveArea("JetData");    
//    ComputeEffectiveArea("QCDMC");    
//   ComputeEffectiveArea("ZMC");
//    ComputeEffectiveArea("DataTagAndProbe");
//    ComputeEffectiveArea("ZMCTagAndProbe");
//    ComputeEffectiveArea("WJetsMC");
    
  ComputeEffectiveAreaElectrons("HWW130");


}
