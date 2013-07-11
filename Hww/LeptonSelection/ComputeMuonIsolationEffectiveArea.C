 
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


void ComputeEffectiveAreaMuons(string Label = "") {

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
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.muons.root", "READ");
  
  TLegend *legend = 0;
  
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



  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_01Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_01Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_05Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_05Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_10Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_10Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_15Threshold_MuonBarrel;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonEndcap;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_15Threshold_MuonEndcap;




  vector<Double_t> RhoDensityMuons_Error;
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

  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_01Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_01Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_05Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_05Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_10Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_10Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_03Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_15Threshold_MuonBarrel_Error;
  vector<Double_t> ChargedIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> ChargedNoPUIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> NeutralIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> TotalPFIsoAverage_04Cone_15Threshold_MuonEndcap_Error;
  vector<Double_t> VertexSelectedPFIsoAverage_04Cone_15Threshold_MuonEndcap_Error;



  for (int n=0; n < 20; ++n) {     
    TCanvas *cv = new TCanvas("cv","cv", 800,600);

    TH1F *tmpRhoMuon = 0;
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



    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel = 0;
    TH1F *tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = 0;
    TH1F *tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap = 0;




 


    
    cout << "here " << n << endl;


 


    if (Label == "ZMC") {
      cout << label << " " << n << endl;

      tmpRhoMuon = (TH1F*)f->Get((string("RhoMuon_") + IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_caloIso_Barrel = (TH1F*)f->Get((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ecalIso_Barrel = (TH1F*)f->Get((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_hcalIso_Barrel = (TH1F*)f->Get((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_trkIso_Barrel = (TH1F*)f->Get((string("Muon_trkIso_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_caloIso05_Barrel = (TH1F*)f->Get((string("Muon_caloIso05_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ecalIso05_Barrel = (TH1F*)f->Get((string("Muon_ecalIso05_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_hcalIso05_Barrel = (TH1F*)f->Get((string("Muon_hcalIso05_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_trkIso05_Barrel = (TH1F*)f->Get((string("Muon_trkIso05_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_caloIso_Endcap = (TH1F*)f->Get((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ecalIso_Endcap = (TH1F*)f->Get((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_hcalIso_Endcap = (TH1F*)f->Get((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_trkIso_Endcap = (TH1F*)f->Get((string("Muon_trkIso_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_caloIso05_Endcap = (TH1F*)f->Get((string("Muon_caloIso05_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ecalIso05_Endcap = (TH1F*)f->Get((string("Muon_ecalIso05_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_hcalIso05_Endcap = (TH1F*)f->Get((string("Muon_hcalIso05_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_trkIso05_Endcap = (TH1F*)f->Get((string("Muon_trkIso05_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());

      tmpMuon_ChargedIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_01Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_01Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_05Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_05Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_10Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_10Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone03_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_15Threshold_Barrel_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_ChargedNoPUIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_NeutralIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_NeutralIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
      tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFIso_Cone04_15Threshold_Endcap_")+ IntToString(n) + "_Zmm_Pt10To20").c_str());
 


    } else if (Label == "HWW130_Pt10To20") {
      cout << label << " " << n << endl;





    } 

    cout << "load " << n << endl;
    assert( tmpRhoMuon); 
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


    assert(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap);
    assert(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel);
    assert(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap);
    assert(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap);


    RhoDensityMuons.push_back(tmpRhoMuon->GetMean());
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

    ChargedIsoAverage_Cone03_01Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone03_01Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_01Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_01Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone03_01Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_01Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_01Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone04_01Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_01Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_01Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone04_01Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_01Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_05Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone03_05Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_05Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_05Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone03_05Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_05Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_05Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone04_05Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_05Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_05Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone04_05Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_05Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_10Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone03_10Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_10Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_10Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone03_10Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_10Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_10Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone04_10Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_10Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_10Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone04_10Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_10Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone03_15Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone03_15Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone03_15Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone03_15Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone03_15Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone03_15Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap->GetMean());
    ChargedIsoAverage_Cone04_15Threshold_MuonBarrel.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMean());
    NeutralIsoAverage_Cone04_15Threshold_MuonBarrel.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel->GetMean());
    TotalPFIsoAverage_Cone04_15Threshold_MuonBarrel.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel->GetMean());
    VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel->GetMean());
    ChargedIsoAverage_Cone04_15Threshold_MuonEndcap.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap->GetMean());
    ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMean());
    NeutralIsoAverage_Cone04_15Threshold_MuonEndcap.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap->GetMean());
    TotalPFIsoAverage_Cone04_15Threshold_MuonEndcap.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap->GetMean());
    VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap->GetMean());


 

 

    RhoDensityMuons_Error.push_back(tmpRhoMuon->GetMeanError());
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

    ChargedIsoAverage_Cone03_01Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone03_01Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_01Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_01Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_01Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone03_01Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_01Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_01Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_01Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_01Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone04_01Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_01Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_01Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_01Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone04_01Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_01Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_01Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_01Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_05Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone03_05Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_05Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_05Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_05Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone03_05Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_05Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_05Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_05Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_05Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone04_05Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_05Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_05Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_05Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone04_05Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_05Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_05Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_05Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_10Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone03_10Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_10Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_10Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_10Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone03_10Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_10Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_10Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_10Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_10Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone04_10Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_10Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_10Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_10Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone04_10Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_10Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_10Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_10Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone03_15Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone03_15Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone03_15Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone03_15Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone03_15Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone03_15Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone03_15Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone03_15Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone03_15Threshold_Endcap->GetMeanError());
    ChargedIsoAverage_Cone04_15Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Barrel->GetMeanError());
    NeutralIsoAverage_Cone04_15Threshold_MuonBarrel_Error.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Barrel->GetMeanError());
    TotalPFIsoAverage_Cone04_15Threshold_MuonBarrel_Error.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Barrel->GetMeanError());
    ChargedIsoAverage_Cone04_15Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedIso_Cone04_15Threshold_Endcap->GetMeanError());
    ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap_Error.push_back(tmpMuon_ChargedNoPUIso_Cone04_15Threshold_Endcap->GetMeanError());
    NeutralIsoAverage_Cone04_15Threshold_MuonEndcap_Error.push_back(tmpMuon_NeutralIso_Cone04_15Threshold_Endcap->GetMeanError());
    TotalPFIsoAverage_Cone04_15Threshold_MuonEndcap_Error.push_back(tmpMuon_TotalPFIso_Cone04_15Threshold_Endcap->GetMeanError());
    VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap_Error.push_back(tmpMuon_VertexSelectedPFIso_Cone04_15Threshold_Endcap->GetMeanError());
  





  }


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  ofstream outputfile("MuonEffectiveArea.txt");
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

  Double_t MuonRho = 0;

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
  MuonRho = f1->GetParameter(1);
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphRhoDensityMuons" + label + ".gif").c_str());
  outputfile << "MuonRho: " << MuonRho << endl;


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
  outputfile << "GraphCaloIsoAverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphEcalIsoAverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphHcalIsoAverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


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
  outputfile << "GraphTrkIsoAverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;






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
  outputfile << "GraphCaloIsoAverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphEcalIsoAverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphHcalIsoAverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


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
  outputfile << "GraphTrkIsoAverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



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
  outputfile << "GraphCaloIso05AverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphEcalIso05AverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphHcalIso05AverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


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
  outputfile << "GraphTrkIso05AverageMuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;






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
  outputfile << "GraphCaloIso05AverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphEcalIso05AverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

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
  outputfile << "GraphHcalIso05AverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


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
  outputfile << "GraphTrkIso05AverageMuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;




  //*******************************************************************************************
  // 01Threshold
  //*******************************************************************************************

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_01Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_01Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_01Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_01Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_01Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_01Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;





































  //*******************************************************************************************
  // 05Threshold
  //*******************************************************************************************

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_05Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_05Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_05Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_05Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_05Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_05Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;































  //*******************************************************************************************
  // 10Threshold
  //*******************************************************************************************

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_10Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_10Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_10Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_10Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_10Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_10Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;

































  //*******************************************************************************************
  // 15Threshold
  //*******************************************************************************************

  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_15Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone03_15Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone03_15Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->Fit(("ChargedIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_15Threshold_MuonBarrel[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->SetName(("GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel"+label).c_str());
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->SetTitle("");
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->Fit(("NeutralIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonBarrel[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->Fit(("TotalPFIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrelError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->Fit(("VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrelFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************



  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedIsoAverage_Cone04_15Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->SetName(("GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap"+label).c_str());
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->SetTitle("");
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->Fit(("ChargedIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap[i];
    yErrLow[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->SetName(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap"+label).c_str());
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->SetTitle("");
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->Fit(("ChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = NeutralIsoAverage_Cone04_15Threshold_MuonEndcap[i];
    yErrLow[i] = NeutralIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = NeutralIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->SetName(("GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap"+label).c_str());
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->SetTitle("");
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("NeutralIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->Fit(("NeutralIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;


  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonEndcap[i];
    yErrLow[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = TotalPFIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetName(("GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap"+label).c_str());
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetTitle("");
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("TotalPFIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->Fit(("TotalPFIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;



  //*******************************************************************************************

  for (UInt_t i=0; i<nPoints; ++i) {    
    y[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap[i];
    yErrLow[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    yErrHigh[i] = VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcapError[i];
    NPileup[i] = i;
  }

  TGraphAsymmErrors *GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap = new TGraphAsymmErrors (nPoints,  NPileup, y, NPileupError, NPileupError, yErrLow, yErrHigh);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetName(("GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap"+label).c_str());
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetTitle("");
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->SetMarkerColor(kRed);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetXaxis()->SetTitleOffset(1.02);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetYaxis()->SetTitleOffset(1.05);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->Draw("AP");
  
  f1 = new TF1(("VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(), "pol1", 0.5, 9);
  GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->Fit(("VertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcapFit"+label).c_str(),"R");
  tmpLabel = "Slope : " + DoubleToString(f1->GetParameter(1));
  LabelText = new TPaveLabel(0.7,0.3,0.9,0.5, tmpLabel.c_str(), "NDC");
  LabelText->SetBorderSize(0);
  LabelText->Draw();
  cv1->SaveAs(("GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap" + label + ".gif").c_str());
  outputfile << "GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap" + label << " : " << f1->GetParameter(1) / MuonRho << endl;









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


  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone03_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone04_01Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone03_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone04_05Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone03_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone04_10Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel,GraphChargedIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap,GraphChargedIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone03_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel,GraphChargedIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel,GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel,GraphNeutralIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel,GraphTotalPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel,GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonBarrel->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap,GraphChargedIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap,GraphChargedNoPUIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap,GraphNeutralIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap,GraphTotalPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");
  file->WriteTObject(GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap,GraphVertexSelectedPFIsoAverage_Cone04_15Threshold_MuonEndcap->GetName(), "WriteDelete");





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
