//root -l EWKAna/Hww/DYBkg/PlotDYBkgSystematics.C+

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


void PlotRInOutRatio() {

  string label = "_AfterMinMassCutJetVeto";

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);

  TStyle *MITStyle = new TStyle("MIT-Style","The Perfect Style for Plots ;-)");
  // gStyle = MITStyle;
  MITStyle->SetPalette(1);
//   TFile *f = new TFile("HwwSelectionPlots_ZllSelectionWithMassVetoJetVeto.root", "READ");
//   TFile *f = new TFile("HwwSelectionPlots_AfterJetVeto.root", "READ");
  TFile *f = new TFile("HwwDYBkgSystematicsPlots.root", "READ");
  
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;


  //**********************
  //BASELINE 2010 PU
  //**********************
  TGraphAsymmErrors *RInOutAfterMetdeltaPhilEtCut_mm = (TGraphAsymmErrors*)f->Get("RInOutAfterTCMetdeltaPhilEtCut_Zmm_2010PU_mm");
  TGraphAsymmErrors *RInOutAfterMetdeltaPhilEtCut_ee = (TGraphAsymmErrors*)f->Get("RInOutAfterTCMetdeltaPhilEtCut_Zee_2010PU_ee");
//   TGraphAsymmErrors *RInOutAfterMetdeltaPhilEtCut_mm = (TGraphAsymmErrors*)f->Get("RInOutAfterTCMetdeltaPhilEtCut_Zmm_2010PU_2020Selection_mm");
//   TGraphAsymmErrors *RInOutAfterMetdeltaPhilEtCut_ee = (TGraphAsymmErrors*)f->Get("RInOutAfterTCMetdeltaPhilEtCut_Zee_2010PU_2020Selection_ee");
  assert(RInOutAfterMetdeltaPhilEtCut_mm);
  assert(RInOutAfterMetdeltaPhilEtCut_ee);

  RInOutAfterMetdeltaPhilEtCut_mm->SetMarkerColor(colors[0]); 
  RInOutAfterMetdeltaPhilEtCut_ee->SetMarkerColor(colors[1]); 
  RInOutAfterMetdeltaPhilEtCut_mm->SetLineColor(colors[0]); 
  RInOutAfterMetdeltaPhilEtCut_ee->SetLineColor(colors[1]); 
  RInOutAfterMetdeltaPhilEtCut_mm->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMetdeltaPhilEtCut_ee->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMetdeltaPhilEtCut_mm->GetXaxis()->SetRangeUser(0, 40);
  RInOutAfterMetdeltaPhilEtCut_ee->GetXaxis()->SetRangeUser(0, 40);

  legend = new TLegend(0.43,0.70,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterMetdeltaPhilEtCut_mm, "#mu#mu final state", "LP");
  legend->AddEntry(RInOutAfterMetdeltaPhilEtCut_ee, "ee final state", "LP");
  

  RInOutAfterMetdeltaPhilEtCut_mm->SetTitle("");
  RInOutAfterMetdeltaPhilEtCut_mm->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterMetdeltaPhilEtCut_mm->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterMetdeltaPhilEtCut_mm->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterMetdeltaPhilEtCut_mm->Draw("AP");
  RInOutAfterMetdeltaPhilEtCut_ee->Draw("P");
  legend->Draw();
  cv->SaveAs(("RInOut" + label + "_MetCut_2010PU.gif").c_str());


  //**********************
  //**********************




  return;



} 
  





void PlotDYBkgSystematics() {

  PlotRInOutRatio();

}
