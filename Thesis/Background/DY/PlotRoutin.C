#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TAxis.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TGraphErrors.h>           // Graph class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O

#endif

void PlotRoutin() {

  TFile *file = new TFile("RVsMetPlots.root","READ");

  TGraph *Routin_0Jet = (TGraph*)file->Get("RVsMetGraph_ll_mH0_jetBin0");
  TGraph *Routin_1Jet = (TGraph*)file->Get("RVsMetGraph_ll_mH0_jetBin1");
  TGraph *Routin_2Jet = (TGraph*)file->Get("RVsMetGraph_ll_mH0_jetBin2");

  assert(Routin_0Jet);
  assert(Routin_1Jet);
  assert(Routin_2Jet);

  TCanvas *cv = new TCanvas("cv","cv",800,600);

  Routin_0Jet->Draw("AP");
  Routin_0Jet->SetTitle("");
  Routin_0Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_0Jet->GetYaxis()->SetTitleOffset(1.2);
  Routin_0Jet->GetYaxis()->SetRangeUser(0,0.4);
  Routin_0Jet->SetMarkerStyle(20);
  cv->SaveAs("Routin_0Jet.png");
  cv->SaveAs("Routin_0Jet.eps");


  Routin_1Jet->Draw("AP");
  Routin_1Jet->SetTitle("");
  Routin_1Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_1Jet->GetYaxis()->SetTitleOffset(1.2);
  Routin_1Jet->GetYaxis()->SetRangeUser(0,0.4);
  Routin_1Jet->SetMarkerStyle(20);
  cv->SaveAs("Routin_1Jet.png");
  cv->SaveAs("Routin_1Jet.eps");


  Routin_2Jet->Draw("AP");
  Routin_2Jet->SetTitle("");
  Routin_2Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_2Jet->GetYaxis()->SetTitleOffset(1.2);
  Routin_2Jet->GetYaxis()->SetRangeUser(0,0.4);
  Routin_2Jet->SetMarkerStyle(20);
  cv->SaveAs("Routin_2Jet.png");
  cv->SaveAs("Routin_2Jet.eps");



}
