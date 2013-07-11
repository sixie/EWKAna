#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TAxis.h>                  // class to access ntuples
#include <TH1F.h>                   
#include <TLegend.h>                
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


void DoPlotROutIn() {


  TFile *file = new TFile("RVsMetPlots.root","READ");

  TGraphErrors *Routin_0Jet = (TGraphErrors*)file->Get("RVsMetGraph_ll_mH0_jetBin0");
  TGraphErrors *Routin_1Jet = (TGraphErrors*)file->Get("RVsMetGraph_ll_mH0_jetBin1");
  TGraphErrors *Routin_2Jet = (TGraphErrors*)file->Get("RVsMetGraph_ll_mH0_jetBin2");

  assert(Routin_0Jet);
  assert(Routin_1Jet);
  assert(Routin_2Jet);

  TCanvas *cv = 0;


  cv = new TCanvas("cv","cv",800,600);         
  Routin_0Jet->SetTitle("");
  Routin_0Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_0Jet->GetXaxis()->SetRangeUser(20,45);
  Routin_0Jet->GetXaxis()->SetTitleOffset(1.08);
  Routin_0Jet->GetYaxis()->SetTitleOffset(1.25);
  Routin_0Jet->GetYaxis()->SetTitle("R_{out/in}");
  Routin_0Jet->GetYaxis()->SetRangeUser(0.0,0.40);
  Routin_0Jet->Draw("AP");
  cv->SaveAs("Routin_0Jet.eps");

  cv = new TCanvas("cv","cv",800,600);         
  Routin_1Jet->SetTitle("");
  Routin_1Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_1Jet->GetXaxis()->SetRangeUser(20,45);
  Routin_1Jet->GetYaxis()->SetTitle("R_{out/in}");
  Routin_1Jet->GetYaxis()->SetRangeUser(0.0,0.40);
  Routin_1Jet->GetXaxis()->SetTitleOffset(1.08);
  Routin_1Jet->GetYaxis()->SetTitleOffset(1.25);
  Routin_1Jet->Draw("AP");
  cv->SaveAs("Routin_1Jet.eps");

  cv = new TCanvas("cv","cv",800,600);         
  Routin_2Jet->SetTitle("");
  Routin_2Jet->GetXaxis()->SetTitle("Minimum Projected Missing Transverse Energy [GeV]");
  Routin_2Jet->GetXaxis()->SetRangeUser(20,45);
  Routin_2Jet->GetYaxis()->SetTitle("R_{out/in}");
  Routin_2Jet->GetYaxis()->SetRangeUser(0.0,0.40);
  Routin_2Jet->GetXaxis()->SetTitleOffset(1.08);
  Routin_2Jet->GetYaxis()->SetTitleOffset(1.25);
  Routin_2Jet->Draw("AP");
  cv->SaveAs("Routin_2Jet.eps");



}


void PlotROutInVsMHiggs() {

  Double_t x[80];

  for( UInt_t i=0; i < 80; ++i) x[i] = i;

  TH1F *Routin_0Jet = new TH1F("", "", 79, x);
  Routin_0Jet->GetYaxis()->SetTitle("R_{out/in} Ratio");
  Routin_0Jet->GetYaxis()->SetTitleOffset(1.2);
  Routin_0Jet->SetBinContent(4,0.20); Routin_0Jet->SetBinError(4,TMath::Sqrt(pow(0.02,2) + pow(0.01,2)));
  Routin_0Jet->SetBinContent(8, 0.15); Routin_0Jet->SetBinError(8,TMath::Sqrt(pow(0.02,2) + pow(0.08,2)));
  Routin_0Jet->SetBinContent(12, 0.15); Routin_0Jet->SetBinError(12,TMath::Sqrt(pow(0.02,2) + pow(0.08,2)));
  Routin_0Jet->SetBinContent(16, 0.15); Routin_0Jet->SetBinError(16,TMath::Sqrt(pow(0.02,2) + pow(0.08,2)));
  Routin_0Jet->SetBinContent(20, 0.20); Routin_0Jet->SetBinError(20,TMath::Sqrt(pow(0.03,2) + pow(0.09,2)));
  Routin_0Jet->SetBinContent(24, 0.25); Routin_0Jet->SetBinError(24,TMath::Sqrt(pow(0.04,2) + pow(0.12,2)));
  Routin_0Jet->SetBinContent(28, 0.32); Routin_0Jet->SetBinError(28,TMath::Sqrt(pow(0.05,2) + pow(0.21,2)));
  Routin_0Jet->SetBinContent(32, 0.41); Routin_0Jet->SetBinError(32,TMath::Sqrt(pow(0.07,2) + pow(0.17,2)));
  Routin_0Jet->SetBinContent(36, 0.51); Routin_0Jet->SetBinError(36,TMath::Sqrt(pow(0.09,2) + pow(0.20,2)));
  Routin_0Jet->SetBinContent(40, 0.51); Routin_0Jet->SetBinError(40,TMath::Sqrt(pow(0.09,2) + pow(0.20,2)));
  Routin_0Jet->SetBinContent(44, 0.51); Routin_0Jet->SetBinError(44,TMath::Sqrt(pow(0.09,2) + pow(0.20,2)));
  Routin_0Jet->SetBinContent(48, 0.33); Routin_0Jet->SetBinError(48,TMath::Sqrt(pow(0.09,2) + pow(0.34,2)));
  Routin_0Jet->SetBinContent(52, 0.65); Routin_0Jet->SetBinError(52,TMath::Sqrt(pow(0.27,2) + pow(0.25,2)));
  Routin_0Jet->SetBinContent(56, 0.68); Routin_0Jet->SetBinError(56,TMath::Sqrt(pow(0.30,2) + pow(0.09,2)));
  Routin_0Jet->SetBinContent(60, 0.59); Routin_0Jet->SetBinError(60,TMath::Sqrt(pow(0.27,2) + pow(0.19,2)));
  Routin_0Jet->SetBinContent(64, 0.39); Routin_0Jet->SetBinError(64,TMath::Sqrt(pow(0.13,2) + pow(0.12,2)));
  Routin_0Jet->SetBinContent(68, 0.18); Routin_0Jet->SetBinError(68,TMath::Sqrt(pow(0.06,2) + pow(0.30,2)));
  Routin_0Jet->SetBinContent(72, 0.05); Routin_0Jet->SetBinError(72,TMath::Sqrt(pow(0.01,2) + pow(0.01,2)));
  Routin_0Jet->SetBinContent(76, 0.26); Routin_0Jet->SetBinError(76,TMath::Sqrt(pow(0.07,2) + pow(0.14,2)));
  Routin_0Jet->GetXaxis()->SetTitle("Higgs Mass Hypothesis [GeV/c^{2}]");
  Routin_0Jet->GetXaxis()->SetTitleSize(0.045);
  Routin_0Jet->GetXaxis()->SetTitleOffset(1.25);
  Routin_0Jet->GetXaxis()->SetTickLength(0);
  Routin_0Jet->GetXaxis()->SetBinLabel(4,"0");
  Routin_0Jet->GetXaxis()->SetBinLabel(8,"115");
  Routin_0Jet->GetXaxis()->SetBinLabel(12,"118");
  Routin_0Jet->GetXaxis()->SetBinLabel(16,"120");
  Routin_0Jet->GetXaxis()->SetBinLabel(20,"122");
  Routin_0Jet->GetXaxis()->SetBinLabel(24,"124");
  Routin_0Jet->GetXaxis()->SetBinLabel(28,"126");
  Routin_0Jet->GetXaxis()->SetBinLabel(32,"128");
  Routin_0Jet->GetXaxis()->SetBinLabel(36,"130");
  Routin_0Jet->GetXaxis()->SetBinLabel(40,"135");
  Routin_0Jet->GetXaxis()->SetBinLabel(44,"140");
  Routin_0Jet->GetXaxis()->SetBinLabel(48,"150");
  Routin_0Jet->GetXaxis()->SetBinLabel(52,"160");
  Routin_0Jet->GetXaxis()->SetBinLabel(56,"170");
  Routin_0Jet->GetXaxis()->SetBinLabel(60,"180");
  Routin_0Jet->GetXaxis()->SetBinLabel(64,"190");
  Routin_0Jet->GetXaxis()->SetBinLabel(68,"200");
  Routin_0Jet->GetXaxis()->SetBinLabel(72,"250");
  Routin_0Jet->GetXaxis()->SetBinLabel(76,"300");


  TH1F *Routin_1Jet = new TH1F("", "", 79, x);
  Routin_1Jet->GetYaxis()->SetTitle("R_{out/in} Ratio");
  Routin_1Jet->GetYaxis()->SetTitleOffset(1.2);
  Routin_1Jet->SetBinContent(4,0.16); Routin_1Jet->SetBinError(4,TMath::Sqrt(pow(0.01,2) + pow(0.02,2)));
  Routin_1Jet->SetBinContent(8, 0.08); Routin_1Jet->SetBinError(8,TMath::Sqrt(pow(0.01,2) + pow(0.01,2)));
  Routin_1Jet->SetBinContent(12, 0.08); Routin_1Jet->SetBinError(12,TMath::Sqrt(pow(0.01,2) + pow(0.01,2)));
  Routin_1Jet->SetBinContent(16, 0.08); Routin_1Jet->SetBinError(16,TMath::Sqrt(pow(0.01,2) + pow(0.01,2)));
  Routin_1Jet->SetBinContent(20, 0.09); Routin_1Jet->SetBinError(20,TMath::Sqrt(pow(0.01,2) + pow(0.02,2)));
  Routin_1Jet->SetBinContent(24, 0.11); Routin_1Jet->SetBinError(24,TMath::Sqrt(pow(0.01,2) + pow(0.03,2)));
  Routin_1Jet->SetBinContent(28, 0.13); Routin_1Jet->SetBinError(28,TMath::Sqrt(pow(0.01,2) + pow(0.04,2)));
  Routin_1Jet->SetBinContent(32, 0.16); Routin_1Jet->SetBinError(32,TMath::Sqrt(pow(0.02,2) + pow(0.05,2)));
  Routin_1Jet->SetBinContent(36, 0.19); Routin_1Jet->SetBinError(36,TMath::Sqrt(pow(0.02,2) + pow(0.05,2)));
  Routin_1Jet->SetBinContent(40, 0.19); Routin_1Jet->SetBinError(40,TMath::Sqrt(pow(0.02,2) + pow(0.05,2)));
  Routin_1Jet->SetBinContent(44, 0.19); Routin_1Jet->SetBinError(44,TMath::Sqrt(pow(0.02,2) + pow(0.05,2)));
  Routin_1Jet->SetBinContent(48, 0.11); Routin_1Jet->SetBinError(48,TMath::Sqrt(pow(0.01,2) + pow(0.05,2)));
  Routin_1Jet->SetBinContent(52, 0.23); Routin_1Jet->SetBinError(52,TMath::Sqrt(pow(0.04,2) + pow(0.07,2)));
  Routin_1Jet->SetBinContent(56, 0.21); Routin_1Jet->SetBinError(56,TMath::Sqrt(pow(0.04,2) + pow(0.06,2)));
  Routin_1Jet->SetBinContent(60, 0.19); Routin_1Jet->SetBinError(60,TMath::Sqrt(pow(0.03,2) + pow(0.06,2)));
  Routin_1Jet->SetBinContent(64, 0.17); Routin_1Jet->SetBinError(64,TMath::Sqrt(pow(0.02,2) + pow(0.03,2)));
  Routin_1Jet->SetBinContent(68, 0.12); Routin_1Jet->SetBinError(68,TMath::Sqrt(pow(0.01,2) + pow(0.01,2)));
  Routin_1Jet->SetBinContent(72, 0.06); Routin_1Jet->SetBinError(72,TMath::Sqrt(pow(0.01,2) + pow(0.02,2)));
  Routin_1Jet->SetBinContent(76, 0.08); Routin_1Jet->SetBinError(76,TMath::Sqrt(pow(0.01,2) + pow(0.04,2)));
  Routin_1Jet->GetXaxis()->SetTitle("Higgs Mass Hypothesis [GeV/c^{2}]");
  Routin_1Jet->GetXaxis()->SetTitleSize(0.045);
  Routin_1Jet->GetXaxis()->SetTitleOffset(1.25);
  Routin_1Jet->GetXaxis()->SetTickLength(0);
  Routin_1Jet->GetXaxis()->SetBinLabel(4,"0");
  Routin_1Jet->GetXaxis()->SetBinLabel(8,"115");
  Routin_1Jet->GetXaxis()->SetBinLabel(12,"118");
  Routin_1Jet->GetXaxis()->SetBinLabel(16,"120");
  Routin_1Jet->GetXaxis()->SetBinLabel(20,"122");
  Routin_1Jet->GetXaxis()->SetBinLabel(24,"124");
  Routin_1Jet->GetXaxis()->SetBinLabel(28,"126");
  Routin_1Jet->GetXaxis()->SetBinLabel(32,"128");
  Routin_1Jet->GetXaxis()->SetBinLabel(36,"130");
  Routin_1Jet->GetXaxis()->SetBinLabel(40,"135");
  Routin_1Jet->GetXaxis()->SetBinLabel(44,"140");
  Routin_1Jet->GetXaxis()->SetBinLabel(48,"150");
  Routin_1Jet->GetXaxis()->SetBinLabel(52,"160");
  Routin_1Jet->GetXaxis()->SetBinLabel(56,"170");
  Routin_1Jet->GetXaxis()->SetBinLabel(60,"180");
  Routin_1Jet->GetXaxis()->SetBinLabel(64,"190");
  Routin_1Jet->GetXaxis()->SetBinLabel(68,"200");
  Routin_1Jet->GetXaxis()->SetBinLabel(72,"250");
  Routin_1Jet->GetXaxis()->SetBinLabel(76,"300");


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  Routin_0Jet->Draw("E1");
  cv->SaveAs("Routin_0Jet.eps");


  Routin_1Jet->Draw("E1");
  cv->SaveAs("Routin_1Jet.eps");


}

void PlotROutIn() {

  PlotROutInVsMHiggs();

}
