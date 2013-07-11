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

//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

void MakeDataMCPredictionComparison() {

  TFile *file = 0;
  TH1F *DataPrediction = 0;
  TH1F *MCPrediction = 0;
  TCanvas *cv = 0;
  TLegend *legend = 0;

  file = new TFile("FakeLeptonBkgPlots.root","READ");


  //*********************************************************
  //0 Jet Bin - MVA
  //*********************************************************
  DataPrediction = (TH1F*)file->Get("hMVA_FakeLeptonBkg_Data_FakeLepton_AllFinalStates_ZeroJetBin");
  MCPrediction = (TH1F*)file->Get("hMVA_FakeLeptonBkg_MC_FakeLepton_AllFinalStates_ZeroJetBin");

  assert(DataPrediction);
  assert(MCPrediction);

  NormalizeHist(DataPrediction);
  NormalizeHist(MCPrediction);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.5, 0.75, 0.9, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DataPrediction, "Data-Driven Prediction", "L");
  legend->AddEntry(MCPrediction, "MC Simulation Prediction", "FP");

  DataPrediction->SetLineColor(kBlue);
  MCPrediction->SetLineColor(kRed);
  DataPrediction->SetLineWidth(2);
  MCPrediction->SetLineWidth(2);
  DataPrediction->SetTitle("");
  DataPrediction->GetXaxis()->SetTitle("M_{ll} [GeV/c^{2}]");
  DataPrediction->GetYaxis()->SetTitle("Fraction of Events");
  DataPrediction->GetYaxis()->SetTitleOffset(1.4);
//  DataPrediction->GetXaxis()->SetRangeUser(-1.0,1.0);
  DataPrediction->GetYaxis()->SetRangeUser(0.0,0.30);


  DataPrediction->Draw("hist");
  MCPrediction->SetMarkerSize(0);
  MCPrediction->SetMarkerColor(kRed);
  MCPrediction->SetFillColor(kRed);
  MCPrediction->SetFillStyle(3001);
  MCPrediction->Draw("E2,same");
  DataPrediction->Draw("hist,same");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkg_MVA_DataDrivenVsMC.png");
  cv->SaveAs("FakeLeptonBkg_MVA_DataDrivenVsMC.eps");


  return;


  //*********************************************************
  //0 Jet Bin - DileptonMass
  //*********************************************************
  DataPrediction = (TH1F*)file->Get("hDileptonMass_FakeLeptonBkg_Data_FakeLepton_AllFinalStates_ZeroJetBin");
  MCPrediction = (TH1F*)file->Get("hDileptonMass_FakeLeptonBkg_MC_FakeLepton_AllFinalStates_ZeroJetBin");

  assert(DataPrediction);
  assert(MCPrediction);

  NormalizeHist(DataPrediction);
  NormalizeHist(MCPrediction);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.5, 0.75, 0.9, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DataPrediction, "Data-Driven Prediction", "L");
  legend->AddEntry(MCPrediction, "MC Simulation Prediction", "FP");

  DataPrediction->SetLineColor(kBlue);
  MCPrediction->SetLineColor(kRed);
  DataPrediction->SetLineWidth(2);
  MCPrediction->SetLineWidth(2);
  DataPrediction->SetTitle("");
  DataPrediction->GetXaxis()->SetTitle("M_{ll} [GeV/c^{2}]");
  DataPrediction->GetYaxis()->SetTitle("Fraction of Events");
  DataPrediction->GetYaxis()->SetTitleOffset(1.4);
//  DataPrediction->GetXaxis()->SetRangeUser(-1.0,1.0);
  DataPrediction->GetYaxis()->SetRangeUser(0.0,0.20);


  DataPrediction->Draw("hist");
  MCPrediction->SetMarkerSize(0);
  MCPrediction->SetMarkerColor(kRed);
  MCPrediction->SetFillColor(kRed);
  MCPrediction->SetFillStyle(3001);
  MCPrediction->Draw("E2,same");
  DataPrediction->Draw("hist,same");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkg_DileptonMass_DataDrivenVsMC.png");
  cv->SaveAs("FakeLeptonBkg_DileptonMass_DataDrivenVsMC.eps");




  //*********************************************************
  //0 Jet Bin - DeltaPhi
  //*********************************************************
  DataPrediction = (TH1F*)file->Get("hDeltaPhi_FakeLeptonBkg_Data_FakeLepton_AllFinalStates_ZeroJetBin");
  MCPrediction = (TH1F*)file->Get("hDeltaPhi_FakeLeptonBkg_MC_FakeLepton_AllFinalStates_ZeroJetBin");

  assert(DataPrediction);
  assert(MCPrediction);

  NormalizeHist(DataPrediction);
  NormalizeHist(MCPrediction);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.5, 0.75, 0.9, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DataPrediction, "Data-Driven Prediction", "L");
  legend->AddEntry(MCPrediction, "MC Simulation Prediction", "FP");

  DataPrediction->SetLineColor(kBlue);
  MCPrediction->SetLineColor(kRed);
  DataPrediction->SetLineWidth(2);
  MCPrediction->SetLineWidth(2);
  DataPrediction->SetTitle("");
  DataPrediction->GetXaxis()->SetTitle("#Delta#phi(l,l)");
  DataPrediction->GetYaxis()->SetTitle("Fraction of Events");
  DataPrediction->GetYaxis()->SetTitleOffset(1.4);
//  DataPrediction->GetXaxis()->SetRangeUser(-1.0,1.0);
  DataPrediction->GetYaxis()->SetRangeUser(0.0,0.20);


  DataPrediction->Draw("hist");
  MCPrediction->SetMarkerSize(0);
  MCPrediction->SetMarkerColor(kRed);
  MCPrediction->SetFillColor(kRed);
  MCPrediction->SetFillStyle(3001);
  MCPrediction->Draw("E2,same");
  DataPrediction->Draw("hist,same");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkg_DeltaPhi_DataDrivenVsMC.png");
  cv->SaveAs("FakeLeptonBkg_DeltaPhi_DataDrivenVsMC.eps");




 //*********************************************************
  //0 Jet Bin - MtHiggs
  //*********************************************************
  DataPrediction = (TH1F*)file->Get("hMtHiggs_FakeLeptonBkg_Data_FakeLepton_AllFinalStates_ZeroJetBin");
  MCPrediction = (TH1F*)file->Get("hMtHiggs_FakeLeptonBkg_MC_FakeLepton_AllFinalStates_ZeroJetBin");

  assert(DataPrediction);
  assert(MCPrediction);

  NormalizeHist(DataPrediction);
  NormalizeHist(MCPrediction);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.5, 0.75, 0.9, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DataPrediction, "Data-Driven Prediction", "L");
  legend->AddEntry(MCPrediction, "MC Simulation Prediction", "FP");

  DataPrediction->SetLineColor(kBlue);
  MCPrediction->SetLineColor(kRed);
  DataPrediction->SetLineWidth(2);
  MCPrediction->SetLineWidth(2);
  DataPrediction->SetTitle("");
  DataPrediction->GetXaxis()->SetTitle("Higgs Transverse Mass [GeV/c^{2}]");
  DataPrediction->GetYaxis()->SetTitle("Fraction of Events");
  DataPrediction->GetYaxis()->SetTitleOffset(1.4);
//  DataPrediction->GetXaxis()->SetRangeUser(-1.0,1.0);
  DataPrediction->GetYaxis()->SetRangeUser(0.0,0.30);


  DataPrediction->Draw("hist");
  MCPrediction->SetMarkerSize(0);
  MCPrediction->SetMarkerColor(kRed);
  MCPrediction->SetFillColor(kRed);
  MCPrediction->SetFillStyle(3001);
  MCPrediction->Draw("E2,same");
  DataPrediction->Draw("hist,same");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkg_MtHiggs_DataDrivenVsMC.png");
  cv->SaveAs("FakeLeptonBkg_MtHiggs_DataDrivenVsMC.eps");







}
