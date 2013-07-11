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


void PlotElectronPtSpectra() {

  TFile *file = 0;
  TH1F *Source = 0;
  TH1F *Target = 0;
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*********************************************************
  //Background Source and Target
  //*********************************************************

  file = new TFile("ElectronPtSpectrum.root","READ");

  Source = (TH1F*)file->Get("FakeElePtSource");
  Target = (TH1F*)file->Get("FakeElePtTarget");

  assert(Source);
  assert(Target);

  NormalizeHist(Source);
  NormalizeHist(Target);


  cv = new TCanvas("cv","cv",800,600);
  Source->SetTitle("");
  Source->GetXaxis()->SetTitle("p_{T} [GeV]");
  Source->GetYaxis()->SetTitle("Fraction of Events");
  Source->GetYaxis()->SetTitleOffset(1.3);
  Source->GetXaxis()->SetTitleOffset(1.05);
  Source->Draw("hist");
  cv->SaveAs("ElectronMVATraining_Fakes_Pt_Source.eps");

  cv = new TCanvas("cv","cv",800,600);
  Target->Smooth(10);
  Target->SetTitle("");
  Target->GetXaxis()->SetTitle("p_{T} [GeV]");
  Target->GetYaxis()->SetTitle("Fraction of Events");
  Target->GetYaxis()->SetTitleOffset(1.3);
  Target->GetXaxis()->SetTitleOffset(1.05);
  Target->Draw("hist");
  cv->SaveAs("ElectronMVATraining_Fakes_Pt_Target.eps");


  //*********************************************************
  //Signal Source and Target
  //*********************************************************

  Source = (TH1F*)file->Get("RealElePtSource");
  Target = (TH1F*)file->Get("RealElePtTarget");

  assert(Source);
  assert(Target);

  NormalizeHist(Source);
  NormalizeHist(Target);


  cv = new TCanvas("cv","cv",800,600);
  Source->SetTitle("");
  Source->GetXaxis()->SetTitle("p_{T} [GeV]");
  Source->GetYaxis()->SetTitle("Fraction of Events");
  Source->GetYaxis()->SetTitleOffset(1.3);
  Source->GetXaxis()->SetTitleOffset(1.05);
  Source->Draw("hist");
  cv->SaveAs("ElectronMVATraining_Real_Pt_Source.eps");

  cv = new TCanvas("cv","cv",800,600);
  Target->SetTitle("");
  Target->GetXaxis()->SetTitle("p_{T} [GeV]");
  Target->GetYaxis()->SetTitle("Fraction of Events");
  Target->GetYaxis()->SetTitleOffset(1.3);
  Target->GetXaxis()->SetTitleOffset(1.05);
  Target->Draw("hist");
  cv->SaveAs("ElectronMVATraining_Real_Pt_Target.eps");



}

void PlotMuonPtSpectra() {

  TFile *file = 0;
  TH1F *Source = 0;
  TH1F *Target = 0;
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*********************************************************
  //Background Source and Target
  //*********************************************************

  file = new TFile("MuonPtSpectrum.root","READ");

  Source = (TH1F*)file->Get("FakeMuPtSource");
  Target = (TH1F*)file->Get("FakeMuPtTarget");

  assert(Source);
  assert(Target);

  NormalizeHist(Source);
  NormalizeHist(Target);


  cv = new TCanvas("cv","cv",800,600);
  Source->SetTitle("");
  Source->GetXaxis()->SetTitle("p_{T} [GeV]");
  Source->GetYaxis()->SetTitle("Fraction of Events");
  Source->GetYaxis()->SetTitleOffset(1.3);
  Source->GetXaxis()->SetTitleOffset(1.05);
  Source->Draw("hist");
  cv->SaveAs("MuonMVATraining_Fakes_Pt_Source.eps");

  cv = new TCanvas("cv","cv",800,600);
  Target->Smooth(10);
  Target->SetTitle("");
  Target->GetXaxis()->SetTitle("p_{T} [GeV]");
  Target->GetYaxis()->SetTitle("Fraction of Events");
  Target->GetYaxis()->SetTitleOffset(1.3);
  Target->GetXaxis()->SetTitleOffset(1.05);
  Target->Draw("hist");
  cv->SaveAs("MuonMVATraining_Fakes_Pt_Target.eps");


  //*********************************************************
  //Signal Source and Target
  //*********************************************************

  Source = (TH1F*)file->Get("RealMuPtSource");
  Target = (TH1F*)file->Get("RealMuPtTarget");

  assert(Source);
  assert(Target);

  NormalizeHist(Source);
  NormalizeHist(Target);


  cv = new TCanvas("cv","cv",800,600);
  Source->SetTitle("");
  Source->GetXaxis()->SetTitle("p_{T} [GeV]");
  Source->GetYaxis()->SetTitle("Fraction of Events");
  Source->GetYaxis()->SetTitleOffset(1.3);
  Source->GetXaxis()->SetTitleOffset(1.05);
  Source->Draw("hist");
  cv->SaveAs("MuonMVATraining_Real_Pt_Source.eps");

  cv = new TCanvas("cv","cv",800,600);
  Target->SetTitle("");
  Target->GetXaxis()->SetTitle("p_{T} [GeV]");
  Target->GetYaxis()->SetTitle("Fraction of Events");
  Target->GetYaxis()->SetTitleOffset(1.3);
  Target->GetXaxis()->SetTitleOffset(1.05);
  Target->Draw("hist");
  cv->SaveAs("MuonMVATraining_Real_Pt_Target.eps");



}

void PlotMVATrainingSpectra() {

  PlotElectronPtSpectra();
  PlotMuonPtSpectra();

}
