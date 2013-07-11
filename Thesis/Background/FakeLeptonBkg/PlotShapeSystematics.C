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

void PlotShapeSystematics() {

  TFile *file = 0;
  TH1F *DefaultShape = 0;
  TH1F *UpShape = 0;
  TH1F *DownShape = 0;
  TCanvas *cv = 0;
  TLegend *legend = 0;

  //*********************************************************
  //0 Jet Bin - Jet Pt Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_of_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.25);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_OF0Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_OF0Jet.eps");


  //*********************************************************
  //0 Jet Bin - Jet Pt Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.25);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_SF0Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_SF0Jet.eps");





  //*********************************************************
  //1 Jet Bin - Jet Pt Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_of_1j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.25);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_OF1Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_OF1Jet.eps");





  //*********************************************************
  //1 Jet Bin - Jet Pt Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_sf_1j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.25);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_SF1Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_JetSpectrum_SF1Jet.eps");






  //*********************************************************
  //0 Jet Bin - MC Closure Test Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_of_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.4);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_OF0Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_OF0Jet.eps");



  //*********************************************************
  //0 Jet Bin - MC Closure Test Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.4);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_SF0Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_SF0Jet.eps");



  //*********************************************************
  //1 Jet Bin - MC Closure Test Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_of_1j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.4);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_OF1Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_OF1Jet.eps");














  //*********************************************************
  //1 Jet Bin - MC Closure Test Systematic
  //*********************************************************

  file = new TFile("/build/sixie/Thesis/cards/MVAIDIsoCombined2011/130/hww130_sf_1j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_Wjets");
  UpShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwUp");
  DownShape = (TH1F*)file->Get("histo_Wjets_CMS_MVAWMCBounding_hwwDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  NormalizeHist(DefaultShape);
  NormalizeHist(UpShape);
  NormalizeHist(DownShape);


  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Bounding Shape (Up)", "L");
  legend->AddEntry(DownShape, "Bounding Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Fraction of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.4);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,0.4);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_SF1Jet.png");
  cv->SaveAs("FakeLeptonBkgShapeVariation_SampleDependence_SF1Jet.eps");




}
