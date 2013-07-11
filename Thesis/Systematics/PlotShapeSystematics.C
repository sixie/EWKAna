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
#include <MitStyle.C>             
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

TH1F* MakeRelative(TH1F *hist, TH1F *reference) {

  TH1F *relative = (TH1F*)hist->Clone((string(hist->GetName())+"_relative").c_str());
  relative->SetTitle("");
  assert(hist->GetXaxis()->GetNbins() == reference->GetXaxis()->GetNbins());

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {        
    if (reference->GetBinContent(b) > 0) {
      relative->SetBinContent(b,100 * ( hist->GetBinContent(b) - reference->GetBinContent(b) )/reference->GetBinContent(b));
//       relative->SetBinContent(b,( hist->GetBinContent(b) - reference->GetBinContent(b) )/reference->GetBinError(b));
    } else {
      relative->SetBinContent(b,0);
    }    
  }

  return relative;
}

void PlotLeptonEfficiencyShapeSystematics() {


  TFile *file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_0j.input.root","READ");
  TH1F *DefaultShape;
  TH1F *UpShape;
  TH1F *DownShape;
  TCanvas *cv ;
  TLegend *legend;



  //*************************
  //ggH signal - DF 0Jet
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-5,5);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");

  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_DF.eps");

  //*************************
  //ggH signal - SF 0Jet
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,5.5);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-5,5);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");

  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_SF.eps");

  return;

  //*************************
  //ggH signal - DF
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.2);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_DF.eps");



  //*************************
  //ggH signal - SF
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepEffBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.2);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,6);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("LeptonEfficiencyShapeVariation_ggH_0Jet_SF.eps");



}



void PlotLeptonScaleAndResolutionShapeSystematics() {


  TFile *file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_0j.input.root","READ");
  TH1F *DefaultShape;
  TH1F *UpShape;
  TH1F *DownShape;
  TCanvas *cv ;
  TLegend *legend;

  //*************************
  //ggH signal - DF 0Jet
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-120,120);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_DF.eps");

  //*************************
  //ggH signal - SF
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,5.5);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-120,120);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");

  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_SF.eps");


  return;

  //*************************
  //ggH signal - DF
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.2);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_DF.eps");


  //*************************
  //ggH signal - SF
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVALepResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);

  cv = new TCanvas("cv","cv",800,600);

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.2);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,6);


  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("LeptonScaleAndResolutionShapeVariation_ggH_0Jet_SF.eps");



}


void PlotMETResolutionShapeSystematics() {


  TFile *file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_0j.input.root","READ");
  TH1F *DefaultShape;
  TH1F *UpShape;
  TH1F *DownShape;
  TCanvas *cv ;
  TLegend *legend;



  //*************************
  //ggH signal - DF 0Jet
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAMETResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAMETResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-50,50);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("METResolutionShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("METResolutionShapeVariation_ggH_0Jet_DF.eps");

  //*************************
  //ggH signal - SF
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAMETResBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAMETResBoundingDown");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAMETResBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,5.5);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-100,100);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("METResolutionShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("METResolutionShapeVariation_ggH_0Jet_SF.eps");

}

void PlotJESShapeSystematics() {


  TFile *file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_0j.input.root","READ");
  TH1F *DefaultShape;
  TH1F *UpShape;
  TH1F *DownShape;
  TCanvas *cv ;
  TLegend *legend;



  //*************************
  //ggH signal - DF 0Jet
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAJESBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAJESBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-20,20);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("JESShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("JESShapeVariation_ggH_0Jet_DF.eps");

  //*************************
  //ggH signal - DF 1Jet
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_1j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAJESBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAJESBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,7.75);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetNdivisions(306);
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-70,70);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("JESShapeVariation_ggH_1Jet_DF.png");
  cv->SaveAs("JESShapeVariation_ggH_1Jet_DF.eps");

  return;


}


void PlotGluonFusionHiggsShapeSystematics() {


  TFile *file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_of_0j.input.root","READ");
  TH1F *DefaultShape;
  TH1F *UpShape;
  TH1F *DownShape;
  TCanvas *cv ;
  TLegend *legend;



  //*************************
  //ggH signal - DF 0Jet
  //*************************

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAggHBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAggHBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,13);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-3,3);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("ggHShapeVariation_ggH_0Jet_DF.png");
  cv->SaveAs("ggHShapeVariation_ggH_0Jet_DF.eps");

  //*************************
  //ggH signal - SF 0Jet
  //*************************
  file = new TFile("/data/blue/sixie/Thesis/Limits/MVAIDIsoCombinedDetIsoSameSigWP/130/hww130_sf_0j.input.root","READ");

  DefaultShape = (TH1F*)file->Get("histo_ggH");
  UpShape = (TH1F*)file->Get("histo_ggH_CMS_MVAggHBoundingUp");
  DownShape = (TH1F*)file->Get("histo_ggH_CMS_MVAggHBoundingDown");

  assert(DefaultShape);
  assert(UpShape);
  assert(DownShape);


  //Make relative histograms
  UpShapeRelative = MakeRelative(UpShape,DefaultShape);
  DownShapeRelative = MakeRelative(DownShape,DefaultShape);


  cv = new TCanvas("cv","cv",800,600);

  pad1 = new TPad("pad1","pad1", 0,0.2,1,1);
  pad1->SetBottomMargin(0.125);
  pad1->Draw();
  pad1->cd();

  legend = new TLegend(0.2, 0.6, 0.5, 0.9);
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillColor(kWhite);
  legend->AddEntry(DefaultShape, "Default Shape", "L");
  legend->AddEntry(UpShape, "Systematics Shape (Up)", "L");
  legend->AddEntry(DownShape, "Systematics Shape (Down)", "L");

  DefaultShape->SetLineColor(kBlack);
  UpShape->SetLineColor(kBlue);
  DownShape->SetLineColor(kRed);
  DefaultShape->SetLineWidth(2);
  UpShape->SetLineWidth(2);
  DownShape->SetLineWidth(2);
  DefaultShape->SetTitle("");
  DefaultShape->GetXaxis()->SetTitle("MVA discriminator");
  DefaultShape->GetYaxis()->SetTitle("Number of Events");
  DefaultShape->GetYaxis()->SetTitleOffset(1.1);
  DefaultShape->GetXaxis()->SetTitleOffset(1.05);
  DefaultShape->GetXaxis()->SetRangeUser(-1.0,1.0);
  DefaultShape->GetYaxis()->SetRangeUser(0.0,5.5);

  DefaultShape->Draw("hist");
  UpShape->Draw("same,hist");
  DownShape->Draw("same,hist");
  legend->Draw();

  cv->cd();
  pad2 = new TPad("pad2","pad2", 0,0,1,0.2);
  pad1->SetTopMargin(0.01);
  pad2->Draw();
  pad2->cd();

  UpShapeRelative->GetYaxis()->SetTitle("% Difference");
  UpShapeRelative->GetYaxis()->SetTitleSize(0.15);
  UpShapeRelative->GetYaxis()->SetTitleOffset(0.3);
  UpShapeRelative->GetYaxis()->SetRangeUser(-3,3);
  UpShapeRelative->GetYaxis()->SetLabelSize(0.15);
  UpShapeRelative->GetXaxis()->SetLabelSize(0.0);
  UpShapeRelative->SetLineColor(kBlue);
  UpShapeRelative->SetMarkerColor(kBlue);
  UpShapeRelative->Draw("hist");
  DownShapeRelative->SetLineColor(kRed);
  DownShapeRelative->SetMarkerColor(kRed);
  DownShapeRelative->Draw("hist,same");


  cv->SaveAs("ggHShapeVariation_ggH_0Jet_SF.png");
  cv->SaveAs("ggHShapeVariation_ggH_0Jet_SF.eps");



}


void PlotShapeSystematics() {

   PlotLeptonEfficiencyShapeSystematics();
//   PlotLeptonScaleAndResolutionShapeSystematics();
//   PlotMETResolutionShapeSystematics();
//    PlotJESShapeSystematics();
//   PlotGluonFusionHiggsShapeSystematics();

}
