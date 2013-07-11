//root -l EWKAna/Hww/FakeRate/plotFakeRateComparison.C+\(\)
//================================================================================================
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TCanvas.h>                // class for drawing
#include <TBenchmark.h>             // class to track macro running statistics
#include <iostream>                 // standard I/O
#include <iomanip>
#include <fstream>


// RooFit headers

#include "TFile.h"
#include "TH1D.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"

#endif

#define LUMINOSITY 2.88 //(in pb^-1)
#define NBINSPASS 60
#define NBINSFAIL 24




TH1F* NormalizeHist(TH1F *hist) {
  Double_t norm = 0;

  TH1F *histCopy = (TH1F*)hist->Clone((string(hist->GetName()) + "_copy").c_str());

  for (UInt_t b=0; b<histCopy->GetXaxis()->GetNbins()+2; ++b) {
    norm += histCopy->GetBinContent(b);
  }
  for (UInt_t b=0; b<histCopy->GetXaxis()->GetNbins()+2; ++b) {
    histCopy->SetBinContent(b,histCopy->GetBinContent(b) / norm);
    histCopy->SetBinError(b,histCopy->GetBinError(b) / norm);
  }
//  histCopy->Rebin(5);
  histCopy->Rebin(2);

  return histCopy;
}

void plotElectronFakeRateComparison( ) {

  TFile *file = 0;

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);

  TGraphAsymmErrors* FakeRate = 0;

  //Compare CutBased Vs LH Vs MVA  PT
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  FakeRate = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");
  
  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRate);
  labels.push_back("BDTG (IDIsoCombined) WP");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.2);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(10,35);
    plots[i]->GetYaxis()->SetRangeUser(0,0.10);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
 
  cv->SaveAs("ElectronFakeRate_Pt.gif");
  cv->SaveAs("ElectronFakeRate_Pt.eps");



  //Compare CutBased Vs LH Vs MVA  Eta
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  FakeRate = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Eta");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRate);
  labels.push_back("BDTG (IDIsoCombined) WP");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.2);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,2.5);
    plots[i]->GetYaxis()->SetRangeUser(0,0.10);
  }
  
  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }
  }
 
  cv->SaveAs("ElectronFakeRate_Eta.gif");
  cv->SaveAs("ElectronFakeRate_Eta.eps");





}






void plotMuonFakeRateComparison( ) {

  TFile *file = 0;

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);

  TGraphAsymmErrors* FakeRate = 0;

  //Compare CutBased Vs LH Vs MVA  PT
  file = new TFile("MuonFakeRate.BDTGIDIsoCombinedDetIso.root", "READ");
  FakeRate = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt");
  
  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRate);
  labels.push_back("BDTG (IDIsoCombined) WP");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.3);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(10,35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
 
  cv->SaveAs("MuonFakeRate_Pt.gif");
  cv->SaveAs("MuonFakeRate_Pt.eps");



  //Compare CutBased Vs LH Vs MVA  Eta
  file = new TFile("MuonFakeRate.BDTGIDIsoCombinedDetIso.root", "READ");
  FakeRate = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Eta");
 
  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRate);
  labels.push_back("BDTG (IDIsoCombined) WP");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.3);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,2.5);
  }
  
  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }
  }
 
  cv->SaveAs("MuonFakeRate_Eta.gif");
  cv->SaveAs("MuonFakeRate_Eta.eps");





}






void plotMuonFakeRateVsPileup( ) {

  TFile *file = new TFile("MuonFakeRate.BDTGIDIsoCombinedDetIso.root", "READ");

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = 0;


  //********************************************************************************
  //NVTX
  //********************************************************************************

  tmpLegend = new TLegend(0.23,0.65,0.53,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt10To20_Barrel_NVtx"));
  labels.push_back("Barrel, 10 < p_{T} < 20");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt10To20_Endcap_NVtx"));
  labels.push_back("Endcap, 10 < p_{T} < 20");
  colors.push_back(kGreen+2);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt20ToInf_Barrel_NVtx"));
  labels.push_back("Barrel, 20 < p_{T} < 35");
  colors.push_back(kBlue);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt20ToInf_Endcap_NVtx"));
  labels.push_back("Endcap, 20< p_{T} < 35");
  colors.push_back(kRed);

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.4);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,14);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("MuonFakeRate_NVtx.gif");
  cv->SaveAs("MuonFakeRate_NVtx.eps");


  //********************************************************************************
  //Rho
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt10To20_Barrel_Rho"));
  labels.push_back("Barrel, 10 < p_{T} < 20");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt10To20_Endcap_Rho"));
  labels.push_back("Endcap, 10 < p_{T} < 20");
  colors.push_back(kGreen+2);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt20ToInf_Barrel_Rho"));
  labels.push_back("Barrel, 20 < p_{T} < 35");
  colors.push_back(kBlue);
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt20ToInf_Endcap_Rho"));
  labels.push_back("Endcap, 20< p_{T} < 35");
  colors.push_back(kRed);

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.8);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,14);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("MuonFakeRate_Rho.gif");

}







void plotElectronFakeRateVsPileup( ) {

   TFile *file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");


  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);






  //********************************************************************************
  //NVtx
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_NVtx"));
  labels.push_back("Barrel, 10 < p_{T} < 20");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_NVtx"));
  labels.push_back("Endcap, 10 < p_{T} < 20");
  colors.push_back(kGreen+2);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_NVtx"));
  labels.push_back("Barrel, 20 < p_{T} < 35");
  colors.push_back(kBlue);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_NVtx"));
  labels.push_back("Endcap, 20< p_{T} < 35");
  colors.push_back(kRed);

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.1);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,10);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_NVtx.gif");
  cv->SaveAs("ElectronFakeRate_NVtx.eps");


  //********************************************************************************
  //Rho
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_Rho"));
  labels.push_back("Barrel, 10 < p_{T} < 20");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_Rho"));
  labels.push_back("Endcap, 10 < p_{T} < 20");
  colors.push_back(kGreen+2);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_Rho"));
  labels.push_back("Barrel, 20 < p_{T} < 35");
  colors.push_back(kBlue);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_Rho"));
  labels.push_back("Endcap, 20< p_{T} < 35");
  colors.push_back(kRed);

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.1);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,10);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Rho.gif");


}





void PlotLeptonJetPt(string Label = "" ) {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kBlue);
  colors.push_back(kRed);
  colors.push_back(kMagenta);
  colors.push_back(kGreen+3);
//   colors.push_back(kBlack);
//   colors.push_back(kGreen);
//   colors.push_back(kBlue);
//   colors.push_back(kRed);
//   colors.push_back(kMagenta);
//   colors.push_back(kCyan);
//   colors.push_back(kBlack);
//   colors.push_back(kGreen);
//   colors.push_back(kBlue);
//   colors.push_back(kRed);
//   colors.push_back(kMagenta);
//   colors.push_back(kCyan);
//   colors.push_back(kBlack);
//   colors.push_back(kGreen);
  
  TFile *fileFR = new TFile("MuonFakeRate.SmurfV5.root", "READ");

  vector<TH1F*> hists;
  vector<string> histLabels;
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
  TPaveLabel *LabelText;
  TLegend *legend = new TLegend(0.55,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);





 

  //****************************************************************************************************

  //****************************************************************************************************
  hists.clear();
  histLabels.clear();
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorM2_WJetsMCFiltered_ptThreshold15")); 
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorM2_Mu8Sample_ptThreshold0"));   
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorM2_Mu8Sample_ptThreshold15"));   
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorM2_Mu8Sample_ptThreshold30"));   
  
  histLabels.push_back("WJetsMC");
  histLabels.push_back("Data Jet0");
  histLabels.push_back("Data Jet15");
  histLabels.push_back("Data Jet30");  
  
  legend->Clear();
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

   for(UInt_t i=0; i<hists.size() ; ++i) {
     cout << "Load " << i << endl;
    assert(hists[i]);
    hists[i]->Rebin(2);
    hists[i] = NormalizeHist(hists[i]);
    hists[i]->SetTitle("");
    hists[i]->GetXaxis()->SetTitleOffset(1.05);
    hists[i]->GetYaxis()->SetTitleOffset(1.2);
    hists[i]->SetLineColor(colors[i]);
    hists[i]->SetMarkerColor(colors[i]); 
    hists[i]->GetYaxis()->SetRangeUser(0, 0.30); 
    hists[i]->GetXaxis()->SetTitle( "Jet p_{T} [GeV/c]"); 
    hists[i]->GetYaxis()->SetTitle( "Fraction of Events"); 
    
    if (i==0) legend->AddEntry(hists[i],histLabels[i].c_str(),"LP");
    else legend->AddEntry(hists[i],histLabels[i].c_str(),"L");
  }


  for(UInt_t i=0; i<hists.size() ; ++i) {
    if (i==0) {
      hists[i]->Draw("E1");
    } 
    else {
      hists[i]->Draw("hist,same");
    }
  }

  legend->Draw();
  cv->SaveAs("LeptonJetPt_MuonM2_0To30.gif");
  cv->SaveAs("LeptonJetPt_MuonM2_0To30.eps");

  //****************************************************************************************************

  //****************************************************************************************************
  fileFR = new TFile("ElectronFakeRate.SmurfV5.root", "READ");

  hists.clear();
  histLabels.clear();
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorV4_WJetsMCFiltered_ptThreshold15")); 
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorV4_Ele8CaloIdLCaloIsoVLSample_ptThreshold20"));   
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorV4_Ele8CaloIdLCaloIsoVLSample_ptThreshold35"));   
  hists.push_back((TH1F*)fileFR->Get("histLeptonJetPt_DenominatorV4_Ele8CaloIdLCaloIsoVLSample_ptThreshold50"));   


  histLabels.push_back("WJetsMC");
  histLabels.push_back("Data Jet20");
  histLabels.push_back("Data Jet35");
  histLabels.push_back("Data Jet50");  
  
  legend->Clear();
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);

   for(UInt_t i=0; i<hists.size() ; ++i) {
     cout << "Load " << i << endl;

    assert(hists[i]);
    hists[i]->Rebin(2);
    hists[i] = NormalizeHist(hists[i]);
    hists[i]->SetTitle("");
    hists[i]->GetXaxis()->SetTitleOffset(1.05);
    hists[i]->GetYaxis()->SetTitleOffset(1.2);
    hists[i]->SetLineColor(colors[i]);
    hists[i]->SetMarkerColor(colors[i]); 
    hists[i]->GetYaxis()->SetRangeUser(0, 0.20); 
    hists[i]->GetXaxis()->SetTitle( "Jet p_{T} [GeV/c]"); 
    hists[i]->GetYaxis()->SetTitle( "Fraction of Events"); 

    if (i==0) legend->AddEntry(hists[i],histLabels[i].c_str(),"LP");
    else legend->AddEntry(hists[i],histLabels[i].c_str(),"L");
  }


  for(UInt_t i=0; i<hists.size() ; ++i) {
    if (i==0) {
      hists[i]->Draw("E1");
    } 
    else {
      hists[i]->Draw("hist,same");
    }
  }

  legend->Draw();
  cv->SaveAs("LeptonJetPt_ElectronV4_20To50.gif");
  cv->SaveAs("LeptonJetPt_ElectronV4_20To50.eps");




//   //****************************************************************************************************

//   //****************************************************************************************************
//   hists.clear();
//   histLabels.clear();
//   hists.push_back((TH1F*)(fileFR->Get("histLeptonJetPt_TightPlusFailSample")->Clone("1"))); 
//   hists.push_back((TH1F*)(fileFR->Get("histLeptonJetPt_TightPlusFailSample")->Clone("2")));   
  
//   histLabels.push_back("Simulation");
//   histLabels.push_back("Predicted");

  
//   legend->Clear();
//   legend->SetBorderSize(0);
//   legend->SetFillStyle(0);


//    for(UInt_t i=0; i<hists.size() ; ++i) {
//      cout << "Load " << i << endl;

//     assert(hists[i]);
//     hists[i] = NormalizeHist(hists[i]);
//     hists[i]->SetTitle("");
//     hists[i]->GetXaxis()->SetTitleOffset(1.05);
//     hists[i]->GetYaxis()->SetTitleOffset(1.2);
//     hists[i]->SetLineColor(colors[i]);
//     hists[i]->SetMarkerColor(colors[i]); 
//     hists[i]->GetYaxis()->SetRangeUser(0, 0.15); 
//     hists[i]->GetXaxis()->SetTitle( "Jet p_{T} [GeV/c]"); 
    
//     if (i==0) legend->AddEntry(hists[i],histLabels[i].c_str(),"LP");
//     else legend->AddEntry(hists[i],histLabels[i].c_str(),"LP");
//   }
//    hists[0]->SetLineColor(kBlue);
//    hists[0]->SetMarkerColor(kBlue);
//    hists[1]->SetLineColor(kRed);
//    hists[1]->SetMarkerColor(kRed);
  


//   for(UInt_t i=0; i<hists.size() ; ++i) {
//     if (i==0) {
//       hists[i]->Draw("E1");
//     } 
//     else {
//       hists[i]->Draw("hist,same");
//     }
//   }

//   legend->Draw();
//   cv->SaveAs("MCClosureTest.test.gif");


}



void plotFakeRate( ) {

//   plotElectronFakeRateComparison();
//   plotElectronFakeRateVsPileup();

//   plotMuonFakeRateComparison();
//    plotMuonFakeRateVsPileup();;

   PlotLeptonJetPt();


}
