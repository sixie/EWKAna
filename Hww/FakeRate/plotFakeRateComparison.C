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
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"

#endif

#define LUMINOSITY 2.88 //(in pb^-1)
#define NBINSPASS 60
#define NBINSFAIL 24


void compareStdIsoVsPFIso() {

  TFile *file = new TFile("ElectronFakeRate_V4.root", "READ");

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);


  tmpLegend = new TLegend(0.33,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTNoTrackerToWWTight_Ele10Jet30_Pt"));
  labels.push_back("Numerator: WW Cuts");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTNoTrackerToVBTF90Loose_Ele10Jet30_Pt"));
  labels.push_back("Numerator: VBTF90, RelativeCombIso < 0.3");
  colors.push_back(kBlue);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.9);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_HLTNoTracker_VsPt.gif");



}



void plotElectronFakeRateComparison( ) {

  TFile *file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);


  //Compare CutBased Vs LH Vs MVA  PT
  file = new TFile("ElectronFakeRate.SmurfV6.skim.root", "READ");
  TGraphAsymmErrors* FakeRateSmurfV6 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");
  file = new TFile("/data/smurf/sixie/FakeRates/FakeRates_BDTGWithIPInfoElectron.root", "READ");
  TGraphAsymmErrors* FakeRateSmurfMVANoIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");
  file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
  TGraphAsymmErrors* FakeRateSmurfMVAWithIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  TGraphAsymmErrors* FakeRateSmurfMVAIDIsoCombined = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
//   plots.push_back(FakeRateSmurfV6);
//   labels.push_back("Cut-Based");
//   colors.push_back(kGreen+2);
//   plots.push_back(FakeRateSmurfMVANoIPInfo);
//   labels.push_back("BDTG (no IP info) WP");
//   colors.push_back(kBlack);
  plots.push_back(FakeRateSmurfMVAWithIPInfo);
  labels.push_back("BDTG (2011) WP");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAIDIsoCombined);
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
    plots[i]->GetXaxis()->SetRangeUser(0,35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_MVANewVsOld_Pt.gif");



  //Compare CutBased Vs LH Vs MVA  Eta
  file = new TFile("ElectronFakeRate.SmurfV6.skim.root", "READ");
  FakeRateSmurfV6 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Eta");
  file = new TFile("/data/smurf/sixie/FakeRates/FakeRates_BDTGWithIPInfoElectron.root", "READ");
  FakeRateSmurfMVANoIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Eta");
  file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
  FakeRateSmurfMVAWithIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Eta");
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  FakeRateSmurfMVAIDIsoCombined = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Eta");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
//   plots.push_back(FakeRateSmurfV6);
//   labels.push_back("Cut-Based");
//   colors.push_back(kGreen+2);
//   plots.push_back(FakeRateSmurfMVANoIPInfo);
//   labels.push_back("BDTG (no IP info) WP");
//   colors.push_back(kBlack);
  plots.push_back(FakeRateSmurfMVAWithIPInfo);
  labels.push_back("BDTG (2011) WP");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAIDIsoCombined);
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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_MVANewVsOld_Eta.gif");



  //Compare CutBased Vs LH Vs MVA  NVtx
  file = new TFile("ElectronFakeRate.SmurfV6.skim.root", "READ");
  FakeRateSmurfV6 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_NVtx");
  file = new TFile("/data/smurf/sixie/FakeRates/FakeRates_BDTGWithIPInfoElectron.root", "READ");
  FakeRateSmurfMVANoIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_NVtx");
  file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
  FakeRateSmurfMVAWithIPInfo = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_NVtx");
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  FakeRateSmurfMVAIDIsoCombined = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_NVtx");

  assert(FakeRateSmurfMVAWithIPInfo);
  assert(FakeRateSmurfMVAIDIsoCombined);

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
//   plots.push_back(FakeRateSmurfV6);
//   labels.push_back("Cut-Based");
//   colors.push_back(kGreen+2);
//   plots.push_back(FakeRateSmurfMVANoIPInfo);
//   labels.push_back("BDTG (no IP info) WP");
//   colors.push_back(kBlack);
  plots.push_back(FakeRateSmurfMVAWithIPInfo);
  labels.push_back("BDTG (2011) WP");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAIDIsoCombined);
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
    plots[i]->GetXaxis()->SetRangeUser(0,26);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_MVANewVsOld_NVtx.gif");



  return;




//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_v2_Ele10Jet30_Pt"));
//   labels.push_back("Ele10");
//   colors.push_back(kBlack);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_jetcombinedEleTriggerMatched_v2_Jet30_Pt"));
//   labels.push_back("Jet30");
//   colors.push_back(kBlue);
// //   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_r10aeg_v2_Photon30_Pt"));
// //   labels.push_back("Photon30");
// //   colors.push_back(kRed);

//   assert(plots.size() == labels.size());



//   tmpLegend->Clear();
//   for(int i=0; i<labels.size(); ++i) {
//     tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
//     plots[i]->SetMarkerColor(colors[i]);
//     plots[i]->SetLineColor(colors[i]);
//     plots[i]->SetTitle("");
//     plots[i]->GetYaxis()->SetRangeUser(0,0.15);
//     plots[i]->GetYaxis()->SetTitleOffset(1.1);
//    }

//   for(int i=0; i<labels.size(); ++i) {
//     if (i==0) {
//       plots[i]->Draw("AP");
//     } else {
//       plots[i]->Draw("Psame");
//     }    
//   }
//   tmpLegend->Draw();

//   cv->SaveAs("ElectronFakeRateCompareTriggers_VsPt.gif");







  tmpLegend = new TLegend(0.33,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTNoTrackerToWWTight_Ele10Jet30_Pt"));
  labels.push_back("Numerator: WW Cuts");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTNoTrackerToVBTF90Loose_Ele10Jet30_Pt"));
  labels.push_back("Numerator: VBTF90, RelativeCombIso < 0.3");
  colors.push_back(kBlue);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.9);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_HLTNoTracker_VsPt.gif");



  tmpLegend = new TLegend(0.33,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTWithTrackerToWWTight_Ele10Jet30_Pt"));
  labels.push_back("Numerator: WW Cuts");
  colors.push_back(kBlack);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTWithTrackerToVBTF90Loose_Ele10Jet30_Pt"));
  labels.push_back("Numerator: VBTF90, RelativeCombIso < 0.3");
  colors.push_back(kBlue);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.9);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_HLTWithTracker_VsPt.gif");




//   plots.clear();
//   labels.clear();
//   colors.clear();
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EGPt10To20_v2_Ele10Jet30_Eta"));
//   labels.push_back("Ele10");
//   colors.push_back(kBlack);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_jetcombinedPt10To20_v2_Jet30_Eta"));
//   labels.push_back("Jet30");
//   colors.push_back(kBlue);
// //   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_r10aegPt10To20_v2_Photon30_Eta"));
// //   labels.push_back("Photon30");
// //   colors.push_back(kRed);
//   assert(plots.size() == labels.size());

//   tmpLegend->Clear();
//   for(int i=0; i<labels.size(); ++i) {
//     tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
//     plots[i]->SetMarkerColor(colors[i]);
//     plots[i]->SetLineColor(colors[i]);
//     plots[i]->SetTitle("");
//     plots[i]->GetYaxis()->SetRangeUser(0,0.25);
//     plots[i]->GetYaxis()->SetTitleOffset(1.1);
//    }

//   for(int i=0; i<labels.size(); ++i) {
//     if (i==0) {
//       plots[i]->Draw("AP");
//     } else {
//       plots[i]->Draw("Psame");
//     }    
//   }
//   tmpLegend->Draw();

//   cv->SaveAs("ElectronFakeRateCompareTriggers_Pt10To20_VsEta.gif");




//   plots.clear();
//   labels.clear();
//   colors.clear();
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EGPt10To20_v2_Ele10Jet30_Phi"));
//   labels.push_back("Ele10");
//   colors.push_back(kBlack);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_jetcombinedPt10To20_v2_Jet30_Phi"));
//   labels.push_back("Jet30");
//   colors.push_back(kBlue);
// //   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_r10aegPt10To20_v2_Photon30_Phi"));
// //   labels.push_back("Photon30");
// //   colors.push_back(kRed);
//   assert(plots.size() == labels.size());

//   tmpLegend->Clear();
//   for(int i=0; i<labels.size(); ++i) {
//     tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
//     plots[i]->SetMarkerColor(colors[i]);
//     plots[i]->SetLineColor(colors[i]);
//     plots[i]->SetTitle("");
//     plots[i]->GetYaxis()->SetRangeUser(0,0.25);
//     plots[i]->GetYaxis()->SetTitleOffset(1.1);
//    }

//   for(int i=0; i<labels.size(); ++i) {
//     if (i==0) {
//       plots[i]->Draw("AP");
//     } else {
//       plots[i]->Draw("Psame");
//     }    
//   }
//   tmpLegend->Draw();

//   cv->SaveAs("ElectronFakeRateCompareTriggers_Pt10To20_VsPhi.gif");


//   plots.clear();
//   labels.clear();
//   colors.clear();
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3");
//   colors.push_back(kBlack);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombIso03_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso03");
//   colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombIso02_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso02");
//   colors.push_back(kRed);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombIso01_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso01");
//   colors.push_back(kMagenta);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombinedIso01SigmaIEta_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso01SigmaIEta");
//   colors.push_back(kGreen);
//   assert(plots.size() == labels.size());

//   tmpLegend = new TLegend(0.33,0.65,0.93,0.90);   
//   tmpLegend->SetTextSize(0.04);
//   tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
//   for(int i=0; i<labels.size(); ++i) {
//     if (!plots[i]) { cout << "not found " << i << endl; continue;}
//     tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
//     plots[i]->SetMarkerColor(colors[i]);
//     plots[i]->SetLineColor(colors[i]);
//     plots[i]->SetTitle("");
//     plots[i]->GetYaxis()->SetRangeUser(0,0.5);
//     plots[i]->GetYaxis()->SetTitleOffset(1.2);
//     plots[i]->GetXaxis()->SetTitleOffset(1.1);
//    }

//   for(int i=0; i<labels.size(); ++i) {
//     if (i==0) {
//       plots[i]->Draw("AP");
//     } else {
//       plots[i]->Draw("Psame");
//     }    
//   }
//   tmpLegend->Draw();

//   cv->SaveAs("ElectronFakeRate_HLTOption3_TightenIsolation.gif");





//   plots.clear();
//   labels.clear();
//   colors.clear();
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombIso03_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso03");
//   colors.push_back(kBlack);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombinedIso03SigmaIEta_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso03+SigmaIEta");
//   colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRate_EG_HLTOption3PlusCombinedIso03VBTF90_Ele10Jet30_Pt"));
//   labels.push_back("HLTOption3+CombIso03+VBTF90");
//   colors.push_back(kRed);

//   assert(plots.size() == labels.size());

//   tmpLegend->Clear();
//   for(int i=0; i<labels.size(); ++i) {
//     if (!plots[i]) { cout << "not found " << i << endl; continue;}
//     tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
//     plots[i]->SetMarkerColor(colors[i]);
//     plots[i]->SetLineColor(colors[i]);
//     plots[i]->SetTitle("");
//     plots[i]->GetYaxis()->SetRangeUser(0,0.5);
//     plots[i]->GetYaxis()->SetTitleOffset(1.2);
//     plots[i]->GetXaxis()->SetTitleOffset(1.1);
//    }

//   for(int i=0; i<labels.size(); ++i) {
//     if (i==0) {
//       plots[i]->Draw("AP");
//     } else {
//       plots[i]->Draw("Psame");
//     }    
//   }
//   tmpLegend->Draw();

//   cv->SaveAs("ElectronFakeRate_HLTOption3CombIso03_TightenId.gif");



  //EffFile->Close();

}



void plotElectronFakeRateVsJetThresholds( ) {

  TFile *file = 0;
  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);

  TGraphAsymmErrors* FakeRateJet35 = 0;
  TGraphAsymmErrors* FakeRateJet30 = 0;
  TGraphAsymmErrors* FakeRateJet40 = 0;
  TGraphAsymmErrors* FakeRateJet20 = 0;
  TGraphAsymmErrors* FakeRateJet50 = 0;



  //Compare CutBased Vs LH Vs MVA  PT
  file = new TFile("ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "READ");
  FakeRateJet35 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt");
  FakeRateJet30 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold30_Pt");
  FakeRateJet40 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold40_Pt");
  FakeRateJet20 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold20_Pt");
  FakeRateJet50 = (TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold50_Pt");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
//   plots.push_back(FakeRateSmurfV6);
//   labels.push_back("Cut-Based");
//   colors.push_back(kGreen+2);
//   plots.push_back(FakeRateSmurfMVANoIPInfo);
//   labels.push_back("BDTG (no IP info) WP");
//   colors.push_back(kBlack);
  plots.push_back(FakeRateJet35);
  labels.push_back("Nominal (Jet35)");
  colors.push_back(kRed);
  plots.push_back(FakeRateJet30);
  labels.push_back("Jet30");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateJet40);
  labels.push_back("Jet40");
  colors.push_back(kGreen-2);
  plots.push_back(FakeRateJet20);
  labels.push_back("Jet20");
  colors.push_back(kBlue+2);
  plots.push_back(FakeRateJet50);
  labels.push_back("Jet50");
  colors.push_back(kBlue-2);

  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");

    //set nominal    
    if (i==0) {
      plots[i]->SetMarkerStyle(29);
      plots[i]->SetMarkerSize(2.0);
    }
    //erase the error bars
    for(UInt_t p=0; p<plots[i]->GetN(); ++p) {
      plots[i]->SetPointError(p,0,0,0,0);
    }

    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.08);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("ElectronFakeRate_JetSpectrumSystematics_Pt.gif");





  return;


}



void plotMuonFakeRateVsJetThresholds( ) {

  TFile *file = 0;
  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);

  TGraphAsymmErrors* FakeRateJet15 = 0;
  TGraphAsymmErrors* FakeRateJet20 = 0;
  TGraphAsymmErrors* FakeRateJet10 = 0;
  TGraphAsymmErrors* FakeRateJet0 = 0;
  TGraphAsymmErrors* FakeRateJet30 = 0;



  //Compare CutBased Vs LH Vs MVA  PT
   file = new TFile("MuonFakeRate.BDTGIDIsoCombined.root", "READ");
//  file = new TFile("MuonFakeRate.SmurfV7.root", "READ");

  FakeRateJet15 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt");
  FakeRateJet20 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold20_Pt");
  FakeRateJet10 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold10_Pt");
  FakeRateJet0 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold0_Pt");
  FakeRateJet30 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold30_Pt");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
//   plots.push_back(FakeRateSmurfV6);
//   labels.push_back("Cut-Based");
//   colors.push_back(kGreen+2);
//   plots.push_back(FakeRateSmurfMVANoIPInfo);
//   labels.push_back("BDTG (no IP info) WP");
//   colors.push_back(kBlack);
  plots.push_back(FakeRateJet15);
  labels.push_back("Nominal (Jet15)");
  colors.push_back(kRed);
  plots.push_back(FakeRateJet10);
  labels.push_back("Jet10");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateJet20);
  labels.push_back("Jet20");
  colors.push_back(kGreen-2);
  plots.push_back(FakeRateJet0);
  labels.push_back("Jet0");
  colors.push_back(kBlue+2);
  plots.push_back(FakeRateJet30);
  labels.push_back("Jet30");
  colors.push_back(kBlue-2);

  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");

    //set nominal    
    if (i==0) {
      plots[i]->SetMarkerStyle(29);
      plots[i]->SetMarkerSize(2.0);
    }
    //erase the error bars
    for(UInt_t p=0; p<plots[i]->GetN(); ++p) {
      plots[i]->SetPointError(p,0,0,0,0);
    }

    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.18);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
    plots[i]->GetXaxis()->SetRangeUser(0,35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("MuonFakeRate_JetSpectrumSystematics_Pt.gif");





  return;


}






void plotMuonFakeRateComparison( ) {

  TFile *file = new TFile("MuonFakeRate.SmurfV7.root", "READ");

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);



  //
  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);



  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt"));
  labels.push_back("BDTG (with IP info) WP");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

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
    plots[i]->GetXaxis()->SetRangeUser(0,35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
//   tmpLegend->Draw();

  cv->SaveAs("MuonFakeRate_Pt.gif");






  TGraphAsymmErrors* FakeRateSmurfV7 = 0;
  TGraphAsymmErrors* FakeRateSmurfMVAV8 = 0;

  //Compare CutBased Vs MVA V3 PT
  file = new TFile("MuonFakeRate.SmurfV7.root", "READ");
  FakeRateSmurfV7 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt");
  file = new TFile("MuonFakeRate.BDTGIDIsoCombined.root", "READ");
  FakeRateSmurfMVAV8 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Pt");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRateSmurfV7);
  labels.push_back("Cut-Based");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAV8);
  labels.push_back("BDTG IDIsoCombined");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

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
    plots[i]->GetXaxis()->SetRangeUser(0, 35);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("MuonFakeRate_CutBasedVsMVA_Pt.gif");


  //Compare CutBased Vs MVA V3 Eta
  file = new TFile("MuonFakeRate.SmurfV7.root", "READ");
  FakeRateSmurfV7 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Eta");
  file = new TFile("MuonFakeRate.BDTGV8.root", "READ");
  FakeRateSmurfMVAV8 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Eta");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRateSmurfV7);
  labels.push_back("Cut-Based");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAV8);
  labels.push_back("BDTG V8");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.5);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("MuonFakeRate_CutBasedVsMVA_Eta.gif");



  //Compare CutBased Vs MVA V3 NVtx
  file = new TFile("MuonFakeRate.SmurfV7.root", "READ");
  FakeRateSmurfV7 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_NVtx");
  file = new TFile("MuonFakeRate.BDTGV8.root", "READ");
  FakeRateSmurfMVAV8 = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_NVtx");

  tmpLegend = new TLegend(0.63,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetFillStyle(0);
  tmpLegend->SetBorderSize(0);
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back(FakeRateSmurfV7);
  labels.push_back("Cut-Based");
  colors.push_back(kGreen+2);
  plots.push_back(FakeRateSmurfMVAV8);
  labels.push_back("BDTG V8");
  colors.push_back(kRed);
  assert(plots.size() == labels.size());

  tmpLegend->Clear();
  for(int i=0; i<labels.size(); ++i) {
    if (!plots[i]) {cout << "not found " << i << endl; continue;}
    tmpLegend->AddEntry(plots[i], labels[i].c_str(), "LP");
    plots[i]->SetMarkerColor(colors[i]);
    plots[i]->SetLineColor(colors[i]);
    plots[i]->SetTitle("");
    plots[i]->GetYaxis()->SetRangeUser(0,0.5);
    plots[i]->GetYaxis()->SetTitleOffset(1.1);
    plots[i]->GetXaxis()->SetTitleOffset(1.05);
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();

  cv->SaveAs("MuonFakeRate_CutBasedVsMVA_NVtx.gif");



}


void plotMuonFakeRateVsPileup( ) {

  TFile *file = new TFile("MuonFakeRate.SmurfV6.skim.root", "READ");

  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);


  //********************************************************************************
  //NVTX
  //********************************************************************************
  file = new TFile("MuonFakeRate.SmurfV6.skim.root", "READ");
  TGraphAsymmErrors* FakeRate_NVtx = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_NVtx");
  file = new TFile("MuonFakeRate.SmurfV6.skim.root", "READ");
  TGraphAsymmErrors* FakeRate_Rho = (TGraphAsymmErrors*)file->Get("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_Rho");



  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
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
    plots[i]->GetYaxis()->SetRangeUser(0,1.0);
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


  //********************************************************************************
  //Rho
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
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

//   TFile *file = new TFile("ElectronFakeRate.SmurfV6.skim.root", "READ");
//   TFile *file = new TFile("ElectronFakeRate.SmurfLH.skim.root", "READ");
//    TFile *file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
   TFile *fileA = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
   TFile *fileB = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");
   TFile *file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");


  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.6,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);




  //********************************************************************************
  //NVTX
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_NVtx"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_NVtx"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_NVtx"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_NVtx_Pt10to20_Barrel.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_NVtx"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_NVtx"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_NVtx"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_NVtx_Pt10to20_Endcap.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_NVtx"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_NVtx"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_NVtx"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_NVtx_Pt20toInf_Barrel.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_NVtx"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_NVtx"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_NVtx"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_NVtx_Pt20toInf_Endcap.Run2011AVsB.gif");




  //********************************************************************************
  //Rho
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_Rho"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_Rho"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_Rho"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Rho_Pt10to20_Barrel.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_Rho"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_Rho"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Endcap_Rho"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Rho_Pt10to20_Endcap.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_Rho"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_Rho"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Barrel_Rho"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Rho_Pt20toInf_Barrel.Run2011AVsB.gif");

  //********************************************************************************


  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_Rho"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_Rho"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt20ToInf_Endcap_Rho"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Rho_Pt20toInf_Endcap.Run2011AVsB.gif");




  //********************************************************************************
  //NVtx
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
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
    plots[i]->GetYaxis()->SetRangeUser(0,0.3);
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


  //********************************************************************************
  //Rho
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
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
    plots[i]->GetYaxis()->SetRangeUser(0,0.3);
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













  //********************************************************************************
  //PT
  //********************************************************************************

  tmpLegend = new TLegend(0.63,0.65,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)fileA->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt"));
  labels.push_back("Run2011A");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)fileB->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt"));
  labels.push_back("Run2011B");
  colors.push_back(kBlue);
//   plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_Pt10To20_Barrel_Pt"));
//   labels.push_back("Full2011");
//   colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Pt.Run2011AVsB.gif");

  //********************************************************************************









}




void plotElectronFakeRateVsTriggers( ) {

//    TFile *file = new TFile("ElectronFakeRate.SmurfBDTGWithIPInfoElectron.r11b.skim.root", "READ");
  TFile *file = new TFile("ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "READ");


  vector<TGraphAsymmErrors*> plots;
  vector<string> labels;
  vector<Int_t> colors;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.4,0.80,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(1);




  //********************************************************************************
  //
  //********************************************************************************

  tmpLegend = new TLegend(0.4,0.80,0.93,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(0);//   tmpLegend->Clear();
  plots.clear();
  labels.clear();
  colors.clear();
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8Sample_ptThreshold35_Pt"));
  labels.push_back("Ele8");
  colors.push_back(kRed);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLSample_ptThreshold35_Pt"));
  labels.push_back("Ele8_CaloIdL_CaloIsoVL");
  colors.push_back(kBlue);
  plots.push_back((TGraphAsymmErrors*)file->Get("ElectronFakeRateDenominatorV4_Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL_ptThreshold35_Pt"));
  labels.push_back("Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL");
  colors.push_back(kGreen+2);

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
   }

  for(int i=0; i<labels.size(); ++i) {
    if (i==0) {
      plots[i]->Draw("AP");
    } else {
      plots[i]->Draw("Psame");
    }    
  }
  tmpLegend->Draw();
  cv->SaveAs("ElectronFakeRate_Triggers.gif");


}



void plotFakeRateComparison( ) {

//   plotElectronFakeRateComparison();
//   plotElectronFakeRateVsJetThresholds();
//   plotElectronFakeRateVsTriggers();
//   plotElectronFakeRateVsPileup();

//   plotMuonFakeRateComparison();
   plotMuonFakeRateVsJetThresholds();
//   plotMuonFakeRateVsPileup();
//   plotMuonFakeRateComparison();


}
