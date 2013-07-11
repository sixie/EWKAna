//root -l EWKAna/Hww/Selection/PlotMetStudies.C+

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
  TFile *f = new TFile("HwwSelectionPlots.root", "READ");
  
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;


  //**********************
  //BASELINE 2010 PU
  //**********************

  TGraphAsymmErrors *RInOutAfterPFMetdeltaPhilEtCut_mm = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2010PU_mm");
  TGraphAsymmErrors *RInOutAfterPFMetdeltaPhilEtCut_ee = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterPFMetdeltaPhilEtCut_mm->SetMarkerColor(colors[0]); 
  RInOutAfterPFMetdeltaPhilEtCut_ee->SetMarkerColor(colors[1]); 
  RInOutAfterPFMetdeltaPhilEtCut_mm->SetLineColor(colors[0]); 
  RInOutAfterPFMetdeltaPhilEtCut_ee->SetLineColor(colors[1]); 
  RInOutAfterPFMetdeltaPhilEtCut_mm->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut_ee->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut_mm->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFMetdeltaPhilEtCut_ee->GetXaxis()->SetRangeUser(0, 60);

  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut_mm, "#mu#mu final state", "LP");
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut_ee, "ee final state", "LP");
  

  RInOutAfterPFMetdeltaPhilEtCut_mm->SetTitle("");
  RInOutAfterPFMetdeltaPhilEtCut_mm->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterPFMetdeltaPhilEtCut_mm->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetdeltaPhilEtCut_mm->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetdeltaPhilEtCut_mm->Draw("AP");
  RInOutAfterPFMetdeltaPhilEtCut_ee->Draw("P");
  legend->Draw();
  cv->SaveAs(("RInOut" + label + "_MetCut_2010PU.gif").c_str());
  //**********************
  //**********************




  return;



  TGraphAsymmErrors *RInOutAfterPFMetCut = 0;
  TGraphAsymmErrors *RInOutAfterPFTrackMetCut = 0;
  TGraphAsymmErrors *RInOutAfterPFNoFwdMetCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFTrackMetPFMetCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFNoFwdMetPFMetCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = 0;
  TGraphAsymmErrors *RInOutAfterPFMetdeltaPhilEtCut = 0;
  TGraphAsymmErrors *RInOutAfterPFTrackMetdeltaPhilEtCut = 0;
  TGraphAsymmErrors *RInOutAfterPFNoFwdMetdeltaPhilEtCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = 0;
  TGraphAsymmErrors *RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = 0;








  RInOutAfterPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetCut_Z_2011PU_mm");
  RInOutAfterPFTrackMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetCut_Z_2011PU_mm");
  RInOutAfterPFNoFwdMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetCut_Z_2011PU_mm");
  RInOutAfterMinPFTrackMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetCut_Z_2011PU_mm");
  RInOutAfterMinPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetCut_Z_2011PU_mm");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut_Z_2011PU_mm");
  RInOutAfterPFMetCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetCut->SetTitle("");
  RInOutAfterPFMetCut->GetXaxis()->SetTitle("Met Cut Value [GeV]");
  RInOutAfterPFMetCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetCut->Draw("AP");
  RInOutAfterPFTrackMetCut->Draw("P");
  RInOutAfterPFNoFwdMetCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->Draw("P");
  legend->Draw();
  cv->SaveAs(("RInOut" + label + "_MetCut_Zmm_2011PU.gif").c_str());





  RInOutAfterPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetCut_Z_2011PU_ee");
  RInOutAfterPFTrackMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetCut_Z_2011PU_ee");
  RInOutAfterPFNoFwdMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetCut_Z_2011PU_ee");
  RInOutAfterMinPFTrackMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetCut_Z_2011PU_ee");
  RInOutAfterMinPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetCut_Z_2011PU_ee");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut_Z_2011PU_ee");
  RInOutAfterPFMetCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetCut->SetTitle("");
  RInOutAfterPFMetCut->GetXaxis()->SetTitle("Met Cut Value [GeV]");
  RInOutAfterPFMetCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetCut->Draw("AP");
  RInOutAfterPFTrackMetCut->Draw("P");
  RInOutAfterPFNoFwdMetCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetCut_Zee_2011PU.gif").c_str());






  RInOutAfterPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetCut_Z_2010PU_mm");
  RInOutAfterPFTrackMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetCut_Z_2010PU_mm");
  RInOutAfterPFNoFwdMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetCut_Z_2010PU_mm");
  RInOutAfterMinPFTrackMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetCut_Z_2010PU_mm");
  RInOutAfterMinPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetCut_Z_2010PU_mm");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut_Z_2010PU_mm");
  RInOutAfterPFMetCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetCut->SetTitle("");
  RInOutAfterPFMetCut->GetXaxis()->SetTitle("Met Cut Value [GeV]");
  RInOutAfterPFMetCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetCut->Draw("AP");
  RInOutAfterPFTrackMetCut->Draw("P");
  RInOutAfterPFNoFwdMetCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetCut_Zmm_2010PU.gif").c_str());



  RInOutAfterPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetCut_Z_2010PU_ee");
  RInOutAfterPFTrackMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetCut_Z_2010PU_ee");
  RInOutAfterPFNoFwdMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetCut_Z_2010PU_ee");
  RInOutAfterMinPFTrackMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetCut_Z_2010PU_ee");
  RInOutAfterMinPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetCut_Z_2010PU_ee");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut_Z_2010PU_ee");
  RInOutAfterPFMetCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetCut->SetTitle("");
  RInOutAfterPFMetCut->GetXaxis()->SetTitle("Met Cut Value [GeV]");
  RInOutAfterPFMetCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetCut->Draw("AP");
  RInOutAfterPFTrackMetCut->Draw("P");
  RInOutAfterPFNoFwdMetCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetCut_Zee_2010PU.gif").c_str());










  RInOutAfterPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterPFTrackMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut_Z_2011PU_mm");
  RInOutAfterPFMetdeltaPhilEtCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetdeltaPhilEtCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetdeltaPhilEtCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetdeltaPhilEtCut->SetTitle("");
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetdeltaPhilEtCut->Draw("AP");
  RInOutAfterPFTrackMetdeltaPhilEtCut->Draw("P");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetdeltaPhilEtCut_Zmm_2011PU.gif").c_str());



  RInOutAfterPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterPFTrackMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut_Z_2011PU_ee");
  RInOutAfterPFMetdeltaPhilEtCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetdeltaPhilEtCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetdeltaPhilEtCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetdeltaPhilEtCut->SetTitle("");
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetdeltaPhilEtCut->Draw("AP");
  RInOutAfterPFTrackMetdeltaPhilEtCut->Draw("P");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetdeltaPhilEtCut_Zee_2011PU.gif").c_str());





  RInOutAfterPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterPFTrackMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut_Z_2010PU_mm");
  RInOutAfterPFMetdeltaPhilEtCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetdeltaPhilEtCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetdeltaPhilEtCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetdeltaPhilEtCut->SetTitle("");
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetdeltaPhilEtCut->Draw("AP");
  RInOutAfterPFTrackMetdeltaPhilEtCut->Draw("P");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
  legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetdeltaPhilEtCut_Zmm_2010PU.gif").c_str());



  RInOutAfterPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterPFTrackMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFTrackMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterPFNoFwdMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = (TGraphAsymmErrors*)f->Get("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut_Z_2010PU_ee");
  RInOutAfterPFMetdeltaPhilEtCut->SetMarkerColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetMarkerColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetMarkerColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetMarkerColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->SetLineColor(colors[0]); 
  RInOutAfterPFTrackMetdeltaPhilEtCut->SetLineColor(colors[1]); 
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->SetLineColor(colors[2]); 
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetLineColor(colors[3]); 
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[4]); 
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetLineColor(colors[5]); 
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetYaxis()->SetRangeUser(0, 0.5);
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFTrackMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);
  RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->GetXaxis()->SetRangeUser(0, 60);


  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(RInOutAfterPFMetdeltaPhilEtCut, "PFMet", "LP");
  legend->AddEntry(RInOutAfterPFTrackMetdeltaPhilEtCut, "PFTrackMet", "LP");
  legend->AddEntry(RInOutAfterPFNoFwdMetdeltaPhilEtCut, "PFNoFwd", "LP");
  legend->AddEntry(RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFMet)", "LP");
  legend->AddEntry(RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut , "Min(PFNoFwdMet,PFMet)", "LP");
//   legend->AddEntry(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut, "Min(PFTrackMet,PFNoFwdMet,PFMet)", "LP");
  

  RInOutAfterPFMetdeltaPhilEtCut->SetTitle("");
  RInOutAfterPFMetdeltaPhilEtCut->GetXaxis()->SetTitle("Projected Met Cut Value [GeV]");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitle("R_{out/in}");
  RInOutAfterPFMetdeltaPhilEtCut->GetYaxis()->SetTitleOffset(1.2);
  RInOutAfterPFMetdeltaPhilEtCut->Draw("AP");
  RInOutAfterPFTrackMetdeltaPhilEtCut->Draw("P");
  RInOutAfterPFNoFwdMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->Draw("P");
  RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
//   RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->Draw("P");
legend->Draw();
cv->SaveAs(("RInOut" + label + "_MetdeltaPhilEtCut_Zee_2010PU.gif").c_str());










} 
  


void PlotPFTrackMetVsPFMet() {

  TStyle *MITStyle = new TStyle("MIT-Style","The Perfect Style for Plots ;-)");
  // gStyle = MITStyle;
  MITStyle->SetPalette(1);
//    TFile *f = new TFile("HwwSelectionPlots_ZllSelectionWithMassVetoJetVeto.root", "READ");
  TFile *f = new TFile("HwwSelectionPlots.root", "READ");
  
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;
  TH2F *tmpBkg = 0;
  TH2F *tmpSig = 0;

  tmpBkg = (TH2F*)f->Get("hPFTrackMetVsPFMet_Z_2011PU");
  tmpBkg->GetYaxis()->SetTitleOffset(1.1);
  tmpBkg->GetXaxis()->SetTitleOffset(1.05);

  tmpBkg->Draw("colz");
  cv->SaveAs("PFTrackMetVsPFMet_Z_2011PU.gif");
  tmpBkg->GetXaxis()->SetRangeUser(20,100);
  cv->SaveAs("PFTrackMetVsPFMet_Z_2011PU_AfterPFTrackMetCut.gif");



  tmpSig = (TH2F*)f->Get("hPFTrackMetVsPFMet_WW_2011PU");
  tmpSig->GetYaxis()->SetTitleOffset(1.1);
  tmpSig->GetXaxis()->SetTitleOffset(1.05);

  tmpSig->Draw("colz");
  cv->SaveAs("PFTrackMetVsPFMet_WW_2011PU.gif");
  tmpSig->GetXaxis()->SetRangeUser(20,100);
  cv->SaveAs("PFTrackMetVsPFMet_WW_2011PU_AfterPFTrackMetCut.gif");


  tmpBkg = (TH2F*)f->Get("hPFTrackMetVsPFMet_Z_2010PU");
  tmpBkg->GetYaxis()->SetTitleOffset(1.1);
  tmpBkg->GetXaxis()->SetTitleOffset(1.05);

  tmpBkg->Draw("colz");
  cv->SaveAs("PFTrackMetVsPFMet_Z_2010PU.gif");
  tmpBkg->GetXaxis()->SetRangeUser(20,100);
  cv->SaveAs("PFTrackMetVsPFMet_Z_2010PU_AfterPFTrackMetCut.gif");



  tmpSig = (TH2F*)f->Get("hPFTrackMetVsPFMet_WW_2010PU");
  tmpSig->GetYaxis()->SetTitleOffset(1.1);
  tmpSig->GetXaxis()->SetTitleOffset(1.05);

  tmpSig->Draw("colz");
  cv->SaveAs("PFTrackMetVsPFMet_WW_2010PU.gif");
  tmpSig->GetXaxis()->SetRangeUser(20,100);
  cv->SaveAs("PFTrackMetVsPFMet_WW_2010PU_AfterPFTrackMetCut.gif");


} 
  

void PlotNSigVsNBkg() {

//   string label = "2011PU_mm_AfterMinMassCutJetVetoZVetoTopVeto";
//   string label = "2011PU_mm_AfterMinMassCutJetVetoZVeto";
   string label = "2011PU_mm_AfterJetVeto";

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
//   TFile *f = new TFile("HwwSelectionPlots_ZllSelectionWithMassVetoJetVeto.root", "READ");
//   TFile *f = new TFile("HwwSelectionPlots_ZllSelectionWithJetVeto.root", "READ");
//   TFile *f = new TFile("HwwSelectionPlots_AfterMinMassCutJetVetoZVetoTopVeto.root", "READ");
//   TFile *f = new TFile("HwwSelectionPlots_AfterMinMassCutJetVetoZVeto.root", "READ");
  TFile *f = new TFile("HwwSelectionPlots_AfterJetVeto.root", "READ");
//    TFile *f = new TFile("HwwSelectionPlots.root", "READ");
  
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;
  
  
  vector<string> MetNames;    
//   MetNames.push_back("hPFMet");
//   MetNames.push_back("hPFTrackMet");   
//   MetNames.push_back("hPFNoFwdMet");
//   MetNames.push_back("hMinPFTrackMetPFMet");
//   MetNames.push_back("hMinPFNoFwdMetPFMet");
//   MetNames.push_back("hMinPFTrackMetPFNoFwdMetPFMet");

  MetNames.push_back("hPFMetDeltaPhilEt");
  MetNames.push_back("hPFTrackMetDeltaPhilEt");   
  MetNames.push_back("hPFNoFwdMetDeltaPhilEt");
  MetNames.push_back("hMinPFTrackMetPFMetDeltaPhilEt");
  MetNames.push_back("hMinPFNoFwdMetPFMetDeltaPhilEt");
//    MetNames.push_back("hMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt");

  vector<string> labels;
  labels.push_back("PFMet");
  labels.push_back("PFTrackMet");   
  labels.push_back("PFNoFwdMet");
  labels.push_back("Min(PFTrackMet,PFMet)");
  labels.push_back("Min(PFNoFwdMet,PFMet)");
//    labels.push_back("Min(PFTrackMet,PFNoFwdMet,PFMet)");


// //Projected Met Study
//   vector<string> MetNames;    
//   MetNames.push_back("hPFMet");
//   MetNames.push_back("hPFMetDeltaPhilEt");
//   MetNames.push_back("hPFTrackMet");   
//   MetNames.push_back("hPFTrackMetDeltaPhilEt");   
//   MetNames.push_back("hMinPFTrackMetPFMet");
//   MetNames.push_back("hMinPFTrackMetPFMetDeltaPhilEt");

//   vector<string> labels;
//   labels.push_back("PFMet");
//   labels.push_back("PF ProjMet ");
//   labels.push_back("PFTrackMet ");   
//   labels.push_back("PFTrack ProjMet ");   
//   labels.push_back("Min(PFMet,PFTrackMet)");
//   labels.push_back("Min(PF ProjMet,PFTrack ProjMet)");
// //   labels.push_back("Proj Min(PFMet,PFTrackMet)");


  vector<string> ProjectedMetNames;
  vector<TH1F*> signalMet;
  vector<TH1F*> bkgMet;
  for (UInt_t i=0; i<MetNames.size(); ++i) {
    TH1F *tmpSignal = (TH1F*)f->Get((MetNames[i]+"_WW_2011PU_mm").c_str());
    signalMet.push_back(tmpSignal);
    TH1F *tmpBkg = (TH1F*)f->Get((MetNames[i]+"_Z_2011PU_mm").c_str());
    bkgMet.push_back(tmpBkg);
  }

    

//     Double_t maxy = Signal->GetMaximum();
//     if (Bkg->GetMaximum() > maxy) maxy = Bkg->GetMaximum(); maxy = maxy * 1.2;
//     Signal->SetMaximum(maxy);
//     Signal->Draw("hist");
//     Bkg->Draw("histsame");
    
  vector<TGraphAsymmErrors*> Met_NSigVsCut;
  vector<TGraphAsymmErrors*> Met_NSigVsNBkg;
  vector<TGraphAsymmErrors*> Met_SigEffVsBkgEff;
    

  for (UInt_t i=0; i<MetNames.size(); ++i) {

    //Make Met Plots
    const int nPoints = signalMet[i]->GetXaxis()->GetNbins();
    double cutValue[nPoints];
    double cutValueErr[nPoints];
    double NSig[nPoints];
    double NSigErr[nPoints];
    double NBkg[nPoints];
    double NBkgErr[nPoints];
    double SigEff[nPoints];
    double BkgEff[nPoints];

    double NSigTotal = 0;
    double NBkgTotal = 0;
    for (UInt_t q=0; q < nPoints+2; ++q) {
      NSigTotal += signalMet[i]->GetBinContent(q);
      NBkgTotal += bkgMet[i]->GetBinContent(q);
    }   

    for(UInt_t b=0; b < nPoints; ++b) {
      cutValue[b] = signalMet[i]->GetXaxis()->GetBinCenter(b);
      cutValueErr[b] = 0;
      Double_t nsig = 0;
      Double_t nsigErrSqr = 0;
      Double_t nbkg = 0;
      Double_t nbkgErrSqr = 0;
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalMet[i]->GetBinContent(q);
        nsigErrSqr += pow(signalMet[i]->GetBinError(q),2);
        nbkg += bkgMet[i]->GetBinContent(q);
        nbkgErrSqr += pow(bkgMet[i]->GetBinError(q),2);
      }
      NSig[b] = nsig;
      NSigErr[b] = 0; TMath::Sqrt(nsigErrSqr);
      NBkg[b] = nbkg;
      NBkgErr[b] = 0; TMath::Sqrt(nbkgErrSqr);
      SigEff[b] = nsig / NSigTotal;
      BkgEff[b] = nbkg / NBkgTotal;
      cout << b << " : " << nsig << " , " << nbkg << endl;
    }
    TGraphAsymmErrors *tmpNSigVsCut = new TGraphAsymmErrors (nPoints, cutValue, NSig, cutValueErr, cutValueErr, NSigErr, NSigErr);
    tmpNSigVsCut->SetTitle("");
    tmpNSigVsCut->GetXaxis()->SetTitle("Met Cut [GeV]");
    tmpNSigVsCut->GetYaxis()->SetTitle("Number of WW Signal / fb^{-1}");
    tmpNSigVsCut->GetYaxis()->SetTitleOffset(1.1);
    tmpNSigVsCut->GetXaxis()->SetTitleOffset(1.05);



    TGraphAsymmErrors *tmpNSigVsNBkg = new TGraphAsymmErrors (nPoints, NBkg, NSig,NBkgErr, NBkgErr, NSigErr, NSigErr );
    tmpNSigVsNBkg->SetTitle("");
    tmpNSigVsNBkg->GetXaxis()->SetTitle("Number of DY Bkg / fb^{-1}");
    tmpNSigVsNBkg->GetYaxis()->SetTitle("Number of WW Signal / fb^{-1}");
    tmpNSigVsNBkg->GetYaxis()->SetTitleOffset(1.1);
    tmpNSigVsNBkg->GetXaxis()->SetTitleOffset(1.05);


    TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff,NBkgErr, NBkgErr, NSigErr, NSigErr );
    
    tmpSigEffVsBkgEff->SetTitle("");
    tmpSigEffVsBkgEff->GetXaxis()->SetTitle("DY Bkg Eff");
    tmpSigEffVsBkgEff->GetYaxis()->SetTitle("WW Signal Eff");
    tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
    tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);

    Met_NSigVsCut.push_back(tmpNSigVsCut);
    Met_NSigVsNBkg.push_back(tmpNSigVsNBkg);   
    Met_SigEffVsBkgEff.push_back(tmpSigEffVsBkgEff);
  }

  //     NSigVsCut->Draw("AP");
//     cv->SaveAs("NSigVsCut_AfterProjMetCut_PFMet.gif");
  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<MetNames.size(); ++i) {
    legend->AddEntry(Met_NSigVsCut[i],labels[i].c_str(), "LP");
    Met_NSigVsCut[i]->SetMarkerColor(colors[i]);
    Met_NSigVsCut[i]->SetMarkerSize(0.5);
    Met_NSigVsCut[i]->SetLineColor(colors[i]);
    
    Met_NSigVsCut[i]->GetYaxis()->SetRangeUser(0,150);    
    Met_NSigVsCut[i]->GetXaxis()->SetRangeUser(0,80);    
    if (i==0) {
      Met_NSigVsCut[i]->Draw("AP");
    } else {
      Met_NSigVsCut[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("NSigVsCut_CompareDifferentMetCut" + label + ".gif").c_str());



  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<MetNames.size(); ++i) {
    legend->AddEntry(Met_NSigVsNBkg[i],labels[i].c_str(), "LP");
    Met_NSigVsNBkg[i]->SetMarkerColor(colors[i]);
    Met_NSigVsNBkg[i]->SetLineColor(colors[i]);
    Met_NSigVsNBkg[i]->SetMarkerSize(0.5);
   
    Met_NSigVsNBkg[i]->GetYaxis()->SetRangeUser(0,140);    
    Met_NSigVsNBkg[i]->GetXaxis()->SetRangeUser(0,200);    
    if (i==0) {
      Met_NSigVsNBkg[i]->Draw("AP");
    } else {
      Met_NSigVsNBkg[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("NSigVsNBkg_CompareDifferentMetCut" + label + ".gif").c_str());
    



  legend = new TLegend(0.43,0.70,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<MetNames.size(); ++i) {
    legend->AddEntry(Met_SigEffVsBkgEff[i],labels[i].c_str(), "LP");
    Met_SigEffVsBkgEff[i]->SetMarkerColor(colors[i]);
    Met_SigEffVsBkgEff[i]->SetLineColor(colors[i]);
    Met_SigEffVsBkgEff[i]->SetMarkerSize(0.5);
   
    Met_SigEffVsBkgEff[i]->GetYaxis()->SetRangeUser(0,1.0);    
    Met_SigEffVsBkgEff[i]->GetXaxis()->SetRangeUser(0,0.01);    
    if (i==0) {
      Met_SigEffVsBkgEff[i]->Draw("AP");
    } else {
      Met_SigEffVsBkgEff[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_CompareDifferentMetCut" + label + ".gif").c_str());
    


  }

void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  for (int b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) { norm += hist->GetBinContent(b); }
  for (int b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }  
}



void CompareZllMet() {

    TFile *f = new TFile("HwwSelectionPlots_ZllSelection.root", "READ");

    TCanvas *cv = new TCanvas("cv","cv", 800,600);
    TLegend *legend = 0;
    legend = new TLegend(0.43,0.70,0.93,0.90);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(1);
    
    string filename;
    vector<string> histNames;
    
//     filename = "pfMetComparison.gif";
//     histNames.push_back("hPFMet_Z_2010PULessThan2");
//     histNames.push_back("hPFMet_Z_2010PUMoreThan4");
//     histNames.push_back("hPFMet_Z_2011PULessThan2");
//     histNames.push_back("hPFMet_Z_2011PU5To7");
//     histNames.push_back("hPFMet_Z_2011MoreThan15");

//     filename = "pfMetComparison_WW.gif";
//     histNames.push_back("hPFMet_WW_2010PULessThan2");
//     histNames.push_back("hPFMet_WW_2010PUMoreThan4");
//     histNames.push_back("hPFMet_WW_2011PULessThan2");
//     histNames.push_back("hPFMet_WW_2011PU5To7");
//     histNames.push_back("hPFMet_WW_2011MoreThan15");



//     filename = "pfTrackMetComparison.gif";
//     histNames.push_back("hPFTrackMet_Z_2010PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2010PUMoreThan4");
//     histNames.push_back("hPFTrackMet_Z_2011PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2011PU5To7");
//     histNames.push_back("hPFTrackMet_Z_2011MoreThan15");

//     filename = "pfNoFwdMetComparison.gif";
//     histNames.push_back("hPFNoFwdMet_Z_2010PULessThan2");
//     histNames.push_back("hPFNoFwdMet_Z_2010PUMoreThan4");
//     histNames.push_back("hPFNoFwdMet_Z_2011PULessThan2");
//     histNames.push_back("hPFNoFwdMet_Z_2011PU5To7");
//     histNames.push_back("hPFNoFwdMet_Z_2011MoreThan15");

//     filename = "MinPFTrackMetPFMetComparison.gif";
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2010PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2010PUMoreThan4");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011PU5To7");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011MoreThan15");

//     filename = "DifferentMetOptions_2010PU.gif";
//     histNames.push_back("hPFMet_Z_2010PULessThan2");
//     histNames.push_back("hPFMet_Z_2010PUMoreThan4");
//     histNames.push_back("hPFTrackMet_Z_2010PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2010PUMoreThan4");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2010PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2010PUMoreThan4");

//     filename = "DifferentMetOptions_OutOfTimePUEffect.gif";
//     histNames.push_back("hPFMet_Z_2010PULessThan2");
//     histNames.push_back("hPFMet_Z_2011PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2010PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2011PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2010PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011PULessThan2");

//     filename = "DifferentMetOptions_2011PU.gif";
//     histNames.push_back("hPFMet_Z_2011PULessThan2");
//     histNames.push_back("hPFMet_Z_2011MoreThan15");
//     histNames.push_back("hPFTrackMet_Z_2011PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2011MoreThan15");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011PULessThan2");
//     histNames.push_back("hMinPFTrackMetPFMet_Z_2011MoreThan15");


//     filename = "pfMet_SignalVsBkg.gif";
//     histNames.push_back("hPFMet_Z_2010PULessThan2");
//     histNames.push_back("hPFMet_WW_2010PULessThan2");
//     histNames.push_back("hPFMet_Z_2011MoreThan15");
//     histNames.push_back("hPFMet_WW_2011MoreThan15");

//     filename = "pfTrackMet_SignalVsBkg.gif";
//     histNames.push_back("hPFTrackMet_Z_2010PULessThan2");
//     histNames.push_back("hPFTrackMet_WW_2010PULessThan2");
//     histNames.push_back("hPFTrackMet_Z_2011MoreThan15");
//     histNames.push_back("hPFTrackMet_WW_2011MoreThan15");

//     filename = "pfNoFwdMet_SignalVsBkg.gif";
//     histNames.push_back("hPFNoFwdMet_Z_2010PULessThan2");
//     histNames.push_back("hPFNoFwdMet_WW_2010PULessThan2");
//     histNames.push_back("hPFNoFwdMet_Z_2011MoreThan15");
//     histNames.push_back("hPFNoFwdMet_WW_2011MoreThan15");

    filename = "MinPFTrackMetPFMet_SignalVsBkg.gif";
    histNames.push_back("hMinPFTrackMetPFMet_Z_2010PULessThan2");
    histNames.push_back("hMinPFTrackMetPFMet_WW_2010PULessThan2");
    histNames.push_back("hMinPFTrackMetPFMet_Z_2011MoreThan15");
    histNames.push_back("hMinPFTrackMetPFMet_WW_2011MoreThan15");

    vector<Int_t> colors;
    colors.push_back(kRed);
    colors.push_back(kBlue);
    colors.push_back(kMagenta);
    colors.push_back(kCyan);
    colors.push_back(kBlack);
    colors.push_back(kGreen);
    vector<TH1F*> hists;

    for(UInt_t i=0; i<histNames.size(); ++i) {
      TH1F *tmp1 = (TH1F*)f->Get(histNames[i].c_str()); cout << i << endl; assert(tmp1);
      TH1F *tmp = (TH1F*)tmp1->Clone(histNames[i].c_str());
      NormalizeHist(tmp);
      tmp->GetYaxis()->SetTitleOffset(1.1);
      tmp->GetXaxis()->SetTitleOffset(1.05);
      tmp->SetLineColor(colors[i]);
      tmp->SetMarkerColor(colors[i]);
      hists.push_back(tmp);
      legend->AddEntry(tmp, histNames[i].c_str(), "LP");
    }

    Double_t maxy = hists[0]->GetMaximum();
    for(UInt_t i=0; i<histNames.size(); ++i) {
      if (hists[i]->GetMaximum() > maxy) maxy = hists[i]->GetMaximum();
    }
    maxy = maxy*1.2;
    
    
    maxy = maxy*5;
    for(UInt_t i=0; i<histNames.size(); ++i) {
      hists[i]->SetMaximum(maxy);
      if (i==0) {        
        hists[i]->Draw("hist");
      } else {
        hists[i]->Draw("histsame");
      }
    }  
    legend->Draw();

    cv->SetLogy();
    cv->SaveAs(filename.c_str());
}


void PlotMetStudies() {

//   CompareZllMet();
    PlotNSigVsNBkg();
//   PlotPFTrackMetVsPFMet();
//         PlotRInOutRatio();

}
