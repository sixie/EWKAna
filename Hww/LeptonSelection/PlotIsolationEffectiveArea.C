//root -l EWKAna/Hww/LeptonSelection/PlotIsolationEffectiveArea.C+
 
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
#include <TF1.h>                   // 1D histograms
#include <TPaveLabel.h>                   // 1D histograms
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




//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string DoubleToString(double i) {
  char temp[100];
  sprintf(temp, "%.4f", i);
  string str = temp;
  return str;
}


void PlotIsolationEffectiveArea(string Label = "") {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
//   TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.bkg.root", "READ");
  TFile *f = new TFile("RhoPlots.root", "READ");
  
  TLegend *legend = 0;
  
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  
  TGraphAsymmErrors *tmpIso_ZMC = 0;
  TGraphAsymmErrors *tmpIso_TagAndProbe = 0;
  TGraphAsymmErrors *tmpIso_QCDMC = 0;
  TGraphAsymmErrors *tmpIso_WJetsMC = 0;
  TGraphAsymmErrors *tmpIso_JetData = 0;
  string tmpLabel;
  TPaveLabel *SlopeText_TagAndProbe;
  TPaveLabel *SlopeText_ZMC;
  TPaveLabel *SlopeText_WJetsMC;
  TPaveLabel *SlopeText_QCDMC;
  TPaveLabel *SlopeText_JetData;
  TF1 *function_TagAndProbe = 0;
  TF1 *function_ZMC = 0;
  TF1 *function_WJetsMC = 0;
  TF1 *function_QCDMC = 0;
  TF1 *function_JetData = 0;




  //*************************************************************************************************8
  // Electron CaloIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("CaloIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("CaloIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("CaloIsoElectronBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Electron ECalIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


   function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("EcalIsoAverageElectronBarrelFit_DataTagAndProbe");
   assert(function_TagAndProbe);
   tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
   SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
   SlopeText_TagAndProbe->SetBorderSize(0);
   SlopeText_TagAndProbe->SetTextSize(0.5);
   delete function_TagAndProbe;

   function_ZMC = tmpIso_ZMC->GetFunction("EcalIsoAverageElectronBarrelFit_ZMC");
   assert(function_ZMC);
   tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
   SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
   SlopeText_ZMC->SetBorderSize(0);
   SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("EcalIsoElectronBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Electron HCalIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("HcalIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("HcalIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("HcalIsoElectronBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Electron TrkIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TrkIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TrkIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TrkIsoElectronBarrel_Z.gif");
 












  //*************************************************************************************************8
  // Electron CaloIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("CaloIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("CaloIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("CaloIsoElectronEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Electron ECalIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


   function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("EcalIsoAverageElectronEndcapFit_DataTagAndProbe");
   assert(function_TagAndProbe);
   tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
   SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
   SlopeText_TagAndProbe->SetBorderSize(0);
   SlopeText_TagAndProbe->SetTextSize(0.5);
   delete function_TagAndProbe;

   function_ZMC = tmpIso_ZMC->GetFunction("EcalIsoAverageElectronEndcapFit_ZMC");
   assert(function_ZMC);
   tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
   SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
   SlopeText_ZMC->SetBorderSize(0);
   SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("EcalIsoElectronEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Electron HCalIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("HcalIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("HcalIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("HcalIsoElectronEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Electron TrkIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TrkIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TrkIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TrkIsoElectronEndcap_Z.gif");
 













  //*************************************************************************************************8
  // Electron ChargedIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoElectronBarrel_Z.gif");
 

 //*************************************************************************************************8
  // Electron ChargedNoPUIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoNoPUAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;

  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoNoPUAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoNoPUElectronBarrel_Z.gif");
 



 //*************************************************************************************************8
  // Electron NeutralHadronIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("NeutralHadronIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("NeutralHadronIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("NeutralHadronIsoElectronBarrel_Z.gif");
 


 //*************************************************************************************************8
  // Electron GammaIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("GammaIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("GammaIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("GammaIsoElectronBarrel_Z.gif");
 

 //*************************************************************************************************8
  // Electron TotalPFIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TotalPFIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TotalPFIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TotalPFIsoElectronBarrel_Z.gif");
 

 //*************************************************************************************************8
  // Electron FootprintRemovedPFIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("FootprintRemovedPFIsoAverageElectronBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("FootprintRemovedPFIsoAverageElectronBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);


  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("FootprintRemovedPFIsoElectronBarrel_Z.gif");
 


  //*************************************************************************************************8
  // Electron ChargedIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoElectronEndcap_Z.gif");
 

 //*************************************************************************************************8
  // Electron ChargedNoPUIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoNoPUAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;

  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoNoPUAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoNoPUElectronEndcap_Z.gif");
 



 //*************************************************************************************************8
  // Electron NeutralHadronIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphNeutralHadronIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("NeutralHadronIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("NeutralHadronIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("NeutralHadronIsoElectronEndcap_Z.gif");
 


 //*************************************************************************************************8
  // Electron GammaIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphGammaIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("GammaIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("GammaIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("GammaIsoElectronEndcap_Z.gif");
 

 //*************************************************************************************************8
  // Electron TotalPFIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TotalPFIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TotalPFIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TotalPFIsoElectronEndcap_Z.gif");
 

 //*************************************************************************************************8
  // Electron FootprintRemovedPFIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphFootprintRemovedPFIsoAverageElectronEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("FootprintRemovedPFIsoAverageElectronEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("FootprintRemovedPFIsoAverageElectronEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);


  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("FootprintRemovedPFIsoElectronEndcap_Z.gif");
 













  //*************************************************************************************************8
  // Muon CaloIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("CaloIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("CaloIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("CaloIsoMuonBarrel_Z.gif");

  //*************************************************************************************************8
  // Muon ECalIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


   function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("EcalIsoAverageMuonBarrelFit_DataTagAndProbe");
   assert(function_TagAndProbe);
   tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
   SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
   SlopeText_TagAndProbe->SetBorderSize(0);
   SlopeText_TagAndProbe->SetTextSize(0.5);
   delete function_TagAndProbe;

   function_ZMC = tmpIso_ZMC->GetFunction("EcalIsoAverageMuonBarrelFit_ZMC");
   assert(function_ZMC);
   tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
   SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
   SlopeText_ZMC->SetBorderSize(0);
   SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("EcalIsoMuonBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Muon HCalIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("HcalIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("HcalIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("HcalIsoMuonBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Muon TrkIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TrkIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TrkIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TrkIsoMuonBarrel_Z.gif");
 


  //*************************************************************************************************8
  // Muon CaloIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphCaloIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("CaloIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("CaloIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
  
  
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("CaloIsoMuonEndcap_Z.gif");

  //*************************************************************************************************8
  // Muon ECalIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphEcalIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


   function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("EcalIsoAverageMuonEndcapFit_DataTagAndProbe");
   assert(function_TagAndProbe);
   tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
   SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
   SlopeText_TagAndProbe->SetBorderSize(0);
   SlopeText_TagAndProbe->SetTextSize(0.5);
   delete function_TagAndProbe;

   function_ZMC = tmpIso_ZMC->GetFunction("EcalIsoAverageMuonEndcapFit_ZMC");
   assert(function_ZMC);
   tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
   SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
   SlopeText_ZMC->SetBorderSize(0);
   SlopeText_ZMC->SetTextSize(0.5);

  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("EcalIsoMuonEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Muon HCalIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphHcalIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("HcalIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("HcalIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("HcalIsoMuonEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Muon TrkIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTrkIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TrkIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TrkIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TrkIsoMuonEndcap_Z.gif");
 




  //*************************************************************************************************8
  // Muon ChargedIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoMuonBarrel_Z.gif");
 

  //*************************************************************************************************8
  // Muon ChargedIsoNoPU Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoNoPUAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoNoPUAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoNoPUMuonBarrel_Z.gif");
 


  //*************************************************************************************************8
  // Muon NeutralIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("NeutralIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("NeutralIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("NeutralIsoMuonBarrel_Z.gif");
 
  //*************************************************************************************************8
  // Muon TotalPFIso Barrel
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonBarrel_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonBarrel_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonBarrel_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonBarrel_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonBarrel_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TotalPFIsoAverageMuonBarrelFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TotalPFIsoAverageMuonBarrelFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TotalPFIsoMuonBarrel_Z.gif");
 

  //*************************************************************************************************8
  // Muon ChargedIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoMuonEndcap_Z.gif");
 

  //*************************************************************************************************8
  // Muon ChargedIsoNoPU Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphChargedIsoNoPUAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("ChargedIsoNoPUAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("ChargedIsoNoPUAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("ChargedIsoNoPUMuonEndcap_Z.gif");
 


  //*************************************************************************************************8
  // Muon NeutralIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphNeutralIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("NeutralIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("NeutralIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("NeutralIsoMuonEndcap_Z.gif");
 
  //*************************************************************************************************8
  // Muon TotalPFIso Endcap
  //*************************************************************************************************8
  tmpIso_ZMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonEndcap_ZMC");
  tmpIso_TagAndProbe = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonEndcap_DataTagAndProbe");
  tmpIso_QCDMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonEndcap_QCDMC");
  tmpIso_WJetsMC = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonEndcap_WJetsMC");
  tmpIso_JetData = (TGraphAsymmErrors*)f->Get("GraphTotalPFIsoAverageMuonEndcap_JetData");
  tmpIso_ZMC->GetXaxis()->SetTitle("Number of Primary Vertices");
  tmpIso_ZMC->GetYaxis()->SetTitle("<Isolation> [GeV]");

  tmpIso_TagAndProbe->SetFillColor(kBlue);
  tmpIso_TagAndProbe->SetFillStyle(3001);
  tmpIso_TagAndProbe->SetMarkerColor(kBlue);
  tmpIso_TagAndProbe->SetLineColor(kBlue);


  function_TagAndProbe = tmpIso_TagAndProbe->GetFunction("TotalPFIsoAverageMuonEndcapFit_DataTagAndProbe");
  assert(function_TagAndProbe);
  tmpLabel = "Data Slope : " + DoubleToString(function_TagAndProbe->GetParameter(1)) + " +- " + DoubleToString(function_TagAndProbe->GetParError(1));
  SlopeText_TagAndProbe = new TPaveLabel(0.25,0.67,0.43,0.75, tmpLabel.c_str(), "NDC");
  SlopeText_TagAndProbe->SetBorderSize(0);
  SlopeText_TagAndProbe->SetTextSize(0.5);
  delete function_TagAndProbe;
  
  function_ZMC = tmpIso_ZMC->GetFunction("TotalPFIsoAverageMuonEndcapFit_ZMC");
  assert(function_ZMC);
  tmpLabel = "MC Slope : " + DoubleToString(function_ZMC->GetParameter(1)) + " +- " + DoubleToString(function_ZMC->GetParError(1));
  SlopeText_ZMC = new TPaveLabel(0.25,0.59,0.43,0.67, tmpLabel.c_str(), "NDC");
  SlopeText_ZMC->SetBorderSize(0);
  SlopeText_ZMC->SetTextSize(0.5);
 
  assert( tmpIso_ZMC );
  assert( tmpIso_TagAndProbe ); 

  //    assert( tmpIso_QCDMC ); 
//     assert( tmpIso_WJetsMC); 
//     assert( tmpIso_JetData); 

  legend = new TLegend(0.20,0.75,0.43,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  legend->AddEntry(tmpIso_ZMC, "Z->ll MC", "P");
  legend->AddEntry(tmpIso_TagAndProbe, "Data T&P", "F");


  tmpIso_ZMC->Draw("ap");
  tmpIso_TagAndProbe->Draw("2same");
  legend->Draw();
  SlopeText_TagAndProbe->Draw();
  SlopeText_ZMC->Draw();

  cv->SaveAs("TotalPFIsoMuonEndcap_Z.gif");
 


}

