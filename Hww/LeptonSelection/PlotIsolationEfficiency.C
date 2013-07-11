//root -l -b EWKAna/Hww/LeptonSelection/PlotIsolationEfficiency.C+
 
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
#include "TLine.h"
#include "TStyle.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
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


void PlotIsolationEfficiency(string Label = "" ) {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kBlue);
  colors.push_back(kRed);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.root", "READ");
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
  TPaveLabel *LabelText;
  TLegend *legend = 0;
  

//   //----------------------------------------------------------------------------------------------------
//   // 
//   //====================================================================================================
//   TGraphAsymmErrors* EfficiencyDataTagAndProbe;
//   TGraphAsymmErrors* EfficiencyMCTagAndProbe;
  
//   EfficiencyDataTagAndProbe = (TGraphAsymmErrors*)f->Get("ElectronIsolationBarrelEffVsNVertices_DataTagAndProbe");
//   EfficiencyMCTagAndProbe = (TGraphAsymmErrors*)f->Get("ElectronIsolationBarrelEffVsNVertices_ZMCTagAndProbe");

//   legend = new TLegend(0.20,0.25,0.60,0.40);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(1);
  
//   legend->AddEntry(EfficiencyDataTagAndProbe, "Data", "LP");
//   legend->AddEntry(EfficiencyMCTagAndProbe, "MC", "LP");
        
//   EfficiencyDataTagAndProbe->SetFillStyle(3001);
//   EfficiencyDataTagAndProbe->SetFillColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetLineColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyDataTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyDataTagAndProbe->GetYaxis()->SetTitleOffset(1.2);
//   EfficiencyMCTagAndProbe->SetFillStyle(3001);
//   EfficiencyMCTagAndProbe->SetFillColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetLineColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyMCTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyMCTagAndProbe->GetYaxis()->SetTitleOffset(1.2);

//   EfficiencyDataTagAndProbe->Draw("AP");  
//   EfficiencyMCTagAndProbe->Draw("Psame");
//   legend->Draw();
  
//   cv->SaveAs("ElectronIsolationBarrelEffVsNVertices_TagAndProbe.gif");
  



//   EfficiencyDataTagAndProbe = (TGraphAsymmErrors*)f->Get("ElectronIsolationEndcapEffVsNVertices_DataTagAndProbe");
//   EfficiencyMCTagAndProbe = (TGraphAsymmErrors*)f->Get("ElectronIsolationEndcapEffVsNVertices_ZMCTagAndProbe");

//   legend = new TLegend(0.20,0.25,0.60,0.40);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(1);
  
//   legend->AddEntry(EfficiencyDataTagAndProbe, "Data", "LP");
//   legend->AddEntry(EfficiencyMCTagAndProbe, "MC", "LP");
        
//   EfficiencyDataTagAndProbe->SetFillStyle(3001);
//   EfficiencyDataTagAndProbe->SetFillColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetLineColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyDataTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyDataTagAndProbe->GetYaxis()->SetTitleOffset(1.2);
//   EfficiencyMCTagAndProbe->SetFillStyle(3001);
//   EfficiencyMCTagAndProbe->SetFillColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetLineColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyMCTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyMCTagAndProbe->GetYaxis()->SetTitleOffset(1.2);

//   EfficiencyDataTagAndProbe->Draw("AP");  
//   EfficiencyMCTagAndProbe->Draw("Psame");
//   legend->Draw();
  
//   cv->SaveAs("ElectronIsolationEndcapEffVsNVertices_TagAndProbe.gif");
  


//   EfficiencyDataTagAndProbe = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_DataTagAndProbe");
//   EfficiencyMCTagAndProbe = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_ZMCTagAndProbe");

//   assert(EfficiencyDataTagAndProbe);
//   assert(EfficiencyMCTagAndProbe);

//   legend = new TLegend(0.20,0.25,0.60,0.40);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(1);
  
//   legend->AddEntry(EfficiencyDataTagAndProbe, "Data", "LP");
//   legend->AddEntry(EfficiencyMCTagAndProbe, "MC", "LP");
        
//   EfficiencyDataTagAndProbe->SetFillStyle(3001);
//   EfficiencyDataTagAndProbe->SetFillColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetLineColor(colors[0]);
//   EfficiencyDataTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyDataTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyDataTagAndProbe->GetYaxis()->SetTitleOffset(1.2);
//   EfficiencyDataTagAndProbe->GetYaxis()->SetRangeUser(0.85, 1.0);
//   EfficiencyMCTagAndProbe->SetFillStyle(3001);
//   EfficiencyMCTagAndProbe->SetFillColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetLineColor(colors[1]);
//   EfficiencyMCTagAndProbe->SetMarkerSize(1.0);
//   EfficiencyMCTagAndProbe->GetXaxis()->SetTitleOffset(1.05);
//   EfficiencyMCTagAndProbe->GetYaxis()->SetTitleOffset(1.2);
//   EfficiencyMCTagAndProbe->GetYaxis()->SetRangeUser(0.85, 1.0);

//   EfficiencyDataTagAndProbe->Draw("AP");  
//   EfficiencyMCTagAndProbe->Draw("Psame");
//   legend->Draw();
  
//   cv->SaveAs("MuonIsolationEffVsNVertices_TagAndProbe.gif");
  


  //----------------------------------------------------------------------------------------------------
  // 
  //====================================================================================================
  TGraphAsymmErrors* ElectronIsolationEfficiencyHWWPt20;
  TGraphAsymmErrors* ElectronIsolationEfficiencyHWWPt10To20;
  TGraphAsymmErrors* ElectronIsolationEfficiencyHWWPt35To50;
  TLine *lowerline;
  TLine *upperline;
  
  ElectronIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("ElectronIsolationEffVsNVertices_HWW130_PtSelected");
  ElectronIsolationEfficiencyHWWPt10To20 = (TGraphAsymmErrors*)f->Get("ElectronIsolationEffVsNVertices_HWW130_Pt15To20");
  ElectronIsolationEfficiencyHWWPt35To50 = (TGraphAsymmErrors*)f->Get("ElectronIsolationEffVsNVertices_HWW130_Pt30To50");

  assert(ElectronIsolationEfficiencyHWWPt20);
  assert(ElectronIsolationEfficiencyHWWPt10To20);
  assert(ElectronIsolationEfficiencyHWWPt35To50);

  legend = new TLegend(0.20,0.25,0.60,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  
  legend->AddEntry(ElectronIsolationEfficiencyHWWPt20, "All Selected", "LP");
  legend->AddEntry(ElectronIsolationEfficiencyHWWPt10To20, "Pt in [15,20] GeV", "LP");
  legend->AddEntry(ElectronIsolationEfficiencyHWWPt35To50, "Pt in [30,50] GeV", "LP");
       
  ElectronIsolationEfficiencyHWWPt20->SetFillStyle(3001);
  ElectronIsolationEfficiencyHWWPt20->SetFillColor(colors[0]);
  ElectronIsolationEfficiencyHWWPt20->SetMarkerColor(colors[0]);
  ElectronIsolationEfficiencyHWWPt20->SetLineColor(colors[0]);
  ElectronIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
  ElectronIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);
  ElectronIsolationEfficiencyHWWPt20->GetYaxis()->SetRangeUser(0.30, 1.0);

  ElectronIsolationEfficiencyHWWPt10To20->SetFillStyle(3001);
  ElectronIsolationEfficiencyHWWPt10To20->SetFillColor(colors[1]);
  ElectronIsolationEfficiencyHWWPt10To20->SetMarkerColor(colors[1]);
  ElectronIsolationEfficiencyHWWPt10To20->SetLineColor(colors[1]);
  ElectronIsolationEfficiencyHWWPt10To20->SetMarkerSize(1.0);
  ElectronIsolationEfficiencyHWWPt10To20->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEfficiencyHWWPt10To20->GetYaxis()->SetTitleOffset(1.2);
  ElectronIsolationEfficiencyHWWPt10To20->GetYaxis()->SetRangeUser(0.30, 1.0);

  ElectronIsolationEfficiencyHWWPt35To50->SetFillStyle(3001);
  ElectronIsolationEfficiencyHWWPt35To50->SetFillColor(colors[2]);
  ElectronIsolationEfficiencyHWWPt35To50->SetMarkerColor(colors[2]);
  ElectronIsolationEfficiencyHWWPt35To50->SetLineColor(colors[2]);
  ElectronIsolationEfficiencyHWWPt35To50->SetMarkerSize(1.0);
  ElectronIsolationEfficiencyHWWPt35To50->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEfficiencyHWWPt35To50->GetYaxis()->SetTitleOffset(1.2);
  ElectronIsolationEfficiencyHWWPt35To50->GetYaxis()->SetRangeUser(0.30, 1.0);


  ElectronIsolationEfficiencyHWWPt20->Draw("AP");  
  ElectronIsolationEfficiencyHWWPt10To20->Draw("Psame");
  ElectronIsolationEfficiencyHWWPt35To50->Draw("Psame");
  legend->Draw();
  lowerline = new TLine(0.15,0.77,0.95,0.77); 
  lowerline->SetLineWidth(3.0);
  lowerline->SetLineStyle(2);
  lowerline->DrawLineNDC(0.15,0.77,0.95,0.77);
  upperline = new TLine(0.15,0.835,0.95,0.835); 
  upperline->SetLineWidth(3.0);
  upperline->SetLineStyle(2);
  upperline->DrawLineNDC(0.15,0.835,0.95,0.835);
//   TPaveLabel *LabelTextUpper = new TPaveLabel(0.15,0.86,0.95,0.86, "0.92" );
//   LabelTextUpper->Draw();
//   TPaveLabel *LabelTextLower = new TPaveLabel(0.15,0.78,0.95,0.78, "0.83" );
//   LabelTextLower->Draw();
  cv->SaveAs("ElectronIsolationVsNVertices_HWW130.gif");
  cv->SaveAs("ElectronIsolationVsNVertices_HWW130.pdf");
  

  TGraphAsymmErrors* MuonIsolationEfficiencyHWWPt20;
  TGraphAsymmErrors* MuonIsolationEfficiencyHWWPt10To20;
  TGraphAsymmErrors* MuonIsolationEfficiencyHWWPt35To50;
  
  MuonIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_HWW130_PtSelected");
  MuonIsolationEfficiencyHWWPt10To20 = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_HWW130_Pt10To20");
  MuonIsolationEfficiencyHWWPt35To50 = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_HWW130_Pt30To50");
  assert(MuonIsolationEfficiencyHWWPt20);
  assert(MuonIsolationEfficiencyHWWPt10To20);
  assert(MuonIsolationEfficiencyHWWPt35To50);

  legend = new TLegend(0.20,0.25,0.60,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  
  legend->AddEntry(MuonIsolationEfficiencyHWWPt20, "All Selected", "LP");
  legend->AddEntry(MuonIsolationEfficiencyHWWPt10To20, "Pt in [10,20] GeV", "LP");
  legend->AddEntry(MuonIsolationEfficiencyHWWPt35To50, "Pt in [30,50] GeV", "LP");
       
  MuonIsolationEfficiencyHWWPt20->SetFillStyle(3001);
  MuonIsolationEfficiencyHWWPt20->SetFillColor(colors[0]);
  MuonIsolationEfficiencyHWWPt20->SetMarkerColor(colors[0]);
  MuonIsolationEfficiencyHWWPt20->SetLineColor(colors[0]);
  MuonIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
  MuonIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);

  MuonIsolationEfficiencyHWWPt10To20->SetFillStyle(3001);
  MuonIsolationEfficiencyHWWPt10To20->SetFillColor(colors[1]);
  MuonIsolationEfficiencyHWWPt10To20->SetMarkerColor(colors[1]);
  MuonIsolationEfficiencyHWWPt10To20->SetLineColor(colors[1]);
  MuonIsolationEfficiencyHWWPt10To20->SetMarkerSize(1.0);
  MuonIsolationEfficiencyHWWPt10To20->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEfficiencyHWWPt10To20->GetYaxis()->SetTitleOffset(1.2);

  MuonIsolationEfficiencyHWWPt35To50->SetFillStyle(3001);
  MuonIsolationEfficiencyHWWPt35To50->SetFillColor(colors[2]);
  MuonIsolationEfficiencyHWWPt35To50->SetMarkerColor(colors[2]);
  MuonIsolationEfficiencyHWWPt35To50->SetLineColor(colors[2]);
  MuonIsolationEfficiencyHWWPt35To50->SetMarkerSize(1.0);
  MuonIsolationEfficiencyHWWPt35To50->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEfficiencyHWWPt35To50->GetYaxis()->SetTitleOffset(1.2);


  MuonIsolationEfficiencyHWWPt20->Draw("AP");  
  MuonIsolationEfficiencyHWWPt10To20->Draw("Psame");
  MuonIsolationEfficiencyHWWPt35To50->Draw("Psame");
  legend->Draw();
  
  lowerline = new TLine(0.15,0.82,0.95,0.82); 
  lowerline->SetLineWidth(3.0);
  lowerline->SetLineStyle(2);
  lowerline->DrawLineNDC(0.15,0.82,0.95,0.82);
  upperline = new TLine(0.15,0.86,0.95,0.86); 
  upperline->SetLineWidth(3.0);
  upperline->SetLineStyle(2);
  upperline->DrawLineNDC(0.15,0.86,0.95,0.86);

  cv->SaveAs("MuonIsolationVsNVertices_HWW130.gif");
  cv->SaveAs("MuonIsolationVsNVertices_HWW130.pdf");
  




//   //----------------------------------------------------------------------------------------------------
//   // Corrections
//   //====================================================================================================
// //   TGraphAsymmErrors* ElectronIsolationEfficiencyHWWPt20;
//   TGraphAsymmErrors* ElectronFastjetCorrectedIsolationEfficiencyHWWPt20;
//   TGraphAsymmErrors* ElectronIsolationEfficiencyHWWPt35To50;
//   TLine *lowerline;
//   TLine *upperline;
  
//   ElectronIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("ElectronIsolationEffVsNVertices_HWW130_Pt20");
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("ElectronFastjetCorrectedIsolationEffVsNVertices_HWW130_Pt20");

//   legend = new TLegend(0.20,0.25,0.60,0.40);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(1);
  
//   legend->AddEntry(ElectronIsolationEfficiencyHWWPt20, "No corrections", "LP");
//   legend->AddEntry(ElectronFastjetCorrectedIsolationEfficiencyHWWPt20, "Fastjet corrected", "LP");
       
//   ElectronIsolationEfficiencyHWWPt20->SetFillStyle(3001);
//   ElectronIsolationEfficiencyHWWPt20->SetFillColor(colors[0]);
//   ElectronIsolationEfficiencyHWWPt20->SetMarkerColor(colors[0]);
//   ElectronIsolationEfficiencyHWWPt20->SetLineColor(colors[0]);
//   ElectronIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
//   ElectronIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
//   ElectronIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);
//   ElectronIsolationEfficiencyHWWPt20->GetYaxis()->SetRangeUser(0.75, 1.0);

//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->SetFillStyle(3001);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->SetFillColor(colors[1]);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->SetMarkerColor(colors[1]);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->SetLineColor(colors[1]);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->GetYaxis()->SetRangeUser(0.75, 1.0);
 


//   ElectronIsolationEfficiencyHWWPt20->Draw("AP");  
//   ElectronFastjetCorrectedIsolationEfficiencyHWWPt20->Draw("Psame");
//   legend->Draw();
//   lowerline = new TLine(0.15,0.56,0.95,0.56); 
//   lowerline->SetLineWidth(3.0);
//   lowerline->SetLineStyle(2);
//   lowerline->DrawLineNDC(0.15,0.56,0.95,0.56);
//   upperline = new TLine(0.15,0.72,0.95,0.72); 
//   upperline->SetLineWidth(3.0);
//   upperline->SetLineStyle(2);
//   upperline->DrawLineNDC(0.15,0.72,0.95,0.72);
// //   TPaveLabel *LabelTextUpper = new TPaveLabel(0.15,0.86,0.95,0.86, "0.92" );
// //   LabelTextUpper->Draw();
// //   TPaveLabel *LabelTextLower = new TPaveLabel(0.15,0.78,0.95,0.78, "0.83" );
// //   LabelTextLower->Draw();
//   cv->SaveAs("ElectronIsolationVsNVertices_HWW130_FastjetCorrection.gif");
  





// //   TGraphAsymmErrors* MuonIsolationEfficiencyHWWPt20;
//   TGraphAsymmErrors* MuonFastjetCorrectedIsolationEfficiencyHWWPt20;
  
//   MuonIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("MuonIsolationEffVsNVertices_HWW130_Pt20");
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20 = (TGraphAsymmErrors*)f->Get("MuonFastjetCorrectedIsolationEffVsNVertices_HWW130_Pt20");

//   legend = new TLegend(0.20,0.25,0.60,0.40);
//   legend->SetTextSize(0.03);
//   legend->SetBorderSize(1);
  
//   legend->AddEntry(MuonIsolationEfficiencyHWWPt20, "No correction", "LP");
//   legend->AddEntry(MuonFastjetCorrectedIsolationEfficiencyHWWPt20, "Fastjet corrected", "LP");
       
//   MuonIsolationEfficiencyHWWPt20->SetFillStyle(3001);
//   MuonIsolationEfficiencyHWWPt20->SetFillColor(colors[0]);
//   MuonIsolationEfficiencyHWWPt20->SetMarkerColor(colors[0]);
//   MuonIsolationEfficiencyHWWPt20->SetLineColor(colors[0]);
//   MuonIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
//   MuonIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
//   MuonIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);
//   MuonIsolationEfficiencyHWWPt20->GetYaxis()->SetRangeUser(0.80, 1.0);

//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->SetFillStyle(3001);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->SetFillColor(colors[1]);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->SetMarkerColor(colors[1]);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->SetLineColor(colors[1]);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->SetMarkerSize(1.0);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->GetXaxis()->SetTitleOffset(1.05);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->GetYaxis()->SetTitleOffset(1.2);
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->GetYaxis()->SetRangeUser(0.80, 1.0);

//   MuonIsolationEfficiencyHWWPt20->Draw("AP");  
//   MuonFastjetCorrectedIsolationEfficiencyHWWPt20->Draw("Psame");
//   legend->Draw();
  
//   lowerline = new TLine(0.15,0.59,0.95,0.59); 
//   lowerline->SetLineWidth(3.0);
//   lowerline->SetLineStyle(2);
//   lowerline->DrawLineNDC(0.15,0.59,0.95,0.59);
//   upperline = new TLine(0.15,0.680,0.95,0.680); 
//   upperline->SetLineWidth(3.0);
//   upperline->SetLineStyle(2);
//   upperline->DrawLineNDC(0.15,0.680,0.95,0.680);

//   cv->SaveAs("MuonIsolationVsNVertices_HWW130_FastjetCorrection.gif");
  



}


