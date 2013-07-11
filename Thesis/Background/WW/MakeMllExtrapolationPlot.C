#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TLine.h"
#include "TArrow.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "HWWCuts.h"
#include "factors.h"
#include "OtherBkgScaleFactors.h"
#include "DYBkgScaleFactors.h"
#include "TopBkgScaleFactors.h"


void MakeMllExtrapolationPlot() {

  TCanvas *c1 = new TCanvas("cv","cv",800,600);
  
  TFile *_file130 = TFile::Open("/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/hww130.root");
  TTree *t130 = (TTree*) _file130->Get("tree");
  TH1F* h130 = new TH1F("h130","h130",50,0,300);
  t130->Draw("dilep.mass()>>h130");
  h130->SetLineColor(kRed);
  h130->SetLineWidth(4);
  h130->SetLineStyle(9);
  //h130->Sumw2();


  TFile *_file200 = TFile::Open("/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/hww200.root");
  TTree *t200 = (TTree*) _file200->Get("tree");
  TH1F* h200 = new TH1F("h200","h200",50,0,300);
  t200->Draw("dilep.mass()>>h200");
  h200->SetLineColor(kBlue);
  h200->SetLineWidth(4);
  h200->SetLineStyle(1);
  //h200->Sumw2();


  TFile *_file300 = TFile::Open("/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/hww300.root");
  TTree *t300 = (TTree*) _file300->Get("tree");
  TH1F* h300 = new TH1F("h300","h300",50,0,300);
  t300->Draw("dilep.mass()>>h300");
  h300->SetLineColor(kGreen+3);
  h300->SetLineWidth(4);
  h300->SetLineStyle(7);
  //h300->Sumw2();


  TFile *_fileww = TFile::Open("/data/smurf/sixie/data/Thesis/Run2011_Summer11_SmurfV7_42X/mitf-alljets/qqww.root");
  TTree *tww = (TTree*) _fileww->Get("tree");
  TH1F* hww = new TH1F("hww","hww",50,0,300);
  tww->Draw("dilep.mass()>>hww");
  hww->SetFillColor(kAzure-9);
  hww->SetFillStyle(1001);
  hww->SetLineWidth(1);
  //hww->Sumw2();


  h130->SetTitle("");
  h130->GetYaxis()->SetTitle("Fraction of Events");
  h130->GetYaxis()->SetTitleOffset(1.2);
  h130->GetXaxis()->SetTitle("m_{ll} [GeV/c^{2}]");
  h130->GetXaxis()->SetTitleOffset(1.05);
  h130->DrawNormalized("hist");
  hww->DrawNormalized("same,hist");
  h130->DrawNormalized("same,hist");
  h200->DrawNormalized("same,hist");
  h300->DrawNormalized("same,hist");
  h130->GetYaxis()->SetRangeUser(0.0,0.15);
  hww->GetYaxis()->SetRangeUser(0.0,0.15);  
  h200->GetYaxis()->SetRangeUser(0.0,0.15);
  h300->GetYaxis()->SetRangeUser(0.0,0.15);



  TLegend *l1 = new TLegend(0.65, 0.5, 0.9, 0.9);
  l1->SetTextSize(0.04);
  l1->SetBorderSize(0);
  l1->SetFillColor(kWhite);
  l1->AddEntry(h130, "m_{H}=130 GeV", "l");
  l1->AddEntry(h300, "m_{H}=300 GeV", "l");
  l1->AddEntry(h200, "m_{H}=200 GeV", "l");
  l1->AddEntry(hww, "qq#rightarrow WW", "f");
  l1->Draw();

  TLine *cutline = new TLine(100,0.12,100,0.0);
  cutline->SetLineColor(kBlack);
  cutline->SetLineWidth(4);
  cutline->Draw();

  TArrow *cutarrow = new TArrow(100,0.09,125,0.09,0.03,"|>");
  cutarrow->SetLineColor(kBlack);
  cutarrow->SetLineWidth(4);
  cutarrow->SetArrowSize(0.03);
  cutarrow->SetAngle(30);
  cutarrow->Draw();

  TPaveLabel *label = new TPaveLabel(115,0.098,165,0.108,"WW Control Region");
  label->SetFillStyle(0);
  label->SetBorderSize(0);
  label->SetTextSize(0.5);
  label->Draw();


  c1->SaveAs("DileptonMassWWControlRegion.eps");
  c1->SaveAs("DileptonMassWWControlRegion.png");


}
