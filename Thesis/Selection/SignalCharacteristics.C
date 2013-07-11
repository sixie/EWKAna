//root -l EWKAna/Hww/Thesis/SignalCharacteristics/SignalCharacteristics.C+\(\)
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
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree


// RooFit headers

#include "TFile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TLegend.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================
//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwGenNtuplerMod");
    if (!dir) dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

void NormalizeHiggsPt(TH1D *hist) {
  Double_t norm = 0;

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / (hist->GetXaxis()->GetBinUpEdge(b) - hist->GetXaxis()->GetBinLowEdge(b)));
    hist->SetBinError(b,hist->GetBinError(b) / (hist->GetXaxis()->GetBinUpEdge(b) - hist->GetXaxis()->GetBinLowEdge(b)));
  }

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

void NormalizeHist(TH1D *hist) {
  Double_t norm = 0;

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


void SignalCharacteristics() {


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *HiggsPt_HWW160 = new TH1D("HiggsPt_HWW160", ";Higgs p_{T} [GeV/c];Fraction of Events", 80, 0, 200);
  TH1D *DeltaPhi_HWW160 = new TH1D("DeltaPhi_HWW160", ";#Delta#phi(l,l)[Degrees];Fraction of Events", 25, 0, 180);
  TH1D *Met_HWW160 = new TH1D("Met_HWW160", ";Missing Transverse Energy [GeV/c];Fraction of Events", 60, 0, 120);
  TH1D *PtMin_HWW120 = new TH1D("PtMin_HWW120", ";Trailing Lepton p_{T} [GeV/c];Fraction of Events", 50, 0, 100);
  TH1D *PtMax_HWW120= new TH1D("PtMax_HWW120", ";Leading Lepton p_{T} [GeV/c];Fraction of Events", 100, 0, 100);
  TH1D *DileptonMass = new TH1D("DileptonMass", ";M_{ll} [GeV/c^{2}];Fraction of Events", 100, 0, 100);


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *genInfo    = new mithep::TGenInfo();

  
  //********************************************************
  // mH 160
  //********************************************************
  eventTree = getTreeFromFile("/home/sixie/hist/HwwAcceptance/MCAtNLOHerwig_ggHToWW_mH160.root","Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info", &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &genInfo);        TBranch *genInfoBr = eventTree->GetBranch("Gen");


  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);

    Double_t eventweight = info->eventweight;
 
    HiggsPt_HWW160->Fill(genInfo->ptBosonSystem);
    DeltaPhi_HWW160->Fill((180.0/TMath::Pi())*acos(cos(genInfo->phi_1 - genInfo->phi_2)));
    Met_HWW160->Fill(genInfo->met);

  } //end loop over data     


  //********************************************************
  // mH 120
  //********************************************************
  eventTree = getTreeFromFile("/home/sixie/hist/HwwAcceptance/MCAtNLOHerwig_ggHToWW_mH120.root","Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info", &info);          infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &genInfo);        genInfoBr = eventTree->GetBranch("Gen");


  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);

    Double_t eventweight = info->eventweight;
 
    PtMax_HWW120->Fill(genInfo->pt_1);
    PtMin_HWW120->Fill(genInfo->pt_2);

  } //end loop over data     



  //********************************************************
  // Normalize Hists
  //********************************************************
  NormalizeHiggsPt(HiggsPt_HWW160);
  NormalizeHist(DeltaPhi_HWW160);
  NormalizeHist(Met_HWW160);
  NormalizeHist(PtMax_HWW120);
  NormalizeHist(PtMin_HWW120);

  //********************************************************
  // Plots
  //********************************************************
  TFile *file = new TFile("ThesisPlots.SignalCharacteristics.root","UPDATE");
  file->WriteTObject(HiggsPt_HWW160, HiggsPt_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(DeltaPhi_HWW160, DeltaPhi_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(Met_HWW160, Met_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(PtMax_HWW120, PtMax_HWW120->GetName(), "WriteDelete");
  file->WriteTObject(PtMin_HWW120, PtMin_HWW120->GetName(), "WriteDelete");
  file->Close();
   
  //********************************************************
  // Draw
  //********************************************************
  TCanvas *cv = 0;
  
  cv = new TCanvas("cv", "cv", 800, 600);
//   cv->SetLogy();
  HiggsPt_HWW160->GetYaxis()->SetTitleOffset(1.4);
  HiggsPt_HWW160->GetXaxis()->SetTitleOffset(1.05);
  HiggsPt_HWW160->Draw("hist");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_HiggsPt.gif");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_HiggsPt.eps");

  cv = new TCanvas("cv", "cv", 800, 600);
  Met_HWW160->GetYaxis()->SetTitleOffset(1.4);
  Met_HWW160->GetXaxis()->SetTitleOffset(1.05);
  Met_HWW160->Draw("hist");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_Met.gif");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_Met.eps");

  cv = new TCanvas("cv", "cv", 800, 600);
  DeltaPhi_HWW160->GetYaxis()->SetTitleOffset(1.4);
  DeltaPhi_HWW160->GetXaxis()->SetTitleOffset(1.05);
  DeltaPhi_HWW160->Draw("hist");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_DeltaPhi.gif");
  cv->SaveAs("Selection_SignalCharacteristics_HWW160_DeltaPhi.eps");


  cv = new TCanvas("cv", "cv", 800, 600);
  PtMin_HWW120->GetYaxis()->SetTitleOffset(1.4);
  PtMin_HWW120->GetXaxis()->SetTitleOffset(1.05);
  PtMin_HWW120->Draw("hist");
  cv->SaveAs("Selection_SignalCharacteristics_HWW120_PtMin.gif");
  cv->SaveAs("Selection_SignalCharacteristics_HWW120_PtMin.eps");

  cv = new TCanvas("cv", "cv", 800, 600);
  PtMax_HWW120->GetYaxis()->SetTitleOffset(1.4);
  PtMax_HWW120->GetXaxis()->SetTitleOffset(1.05);
  PtMax_HWW120->Draw("hist");
  cv->SaveAs("Selection_SignalCharacteristics_HWW120_PtMax.gif");
  cv->SaveAs("Selection_SignalCharacteristics_HWW120_PtMax.eps");



}



