//root -l EWKAna/Hww/Acceptance/HiggsPtAnalysis.C+\(\)
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
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TLegend.h"

#endif

#define LUMINOSITY 2.88 //(in pb^-1)
#define NBINSPASS 60
#define NBINSFAIL 24

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
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
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


void PlotKFactor() {

  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  TCanvas *cv = 0;
  cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.23,0.55,0.43,0.70);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);


  TH1D* PowhegKFactor160 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww160_PowhegToNNLL");
  TH1D* MCAtNLOKFactor160 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww160_MCAtNLOToNNLL");
  tmpLegend->Clear();
  tmpLegend->AddEntry(PowhegKFactor160, "Powheg", "LP");   
  tmpLegend->AddEntry(MCAtNLOKFactor160, "MC@NLO", "LP");   


  MCAtNLOKFactor160->SetLineColor(kBlack);
  MCAtNLOKFactor160->Draw("hist");
  MCAtNLOKFactor160->GetYaxis()->SetTitleOffset(1.2);
  MCAtNLOKFactor160->GetYaxis()->SetTitle("KFactor");
  MCAtNLOKFactor160->GetXaxis()->SetTitleOffset(1.05);
  MCAtNLOKFactor160->GetYaxis()->SetRangeUser(0.3,3.2);

  PowhegKFactor160->SetLineColor(kRed);
  PowhegKFactor160->Draw("hist,same");
  tmpLegend->Draw();
  cv->SaveAs("KFactorPowhegVsMCAtNLO.gif");




//   TH1D* KFactor160 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww160_PowhegToNNLL");
//   TH1D* KFactor200 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww200_PowhegToNNLL");
//   TH1D* KFactor250 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww250_PowhegToNNLL");

//   TCanvas *cv = new TCanvas("cv","cv", 800,600);

//   TLegend *tmpLegend = new TLegend(0.73,0.55,0.93,0.70);   
//   tmpLegend->SetTextSize(0.03);
//   tmpLegend->SetBorderSize(1);
//   tmpLegend->AddEntry(KFactor160, "m_{H} = 160", "LP");   
//   tmpLegend->AddEntry(KFactor200, "m_{H} = 200", "LP");   
//   tmpLegend->AddEntry(KFactor250, "m_{H} = 250", "LP");   

//   KFactor160->SetLineColor(kBlack);
//   KFactor160->Draw("hist");
//   KFactor160->GetYaxis()->SetTitleOffset(1.2);
//   KFactor160->GetYaxis()->SetTitle("KFactor");
//   KFactor160->GetXaxis()->SetTitleOffset(1.05);
//   KFactor200->SetLineColor(kRed);
//   KFactor200->Draw("hist,same");
//   KFactor250->SetLineColor(kBlue);
//   KFactor250->Draw("hist,same");
//   tmpLegend->Draw();
//   cv->SaveAs("KFactorVsMass.gif");


 
}

void MakeHQTHiggsPtDistribution( string inputFilename ,  const string Label) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
  string label = "";
  if (Label != "") label = "_" + Label;
 


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  ifstream infile(inputFilename.c_str());
  vector<double> lowEdges;
  vector<double> binContent;
  double x_min;
  double x_max;
  double y;
  Int_t binNumber = 1;
  while(infile >> x_min) {
    infile >> x_max >> y;
    lowEdges.push_back(x_min);
    binContent.push_back(y);
  }

  Double_t *xLowEdges = new Double_t[lowEdges.size()];
  for (UInt_t i=0; i < lowEdges.size();i++) {
    xLowEdges[i] = lowEdges[i];
  }


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *bosonSystemPt = new TH1D((string("HiggsPt")+ label).c_str(), ";  Higgs p_{T} [GeV/c]; Number of Events", lowEdges.size()-1, xLowEdges);

  for (UInt_t i=1;i<lowEdges.size();i++) {
    bosonSystemPt->SetBinContent(i,binContent[i]);
  }

  TFile *file = new TFile("HwwAcceptanceSystematics.root", "UPDATE");
  file->WriteTObject(bosonSystemPt, bosonSystemPt->GetName(), "WriteDelete");
  file->Close();

}


void MakeMCHiggsPtDistribution( string inputFilename ,  const string Label) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
  string label = "";
  if (Label != "") label = "_" + Label;
 

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *genMet = new TH1D("genMet", ";  Met [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonPtMin = new TH1D("leptonPtMin", ";  Lepton p_{T} [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonEta = new TH1D("leptonEta", "; Lepton #eta; Number of Events", 50, -5, 5);
  TH1D *bosonSystemPt = new TH1D((string("bosonSystemPt")+ label).c_str(), ";  Higgs p_{T} [GeV/c]; Number of Events", 80, 1.25, 201.25);
  TH1D *bosonSystemPtFineBinning = new TH1D((string("bosonSystemPtFineBinning")+ label).c_str(), ";  Higgs p_{T} [GeV/c]; Number of Events", 4000, 1.225, 201.225);
  TH1D *leadingGenJetPt = new TH1D((string("leadingGenJetPt")+ label).c_str() , "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *leadingGenJetPt_reweighted = new TH1D((string("leadingGenJetPt_reweighted")+ label).c_str() , "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *nGenJets = new TH1D("nGenJets", "; Mass [GeV/c^{2}]; Number of Events", 10, -0.5, 9.5);
  TH1D *leptonDeltaPhi = new TH1D("leptonDeltaPhi", "; #Delta #phi (l,l); Number of Events", 50, 0, 180);


  Double_t NEventsGenerated = 0;
  Double_t NEventsAccepted = 0;

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
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &genInfo); TBranch *genInfoBr = eventTree->GetBranch("Gen");

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);
		
    NEventsGenerated++;

    mithep::FourVectorM lepton1;
    mithep::FourVectorM lepton2;
    lepton1.SetCoordinates(genInfo->pt_1, genInfo->eta_1, genInfo->phi_1, 0.51099892e-3 );
    lepton2.SetCoordinates(genInfo->pt_2, genInfo->eta_2, genInfo->phi_2, 0.51099892e-3 );
    mithep::FourVectorM dilepton = lepton1+lepton2;

    if (info->eventweight >= 0 && genInfo->ptBosonSystem <= 300) {
      bosonSystemPt->Fill(genInfo->ptBosonSystem);
      bosonSystemPtFineBinning->Fill(genInfo->ptBosonSystem);
    } else {
      bosonSystemPt->Fill(genInfo->ptBosonSystem, -1);
      bosonSystemPtFineBinning->Fill(genInfo->ptBosonSystem,-1);
    }


  } //end loop over data     

  delete info;
  delete genInfo;

  //--------------------------------------------------------------------------------------------------------------
  // Normalize Higgs Pt spectrum
  //==============================================================================================================
  Double_t norm = 0 ;
  for(int i=0; i<bosonSystemPt->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += bosonSystemPt->GetBinContent(i);
  }
  for(int i=0; i<bosonSystemPt->GetXaxis()->GetNbins()+2 ; ++i) {
    bosonSystemPt->SetBinContent(i, bosonSystemPt->GetBinContent(i) / norm);
  }
//   for(int i=0; i<bosonSystemPt->GetXaxis()->GetNbins()+2 ; ++i) {
//     cout << "BosonPt " << i << " " << bosonSystemPt->GetBinContent(i) << endl;
//   }

  //smoothing - smooth to 1GeV bins
  double total = 0;
  for(int i=0; i<25 ; ++i) {
    total = 0;
    for (int j=0; j<20;j++) {
      total += bosonSystemPtFineBinning->GetBinContent(i*20 + j);
    }
    for (int j=0; j<20;j++) {
      bosonSystemPtFineBinning->SetBinContent(i*20 + j, total / 20);
    }
  }
  //smoothing - first big bin
  total = 0;
  for(int i=477; i<502 ; ++i) {
    total += bosonSystemPtFineBinning->GetBinContent(i);
  }
  for(int i=477; i<502 ; ++i) {
    bosonSystemPtFineBinning->SetBinContent(i,total / 26);
  }
  //smoothing - remaining big bins
  for (int i=10; i<80; ++i) {
    total = 0;
    for (int j=0; j<50;j++) {
      if (i*50 + 2 + j < 4001) {
        total += bosonSystemPtFineBinning->GetBinContent(i*50 + 2 + j);
      }
    }
    for (int j=0; j<50;j++) {
      if (i < 79) {        
        bosonSystemPtFineBinning->SetBinContent(i*50 + 2 + j, total / 50);
      } else {
        bosonSystemPtFineBinning->SetBinContent(i*50 + 2 + j, total / 49);
      }
    }
  }
  //normalize;
  norm = 0 ;
  for(int i=0; i<bosonSystemPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += bosonSystemPtFineBinning->GetBinContent(i);
  }
  for(int i=0; i<bosonSystemPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    bosonSystemPtFineBinning->SetBinContent(i, bosonSystemPtFineBinning->GetBinContent(i) / norm);
  }

  TFile *file = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");
  file->WriteTObject(bosonSystemPt, bosonSystemPt->GetName(), "WriteDelete");
  file->WriteTObject(bosonSystemPtFineBinning, bosonSystemPtFineBinning->GetName(), "WriteDelete");
//   file->Close();

}


void ComputeHiggsPtKFactor( ) {
  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  TH1D* MCAtNLOHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_MCAtNLO_default");
  TH1D* PowhegHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_default");
  PowhegHiggsPt->SetLineColor(kRed);
  TH1D* NNLLHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_NNLL_default");
  NNLLHiggsPt->SetLineColor(kBlue);

  TH1D* MCAtNLOHiggsPtFineBinning = (TH1D*)EffFile->Get("bosonSystemPtFineBinning_ggHww160_MCAtNLO_default");
  TH1D* PowhegHiggsPtFineBinning = (TH1D*)EffFile->Get("bosonSystemPtFineBinning_ggHww160_default");
  TH1D* NNLLHiggsPtFineBinning = (TH1D*)EffFile->Get("bosonSystemPtFineBinning_ggHww160_NNLL_default");

  assert(MCAtNLOHiggsPt);
  assert(PowhegHiggsPt);
  assert(NNLLHiggsPt);
  assert(MCAtNLOHiggsPtFineBinning);
  assert(PowhegHiggsPtFineBinning);
  assert(NNLLHiggsPtFineBinning);


  Double_t norm = 0 ;

  for(int i=0; i<PowhegHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += PowhegHiggsPt->GetBinContent(i);
  }
  for(int i=0; i<PowhegHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    PowhegHiggsPt->SetBinContent(i, PowhegHiggsPt->GetBinContent(i) / norm);
//     cout << i << " " << PowhegHiggsPt->GetBinContent(i) << endl;
  }

  norm = 0;
  for(int i=0; i<MCAtNLOHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += MCAtNLOHiggsPt->GetBinContent(i);
  }
  for(int i=0; i<MCAtNLOHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    MCAtNLOHiggsPt->SetBinContent(i, MCAtNLOHiggsPt->GetBinContent(i) / norm);
    // cout << i << " " << MCAtNLOHiggsPt->GetBinContent(i) << endl;
  }
  
  norm = 0 ;
  for(int i=0; i<PowhegHiggsPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += PowhegHiggsPtFineBinning->GetBinContent(i);
  }
  for(int i=0; i<PowhegHiggsPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    PowhegHiggsPtFineBinning->SetBinContent(i, PowhegHiggsPtFineBinning->GetBinContent(i) / norm);
  }

  norm = 0;
  for(int i=0; i<MCAtNLOHiggsPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += MCAtNLOHiggsPtFineBinning->GetBinContent(i);
  }
  for(int i=0; i<MCAtNLOHiggsPtFineBinning->GetXaxis()->GetNbins()+2 ; ++i) {
    MCAtNLOHiggsPtFineBinning->SetBinContent(i, MCAtNLOHiggsPtFineBinning->GetBinContent(i) / norm);
    // cout << i << " " << MCAtNLOHiggsPtFineBinning->GetBinContent(i) << endl;
  }


  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.73,0.55,0.93,0.70);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(MCAtNLOHiggsPt, "MC@NLO", "LP");   
  tmpLegend->AddEntry(PowhegHiggsPt, "Powheg", "LP");   
  tmpLegend->AddEntry(NNLLHiggsPt, "NNLO+NNLL", "LP");   

  MCAtNLOHiggsPt->SetMaximum(0.07);
  MCAtNLOHiggsPt->Draw("hist");
  MCAtNLOHiggsPt->GetYaxis()->SetTitleOffset(1.2);
  MCAtNLOHiggsPt->GetYaxis()->SetTitle("Fraction of Events");
  MCAtNLOHiggsPt->GetXaxis()->SetTitleOffset(1.05);
  PowhegHiggsPt->Draw("hist,same");
  NNLLHiggsPt->Draw("hist,same");
  tmpLegend->Draw();
  cv->SaveAs("HiggsPtComparison.gif");

  cv->SetLogy();
  cv->SaveAs("HiggsPtComparison_logY.gif");




  TH1D* kFactorHiggsPt_PowhegToNNLL = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_PowhegToNNLL");
  for(int i=0; i<kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    if (i==0) {
      kFactorHiggsPt_PowhegToNNLL->SetBinContent(i,1.0);
    } else if (i == kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+1) {
      kFactorHiggsPt_PowhegToNNLL->SetBinContent(i,0.6);
    } else {
      kFactorHiggsPt_PowhegToNNLL->SetBinContent(i,NNLLHiggsPt->GetBinContent(i) / PowhegHiggsPt->GetBinContent(i));    
    }
//     cout << i << " " << NNLLHiggsPt->GetBinContent(i) << " / " << PowhegHiggsPt->GetBinContent(i) << " " << kFactorHiggsPt_PowhegToNNLL->GetBinContent(i) << endl;
  }
  TH1D* kFactorHiggsPt_MCAtNLOToNNLL = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_MCAtNLOToNNLL");
  for(int i=0; i<kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    if (i==0) {
      kFactorHiggsPt_MCAtNLOToNNLL->SetBinContent(i,1.0);
    } else if (i == kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+1) {
      kFactorHiggsPt_MCAtNLOToNNLL->SetBinContent(i,1.5);
    } else {
      kFactorHiggsPt_MCAtNLOToNNLL->SetBinContent(i,NNLLHiggsPt->GetBinContent(i) / MCAtNLOHiggsPt->GetBinContent(i));
    }
  }

  TH1D* kFactorHiggsPtFineBinning_PowhegToNNLL = (TH1D*)MCAtNLOHiggsPtFineBinning->Clone("kFactorHiggsPtFineBinning_ggHww160_PowhegToNNLL");
  for(int i=0; i<kFactorHiggsPtFineBinning_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    if (i==0) {
    kFactorHiggsPtFineBinning_PowhegToNNLL->SetBinContent(i,1.0);
    } else if (i == kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+1) {
    kFactorHiggsPtFineBinning_PowhegToNNLL->SetBinContent(i,0.6);
    } else {
    kFactorHiggsPtFineBinning_PowhegToNNLL->SetBinContent(i,NNLLHiggsPtFineBinning->GetBinContent(i) / PowhegHiggsPtFineBinning->GetBinContent(i));
    }
//     cout << i << " " << NNLLHiggsPtFineBinning->GetBinContent(i) << " / " << PowhegHiggsPtFineBinning->GetBinContent(i) << " " << kFactorHiggsPtFineBinning_PowhegToNNLL->GetBinContent(i) << endl;
  }
  TH1D* kFactorHiggsPtFineBinning_MCAtNLOToNNLL = (TH1D*)MCAtNLOHiggsPtFineBinning->Clone("kFactorHiggsPtFineBinning_ggHww160_MCAtNLOToNNLL");
  for(int i=0; i<kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    if (i==0) {
      kFactorHiggsPtFineBinning_MCAtNLOToNNLL->SetBinContent(i,1.0);
    } else if (i == kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+1) {
      kFactorHiggsPtFineBinning_MCAtNLOToNNLL->SetBinContent(i,1.5);
    } else {
      kFactorHiggsPtFineBinning_MCAtNLOToNNLL->SetBinContent(i,NNLLHiggsPtFineBinning->GetBinContent(i) / MCAtNLOHiggsPtFineBinning->GetBinContent(i));
    }
//    cout << i << " " << NNLLHiggsPtFineBinning->GetBinContent(i) << " / " << MCAtNLOHiggsPtFineBinning->GetBinContent(i) << " " << kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetBinContent(i) << endl;
  }

  cv->SetLogy(0);
  tmpLegend = new TLegend(0.20,0.75,0.55,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(kFactorHiggsPt_PowhegToNNLL, "NNLO+NNLL/Powheg", "LP");   
   tmpLegend->AddEntry(kFactorHiggsPt_MCAtNLOToNNLL, "NNLO+NNLL/MC@NLO", "LP");   
//   kFactorHiggsPt_MCAtNLOToNNLL->GetYaxis()->SetTitle("Ratio");
//    kFactorHiggsPt_MCAtNLOToNNLL->Draw("hist");
//   kFactorHiggsPt_MCAtNLOToNNLL->SetMinimum(0.3);
//   kFactorHiggsPt_MCAtNLOToNNLL->SetMaximum(3.2);
  kFactorHiggsPt_PowhegToNNLL->SetLineColor(kRed);
  kFactorHiggsPt_PowhegToNNLL->Draw("hist");
  kFactorHiggsPt_PowhegToNNLL->SetMinimum(0.0);
  kFactorHiggsPt_PowhegToNNLL->SetMaximum(3.2);
  kFactorHiggsPt_MCAtNLOToNNLL->Draw("hist,same");  
  tmpLegend->Draw(); 
//   cv->SetLogy(0);
  cv->SaveAs("HiggsPtKFactors.gif");



  //*****************************************************************************************
  //Print KFactors
  //*****************************************************************************************
//   //To Produce Tex
//   for(int i=0; i<kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
//     char kfactor[10];
//     sprintf(kfactor,"%.4f ",kFactorHiggsPt_PowhegToNNLL->GetBinContent(i));

//     cout << kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetBinLowEdge(i) << " - " << kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetBinUpEdge(i) << " &  " << kfactor <<  " \\\\" << endl;
//   }

//   for(int i=0; i<kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
//     char kfactor[10];
//     sprintf(kfactor,"%.4f ",kFactorHiggsPt_MCAtNLOToNNLL->GetBinContent(i));

//     cout << kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetBinLowEdge(i) << " - " << kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetBinUpEdge(i) << " &  " << kfactor <<  " \\\\" << endl;
//   }


//   //To Produce twiki table
//   for(int i=0; i<kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
//     char kfactor[10];
//     sprintf(kfactor,"%.4f ",kFactorHiggsPt_PowhegToNNLL->GetBinContent(i));

//     cout << "| " << kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetBinLowEdge(i) << " - " << kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetBinUpEdge(i) << " | " << kfactor <<  " | " << endl;
//   }

//   for(int i=0; i<kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
//     char kfactor[10];
//     sprintf(kfactor,"%.4f ",kFactorHiggsPt_MCAtNLOToNNLL->GetBinContent(i));

//     cout << "| " << kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetBinLowEdge(i) << " - " << kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetBinUpEdge(i) << " | " << kfactor <<  " |" << endl;
//   }


  //To produce the KFactor file
  for(int i=0; i<kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    char kfactor[10];
    sprintf(kfactor,"%.4f ",kFactorHiggsPt_PowhegToNNLL->GetBinContent(i));    
    cout << i << " " << kfactor <<  endl;
  }

  for(int i=0; i<kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    char kfactor[10];
    sprintf(kfactor,"%.4f ",kFactorHiggsPt_MCAtNLOToNNLL->GetBinContent(i));
    cout << i << " " << kfactor <<  endl;
  }




  EffFile->WriteTObject(kFactorHiggsPt_PowhegToNNLL, kFactorHiggsPt_PowhegToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPt_MCAtNLOToNNLL, kFactorHiggsPt_MCAtNLOToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPtFineBinning_PowhegToNNLL, kFactorHiggsPtFineBinning_PowhegToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPtFineBinning_MCAtNLOToNNLL, kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetName(), "WriteDelete");
  EffFile->Close();

}


void HiggsPtAnalysis() {


  //*****************************************************************************************
  //Produce Higgs Pt Distribution From HQT files
  //*****************************************************************************************
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh115.txt","ggHww_mH115_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh120.txt","ggHww_mH120_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh130.txt","ggHww_mH130_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh140.txt","ggHww_mH140_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh150.txt","ggHww_mH150_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh160.txt","ggHww_mH160_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh170.txt","ggHww_mH170_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh180.txt","ggHww_mH180_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh190.txt","ggHww_mH190_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh200.txt","ggHww_mH200_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh210.txt","ggHww_mH210_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh220.txt","ggHww_mH220_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh230.txt","ggHww_mH230_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh250.txt","ggHww_mH250_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh300.txt","ggHww_mH300_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh350.txt","ggHww_mH350_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh400.txt","ggHww_mH400_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh450.txt","ggHww_mH450_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh500.txt","ggHww_mH500_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh550.txt","ggHww_mH550_default_HQT");
  MakeHQTHiggsPtDistribution("MitPhysics/data/HiggsPt/HqT-mT172p5.mh600.txt","ggHww_mH600_default_HQT");


  //*****************************************************************************************
  //Produce Higgs Pt Distribution for MC's
  //*****************************************************************************************
//   MakeMCHiggsPtDistribution("/home/sixie/hist/HwwAcceptance/ggHWW160GenOnly_default.root","ggHww160_default");
//   MakeMCHiggsPtDistribution("/home/sixie/hist/HwwAcceptance/ggHWW160GenOnly_MCAtNLO_default.root","ggHww160_MCAtNLO_default");


  //*****************************************************************************************
  //Compute KFactor
  //*****************************************************************************************
//   ComputeHiggsPtKFactor();

  //*****************************************************************************************
  //Plot & Print KFactors
  //*****************************************************************************************
//   PlotKFactor();


}
