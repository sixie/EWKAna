//root -l EWKAna/Hww/Acceptance/CompareHiggsPt.C+\(\)
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

void PlotKFactor() {

  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  TH1D* KFactor160 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww160_PowhegToNNLL");
  TH1D* KFactor200 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww200_PowhegToNNLL");
  TH1D* KFactor250 = (TH1D*)EffFile->Get("kFactorHiggsPt_ggHww250_PowhegToNNLL");

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.73,0.55,0.93,0.70);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(KFactor160, "m_{H} = 160", "LP");   
  tmpLegend->AddEntry(KFactor200, "m_{H} = 200", "LP");   
  tmpLegend->AddEntry(KFactor250, "m_{H} = 250", "LP");   

  KFactor160->SetLineColor(kBlack);
  KFactor160->Draw("hist");
  KFactor160->GetYaxis()->SetTitleOffset(1.2);
  KFactor160->GetYaxis()->SetTitle("KFactor");
  KFactor160->GetXaxis()->SetTitleOffset(1.05);
  KFactor200->SetLineColor(kRed);
  KFactor200->Draw("hist,same");
  KFactor250->SetLineColor(kBlue);
  KFactor250->Draw("hist,same");
  tmpLegend->Draw();
  cv->SaveAs("KFactorVsMass.gif");


 
}

void CompareHiggsPt( ) {
  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  TH1D* MCAtNLOHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_MCAtNLO_default");
  TH1D* PowhegHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_default");
  PowhegHiggsPt->SetLineColor(kRed);
//   TH1D* ResbosHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_Resbos");
//   ResbosHiggsPt->SetLineColor(kBlue);
  TH1D* NNLLHiggsPt = (TH1D*)EffFile->Get("bosonSystemPt_ggHww160_NNLL_default");
  NNLLHiggsPt->SetLineColor(kGreen);

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

//   norm = 0;
//   for(int i=0; i<ResbosHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
//     norm += ResbosHiggsPt->GetBinContent(i);
//   }
//   for(int i=0; i<ResbosHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
//     ResbosHiggsPt->SetBinContent(i, ResbosHiggsPt->GetBinContent(i) / norm);
//     cout << i << " " << ResbosHiggsPt->GetBinContent(i) << endl;
//   }

  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.73,0.55,0.93,0.70);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(MCAtNLOHiggsPt, "MC@NLO", "LP");   
  tmpLegend->AddEntry(PowhegHiggsPt, "Powheg", "LP");   
//   tmpLegend->AddEntry(ResbosHiggsPt, "Resbos", "LP");   
  tmpLegend->AddEntry(NNLLHiggsPt, "NNLO+NNLL", "LP");   

  MCAtNLOHiggsPt->SetMaximum(0.07);
  MCAtNLOHiggsPt->Draw("hist");
  MCAtNLOHiggsPt->GetYaxis()->SetTitleOffset(1.2);
  MCAtNLOHiggsPt->GetYaxis()->SetTitle("Fraction of Events");
  MCAtNLOHiggsPt->GetXaxis()->SetTitleOffset(1.05);
  PowhegHiggsPt->Draw("hist,same");
//   ResbosHiggsPt->Draw("hist,same");
  NNLLHiggsPt->Draw("hist,same");
  tmpLegend->Draw();
  cv->SaveAs("HiggsPtComparison.gif");

  cv->SetLogy();
  cv->SaveAs("HiggsPtComparison_logY.gif");


//   TH1D* leadingGenJet_MCAtNLO = (TH1D*)EffFile->Get("leadingGenJetPt_ggHww160_MCAtNLO_default");
//   norm = 0;
//   for(int i=0; i<leadingGenJet_MCAtNLO->GetXaxis()->GetNbins()+2 ; ++i) {
//     norm += leadingGenJet_MCAtNLO->GetBinContent(i);
//   }
//   for(int i=0; i<leadingGenJet_MCAtNLO->GetXaxis()->GetNbins()+2 ; ++i) {
//     leadingGenJet_MCAtNLO->SetBinContent(i, leadingGenJet_MCAtNLO->GetBinContent(i) / norm);
//   }
//   leadingGenJet_MCAtNLO->SetLineColor(kBlack);
//   leadingGenJet_MCAtNLO->GetYaxis()->SetTitleOffset(1.2);
//   TH1D* leadingGenJet_Powheg = (TH1D*)EffFile->Get("leadingGenJetPt_ggHww160_default");
//   norm = 0;
//   for(int i=0; i<leadingGenJet_Powheg->GetXaxis()->GetNbins()+2 ; ++i) {
//     norm += leadingGenJet_Powheg->GetBinContent(i);
//   }
//   for(int i=0; i<leadingGenJet_Powheg->GetXaxis()->GetNbins()+2 ; ++i) {
//     leadingGenJet_Powheg->SetBinContent(i, leadingGenJet_Powheg->GetBinContent(i) / norm);
//   }

//   tmpLegend->Clear();
//   tmpLegend->AddEntry(leadingGenJet_MCAtNLO, "MC@NLO", "LP");    
//   tmpLegend->AddEntry(leadingGenJet_Powheg, "Powheg", "LP");   
  
//   leadingGenJet_MCAtNLO->Draw("hist");
//   leadingGenJet_Powheg->Draw("hist,same");
//   tmpLegend->Draw();
//   cv->SaveAs("LeadingGenJetPtComparison_logY.gif");



  TH1D* kFactorHiggsPt_PowhegToMCAtNLO = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_PowhegToMCAtNLO");
  for(int i=0; i<kFactorHiggsPt_PowhegToMCAtNLO->GetXaxis()->GetNbins()+2 ; ++i) {
    kFactorHiggsPt_PowhegToMCAtNLO->SetBinContent(i,MCAtNLOHiggsPt->GetBinContent(i) / PowhegHiggsPt->GetBinContent(i));
  }
//   TH1D* kFactorHiggsPt_PowhegToResbos = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_PowhegToResbos");
//   for(int i=0; i<kFactorHiggsPt_PowhegToResbos->GetXaxis()->GetNbins()+2 ; ++i) {
//     kFactorHiggsPt_PowhegToResbos->SetBinContent(i,ResbosHiggsPt->GetBinContent(i) / PowhegHiggsPt->GetBinContent(i));
//   }
//   TH1D* kFactorHiggsPt_MCAtNLOToResbos = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_MCAtNLOToResbos");
//   for(int i=0; i<kFactorHiggsPt_MCAtNLOToResbos->GetXaxis()->GetNbins()+2 ; ++i) {
//     kFactorHiggsPt_MCAtNLOToResbos->SetBinContent(i,ResbosHiggsPt->GetBinContent(i) / MCAtNLOHiggsPt->GetBinContent(i));
//   }
  TH1D* kFactorHiggsPt_PowhegToNNLL = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_PowhegToNNLL");
  for(int i=0; i<kFactorHiggsPt_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    kFactorHiggsPt_PowhegToNNLL->SetBinContent(i,NNLLHiggsPt->GetBinContent(i) / PowhegHiggsPt->GetBinContent(i));
  }
  TH1D* kFactorHiggsPt_MCAtNLOToNNLL = (TH1D*)MCAtNLOHiggsPt->Clone("kFactorHiggsPt_ggHww160_MCAtNLOToNNLL");
  for(int i=0; i<kFactorHiggsPt_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    kFactorHiggsPt_MCAtNLOToNNLL->SetBinContent(i,NNLLHiggsPt->GetBinContent(i) / MCAtNLOHiggsPt->GetBinContent(i));
  }
  TH1D* kFactorHiggsPtFineBinning_PowhegToNNLL = (TH1D*)MCAtNLOHiggsPtFineBinning->Clone("kFactorHiggsPtFineBinning_ggHww160_PowhegToNNLL");
  for(int i=0; i<kFactorHiggsPtFineBinning_PowhegToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    kFactorHiggsPtFineBinning_PowhegToNNLL->SetBinContent(i,NNLLHiggsPtFineBinning->GetBinContent(i) / PowhegHiggsPtFineBinning->GetBinContent(i));
    cout << i << " " << NNLLHiggsPtFineBinning->GetBinContent(i) << " / " << PowhegHiggsPtFineBinning->GetBinContent(i) << " " << kFactorHiggsPtFineBinning_PowhegToNNLL->GetBinContent(i) << endl;
  }
  TH1D* kFactorHiggsPtFineBinning_MCAtNLOToNNLL = (TH1D*)MCAtNLOHiggsPtFineBinning->Clone("kFactorHiggsPtFineBinning_ggHww160_MCAtNLOToNNLL");
  for(int i=0; i<kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetXaxis()->GetNbins()+2 ; ++i) {
    kFactorHiggsPtFineBinning_MCAtNLOToNNLL->SetBinContent(i,NNLLHiggsPtFineBinning->GetBinContent(i) / MCAtNLOHiggsPtFineBinning->GetBinContent(i));
//    cout << i << " " << NNLLHiggsPtFineBinning->GetBinContent(i) << " / " << MCAtNLOHiggsPtFineBinning->GetBinContent(i) << " " << kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetBinContent(i) << endl;
  }

  tmpLegend = new TLegend(0.20,0.75,0.55,0.90);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(kFactorHiggsPt_PowhegToNNLL, "NNLO+NNLL/Powheg", "LP");   
  tmpLegend->AddEntry(kFactorHiggsPt_MCAtNLOToNNLL, "NNLO+NNLL/MC@NLO", "LP");   
  kFactorHiggsPt_MCAtNLOToNNLL->GetYaxis()->SetTitle("Ratio");
  kFactorHiggsPt_MCAtNLOToNNLL->Draw("hist");
  kFactorHiggsPt_MCAtNLOToNNLL->SetMinimum(0.3);
  kFactorHiggsPt_MCAtNLOToNNLL->SetMaximum(3.2);
  kFactorHiggsPt_PowhegToNNLL->SetLineColor(kRed);
  kFactorHiggsPt_PowhegToNNLL->Draw("hist,same");
  tmpLegend->Draw();
  cv->SetLogy(0);
  cv->SaveAs("HiggsPtKFactors.gif");


//   EffFile->WriteTObject(kFactorHiggsPt_PowhegToMCAtNLO, kFactorHiggsPt_PowhegToMCAtNLO->GetName(), "WriteDelete");
//   EffFile->WriteTObject(kFactorHiggsPt_PowhegToResbos, kFactorHiggsPt_PowhegToResbos->GetName(), "WriteDelete");
//   EffFile->WriteTObject(kFactorHiggsPt_MCAtNLOToResbos, kFactorHiggsPt_MCAtNLOToResbos->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPt_PowhegToNNLL, kFactorHiggsPt_PowhegToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPt_MCAtNLOToNNLL, kFactorHiggsPt_MCAtNLOToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPtFineBinning_PowhegToNNLL, kFactorHiggsPtFineBinning_PowhegToNNLL->GetName(), "WriteDelete");
  EffFile->WriteTObject(kFactorHiggsPtFineBinning_MCAtNLOToNNLL, kFactorHiggsPtFineBinning_MCAtNLOToNNLL->GetName(), "WriteDelete");
  EffFile->Close();

  PlotKFactor();

}


