//root -l EWKAna/Hww/Acceptance/MakeResbosHiggsPt.C+\(\)
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
#include "TTree.h"

#endif

void MakeResbosHiggsPt( ) {
  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  TFile *inputFile = new TFile("/home/sixie/hist/HwwAcceptance/Resbos/Resbos_H_NNLO_mh160.root", "UPDATE");
  TTree *tree = (TTree*)inputFile->Get("h10");
  TH1F* ResbosHiggsPt = new TH1F("bosonSystemPt_ggHww160_Resbos", "", 80, 0, 200);

  tree->Draw("pT_H>>bosonSystemPt_ggHww160_Resbos","pT_H<=1000");
 
  Double_t norm = 0 ;
  for(int i=0; i<ResbosHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    norm += ResbosHiggsPt->GetBinContent(i);
  }
  for(int i=0; i<ResbosHiggsPt->GetXaxis()->GetNbins()+2 ; ++i) {
    ResbosHiggsPt->SetBinContent(i, ResbosHiggsPt->GetBinContent(i) / norm);
    cout << i << " " << ResbosHiggsPt->GetBinContent(i) << endl;
  }

  EffFile->WriteTObject(ResbosHiggsPt, ResbosHiggsPt->GetName(), "WriteDelete");
  EffFile->Close();





}


