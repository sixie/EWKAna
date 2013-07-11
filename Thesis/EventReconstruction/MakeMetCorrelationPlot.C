#include "SmurfTree.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TStyle.h"
#include "MitStyle.h"

//------------------------------------------------------------------------------
// WW control region macro
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void MakeMetCorrelationPlot ()
{


  TH2F *signalMet = new TH2F("signalMet",";Projected MET [GeV/c] ;Projected TrackMET [GeV/c];", 50,0,100,50,0,100);
  TH2F *bkgMet = new TH2F("bkgMet",";Projected MET [GeV/c] ;Projected TrackMET [GeV/c];", 50,0,100,50,0,100);

// ***********************************************************************************************
// Full Run2011A Data 
// ***********************************************************************************************
  TChain *chsignal = new TChain("tree");
  chsignal->Add("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/hww130.root");
  TTree *signal = (TTree*) chsignal;

  TChain *chbackground = new TChain("tree");
  chbackground->Add("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/dymm.root");
  TTree *background = (TTree*) chbackground;


  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          npu;
  UInt_t          njets;
  UInt_t          run;
  UInt_t          event;
  Float_t         scale1fb;
  UInt_t          type;
  Float_t         met;
  Float_t         trackMet;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           processId;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;
  Float_t         jet2Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Float_t         higgsPt = -999;

  //***********************************************************************************************
  //Signal
  //***********************************************************************************************

  signal->SetBranchAddress( "pmet"         , &pmet         );
  signal->SetBranchAddress( "pTrackMet"    , &pTrackMet    );

  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    signal->GetEntry(i);
    signalMet->Fill(pmet,pTrackMet);
  }

//   //***********************************************************************************************
//   //Background
//   //***********************************************************************************************
//   background->SetBranchAddress( "pmet"         , &pmet         );
//   background->SetBranchAddress( "pTrackMet"    , &pTrackMet    );

//   cout << "Total Bkg: " << background->GetEntries() << endl;
//   for (UInt_t i=0; i<background->GetEntries(); i++) {
//     if (i%100000==0) cout << "Event " << i << endl;
//     background->GetEntry(i);
//     bkgMet->Fill(pmet,pTrackMet);
//   }
  
  TCanvas *cv = 0;
  cv = new TCanvas("cv","cv",800,600);
  signalMet->SetTitle("");
  signalMet->GetXaxis()->SetTitle("Projected MET [GeV/c]");
  signalMet->GetYaxis()->SetTitle("Projected TrackMET [GeV/c]");
  signalMet->GetYaxis()->SetTitleOffset(1.3);
  signalMet->GetXaxis()->SetTitleOffset(1.05);
  signalMet->Draw("colz");
  cv->SaveAs("MetCorrelations_HWW130.eps");

//   cv = new TCanvas("cv","cv",800,600);
//   bkgMet->SetTitle("");
//   bkgMet->GetXaxis()->SetTitle("Projected MET [GeV/c]");
//   bkgMet->GetYaxis()->SetTitle("Projected TrackMET [GeV/c]");
//   bkgMet->GetYaxis()->SetTitleOffset(1.3);
//   bkgMet->GetXaxis()->SetTitleOffset(1.05);
//   bkgMet->Draw("colz");
//   cv->SaveAs("MetCorrelations_Bkg.eps");




  return;

}


