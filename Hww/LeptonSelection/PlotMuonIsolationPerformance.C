 
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
#include "MitCommon/MathTools/interface/MathUtils.h"
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

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
void AddIsoHistograms(TH1F* summedHist, vector<TH1F*> histVector, UInt_t StartIndex, UInt_t EndIndex) {

  assert (StartIndex < histVector.size());
  assert (EndIndex < histVector.size());

  cout << histVector.size() << endl;

  for (UInt_t b=0; b < UInt_t(summedHist->GetXaxis()->GetNbins()+2); ++b) {
    Double_t binContent = 0;
    Double_t binErrorSquared = 0;

    for( UInt_t i = StartIndex; i <= EndIndex; ++i) {
      binContent += histVector[i]->GetBinContent(b);
      binErrorSquared += pow(histVector[i]->GetBinError(b),2);
    }

    summedHist->SetBinContent(b,binContent);
    summedHist->SetBinError(b,TMath::Sqrt(binErrorSquared));
  }

  return;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double SigEff[nPoints];
  double BkgEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double BkgEffErrLow[nPoints];
  double BkgEffErrHigh[nPoints];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;
    for (UInt_t q=0; q <= b; ++q) {
      nsig += signalHist->GetBinContent(q);
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

    n1 = TMath::Nint(nbkg);
    n2 = TMath::Nint(NBkgTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    BkgEff[b] = ratio;
    BkgEffErrLow[b] = 0;
    BkgEffErrHigh[b] = 0;
//     BkgEffErrLow[b] = errLow;
//     BkgEffErrHigh[b] = errHigh;
  }

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh, SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);

  return tmpSigEffVsBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name, Double_t myCutValue ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double SigEff[1];
  double BkgEff[1];
  double SigEffErrLow[1];
  double SigEffErrHigh[1];
  double BkgEffErrLow[1];
  double BkgEffErrHigh[1];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  double cutValue;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }
  
  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;
    for (UInt_t q=0; q <= b; ++q) {
      nsig += signalHist->GetBinContent(q);
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      cout << "Cut: " << myCutValue << " : " << cutValue << " " << signalHist->GetXaxis()->GetBinCenter(b) << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue)) {
      cutValue = signalHist->GetXaxis()->GetBinCenter(b);
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
//     SigEffErrLow[0] = errLow;
//     SigEffErrHigh[0] = errHigh;

      n1 = TMath::Nint(nbkg);
      n2 = TMath::Nint(NBkgTotal);
      mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
      BkgEff[0] = ratio;
      BkgEffErrLow[0] = 0;
      BkgEffErrHigh[0] = 0;
//     BkgEffErrLow[0] = errLow;
//     BkgEffErrHigh[0] = errHigh;
    }
  }
    cout << "Final value : " << cutValue << " " << SigEff[0] << " " << BkgEff[0] << endl;
    TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (1, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh , SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerColor(kBlack);
  tmpSigEffVsBkgEff->SetLineColor(kBlack);
  tmpSigEffVsBkgEff->SetMarkerSize(1.5);

  return tmpSigEffVsBkgEff;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeSigEffVsCutValueGraph(TH1F* signalHist, string name ) {

  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[nPoints];
  double cutValueErr[nPoints];
  double SigEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double NSigTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    cutValue[b] = signalHist->GetXaxis()->GetBinCenter(b);
    cutValueErr[b] = 0;
    Double_t nsig = 0;
    for (UInt_t q=0; q <= b; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

  }

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (nPoints, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);

  return tmpSigEffVsCut;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, string name, Double_t myCutValue ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double cutValue[1] = {0};
  double cutValueErr[1] = {0};
  double SigEff[1] = {0};
  double SigEffErrLow[1] = {0};
  double SigEffErrHigh[1] = {0};
  double NSigTotal = 0;
  
  Double_t effDiff = 9999;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=0; q <= b; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    
      cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue[0])) {
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
//     SigEffErrLow[0] = errLow;
//     SigEffErrHigh[0] = errHigh;
      cutValue[0] = signalHist->GetXaxis()->GetBinCenter(b);
      cutValueErr[0] = 0;
    }
  }

   cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (1, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsCut->SetMarkerColor(kBlack);
  tmpSigEffVsCut->SetLineColor(kBlack);
  tmpSigEffVsCut->SetMarkerSize(1.5);

  return tmpSigEffVsCut;
}






//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeEffVsNVertex(vector<TH1F*> histVector, string name, Double_t cutValue, Int_t CutAboveOrBelow ) {
  //Make Met Plots
  const UInt_t nPoints = histVector.size();
  double NPileup[nPoints];
  double NPileupError[nPoints];
  double Eff[nPoints];
  double EffErrLow[nPoints];
  double EffErrHigh[nPoints];
  double NTotal = 0;
  

  for (UInt_t i=0; i<nPoints; ++i) {
    NTotal = 0;
    for (UInt_t q=0; q < UInt_t(histVector[i]->GetXaxis()->GetNbins()+2); ++q) {
      NTotal += histVector[i]->GetBinContent(q);
    }

    Double_t NPass = 0;
    if (CutAboveOrBelow == 1) {
      for (UInt_t q=0; q < UInt_t(histVector[i]->GetXaxis()->FindFixBin(cutValue)); ++q) {
        NPass += histVector[i]->GetBinContent(q);
      }
    } else if (CutAboveOrBelow == -1) {
      for (UInt_t q=histVector[i]->GetXaxis()->FindFixBin(0.15); q < UInt_t(histVector[i]->GetXaxis()->GetNbins()+2); ++q) {
        NPass += histVector[i]->GetBinContent(q);
      }
    }
      
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = TMath::Nint(NPass);
    Double_t n2 = TMath::Nint(NTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    Eff[i] = ratio;
    EffErrLow[i] = errLow;
    EffErrHigh[i] = errHigh;
    NPileup[i] = i;
//    cout << "MuonIsolationEff " << i << " : " << MuonIsolationEff[i] << endl;    
  }

  TGraphAsymmErrors *tmpEffVsNVertex = new TGraphAsymmErrors (nPoints,  NPileup, Eff, NPileupError, NPileupError, EffErrLow, EffErrHigh);
  tmpEffVsNVertex->SetName(name.c_str());
  tmpEffVsNVertex->SetTitle("");
  tmpEffVsNVertex->SetMarkerColor(kRed);
  tmpEffVsNVertex->GetXaxis()->SetTitleOffset(1.02);
  tmpEffVsNVertex->GetYaxis()->SetTitleOffset(1.05);

  return tmpEffVsNVertex;
}






void PlotMuonIsoVsNVertices(string Label = "" , string signalSample = "Zmm", string bkgSample = "WJetsMC") {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.root", "READ");
  
  vector<TH1F*>   SignalMuon_RelIso_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_Endcap;
  vector<TH1F*>   SignalMuon_RelIso_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;

  vector<TH1F*>   BkgMuon_RelIso_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_Endcap;
  vector<TH1F*>   BkgMuon_RelIso_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;



  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalMuon_RelIso_Barrel;
    TH1F *tmpSignalMuon_RelIso_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso_Endcap;
    TH1F *tmpSignalMuon_RelIso_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;

    TH1F *tmpBkgMuon_RelIso_Barrel;
    TH1F *tmpBkgMuon_RelIso_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso_Endcap;
    TH1F *tmpBkgMuon_RelIso_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;


    tmpSignalMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());

    tmpBkgMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());


    cout << "load " << n << endl;
    assert( tmpSignalMuon_RelIso_Barrel );
    assert( tmpSignalMuon_RelIso_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso_Endcap );
    assert( tmpSignalMuon_RelIso_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );

    assert( tmpBkgMuon_RelIso_Barrel );
    assert( tmpBkgMuon_RelIso_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso_Endcap );
    assert( tmpBkgMuon_RelIso_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );


    SignalMuon_RelIso_Barrel.push_back(tmpSignalMuon_RelIso_Barrel);
    SignalMuon_RelIso_RhoCorrected_Barrel.push_back(tmpSignalMuon_RelIso_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalMuon_RelIso_Endcap.push_back(tmpSignalMuon_RelIso_Endcap);
    SignalMuon_RelIso_RhoCorrected_Endcap.push_back(tmpSignalMuon_RelIso_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);

    BkgMuon_RelIso_Barrel.push_back(tmpBkgMuon_RelIso_Barrel);
    BkgMuon_RelIso_RhoCorrected_Barrel.push_back(tmpBkgMuon_RelIso_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgMuon_RelIso_Endcap.push_back(tmpBkgMuon_RelIso_Endcap);
    BkgMuon_RelIso_RhoCorrected_Endcap.push_back(tmpBkgMuon_RelIso_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);


  }

  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================

  vector<TGraphAsymmErrors*> EffVsNVertices;
  vector<string> GraphLabels;
  string plotname;

  //Make Met Plots
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
//   TPaveLabel *LabelText;
  TLegend * legend;
  TGraphAsymmErrors *tmpEffVsNVertex = 0;


  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  EffVsNVertices.clear();
  GraphLabels.clear();
  plotname = "MuonIsolationEfficiency_Barrel_Signal";
  //*************************************************************************************

  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_TotalPFRelIso_Barrel, "EffVsNVertex_TotalPFRelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_VertexSelectedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedPFRelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_TotalPFRelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF Corr");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedPFRelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU Corr");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(EffVsNVertices[i],GraphLabels[i].c_str(), "LP");

    EffVsNVertices[i]->SetMarkerColor(colors[i]);
    EffVsNVertices[i]->SetLineColor(colors[i]);
    EffVsNVertices[i]->SetMarkerSize(0.5);
   
    EffVsNVertices[i]->GetYaxis()->SetRangeUser(0.0,1.2);    
//    EffVsNVertices[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      EffVsNVertices[i]->Draw("AP");
    } else {
      EffVsNVertices[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("EffVsNVertices_" + plotname + label + ".gif").c_str());




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  EffVsNVertices.clear();
  GraphLabels.clear();
  plotname = "MuonIsolationEfficiency_Barrel_Bkg";
  //*************************************************************************************

  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_TotalPFRelIso_Barrel, "EffVsNVertex_TotalPFRelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_VertexSelectedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedPFRelIso_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_TotalPFRelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF Corr");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedPFRelIsoRhoCorrected_Barrel", 0.15, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU Corr");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(EffVsNVertices[i],GraphLabels[i].c_str(), "LP");

    EffVsNVertices[i]->SetMarkerColor(colors[i]);
    EffVsNVertices[i]->SetLineColor(colors[i]);
    EffVsNVertices[i]->SetMarkerSize(0.5);
   
    EffVsNVertices[i]->GetYaxis()->SetRangeUser(0.0,1.2);    
//    EffVsNVertices[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      EffVsNVertices[i]->Draw("AP");
    } else {
      EffVsNVertices[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("EffVsNVertices_" + plotname + label + ".gif").c_str());

}

















void MakeIsolationPerformancePlots(string Label = "", string signalSample = "Zmm", string bkgSample = "WJetsMC", Double_t currentWPCut = 0.15) {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  colors.push_back(kRed+3);
  colors.push_back(kGreen+3);
  colors.push_back(kAzure);
  
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.muons.root", "READ");
    
  vector<TH1F*>   SignalMuon_RelIso_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso05_Barrel;
  vector<TH1F*>   SignalMuon_RelIso05_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_Endcap;
  vector<TH1F*>   SignalMuon_RelIso_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_RelIso05_Endcap;
  vector<TH1F*>   SignalMuon_RelIso05_RhoCorrected_Endcap;

  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;


  vector<TH1F*>   BkgMuon_RelIso_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso05_Barrel;
  vector<TH1F*>   BkgMuon_RelIso05_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_Endcap;
  vector<TH1F*>   BkgMuon_RelIso_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_RelIso05_Endcap;
  vector<TH1F*>   BkgMuon_RelIso05_RhoCorrected_Endcap;

  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;






  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalMuon_RelIso_Barrel;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso05_Barrel;
    TH1F *tmpSignalMuon_RelIso05RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso_Endcap;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_RelIso05_Endcap;
    TH1F *tmpSignalMuon_RelIso05RhoCorrected_Endcap;

    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;



    TH1F *tmpBkgMuon_RelIso_Barrel;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso05_Barrel;
    TH1F *tmpBkgMuon_RelIso05RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso_Endcap;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_RelIso05_Endcap;
    TH1F *tmpBkgMuon_RelIso05RhoCorrected_Endcap;

    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;


    tmpSignalMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05_Barrel = (TH1F*)f->Get((string("Muon_RelIso05_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05_Endcap = (TH1F*)f->Get((string("Muon_RelIso05_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());

    tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());





    tmpBkgMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05_Barrel = (TH1F*)f->Get((string("Muon_RelIso05_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05_Endcap = (TH1F*)f->Get((string("Muon_RelIso05_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());

    tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());










    cout << "load " << n << endl;
    assert( tmpSignalMuon_RelIso_Barrel );
    assert( tmpSignalMuon_RelIso_RhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso05_Barrel );
    assert( tmpSignalMuon_RelIso05_RhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso_Endcap );
    assert( tmpSignalMuon_RelIso_RhoCorrected_Endcap );
    assert( tmpSignalMuon_RelIso05_Endcap );
    assert( tmpSignalMuon_RelIso05_RhoCorrected_Endcap );

    assert( tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );






    assert( tmpBkgMuon_RelIso_Barrel );
    assert( tmpBkgMuon_RelIso_RhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso05_Barrel );
    assert( tmpBkgMuon_RelIso05_RhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso_Endcap );
    assert( tmpBkgMuon_RelIso_RhoCorrected_Endcap );
    assert( tmpBkgMuon_RelIso05_Endcap );
    assert( tmpBkgMuon_RelIso05_RhoCorrected_Endcap );

    assert( tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );


    SignalMuon_RelIso_Barrel.push_back(tmpSignalMuon_RelIso_Barrel);
    SignalMuon_RelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_RelIsoRhoCorrected_Barrel);
    SignalMuon_RelIso05_Barrel.push_back(tmpSignalMuon_RelIso05_Barrel);
    SignalMuon_RelIso05RhoCorrected_Barrel.push_back(tmpSignalMuon_RelIso05RhoCorrected_Barrel);
    SignalMuon_RelIso_Endcap.push_back(tmpSignalMuon_RelIso_Endcap);
    SignalMuon_RelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_RelIsoRhoCorrected_Endcap);
    SignalMuon_RelIso05_Endcap.push_back(tmpSignalMuon_RelIso05_Endcap);
    SignalMuon_RelIso05RhoCorrected_Endcap.push_back(tmpSignalMuon_RelIso05RhoCorrected_Endcap);


    SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);




    BkgMuon_RelIso_Barrel.push_back(tmpBkgMuon_RelIso_Barrel);
    BkgMuon_RelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_RelIsoRhoCorrected_Barrel);
    BkgMuon_RelIso05_Barrel.push_back(tmpBkgMuon_RelIso05_Barrel);
    BkgMuon_RelIso05RhoCorrected_Barrel.push_back(tmpBkgMuon_RelIso05RhoCorrected_Barrel);
    BkgMuon_RelIso_Endcap.push_back(tmpBkgMuon_RelIso_Endcap);
    BkgMuon_RelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_RelIsoRhoCorrected_Endcap);
    BkgMuon_RelIso05_Endcap.push_back(tmpBkgMuon_RelIso05_Endcap);
    BkgMuon_RelIso05RhoCorrected_Endcap.push_back(tmpBkgMuon_RelIso05RhoCorrected_Endcap);


    BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);


  }

  
  TH1F *NVertex0To4_SignalMuon_RelIso_Barrel = (TH1F*)SignalMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05_Barrel = (TH1F*)SignalMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)SignalMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso_Endcap = (TH1F*)SignalMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05_Endcap = (TH1F*)SignalMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)SignalMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());


  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());





  TH1F *NVertex0To4_BkgMuon_RelIso_Barrel = (TH1F*)BkgMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05_Barrel = (TH1F*)BkgMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)BkgMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso_Endcap = (TH1F*)BkgMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05_Endcap = (TH1F*)BkgMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)BkgMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());

  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());





  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso_Barrel, SignalMuon_RelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIsoRhoCorrected_Barrel, SignalMuon_RelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05_Barrel, SignalMuon_RelIso05_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05RhoCorrected_Barrel, SignalMuon_RelIso05RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso_Endcap, SignalMuon_RelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIsoRhoCorrected_Endcap, SignalMuon_RelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05_Endcap, SignalMuon_RelIso05_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05RhoCorrected_Endcap, SignalMuon_RelIso05RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);






  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso_Barrel, BkgMuon_RelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel, BkgMuon_RelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05_Barrel, BkgMuon_RelIso05_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05RhoCorrected_Barrel, BkgMuon_RelIso05RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso_Endcap, BkgMuon_RelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap, BkgMuon_RelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05_Endcap, BkgMuon_RelIso05_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05RhoCorrected_Endcap, BkgMuon_RelIso05RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);













  TH1F *NVertex7To15_SignalMuon_RelIso_Barrel = (TH1F*)SignalMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05_Barrel = (TH1F*)SignalMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)SignalMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso_Endcap = (TH1F*)SignalMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05_Endcap = (TH1F*)SignalMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)SignalMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());


  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());





  TH1F *NVertex7To15_BkgMuon_RelIso_Barrel = (TH1F*)BkgMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05_Barrel = (TH1F*)BkgMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)BkgMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso_Endcap = (TH1F*)BkgMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05_Endcap = (TH1F*)BkgMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)BkgMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());

  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());





  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso_Barrel, SignalMuon_RelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, SignalMuon_RelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05_Barrel, SignalMuon_RelIso05_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel, SignalMuon_RelIso05RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso_Endcap, SignalMuon_RelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, SignalMuon_RelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05_Endcap, SignalMuon_RelIso05_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap, SignalMuon_RelIso05RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7, 15);






  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso_Barrel, BkgMuon_RelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIsoRhoCorrected_Barrel, BkgMuon_RelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05_Barrel, BkgMuon_RelIso05_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05RhoCorrected_Barrel, BkgMuon_RelIso05RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso_Endcap, BkgMuon_RelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIsoRhoCorrected_Endcap, BkgMuon_RelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05_Endcap, BkgMuon_RelIso05_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05RhoCorrected_Endcap, BkgMuon_RelIso05RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone03_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone03_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7, 15);










  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  TFile *canvasFile = new TFile("IsolationPerformancePlots.muons.root","UPDATE");

  vector<TGraphAsymmErrors*> SignalEffVsCutGraphs;
  vector<TGraphAsymmErrors*> SignalEffVsBkgEffGraphs;
  vector<string> GraphLabels;
  string plotname;


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  string tmpLabel;
//   TPaveLabel *LabelText;
  TLegend * legend;
  TGraphAsymmErrors *tmpSigEffVsCut = 0;
  TGraphAsymmErrors *tmpSigEffVsBkgEff = 0;


  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex0To4_PFConesize_01Threshold";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso05_Barrel, "SigEffVsCut_Muon_RelIso05_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso05_Barrel, NVertex0To4_BkgMuon_RelIso05_Barrel, "SigEffVsBkgEff_Muon_RelIso05_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.5 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.4 Cone");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
 

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
  




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex0To4_dZCut_01Threshold";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso No dZ Cut");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso |dZ|<0.1cm");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 







  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_PFConesize_01Threshold";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05_Barrel, "SigEffVsCut_Muon_RelIso05_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05_Barrel, NVertex7To15_BkgMuon_RelIso05_Barrel, "SigEffVsBkgEff_Muon_RelIso05_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.5 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.4 Cone");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
 

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
  




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_dZCut_01Threshold";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso No dZ Cut");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso |dZ|<0.1cm");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 





  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_RhoCorrection";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Barrel, "SigEffVsCut_Muon_RelIso_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Barrel, NVertex7To15_BkgMuon_RelIso_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso ");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_PtThresholds";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>0.1");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>0.5");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5");





 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 





  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_RhoCorrectionWithPtThresholds";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0 EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5 EA Corr");




 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 



  //*************************************************************************************
  //Best Choices Barrel
  //*************************************************************************************


  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronBarrel_BestChoices_NVertex0To4";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_RhoCorrected_Barrel, "SigEffVsCut_Muon_RelIso_RhoCorrected_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_RhoCorrected_Barrel, NVertex0To4_BkgMuon_RelIso_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "P");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    

  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronBarrel_BestChoices_NVertex7To15";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex7To15_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Barrel, "SigEffVsCut_Muon_RelIso_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Barrel, NVertex7To15_BkgMuon_RelIso_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "P");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    



























  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex0To4_PFConesize_01Threshold";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso05_Endcap, "SigEffVsCut_Muon_RelIso05_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso05_Endcap, NVertex0To4_BkgMuon_RelIso05_Endcap, "SigEffVsBkgEff_Muon_RelIso05_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.5 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.4 Cone");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
 

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
  




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex0To4_dZCut_01Threshold";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso No dZ Cut");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso |dZ|<0.1cm");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 







  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_PFConesize_01Threshold";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05_Endcap, "SigEffVsCut_Muon_RelIso05_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05_Endcap, NVertex7To15_BkgMuon_RelIso05_Endcap, "SigEffVsBkgEff_Muon_RelIso05_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso 0.5 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.3 Cone");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso 0.4 Cone");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
 

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
  




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_dZCut_01Threshold";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso No dZ Cut");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso |dZ|<0.1cm");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 





  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_RhoCorrection";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Endcap, "SigEffVsCut_Muon_RelIso_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Endcap, NVertex7To15_BkgMuon_RelIso_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso ");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_PtThresholds";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>0.1");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_05Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>0.5");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5");





 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 





  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_RhoCorrectionWithPtThresholds";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0 EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.5 EA Corr");




 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
 canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->AddEntry(SignalEffVsCutGraphs[0],GraphLabels[0].c_str(), "LP");
  SignalEffVsCutGraphs[0]->Draw("Psame");
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsCut_" + plotname + label).c_str(), "WriteDelete") ;
 



  //*************************************************************************************
  //Best Choices Endcap
  //*************************************************************************************


  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronEndcap_BestChoices_NVertex0To4";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso_RhoCorrected_Endcap, "SigEffVsCut_Muon_RelIso_RhoCorrected_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso_RhoCorrected_Endcap, NVertex0To4_BkgMuon_RelIso_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "P");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    

  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronEndcap_BestChoices_NVertex7To15";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex7To15_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Endcap, "SigEffVsCut_Muon_RelIso_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_RhoCorrected_Endcap, NVertex7To15_BkgMuon_RelIso_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIso_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Std Iso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso EA Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, NVertex7To15_BkgMuon_TotalPFRelIso_Cone04_10Threshold_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFIso NeuPt>1.0");

 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.75);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "P");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    








  canvasFile->Close();
  return;






}




void PlotMuonIsolationPerformance(Int_t i) {



  if (i==0) {
    MakeIsolationPerformancePlots("HWW130SignalDataEle10Jet30Bkg_Pt10To20", "HWW130_Pt10To20", "Ele10Jet30_Pt10To20", 0.10);
    PlotMuonIsoVsNVertices("HWW130SignalDataEle10Jet30Bkg_Pt10To20", "HWW130_Pt10To20", "Ele10Jet30_Pt10To20");
  }
  
  if (i==1) {
    MakeIsolationPerformancePlots("HWW130SignalDataEle10Jet30Bkg_Pt20To30", "HWW130_Pt20To30", "Ele10Jet30_Pt20To30", 0.10);
    PlotMuonIsoVsNVertices("HWW130SignalDataEle10Jet30Bkg_Pt20To30", "HWW130_Pt20To30", "Ele10Jet30_Pt20To30");
  }

  if (i==10) {
    MakeIsolationPerformancePlots("ZeeSignalDataEle10Jet30Bkg_Pt10To20", "Zee_Pt10To20", "Ele10Jet30_Pt10To20", 0.10);
    PlotMuonIsoVsNVertices("ZeeSignalDataEle10Jet30Bkg_Pt10To20", "Zee_Pt10To20", "Ele10Jet30_Pt10To20");
  }
  if (i==11) {
    MakeIsolationPerformancePlots("ZeeSignalDataEle10Jet30Bkg_Pt20To30", "Zee_Pt20To30", "Ele10Jet30_Pt20To30", 0.10);
    PlotMuonIsoVsNVertices("ZeeSignalDataEle10Jet30Bkg_Pt20To30", "Zee_Pt20To30", "Ele10Jet30_Pt20To30");
  }



}
