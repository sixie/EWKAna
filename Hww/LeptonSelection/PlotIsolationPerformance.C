//root -l EWKAna/Hww/LeptonSelection/PlotIsolationPerformance.C+
 
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





void PlotMuonIsolationPerformance(string Label = "", string signalSample = "Zmm", string bkgSample = "WJetsMC" , Double_t currentWPCut = 0.15) {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kGreen);
  colors.push_back(kRed+3);  

  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.muons.root", "READ");
    
  vector<TH1F*>   SignalMuon_RelIso_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   SignalMuon_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   SignalMuon_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_RelIso05_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso04_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso04_Barrel;
  vector<TH1F*>   SignalMuon_RelIso05RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIso04RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso05_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso04_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso04_Endcap;
  vector<TH1F*>   SignalMuon_RelIso05RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso04RhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap;

  vector<TH1F*>   BkgMuon_RelIso_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   BkgMuon_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   BkgMuon_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_RelIso05_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso04_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso04_Barrel;
  vector<TH1F*>   BkgMuon_RelIso05RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso04RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso05_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso04_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso04_Endcap;
  vector<TH1F*>   BkgMuon_RelIso05RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso04RhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap;

  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalMuon_RelIso_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_RelIso05_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso04_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso04_Barrel;
    TH1F *tmpSignalMuon_RelIso05RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso04RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso05_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso04_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso04_Endcap;
    TH1F *tmpSignalMuon_RelIso05RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso04RhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap;

    TH1F *tmpBkgMuon_RelIso_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_RelIso05_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso04_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso04_Barrel;
    TH1F *tmpBkgMuon_RelIso05RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso04RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso05_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso04_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso04_Endcap;
    TH1F *tmpBkgMuon_RelIso05RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso04RhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap;

    tmpSignalMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05_Barrel = (TH1F*)f->Get((string("Muon_RelIso05_Barrel_") + IntToString(n) + "_" + signalSample).c_str());

    cout << (string("Muon_RelIso05_Barrel_") + IntToString(n) + "_" + signalSample).c_str() << endl;

    tmpSignalMuon_TotalPFRelIso04_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso04_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05_Endcap = (TH1F*)f->Get((string("Muon_RelIso05_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso04_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso04_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
 
    tmpBkgMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05_Barrel = (TH1F*)f->Get((string("Muon_RelIso05_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso04_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso04_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05_Endcap = (TH1F*)f->Get((string("Muon_RelIso05_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso04_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso04_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIso05RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
 
    cout << "load " << n << endl;
    assert( tmpSignalMuon_RelIso_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Barrel );
    assert( tmpSignalMuon_RelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Endcap );
    assert( tmpSignalMuon_RelIsoRhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalMuon_RelIso05_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso04_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso04_Barrel );
    assert( tmpSignalMuon_RelIso05RhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso04RhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso05_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso04_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso04_Endcap );
    assert( tmpSignalMuon_RelIso05RhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso04RhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap );

    assert( tmpBkgMuon_RelIso_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Barrel );
    assert( tmpBkgMuon_RelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Endcap );
    assert( tmpBkgMuon_RelIsoRhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgMuon_RelIso05_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso04_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso04_Barrel );
    assert( tmpBkgMuon_RelIso05RhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso04RhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso05_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso04_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso04_Endcap );
    assert( tmpBkgMuon_RelIso05RhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso04RhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap );

    SignalMuon_RelIso_Barrel.push_back(tmpSignalMuon_RelIso_Barrel);
    SignalMuon_TotalPFRelIso_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Barrel);
    SignalMuon_RelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_RelIsoRhoCorrected_Barrel);
    SignalMuon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    SignalMuon_RelIso_Endcap.push_back(tmpSignalMuon_RelIso_Endcap);
    SignalMuon_TotalPFRelIso_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Endcap);
    SignalMuon_RelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_RelIsoRhoCorrected_Endcap);
    SignalMuon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    SignalMuon_RelIso05_Barrel.push_back(tmpSignalMuon_RelIso05_Barrel);
    SignalMuon_TotalPFRelIso04_Barrel.push_back(tmpSignalMuon_TotalPFRelIso04_Barrel);
    SignalMuon_VertexSelectedPFRelIso04_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso04_Barrel);
    SignalMuon_RelIso05RhoCorrected_Barrel.push_back(tmpSignalMuon_RelIso05RhoCorrected_Barrel);
    SignalMuon_TotalPFRelIso04RhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIso04RhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel);
    SignalMuon_RelIso05_Endcap.push_back(tmpSignalMuon_RelIso05_Endcap);
    SignalMuon_TotalPFRelIso04_Endcap.push_back(tmpSignalMuon_TotalPFRelIso04_Endcap);
    SignalMuon_VertexSelectedPFRelIso04_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso04_Endcap);
    SignalMuon_RelIso05RhoCorrected_Endcap.push_back(tmpSignalMuon_RelIso05RhoCorrected_Endcap);
    SignalMuon_TotalPFRelIso04RhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIso04RhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap);
    BkgMuon_RelIso_Barrel.push_back(tmpBkgMuon_RelIso_Barrel);
    BkgMuon_TotalPFRelIso_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Barrel);
    BkgMuon_RelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_RelIsoRhoCorrected_Barrel);
    BkgMuon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    BkgMuon_RelIso_Endcap.push_back(tmpBkgMuon_RelIso_Endcap);
    BkgMuon_TotalPFRelIso_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Endcap);
    BkgMuon_RelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_RelIsoRhoCorrected_Endcap);
    BkgMuon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    BkgMuon_RelIso05_Barrel.push_back(tmpBkgMuon_RelIso05_Barrel);
    BkgMuon_TotalPFRelIso04_Barrel.push_back(tmpBkgMuon_TotalPFRelIso04_Barrel);
    BkgMuon_VertexSelectedPFRelIso04_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso04_Barrel);
    BkgMuon_RelIso05RhoCorrected_Barrel.push_back(tmpBkgMuon_RelIso05RhoCorrected_Barrel);
    BkgMuon_TotalPFRelIso04RhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIso04RhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel);
    BkgMuon_RelIso05_Endcap.push_back(tmpBkgMuon_RelIso05_Endcap);
    BkgMuon_TotalPFRelIso04_Endcap.push_back(tmpBkgMuon_TotalPFRelIso04_Endcap);
    BkgMuon_VertexSelectedPFRelIso04_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso04_Endcap);
    BkgMuon_RelIso05RhoCorrected_Endcap.push_back(tmpBkgMuon_RelIso05RhoCorrected_Endcap);
    BkgMuon_TotalPFRelIso04RhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIso04RhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap);
  }


  TH1F *NVertex0To4_SignalMuon_RelIso_Barrel = (TH1F*)SignalMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso_Endcap = (TH1F*)SignalMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05_Barrel = (TH1F*)SignalMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso04_Barrel = (TH1F*)SignalMuon_TotalPFRelIso04_Barrel[0]->Clone((string("Muon_TotalPFRelIso04_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso04_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)SignalMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05_Endcap = (TH1F*)SignalMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso04_Endcap = (TH1F*)SignalMuon_TotalPFRelIso04_Endcap[0]->Clone((string("Muon_TotalPFRelIso04_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso04_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)SignalMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
  TH1F *NVertex0To4_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());

  TH1F *NVertex0To4_BkgMuon_RelIso_Barrel = (TH1F*)BkgMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso_Endcap = (TH1F*)BkgMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05_Barrel = (TH1F*)BkgMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso04_Barrel = (TH1F*)BkgMuon_TotalPFRelIso04_Barrel[0]->Clone((string("Muon_TotalPFRelIso04_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso04_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)BkgMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05_Endcap = (TH1F*)BkgMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso04_Endcap = (TH1F*)BkgMuon_TotalPFRelIso04_Endcap[0]->Clone((string("Muon_TotalPFRelIso04_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso04_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)BkgMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
  TH1F *NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex0To4_") + bkgSample).c_str());
 
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso_Barrel, SignalMuon_RelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Barrel, SignalMuon_TotalPFRelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Barrel, SignalMuon_VertexSelectedPFRelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIsoRhoCorrected_Barrel, SignalMuon_RelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso_Endcap, SignalMuon_RelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso_Endcap, SignalMuon_TotalPFRelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Endcap, SignalMuon_VertexSelectedPFRelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIsoRhoCorrected_Endcap, SignalMuon_RelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05_Barrel, SignalMuon_RelIso05_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso04_Barrel, SignalMuon_TotalPFRelIso04_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Barrel, SignalMuon_VertexSelectedPFRelIso04_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05RhoCorrected_Barrel, SignalMuon_RelIso05RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05_Endcap, SignalMuon_RelIso05_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso04_Endcap, SignalMuon_TotalPFRelIso04_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Endcap, SignalMuon_VertexSelectedPFRelIso04_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_RelIso05RhoCorrected_Endcap, SignalMuon_RelIso05RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, 0, 4);

  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso_Barrel, BkgMuon_RelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Barrel, BkgMuon_TotalPFRelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Barrel, BkgMuon_VertexSelectedPFRelIso_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel, BkgMuon_RelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso_Endcap, BkgMuon_RelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso_Endcap, BkgMuon_TotalPFRelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Endcap, BkgMuon_VertexSelectedPFRelIso_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap, BkgMuon_RelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05_Barrel, BkgMuon_RelIso05_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso04_Barrel, BkgMuon_TotalPFRelIso04_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Barrel, BkgMuon_VertexSelectedPFRelIso04_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05RhoCorrected_Barrel, BkgMuon_RelIso05RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Barrel, BkgMuon_TotalPFRelIso04RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05_Endcap, BkgMuon_RelIso05_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso04_Endcap, BkgMuon_TotalPFRelIso04_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Endcap, BkgMuon_VertexSelectedPFRelIso04_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_RelIso05RhoCorrected_Endcap, BkgMuon_RelIso05RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Endcap, BkgMuon_TotalPFRelIso04RhoCorrected_Endcap, 0, 4);
  AddIsoHistograms(NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, 0, 4);


  TH1F *NVertex7To15_SignalMuon_RelIso_Barrel = (TH1F*)SignalMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Barrel = (TH1F*)SignalMuon_TotalPFRelIso_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso_Endcap = (TH1F*)SignalMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso_Endcap = (TH1F*)SignalMuon_TotalPFRelIso_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05_Barrel = (TH1F*)SignalMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso04_Barrel = (TH1F*)SignalMuon_TotalPFRelIso04_Barrel[0]->Clone((string("Muon_TotalPFRelIso04_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso04_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel = (TH1F*)SignalMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)SignalMuon_TotalPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05_Endcap = (TH1F*)SignalMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso04_Endcap = (TH1F*)SignalMuon_TotalPFRelIso04_Endcap[0]->Clone((string("Muon_TotalPFRelIso04_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso04_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap = (TH1F*)SignalMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)SignalMuon_TotalPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  TH1F *NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
  
  TH1F *NVertex7To15_BkgMuon_RelIso_Barrel = (TH1F*)BkgMuon_RelIso_Barrel[0]->Clone((string("Muon_RelIso_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Barrel = (TH1F*)BkgMuon_TotalPFRelIso_Barrel[0]->Clone((string("Muon_TotalPFRelIso_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_RelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_RelIsoRhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso_Endcap = (TH1F*)BkgMuon_RelIso_Endcap[0]->Clone((string("Muon_RelIso_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso_Endcap = (TH1F*)BkgMuon_TotalPFRelIso_Endcap[0]->Clone((string("Muon_TotalPFRelIso_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_RelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_RelIsoRhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05_Barrel = (TH1F*)BkgMuon_RelIso05_Barrel[0]->Clone((string("Muon_RelIso05_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso04_Barrel = (TH1F*)BkgMuon_TotalPFRelIso04_Barrel[0]->Clone((string("Muon_TotalPFRelIso04_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso04_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso04_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05RhoCorrected_Barrel = (TH1F*)BkgMuon_RelIso05RhoCorrected_Barrel[0]->Clone((string("Muon_RelIso05RhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso04RhoCorrected_Barrel = (TH1F*)BkgMuon_TotalPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel = (TH1F*)BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05_Endcap = (TH1F*)BkgMuon_RelIso05_Endcap[0]->Clone((string("Muon_RelIso05_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso04_Endcap = (TH1F*)BkgMuon_TotalPFRelIso04_Endcap[0]->Clone((string("Muon_TotalPFRelIso04_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso04_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso04_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_RelIso05RhoCorrected_Endcap = (TH1F*)BkgMuon_RelIso05RhoCorrected_Endcap[0]->Clone((string("Muon_RelIso05RhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_TotalPFRelIso04RhoCorrected_Endcap = (TH1F*)BkgMuon_TotalPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
  TH1F *NVertex7To15_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap = (TH1F*)BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap[0]->Clone((string("Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex7To15_") + bkgSample).c_str());
 
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso_Barrel, SignalMuon_RelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Barrel, SignalMuon_TotalPFRelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel, SignalMuon_VertexSelectedPFRelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, SignalMuon_RelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso_Endcap, SignalMuon_RelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso_Endcap, SignalMuon_TotalPFRelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap, SignalMuon_VertexSelectedPFRelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, SignalMuon_RelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05_Barrel, SignalMuon_RelIso05_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso04_Barrel, SignalMuon_TotalPFRelIso04_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Barrel, SignalMuon_VertexSelectedPFRelIso04_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel, SignalMuon_RelIso05RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05_Endcap, SignalMuon_RelIso05_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso04_Endcap, SignalMuon_TotalPFRelIso04_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Endcap, SignalMuon_VertexSelectedPFRelIso04_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap, SignalMuon_RelIso05RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, 7, 15);

  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso_Barrel, BkgMuon_RelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Barrel, BkgMuon_TotalPFRelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Barrel, BkgMuon_VertexSelectedPFRelIso_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIsoRhoCorrected_Barrel, BkgMuon_RelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso_Endcap, BkgMuon_RelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso_Endcap, BkgMuon_TotalPFRelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso_Endcap, BkgMuon_VertexSelectedPFRelIso_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIsoRhoCorrected_Endcap, BkgMuon_RelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05_Barrel, BkgMuon_RelIso05_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso04_Barrel, BkgMuon_TotalPFRelIso04_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso04_Barrel, BkgMuon_VertexSelectedPFRelIso04_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05RhoCorrected_Barrel, BkgMuon_RelIso05RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso04RhoCorrected_Barrel, BkgMuon_TotalPFRelIso04RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05_Endcap, BkgMuon_RelIso05_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso04_Endcap, BkgMuon_TotalPFRelIso04_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso04_Endcap, BkgMuon_VertexSelectedPFRelIso04_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_RelIso05RhoCorrected_Endcap, BkgMuon_RelIso05RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_TotalPFRelIso04RhoCorrected_Endcap, BkgMuon_TotalPFRelIso04RhoCorrected_Endcap, 7, 15);
  AddIsoHistograms(NVertex7To15_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, 7, 15);






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
  TGraphAsymmErrors* tmpSigEffVsCut = 0;
  TGraphAsymmErrors* tmpSigEffVsBkgEff = 0;

  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex0To4";
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
  GraphLabels.push_back("Standard RelIso");

 //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Barrel_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Barrel_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");

  //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso05_Barrel, "SigEffVsCut_Muon_RelIso05_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso05_Barrel, NVertex0To4_BkgMuon_RelIso05_Barrel, "SigEffVsBkgEff_Muon_RelIso05_Barrel_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso05");

 //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso04_Barrel, "SigEffVsCut_Muon_TotalPFRelIso04_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso04_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso04_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso04_Barrel_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso04");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso04_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04_Barrel_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso04");

 //*******************************************************************************************
 
  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {

    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");
    
    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
    
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
    
    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i-1]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i-1]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
    
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.0,1.0);    
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
  plotname = "MuonBarrel_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");




  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05_Barrel, "SigEffVsCut_Muon_RelIso05_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05_Barrel, NVertex0To4_BkgMuon_RelIso05_Barrel, "SigEffVsBkgEff_Muon_RelIso05_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso05");

 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso04_Barrel, "SigEffVsCut_Muon_TotalPFRelIso04_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso04_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso04_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso04_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso04");




  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso04_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso04");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    



  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex0To4";
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
  GraphLabels.push_back("Standard RelIso");

 //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Endcap_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Endcap_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");

  //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_RelIso05_Endcap, "SigEffVsCut_Muon_RelIso05_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_RelIso05_Endcap, NVertex0To4_BkgMuon_RelIso05_Endcap, "SigEffVsBkgEff_Muon_RelIso05_Endcap_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso05");

 //*************************************************************************************

  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_TotalPFRelIso04_Endcap, "SigEffVsCut_Muon_TotalPFRelIso04_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_TotalPFRelIso04_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso04_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso04_Endcap_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso04");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso04_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalMuon_VertexSelectedPFRelIso04_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04_Endcap_NVertex0To4");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso04");

 //*******************************************************************************************
 
  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {

    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");
    
    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
    
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
    
    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i-1]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i-1]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
    
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.0,1.0);    
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
  plotname = "MuonEndcap_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");




  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05_Endcap, "SigEffVsCut_Muon_RelIso05_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05_Endcap, NVertex0To4_BkgMuon_RelIso05_Endcap, "SigEffVsBkgEff_Muon_RelIso05_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso05");

 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso04_Endcap, "SigEffVsCut_Muon_TotalPFRelIso04_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso04_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso04_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso04_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso04");




  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso04_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso04");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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






  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_Corrections";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Barrel, "SigEffVsCut_Muon_TotalPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Barrel_NVertex7To15");


  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_RelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
   



  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonBarrel_NVertex7To15_ConeSizeWithCorrections";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Barrel, "CurrentWP_SigEffVsCut_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Barrel, NVertex0To4_BkgMuon_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Barrel_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_RelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_RelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel, "SigEffVsCut_Muon_RelIso05RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Barrel, NVertex0To4_BkgMuon_RelIso05RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_RelIso05RhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso05 Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, "SigEffVsCut_Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Barrel, NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_TotalPFRelIso04RhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso04 Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, "SigEffVsCut_Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Barrel, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04RhoCorrected_Barrel_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
   




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_Corrections";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Endcap, "SigEffVsCut_Muon_TotalPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso_Endcap_NVertex7To15");


  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_RelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonEndcap_NVertex7To15_ConeSizeWithCorrections";
  //*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso_Endcap, "CurrentWP_SigEffVsCut_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso_Endcap, NVertex0To4_BkgMuon_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Muon_RelIso_Endcap_NVertex7To15", currentWPCut);

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_RelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_RelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap, "SigEffVsCut_Muon_RelIso05RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_RelIso05RhoCorrected_Endcap, NVertex0To4_BkgMuon_RelIso05RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_RelIso05RhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso05 Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, "SigEffVsCut_Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_TotalPFRelIso04RhoCorrected_Endcap, NVertex0To4_BkgMuon_TotalPFRelIso04RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_TotalPFRelIso04RhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso04 Corr");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, "SigEffVsCut_Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, NVertex0To4_BkgMuon_VertexSelectedPFRelIso04RhoCorrected_Endcap, "SigEffVsBkgEff_Muon_VertexSelectedPFRelIso04RhoCorrected_Endcap_NVertex7To15");

  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==1) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }

  SignalEffVsBkgEffGraphs[0]->Draw("Psame");
  legend->AddEntry(SignalEffVsBkgEffGraphs[0],GraphLabels[0].c_str(), "LP");
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
  canvasFile->WriteTObject(cv,("SigEffVsBkgEff_" + plotname + label).c_str(), "WriteDelete") ;
   




  canvasFile->Close();







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
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   SignalMuon_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalMuon_RelIso_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIso_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   SignalMuon_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;

  vector<TH1F*>   BkgMuon_RelIso_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   BkgMuon_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgMuon_RelIso_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIso_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   BkgMuon_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;

  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalMuon_RelIso_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIso_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalMuon_RelIso_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIso_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpSignalMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;

    TH1F *tmpBkgMuon_RelIso_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIso_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgMuon_RelIso_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIso_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpBkgMuon_RelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;


    tmpSignalMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
 
    tmpBkgMuon_RelIso_Barrel = (TH1F*)f->Get((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIso_Endcap = (TH1F*)f->Get((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
 
    cout << "load " << n << endl;
    assert( tmpSignalMuon_RelIso_Barrel );
    assert( tmpSignalMuon_TotalPFRelIso_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Barrel );
    assert( tmpSignalMuon_RelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalMuon_RelIso_Endcap );
    assert( tmpSignalMuon_TotalPFRelIso_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIso_Endcap );
    assert( tmpSignalMuon_RelIsoRhoCorrected_Endcap );
    assert( tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap );

    assert( tmpBkgMuon_RelIso_Barrel );
    assert( tmpBkgMuon_TotalPFRelIso_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Barrel );
    assert( tmpBkgMuon_RelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgMuon_RelIso_Endcap );
    assert( tmpBkgMuon_TotalPFRelIso_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIso_Endcap );
    assert( tmpBkgMuon_RelIsoRhoCorrected_Endcap );
    assert( tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap );

    SignalMuon_RelIso_Barrel.push_back(tmpSignalMuon_RelIso_Barrel);
    SignalMuon_TotalPFRelIso_Barrel.push_back(tmpSignalMuon_TotalPFRelIso_Barrel);
    SignalMuon_VertexSelectedPFRelIso_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Barrel);
    SignalMuon_RelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_RelIsoRhoCorrected_Barrel);
    SignalMuon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_TotalPFRelIsoRhoCorrected_Barrel);
    SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    SignalMuon_RelIso_Endcap.push_back(tmpSignalMuon_RelIso_Endcap);
    SignalMuon_TotalPFRelIso_Endcap.push_back(tmpSignalMuon_TotalPFRelIso_Endcap);
    SignalMuon_VertexSelectedPFRelIso_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIso_Endcap);
    SignalMuon_RelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_RelIsoRhoCorrected_Endcap);
    SignalMuon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_TotalPFRelIsoRhoCorrected_Endcap);
    SignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    BkgMuon_RelIso_Barrel.push_back(tmpBkgMuon_RelIso_Barrel);
    BkgMuon_TotalPFRelIso_Barrel.push_back(tmpBkgMuon_TotalPFRelIso_Barrel);
    BkgMuon_VertexSelectedPFRelIso_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Barrel);
    BkgMuon_RelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_RelIsoRhoCorrected_Barrel);
    BkgMuon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_TotalPFRelIsoRhoCorrected_Barrel);
    BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    BkgMuon_RelIso_Endcap.push_back(tmpBkgMuon_RelIso_Endcap);
    BkgMuon_TotalPFRelIso_Endcap.push_back(tmpBkgMuon_TotalPFRelIso_Endcap);
    BkgMuon_VertexSelectedPFRelIso_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIso_Endcap);
    BkgMuon_RelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_RelIsoRhoCorrected_Endcap);
    BkgMuon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_TotalPFRelIsoRhoCorrected_Endcap);
    BkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);

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




void PlotElectronIsoVsNVertices(string Label = "", string signalSample = "Zmm", string bkgSample = "WJetsMC") {

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
  
  vector<TH1F*>   SignalElectron_RelIso_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
  vector<TH1F*>   SignalElectron_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_RelIso_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
  vector<TH1F*>   SignalElectron_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;
 
  vector<TH1F*>   BkgElectron_RelIso_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
  vector<TH1F*>   BkgElectron_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_RelIso_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
  vector<TH1F*>   BkgElectron_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;
 
  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalElectron_RelIso_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
    TH1F *tmpSignalElectron_RelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_RelIso_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
    TH1F *tmpSignalElectron_RelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;

    TH1F *tmpBkgElectron_RelIso_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
    TH1F *tmpBkgElectron_RelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_RelIso_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
    TH1F *tmpBkgElectron_RelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;

    tmpSignalElectron_RelIso_Barrel = (TH1F*)f->Get((string("Electron_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso_Endcap = (TH1F*)f->Get((string("Electron_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());


    tmpBkgElectron_RelIso_Barrel = (TH1F*)f->Get((string("Electron_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso_Endcap = (TH1F*)f->Get((string("Electron_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());


    cout << "load " << n << endl;
    assert( tmpSignalElectron_RelIso_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel );
    assert( tmpSignalElectron_RelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_RelIso_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap );
    assert( tmpSignalElectron_RelIsoRhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_RelIso_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel );
    assert( tmpBkgElectron_RelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_RelIso_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap );
    assert( tmpBkgElectron_RelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap );

    SignalElectron_RelIso_Barrel.push_back(tmpSignalElectron_RelIso_Barrel);
    SignalElectron_TotalPFRelIso_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel);
    SignalElectron_RelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_RelIsoRhoCorrected_Barrel);
    SignalElectron_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIsoRhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel);
    SignalElectron_RelIso_Endcap.push_back(tmpSignalElectron_RelIso_Endcap);
    SignalElectron_TotalPFRelIso_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap);
    SignalElectron_RelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_RelIsoRhoCorrected_Endcap);
    SignalElectron_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIsoRhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap);

    BkgElectron_RelIso_Barrel.push_back(tmpBkgElectron_RelIso_Barrel);
    BkgElectron_TotalPFRelIso_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel);
    BkgElectron_RelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_RelIsoRhoCorrected_Barrel);
    BkgElectron_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIsoRhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel);
    BkgElectron_RelIso_Endcap.push_back(tmpBkgElectron_RelIso_Endcap);
    BkgElectron_TotalPFRelIso_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap);
    BkgElectron_RelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_RelIsoRhoCorrected_Endcap);
    BkgElectron_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIsoRhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap);

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
  plotname = "ElectronIsolationEfficiency_Barrel_Signal";
  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_TotalPFRelIso_Barrel, "EffVsNVertex_TotalPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_VertexSelectedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_TotalPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_TotalPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF Corr");

  //*************************************************************************************
   tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
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
  plotname = "ElectronIsolationEfficiency_FootprintRemoved_Barrel_Signal";
  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_FootprintRemovedPFRelIso_Barrel, "EffVsNVertex_FootprintRemovedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedFootprintRemovedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_FootprintRemovedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed Corr");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PFNoPU Corr");


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
  plotname = "ElectronIsolationEfficiency_Barrel_Bkg";
  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_TotalPFRelIso_Barrel, "EffVsNVertex_TotalPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_VertexSelectedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_TotalPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_TotalPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("PF Corr");

  //*************************************************************************************
   tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
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
  plotname = "ElectronIsolationEfficiency_FootprintRemoved_Barrel_Bkg";
  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_RelIso_Barrel, "EffVsNVertex_RelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_FootprintRemovedPFRelIso_Barrel, "EffVsNVertex_FootprintRemovedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PF");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel, "EffVsNVertex_VertexSelectedFootprintRemovedPFRelIso_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PFNoPU");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_RelIsoRhoCorrected_Barrel, "EffVsNVertex_RelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Standard Corr");


  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_FootprintRemovedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed Corr");

  //*************************************************************************************
  tmpEffVsNVertex = MakeEffVsNVertex(BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel, "EffVsNVertex_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel", 0.10, 1);
  EffVsNVertices.push_back(tmpEffVsNVertex);
  GraphLabels.push_back("Footprint Removed PFNoPU Corr");


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


















void PlotElectronIsolationPerformance(string Label = "", string signalSample = "Zmm", string bkgSample = "WJetsMC", Double_t currentWPCut = 0.15) {

  string label = Label;
  if (Label != "") label = "_" + Label;

  vector<Int_t> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kMagenta);
  colors.push_back(kCyan);
  colors.push_back(kBlack);
  colors.push_back(kGreen);
  
  TFile *f = new TFile("HwwSelectionPlots_LeptonEfficiency.electrons.root", "READ");
    
  vector<TH1F*>   SignalElectron_RelIso_Barrel;
  vector<TH1F*>   SignalElectron_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_RelIso_Endcap;
  vector<TH1F*>   SignalElectron_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_RelIso04_Barrel;
  vector<TH1F*>   SignalElectron_RelIso04RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_RelIso04_Endcap;
  vector<TH1F*>   SignalElectron_RelIso04RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;





  vector<TH1F*>   BkgElectron_RelIso_Barrel;
  vector<TH1F*>   BkgElectron_RelIsoRhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_RelIso_Endcap;
  vector<TH1F*>   BkgElectron_RelIsoRhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_RelIso04_Barrel;
  vector<TH1F*>   BkgElectron_RelIso04RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_RelIso04_Endcap;
  vector<TH1F*>   BkgElectron_RelIso04RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
  vector<TH1F*>   BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
  vector<TH1F*>   BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;





  for (int n=0; n < 20; ++n) {     

    TH1F *tmpSignalElectron_RelIso_Barrel;
    TH1F *tmpSignalElectron_RelIsoRhoCorrected_Barrel;
    TH1F *tmpSignalElectron_RelIso_Endcap;
    TH1F *tmpSignalElectron_RelIsoRhoCorrected_Endcap;
    TH1F *tmpSignalElectron_RelIso04_Barrel;
    TH1F *tmpSignalElectron_RelIso04RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_RelIso04_Endcap;
    TH1F *tmpSignalElectron_RelIso04RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;



    TH1F *tmpBkgElectron_RelIso_Barrel;
    TH1F *tmpBkgElectron_RelIsoRhoCorrected_Barrel;
    TH1F *tmpBkgElectron_RelIso_Endcap;
    TH1F *tmpBkgElectron_RelIsoRhoCorrected_Endcap;
    TH1F *tmpBkgElectron_RelIso04_Barrel;
    TH1F *tmpBkgElectron_RelIso04RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_RelIso04_Endcap;
    TH1F *tmpBkgElectron_RelIso04RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap;
    TH1F *tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;
    TH1F *tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap;





    tmpSignalElectron_RelIso_Barrel = (TH1F*)f->Get((string("Electron_RelIso_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso_Endcap = (TH1F*)f->Get((string("Electron_RelIso_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso04_Barrel = (TH1F*)f->Get((string("Electron_RelIso04_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso04_Endcap = (TH1F*)f->Get((string("Electron_RelIso04_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_RelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
     tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());
    tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + signalSample).c_str());





    tmpBkgElectron_RelIso_Barrel = (TH1F*)f->Get((string("Electron_RelIso_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIsoRhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso_Endcap = (TH1F*)f->Get((string("Electron_RelIso_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIsoRhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso04_Barrel = (TH1F*)f->Get((string("Electron_RelIso04_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso04RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_RelIso04RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso04_Endcap = (TH1F*)f->Get((string("Electron_RelIso04_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_RelIso04RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_RelIso04RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
     tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
    tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)f->Get((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_") + IntToString(n) + "_" + bkgSample).c_str());
 


    cout << "load " << n << endl;
    assert( tmpSignalElectron_RelIso_Barrel );
    assert( tmpSignalElectron_RelIsoRhoCorrected_Barrel );
    assert( tmpSignalElectron_RelIso_Endcap );
    assert( tmpSignalElectron_RelIsoRhoCorrected_Endcap );
    assert( tmpSignalElectron_RelIso04_Barrel );
    assert( tmpSignalElectron_RelIso04RhoCorrected_Barrel );
    assert( tmpSignalElectron_RelIso04_Endcap );
    assert( tmpSignalElectron_RelIso04RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );




    assert( tmpBkgElectron_RelIso_Barrel );
    assert( tmpBkgElectron_RelIsoRhoCorrected_Barrel );
    assert( tmpBkgElectron_RelIso_Endcap );
    assert( tmpBkgElectron_RelIsoRhoCorrected_Endcap );
    assert( tmpBkgElectron_RelIso04_Barrel );
    assert( tmpBkgElectron_RelIso04RhoCorrected_Barrel );
    assert( tmpBkgElectron_RelIso04_Endcap );
    assert( tmpBkgElectron_RelIso04RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap );
    assert( tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );
    assert( tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap );




    SignalElectron_RelIso_Barrel.push_back(tmpSignalElectron_RelIso_Barrel);
    SignalElectron_RelIsoRhoCorrected_Barrel.push_back(tmpSignalElectron_RelIsoRhoCorrected_Barrel);
    SignalElectron_RelIso_Endcap.push_back(tmpSignalElectron_RelIso_Endcap);
    SignalElectron_RelIsoRhoCorrected_Endcap.push_back(tmpSignalElectron_RelIsoRhoCorrected_Endcap);
    SignalElectron_RelIso04_Barrel.push_back(tmpSignalElectron_RelIso04_Barrel);
    SignalElectron_RelIso04RhoCorrected_Barrel.push_back(tmpSignalElectron_RelIso04RhoCorrected_Barrel);
    SignalElectron_RelIso04_Endcap.push_back(tmpSignalElectron_RelIso04_Endcap);
    SignalElectron_RelIso04RhoCorrected_Endcap.push_back(tmpSignalElectron_RelIso04RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap);
    SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpSignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);






    BkgElectron_RelIso_Barrel.push_back(tmpBkgElectron_RelIso_Barrel);
    BkgElectron_RelIsoRhoCorrected_Barrel.push_back(tmpBkgElectron_RelIsoRhoCorrected_Barrel);
    BkgElectron_RelIso_Endcap.push_back(tmpBkgElectron_RelIso_Endcap);
    BkgElectron_RelIsoRhoCorrected_Endcap.push_back(tmpBkgElectron_RelIsoRhoCorrected_Endcap);
    BkgElectron_RelIso04_Barrel.push_back(tmpBkgElectron_RelIso04_Barrel);
    BkgElectron_RelIso04RhoCorrected_Barrel.push_back(tmpBkgElectron_RelIso04RhoCorrected_Barrel);
    BkgElectron_RelIso04_Endcap.push_back(tmpBkgElectron_RelIso04_Endcap);
    BkgElectron_RelIso04RhoCorrected_Endcap.push_back(tmpBkgElectron_RelIso04RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel);
    BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap);
    BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);
    BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap.push_back(tmpBkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap);


  }





    TH1F *NVertex0To4_SignalElectron_RelIso_Barrel = (TH1F*)SignalElectron_RelIso_Barrel[0]->Clone((string("Electron_RelIso_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIsoRhoCorrected_Barrel = (TH1F*)SignalElectron_RelIsoRhoCorrected_Barrel[0]->Clone((string("Electron_RelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIso_Endcap = (TH1F*)SignalElectron_RelIso_Endcap[0]->Clone((string("Electron_RelIso_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIsoRhoCorrected_Endcap = (TH1F*)SignalElectron_RelIsoRhoCorrected_Endcap[0]->Clone((string("Electron_RelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIso04_Barrel = (TH1F*)SignalElectron_RelIso04_Barrel[0]->Clone((string("Electron_RelIso04_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIso04RhoCorrected_Barrel = (TH1F*)SignalElectron_RelIso04RhoCorrected_Barrel[0]->Clone((string("Electron_RelIso04RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIso04_Endcap = (TH1F*)SignalElectron_RelIso04_Endcap[0]->Clone((string("Electron_RelIso04_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_RelIso04RhoCorrected_Endcap = (TH1F*)SignalElectron_RelIso04RhoCorrected_Endcap[0]->Clone((string("Electron_RelIso04RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());


    TH1F *NVertex0To4_BkgElectron_RelIso_Barrel = (TH1F*)BkgElectron_RelIso_Barrel[0]->Clone((string("Electron_RelIso_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIsoRhoCorrected_Barrel = (TH1F*)BkgElectron_RelIsoRhoCorrected_Barrel[0]->Clone((string("Electron_RelIsoRhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIso_Endcap = (TH1F*)BkgElectron_RelIso_Endcap[0]->Clone((string("Electron_RelIso_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIsoRhoCorrected_Endcap = (TH1F*)BkgElectron_RelIsoRhoCorrected_Endcap[0]->Clone((string("Electron_RelIsoRhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIso04_Barrel = (TH1F*)BkgElectron_RelIso04_Barrel[0]->Clone((string("Electron_RelIso04_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIso04RhoCorrected_Barrel = (TH1F*)BkgElectron_RelIso04RhoCorrected_Barrel[0]->Clone((string("Electron_RelIso04RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIso04_Endcap = (TH1F*)BkgElectron_RelIso04_Endcap[0]->Clone((string("Electron_RelIso04_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_RelIso04RhoCorrected_Endcap = (TH1F*)BkgElectron_RelIso04RhoCorrected_Endcap[0]->Clone((string("Electron_RelIso04RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());
    TH1F *NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex0To4_") + signalSample).c_str());





    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso_Barrel, SignalElectron_RelIso_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIsoRhoCorrected_Barrel, SignalElectron_RelIsoRhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso_Endcap, SignalElectron_RelIso_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIsoRhoCorrected_Endcap, SignalElectron_RelIsoRhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso04_Barrel, SignalElectron_RelIso04_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso04RhoCorrected_Barrel, SignalElectron_RelIso04RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso04_Endcap, SignalElectron_RelIso04_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_RelIso04RhoCorrected_Endcap, SignalElectron_RelIso04RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);




    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso_Barrel, BkgElectron_RelIso_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIsoRhoCorrected_Barrel, BkgElectron_RelIsoRhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso_Endcap, BkgElectron_RelIso_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIsoRhoCorrected_Endcap, BkgElectron_RelIsoRhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso04_Barrel, BkgElectron_RelIso04_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso04RhoCorrected_Barrel, BkgElectron_RelIso04RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso04_Endcap, BkgElectron_RelIso04_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_RelIso04RhoCorrected_Endcap, BkgElectron_RelIso04RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 0, 4);
    AddIsoHistograms(NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 0, 4);








    TH1F *NVertex7To15_SignalElectron_RelIso_Barrel = (TH1F*)SignalElectron_RelIso_Barrel[0]->Clone((string("Electron_RelIso_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIsoRhoCorrected_Barrel = (TH1F*)SignalElectron_RelIsoRhoCorrected_Barrel[0]->Clone((string("Electron_RelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIso_Endcap = (TH1F*)SignalElectron_RelIso_Endcap[0]->Clone((string("Electron_RelIso_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIsoRhoCorrected_Endcap = (TH1F*)SignalElectron_RelIsoRhoCorrected_Endcap[0]->Clone((string("Electron_RelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIso04_Barrel = (TH1F*)SignalElectron_RelIso04_Barrel[0]->Clone((string("Electron_RelIso04_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIso04RhoCorrected_Barrel = (TH1F*)SignalElectron_RelIso04RhoCorrected_Barrel[0]->Clone((string("Electron_RelIso04RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIso04_Endcap = (TH1F*)SignalElectron_RelIso04_Endcap[0]->Clone((string("Electron_RelIso04_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_RelIso04RhoCorrected_Endcap = (TH1F*)SignalElectron_RelIso04RhoCorrected_Endcap[0]->Clone((string("Electron_RelIso04RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());


    TH1F *NVertex7To15_BkgElectron_RelIso_Barrel = (TH1F*)BkgElectron_RelIso_Barrel[0]->Clone((string("Electron_RelIso_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIsoRhoCorrected_Barrel = (TH1F*)BkgElectron_RelIsoRhoCorrected_Barrel[0]->Clone((string("Electron_RelIsoRhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIso_Endcap = (TH1F*)BkgElectron_RelIso_Endcap[0]->Clone((string("Electron_RelIso_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIsoRhoCorrected_Endcap = (TH1F*)BkgElectron_RelIsoRhoCorrected_Endcap[0]->Clone((string("Electron_RelIsoRhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIso04_Barrel = (TH1F*)BkgElectron_RelIso04_Barrel[0]->Clone((string("Electron_RelIso04_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIso04RhoCorrected_Barrel = (TH1F*)BkgElectron_RelIso04RhoCorrected_Barrel[0]->Clone((string("Electron_RelIso04RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIso04_Endcap = (TH1F*)BkgElectron_RelIso04_Endcap[0]->Clone((string("Electron_RelIso04_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_RelIso04RhoCorrected_Endcap = (TH1F*)BkgElectron_RelIso04RhoCorrected_Endcap[0]->Clone((string("Electron_RelIso04RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());
    TH1F *NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap = (TH1F*)BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap[0]->Clone((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap_NVertex7To15_") + signalSample).c_str());





    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso_Barrel, SignalElectron_RelIso_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Barrel, SignalElectron_RelIsoRhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso_Endcap, SignalElectron_RelIso_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Endcap, SignalElectron_RelIsoRhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso04_Barrel, SignalElectron_RelIso04_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso04RhoCorrected_Barrel, SignalElectron_RelIso04RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso04_Endcap, SignalElectron_RelIso04_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_RelIso04RhoCorrected_Endcap, SignalElectron_RelIso04RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);




    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso_Barrel, BkgElectron_RelIso_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIsoRhoCorrected_Barrel, BkgElectron_RelIsoRhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso_Endcap, BkgElectron_RelIso_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIsoRhoCorrected_Endcap, BkgElectron_RelIsoRhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso04_Barrel, BkgElectron_RelIso04_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso04RhoCorrected_Barrel, BkgElectron_RelIso04RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso04_Endcap, BkgElectron_RelIso04_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_RelIso04RhoCorrected_Endcap, BkgElectron_RelIso04RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_05Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_05Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_10Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_10Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_TotalPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Barrel, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_TotalPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone04_15Threshold_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_TotalPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_FootprintRemovedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedPFRelIso_Cone04_15Threshold_RhoCorrected_Endcap, 7,15);
    AddIsoHistograms(NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_15Threshold_RhoCorrected_Endcap, 7,15);



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
  plotname = "ElectronBarrel_NVertex0To4_PFConesize";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Barrel, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Barrel, NVertex0To4_BkgElectron_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone04_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone04_01Threshold");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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
  plotname = "ElectronBarrel_NVertex0To4";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Barrel, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Barrel, NVertex0To4_BkgElectron_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Barrel, "SigEffVsCut_Electron_RelIso_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Barrel, NVertex0To4_BkgElectron_RelIso_Barrel, "SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso_Cone03_01Threshold");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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
  plotname = "ElectronBarrel_NVertex7To15";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Barrel, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Barrel_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Barrel, NVertex7To15_BkgElectron_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Barrel, "SigEffVsCut_Electron_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Barrel, NVertex7To15_BkgElectron_RelIso_Barrel, "SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso_Cone03_01Threshold");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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
  plotname = "ElectronEndcap_NVertex0To4_PFConesize";

//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Endcap, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Endcap, NVertex0To4_BkgElectron_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, NVertex0To4_BkgElectron_TotalPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone04_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone04_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone04_01Threshold");



 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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
  plotname = "ElectronEndcap_NVertex0To4";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Endcap, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Endcap, NVertex0To4_BkgElectron_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Endcap, "SigEffVsCut_Electron_RelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Endcap, NVertex0To4_BkgElectron_RelIso_Endcap, "SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso_Cone03_01Threshold");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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
  plotname = "ElectronEndcap_NVertex7To15";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Endcap, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex7To15", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Endcap, NVertex7To15_BkgElectron_RelIso_Endcap, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex7To15", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Endcap, "SigEffVsCut_Electron_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Endcap, NVertex7To15_BkgElectron_RelIso_Endcap, "SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, NVertex7To15_BkgElectron_TotalPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso_Cone03_01Threshold");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso_Cone03_01Threshold");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Cone03_01Threshold_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso_Cone03_01Threshold");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=1; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
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
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
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









  canvasFile->Close();
  return;






  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronEndcap_RhoCorrection_NVertex7To15";
//*************************************************************************************
  tmpSigEffVsCut = MakeCurrentWPSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Barrel, "CurrentWP_SigEffVsCut_Electron_TotalPFRelIso_Barrel_NVertex0To4", currentWPCut);
  tmpSigEffVsBkgEff = MakeCurrentWPSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Barrel, NVertex0To4_BkgElectron_RelIso_Barrel, "CurrentWP_SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex0To4", currentWPCut);
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Current WP");

  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Barrel, "SigEffVsCut_Electron_RelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Barrel, NVertex7To15_BkgElectron_RelIso_Barrel, "SigEffVsBkgEff_Electron_RelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Barrel, "SigEffVsCut_Electron_TotalPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Barrel, NVertex7To15_BkgElectron_TotalPFRelIso_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Barrel, "SigEffVsCut_Electron_RelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Barrel, NVertex7To15_BkgElectron_RelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Electron_RelIsoRhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Electron_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIsoRhoCorrected_Barrel, NVertex7To15_BkgElectron_TotalPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Electron_TotalPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
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
  plotname = "ElectronBarrel_VertexSelectedRhoCorrection_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Barrel, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Barrel, NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel, NVertex7To15_BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Footprint Removed PFNoPU");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Footprint Removed PFNoPU Corr");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
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
  plotname = "ElectronEndcap_NVertex0To4";
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_RelIso_Endcap, "SigEffVsCut_Electron_RelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_RelIso_Endcap, NVertex0To4_BkgElectron_RelIso_Endcap, "SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_TotalPFRelIso_Endcap, NVertex0To4_BkgElectron_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_FootprintRemovedPFRelIso_Endcap, NVertex0To4_BkgElectron_FootprintRemovedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedPFRelIso_Endcap, NVertex0To4_BkgElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex0To4");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex0To4_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, NVertex0To4_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex0To4");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
    




  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronEndcap_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Endcap, "SigEffVsCut_Electron_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Endcap, NVertex7To15_BkgElectron_RelIso_Endcap, "SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");

//*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Endcap, NVertex7To15_BkgElectron_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("TotalPFRelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Endcap, "SigEffVsCut_Electron_FootprintRemovedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_FootprintRemovedPFRelIso_Endcap, NVertex7To15_BkgElectron_FootprintRemovedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_FootprintRemovedPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFRelIso");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Endcap, NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPURelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("FootprintRemovedPFNoPURelIso");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());
    

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsCutGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsCutGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsCutGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsCutGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsCutGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsCutGraphs[i]->Draw("AP");
    } else {
      SignalEffVsCutGraphs[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
//   cv->SaveAs(("SigEffVsCut_" + plotname + label + ".gif").c_str());
    






  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  //*************************************************************************************
  SignalEffVsCutGraphs.clear();
  SignalEffVsBkgEffGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronEndcap_RhoCorrection_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIso_Endcap, "SigEffVsCut_Electron_RelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIso_Endcap, NVertex7To15_BkgElectron_RelIso_Endcap, "SigEffVsBkgEff_Electron_RelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Standard RelIso");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Endcap, "SigEffVsCut_Electron_TotalPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIso_Endcap, NVertex7To15_BkgElectron_TotalPFRelIso_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Endcap, "SigEffVsCut_Electron_RelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_RelIsoRhoCorrected_Endcap, NVertex7To15_BkgElectron_RelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Electron_RelIsoRhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("RelIso Corr");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Electron_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_TotalPFRelIsoRhoCorrected_Endcap, NVertex7To15_BkgElectron_TotalPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Electron_TotalPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFRelIso Corr");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
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
  plotname = "ElectronEndcap_VertexSelectedRhoCorrection_NVertex7To15";
  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIso_Endcap, NVertex7To15_BkgElectron_VertexSelectedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU");


 //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap, NVertex7To15_BkgElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("PFNoPU Corr");



  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Footprint Removed PFNoPU");


  //*************************************************************************************
  tmpSigEffVsCut = MakeSigEffVsCutValueGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap, "SigEffVsCut_Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  tmpSigEffVsBkgEff = MakeSigEffVsBkgEffGraph(NVertex7To15_SignalElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap, NVertex7To15_BkgElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap, "SigEffVsBkgEff_Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap_NVertex7To15");
  SignalEffVsCutGraphs.push_back(tmpSigEffVsCut);
  SignalEffVsBkgEffGraphs.push_back(tmpSigEffVsBkgEff);
  GraphLabels.push_back("Footprint Removed PFNoPU Corr");


 //*******************************************************************************************

  legend = new TLegend(0.43,0.20,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(SignalEffVsBkgEffGraphs[i],GraphLabels[i].c_str(), "LP");

    SignalEffVsBkgEffGraphs[i]->SetMarkerColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetLineColor(colors[i]);
    SignalEffVsBkgEffGraphs[i]->SetMarkerSize(0.5);
   
    SignalEffVsBkgEffGraphs[i]->GetYaxis()->SetRangeUser(0.5,1.0);    
    SignalEffVsBkgEffGraphs[i]->GetXaxis()->SetRangeUser(0,0.2);    
    if (i==0) {
      SignalEffVsBkgEffGraphs[i]->Draw("AP");
    } else {
      SignalEffVsBkgEffGraphs[i]->Draw("Psame");
    }
  }
  legend->Draw();
  
  cv->SaveAs(("SigEffVsBkgEff_" + plotname + label + ".gif").c_str());








  canvasFile->Close();
  return;






}




void PlotElectronIsolationPerformance() {



//   PlotElectronIsolationPerformance("HWW130SignalWJetsMCBkg_Pt10To20", "HWW130_Pt10To20", "WJetsMC_Pt10To20", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalWJetsMCBkg_Pt20To35", "HWW130_Pt20To35", "WJetsMC_Pt20To35", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalWJetsMCBkg_Pt35To50", "HWW130_Pt35To50", "WJetsMC_Pt35To50", 0.10);

//   PlotElectronIsolationPerformance("HWW130SignalQCD15To3000MCBkg_Pt10To20", "HWW130_Pt10To20", "QCD15To3000MC_Ele10Jet30_Pt10To20", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalQCD15To3000MCBkg_Pt20To35", "HWW130_Pt20To35", "QCD15To3000MC_Ele10Jet30_Pt20To35", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalQCD15To3000MCBkg_Pt35To50", "HWW130_Pt35To50", "QCD15To3000MC_Ele10Jet30_Pt35To50", 0.10);

//   PlotElectronIsolationPerformance("HWW130SignalQCD30To50MCBkg_Pt10To20", "HWW130_Pt10To20", "QCD30To50MC_Ele10Jet30_Pt10To20", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalQCD30To50MCBkg_Pt20To35", "HWW130_Pt20To35", "QCD30To50MC_Ele10Jet30_Pt20To35", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalQCD30To50MCBkg_Pt35To50", "HWW130_Pt35To50", "QCD30To50MC_Ele10Jet30_Pt35To50", 0.10);

  PlotElectronIsolationPerformance("HWW130SignalDataEle10Jet30Bkg_Pt10To20", "HWW130_Pt10To20", "Ele10Jet30_Pt10To20", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalDataEle10Jet30Bkg_Pt20To35", "HWW130_Pt20To35", "Ele10Jet30_Pt20To35", 0.10);
//   PlotElectronIsolationPerformance("HWW130SignalDataEle10Jet30Bkg_Pt35To50", "HWW130_Pt35To50", "Ele10Jet30_Pt35To50", 0.10);

 



//   PlotElectronIsolationPerformance("HWW130MCSignalWJetsMCBkg_Pt10To20", "HWW130_Pt10To20", "WJetsMC_Pt10To20");
//   PlotElectronIsolationPerformance("TPSignalWJetsMCBkg_Pt10To20", "DataTagAndProbe_Pt10To20", "WJetsMC_Pt10To20");

//   PlotElectronIsolationPerformance("ZMCSignalQCDMCEle10Jet30Bkg", "Zee", "QCDMC_Ele10Jet30");

//   PlotElectronIsolationPerformance("ZMCSignalDataEle10Jet30Bkg", "Zee", "Ele10Jet30");

//   PlotElectronIsolationPerformance("ZMCSignalDataEle10Bkg", "Zee", "Ele10");




//   PlotElectronIsoVsNVertices("DataTPSignalWJetsMCBkg", "DataTagAndProbe", "WJetsMC");
//   PlotElectronIsoVsNVertices("HWW130SignalWJetsMCBkg", "HWW130_Pt10To20", "WJetsMC_Pt10To20");

//    PlotElectronIsoVsNVertices("QCDMCEle10Jet30Bkg", "Zee", "QCDMC_Ele10Jet30");

//    PlotElectronIsoVsNVertices("DataEle10Jet30Bkg", "Zee", "Ele10Jet30");

//   PlotElectronIsoVsNVertices("DataEle10Bkg", "Zee", "Ele10");



}
