//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeMuonIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/MuonMVA/output/MuonNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/MuonMVA/output/MuonNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <MitStyle.h>
#include "TLegend.h"

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"
#include "EWKAna/Utils/LeptonIDCuts.hh"

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
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
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
    for (UInt_t q=b; q < nPoints+2; ++q) {
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
  tmpSigEffVsBkgEff->SetMarkerSize(0.5);
  tmpSigEffVsBkgEff->SetMarkerStyle(20);

  return tmpSigEffVsBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, string name ) {
  //Make Met Plots
  double SigEff[1];
  double BkgEff[1];
  double SigEffErrLow[1];
  double SigEffErrHigh[1];
  double BkgEffErrLow[1];
  double BkgEffErrHigh[1];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  double cutValue;

  SigEff[0] = signalEff;
  SigEffErrLow[0] = 0;
  SigEffErrHigh[0] = 0;
  BkgEff[0] = bkgEff;
  BkgEffErrLow[0] = 0;
  BkgEffErrHigh[0] = 0;

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
    for (UInt_t q=b; q < nPoints+2; ++q) {
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
    for (UInt_t q=b; q < nPoints+2; ++q) {
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

//   cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

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
Double_t FindCutValueAtFixedEfficiency(TH1F* hist, Double_t targetEff ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentEff = 0;
  const UInt_t nPoints = hist->GetXaxis()->GetNbins();
  double NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += hist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += hist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
//     cout << targetEff << " : " << ratio << " , " << Hist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetEff - ratio) < fabs(targetEff - bestCurrentEff)) {
      targetCutValue = hist->GetXaxis()->GetBinCenter(b);
      bestCurrentEff = ratio;
    }
  }

  return targetCutValue;
}





//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetSignalEff ) {
  //Make Met Plots


  Double_t targetBkgEff = 0;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      bestCurrentSignalEff = ratio;
      targetBkgEff = bkgEff;
    }
  }

  return targetBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
Double_t FindSigEffAtFixedBkgEfficiency(TH1F* signalHist, TH1F* bkgHist, Double_t targetBkgEff ) {
  //Make Met Plots

  Double_t targetSignalEff = 0;
  Double_t bestCurrentBkgEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;
  double NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t nbkg = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nbkg += bkgHist->GetBinContent(q);
    }

    Double_t sigEff = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetBkgEff - bkgEff) < fabs(targetBkgEff - bestCurrentBkgEff)) {
      bestCurrentBkgEff = bkgEff;
      targetSignalEff = sigEff;
    }
  }

  return targetSignalEff;
}

Bool_t passCutBasedIsoOnly(Double_t fMuPt, 
                    Double_t fMuEta, 
                    Double_t fMuPFIso ) {

  Bool_t pass = kTRUE;
  
  if (fMuPt < 10) pass = kFALSE;
  if (fabs(fMuEta) > 2.4) pass = kFALSE;

  Double_t iso = fMuPFIso;
  Double_t isoCutValue = 0;

  if (fabs(fMuEta) < 1.479) {
    if (fMuPt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (fMuPt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  } 

  if (! (fMuPFIso / fMuPt < isoCutValue)
    )
    pass = kFALSE;

  return pass;
}



Bool_t passCutBased(Double_t fMuPt, 
                    Double_t fMuEta, 
                    Double_t fMuPFIso,   
                    Double_t fMuTkNchi2,
                    Double_t fMuGlobalNchi2,
                    Double_t fMuNValidHits,
                    Double_t fMuNTrackerHits,
                    Double_t fMuNPixelHits,
                    Double_t fMuNMatches,
                    Double_t fMuD0 ) {

  Bool_t pass = kTRUE;
  
  if (fMuPt < 10) pass = kFALSE;
  if (fabs(fMuEta) > 2.4) pass = kFALSE;

  Double_t iso = fMuPFIso;
  Double_t isoCutValue = 0;

  if (fabs(fMuEta) < 1.479) {
    if (fMuPt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (fMuPt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  } 

  if (! 
      ( ( fMuGlobalNchi2 < 10.0
        && (fMuNValidHits > 0)
        && (fMuNMatches > 1 )
        )
        && fMuNTrackerHits > 10
        && (fMuNPixelHits > 0)
        && fabs(fMuD0) < 0.02
        && fMuPFIso / fMuPt < isoCutValue
        )
    )
    pass = kFALSE;


  if (fMuPt < 20) {
    if (!
        ( fabs(fMuD0) < 0.01
          )
      ) {
      pass = kFALSE;
    }    
  }
  return pass;
}



Bool_t passMVASameCutBasedSig( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 14.5) ptBin = 1;
  if (fMuPt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  //Cut Values for MVA-V8 (PFIso)
//   if (MVABin == 0) MVACut = -0.5514;
//   if (MVABin == 1) MVACut = -0.303;
//   if (MVABin == 2) MVACut = -0.4562;
//   if (MVABin == 3) MVACut = -0.269;
//   if (MVABin == 4) MVACut = 0.1726;
//   if (MVABin == 5) MVACut = 0.801;

  //Cut Values for MVA-V10 (Detector Based Iso)
  if (MVABin == 0) MVACut = -0.5618;
  if (MVABin == 1) MVACut = -0.3002;
  if (MVABin == 2) MVACut = -0.4642;
  if (MVABin == 3) MVACut = -0.2478;
  if (MVABin == 4) MVACut = 0.1706;
  if (MVABin == 5) MVACut = 0.8146;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passMVAHalfCutBasedBkg( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 14.5) ptBin = 1;
  if (fMuPt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.377;
  if (MVABin == 1) MVACut = -0.0902;
  if (MVABin == 2) MVACut = -0.221;
  if (MVABin == 3) MVACut = -0.0154;
  if (MVABin == 4) MVACut = 0.459;
  if (MVABin == 5) MVACut = 0.9158;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}


Bool_t passMVAOneThirdCutBasedBkg( Double_t fMuPt, Double_t fMuEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fMuEta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (fMuPt > 14.5) ptBin = 1;
  if (fMuPt > 20.0) ptBin = 2;
  
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = -0.2362;
  if (MVABin == 1) MVACut = 0.0058;
  if (MVABin == 2) MVACut = -0.0134;
  if (MVABin == 3) MVACut = 0.1526;
  if (MVABin == 4) MVACut = 0.7362;
  if (MVABin == 5) MVACut = 0.9474;

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}




//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeMuonIDMVAPerformancePlots(string RealMuonFile, string FakeMuonFile, string Label, Int_t Option, Int_t NVtxBin = -1)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *MuIDTMVALikelihood_Real = new TH1F(("MuIDTMVALikelihood_Real"+label).c_str(), "; TMVA Likelihood ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDKNN_Real = new TH1F(("MuIDKNN_Real"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLP_Real = new TH1F(("MuIDMLP_Real"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLPBNN_Real = new TH1F(("MuIDMLPBNN_Real"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDT_Real = new TH1F(("MuIDBDT_Real"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTG_Real = new TH1F(("MuIDBDTG_Real"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDCombinedMVA_Real = new TH1F(("MuIDCombinedMVA_Real"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDKNN_C_Real = new TH1F(("MuIDKNN_C_Real"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLP_C_Real = new TH1F(("MuIDMLP_C_Real"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLPBNN_C_Real = new TH1F(("MuIDMLPBNN_C_Real"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDT_C_Real = new TH1F(("MuIDBDT_C_Real"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTG_C_Real = new TH1F(("MuIDBDTG_C_Real"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDCombinedMVA_C_Real = new TH1F(("MuIDCombinedMVA_C_Real"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);

  TH1F *MuIDTMVALikelihood_Fake = new TH1F(("MuIDTMVALikelihood_Fake"+label).c_str(), "; TMVA Likelihood ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDKNN_Fake = new TH1F(("MuIDKNN_Fake"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLP_Fake = new TH1F(("MuIDMLP_Fake"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLPBNN_Fake = new TH1F(("MuIDMLPBNN_Fake"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDT_Fake = new TH1F(("MuIDBDT_Fake"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTG_Fake = new TH1F(("MuIDBDTG_Fake"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDCombinedMVA_Fake = new TH1F(("MuIDCombinedMVA_Fake"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDKNN_C_Fake = new TH1F(("MuIDKNN_C_Fake"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLP_C_Fake = new TH1F(("MuIDMLP_C_Fake"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDMLPBNN_C_Fake = new TH1F(("MuIDMLPBNN_C_Fake"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDT_C_Fake = new TH1F(("MuIDBDT_C_Fake"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTG_C_Fake = new TH1F(("MuIDBDTG_C_Fake"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDCombinedMVA_C_Fake = new TH1F(("MuIDCombinedMVA_C_Fake"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);


  TH1F *MuIDBDTGV1_Real = new TH1F(("MuIDBDTGV1_Real"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTV2_Real = new TH1F(("MuIDBDTV2_Real"+label).c_str(), "; BDT V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV2_Real = new TH1F(("MuIDBDTGV2_Real"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDNNV2_Real = new TH1F(("MuIDNNV2_Real"+label).c_str(), "; NN V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodV2_Real = new TH1F(("MuIDTMVALikelihoodV2_Real"+label).c_str(), "; TMVA Likelihood V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodDV2_Real = new TH1F(("MuIDTMVALikelihoodDV2_Real"+label).c_str(), "; TMVA LikelihoodD V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV1_Fake = new TH1F(("MuIDBDTGV1_Fake"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTV2_Fake = new TH1F(("MuIDBDTV2_Fake"+label).c_str(), "; BDT V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV2_Fake = new TH1F(("MuIDBDTGV2_Fake"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDNNV2_Fake = new TH1F(("MuIDNNV2_Fake"+label).c_str(), "; NN V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodV2_Fake = new TH1F(("MuIDTMVALikelihoodV2_Fake"+label).c_str(), "; TMVA Likelihood V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodDV2_Fake = new TH1F(("MuIDTMVALikelihoodDV2_Fake"+label).c_str(), "; TMVA LikelihoodD V2 ; Number of Events ",  10000, -2 , 2);

  TH1F *MuIDBDTV3_Real = new TH1F(("MuIDBDTV3_Real"+label).c_str(), "; BDT V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV3_Real = new TH1F(("MuIDBDTGV3_Real"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDNNV3_Real = new TH1F(("MuIDNNV3_Real"+label).c_str(), "; NN V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodV3_Real = new TH1F(("MuIDTMVALikelihoodV3_Real"+label).c_str(), "; TMVA Likelihood V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodDV3_Real = new TH1F(("MuIDTMVALikelihoodDV3_Real"+label).c_str(), "; TMVA LikelihoodD V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTV3_Fake = new TH1F(("MuIDBDTV3_Fake"+label).c_str(), "; BDT V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV3_Fake = new TH1F(("MuIDBDTGV3_Fake"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDNNV3_Fake = new TH1F(("MuIDNNV3_Fake"+label).c_str(), "; NN V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodV3_Fake = new TH1F(("MuIDTMVALikelihoodV3_Fake"+label).c_str(), "; TMVA Likelihood V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDTMVALikelihoodDV3_Fake = new TH1F(("MuIDTMVALikelihoodDV3_Fake"+label).c_str(), "; TMVA LikelihoodD V3 ; Number of Events ",  10000, -2 , 2);

  TH1F *MuIDBDTGV4_Real = new TH1F(("MuIDBDTGV4_Real"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV4_Fake = new TH1F(("MuIDBDTGV4_Fake"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV5_Real = new TH1F(("MuIDBDTGV5_Real"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV5_Fake = new TH1F(("MuIDBDTGV5_Fake"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV6_Real = new TH1F(("MuIDBDTGV6_Real"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV6_Fake = new TH1F(("MuIDBDTGV6_Fake"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV7_Real = new TH1F(("MuIDBDTGV7_Real"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV7_Fake = new TH1F(("MuIDBDTGV7_Fake"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV8_Real = new TH1F(("MuIDBDTGV8_Real"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV8_Fake = new TH1F(("MuIDBDTGV8_Fake"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV9_Real = new TH1F(("MuIDBDTGV9_Real"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV9_Fake = new TH1F(("MuIDBDTGV9_Fake"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV10_Real = new TH1F(("MuIDBDTGV10_Real"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV10_Fake = new TH1F(("MuIDBDTGV10_Fake"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV11_Real = new TH1F(("MuIDBDTGV11_Real"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV11_Fake = new TH1F(("MuIDBDTGV11_Fake"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);


  TH1F *BDTGSignalEfficiencyNumerator_Pt = new TH1F("BDTGSignalEfficiencyNumerator_Pt" , "; p_{T} [GeV^{c^}];Number of Events",  100, 0, 100);
  TH1F *BDTGSignalEfficiencyNumerator_Eta = new TH1F("BDTGSignalEfficiencyNumerator_Eta" , "; #eta;Number of Events",  100, -2.5, 2.5);
  TH1F *BDTGSignalEfficiencyNumerator_Phi = new TH1F("BDTGSignalEfficiencyNumerator_Phi" , "; #phi;Number of Events",  100, -3.2, 3.2);
  TH1F *BDTGSignalEfficiencyDenominator_Pt = new TH1F("BDTGSignalEfficiencyDenominator_Pt" , "; p_{T} [GeV^{c^}];Number of Events",  100, 0, 100);
  TH1F *BDTGSignalEfficiencyDenominator_Eta = new TH1F("BDTGSignalEfficiencyDenominator_Eta" , "; #eta;Number of Events",  100, -2.5, 2.5);
  TH1F *BDTGSignalEfficiencyDenominator_Phi = new TH1F("BDTGSignalEfficiencyDenominator_Phi" , "; #phi;Number of Events",  100, -3.2, 3.2);


  Double_t RealMuons = 0;
  Double_t FakeMuons = 0;
  Double_t RealMuonPassCutBased = 0;
  Double_t FakeMuonPassCutBased = 0;

  Double_t RealMuonPassBDTGV10_SameCutBasedSig = 0;
  Double_t FakeMuonPassBDTGV10_SameCutBasedSig = 0;
  Double_t RealMuonPassBDTGV10_HalfCutBasedBkg = 0;
  Double_t FakeMuonPassBDTGV10_HalfCutBasedBkg = 0;
  Double_t RealMuonPassBDTGV10_OneThirdCutBasedBkg = 0;
  Double_t FakeMuonPassBDTGV10_OneThirdCutBasedBkg = 0;



  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 



  Float_t                 fMuLikelihood; 
  Float_t                 fMuKNN; 
  Float_t                 fMuMLP; 
  Float_t                 fMuMLPBNN; 
  Float_t                 fMuBDTG; 
  Float_t                 fMuBDT; 
  Float_t                 fMuCombinedMVA; 
  Float_t                 fMuKNN_C; 
  Float_t                 fMuMLP_C; 
  Float_t                 fMuMLPBNN_C; 
  Float_t                 fMuBDTG_C; 
  Float_t                 fMuBDT_C; 
  Float_t                 fMuCombinedMVA_C; 

  Float_t                 fMuBDTGV1; 
  Float_t                 fMuBDTV2; 
  Float_t                 fMuBDTGV2; 
  Float_t                 fMuNNV2; 
  Float_t                 fMuLikelihoodV2; 
  Float_t                 fMuLikelihoodDV2; 
  Float_t                 fMuBDTV3; 
  Float_t                 fMuBDTGV3; 
  Float_t                 fMuNNV3; 
  Float_t                 fMuLikelihoodV3; 
  Float_t                 fMuLikelihoodDV3; 

  Float_t                 fMuBDTGV4; 
  Float_t                 fMuBDTGV5; 
  Float_t                 fMuBDTGV6; 
  Float_t                 fMuBDTGV7; 
  Float_t                 fMuBDTGV8; 
  Float_t                 fMuBDTGV9; 
  Float_t                 fMuBDTGV10; 
  Float_t                 fMuBDTGV11; 


  //*****************************************************************************************
  //RealMuTree
  //*****************************************************************************************
  TFile *RealMuFile = new TFile(RealMuonFile.c_str(), "READ");
  TTree *RealMuTree = (TTree*)RealMuFile->Get("Muons");
  RealMuTree->SetBranchAddress( "weight", &fWeight);
  RealMuTree->SetBranchAddress( "run", &fRunNumber);
  RealMuTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealMuTree->SetBranchAddress( "event", &fEventNumber);
  RealMuTree->SetBranchAddress( "pt", &fMuPt); 
  RealMuTree->SetBranchAddress( "eta", &fMuEta); 
  RealMuTree->SetBranchAddress( "phi", &fMuPhi); 
  RealMuTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  RealMuTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  RealMuTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  RealMuTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  RealMuTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  RealMuTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  RealMuTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  RealMuTree->SetBranchAddress( "D0", &fMuD0); 
  RealMuTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  RealMuTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  RealMuTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  RealMuTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  RealMuTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  RealMuTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  RealMuTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  RealMuTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  RealMuTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  RealMuTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  RealMuTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  RealMuTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  RealMuTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  RealMuTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  RealMuTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  RealMuTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  RealMuTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  RealMuTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  RealMuTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  RealMuTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  RealMuTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  RealMuTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  RealMuTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  RealMuTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  RealMuTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  RealMuTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  RealMuTree->SetBranchAddress( "Rho", &fRho); 
  RealMuTree->SetBranchAddress( "NVertices", &fNVertices); 



  RealMuTree->SetBranchAddress( "LikelihoodV0", &fMuLikelihood); 
  RealMuTree->SetBranchAddress( "KNNV0", &fMuKNN); 
  RealMuTree->SetBranchAddress( "MLPV0", &fMuMLP); 
  RealMuTree->SetBranchAddress( "MLPBNNV0", &fMuMLPBNN); 
  RealMuTree->SetBranchAddress( "BDTGV0", &fMuBDTG); 
  RealMuTree->SetBranchAddress( "BDTV0", &fMuBDT); 
  RealMuTree->SetBranchAddress( "CombinedMVAV0", &fMuCombinedMVA); 
  RealMuTree->SetBranchAddress( "KNNV0_C", &fMuKNN_C); 
  RealMuTree->SetBranchAddress( "MLPV0_C", &fMuMLP_C); 
  RealMuTree->SetBranchAddress( "MLPBNNV0_C", &fMuMLPBNN_C); 
  RealMuTree->SetBranchAddress( "BDTGV0_C", &fMuBDTG_C); 
  RealMuTree->SetBranchAddress( "BDTV0_C", &fMuBDT_C); 
  RealMuTree->SetBranchAddress( "CombinedMVAV0_C", &fMuCombinedMVA_C); 

  RealMuTree->SetBranchAddress( "BDTGV1", &fMuBDTGV1); 
  RealMuTree->SetBranchAddress( "BDTV2", &fMuBDTV2); 
  RealMuTree->SetBranchAddress( "BDTGV2", &fMuBDTGV2); 
  RealMuTree->SetBranchAddress( "NNV2", &fMuNNV2); 
  RealMuTree->SetBranchAddress( "LikelihoodV2", &fMuLikelihoodV2); 
  RealMuTree->SetBranchAddress( "LikelihoodDV2", &fMuLikelihoodDV2); 
  RealMuTree->SetBranchAddress( "BDTV3", &fMuBDTV3); 
  RealMuTree->SetBranchAddress( "BDTGV3", &fMuBDTGV3); 
  RealMuTree->SetBranchAddress( "NNV3", &fMuNNV3); 
  RealMuTree->SetBranchAddress( "LikelihoodV3", &fMuLikelihoodV3); 
  RealMuTree->SetBranchAddress( "LikelihoodDV3", &fMuLikelihoodDV3); 

  RealMuTree->SetBranchAddress( "BDTGV4", &fMuBDTGV4); 
  RealMuTree->SetBranchAddress( "BDTGV5", &fMuBDTGV5); 
  RealMuTree->SetBranchAddress( "BDTGV6", &fMuBDTGV6); 
  RealMuTree->SetBranchAddress( "BDTGV7", &fMuBDTGV7); 
  RealMuTree->SetBranchAddress( "BDTGV8", &fMuBDTGV8); 
  RealMuTree->SetBranchAddress( "BDTGV9", &fMuBDTGV9); 
  RealMuTree->SetBranchAddress( "BDTGV10", &fMuBDTGV10); 
  RealMuTree->SetBranchAddress( "BDTGV11", &fMuBDTGV11); 


  for(UInt_t ientry=0; ientry < RealMuTree->GetEntries(); ientry++) {       	
    RealMuTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //don't evaluate performance using training events
    if (fEventNumber % 2 == 0) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fMuEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fMuPt > 14.5) ptBin = 1;
    if (fMuPt > 20.0) ptBin = 2;
//     if (fMuPt > 35.0) ptBin = 3;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (!passCuts) continue;    

    if (NVtxBin == 0) if (!(fNVertices >= 0 && fNVertices <= 5)) continue;
    if (NVtxBin == 1) if (!(fNVertices >= 6 && fNVertices <= 10)) continue;
    if (NVtxBin == 2) if (!(fNVertices >= 11 && fNVertices <= 15)) continue;
    if (NVtxBin == 3) if (!(fNVertices >= 16 && fNVertices <= 20)) continue;
    if (NVtxBin == 4) if (!(fNVertices >= 21 && fNVertices <= 25)) continue;

    //apply denominator iso cuts
    if( (fMuChargedIso03 + fMuNeutralIso03_10Threshold)/fMuPt > 0.4 ) continue;
    if( (fMuChargedIso03 + fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03, fMuEta))/fMuPt > 0.4 ) continue;
    if( (fMuTrkIso03 + fMuEMIso03 + fMuHadIso03 - rho*MuonEffectiveArea(kMuEMIso03, fMuEta) - rho*MuonEffectiveArea(kMuHadIso03, fMuEta))/fMuPt > 0.4 ) continue;


    
    RealMuons += fWeight;
    BDTGSignalEfficiencyDenominator_Pt->Fill(fMuPt,fWeight);
    BDTGSignalEfficiencyDenominator_Eta->Fill(fMuEta,fWeight);
    BDTGSignalEfficiencyDenominator_Phi->Fill(fMuPhi,fWeight);

    if (passCutBased(fMuPt, fMuEta, fMuPFIso, fMuTkNchi2, fMuGlobalNchi2, fMuNValidHits, fMuNTrackerHits,fMuNPixelHits,fMuNMatches,fMuD0)) 
      RealMuonPassCutBased += fWeight;
    if (passMVASameCutBasedSig( fMuPt, fMuEta, fMuBDTGV10))
      RealMuonPassBDTGV10_SameCutBasedSig += fWeight;
    if (passMVAHalfCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10)) {
      RealMuonPassBDTGV10_HalfCutBasedBkg += fWeight;
      BDTGSignalEfficiencyNumerator_Pt->Fill(fMuPt,fWeight);
      BDTGSignalEfficiencyNumerator_Eta->Fill(fMuEta,fWeight);
      BDTGSignalEfficiencyNumerator_Phi->Fill(fMuPhi,fWeight);      
    }
    if (passMVAOneThirdCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10))
      RealMuonPassBDTGV10_OneThirdCutBasedBkg += fWeight;


    //Fill Histograms
    if (!passCutBasedIsoOnly(fMuPt, fMuEta, fMuPFIso)) {
      //if muon doesn't pass nominal pfiso cut, artificially make mva value very small
      MuIDBDTG_Real->Fill(-1.99,fWeight);
      MuIDBDTGV1_Real->Fill(-1.99,fWeight);
      MuIDBDTGV2_Real->Fill(-1.99,fWeight);
      MuIDBDTGV3_Real->Fill(-1.99,fWeight);
      MuIDBDTGV4_Real->Fill(-1.99,fWeight);
    } else {
      MuIDBDTG_Real->Fill(fMuBDTG,fWeight);
      MuIDBDTGV1_Real->Fill(fMuBDTGV1,fWeight);
      MuIDBDTGV2_Real->Fill(fMuBDTGV2,fWeight);
      MuIDBDTGV3_Real->Fill(fMuBDTGV3,fWeight);
      MuIDBDTGV4_Real->Fill(fMuBDTGV4,fWeight);
    }

    MuIDBDTGV5_Real->Fill(fMuBDTGV5,fWeight);
    MuIDBDTGV6_Real->Fill(fMuBDTGV6,fWeight);
    MuIDBDTGV7_Real->Fill(fMuBDTGV7,fWeight);
    MuIDBDTGV8_Real->Fill(fMuBDTGV8,fWeight);
    MuIDBDTGV9_Real->Fill(fMuBDTGV9,fWeight);
    MuIDBDTGV10_Real->Fill(fMuBDTGV10,fWeight);
    MuIDBDTGV11_Real->Fill(fMuBDTGV11,fWeight);


  } 
  





  //*****************************************************************************************
  //FakeMuTree
  //*****************************************************************************************
  TFile *FakeMuFile = new TFile(FakeMuonFile.c_str(), "READ");
  TTree *FakeMuTree = (TTree*)FakeMuFile->Get("Muons");
  FakeMuTree->SetBranchAddress( "weight", &fWeight);
  FakeMuTree->SetBranchAddress( "run", &fRunNumber);
  FakeMuTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeMuTree->SetBranchAddress( "event", &fEventNumber);
  FakeMuTree->SetBranchAddress( "pt", &fMuPt); 
  FakeMuTree->SetBranchAddress( "eta", &fMuEta); 
  FakeMuTree->SetBranchAddress( "phi", &fMuPhi); 
  FakeMuTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  FakeMuTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  FakeMuTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  FakeMuTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  FakeMuTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  FakeMuTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  FakeMuTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  FakeMuTree->SetBranchAddress( "D0", &fMuD0); 
  FakeMuTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  FakeMuTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  FakeMuTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  FakeMuTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  FakeMuTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  FakeMuTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  FakeMuTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  FakeMuTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  FakeMuTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  FakeMuTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  FakeMuTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  FakeMuTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  FakeMuTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  FakeMuTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  FakeMuTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  FakeMuTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  FakeMuTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  FakeMuTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  FakeMuTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  FakeMuTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  FakeMuTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  FakeMuTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  FakeMuTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  FakeMuTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  FakeMuTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  FakeMuTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  FakeMuTree->SetBranchAddress( "Rho", &fRho); 
  FakeMuTree->SetBranchAddress( "NVertices", &fNVertices); 


  FakeMuTree->SetBranchAddress( "LikelihoodV0", &fMuLikelihood); 
  FakeMuTree->SetBranchAddress( "KNNV0", &fMuKNN); 
  FakeMuTree->SetBranchAddress( "MLPV0", &fMuMLP); 
  FakeMuTree->SetBranchAddress( "MLPBNNV0", &fMuMLPBNN); 
  FakeMuTree->SetBranchAddress( "BDTGV0", &fMuBDTG); 
  FakeMuTree->SetBranchAddress( "BDTV0", &fMuBDT); 
  FakeMuTree->SetBranchAddress( "CombinedMVAV0", &fMuCombinedMVA); 
  FakeMuTree->SetBranchAddress( "KNNV0_C", &fMuKNN_C); 
  FakeMuTree->SetBranchAddress( "MLPV0_C", &fMuMLP_C); 
  FakeMuTree->SetBranchAddress( "MLPBNNV0_C", &fMuMLPBNN_C); 
  FakeMuTree->SetBranchAddress( "BDTGV0_C", &fMuBDTG_C); 
  FakeMuTree->SetBranchAddress( "BDTV0_C", &fMuBDT_C); 
  FakeMuTree->SetBranchAddress( "CombinedMVAV0_C", &fMuCombinedMVA_C); 

  FakeMuTree->SetBranchAddress( "BDTGV1", &fMuBDTGV1); 
  FakeMuTree->SetBranchAddress( "BDTV2", &fMuBDTV2); 
  FakeMuTree->SetBranchAddress( "BDTGV2", &fMuBDTGV2); 
  FakeMuTree->SetBranchAddress( "NNV2", &fMuNNV2); 
  FakeMuTree->SetBranchAddress( "LikelihoodV2", &fMuLikelihoodV2); 
  FakeMuTree->SetBranchAddress( "LikelihoodDV2", &fMuLikelihoodDV2); 
  FakeMuTree->SetBranchAddress( "BDTV3", &fMuBDTV3); 
  FakeMuTree->SetBranchAddress( "BDTGV3", &fMuBDTGV3); 
  FakeMuTree->SetBranchAddress( "NNV3", &fMuNNV3); 
  FakeMuTree->SetBranchAddress( "LikelihoodV3", &fMuLikelihoodV3); 
  FakeMuTree->SetBranchAddress( "LikelihoodDV3", &fMuLikelihoodDV3); 

  FakeMuTree->SetBranchAddress( "BDTGV4", &fMuBDTGV4); 
  FakeMuTree->SetBranchAddress( "BDTGV5", &fMuBDTGV5); 
  FakeMuTree->SetBranchAddress( "BDTGV6", &fMuBDTGV6); 
  FakeMuTree->SetBranchAddress( "BDTGV7", &fMuBDTGV7); 
  FakeMuTree->SetBranchAddress( "BDTGV8", &fMuBDTGV8); 
  FakeMuTree->SetBranchAddress( "BDTGV9", &fMuBDTGV9); 
  FakeMuTree->SetBranchAddress( "BDTGV10", &fMuBDTGV10); 
  FakeMuTree->SetBranchAddress( "BDTGV11", &fMuBDTGV11); 


  for(UInt_t ientry=0; ientry < FakeMuTree->GetEntries(); ientry++) {       	
    FakeMuTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //don't evaluate performance using training events
    if (fEventNumber % 2 == 0) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fMuEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fMuPt > 14.5) ptBin = 1;
    if (fMuPt > 20.0) ptBin = 2;
//     if (fMuPt > 35.0) ptBin = 3;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (Option == 5) passCuts = (subdet == 1 && ptBin == 2);
    if (!passCuts) continue;    

    if (NVtxBin == 0) if (!(fNVertices >= 0 && fNVertices <= 5)) continue;
    if (NVtxBin == 1) if (!(fNVertices >= 6 && fNVertices <= 10)) continue;
    if (NVtxBin == 2) if (!(fNVertices >= 11 && fNVertices <= 15)) continue;
    if (NVtxBin == 3) if (!(fNVertices >= 16 && fNVertices <= 20)) continue;
    if (NVtxBin == 4) if (!(fNVertices >= 21 && fNVertices <= 25)) continue;

    //apply denominator iso cuts
    if( (fMuChargedIso03 + fMuNeutralIso03_10Threshold)/fMuPt > 0.4 ) continue;
    if( (fMuChargedIso03 + fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03, fMuEta))/fMuPt > 0.4 ) continue;
    if( (fMuTrkIso03 + fMuEMIso03 + fMuHadIso03 - rho*MuonEffectiveArea(kMuEMIso03, fMuEta) - rho*MuonEffectiveArea(kMuHadIso03, fMuEta))/fMuPt > 0.4 ) continue;


    FakeMuons += fWeight;

    if (passCutBased(fMuPt, fMuEta, fMuPFIso, fMuTkNchi2, fMuGlobalNchi2, fMuNValidHits, fMuNTrackerHits,fMuNPixelHits,fMuNMatches,fMuD0))
      FakeMuonPassCutBased += fWeight;
    if (passMVASameCutBasedSig( fMuPt, fMuEta, fMuBDTGV10))
      FakeMuonPassBDTGV10_SameCutBasedSig += fWeight;
    if (passMVAHalfCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10))
      FakeMuonPassBDTGV10_HalfCutBasedBkg += fWeight;
    if (passMVAOneThirdCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10))
      FakeMuonPassBDTGV10_OneThirdCutBasedBkg += fWeight;

 
    //Fill Histograms
    if (!passCutBasedIsoOnly(fMuPt, fMuEta, fMuPFIso)) {
      //if muon doesn't pass nominal pfiso cut, artificially make mva value very small
      MuIDBDTG_Fake->Fill(-1.99,fWeight);
      MuIDBDTGV1_Fake->Fill(-1.99,fWeight);
      MuIDBDTGV2_Fake->Fill(-1.99,fWeight);
      MuIDBDTGV3_Fake->Fill(-1.99,fWeight);
      MuIDBDTGV4_Fake->Fill(-1.99,fWeight);
    } else {
      MuIDBDTG_Fake->Fill(fMuBDTG,fWeight);
      MuIDBDTGV1_Fake->Fill(fMuBDTGV1,fWeight);
      MuIDBDTGV2_Fake->Fill(fMuBDTGV2,fWeight);
      MuIDBDTGV3_Fake->Fill(fMuBDTGV3,fWeight);
      MuIDBDTGV4_Fake->Fill(fMuBDTGV4,fWeight);
    }

    MuIDBDTGV5_Fake->Fill(fMuBDTGV5,fWeight);
    MuIDBDTGV6_Fake->Fill(fMuBDTGV6,fWeight);
    MuIDBDTGV7_Fake->Fill(fMuBDTGV7,fWeight);
    MuIDBDTGV8_Fake->Fill(fMuBDTGV8,fWeight);
    MuIDBDTGV9_Fake->Fill(fMuBDTGV9,fWeight);
    MuIDBDTGV10_Fake->Fill(fMuBDTGV10,fWeight);
    MuIDBDTGV11_Fake->Fill(fMuBDTGV11,fWeight);


  } //loop over electrons
  



  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "Cut-Based Real Muon Efficiency : " << RealMuonPassCutBased << " / " << RealMuons << " = " << RealMuonPassCutBased/RealMuons << endl;
  cout << "Cut-Based Fake Muon Efficiency : " << FakeMuonPassCutBased << " / " << FakeMuons << " = " << FakeMuonPassCutBased/FakeMuons << endl;
  TGraphAsymmErrors* ROC_CutBasedCurrentWP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassCutBased/RealMuons , FakeMuonPassCutBased/FakeMuons, "ROC_CutBasedCurrentWP"+label);

  cout << "BDTGV10 SameCutBasedSig Real Muon Efficiency : " << RealMuonPassBDTGV10_SameCutBasedSig << " / " << RealMuons << " = " << RealMuonPassBDTGV10_SameCutBasedSig/RealMuons << endl;
  cout << "BDTGV10 SameCutBasedSig Fake Muon Efficiency : " << FakeMuonPassBDTGV10_SameCutBasedSig << " / " << FakeMuons << " = " << FakeMuonPassBDTGV10_SameCutBasedSig/FakeMuons << endl;
  TGraphAsymmErrors* ROC_BDTGV10_SameCutBasedSig_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassBDTGV10_SameCutBasedSig/RealMuons , FakeMuonPassBDTGV10_SameCutBasedSig/FakeMuons, "ROC_BDTGV10SameCutBasedSigWP"+label);

  cout << "BDTGV10 HalfCutBasedBkg Real Muon Efficiency : " << RealMuonPassBDTGV10_HalfCutBasedBkg << " / " << RealMuons << " = " << RealMuonPassBDTGV10_HalfCutBasedBkg/RealMuons << endl;
  cout << "BDTGV10 HalfCutBasedBkg Fake Muon Efficiency : " << FakeMuonPassBDTGV10_HalfCutBasedBkg << " / " << FakeMuons << " = " << FakeMuonPassBDTGV10_HalfCutBasedBkg/FakeMuons << endl;
  TGraphAsymmErrors* ROC_BDTGV10_HalfCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassBDTGV10_HalfCutBasedBkg/RealMuons , FakeMuonPassBDTGV10_HalfCutBasedBkg/FakeMuons, "ROC_BDTGV10HalfCutBasedBkgWP"+label);

  cout << "BDTGV10 OneThirdCutBasedBkg Real Muon Efficiency : " << RealMuonPassBDTGV10_OneThirdCutBasedBkg << " / " << RealMuons << " = " << RealMuonPassBDTGV10_OneThirdCutBasedBkg/RealMuons << endl;
  cout << "BDTGV10 OneThirdCutBasedBkg Fake Muon Efficiency : " << FakeMuonPassBDTGV10_OneThirdCutBasedBkg << " / " << FakeMuons << " = " << FakeMuonPassBDTGV10_OneThirdCutBasedBkg/FakeMuons << endl;
  TGraphAsymmErrors* ROC_BDTGV10_OneThirdCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealMuonPassBDTGV10_OneThirdCutBasedBkg/RealMuons , FakeMuonPassBDTGV10_OneThirdCutBasedBkg/FakeMuons, "ROC_BDTGV10OneThirdCutBasedBkgWP"+label);


  cout << "**********************\n";
  cout << "Bkg At LHTight Signal Eff\n";

  Double_t BkgEffCutBased = FakeMuonPassCutBased/FakeMuons;
  Double_t SigEffCutBased = RealMuonPassCutBased/RealMuons;
  Double_t SigEffBDTGV0_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTG_Real, MuIDBDTG_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV1_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV1_Real, MuIDBDTGV1_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV2_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV2_Real, MuIDBDTGV2_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV3_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV3_Real, MuIDBDTGV3_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV4_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV4_Real, MuIDBDTGV4_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV5_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV5_Real, MuIDBDTGV5_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV6_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV6_Real, MuIDBDTGV6_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV7_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV7_Real, MuIDBDTGV7_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV8_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV8_Real, MuIDBDTGV8_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV9_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV9_Real, MuIDBDTGV9_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV10_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV11_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV11_Real, MuIDBDTGV11_Fake, 0.5*BkgEffCutBased);

  Double_t SigEffBDTGV0_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTG_Real, MuIDBDTG_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV1_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV1_Real, MuIDBDTGV1_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV2_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV2_Real, MuIDBDTGV2_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV3_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV3_Real, MuIDBDTGV3_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV4_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV4_Real, MuIDBDTGV4_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV5_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV5_Real, MuIDBDTGV5_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV6_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV6_Real, MuIDBDTGV6_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV7_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV7_Real, MuIDBDTGV7_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV8_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV8_Real, MuIDBDTGV8_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV9_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV9_Real, MuIDBDTGV9_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV10_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, (1.0/3.0)*BkgEffCutBased);
  Double_t SigEffBDTGV11_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV11_Real, MuIDBDTGV11_Fake, (1.0/3.0)*BkgEffCutBased);

  cout << "Signal Efficiency (wrt Cut-based) for : half bkg : one third bkg \n";
  cout << "BDTGV0 : " << SigEffBDTGV0_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV0_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV1 : " << SigEffBDTGV1_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV1_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV2 : " << SigEffBDTGV2_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV2_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV3 : " << SigEffBDTGV3_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV3_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV4 : " << SigEffBDTGV4_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV4_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV5 : " << SigEffBDTGV5_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV5_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV6 : " << SigEffBDTGV6_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV6_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV7 : " << SigEffBDTGV7_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV7_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV8 : " << SigEffBDTGV8_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV8_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV9 : " << SigEffBDTGV9_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV9_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV10 : " << SigEffBDTGV10_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV10_OneThirdBkg/SigEffCutBased << endl;
  cout << "BDTGV11 : " << SigEffBDTGV11_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV11_OneThirdBkg/SigEffCutBased << endl;

  cout << "**********************\n";

  Double_t BkgEffBDTGV0_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTG_Real, MuIDBDTG_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV1_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV1_Real, MuIDBDTGV1_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV2_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV2_Real, MuIDBDTGV2_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV3_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV3_Real, MuIDBDTGV3_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV4_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV4_Real, MuIDBDTGV4_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV5_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV5_Real, MuIDBDTGV5_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV6_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV6_Real, MuIDBDTGV6_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV7_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV7_Real, MuIDBDTGV7_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV8_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV8_Real, MuIDBDTGV8_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV9_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV9_Real, MuIDBDTGV9_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV10_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, SigEffCutBased);
  Double_t BkgEffBDTGV11_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV11_Real, MuIDBDTGV11_Fake, SigEffCutBased);

  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "BDTGV0 : " << BkgEffBDTGV0_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV1 : " << BkgEffBDTGV1_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV2 : " << BkgEffBDTGV2_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV3 : " << BkgEffBDTGV3_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV4 : " << BkgEffBDTGV4_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV5 : " << BkgEffBDTGV5_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV6 : " << BkgEffBDTGV6_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV7 : " << BkgEffBDTGV7_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV8 : " << BkgEffBDTGV8_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV9 : " << BkgEffBDTGV9_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV10 : " << BkgEffBDTGV10_SameSig/BkgEffCutBased << endl;
  cout << "BDTGV11 : " << BkgEffBDTGV11_SameSig/BkgEffCutBased << endl;

  cout << "**********************\n";




  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
//   TGraphAsymmErrors* ROC_TMVALikelihood = MakeSigEffVsBkgEffGraph(MuIDTMVALikelihood_Real, MuIDTMVALikelihood_Fake, "ROC_TMVALikelihood"+label );
//   TGraphAsymmErrors* ROC_KNN = MakeSigEffVsBkgEffGraph(MuIDKNN_Real, MuIDKNN_Fake, "ROC_KNN"+label );
//   TGraphAsymmErrors* ROC_MLP = MakeSigEffVsBkgEffGraph(MuIDMLP_Real, MuIDMLP_Fake, "ROC_MLP"+label );
//   TGraphAsymmErrors* ROC_MLPBNN = MakeSigEffVsBkgEffGraph(MuIDMLPBNN_Real, MuIDMLPBNN_Fake, "ROC_MLPBNN"+label );
//   TGraphAsymmErrors* ROC_BDT = MakeSigEffVsBkgEffGraph(MuIDBDT_Real, MuIDBDT_Fake, "ROC_BDT"+label );
//   TGraphAsymmErrors* ROC_BDTG = MakeSigEffVsBkgEffGraph(MuIDBDTG_Real, MuIDBDTG_Fake, "ROC_BDTG"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA = MakeSigEffVsBkgEffGraph(MuIDCombinedMVA_Real, MuIDCombinedMVA_Fake, "ROC_CombinedMVA"+label );

//   TGraphAsymmErrors* ROC_KNN_C = MakeSigEffVsBkgEffGraph(MuIDKNN_C_Real, MuIDKNN_C_Fake, "ROC_KNN_C"+label );
//    TGraphAsymmErrors* ROC_MLP_C = MakeSigEffVsBkgEffGraph(MuIDMLP_C_Real, MuIDMLP_C_Fake, "ROC_MLP_C"+label );
//   TGraphAsymmErrors* ROC_MLPBNN_C = MakeSigEffVsBkgEffGraph(MuIDMLPBNN_C_Real, MuIDMLPBNN_C_Fake, "ROC_MLPBNN_C"+label );
//   TGraphAsymmErrors* ROC_BDT_C = MakeSigEffVsBkgEffGraph(MuIDBDT_C_Real, MuIDBDT_C_Fake, "ROC_BDT_C"+label );
//   TGraphAsymmErrors* ROC_BDTG_C = MakeSigEffVsBkgEffGraph(MuIDBDTG_C_Real, MuIDBDTG_C_Fake, "ROC_BDTG_C"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA_C = MakeSigEffVsBkgEffGraph(MuIDCombinedMVA_C_Real, MuIDCombinedMVA_C_Fake, "ROC_CombinedMVA_C"+label );


  TGraphAsymmErrors* ROC_BDTGV1 = MakeSigEffVsBkgEffGraph(MuIDBDTGV1_Real, MuIDBDTGV1_Fake, "ROC_BDTGV1"+label );
  TGraphAsymmErrors* ROC_BDTGV2 = MakeSigEffVsBkgEffGraph(MuIDBDTGV2_Real, MuIDBDTGV2_Fake, "ROC_BDTGV2"+label );
  TGraphAsymmErrors* ROC_BDTGV3 = MakeSigEffVsBkgEffGraph(MuIDBDTGV3_Real, MuIDBDTGV3_Fake, "ROC_BDTGV3"+label );
  TGraphAsymmErrors* ROC_BDTGV4 = MakeSigEffVsBkgEffGraph(MuIDBDTGV4_Real, MuIDBDTGV4_Fake, "ROC_BDTGV4"+label );
  TGraphAsymmErrors* ROC_BDTGV5 = MakeSigEffVsBkgEffGraph(MuIDBDTGV5_Real, MuIDBDTGV5_Fake, "ROC_BDTGV5"+label );
  TGraphAsymmErrors* ROC_BDTGV6 = MakeSigEffVsBkgEffGraph(MuIDBDTGV6_Real, MuIDBDTGV6_Fake, "ROC_BDTGV6"+label );
  TGraphAsymmErrors* ROC_BDTGV7 = MakeSigEffVsBkgEffGraph(MuIDBDTGV7_Real, MuIDBDTGV7_Fake, "ROC_BDTGV7"+label );
  TGraphAsymmErrors* ROC_BDTGV8 = MakeSigEffVsBkgEffGraph(MuIDBDTGV8_Real, MuIDBDTGV8_Fake, "ROC_BDTGV8"+label );
  TGraphAsymmErrors* ROC_BDTGV9 = MakeSigEffVsBkgEffGraph(MuIDBDTGV9_Real, MuIDBDTGV9_Fake, "ROC_BDTGV9"+label );
  TGraphAsymmErrors* ROC_BDTGV10 = MakeSigEffVsBkgEffGraph(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, "ROC_BDTGV10"+label );
  TGraphAsymmErrors* ROC_BDTGV11 = MakeSigEffVsBkgEffGraph(MuIDBDTGV11_Real, MuIDBDTGV11_Fake, "ROC_BDTGV11"+label );


  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_BDTGV10_SameSig = FindCutValueAtFixedEfficiency(MuIDBDTGV10_Real, SigEffCutBased );
  Double_t CutValue_BDTGV10_HalfBkg = FindCutValueAtFixedEfficiency(MuIDBDTGV10_Fake, 0.5*BkgEffCutBased );
  Double_t CutValue_BDTGV10_OneThirdBkg = FindCutValueAtFixedEfficiency(MuIDBDTGV10_Fake, (1.0/3.0)*BkgEffCutBased );
  cout << "BDTG V10 Cut Value @ Same Cut-Based Sig: " << CutValue_BDTGV10_SameSig << endl;
  cout << "BDTG V10 Cut Value @ 50% Cut-Based Bkg: " << CutValue_BDTGV10_HalfBkg << endl;
  cout << "BDTG V10 Cut Value @ 33% Cut-Based Bkg: " << CutValue_BDTGV10_OneThirdBkg << endl;


//   TFile *canvasFile = new TFile("MuonIDMVAPerformancePlots.root","UPDATE");
  TLegend* legend;
  TCanvas* cv;
  string plotname;

//   //*****************************************************************************************
//   //Plot Distributions
//   //*****************************************************************************************


//   vector<TH1F*> hists;
//   vector<string> histLabels;




  //*****************************************************************************************
  //Plot ROC Curves
  //*****************************************************************************************
  vector<TGraphAsymmErrors*> ROCGraphs;
  vector<string> GraphLabels;
  vector<Int_t> colors;

//   vector<Int_t> markers;
//   vector<Int_t> colors;
//   colors.push_back(kRed);     markers.push_back(20);
//   colors.push_back(kCyan);    markers.push_back(21);
// //   colors.push_back(kBlue);    markers.push_back(21);
//   colors.push_back(kMagenta); markers.push_back(22);
//   colors.push_back(kCyan);    markers.push_back(34);
//   colors.push_back(kBlack);   markers.push_back(29);
//   colors.push_back(kGreen);   markers.push_back(33);
//   colors.push_back(kRed-2);   markers.push_back(33);
//   colors.push_back(kOrange);   markers.push_back(33);
//   colors.push_back(kBlue-2);   markers.push_back(33);
//   colors.push_back(kGreen-2);   markers.push_back(33);
//   colors.push_back(kMagenta-2);   markers.push_back(33);

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "MuonIDMVA"+label;

  ROCGraphs.push_back(ROC_BDTGV1);
  GraphLabels.push_back("Baseline MVA");
  colors.push_back(kGreen+2);
  
  
//   ROCGraphs.push_back(ROC_BDTGV2);
//   GraphLabels.push_back("V2");
//   colors.push_back(kCyan+2);
  
  
//   ROCGraphs.push_back(ROC_BDTGV3);
//   GraphLabels.push_back("Baseline + SC,CC");
//   colors.push_back(kRed);
  
  
//   ROCGraphs.push_back(ROC_BDTGV4);
//   GraphLabels.push_back("Baseline + SC,CC,CaloEnergy");
//   colors.push_back(kBlue);


//   ROCGraphs.push_back(ROC_BDTGV5);
//   GraphLabels.push_back("V5");
//   colors.push_back(kBlack);

//   ROCGraphs.push_back(ROC_BDTGV6);
//   GraphLabels.push_back("V6");
//   colors.push_back(kOrange);

//    ROCGraphs.push_back(ROC_BDTGV7);
//    GraphLabels.push_back("V7");
//    colors.push_back(kCyan+2);

//    ROCGraphs.push_back(ROC_BDTGV8);
//    GraphLabels.push_back("V8");
//    colors.push_back(kBlue);

//    ROCGraphs.push_back(ROC_BDTGV9);
//    GraphLabels.push_back("ID+Iso03");
//    colors.push_back(kMagenta);

   ROCGraphs.push_back(ROC_BDTGV10);
   GraphLabels.push_back("ID+Iso03+Iso05");
   colors.push_back(kRed+2);
  
//    ROCGraphs.push_back(ROC_BDTGV11);
//    GraphLabels.push_back("ID+Iso03+Iso05+NVtx");
//    colors.push_back(kGreen+2);


  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;
// //   if (Option == 0 )                              { xmin = 0.15; xmax = 0.45; ymin = 0.75; ymax = 0.85; }

//FOr Data
   if (Option == 0 )                              { xmin = 0.05; xmax = 0.26; ymin = 0.35; ymax = 0.85; }
//   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//   if (Option == 2 )                              { xmin = 0.00; xmax = 0.40; ymin = 0.00; ymax = 1.00; }
//   if (Option == 3 )                              { xmin = 0.25; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
  if (Option == 4 )                                 { xmin = 0.02; xmax = 0.28; ymin = 0.75; ymax = 1.00; }
  if (Option == 5 )                                 { xmin = 0.15; xmax = 0.32; ymin = 0.75; ymax = 1.00; }



//   //For MC
//    if (Option == 0 )                              { xmin = 0.05; xmax = 0.26; ymin = 0.35; ymax = 0.85; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 2 )                              { xmin = 0.00; xmax = 0.40; ymin = 0.00; ymax = 1.00; }
// //   if (Option == 3 )                              { xmin = 0.25; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                                 { xmin = 0.02; xmax = 0.23; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                                 { xmin = 0.15; xmax = 0.32; ymin = 0.75; ymax = 1.00; }


//FOr Data
//   if (Option == 0 )                              { xmin = 0.10; xmax = 0.50; ymin = 0.70; ymax = 0.90; }
//   if (Option == 1 )                              { xmin = 0.05; xmax = 0.35; ymin = 0.30; ymax = 0.90; }
//   if (Option == 2 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.30; ymax = 0.80; }
//   if (Option == 3 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }

  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.54,0.14,0.90,0.44);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(0);
  for (UInt_t i=0; i<GraphLabels.size(); ++i) {
    legend->AddEntry(ROCGraphs[i],GraphLabels[i].c_str(), "LP");

    ROCGraphs[i]->SetMarkerColor(colors[i]);
    ROCGraphs[i]->SetLineColor(colors[i]);
    ROCGraphs[i]->SetMarkerSize(0.75);
   
    ROCGraphs[i]->GetXaxis()->SetRangeUser(xmin,xmax);    
    ROCGraphs[i]->GetYaxis()->SetRangeUser(ymin,ymax);    
    if (i==0) {
      ROCGraphs[i]->Draw("AP");
    } else {
      ROCGraphs[i]->Draw("Psame");
    }
  }

//   legend->AddEntry(ROC_LHTightWP, "NEW LH WP", "P");
//   ROC_LHTightWP->SetFillColor(kBlue);
//   ROC_LHTightWP->SetMarkerColor(kBlue);
//   ROC_LHTightWP->SetMarkerStyle(34);
//   ROC_LHTightWP->SetMarkerSize(2.5);
//   ROC_LHTightWP->Draw("Psame");
  
//   legend->AddEntry(ROC_BDTGV1WP, "V1 WP", "P");
//   ROC_BDTGV1WP->SetFillColor(kBlack);
//   ROC_BDTGV1WP->SetMarkerColor(kBlack);
//   ROC_BDTGV1WP->SetMarkerStyle(34);
//   ROC_BDTGV1WP->SetMarkerSize(2.5);
//   ROC_BDTGV1WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV2WP, "V2 WP", "P");
//   ROC_BDTGV2WP->SetFillColor(kRed);
//   ROC_BDTGV2WP->SetMarkerColor(kRed);
//   ROC_BDTGV2WP->SetMarkerStyle(34);
//   ROC_BDTGV2WP->SetMarkerSize(2.5);
//   ROC_BDTGV2WP->Draw("Psame");

  legend->AddEntry(ROC_CutBasedCurrentWP, "2011 Cut-Based", "P");
  ROC_CutBasedCurrentWP->SetFillColor(kGreen+3);
  ROC_CutBasedCurrentWP->SetMarkerColor(kGreen+3);
  ROC_CutBasedCurrentWP->SetMarkerStyle(34);
  ROC_CutBasedCurrentWP->SetMarkerSize(2.5);
  ROC_CutBasedCurrentWP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV8_SameCutBasedSig_WP, "BDT @ Same CutBased Sig", "P");
//   ROC_BDTGV8_SameCutBasedSig_WP->SetFillColor(kMagenta);
//   ROC_BDTGV8_SameCutBasedSig_WP->SetMarkerColor(kMagenta);
//   ROC_BDTGV8_SameCutBasedSig_WP->SetMarkerStyle(34);
//   ROC_BDTGV8_SameCutBasedSig_WP->SetMarkerSize(2.5);
//   ROC_BDTGV8_SameCutBasedSig_WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV8_HalfCutBasedBkg_WP, "BDT @ 50% CutBased Bkg", "P");
//   ROC_BDTGV8_HalfCutBasedBkg_WP->SetFillColor(kRed);
//   ROC_BDTGV8_HalfCutBasedBkg_WP->SetMarkerColor(kRed);
//   ROC_BDTGV8_HalfCutBasedBkg_WP->SetMarkerStyle(34);
//   ROC_BDTGV8_HalfCutBasedBkg_WP->SetMarkerSize(2.5);
//   ROC_BDTGV8_HalfCutBasedBkg_WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV8_OneThirdCutBasedBkg_WP, "BDT @ 33% CutBased Bkg", "P");
//   ROC_BDTGV8_OneThirdCutBasedBkg_WP->SetFillColor(kMagenta);
//   ROC_BDTGV8_OneThirdCutBasedBkg_WP->SetMarkerColor(kMagenta);
//   ROC_BDTGV8_OneThirdCutBasedBkg_WP->SetMarkerStyle(34);
//   ROC_BDTGV8_OneThirdCutBasedBkg_WP->SetMarkerSize(2.5);
//   ROC_BDTGV8_OneThirdCutBasedBkg_WP->Draw("Psame");

  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
//   canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;


  //*****************************************************************************************
  // Make Signal Efficiency Plots
  //*****************************************************************************************
  vector<double> ptbins;
  ptbins.push_back(10);  
  ptbins.push_back(15);  
  ptbins.push_back(20);  
  ptbins.push_back(25);  
  ptbins.push_back(30);  
  ptbins.push_back(35);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
  ptbins.push_back(100);  


  vector<double> etabins;
  for (UInt_t i=0; i< 40; ++i) {
    etabins.push_back(-3.0 + i*(6.0/40));
  }
  vector<double> phibins;
  for (UInt_t i=0; i< 20; ++i) {
    phibins.push_back(-3.2 + i*(6.4/20));
  }


  Int_t ErrorType = 2; //Clopper Pearson errors
  TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Pt, BDTGSignalEfficiencyDenominator_Pt, ("SignalEfficiency_Pt"+label).c_str(), ptbins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Eta, BDTGSignalEfficiencyDenominator_Eta, ("SignalEfficiency_Eta"+label).c_str(), etabins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Phi, BDTGSignalEfficiencyDenominator_Phi, ("SignalEfficiency_Phi"+label).c_str(), phibins, ErrorType, -99, -99, 0, 1);

//   cv = new TCanvas("cv", "cv", 800, 600);
//   efficiency_pt->Draw("AP");
//   cv->SaveAs(("SignalEfficiency"+label+"_Pt.gif").c_str());

//   cv = new TCanvas("cv", "cv", 800, 600);
//   efficiency_eta->Draw("AP");
//   cv->SaveAs(("SignalEfficiency"+label+"_Eta.gif").c_str());

//   cv = new TCanvas("cv", "cv", 800, 600);
//   efficiency_phi->Draw("AP");
//   cv->SaveAs(("SignalEfficiency"+label+"_Phi.gif").c_str());

  //*****************************************************************************************
//   canvasFile->Close();
  //*****************************************************************************************



//   //*****************************************************************************************
//   //Save Histograms in file
//   //*****************************************************************************************
  TFile *file = new TFile("MuonIDMVAResults.root", "UPDATE");
//   file->cd();
//   file->WriteTObject(MuIDBDT_Real, MuIDBDT_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDBDTGV3_Real, MuIDBDTGV3_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDNN_Real, MuIDNN_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDTMVALikelihood_Real, MuIDTMVALikelihood_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDTMVALikelihoodD_Real, MuIDTMVALikelihoodD_Real->GetName(), "WriteDelete");  

//   file->WriteTObject(MuIDBDT_Fake, MuIDBDT_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDBDTGV3_Fake, MuIDBDTGV3_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDNN_Fake, MuIDNN_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDTMVALikelihood_Fake, MuIDTMVALikelihood_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(MuIDTMVALikelihoodD_Fake, MuIDTMVALikelihoodD_Fake->GetName(), "WriteDelete"); 

//   file->WriteTObject(ROC_BDT, ROC_BDT->GetName(), "WriteDelete");  
//    file->WriteTObject(ROC_BDTGV3, ROC_BDTGV3->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_NN, ROC_NN->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihood, ROC_TMVALikelihood->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihoodD, ROC_TMVALikelihoodD->GetName(), "WriteDelete");  
 
   file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");  
   file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");  
   file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");  


   file->Close();
   delete file;




  gBenchmark->Show("WWTemplate");       
} 

