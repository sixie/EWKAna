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
  TH1F *MuIDBDTGV10_Real = new TH1F(("MuIDBDTGV10_Real"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *MuIDBDTGV10_Fake = new TH1F(("MuIDBDTGV10_Fake"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);

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

  Float_t                 fMuBDTGV10; 


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


  RealMuTree->SetBranchAddress( "BDTGV10", &fMuBDTGV10); 


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

    if (passCutBased(fMuPt, fMuEta, fMuPFIso, fMuTkNchi2, fMuGlobalNchi2, fMuNValidHits, fMuNTrackerHits,fMuNPixelHits,fMuNMatches,fMuD0)) 
      RealMuonPassCutBased += fWeight;
    if (passMVASameCutBasedSig( fMuPt, fMuEta, fMuBDTGV10))
      RealMuonPassBDTGV10_SameCutBasedSig += fWeight;
    if (passMVAHalfCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10)) {
      RealMuonPassBDTGV10_HalfCutBasedBkg += fWeight;
    }
    if (passMVAOneThirdCutBasedBkg( fMuPt, fMuEta, fMuBDTGV10))
      RealMuonPassBDTGV10_OneThirdCutBasedBkg += fWeight;


    MuIDBDTGV10_Real->Fill(fMuBDTGV10,fWeight);


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

  FakeMuTree->SetBranchAddress( "BDTGV10", &fMuBDTGV10); 


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

    MuIDBDTGV10_Fake->Fill(fMuBDTGV10,fWeight);


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
  Double_t SigEffBDTGV10_HalfBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, 0.5*BkgEffCutBased);
  Double_t SigEffBDTGV10_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, (1.0/3.0)*BkgEffCutBased);

  cout << "Signal Efficiency (wrt Cut-based) for : half bkg : one third bkg \n";
  cout << "BDTGV10 : " << SigEffBDTGV10_HalfBkg/SigEffCutBased << " : " << SigEffBDTGV10_OneThirdBkg/SigEffCutBased << endl;
 
  cout << "**********************\n";

  Double_t BkgEffBDTGV10_SameSig = FindBkgEffAtFixedSignalEfficiency(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, SigEffCutBased);

  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "BDTGV10 : " << BkgEffBDTGV10_SameSig/BkgEffCutBased << endl;

  cout << "**********************\n";




  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_BDTGV10 = MakeSigEffVsBkgEffGraph(MuIDBDTGV10_Real, MuIDBDTGV10_Fake, "ROC_BDTGV10"+label );


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

   ROCGraphs.push_back(ROC_BDTGV10);
   GraphLabels.push_back("BDT");
   colors.push_back(kRed);
  
  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;


  if (Option == 0 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.0; ymax = 1.0; }
  if (Option == 1 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.0; ymax = 1.0; }
  if (Option == 2 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.0; ymax = 1.0; }
  if (Option == 3 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.0; ymax = 1.0; }
  if (Option == 4 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.65; ymax = 1.0; }
  if (Option == 5 )                              { xmin = 0.0; xmax = 0.5; ymin = 0.65; ymax = 1.0; }

  cv = new TCanvas("cv", "cv", 800, 600);

  legend = new TLegend(0.34,0.14,0.90,0.44);
  legend->SetTextSize(0.07);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
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

  legend->AddEntry(ROC_CutBasedCurrentWP, "Cut-Based", "P");
  ROC_CutBasedCurrentWP->SetFillColor(kBlue);
  ROC_CutBasedCurrentWP->SetMarkerColor(kBlue);
  ROC_CutBasedCurrentWP->SetMarkerStyle(34);
  ROC_CutBasedCurrentWP->SetMarkerSize(2.5);
  ROC_CutBasedCurrentWP->Draw("Psame");

  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
  cv->SaveAs(("ROCGraphs_" + plotname + ".eps").c_str());


  gBenchmark->Show("WWTemplate");       
} 

