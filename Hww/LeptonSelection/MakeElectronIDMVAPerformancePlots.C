//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronIDMVAPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
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
Double_t FindCutValueAtFixedEfficiency(TH1F* signalHist, Double_t targetSignalEff ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  double NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    for (UInt_t q=b; q < nPoints+2; ++q) {
      nsig += signalHist->GetBinContent(q);
    }

    Double_t ratio = nsig / NSigTotal;
//     cout << targetSignalEff << " : " << ratio << " , " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      targetCutValue = signalHist->GetXaxis()->GetBinCenter(b);
      bestCurrentSignalEff = ratio;
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



Bool_t passVBTF( Int_t Option , Double_t fElePt, Double_t fEleSCEta , Double_t fEleEta, Double_t fEleSigmaIEtaIEta, Double_t fEleDEtaIn, Double_t fEleDPhiIn, Double_t fEleHoverE, Double_t fEleD0, Double_t fEleDZ, Double_t fEleFBrem, Double_t fEleEOverP) {
  
  Bool_t passNumerator = kTRUE;

  //VBTF80
  if (Option == 0) {    
    //Barrel 
    if (fabs(fEleSCEta) < 1.479) {
      if (! ( (0==0)
              && fEleSigmaIEtaIEta < 0.01 
              && fabs(fEleDEtaIn) < 0.004
              && fabs(fEleDPhiIn) < 0.06
              && fEleHoverE < 0.04
              && fabs(fEleD0) < 0.02
              && fabs(fEleDZ) < 0.1
            )
        ) {
        passNumerator = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && fEleSigmaIEtaIEta < 0.03
               && fabs(fEleDEtaIn) < 0.007
               && fabs(fEleDPhiIn) < 0.03
               && fEleHoverE < 0.025
               && fabs(fEleD0) < 0.02
               && fabs(fEleDZ) < 0.1
            )
        ) {
        passNumerator = kFALSE;
      }
    } 
  }

  //HWW LP2011 working point
  if (Option == 1) {    
    //Barrel 
    if (fabs(fEleSCEta) < 1.479) {
      if (! ( (0==0)
              && fEleSigmaIEtaIEta < 0.01 
              && fabs(fEleDEtaIn) < 0.004
              && fabs(fEleDPhiIn) < 0.06
              && fEleHoverE < 0.04
              && fabs(fEleD0) < 0.02
              && fabs(fEleDZ) < 0.1
            )
        ) {
        passNumerator = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && fEleSigmaIEtaIEta < 0.03
               && fabs(fEleDEtaIn) < 0.007
               && fabs(fEleDPhiIn) < 0.03
               && fEleHoverE < 0.10
               && fabs(fEleD0) < 0.02
               && fabs(fEleDZ) < 0.1
            )
        ) {
        passNumerator = kFALSE;
      }
    } 

    if (fElePt < 20) {
      //Barrel 
      if (fabs(fEleSCEta) < 1.479) {
        if (! ( (0==0)
                && fabs(fEleDEtaIn) < 0.004
                && fabs(fEleDPhiIn) < 0.03
                && fEleHoverE < 0.025
              )
          ) {
          passNumerator = kFALSE;
        }      
      }
      //Endcap
      else  {
        if (! (  (0==0)
                 && fabs(fEleDEtaIn) < 0.005
                 && fabs(fEleDPhiIn) < 0.02
              )
          ) {
          passNumerator = kFALSE;
        }
      } 

      if (fEleFBrem <= 0.15) {
        if (fabs(fEleEta) > 1.0) {
          passNumerator = kFALSE;
        } else {
          if (!( fEleEOverP > 0.95 )) passNumerator = kFALSE;
        }
      }
    }
  }


  return passNumerator;
}


Bool_t passLH( Double_t fElePt, Double_t fEleSCEta, Double_t fEleNBrem, Double_t fEleStandardLikelihood) {
  Double_t likCut = -999;

  if(fElePt > 20) {
    if(fabs(fEleSCEta) < 1.479){
      if(fEleNBrem == 0)                   likCut = 3.5;
      else                                 likCut = 4.0;
    }
    else  {                                
      if(fEleNBrem == 0)                   likCut = 4.0;
      else                                 likCut = 4.0;
    }
  }
  else {
    if(fabs(fEleSCEta) < 1.479){
      if(fEleNBrem == 0) likCut =  4.0;
      else                                 likCut =  4.5;
    }
    else  {                                
      if(fEleNBrem == 0) likCut =  4.0;
      else                                 likCut =  4.0;
    }
  }
  if (fEleStandardLikelihood > likCut) return kTRUE;
  return kFALSE;
}

Bool_t pass2011MVA( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 1.0) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 20.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.139;
  if (MVABin == 1) MVACut = 0.525;
  if (MVABin == 2) MVACut = 0.543; 
  if (MVABin == 3) MVACut = 0.947;
  if (MVABin == 4) MVACut = 0.950;
  if (MVABin == 5) MVACut = 0.884;
  

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

Bool_t passMVASame2011Sig( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 1.0) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 20.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.4202;
  if (MVABin == 1) MVACut = 0.6206;
  if (MVABin == 2) MVACut = 0.619; 
  if (MVABin == 3) MVACut = 0.959;
  if (MVABin == 4) MVACut = 0.9586;
  if (MVABin == 5) MVACut = 0.9278;
  

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}

Bool_t passMVAHalf2011Bkg( Double_t fElePt, Double_t fEleSCEta, Double_t MVAValue) {

  Int_t subdet = 0;
  if (fabs(fEleSCEta) < 1.0) subdet = 0;
  else if (fabs(fEleSCEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (fElePt > 20.0) ptBin = 1;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;
  if (MVABin == 0) MVACut = 0.5842;
  if (MVABin == 1) MVACut = 0.6998;
  if (MVABin == 2) MVACut = 0.6818; 
  if (MVABin == 3) MVACut = 0.9818;
  if (MVABin == 4) MVACut = 0.9818;
  if (MVABin == 5) MVACut = 0.957;
  

  if (MVAValue > MVACut) return kTRUE;
  return kFALSE;
}



Bool_t passCutBasedIsoOnly(Double_t fElePt, 
                    Double_t fEleEta, 
                    Double_t fElePFIso ) {

  //apply full isolation cut
  Bool_t passIsoCuts = kFALSE;
  if (fElePt >= 10 && fElePt < 20) passIsoCuts = ( fElePFIso  < 0.09 ); 
  if (fElePt >= 20) passIsoCuts = ( fElePFIso  < 0.13 ); 

  return passIsoCuts;
}



//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeElectronIDMVAPerformancePlots(string RealElectronFile, string FakeElectronFile, string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************
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


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *EleIDStandardLikelihood_Real = new TH1F(("EleIDStandardLikelihood_Real"+label).c_str(), "; StandardLikelihood ; Number of Events ",  10000, -20 , 20);
  TH1F *EleIDTMVALikelihood_Real = new TH1F(("EleIDTMVALikelihood_Real"+label).c_str(), "; TMVA Likelihood ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDKNN_Real = new TH1F(("EleIDKNN_Real"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLP_Real = new TH1F(("EleIDMLP_Real"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLPBNN_Real = new TH1F(("EleIDMLPBNN_Real"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDT_Real = new TH1F(("EleIDBDT_Real"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTG_Real = new TH1F(("EleIDBDTG_Real"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDCombinedMVA_Real = new TH1F(("EleIDCombinedMVA_Real"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDKNN_C_Real = new TH1F(("EleIDKNN_C_Real"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLP_C_Real = new TH1F(("EleIDMLP_C_Real"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLPBNN_C_Real = new TH1F(("EleIDMLPBNN_C_Real"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDT_C_Real = new TH1F(("EleIDBDT_C_Real"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTG_C_Real = new TH1F(("EleIDBDTG_C_Real"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDCombinedMVA_C_Real = new TH1F(("EleIDCombinedMVA_C_Real"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);

  TH1F *EleIDStandardLikelihood_Fake = new TH1F(("EleIDStandardLikelihood_Fake"+label).c_str(), "; StandardLikelihood ; Number of Events ",  10000, -20 , 20);
  TH1F *EleIDTMVALikelihood_Fake = new TH1F(("EleIDTMVALikelihood_Fake"+label).c_str(), "; TMVA Likelihood ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDKNN_Fake = new TH1F(("EleIDKNN_Fake"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLP_Fake = new TH1F(("EleIDMLP_Fake"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLPBNN_Fake = new TH1F(("EleIDMLPBNN_Fake"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDT_Fake = new TH1F(("EleIDBDT_Fake"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTG_Fake = new TH1F(("EleIDBDTG_Fake"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDCombinedMVA_Fake = new TH1F(("EleIDCombinedMVA_Fake"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDKNN_C_Fake = new TH1F(("EleIDKNN_C_Fake"+label).c_str(), "; KNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLP_C_Fake = new TH1F(("EleIDMLP_C_Fake"+label).c_str(), "; MLP ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDMLPBNN_C_Fake = new TH1F(("EleIDMLPBNN_C_Fake"+label).c_str(), "; MLPBNN ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDT_C_Fake = new TH1F(("EleIDBDT_C_Fake"+label).c_str(), "; BDT ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTG_C_Fake = new TH1F(("EleIDBDTG_C_Fake"+label).c_str(), "; BDTG ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDCombinedMVA_C_Fake = new TH1F(("EleIDCombinedMVA_C_Fake"+label).c_str(), "; CombinedMVA ; Number of Events ",  10000, -2 , 2);


  TH1F *EleIDBDTGV0_Real = new TH1F(("EleIDBDTGV0_Real"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV0_Fake = new TH1F(("EleIDBDTGV0_Fake"+label).c_str(), "; BDTG V0 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV1_Real = new TH1F(("EleIDBDTGV1_Real"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV1_Fake = new TH1F(("EleIDBDTGV1_Fake"+label).c_str(), "; BDTG V1 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV2_Real = new TH1F(("EleIDBDTGV2_Real"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV2_Fake = new TH1F(("EleIDBDTGV2_Fake"+label).c_str(), "; BDTG V2 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV3_Real = new TH1F(("EleIDBDTGV3_Real"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV3_Fake = new TH1F(("EleIDBDTGV3_Fake"+label).c_str(), "; BDTG V3 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV4_Real = new TH1F(("EleIDBDTGV4_Real"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV4_Fake = new TH1F(("EleIDBDTGV4_Fake"+label).c_str(), "; BDTG V4 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV5_Real = new TH1F(("EleIDBDTGV5_Real"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV5_Fake = new TH1F(("EleIDBDTGV5_Fake"+label).c_str(), "; BDTG V5 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV6_Real = new TH1F(("EleIDBDTGV6_Real"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV6_Fake = new TH1F(("EleIDBDTGV6_Fake"+label).c_str(), "; BDTG V6 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV7_Real = new TH1F(("EleIDBDTGV7_Real"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV7_Fake = new TH1F(("EleIDBDTGV7_Fake"+label).c_str(), "; BDTG V7 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV8_Real = new TH1F(("EleIDBDTGV8_Real"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV8_Fake = new TH1F(("EleIDBDTGV8_Fake"+label).c_str(), "; BDTG V8 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV9_Real = new TH1F(("EleIDBDTGV9_Real"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV9_Fake = new TH1F(("EleIDBDTGV9_Fake"+label).c_str(), "; BDTG V9 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV10_Real = new TH1F(("EleIDBDTGV10_Real"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV10_Fake = new TH1F(("EleIDBDTGV10_Fake"+label).c_str(), "; BDTG V10 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV11_Real = new TH1F(("EleIDBDTGV11_Real"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV11_Fake = new TH1F(("EleIDBDTGV11_Fake"+label).c_str(), "; BDTG V11 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV12_Real = new TH1F(("EleIDBDTGV12_Real"+label).c_str(), "; BDTG V12 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV12_Fake = new TH1F(("EleIDBDTGV12_Fake"+label).c_str(), "; BDTG V12 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV13_Real = new TH1F(("EleIDBDTGV13_Real"+label).c_str(), "; BDTG V13 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV13_Fake = new TH1F(("EleIDBDTGV13_Fake"+label).c_str(), "; BDTG V13 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV14_Real = new TH1F(("EleIDBDTGV14_Real"+label).c_str(), "; BDTG V14 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV14_Fake = new TH1F(("EleIDBDTGV14_Fake"+label).c_str(), "; BDTG V14 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV15_Real = new TH1F(("EleIDBDTGV15_Real"+label).c_str(), "; BDTG V15 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV15_Fake = new TH1F(("EleIDBDTGV15_Fake"+label).c_str(), "; BDTG V15 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV16_Real = new TH1F(("EleIDBDTGV16_Real"+label).c_str(), "; BDTG V16 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV16_Fake = new TH1F(("EleIDBDTGV16_Fake"+label).c_str(), "; BDTG V16 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV17_Real = new TH1F(("EleIDBDTGV17_Real"+label).c_str(), "; BDTG V17 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV17_Fake = new TH1F(("EleIDBDTGV17_Fake"+label).c_str(), "; BDTG V17 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV18_Real = new TH1F(("EleIDBDTGV18_Real"+label).c_str(), "; BDTG V18 ; Number of Events ",  10000, -2 , 2);
  TH1F *EleIDBDTGV18_Fake = new TH1F(("EleIDBDTGV18_Fake"+label).c_str(), "; BDTG V18 ; Number of Events ",  10000, -2 , 2);



  TH1F *BDTGSignalEfficiencyNumerator_Pt = new TH1F("BDTGSignalEfficiencyNumerator_Pt" , "; p_{T} [GeV^{c^}];Number of Events",  100, 0, 100);
  TH1F *BDTGSignalEfficiencyNumerator_Eta = new TH1F("BDTGSignalEfficiencyNumerator_Eta" , "; #eta;Number of Events",  100, -2.5, 2.5);
  TH1F *BDTGSignalEfficiencyNumerator_Phi = new TH1F("BDTGSignalEfficiencyNumerator_Phi" , "; #phi;Number of Events",  100, -3.2, 3.2);
  TH1F *BDTGSignalEfficiencyDenominator_Pt = new TH1F("BDTGSignalEfficiencyDenominator_Pt" , "; p_{T} [GeV^{c^}];Number of Events",  100, 0, 100);
  TH1F *BDTGSignalEfficiencyDenominator_Eta = new TH1F("BDTGSignalEfficiencyDenominator_Eta" , "; #eta;Number of Events",  100, -2.5, 2.5);
  TH1F *BDTGSignalEfficiencyDenominator_Phi = new TH1F("BDTGSignalEfficiencyDenominator_Phi" , "; #phi;Number of Events",  100, -3.2, 3.2);


  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassCutBased = 0;
  Double_t FakeElectronPassCutBased = 0;
  Double_t RealElectronPass2011MVA = 0;
  Double_t FakeElectronPass2011MVA = 0;

  Double_t RealElectronPassBDTGV15_SameCutBasedSig = 0;
  Double_t FakeElectronPassBDTGV15_SameCutBasedSig = 0;
  Double_t RealElectronPassBDTGV15_HalfCutBasedBkg = 0;
  Double_t FakeElectronPassBDTGV15_HalfCutBasedBkg = 0;
  Double_t RealElectronPassBDTGV15_OneThirdCutBasedBkg = 0;
  Double_t FakeElectronPassBDTGV15_OneThirdCutBasedBkg = 0;



  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight = 1.0;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fElePFIso; 
  
  //CutBased Variables
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  //Additional Vars used in Likelihood
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  //Isolation Variables
  Float_t                 fEleStandardLikelihood; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  Float_t                 fEleLikelihood;
  Float_t                 fEleKNN;
  Float_t                 fEleMLP;
  Float_t                 fEleMLPBNN;
  Float_t                 fEleBDTG;
  Float_t                 fEleBDT;
  Float_t                 fEleCombinedMVA;
  Float_t                 fEleKNN_C;
  Float_t                 fEleMLP_C;
  Float_t                 fEleMLPBNN_C;
  Float_t                 fEleBDTG_C;
  Float_t                 fEleBDT_C;
  Float_t                 fEleCombinedMVA_C;


  Float_t                 fEleBDTGV0;
  Float_t                 fEleBDTGV1;
  Float_t                 fEleBDTGV2;
  Float_t                 fEleBDTGV3;
  Float_t                 fEleBDTGV4;
  Float_t                 fEleBDTGV5;
  Float_t                 fEleBDTGV6;
  Float_t                 fEleBDTGV7;
  Float_t                 fEleBDTGV8;
  Float_t                 fEleBDTGV9;
  Float_t                 fEleBDTGV10;
  Float_t                 fEleBDTGV11;
  Float_t                 fEleBDTGV12;
  Float_t                 fEleBDTGV13;
  Float_t                 fEleBDTGV14;
  Float_t                 fEleBDTGV15;
  Float_t                 fEleBDTGV16;
  Float_t                 fEleBDTGV17;
  Float_t                 fEleBDTGV18;




  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  TFile *RealEleFile = new TFile(RealElectronFile.c_str(), "READ");
  TTree *RealEleTree = (TTree*)RealEleFile->Get("Electrons");
  RealEleTree->SetBranchAddress( "weight", &fWeight);
  RealEleTree->SetBranchAddress( "run", &fRunNumber);
  RealEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleTree->SetBranchAddress( "event", &fEventNumber);
  RealEleTree->SetBranchAddress( "pt", &fElePt); 
  RealEleTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  RealEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  RealEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  RealEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  RealEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  RealEleTree->SetBranchAddress( "D0", &fEleD0); 
  RealEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  RealEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  RealEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  RealEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  RealEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  RealEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  RealEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  RealEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  RealEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  RealEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  RealEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "Rho", &fRho); 
  RealEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  RealEleTree->SetBranchAddress( "Likelihood", &fEleLikelihood); 
  RealEleTree->SetBranchAddress( "KNN", &fEleKNN); 
  RealEleTree->SetBranchAddress( "MLP", &fEleMLP); 
  RealEleTree->SetBranchAddress( "MLPBNN", &fEleMLPBNN); 
  RealEleTree->SetBranchAddress( "BDTG", &fEleBDTG); 
  RealEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  RealEleTree->SetBranchAddress( "CombinedMVA", &fEleCombinedMVA); 
  RealEleTree->SetBranchAddress( "KNN_C", &fEleKNN_C); 
  RealEleTree->SetBranchAddress( "MLP_C", &fEleMLP_C); 
  RealEleTree->SetBranchAddress( "MLPBNN_C", &fEleMLPBNN_C); 
  RealEleTree->SetBranchAddress( "BDTG_C", &fEleBDTG_C); 
  RealEleTree->SetBranchAddress( "BDT_C", &fEleBDT_C); 
  RealEleTree->SetBranchAddress( "CombinedMVA_C", &fEleCombinedMVA_C); 

  RealEleTree->SetBranchAddress( "BDTGV0", &fEleBDTGV0); 
  RealEleTree->SetBranchAddress( "BDTGV1", &fEleBDTGV1); 
  RealEleTree->SetBranchAddress( "BDTGV2", &fEleBDTGV2); 
  RealEleTree->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  RealEleTree->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  RealEleTree->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  RealEleTree->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  RealEleTree->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  RealEleTree->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  RealEleTree->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 
  RealEleTree->SetBranchAddress( "BDTGV10", &fEleBDTGV10); 
  RealEleTree->SetBranchAddress( "BDTGV11", &fEleBDTGV11); 
  RealEleTree->SetBranchAddress( "BDTGV12", &fEleBDTGV12); 
  RealEleTree->SetBranchAddress( "BDTGV13", &fEleBDTGV13); 
  RealEleTree->SetBranchAddress( "BDTGV14", &fEleBDTGV14); 
  RealEleTree->SetBranchAddress( "BDTGV15", &fEleBDTGV15); 
  RealEleTree->SetBranchAddress( "BDTGV16", &fEleBDTGV16); 
  RealEleTree->SetBranchAddress( "BDTGV17", &fEleBDTGV17); 
  RealEleTree->SetBranchAddress( "BDTGV18", &fEleBDTGV18); 
 

  for(UInt_t ientry=0; ientry < RealEleTree->GetEntries(); ientry++) {       	
    RealEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    //don't evaluate performance using training events
    if (fEventNumber % 2 == 0) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleSCEta) < 1.0) subdet = 0;
    else if (fabs(fEleSCEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply denominator cuts
    if (!(fabs(fEleDZ) < 0.1 && fabs(fEleD0) < 0.02)) continue;

    
    RealElectrons += fWeight;
    BDTGSignalEfficiencyDenominator_Pt->Fill(fElePt);
    BDTGSignalEfficiencyDenominator_Eta->Fill(fEleEta);
    BDTGSignalEfficiencyDenominator_Phi->Fill(fElePhi);

    if (passVBTF(1, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP) && passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) 
      RealElectronPassCutBased += fWeight;
    if (pass2011MVA( fElePt, fEleSCEta, fEleBDTGV8) && passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) 
      RealElectronPass2011MVA += fWeight;
    if (passMVASame2011Sig( fElePt, fEleSCEta, fEleBDTGV15))
      RealElectronPassBDTGV15_SameCutBasedSig += fWeight;
    if (passMVAHalf2011Bkg( fElePt, fEleSCEta, fEleBDTGV15)) {
      RealElectronPassBDTGV15_HalfCutBasedBkg += fWeight;
      BDTGSignalEfficiencyNumerator_Pt->Fill(fElePt);
      BDTGSignalEfficiencyNumerator_Eta->Fill(fEleEta);
      BDTGSignalEfficiencyNumerator_Phi->Fill(fElePhi);
    }


    //Fill Histograms
    EleIDStandardLikelihood_Real->Fill(fEleStandardLikelihood,fWeight);
    EleIDTMVALikelihood_Real->Fill(fEleLikelihood,fWeight);
    EleIDKNN_Real->Fill(fEleKNN,fWeight);
    EleIDBDT_Real->Fill(fEleBDT,fWeight);
    EleIDBDTG_Real->Fill(fEleBDTG,fWeight);
    EleIDMLP_Real->Fill(fEleMLP,fWeight);
    EleIDMLPBNN_Real->Fill(fEleMLPBNN,fWeight);
    EleIDCombinedMVA_Real->Fill(fEleCombinedMVA,fWeight);

    EleIDKNN_C_Real->Fill(fEleKNN_C,fWeight);
    EleIDBDT_C_Real->Fill(fEleBDT_C,fWeight);
    EleIDBDTG_C_Real->Fill(fEleBDTG_C,fWeight);
    EleIDMLP_C_Real->Fill(fEleMLP_C,fWeight);
    EleIDMLPBNN_C_Real->Fill(fEleMLPBNN_C,fWeight);
    EleIDCombinedMVA_C_Real->Fill(fEleCombinedMVA_C,fWeight);

    
    if (!passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) {
      //if electron doesn't pass nominal pfiso cut, artificially make mva value very small
      EleIDBDTGV0_Real->Fill(-1.99,fWeight);
      EleIDBDTGV1_Real->Fill(-1.99,fWeight);
      EleIDBDTGV2_Real->Fill(-1.99,fWeight);
      EleIDBDTGV3_Real->Fill(-1.99,fWeight);
      EleIDBDTGV4_Real->Fill(-1.99,fWeight);
      EleIDBDTGV5_Real->Fill(-1.99,fWeight);
      EleIDBDTGV6_Real->Fill(-1.99,fWeight);
      EleIDBDTGV7_Real->Fill(-1.99,fWeight);
      EleIDBDTGV8_Real->Fill(-1.99,fWeight);
      EleIDBDTGV9_Real->Fill(-1.99,fWeight);      
    } else {
      EleIDBDTGV0_Real->Fill(fEleBDTGV0,fWeight);
      EleIDBDTGV1_Real->Fill(fEleBDTGV1,fWeight);
      EleIDBDTGV2_Real->Fill(fEleBDTGV2,fWeight);
      EleIDBDTGV3_Real->Fill(fEleBDTGV3,fWeight);
      EleIDBDTGV4_Real->Fill(fEleBDTGV4,fWeight);
      EleIDBDTGV5_Real->Fill(fEleBDTGV5,fWeight);
      EleIDBDTGV6_Real->Fill(fEleBDTGV6,fWeight);
      EleIDBDTGV7_Real->Fill(fEleBDTGV7,fWeight);
      EleIDBDTGV8_Real->Fill(fEleBDTGV8,fWeight);
      EleIDBDTGV9_Real->Fill(fEleBDTGV9,fWeight);      
    }

    EleIDBDTGV10_Real->Fill(fEleBDTGV10,fWeight);
    EleIDBDTGV11_Real->Fill(fEleBDTGV11,fWeight);
    EleIDBDTGV12_Real->Fill(fEleBDTGV12,fWeight);
    EleIDBDTGV13_Real->Fill(fEleBDTGV13,fWeight);
    EleIDBDTGV14_Real->Fill(fEleBDTGV14,fWeight);
    EleIDBDTGV15_Real->Fill(fEleBDTGV15,fWeight);
    EleIDBDTGV16_Real->Fill(fEleBDTGV16,fWeight);
    EleIDBDTGV17_Real->Fill(fEleBDTGV17,fWeight);
    EleIDBDTGV18_Real->Fill(fEleBDTGV18,fWeight);


  } 
  





  //*****************************************************************************************
  //FakeEleTree
  //*****************************************************************************************
  TFile *FakeEleFile = new TFile(FakeElectronFile.c_str(), "READ");
  TTree *FakeEleTree = (TTree*)FakeEleFile->Get("Electrons");
  FakeEleTree->SetBranchAddress( "weight", &fWeight);
  FakeEleTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  FakeEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  FakeEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  FakeEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  FakeEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  FakeEleTree->SetBranchAddress( "D0", &fEleD0); 
  FakeEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  FakeEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  FakeEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  FakeEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  FakeEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  FakeEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  FakeEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  FakeEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  FakeEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  FakeEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  FakeEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  FakeEleTree->SetBranchAddress( "Likelihood", &fEleLikelihood); 
  FakeEleTree->SetBranchAddress( "KNN", &fEleKNN); 
  FakeEleTree->SetBranchAddress( "MLP", &fEleMLP); 
  FakeEleTree->SetBranchAddress( "MLPBNN", &fEleMLPBNN); 
  FakeEleTree->SetBranchAddress( "BDTG", &fEleBDTG); 
  FakeEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  FakeEleTree->SetBranchAddress( "CombinedMVA", &fEleCombinedMVA); 
  FakeEleTree->SetBranchAddress( "KNN_C", &fEleKNN_C); 
  FakeEleTree->SetBranchAddress( "MLP_C", &fEleMLP_C); 
  FakeEleTree->SetBranchAddress( "MLPBNN_C", &fEleMLPBNN_C); 
  FakeEleTree->SetBranchAddress( "BDTG_C", &fEleBDTG_C); 
  FakeEleTree->SetBranchAddress( "BDT_C", &fEleBDT_C); 
  FakeEleTree->SetBranchAddress( "CombinedMVA_C", &fEleCombinedMVA_C); 

  FakeEleTree->SetBranchAddress( "BDTGV0", &fEleBDTGV0); 
  FakeEleTree->SetBranchAddress( "BDTGV1", &fEleBDTGV1); 
  FakeEleTree->SetBranchAddress( "BDTGV2", &fEleBDTGV2); 
  FakeEleTree->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  FakeEleTree->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  FakeEleTree->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  FakeEleTree->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  FakeEleTree->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  FakeEleTree->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  FakeEleTree->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 
  FakeEleTree->SetBranchAddress( "BDTGV10", &fEleBDTGV10); 
  FakeEleTree->SetBranchAddress( "BDTGV11", &fEleBDTGV11); 
  FakeEleTree->SetBranchAddress( "BDTGV12", &fEleBDTGV12); 
  FakeEleTree->SetBranchAddress( "BDTGV13", &fEleBDTGV13); 
  FakeEleTree->SetBranchAddress( "BDTGV14", &fEleBDTGV14); 
  FakeEleTree->SetBranchAddress( "BDTGV15", &fEleBDTGV15); 
  FakeEleTree->SetBranchAddress( "BDTGV16", &fEleBDTGV16); 
  FakeEleTree->SetBranchAddress( "BDTGV17", &fEleBDTGV17); 
  FakeEleTree->SetBranchAddress( "BDTGV18", &fEleBDTGV18); 

  for(UInt_t ientry=0; ientry < FakeEleTree->GetEntries(); ientry++) {       	
    FakeEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //don't evaluate performance using training events
    if (fEventNumber % 2 == 0) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleSCEta) < 1.0) subdet = 0;
    else if (fabs(fEleSCEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply denominator cuts
    if (!(fabs(fEleDZ) < 0.1 && fabs(fEleD0) < 0.02)) continue;


    FakeElectrons += fWeight;

    if (passVBTF(1, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP) && passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) 
      FakeElectronPassCutBased += fWeight;
    if (pass2011MVA( fElePt, fEleSCEta, fEleBDTGV8) && passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) 
      FakeElectronPass2011MVA += fWeight;
    if (passMVASame2011Sig( fElePt, fEleSCEta, fEleBDTGV15))
      FakeElectronPassBDTGV15_SameCutBasedSig += fWeight;
    if (passMVAHalf2011Bkg( fElePt, fEleSCEta, fEleBDTGV15))
      FakeElectronPassBDTGV15_HalfCutBasedBkg += fWeight;



    //Fill Histograms
    EleIDStandardLikelihood_Fake->Fill(fEleStandardLikelihood,fWeight);
    EleIDTMVALikelihood_Fake->Fill(fEleLikelihood,fWeight);
    EleIDKNN_Fake->Fill(fEleKNN,fWeight);
    EleIDBDT_Fake->Fill(fEleBDT,fWeight);
    EleIDBDTG_Fake->Fill(fEleBDTG,fWeight);
    EleIDMLP_Fake->Fill(fEleMLP,fWeight);
    EleIDMLPBNN_Fake->Fill(fEleMLPBNN,fWeight);
    EleIDCombinedMVA_Fake->Fill(fEleCombinedMVA,fWeight);

    EleIDKNN_C_Fake->Fill(fEleKNN_C,fWeight);
    EleIDBDT_C_Fake->Fill(fEleBDT_C,fWeight);
    EleIDBDTG_C_Fake->Fill(fEleBDTG_C,fWeight);
    EleIDMLP_C_Fake->Fill(fEleMLP_C,fWeight);
    EleIDMLPBNN_C_Fake->Fill(fEleMLPBNN_C,fWeight);
    EleIDCombinedMVA_C_Fake->Fill(fEleCombinedMVA_C,fWeight);

    if (!passCutBasedIsoOnly(fElePt, fEleEta, fElePFIso)) {
      EleIDBDTGV0_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV1_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV2_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV3_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV4_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV5_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV6_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV7_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV8_Fake->Fill(-1.99,fWeight);
      EleIDBDTGV9_Fake->Fill(-1.99,fWeight);      
    } else {
      EleIDBDTGV0_Fake->Fill(fEleBDTGV0,fWeight);
      EleIDBDTGV1_Fake->Fill(fEleBDTGV1,fWeight);
      EleIDBDTGV2_Fake->Fill(fEleBDTGV2,fWeight);
      EleIDBDTGV3_Fake->Fill(fEleBDTGV3,fWeight);
      EleIDBDTGV4_Fake->Fill(fEleBDTGV4,fWeight);
      EleIDBDTGV5_Fake->Fill(fEleBDTGV5,fWeight);
      EleIDBDTGV6_Fake->Fill(fEleBDTGV6,fWeight);
      EleIDBDTGV7_Fake->Fill(fEleBDTGV7,fWeight);
      EleIDBDTGV8_Fake->Fill(fEleBDTGV8,fWeight);
      EleIDBDTGV9_Fake->Fill(fEleBDTGV9,fWeight);
    }

    EleIDBDTGV10_Fake->Fill(fEleBDTGV10,fWeight);
    EleIDBDTGV11_Fake->Fill(fEleBDTGV11,fWeight);
    EleIDBDTGV12_Fake->Fill(fEleBDTGV12,fWeight);
    EleIDBDTGV13_Fake->Fill(fEleBDTGV13,fWeight);
    EleIDBDTGV14_Fake->Fill(fEleBDTGV14,fWeight);
    EleIDBDTGV15_Fake->Fill(fEleBDTGV15,fWeight);
    EleIDBDTGV16_Fake->Fill(fEleBDTGV16,fWeight);
    EleIDBDTGV17_Fake->Fill(fEleBDTGV17,fWeight);
    EleIDBDTGV18_Fake->Fill(fEleBDTGV18,fWeight);

  } //loop over electrons
  



  
  //*****************************************************************************************
  //Current Working Points
  //*****************************************************************************************
  cout << "Cut-Based Real Electron Efficiency : " << RealElectronPassCutBased << " / " << RealElectrons << " = " << RealElectronPassCutBased/RealElectrons << endl;
  cout << "Cut-Based Fake Electron Efficiency : " << FakeElectronPassCutBased << " / " << FakeElectrons << " = " << FakeElectronPassCutBased/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_CutBasedCurrentWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCutBased/RealElectrons , FakeElectronPassCutBased/FakeElectrons, "ROC_CutBasedCurrentWP"+label);

  cout << "2011 MVA Real Electron Efficiency : " << RealElectronPass2011MVA << " / " << RealElectrons << " = " << RealElectronPass2011MVA/RealElectrons << endl;
  cout << "2011 MVA Fake Electron Efficiency : " << FakeElectronPass2011MVA << " / " << FakeElectrons << " = " << FakeElectronPass2011MVA/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_2011MVACurrentWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPass2011MVA/RealElectrons , FakeElectronPass2011MVA/FakeElectrons, "ROC_2011MVACurrentWP"+label);


  cout << "BDTGV15 SameCutBasedSig Real Electron Efficiency : " << RealElectronPassBDTGV15_SameCutBasedSig << " / " << RealElectrons << " = " << RealElectronPassBDTGV15_SameCutBasedSig/RealElectrons << endl;
  cout << "BDTGV15 SameCutBasedSig Fake Electron Efficiency : " << FakeElectronPassBDTGV15_SameCutBasedSig << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV15_SameCutBasedSig/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV15_SameCutBasedSig_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV15_SameCutBasedSig/RealElectrons , FakeElectronPassBDTGV15_SameCutBasedSig/FakeElectrons, "ROC_BDTGV15SameCutBasedSigWP"+label);

  cout << "BDTGV15 HalfCutBasedBkg Real Electron Efficiency : " << RealElectronPassBDTGV15_HalfCutBasedBkg << " / " << RealElectrons << " = " << RealElectronPassBDTGV15_HalfCutBasedBkg/RealElectrons << endl;
  cout << "BDTGV15 HalfCutBasedBkg Fake Electron Efficiency : " << FakeElectronPassBDTGV15_HalfCutBasedBkg << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV15_HalfCutBasedBkg/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV15_HalfCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV15_HalfCutBasedBkg/RealElectrons , FakeElectronPassBDTGV15_HalfCutBasedBkg/FakeElectrons, "ROC_BDTGV15HalfCutBasedBkgWP"+label);

  cout << "BDTGV15 OneThirdCutBasedBkg Real Electron Efficiency : " << RealElectronPassBDTGV15_OneThirdCutBasedBkg << " / " << RealElectrons << " = " << RealElectronPassBDTGV15_OneThirdCutBasedBkg/RealElectrons << endl;
  cout << "BDTGV15 OneThirdCutBasedBkg Fake Electron Efficiency : " << FakeElectronPassBDTGV15_OneThirdCutBasedBkg << " / " << FakeElectrons << " = " << FakeElectronPassBDTGV15_OneThirdCutBasedBkg/FakeElectrons << endl;
  TGraphAsymmErrors* ROC_BDTGV15_OneThirdCutBasedBkg_WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassBDTGV15_OneThirdCutBasedBkg/RealElectrons , FakeElectronPassBDTGV15_OneThirdCutBasedBkg/FakeElectrons, "ROC_BDTGV15OneThirdCutBasedBkgWP"+label);

  cout << "**********************\n";
  cout << "Bkg At LHTight Signal Eff\n";

  Double_t BkgEff2011MVA = FakeElectronPass2011MVA/FakeElectrons;
  Double_t SigEff2011MVA = RealElectronPass2011MVA/RealElectrons;
  Double_t SigEffBDTGV0_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV1_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV2_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV3_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV3_Real, EleIDBDTGV3_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV4_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV4_Real, EleIDBDTGV4_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV5_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV5_Real, EleIDBDTGV5_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV6_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV6_Real, EleIDBDTGV6_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV7_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV7_Real, EleIDBDTGV7_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV8_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV8_Real, EleIDBDTGV8_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV9_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV9_Real, EleIDBDTGV9_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV10_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV10_Real, EleIDBDTGV10_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV11_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV11_Real, EleIDBDTGV11_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV12_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV12_Real, EleIDBDTGV12_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV13_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV13_Real, EleIDBDTGV13_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV14_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV14_Real, EleIDBDTGV14_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV15_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV16_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV16_Real, EleIDBDTGV16_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV17_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV17_Real, EleIDBDTGV17_Fake, 0.5*BkgEff2011MVA);
  Double_t SigEffBDTGV18_HalfBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV18_Real, EleIDBDTGV18_Fake, 0.5*BkgEff2011MVA);

  Double_t SigEffBDTGV0_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV1_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV2_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV3_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV3_Real, EleIDBDTGV3_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV4_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV4_Real, EleIDBDTGV4_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV5_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV5_Real, EleIDBDTGV5_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV6_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV6_Real, EleIDBDTGV6_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV7_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV7_Real, EleIDBDTGV7_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV8_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV8_Real, EleIDBDTGV8_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV9_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV9_Real, EleIDBDTGV9_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV10_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV10_Real, EleIDBDTGV10_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV11_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV11_Real, EleIDBDTGV11_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV12_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV12_Real, EleIDBDTGV12_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV13_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV13_Real, EleIDBDTGV13_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV14_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV14_Real, EleIDBDTGV14_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV15_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV16_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV16_Real, EleIDBDTGV16_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV17_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV17_Real, EleIDBDTGV17_Fake, (1.0/3.0)*BkgEff2011MVA);
  Double_t SigEffBDTGV18_OneThirdBkg = FindSigEffAtFixedBkgEfficiency(EleIDBDTGV18_Real, EleIDBDTGV18_Fake, (1.0/3.0)*BkgEff2011MVA);


  cout << "Signal Efficiency (wrt Cut-based) for : half bkg : one third bkg \n";
  cout << "BDTGV0 : " << SigEffBDTGV0_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV0_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV1 : " << SigEffBDTGV1_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV1_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV2 : " << SigEffBDTGV2_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV2_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV3 : " << SigEffBDTGV3_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV3_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV4 : " << SigEffBDTGV4_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV4_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV5 : " << SigEffBDTGV5_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV5_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV6 : " << SigEffBDTGV6_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV6_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV7 : " << SigEffBDTGV7_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV7_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV8 : " << SigEffBDTGV8_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV8_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV9 : " << SigEffBDTGV9_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV9_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV10 : " << SigEffBDTGV10_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV10_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV11 : " << SigEffBDTGV11_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV11_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV12 : " << SigEffBDTGV12_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV12_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV13 : " << SigEffBDTGV13_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV13_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV14 : " << SigEffBDTGV14_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV14_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV15 : " << SigEffBDTGV15_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV15_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV16 : " << SigEffBDTGV16_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV16_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV17 : " << SigEffBDTGV17_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV17_OneThirdBkg/SigEff2011MVA << endl;
  cout << "BDTGV18 : " << SigEffBDTGV18_HalfBkg/SigEff2011MVA << " : " << SigEffBDTGV18_OneThirdBkg/SigEff2011MVA << endl;


  cout << "**********************\n";

  Double_t BkgEffBDTGV0_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV1_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV2_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV3_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV3_Real, EleIDBDTGV3_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV4_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV4_Real, EleIDBDTGV4_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV5_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV5_Real, EleIDBDTGV5_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV6_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV6_Real, EleIDBDTGV6_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV7_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV7_Real, EleIDBDTGV7_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV8_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV8_Real, EleIDBDTGV8_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV9_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV9_Real, EleIDBDTGV9_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV10_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV10_Real, EleIDBDTGV10_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV11_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV11_Real, EleIDBDTGV11_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV12_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV12_Real, EleIDBDTGV12_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV13_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV13_Real, EleIDBDTGV13_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV14_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV14_Real, EleIDBDTGV14_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV15_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV16_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV16_Real, EleIDBDTGV16_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV17_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV17_Real, EleIDBDTGV17_Fake, SigEff2011MVA);
  Double_t BkgEffBDTGV18_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV18_Real, EleIDBDTGV18_Fake, SigEff2011MVA);

  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "BDTGV0 : " << BkgEffBDTGV0_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV1 : " << BkgEffBDTGV1_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV2 : " << BkgEffBDTGV2_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV3 : " << BkgEffBDTGV3_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV4 : " << BkgEffBDTGV4_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV5 : " << BkgEffBDTGV5_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV6 : " << BkgEffBDTGV6_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV7 : " << BkgEffBDTGV7_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV8 : " << BkgEffBDTGV8_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV9 : " << BkgEffBDTGV9_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV10 : " << BkgEffBDTGV10_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV11 : " << BkgEffBDTGV11_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV12 : " << BkgEffBDTGV12_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV13 : " << BkgEffBDTGV13_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV14 : " << BkgEffBDTGV14_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV15 : " << BkgEffBDTGV15_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV16 : " << BkgEffBDTGV16_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV17 : " << BkgEffBDTGV17_SameSig/BkgEff2011MVA << endl;
  cout << "BDTGV18 : " << BkgEffBDTGV18_SameSig/BkgEff2011MVA << endl;

  cout << "**********************\n";



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
//   TGraphAsymmErrors* ROC_StandardLikelihood = MakeSigEffVsBkgEffGraph(EleIDStandardLikelihood_Real, EleIDStandardLikelihood_Fake, "ROC_StandardLikelihood"+label );
//   TGraphAsymmErrors* ROC_TMVALikelihood = MakeSigEffVsBkgEffGraph(EleIDTMVALikelihood_Real, EleIDTMVALikelihood_Fake, "ROC_TMVALikelihood"+label );
//   TGraphAsymmErrors* ROC_KNN = MakeSigEffVsBkgEffGraph(EleIDKNN_Real, EleIDKNN_Fake, "ROC_KNN"+label );
//   TGraphAsymmErrors* ROC_MLP = MakeSigEffVsBkgEffGraph(EleIDMLP_Real, EleIDMLP_Fake, "ROC_MLP"+label );
//   TGraphAsymmErrors* ROC_MLPBNN = MakeSigEffVsBkgEffGraph(EleIDMLPBNN_Real, EleIDMLPBNN_Fake, "ROC_MLPBNN"+label );
//   TGraphAsymmErrors* ROC_BDT = MakeSigEffVsBkgEffGraph(EleIDBDT_Real, EleIDBDT_Fake, "ROC_BDT"+label );
//   TGraphAsymmErrors* ROC_BDTG = MakeSigEffVsBkgEffGraph(EleIDBDTG_Real, EleIDBDTG_Fake, "ROC_BDTG"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA = MakeSigEffVsBkgEffGraph(EleIDCombinedMVA_Real, EleIDCombinedMVA_Fake, "ROC_CombinedMVA"+label );

//   TGraphAsymmErrors* ROC_BDTG_C = MakeSigEffVsBkgEffGraph(EleIDBDTG_C_Real, EleIDBDTG_C_Fake, "ROC_BDTG_C"+label );
//   TGraphAsymmErrors* ROC_CombinedMVA_C = MakeSigEffVsBkgEffGraph(EleIDCombinedMVA_C_Real, EleIDCombinedMVA_C_Fake, "ROC_CombinedMVA_C"+label );


  TGraphAsymmErrors* ROC_BDTGV0 = MakeSigEffVsBkgEffGraph(EleIDBDTGV0_Real, EleIDBDTGV0_Fake, "ROC_BDTGV0"+label );
  TGraphAsymmErrors* ROC_BDTGV1 = MakeSigEffVsBkgEffGraph(EleIDBDTGV1_Real, EleIDBDTGV1_Fake, "ROC_BDTGV1"+label );
  TGraphAsymmErrors* ROC_BDTGV2 = MakeSigEffVsBkgEffGraph(EleIDBDTGV2_Real, EleIDBDTGV2_Fake, "ROC_BDTGV2"+label );
  TGraphAsymmErrors* ROC_BDTGV3 = MakeSigEffVsBkgEffGraph(EleIDBDTGV3_Real, EleIDBDTGV3_Fake, "ROC_BDTGV3"+label );
  TGraphAsymmErrors* ROC_BDTGV4 = MakeSigEffVsBkgEffGraph(EleIDBDTGV4_Real, EleIDBDTGV4_Fake, "ROC_BDTGV4"+label );
  TGraphAsymmErrors* ROC_BDTGV5 = MakeSigEffVsBkgEffGraph(EleIDBDTGV5_Real, EleIDBDTGV5_Fake, "ROC_BDTGV5"+label );
  TGraphAsymmErrors* ROC_BDTGV6 = MakeSigEffVsBkgEffGraph(EleIDBDTGV6_Real, EleIDBDTGV6_Fake, "ROC_BDTGV6"+label );
  TGraphAsymmErrors* ROC_BDTGV7 = MakeSigEffVsBkgEffGraph(EleIDBDTGV7_Real, EleIDBDTGV7_Fake, "ROC_BDTGV7"+label );
  TGraphAsymmErrors* ROC_BDTGV8 = MakeSigEffVsBkgEffGraph(EleIDBDTGV8_Real, EleIDBDTGV8_Fake, "ROC_BDTGV8"+label );
  TGraphAsymmErrors* ROC_BDTGV9 = MakeSigEffVsBkgEffGraph(EleIDBDTGV9_Real, EleIDBDTGV9_Fake, "ROC_BDTGV9"+label );
  TGraphAsymmErrors* ROC_BDTGV10 = MakeSigEffVsBkgEffGraph(EleIDBDTGV10_Real, EleIDBDTGV10_Fake, "ROC_BDTGV10"+label );
  TGraphAsymmErrors* ROC_BDTGV11 = MakeSigEffVsBkgEffGraph(EleIDBDTGV11_Real, EleIDBDTGV11_Fake, "ROC_BDTGV11"+label );
  TGraphAsymmErrors* ROC_BDTGV12 = MakeSigEffVsBkgEffGraph(EleIDBDTGV12_Real, EleIDBDTGV12_Fake, "ROC_BDTGV12"+label );
  TGraphAsymmErrors* ROC_BDTGV13 = MakeSigEffVsBkgEffGraph(EleIDBDTGV13_Real, EleIDBDTGV13_Fake, "ROC_BDTGV13"+label );
  TGraphAsymmErrors* ROC_BDTGV14 = MakeSigEffVsBkgEffGraph(EleIDBDTGV14_Real, EleIDBDTGV14_Fake, "ROC_BDTGV14"+label );
  TGraphAsymmErrors* ROC_BDTGV15 = MakeSigEffVsBkgEffGraph(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, "ROC_BDTGV15"+label );
  TGraphAsymmErrors* ROC_BDTGV16 = MakeSigEffVsBkgEffGraph(EleIDBDTGV16_Real, EleIDBDTGV16_Fake, "ROC_BDTGV16"+label );
  TGraphAsymmErrors* ROC_BDTGV17 = MakeSigEffVsBkgEffGraph(EleIDBDTGV17_Real, EleIDBDTGV17_Fake, "ROC_BDTGV17"+label );
  TGraphAsymmErrors* ROC_BDTGV18 = MakeSigEffVsBkgEffGraph(EleIDBDTGV18_Real, EleIDBDTGV18_Fake, "ROC_BDTGV18"+label );

  //*****************************************************************************************
  //Find Cut with same signal efficiency Make ROC curves
  //*****************************************************************************************
  Double_t CutValue_BDTGV15_SameSig = FindCutValueAtFixedEfficiency(EleIDBDTGV15_Real, SigEff2011MVA );
  Double_t CutValue_BDTGV15_HalfBkg = FindCutValueAtFixedEfficiency(EleIDBDTGV15_Fake, 0.5*BkgEff2011MVA );
  Double_t CutValue_BDTGV15_OneThirdBkg = FindCutValueAtFixedEfficiency(EleIDBDTGV15_Fake, (1.0/3.0)*BkgEff2011MVA );
  cout << "BDTG V15 Cut Value @ Same Cut-Based Sig: " << CutValue_BDTGV15_SameSig << endl;
  cout << "BDTG V15 Cut Value @ 50% Cut-Based Bkg: " << CutValue_BDTGV15_HalfBkg << endl;
  cout << "BDTG V15 Cut Value @ 33% Cut-Based Bkg: " << CutValue_BDTGV15_OneThirdBkg << endl;

//   TFile *canvasFile = new TFile("ElectronIDMVAPerformancePlots.root","UPDATE");
  TLegend* legend;
  TCanvas* cv;
  string plotname;

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
//   colors.push_back(kBlue);    markers.push_back(21);
//   colors.push_back(kMagenta); markers.push_back(22);
//   colors.push_back(kCyan);    markers.push_back(34);
//   colors.push_back(kBlack);   markers.push_back(29);
//   colors.push_back(kGreen);   markers.push_back(33);
//   colors.push_back(kRed-2);   markers.push_back(33);
//   colors.push_back(kOrange);  markers.push_back(33);
//   colors.push_back(kBlue-2);  markers.push_back(33);
//   colors.push_back(kGreen-2);  markers.push_back(33);
//   colors.push_back(kMagenta-2);   markers.push_back(33);

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronIDMVA"+label;


//   ROCGraphs.push_back(ROC_BDTGV8);
//   GraphLabels.push_back("2011MVA (V0)");
//   colors.push_back(kRed);
  
  ROCGraphs.push_back(ROC_BDTGV0);
  GraphLabels.push_back("V0 retrained");
  colors.push_back(kBlue);
  
//   ROCGraphs.push_back(ROC_BDTGV1);
//   GraphLabels.push_back("V0+Track-Related");
//   colors.push_back(kGreen+2);
  
//   ROCGraphs.push_back(ROC_BDTGV2);
//   GraphLabels.push_back("V0+ShowerShape");
//   colors.push_back(kMagenta);
  
//    ROCGraphs.push_back(ROC_BDTGV3);
//    GraphLabels.push_back("V0+Preshower");
//    colors.push_back(kCyan);
   
//   ROCGraphs.push_back(ROC_BDTGV4);
//   GraphLabels.push_back("V0+H/E");
//   colors.push_back(kOrange);

//   ROCGraphs.push_back(ROC_BDTGV5);
//   GraphLabels.push_back("V0+Single Crystal");
//   colors.push_back(kBlack);

//   ROCGraphs.push_back(ROC_BDTGV6);
//   GraphLabels.push_back("All ID Combined");
//   colors.push_back(kAzure+3);

//   ROCGraphs.push_back(ROC_BDTGV7);
//   GraphLabels.push_back("V7");
//   colors.push_back(kRed-2);

//   ROCGraphs.push_back(ROC_BDTGV9);
//   GraphLabels.push_back("All ID Combined (except H/E)");
//   colors.push_back(kRed-2);

//   ROCGraphs.push_back(ROC_BDTGV10);
//   GraphLabels.push_back("V0+PFIso 0.3");
//   colors.push_back(kMagenta+2);

//   ROCGraphs.push_back(ROC_BDTGV11);
//   GraphLabels.push_back("V0+PFIso 0.3 & 0.4");
//   colors.push_back(kCyan+2);

//   ROCGraphs.push_back(ROC_BDTGV12);
//   GraphLabels.push_back("V12");
//   colors.push_back(kGreen);
  
//   ROCGraphs.push_back(ROC_BDTGV13);
//   GraphLabels.push_back("V13");
//   colors.push_back(kMagenta-2);

//   ROCGraphs.push_back(ROC_BDTGV14);
//   GraphLabels.push_back("V0+AllID+AllIso");
//   colors.push_back(kGreen+2);

  ROCGraphs.push_back(ROC_BDTGV15);
  GraphLabels.push_back("V0+AllID+AllRelIso");
  colors.push_back(kMagenta);

//   ROCGraphs.push_back(ROC_BDTGV16);
//   GraphLabels.push_back("V0+AllID+DetIso03");
//   colors.push_back(kGreen);

//   ROCGraphs.push_back(ROC_BDTGV17);
//   GraphLabels.push_back("V0+AllID+DetIso04");
//   colors.push_back(kCyan+2);

//   ROCGraphs.push_back(ROC_BDTGV18);
//   GraphLabels.push_back("V0+AllID+AllDetIso");
//   colors.push_back(kAzure+3);

  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;
// //   if (Option == 0 )                              { xmin = 0.15; xmax = 0.45; ymin = 0.75; ymax = 0.85; }

//   //For MC
//    if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.66; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//    if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//    if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.93; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//    if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.75; ymax = 0.85; }

//FOr Data
//   if (Option == 0 )                              { xmin = 0.10; xmax = 0.50; ymin = 0.70; ymax = 0.90; }
//   if (Option == 1 )                              { xmin = 0.05; xmax = 0.35; ymin = 0.30; ymax = 0.90; }
//   if (Option == 2 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.30; ymax = 0.80; }
//   if (Option == 3 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.80; ymax = 1.00; }
//   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.70; ymax = 1.00; }

//   //**************************************
//   //Data: V0 vs V0+New ID Variables
//   //**************************************
//    if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.66; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//    if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//    if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.93; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//    if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.75; ymax = 0.85; }

//   //**************************************
//   //Data: V0 vs V0+Iso
//   //**************************************
//   if (Option == 0 )                                { xmin = 0.055; xmax = 0.090; ymin = 0.54; ymax = 0.70; }
// //   if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
//   if (Option == 2 )                              { xmin = 0.02; xmax = 0.04; ymin = 0.30; ymax = 0.45; }
//   if (Option == 3 )                                { xmin = 0.11; xmax = 0.145; ymin = 0.84; ymax = 0.95; }
// //   if (Option == 4 )                              { xmin = 0.20; xmax = 0.65; ymin = 0.75; ymax = 1.00; }
//   if (Option == 5 )                              { xmin = 0.07; xmax = 0.12; ymin = 0.77; ymax = 0.88; }
  
  //**************************************
  //Data: V0 vs Best MVAs
  //**************************************
  if (Option == 0 )                                { xmin = 0.040; xmax = 0.090; ymin = 0.40; ymax = 0.70; }
  if (Option == 1 )                              { xmin = 0.00; xmax = 0.30; ymin = 0.00; ymax = 1.00; }
  if (Option == 2 )                              { xmin = 0.015; xmax = 0.035; ymin = 0.30; ymax = 0.45; }
  if (Option == 3 )                                { xmin = 0.09; xmax = 0.15; ymin = 0.75; ymax = 0.95; }
  if (Option == 4 )                              { xmin = 0.00; xmax = 0.2; ymin = 0.75; ymax = 1.00; }
  if (Option == 5 )                              { xmin = 0.05; xmax = 0.12; ymin = 0.65; ymax = 0.90; }


//HWW115
//   if (Option == 3 )                                { xmin = 0.07; xmax = 0.16; ymin = 0.75; ymax = 0.95; }
//Zee
//   if (Option == 3 )                                { xmin = 0.07; xmax = 0.16; ymin = 0.75; ymax = 0.95; }
//WJets vs HWW115
//    if (Option == 3 )                                { xmin = 0.009; xmax = 0.045; ymin = 0.75; ymax = 0.96; }



  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.54,0.14,0.94,0.44);
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
 
//   legend->AddEntry(ROC_CutBasedCurrentWP, "2011 HWW Cut-Based", "P");
//   ROC_CutBasedCurrentWP->SetFillColor(kBlue);
//   ROC_CutBasedCurrentWP->SetMarkerColor(kBlue);
//   ROC_CutBasedCurrentWP->SetMarkerStyle(34);
//   ROC_CutBasedCurrentWP->SetMarkerSize(2.5);
//   ROC_CutBasedCurrentWP->Draw("Psame");


  legend->AddEntry(ROC_2011MVACurrentWP, "2011 MVA", "P");
  ROC_2011MVACurrentWP->SetFillColor(kGreen+3);
  ROC_2011MVACurrentWP->SetMarkerColor(kGreen+3);
  ROC_2011MVACurrentWP->SetMarkerStyle(34);
  ROC_2011MVACurrentWP->SetMarkerSize(2.5);
  ROC_2011MVACurrentWP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV15_SameCutBasedSig_WP, "BDT @ Same 2011 MVA Sig", "P");
//   ROC_BDTGV15_SameCutBasedSig_WP->SetFillColor(kMagenta);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerColor(kMagenta);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerStyle(34);
//   ROC_BDTGV15_SameCutBasedSig_WP->SetMarkerSize(2.5);
//   ROC_BDTGV15_SameCutBasedSig_WP->Draw("Psame");

//   legend->AddEntry(ROC_BDTGV15_HalfCutBasedBkg_WP, "BDT @ 50% 2011 MVA Bkg", "P");
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetFillColor(kRed);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerColor(kRed);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerStyle(34);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->SetMarkerSize(2.5);
//   ROC_BDTGV15_HalfCutBasedBkg_WP->Draw("Psame");

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
  TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Pt, BDTGSignalEfficiencyDenominator_Pt, ("SignalEfficiency"+label).c_str(), ptbins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Eta, BDTGSignalEfficiencyDenominator_Eta, ("SignalEfficiency"+label).c_str(), etabins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(BDTGSignalEfficiencyNumerator_Phi, BDTGSignalEfficiencyDenominator_Phi, ("SignalEfficiency"+label).c_str(), phibins, ErrorType, -99, -99, 0, 1);

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
  TFile *file = new TFile("ElectronIDMVAResults.root", "UPDATE");
//   file->cd();
//   file->WriteTObject(EleIDStandardLikelihood_Real, EleIDStandardLikelihood_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDT_Real, EleIDBDT_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDTG_Real, EleIDBDTG_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDNN_Real, EleIDNN_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihood_Real, EleIDTMVALikelihood_Real->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihoodD_Real, EleIDTMVALikelihoodD_Real->GetName(), "WriteDelete");  

//   file->WriteTObject(EleIDStandardLikelihood_Fake, EleIDStandardLikelihood_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDT_Fake, EleIDBDT_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDTG_Fake, EleIDBDTG_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDNN_Fake, EleIDNN_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihood_Fake, EleIDTMVALikelihood_Fake->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihoodD_Fake, EleIDTMVALikelihoodD_Fake->GetName(), "WriteDelete"); 

//   file->WriteTObject(ROC_StandardLikelihood, ROC_StandardLikelihood->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_BDT, ROC_BDT->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_BDTG, ROC_BDTG->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_NN, ROC_NN->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihood, ROC_TMVALikelihood->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihoodD, ROC_TMVALikelihoodD->GetName(), "WriteDelete");  
 
   file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");  
   file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");  
   file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");  

//   file->Close();
//   delete file;




  gBenchmark->Show("WWTemplate");       
} 

