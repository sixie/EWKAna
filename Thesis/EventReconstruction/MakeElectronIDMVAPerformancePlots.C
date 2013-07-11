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

  Double_t BkgEffCutBased = FakeElectronPassCutBased/FakeElectrons;
  Double_t SigEffCutBased = RealElectronPassCutBased/RealElectrons;
  Double_t BkgEffBDTGV15_SameSig = FindBkgEffAtFixedSignalEfficiency(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, SigEffCutBased);

  cout << "Bkg Efficiency (wrt Cut-based) for same sig eff \n";
  cout << "BDT : " << BkgEffBDTGV15_SameSig/BkgEffCutBased << endl;

  cout << "**********************\n";



  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_BDTGV15 = MakeSigEffVsBkgEffGraph(EleIDBDTGV15_Real, EleIDBDTGV15_Fake, "ROC_BDTGV15"+label );

  TLegend* legend;
  TCanvas* cv;
  string plotname;

  //*****************************************************************************************
  //Plot ROC Curves
  //*****************************************************************************************
  vector<TGraphAsymmErrors*> ROCGraphs;
  vector<string> GraphLabels;
  vector<Int_t> colors;

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronIDMVA"+label;

  ROCGraphs.push_back(ROC_BDTGV15);
  GraphLabels.push_back("BDT");
  colors.push_back(kRed);

  //*****************************************************************************************
  Double_t xmin = 0.0;
  Double_t xmax = 1.0;
  Double_t ymin = 0.0;
  Double_t ymax = 1.0;

  //**************************************
  //Data: V0 vs Best MVAs
  //**************************************
  if (Option == 0 )                              { xmin = 0.0; xmax = 0.2; ymin = 0.0; ymax = 1.0; }
  if (Option == 1 )                              { xmin = 0.0; xmax = 0.2; ymin = 0.0; ymax = 1.0; }
  if (Option == 2 )                              { xmin = 0.0; xmax = 0.2; ymin = 0.0; ymax = 1.0; }
  if (Option == 3 )                              { xmin = 0.05; xmax = 0.25; ymin = 0.75; ymax = 1.0; }
  if (Option == 4 )                              { xmin = 0.05; xmax = 0.25; ymin = 0.65; ymax = 1.0; }
  if (Option == 5 )                              { xmin = 0.0; xmax = 0.3; ymin = 0.0; ymax = 1.0; }

  cv = new TCanvas("cv", "cv", 800, 600);

//    legend = new TLegend(0.45,0.20,0.75,0.50);
  legend = new TLegend(0.44,0.24,0.94,0.44);
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




  TFile *canvasFile = new TFile("ElectronIDMVAPerformancePlots.root","UPDATE");
  canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;
  canvasFile->Close();
  //*****************************************************************************************



  gBenchmark->Show("WWTemplate");       
} 

