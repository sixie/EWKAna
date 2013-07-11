//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/MakeElectronLikelihoodPerformancePlots.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt",0)'
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
TGraphAsymmErrors* MakeFlatROCReference() {

  //Make Met Plots
  const UInt_t nPoints = 400;
  double SigEff[nPoints];
  double BkgEff[nPoints];
  double SigEffErrLow[nPoints];
  double SigEffErrHigh[nPoints];
  double BkgEffErrLow[nPoints];
  double BkgEffErrHigh[nPoints];
  double NSigTotal = 0;
  double NBkgTotal = 0;
  
  for(UInt_t b=0; b < nPoints; ++b) {

    SigEff[b] = -2 + (4.0/400.0) * b;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
    BkgEff[b] = -2 + (4.0/400.0) * b;
    BkgEffErrLow[b] = 0;
    BkgEffErrHigh[b] = 0;
  }

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh, SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName("ReferenceROC");
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
TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name, 
                                           Double_t ReferenceEffReal, Double_t ReferenceEffFake ) {

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
    SigEff[b] = ratio/ReferenceEffReal;
//     SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

    n1 = TMath::Nint(nbkg);
    n2 = TMath::Nint(NBkgTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    BkgEff[b] = ratio/ReferenceEffFake;
//     BkgEff[b] = ratio;
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


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void MakeElectronLikelihoodPerformancePlots(string RealElectronFile, string FakeElectronFile, string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************
  vector<Int_t> markers;
  vector<Int_t> colors;
  colors.push_back(kRed);     markers.push_back(20);
  colors.push_back(kBlue);    markers.push_back(21);
  colors.push_back(kMagenta); markers.push_back(22);
  colors.push_back(kCyan);    markers.push_back(34);
  colors.push_back(kBlack);   markers.push_back(29);
  colors.push_back(kGreen);   markers.push_back(33);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *EleIDStandardLikelihood_Real = new TH1F(("EleIDStandardLikelihood_Real"+label).c_str(), "; StandardLikelihood ; Number of Events ",  1000, -20 , 20);
  TH1F *EleIDPFMVA_Real = new TH1F(("EleIDPFMVA_Real"+label).c_str(), "; PFMVA ; Number of Events ",  1000, -2 , 2);

  TH1F *EleIDStandardLikelihood_Fake = new TH1F(("EleIDStandardLikelihood_Fake"+label).c_str(), "; StandardLikelihood ; Number of Events ",  1000, -20 , 20);
  TH1F *EleIDPFMVA_Fake = new TH1F(("EleIDPFMVA_Fake"+label).c_str(), "; PFMVA ; Number of Events ",  1000, -2 , 2);

  Double_t RealElectrons = 0;
  Double_t FakeElectrons = 0;
  Double_t RealElectronPassCutBased = 0;
  Double_t FakeElectronPassCutBased = 0;
  Double_t RealElectronPassVBTF80 = 0;
  Double_t FakeElectronPassVBTF80 = 0;
  Double_t RealElectronPassLHTight = 0;
  Double_t FakeElectronPassLHTight = 0;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
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
  Float_t                 fElePFMVA; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  Float_t                 fEleBDT;
  Float_t                 fEleBDTG;
  Float_t                 fEleNN;
  Float_t                 fEleLikelihood;
  Float_t                 fEleLikelihoodD;

  
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
  RealEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  RealEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "Rho", &fRho); 
  RealEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  RealEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  RealEleTree->SetBranchAddress( "BDTG", &fEleBDTG); 
  RealEleTree->SetBranchAddress( "NN", &fEleNN); 
  RealEleTree->SetBranchAddress( "Likelihood", &fEleLikelihood); 
  RealEleTree->SetBranchAddress( "LikelihoodD", &fEleLikelihoodD); 

  for(UInt_t ientry=0; ientry < RealEleTree->GetEntries(); ientry++) {       	
    RealEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
        
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
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
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply full isolation cut
    Bool_t passIsoCuts = kFALSE;
    if (ptBin == 0) passIsoCuts = ( fElePFIso / fElePt < 0.09 ); 
    if (ptBin == 1) passIsoCuts = ( fElePFIso / fElePt < 0.13 ); 
    if (!passIsoCuts) continue;

    RealElectrons += fWeight;
    Bool_t passNumerator = kTRUE;

    if (passVBTF(1, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP)) 
      RealElectronPassCutBased += fWeight;
    if (passVBTF(0, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP)) 
      RealElectronPassVBTF80 += fWeight;
    if (passLH( fElePt, fEleSCEta, fEleNBrem, fEleStandardLikelihood))
      RealElectronPassLHTight += fWeight;

    //Fill Histograms
    EleIDStandardLikelihood_Real->Fill(fEleStandardLikelihood,fWeight);
    EleIDPFMVA_Real->Fill(TMath::Max(TMath::Min(Double_t(fElePFMVA),Double_t(1.99)),Double_t(-1.99)),fWeight);

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
  FakeEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  FakeEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  FakeEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  FakeEleTree->SetBranchAddress( "BDTG", &fEleBDTG); 
  FakeEleTree->SetBranchAddress( "NN", &fEleNN); 
  FakeEleTree->SetBranchAddress( "Likelihood", &fEleLikelihood); 
  FakeEleTree->SetBranchAddress( "LikelihoodD", &fEleLikelihoodD); 

  for(UInt_t ientry=0; ientry < FakeEleTree->GetEntries(); ientry++) {       	
    FakeEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    if (!(fabs(fEleEta) < 2.5 && fElePt > 10)) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
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
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply full isolation cut
    Bool_t passIsoCuts = kFALSE;
    if (ptBin == 0) passIsoCuts = ( fElePFIso / fElePt < 0.09 ); 
    if (ptBin == 1) passIsoCuts = ( fElePFIso / fElePt < 0.13 ); 
    if (!passIsoCuts) continue;


    FakeElectrons += fWeight;

    if (passVBTF(1, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP)) 
      FakeElectronPassCutBased += fWeight;
    if (passVBTF(0, fElePt, fEleSCEta,fEleEta, fEleSigmaIEtaIEta, fEleDEtaIn, fEleDPhiIn, fEleHoverE,fEleD0,fEleDZ,fEleFBrem, fEleEOverP)) 
      FakeElectronPassVBTF80 += fWeight;
    if (passLH( fElePt, fEleSCEta, fEleNBrem, fEleStandardLikelihood))
      FakeElectronPassLHTight += fWeight;


    //Fill Histograms
    EleIDStandardLikelihood_Fake->Fill(fEleStandardLikelihood,fWeight);
    EleIDPFMVA_Fake->Fill(TMath::Max(TMath::Min(Double_t(fElePFMVA),Double_t(1.99)),Double_t(-1.99)),fWeight);

  } //loop over electrons
  



  
  //*****************************************************************************************
  //Make ROC curves
  //*****************************************************************************************
  TGraphAsymmErrors* ROC_Reference = MakeFlatROCReference();
  Double_t ReferenceEffReal = RealElectronPassVBTF80/RealElectrons;
  Double_t ReferenceEffFake = FakeElectronPassVBTF80/FakeElectrons;

  cout << "VBTF80 Real Electron Efficiency : " << RealElectronPassVBTF80 << " / " << RealElectrons << " = " << RealElectronPassVBTF80/RealElectrons/ReferenceEffReal << endl;
  cout << "VBTF80 Fake Electron Efficiency : " << FakeElectronPassVBTF80 << " / " << FakeElectrons << " = " << FakeElectronPassVBTF80/FakeElectrons/ReferenceEffFake << endl;
  TGraphAsymmErrors* ROC_VBTF80WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassVBTF80/RealElectrons/ReferenceEffReal , FakeElectronPassVBTF80/FakeElectrons/ReferenceEffFake, "ROC_VBTF80WP"+label);


  cout << "Cut-Based Real Electron Efficiency : " << RealElectronPassCutBased << " / " << RealElectrons << " = " << RealElectronPassCutBased/RealElectrons/ReferenceEffReal << endl;
  cout << "Cut-Based Fake Electron Efficiency : " << FakeElectronPassCutBased << " / " << FakeElectrons << " = " << FakeElectronPassCutBased/FakeElectrons/ReferenceEffFake << endl;
  TGraphAsymmErrors* ROC_CutBasedCurrentWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCutBased/RealElectrons/ReferenceEffReal , FakeElectronPassCutBased/FakeElectrons/ReferenceEffFake, "ROC_CutBasedCurrentWP"+label);
 
  cout << "LHTight Real Electron Efficiency : " << RealElectronPassLHTight << " / " << RealElectrons << " = " << RealElectronPassLHTight/RealElectrons/ReferenceEffReal << endl;
  cout << "LHTight Fake Electron Efficiency : " << FakeElectronPassLHTight << " / " << FakeElectrons << " = " << FakeElectronPassLHTight/FakeElectrons/ReferenceEffFake << endl;
  TGraphAsymmErrors* ROC_LHTightWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassLHTight/RealElectrons/ReferenceEffReal , FakeElectronPassLHTight/FakeElectrons/ReferenceEffFake, "ROC_LHTightCurrentWP"+label);

  TGraphAsymmErrors* ROC_StandardLikelihood = MakeSigEffVsBkgEffGraph(EleIDStandardLikelihood_Real, EleIDStandardLikelihood_Fake, "ROC_StandardLikelihood"+label, ReferenceEffReal, ReferenceEffFake );
  TGraphAsymmErrors* ROC_PFMVA = MakeSigEffVsBkgEffGraph(EleIDPFMVA_Real, EleIDPFMVA_Fake, "ROC_PFMVA"+label, ReferenceEffReal, ReferenceEffFake );


//   //***
//   //Don't Divide by VBTF80 Eff

//   TGraphAsymmErrors* ROC_Reference = MakeFlatROCReference();
//   Double_t ReferenceEffReal = RealElectronPassVBTF80;
//   Double_t ReferenceEffFake = FakeElectronPassVBTF80;

//   cout << "VBTF80 Real Electron Efficiency : " << RealElectronPassVBTF80 << " / " << RealElectrons << " = " << RealElectronPassVBTF80/RealElectrons << endl;
//   cout << "VBTF80 Fake Electron Efficiency : " << FakeElectronPassVBTF80 << " / " << FakeElectrons << " = " << FakeElectronPassVBTF80/FakeElectrons << endl;
//   TGraphAsymmErrors* ROC_VBTF80WP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassVBTF80/RealElectrons , FakeElectronPassVBTF80/FakeElectrons, "ROC_VBTF80WP"+label);


//   cout << "Cut-Based Real Electron Efficiency : " << RealElectronPassCutBased << " / " << RealElectrons << " = " << RealElectronPassCutBased/RealElectrons << endl;
//   cout << "Cut-Based Fake Electron Efficiency : " << FakeElectronPassCutBased << " / " << FakeElectrons << " = " << FakeElectronPassCutBased/FakeElectrons << endl;
//   TGraphAsymmErrors* ROC_CutBasedCurrentWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassCutBased/RealElectrons , FakeElectronPassCutBased/FakeElectrons, "ROC_CutBasedCurrentWP"+label);

//   cout << "LHTight Real Electron Efficiency : " << RealElectronPassLHTight << " / " << RealElectrons << " = " << RealElectronPassLHTight/RealElectrons << endl;
//   cout << "LHTight Fake Electron Efficiency : " << FakeElectronPassLHTight << " / " << FakeElectrons << " = " << FakeElectronPassLHTight/FakeElectrons << endl;
//   TGraphAsymmErrors* ROC_LHTightWP = MakeCurrentWPSigEffVsBkgEffGraph(RealElectronPassLHTight/RealElectrons , FakeElectronPassLHTight/FakeElectrons, "ROC_LHTightCurrentWP"+label);
 

//   TGraphAsymmErrors* ROC_StandardLikelihood = MakeSigEffVsBkgEffGraph(EleIDStandardLikelihood_Real, EleIDStandardLikelihood_Fake, "ROC_StandardLikelihood"+label, ReferenceEffReal, ReferenceEffFake );
//   TGraphAsymmErrors* ROC_PFMVA = MakeSigEffVsBkgEffGraph(EleIDPFMVA_Real, EleIDPFMVA_Fake, "ROC_PFMVA"+label, ReferenceEffReal, ReferenceEffFake );



  //*****************************************************************************************
  //*****************************************************************************************
  TFile *canvasFile = new TFile("ElectronIDMVAPerformancePlots.root","UPDATE");
  TLegend* legend;
  TCanvas* cv;
  string plotname;



  //*****************************************************************************************
  //Plot Distributions
  //*****************************************************************************************
  NormalizeHist(EleIDStandardLikelihood_Real);
  NormalizeHist(EleIDStandardLikelihood_Fake);


  vector<TH1F*> hists;
  vector<string> histLabels;


  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleIDStandardLikelihood_Real, "Real Electrons", "L");
  legend->AddEntry(EleIDStandardLikelihood_Fake, "Fake Electrons", "L");

  EleIDStandardLikelihood_Real->SetLineColor(kBlue);
  EleIDStandardLikelihood_Fake->SetLineColor(kRed);   
  EleIDStandardLikelihood_Real->Draw("hist");
  EleIDStandardLikelihood_Fake->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleIDStandardLikelihood_" + label + ".gif").c_str());


  //*****************************************************************************************
  //Plot ROC Curves
  //*****************************************************************************************
  vector<TGraphAsymmErrors*> ROCGraphs;
  vector<string> GraphLabels;

  //*****************************************************************************************
  //*****************************************************************************************
  ROCGraphs.clear();
  GraphLabels.clear();
  plotname = "ElectronIDMVA"+label;

  ROCGraphs.push_back(ROC_Reference);
  GraphLabels.push_back("Reference Line");

  ROCGraphs.push_back(ROC_StandardLikelihood);
  GraphLabels.push_back("StandardLikelihood");
  ROCGraphs.push_back(ROC_PFMVA);
  GraphLabels.push_back("PFMVA");

  //*****************************************************************************************
  Double_t xmin, xmax, ymin, ymax;
  if (Option == 0 ||Option == 1 ||  Option == 2) { xmin = 0.0; xmax = 2.0; ymin = 0.0; ymax = 2.0; }
  if (Option == 3 )                              { xmin = 0.0; xmax = 2.0; ymin = 0.0; ymax = 2.0; }
  if (Option == 4 ||  Option == 5)               { xmin = 0.0; xmax = 2.0; ymin = 0.0; ymax = 2.0; }

  cv = new TCanvas("cv", "cv", 800, 600);

  legend = new TLegend(0.63,0.20,0.93,0.50);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
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

  legend->AddEntry(ROC_VBTF80WP, "VBTF80", "P");
  ROC_VBTF80WP->SetFillColor(kRed);
  ROC_VBTF80WP->SetMarkerColor(kRed);
  ROC_VBTF80WP->SetMarkerStyle(34);
  ROC_VBTF80WP->SetMarkerSize(2.5);
  ROC_VBTF80WP->Draw("Psame");

  legend->AddEntry(ROC_LHTightWP, "NEW LH WP", "P");
  ROC_LHTightWP->SetFillColor(kBlue);
  ROC_LHTightWP->SetMarkerColor(kBlue);
  ROC_LHTightWP->SetMarkerStyle(34);
  ROC_LHTightWP->SetMarkerSize(2.5);
  ROC_LHTightWP->Draw("Psame");
  
  legend->AddEntry(ROC_CutBasedCurrentWP, "CurrentWP", "P");
  ROC_CutBasedCurrentWP->SetFillColor(kGreen);
  ROC_CutBasedCurrentWP->SetMarkerColor(kGreen);
  ROC_CutBasedCurrentWP->SetMarkerStyle(34);
  ROC_CutBasedCurrentWP->SetMarkerSize(2.5);
  ROC_CutBasedCurrentWP->Draw("Psame");

  legend->Draw();
  
  cv->SaveAs(("ROCGraphs_" + plotname + ".gif").c_str());
  canvasFile->WriteTObject(cv,("ROCGraphs_" + plotname).c_str(), "WriteDelete") ;


  //*****************************************************************************************
  canvasFile->Close();
  //*****************************************************************************************



//   //*****************************************************************************************
//   //Save Histograms in file
//   //*****************************************************************************************
//   TFile *file = new TFile("ElectronIDMVAResults.root", "UPDATE");
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
 
//   file->Close();
//   delete file;




  gBenchmark->Show("WWTemplate");       
} 

