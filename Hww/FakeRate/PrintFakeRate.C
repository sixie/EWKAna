//root -l EWKAna/Hww/FakeRate/PrintFakeRate.C+
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

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
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



void PrintMuonFakeRate() {
  
//   TFile *inputFile = new TFile( "MuonFakeRate.SmurfV6.skim.root", "READ");
  TFile *inputFile = new TFile( "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root", "READ");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Int_t OptimizationVariable = 0;
  Int_t PtBin = 1;
  Int_t EtaBin = 1;

  vector<string> FakeRateHistNames;
  FakeRateHistNames.push_back("MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_PtEta");
 
  mithep::TH2DAsymErr *tmp = (mithep::TH2DAsymErr*)inputFile->Get(FakeRateHistNames[0].c_str());
  assert(tmp);
  


//   //For display
//   for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
//     for(int j=1; j<tmp->GetYaxis()->GetNbins() + 1; ++j) {
//       cout << "Pt Bin " << tmp->GetXaxis()->GetBinCenter(i) << " , Eta Bin " << tmp->GetYaxis()->GetBinCenter(j) << " : ";
//       for (int q=0; q < FakeRateHistNames.size() ; ++q) {        
//         cout << tmp->GetBinContent(i,j) << " + " << tmp->GetBinStatErrorHigh(i,j) << " - " << tmp->GetBinStatErrorLow(i,j) << endl;
//       }
//     }
//   }



  //For tex file
  for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
    cout.width(5);
    cout << " $" << tmp->GetXaxis()->GetBinLowEdge(i) << " < p_{T} <= " << tmp->GetXaxis()->GetBinUpEdge(i) << "$" ;
    cout << " & " ;
    for(int j=1; j<tmp->GetYaxis()->GetNbins()+1 ; ++j) {
      cout.width(8);
      char tempNumber[50];
      char tempNumberErrorHigh[50];
      char tempNumberErrorLow[50];
      sprintf(tempNumber, "%.3f",tmp->GetBinContent(i,j));
      sprintf(tempNumberErrorHigh, "%.3f",tmp->GetBinStatErrorHigh(i,j));
      sprintf(tempNumberErrorLow, "%.3f",tmp->GetBinStatErrorLow(i,j));
      cout << "$" << tempNumber << "^{+" << tempNumberErrorHigh << "}_{-" << tempNumberErrorLow << "}$ ";
      if (j != tmp->GetYaxis()->GetNbins()) cout << "& ";
    }
    cout << " \\\\ " << endl;
    cout << " \\hline" << endl;
  }

  
//   TH2D *tmpHist = (TH2D*)tmp->Clone("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");
//   TH2D *tmpHistErrLow = (TH2D*)tmp->Clone("ElectronFakeRateErrorLowDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");
//   TH2D *tmpHistErrHigh = (TH2D*)tmp->Clone("ElectronFakeRateErrorHighDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");

//   for(int i=0; i<tmp->GetXaxis()->GetNbins() + 2; ++i) {
//     for(int j=0; j<tmp->GetYaxis()->GetNbins() + 2; ++j) {
//       tmpHist->SetBinContent(i,tmp->GetBinContent(i,j));
//       tmpHistErrLow->SetBinContent(i,tmp->GetBinStatErrorLow(i,j));
//       tmpHistErrHigh->SetBinContent(i,tmp->GetBinStatErrorHigh(i,j));
//     }
//   }
//   TFile *f = new TFile("ElectronFakeRate_SmurfV5.root" , "RECREATE");
//   f->WriteTObject(tmpHist, tmpHist->GetName(), "WriteDelete");
//   f->WriteTObject(tmpHistErrLow, tmpHistErrLow->GetName(), "WriteDelete");
//   f->WriteTObject(tmpHistErrHigh, tmpHistErrHigh->GetName(), "WriteDelete");
//   f->Close();
//   delete f;




}








void PrintElectronFakeRate() {
  
  TFile *inputFile = new TFile( "ElectronFakeRate.SmurfV6.skim.root", "READ");
  assert(inputFile);
  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Int_t OptimizationVariable = 0;
  Int_t PtBin = 1;
  Int_t EtaBin = 1;

  vector<string> FakeRateHistNames;
  FakeRateHistNames.push_back("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");
 
  mithep::TH2DAsymErr *tmp = (mithep::TH2DAsymErr*)inputFile->Get(FakeRateHistNames[0].c_str());
  assert(tmp);
  


//   //For display
//   for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
//     for(int j=1; j<tmp->GetYaxis()->GetNbins() + 1; ++j) {
//       cout << "Pt Bin " << tmp->GetXaxis()->GetBinCenter(i) << " , Eta Bin " << tmp->GetYaxis()->GetBinCenter(j) << " : ";
//       for (int q=0; q < FakeRateHistNames.size() ; ++q) {        
//         cout << tmp->GetBinContent(i,j) << " + " << tmp->GetBinStatErrorHigh(i,j) << " - " << tmp->GetBinStatErrorLow(i,j) << endl;
//       }
//     }
//   }



  //For tex file
  for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
    cout.width(5);
    cout << " $" << tmp->GetXaxis()->GetBinLowEdge(i) << " < p_{T} <= " << tmp->GetXaxis()->GetBinUpEdge(i) << "$" ;
    cout << " & " ;
    for(int j=1; j<tmp->GetYaxis()->GetNbins()+1 ; ++j) {
      cout.width(8);
      char tempNumber[50];
      char tempNumberErrorHigh[50];
      char tempNumberErrorLow[50];
      sprintf(tempNumber, "%.3f",tmp->GetBinContent(i,j));
      sprintf(tempNumberErrorHigh, "%.3f",tmp->GetBinStatErrorHigh(i,j));
      sprintf(tempNumberErrorLow, "%.3f",tmp->GetBinStatErrorLow(i,j));
      cout << "$" << tempNumber << "^{+" << tempNumberErrorHigh << "}_{-" << tempNumberErrorLow << "}$ ";
      if (j != tmp->GetYaxis()->GetNbins()) cout << "& ";
    }
    cout << " \\\\ " << endl;
    cout << " \\hline" << endl;
  }

  
//   TH2D *tmpHist = (TH2D*)tmp->Clone("ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");
//   TH2D *tmpHistErrLow = (TH2D*)tmp->Clone("ElectronFakeRateErrorLowDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");
//   TH2D *tmpHistErrHigh = (TH2D*)tmp->Clone("ElectronFakeRateErrorHighDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta");

//   for(int i=0; i<tmp->GetXaxis()->GetNbins() + 2; ++i) {
//     for(int j=0; j<tmp->GetYaxis()->GetNbins() + 2; ++j) {
//       tmpHist->SetBinContent(i,tmp->GetBinContent(i,j));
//       tmpHistErrLow->SetBinContent(i,tmp->GetBinStatErrorLow(i,j));
//       tmpHistErrHigh->SetBinContent(i,tmp->GetBinStatErrorHigh(i,j));
//     }
//   }
//   TFile *f = new TFile("ElectronFakeRate_SmurfV5.root" , "RECREATE");
//   f->WriteTObject(tmpHist, tmpHist->GetName(), "WriteDelete");
//   f->WriteTObject(tmpHistErrLow, tmpHistErrLow->GetName(), "WriteDelete");
//   f->WriteTObject(tmpHistErrHigh, tmpHistErrHigh->GetName(), "WriteDelete");
//   f->Close();
//   delete f;

}

void PrintMuonFakeRateFromSmurfHist() {
  
//   TFile *inputFile = new TFile( "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root", "READ");
  TFile *inputFile = new TFile( "/data/smurf/sixie/FakeRates/FakeRates_SmurfV7Muon_BDTGWithIPInfoElectron.root", "READ");

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Int_t OptimizationVariable = 0;
  Int_t PtBin = 1;
  Int_t EtaBin = 1;

  vector<string> FakeRateHistNames;
  FakeRateHistNames.push_back("MuonFakeRate_M2_ptThreshold15_PtEta");
 
  TH2F *tmp = (TH2F*)inputFile->Get(FakeRateHistNames[0].c_str());
  assert(tmp);
  

  //For tex file
  for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
    cout.width(5);
    cout << " $" << tmp->GetXaxis()->GetBinLowEdge(i) << " < p_{T} <= " << tmp->GetXaxis()->GetBinUpEdge(i) << "$" ;
    cout << " & " ;
    for(int j=1; j<tmp->GetYaxis()->GetNbins()+1 ; ++j) {
      cout.width(8);
      char tempNumber[50];
      char tempNumberError[50];
      sprintf(tempNumber, "%.3f",tmp->GetBinContent(i,j));
      sprintf(tempNumberError, "%.3f",tmp->GetBinError(i,j));
      cout << "$" << tempNumber << " +/- " << tempNumberError << "$ ";
      if (j != tmp->GetYaxis()->GetNbins()) cout << "& ";
    }
    cout << " \\\\ " << endl;
    cout << " \\hline" << endl;
  }


}

void PrintElectronFakeRateFromSmurfHist() {
  
  //TFile *inputFile = new TFile( "/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root", "READ");
//   TFile *inputFile = new TFile( "FakeRates_SmurfV6.root", "READ");
//    TFile *inputFile = new TFile( "FakeRates_SmurfLHElectron.root", "READ");
//     TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo.Run2011A.root", "READ");
//    TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.r11bTest.root", "READ");
//    TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo.root", "READ");
//   TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo.CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root", "READ");

//     TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.Run2011A.root.root", "READ");
//     TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.r11b.root.root", "READ");
//      TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.Run2011A.root.CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root", "READ");
//       TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.r11b.root.CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root", "READ");

//      TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.Run2011A.root", "READ");
//      TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.r11b.root", "READ");
//      TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.Run2011A.CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root", "READ");
//      TFile *inputFile = new TFile( "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.r11b.CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root", "READ");

//   TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo.root", "READ");
//      TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo_NVtx0To2.root", "READ"); 
//    TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo_NVtx3To5.root", "READ");
//     TFile *inputFile = new TFile( "FakeRates_Electron_BDTGWithIPInfo_NVtx6To10.root", "READ");




  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Int_t OptimizationVariable = 0;
  Int_t PtBin = 1;
  Int_t EtaBin = 1;

  vector<string> FakeRateHistNames;
  FakeRateHistNames.push_back("ElectronFakeRate_V4_ptThreshold35_PtEta");
 
  TH2F *tmp = (TH2F*)inputFile->Get(FakeRateHistNames[0].c_str());
  assert(tmp);
  

  //For tex file
  for(int i=1; i<tmp->GetXaxis()->GetNbins() + 1; ++i) {
    cout.width(5);
    cout << " $" << tmp->GetXaxis()->GetBinLowEdge(i) << " < p_{T} <= " << tmp->GetXaxis()->GetBinUpEdge(i) << "$" ;
    cout << " & " ;
    for(int j=1; j<tmp->GetYaxis()->GetNbins()+1 ; ++j) {
      cout.width(8);
      char tempNumber[50];
      char tempNumberError[50];
      sprintf(tempNumber, "%.3f",tmp->GetBinContent(i,j));
      sprintf(tempNumberError, "%.3f",tmp->GetBinError(i,j));
      cout << "$" << tempNumber << " +/- " << tempNumberError << "$ ";
      if (j != tmp->GetYaxis()->GetNbins()) cout << "& ";
    }
    cout << " \\\\ " << endl;
    cout << " \\hline" << endl;
  }


}




void PrintFakeRate() {

//   PrintElectronFakeRate();
//   PrintMuonFakeRate();
//   PrintElectronFakeRateFromSmurfHist();
  PrintMuonFakeRateFromSmurfHist();


}


