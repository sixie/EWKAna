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

#endif


void PrintMuonFakeRateFromSmurfHist() {
  
  TFile *inputFile = new TFile( "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/FakeRates_MVAIDIsoCombinedDetIsoSameSigWP.root", "READ");

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
  
  TFile *inputFile = new TFile( "/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/auxiliar/FakeRates_MVAIDIsoCombinedDetIsoSameSigWP.root", "READ");




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

  PrintElectronFakeRateFromSmurfHist();
  PrintMuonFakeRateFromSmurfHist();
  

}


