//root -l EWKAna/Hww/LeptonSelection/CompareElectronSelectionVariables.C+\(\)
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
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <TLegend.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

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
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2);
void runMakeElectronPlots(const string dataInputFilename,       
                          const string Label, Bool_t isMuonData);
void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname);

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//=== MAIN MACRO =================================================================================================

void CompareElectronSelectionVariables() 
{  
  
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/cern/filefi/018/HwwAnalysis_w10-h130ww2l-gf-z2-v8-pu11_noskim_0000.root","signal_lowPt",kFALSE);


//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12-pu_noskim_normalized.root","signalHww130_Pt10To15",kFALSE);
//    runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12-pu_noskim_normalized.root","signalHww130_Pt15To20",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12-pu_noskim_normalized.root","signalHww130_Pt20ToInf",kFALSE);




//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12-pu_noskim_normalized.root","signalHww130_highPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/HwwAnalysis_w10-zee-powheg-c10-v8-pu11_noskim.root","Zee_MC_highPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/WWAnalysis/WWAnalysis_mu_noskim.root","bkg_lowPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/WWAnalysis/WWAnalysis_mu_noskim.root","bkg_highPt",kFALSE);


//   makePlot("sigmaIPhiIPhi_Barrel_signalHww130_Pt15To20", "Signal", "sigmaIPhiIPhi_Barrel_bkg_Pt15To20", "Bkg", "sigmaIPhiIPhi_Barrel_Pt15To20SignalVsBkg");
//   makePlot("sigmaIPhiIPhi_Endcap_signalHww130_Pt15To20", "Signal", "sigmaIPhiIPhi_Endcap_bkg_Pt15To20", "Bkg", "sigmaIPhiIPhi_Endcap_Pt15To20SignalVsBkg");
  makePlot("sigmaIPhiIPhi_Barrel_signalHww130_Pt20ToInf", "Signal", "sigmaIPhiIPhi_Barrel_bkg_Pt20ToInf", "Bkg", "sigmaIPhiIPhi_Barrel_Pt20ToInfSignalVsBkg");
  makePlot("sigmaIPhiIPhi_Endcap_signalHww130_Pt20ToInf", "Signal", "sigmaIPhiIPhi_Endcap_bkg_Pt20ToInf", "Bkg", "sigmaIPhiIPhi_Endcap_Pt20ToInfSignalVsBkg");




  //*******************************************************************************
  //Compare Pt10To15 Signal Vs Bkg  
  //*******************************************************************************
//   makePlot("relIso_Barrel_signalHww130_Pt10To15", "Signal", "relIso_Barrel_bkg_Pt10To15", "Bkg", "relIso_Barrel_Pt10To15SignalVsBkg");
//   makePlot("relIso04_Barrel_signalHww130_Pt10To15", "Signal", "relIso04_Barrel_bkg_Pt10To15", "Bkg", "relIso04_Barrel_Pt10To15SignalVsBkg");
//   makePlot("relIso_Endcap_signalHww130_Pt10To15", "signal", "relIso_Endcap_bkg_Pt10To15", "bkg", "relIso_Endcap_Pt10To15SignalVsBkg");
//   makePlot("relIso04_Endcap_signalHww130_Pt10To15", "signal", "relIso04_Endcap_bkg_Pt10To15", "bkg", "relIso04_Endcap_Pt10To15SignalVsBkg");

//   makePlot("relIso_Barrel_signalHww130_Pt15To20", "Signal", "relIso_Barrel_bkg_Pt15To20", "Bkg", "relIso_Barrel_Pt15To20SignalVsBkg");
//   makePlot("relIso04_Barrel_signalHww130_Pt15To20", "Signal", "relIso04_Barrel_bkg_Pt15To20", "Bkg", "relIso04_Barrel_Pt15To20SignalVsBkg");
//   makePlot("relIso_Endcap_signalHww130_Pt15To20", "signal", "relIso_Endcap_bkg_Pt15To20", "bkg", "relIso_Endcap_Pt15To20SignalVsBkg");
//   makePlot("relIso04_Endcap_signalHww130_Pt15To20", "signal", "relIso04_Endcap_bkg_Pt15To20", "bkg", "relIso04_Endcap_Pt15To20SignalVsBkg");


//     makePlot("relIso_Barrel_signalHww130_Pt10To15", "Signal", "relIso_Barrel_bkg_Pt10To15", "Bkg", "relIso_Barrel_Pt10To15SignalVsBkg");
//     makePlot("relIso04_Barrel_signalHww130_Pt10To15", "Signal", "relIso04_Barrel_bkg_Pt10To15", "Bkg", "relIso04_Barrel_Pt10To15SignalVsBkg");
//     makePlot("ecalIso03_Barrel_signalHww130_Pt10To15", "Signal", "ecalIso03_Barrel_bkg_Pt10To15", "Bkg", "ecalIso03_Barrel_Pt10To15SignalVsBkg");
//     makePlot("ecalIso04_Barrel_signalHww130_Pt10To15", "Signal", "ecalIso04_Barrel_bkg_Pt10To15", "Bkg", "ecalIso04_Barrel_Pt10To15SignalVsBkg");
//     makePlot("hcalIso03_Barrel_signalHww130_Pt10To15", "Signal", "hcalIso03_Barrel_bkg_Pt10To15", "Bkg", "hcalIso03_Barrel_Pt10To15SignalVsBkg");
//     makePlot("hcalIso04_Barrel_signalHww130_Pt10To15", "Signal", "hcalIso04_Barrel_bkg_Pt10To15", "Bkg", "hcalIso04_Barrel_Pt10To15SignalVsBkg");

//     makePlot("relIso_Barrel_signalHww130_Pt20ToInf", "Signal", "relIso_Barrel_bkg_Pt20ToInf", "Bkg", "relIso_Barrel_Pt20ToInfSignalVsBkg");
//     makePlot("relIso04_Barrel_signalHww130_Pt20ToInf", "Signal", "relIso04_Barrel_bkg_Pt20ToInf", "Bkg", "relIso04_Barrel_Pt20ToInfSignalVsBkg");
//     makePlot("ecalIso03_Barrel_signalHww130_Pt20ToInf", "Signal", "ecalIso03_Barrel_bkg_Pt20ToInf", "Bkg", "ecalIso03_Barrel_Pt20ToInfSignalVsBkg");
//     makePlot("ecalIso04_Barrel_signalHww130_Pt20ToInf", "Signal", "ecalIso04_Barrel_bkg_Pt20ToInf", "Bkg", "ecalIso04_Barrel_Pt20ToInfSignalVsBkg");
//     makePlot("hcalIso03_Barrel_signalHww130_Pt20ToInf", "Signal", "hcalIso03_Barrel_bkg_Pt20ToInf", "Bkg", "hcalIso03_Barrel_Pt20ToInfSignalVsBkg");
//     makePlot("hcalIso04_Barrel_signalHww130_Pt20ToInf", "Signal", "hcalIso04_Barrel_bkg_Pt20ToInf", "Bkg", "hcalIso04_Barrel_Pt20ToInfSignalVsBkg");



//     makePlot("relIso04_Endcap_signalHww130_Pt20ToInf", "signal", "relIso04_Endcap_bkg_Pt20ToInf", "bkg", "relIso04_Endcap_Pt20ToInfSignalVsBkg");

//    makePlot("likelihood_0brem_Barrel_signalHww130_Pt10To15", "signal", "likelihood_0brem_Barrel_bkg_Pt10To15", "bkg", "likelihood_0brem_Barrel_Pt10To15SignalVsBkg");
//    makePlot("likelihood_0brem_Barrel_signalHww130_Pt15To20", "signal", "likelihood_0brem_Barrel_bkg_Pt15To20", "bkg", "likelihood_0brem_Barrel_Pt15To20SignalVsBkg");

//    makePlot("likelihood_0brem_Endcap_signalHww130_Pt10To15", "signal", "likelihood_0brem_Endcap_bkg_Pt10To15", "bkg", "likelihood_0brem_Endcap_Pt10To15SignalVsBkg");
//    makePlot("likelihood_0brem_Endcap_signalHww130_Pt15To20", "signal", "likelihood_0brem_Endcap_bkg_Pt15To20", "bkg", "likelihood_0brem_Endcap_Pt15To20SignalVsBkg");

//    makePlot("likelihood_1brem_Barrel_signalHww130_Pt10To15", "signal", "likelihood_1brem_Barrel_bkg_Pt10To15", "bkg", "likelihood_1brem_Barrel_Pt10To15SignalVsBkg");
//    makePlot("likelihood_1brem_Barrel_signalHww130_Pt15To20", "signal", "likelihood_1brem_Barrel_bkg_Pt15To20", "bkg", "likelihood_1brem_Barrel_Pt15To20SignalVsBkg");

//    makePlot("likelihood_1brem_Endcap_signalHww130_Pt10To15", "signal", "likelihood_1brem_Endcap_bkg_Pt10To15", "bkg", "likelihood_1brem_Endcap_Pt10To15SignalVsBkg");
//    makePlot("likelihood_1brem_Endcap_signalHww130_Pt15To20", "signal", "likelihood_1brem_Endcap_bkg_Pt15To20", "bkg", "likelihood_1brem_Endcap_Pt15To20SignalVsBkg");


  //*******************************************************************************
  //Compare LowPt Signal Vs Bkg  
  //*******************************************************************************


//   makePlot("sigmaIEtaIEta_Barrel_signal_lowPt", "signal", "sigmaIEtaIEta_Barrel_bkg_lowPt", "bkg", "sigmaIEtaIEta_Barrel_LowPtSignalVsBkg");
//   makePlot("DeltaEta_Barrel_signal_lowPt", "signal", "DeltaEta_Barrel_bkg_lowPt", "bkg", "DeltaEta_Barrel_LowPtSignalVsBkg"); 
//   makePlot("DeltaPhi_Barrel_signal_lowPt", "signal", "DeltaPhi_Barrel_bkg_lowPt", "bkg", "DeltaPhi_Barrel_LowPtSignalVsBkg");
//   makePlot("HOverE_Barrel_signal_lowPt", "signal", "HOverE_Barrel_bkg_lowPt", "bkg", "HOverE_Barrel_LowPtSignalVsBkg");
//   makePlot("relIso_Barrel_signal_lowPt", "signal", "relIso_Barrel_bkg_lowPt", "bkg", "relIso_Barrel_LowPtSignalVsBkg");

//   makePlot("sigmaIEtaIEta_Endcap_signal_lowPt", "signal", "sigmaIEtaIEta_Endcap_bkg_lowPt", "bkg", "sigmaIEtaIEta_Endcap_LowPtSignalVsBkg");
//   makePlot("DeltaEta_Endcap_signal_lowPt", "signal", "DeltaEta_Endcap_bkg_lowPt", "bkg", "DeltaEta_Endcap_LowPtSignalVsBkg");
//   makePlot("DeltaPhi_Endcap_signal_lowPt", "signal", "DeltaPhi_Endcap_bkg_lowPt", "bkg", "DeltaPhi_Endcap_LowPtSignalVsBkg");
//   makePlot("HOverE_Endcap_signal_lowPt", "signal", "HOverE_Endcap_bkg_lowPt", "bkg", "HOverE_Endcap_LowPtSignalVsBkg");
//   makePlot("relIso_Endcap_signal_lowPt", "signal", "relIso_Endcap_bkg_lowPt", "bkg", "relIso_Endcap_LowPtSignalVsBkg");

  //*******************************************************************************
  //Compare Signal LowPt Vs HighPt
  //*******************************************************************************


//   makePlot("sigmaIEtaIEta_Barrel_signal_lowPt", "lowPt", "sigmaIEtaIEta_Barrel_signal_highPt", "highPt",  "sigmaIEtaIEta_Barrel_SignalLowPtVsHighPt");
//   makePlot("DeltaEta_Barrel_signal_lowPt", "lowPt", "DeltaEta_Barrel_signal_highPt", "highPt", "DeltaEta_Barrel_SignalLowPtVsHighPt"); 
//   makePlot("DeltaPhi_Barrel_signal_lowPt", "lowPt", "DeltaPhi_Barrel_signal_highPt", "highPt", "DeltaPhi_Barrel_SignalLowPtVsHighPt");
//   makePlot("HOverE_Barrel_signal_lowPt", "lowPt", "HOverE_Barrel_signal_highPt", "highPt", "HOverE_Barrel_SignalLowPtVsHighPt");
//   makePlot("relIso_Barrel_signal_lowPt", "lowPt", "relIso_Barrel_signal_highPt", "highPt", "relIso_Barrel_SignalLowPtVsHighPt");

//   makePlot("sigmaIEtaIEta_Endcap_signal_lowPt", "lowPt", "sigmaIEtaIEta_Endcap_signal_highPt", "highPt", "sigmaIEtaIEta_Endcap_SignalLowPtVsHighPt");
//   makePlot("DeltaEta_Endcap_signal_lowPt", "lowPt", "DeltaEta_Endcap_signal_highPt", "highPt", "DeltaEta_Endcap_SignalLowPtVsHighPt");
//   makePlot("DeltaPhi_Endcap_signal_lowPt", "lowPt", "DeltaPhi_Endcap_signal_highPt", "highPt", "DeltaPhi_Endcap_SignalLowPtVsHighPt");
//   makePlot("HOverE_Endcap_signal_lowPt", "lowPt", "HOverE_Endcap_signal_highPt", "highPt", "HOverE_Endcap_SignalLowPtVsHighPt");
//   makePlot("relIso_Endcap_signal_lowPt", "lowPt", "relIso_Endcap_signal_highPt", "highPt", "relIso_Endcap_SignalLowPtVsHighPt");


  //*******************************************************************************
  //Compare Bkg LowPt Vs HighPt
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_bkg_lowPt", "lowPt", "sigmaIEtaIEta_Barrel_bkg_highPt", "highPt", "sigmaIEtaIEta_Barrel_BkgLowPtVsHighPt");
//   makePlot("DeltaEta_Barrel_bkg_lowPt", "lowPt", "DeltaEta_Barrel_bkg_highPt", "highPt", "DeltaEta_Barrel_BkgLowPtVsHighPt"); 
//   makePlot("DeltaPhi_Barrel_bkg_lowPt", "lowPt", "DeltaPhi_Barrel_bkg_highPt", "highPt", "DeltaPhi_Barrel_BkgLowPtVsHighPt");
//   makePlot("HOverE_Barrel_bkg_lowPt", "lowPt", "HOverE_Barrel_bkg_highPt", "highPt", "HOverE_Barrel_BkgLowPtVsHighPt");
//   makePlot("relIso_Barrel_bkg_lowPt", "lowPt", "relIso_Barrel_bkg_highPt", "highPt", "relIso_Barrel_BkgLowPtVsHighPt");

//   makePlot("sigmaIEtaIEta_Endcap_bkg_lowPt", "lowPt", "sigmaIEtaIEta_Endcap_bkg_highPt", "highPt", "sigmaIEtaIEta_Endcap_BkgLowPtVsHighPt");
//   makePlot("DeltaEta_Endcap_bkg_lowPt", "lowPt", "DeltaEta_Endcap_bkg_highPt", "highPt", "DeltaEta_Endcap_BkgLowPtVsHighPt");
//   makePlot("DeltaPhi_Endcap_bkg_lowPt", "lowPt", "DeltaPhi_Endcap_bkg_highPt", "highPt", "DeltaPhi_Endcap_BkgLowPtVsHighPt");
//   makePlot("HOverE_Endcap_bkg_lowPt", "lowPt", "HOverE_Endcap_bkg_highPt", "highPt", "HOverE_Endcap_BkgLowPtVsHighPt");
//   makePlot("relIso_Endcap_bkg_lowPt", "lowPt", "relIso_Endcap_bkg_highPt", "highPt", "relIso_Endcap_BkgLowPtVsHighPt");



  //*******************************************************************************
  //Compare MC Hww Vs MC Zee
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_signalHww130_highPt", "HwwMC", "sigmaIEtaIEta_Barrel_Zee_MC_highPt", "ZeeMC",  "sigmaIEtaIEta_Barrel_HighPtHwwMCVsZeeMC");
//   makePlot("DeltaEta_Barrel_signalHww130_highPt", "HwwMC", "DeltaEta_Barrel_Zee_MC_highPt", "ZeeMC", "DeltaEta_Barrel_HighPtHwwMCVsZeeMC"); 
//   makePlot("DeltaPhi_Barrel_signalHww130_highPt", "HwwMC", "DeltaPhi_Barrel_Zee_MC_highPt", "ZeeMC", "DeltaPhi_Barrel_HighPtHwwMCVsZeeMC");
//   makePlot("HOverE_Barrel_signalHww130_highPt", "HwwMC", "HOverE_Barrel_Zee_MC_highPt", "ZeeMC", "HOverE_Barrel_HighPtHwwMCVsZeeMC");
//   makePlot("relIso_Barrel_signalHww130_highPt", "HwwMC", "relIso_Barrel_Zee_MC_highPt", "ZeeMC", "relIso_Barrel_HighPtHwwMCVsZeeMC");

//   makePlot("sigmaIEtaIEta_Endcap_signalHww130_highPt", "HwwMC", "sigmaIEtaIEta_Endcap_Zee_MC_highPt", "ZeeMC", "sigmaIEtaIEta_Endcap_HighPtHwwMCVsZeeMC");
//   makePlot("DeltaEta_Endcap_signalHww130_highPt", "HwwMC", "DeltaEta_Endcap_Zee_MC_highPt", "ZeeMC", "DeltaEta_Endcap_HighPtHwwMCVsZeeMC");
//   makePlot("DeltaPhi_Endcap_signalHww130_highPt", "HwwMC", "DeltaPhi_Endcap_Zee_MC_highPt", "ZeeMC", "DeltaPhi_Endcap_HighPtHwwMCVsZeeMC");
//   makePlot("HOverE_Endcap_signalHww130_highPt", "HwwMC", "HOverE_Endcap_Zee_MC_highPt", "ZeeMC", "HOverE_Endcap_HighPtHwwMCVsZeeMC");
//   makePlot("relIso_Endcap_signalHww130_highPt", "HwwMC", "relIso_Endcap_Zee_MC_highPt", "ZeeMC", "relIso_Endcap_HighPtHwwMCVsZeeMC");


  //*******************************************************************************
  //Compare MC Zee Vs MC TagAndProbe Zee
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "sigmaIEtaIEta_Barrel_Zee_MC_highPt", "ZeeMC",  "sigmaIEtaIEta_Barrel_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("DeltaEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaEta_Barrel_Zee_MC_highPt", "ZeeMC", "DeltaEta_Barrel_HighPtZeeTagAndProbeMCVsZeeMC"); 
//   makePlot("DeltaPhi_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaPhi_Barrel_Zee_MC_highPt", "ZeeMC", "DeltaPhi_Barrel_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("HOverE_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "HOverE_Barrel_Zee_MC_highPt", "ZeeMC", "HOverE_Barrel_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("relIso_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "relIso_Barrel_Zee_MC_highPt", "ZeeMC", "relIso_Barrel_HighPtZeeTagAndProbeMCVsZeeMC");

//   makePlot("sigmaIEtaIEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "sigmaIEtaIEta_Endcap_Zee_MC_highPt", "ZeeMC", "sigmaIEtaIEta_Endcap_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("DeltaEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaEta_Endcap_Zee_MC_highPt", "ZeeMC", "DeltaEta_Endcap_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("DeltaPhi_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaPhi_Endcap_Zee_MC_highPt", "ZeeMC", "DeltaPhi_Endcap_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("HOverE_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "HOverE_Endcap_Zee_MC_highPt", "ZeeMC", "HOverE_Endcap_HighPtZeeTagAndProbeMCVsZeeMC");
//   makePlot("relIso_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "relIso_Endcap_Zee_MC_highPt", "ZeeMC", "relIso_Endcap_HighPtZeeTagAndProbeMCVsZeeMC");

  //*******************************************************************************
  //Compare PileUp11 Vs Fall10 Pileup
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "sigmaIEtaIEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU",  "sigmaIEtaIEta_Barrel_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("DeltaEta_Barrel_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "DeltaEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "DeltaEta_Barrel_HighPtZeeMCTPFall10PUVsPU11"); 
//   makePlot("DeltaPhi_Barrel_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "DeltaPhi_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "DeltaPhi_Barrel_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("HOverE_Barrel_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "HOverE_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "HOverE_Barrel_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("relIso_Barrel_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "relIso_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "relIso_Barrel_HighPtZeeMCTPFall10PUVsPU11");

//   makePlot("sigmaIEtaIEta_Endcap_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "sigmaIEtaIEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "sigmaIEtaIEta_Endcap_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("DeltaEta_Endcap_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "DeltaEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "DeltaEta_Endcap_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("DeltaPhi_Endcap_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "DeltaPhi_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "DeltaPhi_Endcap_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("HOverE_Endcap_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "HOverE_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "HOverE_Endcap_HighPtZeeMCTPFall10PUVsPU11");
//   makePlot("relIso_Endcap_Zee_MCTagAndProbePU11_highPt", "ZeeMCTP_PU11", "relIso_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTPFall10PU", "relIso_Endcap_HighPtZeeMCTPFall10PUVsPU11");


  //*******************************************************************************
  //Compare DataZee Vs HwwMC
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_DataZee_highPt", "DataZee", "sigmaIEtaIEta_Barrel_signal_highPt", "HwwMC",  "sigmaIEtaIEta_Barrel_HighPtDataZeeVsHwwMC");
//   makePlot("DeltaEta_Barrel_DataZee_highPt", "DataZee", "DeltaEta_Barrel_signal_highPt", "HwwMC", "DeltaEta_Barrel_HighPtDataZeeVsHwwMC"); 
//   makePlot("DeltaPhi_Barrel_DataZee_highPt", "DataZee", "DeltaPhi_Barrel_signal_highPt", "HwwMC", "DeltaPhi_Barrel_HighPtDataZeeVsHwwMC");
//   makePlot("HOverE_Barrel_DataZee_highPt", "DataZee", "HOverE_Barrel_signal_highPt", "HwwMC", "HOverE_Barrel_HighPtDataZeeVsHwwMC");
//   makePlot("relIso_Barrel_DataZee_highPt", "DataZee", "relIso_Barrel_signal_highPt", "HwwMC", "relIso_Barrel_HighPtDataZeeVsHwwMC");

//   makePlot("sigmaIEtaIEta_Endcap_DataZee_highPt", "DataZee", "sigmaIEtaIEta_Endcap_signal_highPt", "HwwMC", "sigmaIEtaIEta_Endcap_HighPtDataZeeVsHwwMC");
//   makePlot("DeltaEta_Endcap_DataZee_highPt", "DataZee", "DeltaEta_Endcap_signal_highPt", "HwwMC", "DeltaEta_Endcap_HighPtDataZeeVsHwwMC");
//   makePlot("DeltaPhi_Endcap_DataZee_highPt", "DataZee", "DeltaPhi_Endcap_signal_highPt", "HwwMC", "DeltaPhi_Endcap_HighPtDataZeeVsHwwMC");
//   makePlot("HOverE_Endcap_DataZee_highPt", "DataZee", "HOverE_Endcap_signal_highPt", "HwwMC", "HOverE_Endcap_HighPtDataZeeVsHwwMC");
//   makePlot("relIso_Endcap_DataZee_highPt", "DataZee", "relIso_Endcap_signal_highPt", "HwwMC", "relIso_Endcap_HighPtDataZeeVsHwwMC");

  //*******************************************************************************
  //Compare DataZee Vs MC Zee
  //*******************************************************************************

//   makePlot("sigmaIEtaIEta_Barrel_Zee_Data_highPt", "DataZee", "sigmaIEtaIEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP",  "sigmaIEtaIEta_Barrel_HighPtDataZeeVsZeeMCTP");
//   makePlot("DeltaEta_Barrel_Zee_Data_highPt", "DataZee", "DeltaEta_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaEta_Barrel_HighPtDataZeeVsZeeMCTP"); 
//   makePlot("DeltaPhi_Barrel_Zee_Data_highPt", "DataZee", "DeltaPhi_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaPhi_Barrel_HighPtDataZeeVsZeeMCTP");
//   makePlot("HOverE_Barrel_Zee_Data_highPt", "DataZee", "HOverE_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "HOverE_Barrel_HighPtDataZeeVsZeeMCTP");
//   makePlot("relIso_Barrel_Zee_Data_highPt", "DataZee", "relIso_Barrel_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "relIso_Barrel_HighPtDataZeeVsZeeMCTP");

//   makePlot("sigmaIEtaIEta_Endcap_Zee_Data_highPt", "DataZee", "sigmaIEtaIEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "sigmaIEtaIEta_Endcap_HighPtDataZeeVsZeeMCTP");
//   makePlot("DeltaEta_Endcap_Zee_Data_highPt", "DataZee", "DeltaEta_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaEta_Endcap_HighPtDataZeeVsZeeMCTP");
//   makePlot("DeltaPhi_Endcap_Zee_Data_highPt", "DataZee", "DeltaPhi_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "DeltaPhi_Endcap_HighPtDataZeeVsZeeMCTP");
//   makePlot("HOverE_Endcap_Zee_Data_highPt", "DataZee", "HOverE_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "HOverE_Endcap_HighPtDataZeeVsZeeMCTP");
//   makePlot("relIso_Endcap_Zee_Data_highPt", "DataZee", "relIso_Endcap_Zee_MCTagAndProbe_highPt", "ZeeMCTP", "relIso_Endcap_HighPtDataZeeVsZeeMCTP");

//   makePlot("sigmaIEtaIEta_Barrel_Zee_Data_lowPt", "DataZee", "sigmaIEtaIEta_Barrel_Zee_MCTagAndProbe_lowPt", "ZeeMCTP",  "sigmaIEtaIEta_Barrel_LowPtDataZeeVsZeeMCTP");
//   makePlot("DeltaEta_Barrel_Zee_Data_lowPt", "DataZee", "DeltaEta_Barrel_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "DeltaEta_Barrel_LowPtDataZeeVsZeeMCTP"); 
//   makePlot("DeltaPhi_Barrel_Zee_Data_lowPt", "DataZee", "DeltaPhi_Barrel_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "DeltaPhi_Barrel_LowPtDataZeeVsZeeMCTP");
//   makePlot("HOverE_Barrel_Zee_Data_lowPt", "DataZee", "HOverE_Barrel_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "HOverE_Barrel_LowPtDataZeeVsZeeMCTP");
//   makePlot("relIso_Barrel_Zee_Data_lowPt", "DataZee", "relIso_Barrel_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "relIso_Barrel_LowPtDataZeeVsZeeMCTP");

//   makePlot("sigmaIEtaIEta_Endcap_Zee_Data_lowPt", "DataZee", "sigmaIEtaIEta_Endcap_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "sigmaIEtaIEta_Endcap_LowPtDataZeeVsZeeMCTP");
//   makePlot("DeltaEta_Endcap_Zee_Data_lowPt", "DataZee", "DeltaEta_Endcap_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "DeltaEta_Endcap_LowPtDataZeeVsZeeMCTP");
//   makePlot("DeltaPhi_Endcap_Zee_Data_lowPt", "DataZee", "DeltaPhi_Endcap_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "DeltaPhi_Endcap_LowPtDataZeeVsZeeMCTP");
//   makePlot("HOverE_Endcap_Zee_Data_lowPt", "DataZee", "HOverE_Endcap_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "HOverE_Endcap_LowPtDataZeeVsZeeMCTP");
//   makePlot("relIso_Endcap_Zee_Data_lowPt", "DataZee", "relIso_Endcap_Zee_MCTagAndProbe_lowPt", "ZeeMCTP", "relIso_Endcap_LowPtDataZeeVsZeeMCTP");


}



void runMakeElectronPlots(const string dataInputFilename,       
                const string Label, Bool_t isMuonData) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    
 // configuring electron likelihood
  TFile *fileLH = TFile::Open("MitPhysics/data/ElectronLikelihoodPdfs_MC.root");
  TDirectory *EBlt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EBgt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");
  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useSigmaPhiPhi = true;
  ElectronLikelihood *LH = new ElectronLikelihood(&(*EBlt15dir), &(*EElt15dir), &(*EBgt15dir), &(*EEgt15dir),
                              defaultSwitches, std::string("class"),std::string("class"),true,true);



  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_ee = new TH1D("dileptonMass_ee", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_emu = new TH1D("dileptonMass_emu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_mumu = new TH1D("dileptonMass_mumu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  Int_t Count_ee_BB = 0;
  Int_t Count_ee_BE = 0;
  Int_t Count_ee_EE = 0;
  Int_t Count_mm_BB = 0;
  Int_t Count_mm_BE = 0;
  Int_t Count_mm_EE = 0;
  Int_t Count_em_BB = 0;
  Int_t Count_em_BE = 0;
  Int_t Count_em_EE = 0;
  ofstream eventListFile("eventList.txt");

  TH1D *ptMin = new TH1D("ptMin", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  TH1D *ptMax = new TH1D("ptMax", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  string tmplabel = Label;
  if (tmplabel != "") tmplabel = "_" + Label;
  //electron variables
  TH1F *sigmaIEtaIEta_Barrel = new TH1F((string("sigmaIEtaIEta_Barrel")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.05);
  TH1F *sigmaIPhiIPhi_Barrel = new TH1F((string("sigmaIPhiIPhi_Barrel")+tmplabel).c_str(), "; sigma iphi iphi; Fraction of Events ", 100, 0, 0.05);
  TH1F *DeltaEta_Barrel = new TH1F((string("DeltaEta_Barrel")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Barrel = new TH1F((string("DeltaPhi_Barrel")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.03, 0.03);
  TH1F *HOverE_Barrel = new TH1F((string("HOverE_Barrel")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *relIso_Barrel = new TH1F((string("relIso_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *ecalIso03_Barrel = new TH1F((string("ecalIso03_Barrel")+tmplabel).c_str(), "; ecalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *ecalIso04_Barrel = new TH1F((string("ecalIso04_Barrel")+tmplabel).c_str(), "; ecalIso04; Fraction of Events ", 100, 0, 1.0);
  TH1F *hcalIso03_Barrel = new TH1F((string("hcalIso03_Barrel")+tmplabel).c_str(), "; hcalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *hcalIso04_Barrel = new TH1F((string("hcalIso04_Barrel")+tmplabel).c_str(), "; hcalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *relIso04_Barrel = new TH1F((string("relIso04_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Barrel = new TH1F((string("fBrem_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *likelihood_0brem_Barrel = new TH1F((string("likelihood_0brem_Barrel")+tmplabel).c_str(), "; likelihood; Fraction of Events ", 100, 0, 1.0);
  TH1F *likelihood_1brem_Barrel = new TH1F((string("likelihood_1brem_Barrel")+tmplabel).c_str(), "; likelihood; Fraction of Events ", 100, 0, 1.0);

  TH1F *sigmaIEtaIEta_Endcap = new TH1F((string("sigmaIEtaIEta_Endcap")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.1);
  TH1F *sigmaIPhiIPhi_Endcap = new TH1F((string("sigmaIPhiIPhi_Endcap")+tmplabel).c_str(), "; sigma iphi iphi; Fraction of Events ", 100, 0, 0.1);  TH1F *DeltaEta_Endcap = new TH1F((string("DeltaEta_Endcap")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Endcap = new TH1F((string("DeltaPhi_Endcap")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.03, 0.03);
  TH1F *HOverE_Endcap = new TH1F((string("HOverE_Endcap")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *ecalIso03_Endcap = new TH1F((string("ecalIso03_Endcap")+tmplabel).c_str(), "; ecalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *ecalIso04_Endcap = new TH1F((string("ecalIso04_Endcap")+tmplabel).c_str(), "; ecalIso04; Fraction of Events ", 100, 0, 1.0);
  TH1F *hcalIso03_Endcap = new TH1F((string("hcalIso03_Endcap")+tmplabel).c_str(), "; hcalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *hcalIso04_Endcap = new TH1F((string("hcalIso04_Endcap")+tmplabel).c_str(), "; hcalIso03; Fraction of Events ", 100, 0, 1.0);
  TH1F *relIso_Endcap = new TH1F((string("relIso_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *relIso04_Endcap = new TH1F((string("relIso04_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Endcap = new TH1F((string("fBrem_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *likelihood_0brem_Endcap = new TH1F((string("likelihood_0brem_Endcap")+tmplabel).c_str(), "; likelihood; Fraction of Events ", 100, 0, 1.0);
  TH1F *likelihood_1brem_Endcap = new TH1F((string("likelihood_1brem_Endcap")+tmplabel).c_str(), "; likelihood; Fraction of Events ", 100, 0, 1.0);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  

  infile = new TFile(dataInputFilename.c_str());
  assert(infile);

    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kFALSE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("merged_JsonReRecoSep17_JsonStreamExpressV2.txt"); 
  

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(), "Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
	
//     mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//     if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
//     if (!passHLT(info->triggerBits, info->runNum, isMuonData)) continue;


    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear(); 
    jetArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    jetBr->GetEntry(ientry);


    //********************************************************
    // TcMet
    //********************************************************
    TVector3 met;        
    if(info->tcMEx!=0 || info->tcMEy!=0) {       
      met.SetXYZ(info->tcMEx, info->tcMEy, 0);
    }
	
    //********************************************************
    // TcMet
    //********************************************************

    Int_t NLeptons = 0;
    vector<Int_t> leptonType;
    vector<Int_t> leptonIndex;
    vector<Double_t> leptonPt;
    vector<Double_t> leptonEta;
    vector<Double_t> leptonPhi;

    Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      if ( (0==0)
//            && 
//            passElectronCuts(ele)
           &&
           fabs(ele->scEta) < 2.5
           && 
           ele->pt > 10.0
        ) {
        leptonPt.push_back(ele->pt);
        leptonEta.push_back(ele->eta);
        leptonPhi.push_back(ele->phi);
        leptonType.push_back(11);
        leptonIndex.push_back(i);
      }
    }
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
      if ( (0==0)
           &&
           passMuonCuts(mu)
           &&
           fabs(mu->eta) < 2.1
           && 
           mu->pt > 10.0
        ) {
        leptonPt.push_back(mu->pt);
        leptonEta.push_back(mu->eta);
        leptonPhi.push_back(mu->phi);
        leptonType.push_back(13);
        leptonIndex.push_back(i);  
      }
    }
    //sort leptons
    Int_t tempType;
    Int_t tempIndex;
    Double_t tempPt;
    Double_t tempEta;
    Double_t tempPhi;
    for (int l=0; l<leptonIndex.size(); l++) {
      for (int k=0; k < leptonIndex.size() - 1; k++) {
        if (leptonPt[k+1] > leptonPt[k]) {
          tempType = leptonType[k];
          tempIndex = leptonIndex[k];
          tempPt = leptonPt[k];
          tempEta = leptonEta[k];
          tempPhi = leptonPhi[k];
          
          leptonType[k] = leptonType[k+1];
          leptonIndex[k] = leptonIndex[k+1];
          leptonPt[k] = leptonPt[k+1];
          leptonEta[k] = leptonEta[k+1];
          leptonPhi[k] = leptonPhi[k+1];

          leptonType[k+1] = tempType;
          leptonIndex[k+1] = tempIndex;
          leptonPt[k+1] = tempPt;
          leptonEta[k+1] = tempEta;
          leptonPhi[k+1] = tempPhi;
          
        }
      }
    }
    
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

      if (!(jet->pt > 25 && fabs(jet->eta) < 5.0)) continue;
      if (!leadingJet || jet->pt > leadingJet->pt) {
        leadingJet = jet;
      }
      NJets++;
    }


    //******************************************************************************
    //dilepton preselection
    //******************************************************************************
    if (leptonPt.size() < 2) continue;
    if (!(leptonPt[0] > 20.0 && leptonPt[1] > 10.0)) continue;
    if (leptonPt.size() >= 3 || leptonPt[2]>10.0) continue;



    Int_t finalState = -1;
    if (leptonType[0] == 11 && leptonType[1] == 11) {
      finalState = 0;
    } else if (leptonType[0] == 13 && leptonType[1] == 13) {
      finalState = 1;
    } else if (leptonType[0] == 11 && leptonType[1] == 13) {
      finalState = 2;
    } else if (leptonType[0] == 13 && leptonType[1] == 11) {
      finalState = 3;
    }

    //******************************************************************************
    //construct event variables
    //******************************************************************************
    mithep::FourVectorM lepton1;
    mithep::FourVectorM lepton2;
    if (leptonType[0] == 11) {
      lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
    } else {
      lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
    }
    if (leptonType[1] == 11) {
      lepton2.SetCoordinates(leptonPt[1], leptonEta[1], leptonPhi[1], 0.51099892e-3 );
    } else {
      lepton2.SetCoordinates(leptonPt[1], leptonEta[1], leptonPhi[1], 105.658369e-3 );
    }
    mithep::FourVectorM dilepton = lepton1+lepton2;

    double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(leptonPhi[0], 
                                                         leptonPhi[1])* 180.0 / TMath::Pi();    
    double deltaPhiDileptonMet = mithep::MathUtils::DeltaPhi(met.Phi(), 
                                                     dilepton.Phi())*180.0 / TMath::Pi();    
    double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt() * met.Phi()*
                                 (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

    //angle between MET and closest lepton
    double deltaPhiMetLepton[2] = {mithep::MathUtils::DeltaPhi(met.Phi(), lepton1.Phi()),
                                   mithep::MathUtils::DeltaPhi(met.Phi(), lepton2.Phi())};
  
    double mTW[2] = {TMath::Sqrt(2.0*lepton1.Pt()*met.Pt()*
                                 (1.0 - cos(deltaPhiMetLepton[0]))),
                     TMath::Sqrt(2.0*lepton2.Pt()*met.Pt()*
                                 (1.0 - cos(deltaPhiMetLepton[1])))};

    double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
      deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

    double METdeltaPhilEt = met.Pt();
    if(minDeltaPhiMetLepton < TMath::Pi()/2.)
      METdeltaPhilEt = METdeltaPhilEt * sin(minDeltaPhiMetLepton);

    
    //*********************************************************************************************
    //Define Cuts
    //*********************************************************************************************
    const int nCuts = 7;
    bool passCut[nCuts] = {false, false, false, false, false, false, false};
  
    if(lepton1.Pt() >  20.0 &&
       lepton2.Pt() >= 20.0) passCut[0] = true;
    
    if(met.Pt()    > 20.0)               passCut[1] = true;
  
    if(dilepton.M() > 12.0)            passCut[2] = true;
   
    if(NJets     < 1)              passCut[5] = true;

//     if(CleanLeptons->GetEntries() == 2 &&
//        SoftMuons->GetEntries() == 0)      passCut[6] = true;

    if (finalState == 0 || finalState == 1){ // mumu/ee
      if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[3] = true;
      if(METdeltaPhilEt > 35) passCut[4] = true;
    }
    else if(finalState == 2 ||finalState == 3 ) { // emu
      passCut[3] = true;
      if(METdeltaPhilEt > 20) passCut[4] = true;
    }
 
    //*********************************************************************************************
    //Fill Histograms
    //*********************************************************************************************
    if (
//       passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] 
//         && 
      (finalState == 0 || finalState == 3 ) 
      ) {
      ptMin->Fill(lepton2.Pt());
      ptMax->Fill(lepton1.Pt());

      //fill electron variables for 2nd electron
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[1]]);

      //likelihood
      LikelihoodMeasurements measurements;
      measurements.pt = ele->pt;
      measurements.subdet = (fabs(ele->eta)<1.479) ? 0 : 1;
      measurements.deltaPhi = ele->deltaPhiIn;
      measurements.deltaEta = ele->deltaEtaIn;
      measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout;
      measurements.eSuperClusterOverP = ele->EOverP;
      measurements.hadronicOverEm = ele->HoverE;
      measurements.sigmaIEtaIEta = ele->sigiEtaiEta;
      measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->sigiPhiiPhi);
      measurements.fBrem = ele->fBrem;
      measurements.nBremClusters = ele->nBrem;
      double likelihood = LH->result(measurements);




      if ( passElectronCuts(ele) &&
           ele->pt >= 20
//            ele->pt >= 15 && ele->pt < 20
//         ele->pt >= 10 && ele->pt < 15
          && ele->isUsed == kTRUE
        ) {
        if (fabs(ele->eta) < 1.5) {
//           if ((ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1) {
            sigmaIEtaIEta_Barrel->Fill( ele->sigiEtaiEta );
            DeltaEta_Barrel->Fill( ele->deltaEtaIn );
            DeltaPhi_Barrel->Fill( ele->deltaPhiIn );
            HOverE_Barrel->Fill( ele->HoverE );
            if (ele->nBrem == 0) {
              sigmaIPhiIPhi_Barrel->Fill( TMath::Sqrt(ele->sigiPhiiPhi) );
              likelihood_0brem_Barrel->Fill( likelihood );
//               cout << "likelihood : " << likelihood << " " << measurements.pt << " " << measurements.subdet << " " << measurements.deltaPhi << " " << measurements.deltaEta << " " << measurements.eSeedClusterOverPout << " " << measurements.eSuperClusterOverP << " " << measurements.hadronicOverEm << " " << measurements.sigmaIEtaIEta << " " << measurements.sigmaIPhiIPhi << " " << measurements.fBrem << " " << measurements.nBremClusters << endl;
            } else {
              likelihood_1brem_Barrel->Fill( likelihood );
            }
//           }
//           if (ele->sigiEtaiEta < 0.012) {
            ecalIso03_Barrel->Fill( TMath::Max(ele->emIso03 - 1.0, 0.0)  / ele->pt );
            ecalIso04_Barrel->Fill( TMath::Max(ele->emIso04 - 1.0, 0.0)  / ele->pt );
            hcalIso03_Barrel->Fill( ele->hadIso03  / ele->pt );
            hcalIso04_Barrel->Fill( ele->hadIso04  / ele->pt );
            relIso_Barrel->Fill( (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt );

            relIso04_Barrel->Fill( (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt );
//             if ((ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < 0) {
//               cout << ele->trkIso03 << " " << TMath::Max(ele->emIso04 - 1.0, 0.0) << " " << ele->hadIso04 << endl;
//             }
//           }
        } else {
//           if ((ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1) {
            sigmaIEtaIEta_Endcap->Fill( ele->sigiEtaiEta );
            DeltaEta_Endcap->Fill( ele->deltaEtaIn );
            DeltaPhi_Endcap->Fill( ele->deltaPhiIn );
            HOverE_Endcap->Fill( ele->HoverE );
             if (ele->nBrem == 0) {
              sigmaIPhiIPhi_Endcap->Fill( TMath::Sqrt(ele->sigiPhiiPhi) );
              likelihood_0brem_Endcap->Fill( likelihood );
            } else {
              likelihood_1brem_Endcap->Fill( likelihood );
            }
 //          } 
 //          if (ele->sigiEtaiEta < 0.035) {
            ecalIso03_Endcap->Fill( TMath::Max(ele->emIso03 - 1.0, 0.0)  / ele->pt );
            ecalIso04_Endcap->Fill( TMath::Max(ele->emIso04 - 1.0, 0.0)  / ele->pt );
            hcalIso03_Endcap->Fill( ele->hadIso03  / ele->pt );
            hcalIso04_Endcap->Fill( ele->hadIso04  / ele->pt );
            relIso_Endcap->Fill( (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt );
            relIso04_Endcap->Fill( (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt );
 //          }
        }
      }

    }




  } //end loop over data     

  infile->Close(); delete infile;
  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwElectronSelectionVariables.root", "UPDATE");

  file->WriteTObject(ptMin,ptMin->GetName(), "WriteDelete");
  file->WriteTObject(ptMax,ptMax->GetName(), "WriteDelete");


//   file->WriteTObject(sigmaIEtaIEta_Barrel,sigmaIEtaIEta_Barrel->GetName(), "WriteDelete");
   file->WriteTObject(sigmaIPhiIPhi_Barrel,sigmaIPhiIPhi_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(DeltaEta_Barrel,DeltaEta_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(DeltaPhi_Barrel,DeltaPhi_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(HOverE_Barrel,HOverE_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(likelihood_0brem_Barrel,likelihood_0brem_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(likelihood_1brem_Barrel,likelihood_1brem_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(ecalIso03_Barrel,ecalIso03_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(ecalIso04_Barrel,ecalIso04_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(hcalIso03_Barrel,hcalIso03_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(hcalIso04_Barrel,hcalIso04_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(relIso_Barrel,relIso_Barrel->GetName(), "WriteDelete");
//   file->WriteTObject(relIso04_Barrel,relIso04_Barrel->GetName(), "WriteDelete");
  
//   file->WriteTObject(sigmaIEtaIEta_Endcap,sigmaIEtaIEta_Endcap->GetName(), "WriteDelete");
   file->WriteTObject(sigmaIPhiIPhi_Endcap,sigmaIPhiIPhi_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(DeltaEta_Endcap,DeltaEta_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(DeltaPhi_Endcap,DeltaPhi_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(HOverE_Endcap,HOverE_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(likelihood_0brem_Endcap,likelihood_0brem_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(likelihood_1brem_Endcap,likelihood_1brem_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(ecalIso03_Endcap,ecalIso03_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(ecalIso04_Endcap,ecalIso04_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(hcalIso03_Endcap,hcalIso03_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(hcalIso04_Endcap,hcalIso04_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(relIso_Endcap,relIso_Endcap->GetName(), "WriteDelete");
//   file->WriteTObject(relIso04_Endcap,relIso04_Endcap->GetName(), "WriteDelete");
  file->Close();

//   //--------------------------------------------------------------------------------------------------------------
//   // Make plots
//   //==============================================================================================================
   TCanvas *cv = new TCanvas("cv","cv", 800,600);
//    ptMin->Draw();
//    cv->SaveAs("ptMin.gif");
//   ptMax->Draw();
//   cv->SaveAs("ptMax.gif");

//   sigmaIEtaIEta_Barrel_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Barrel_signal.gif");
   sigmaIPhiIPhi_Barrel->Draw();
   cv->SaveAs("sigmaIPhiIPhi_Barrel_signal.gif");
//   DeltaEta_Barrel_signal->Draw();
//   cv->SaveAs("DeltaEta_Barrel_signal.gif");
//   DeltaPhi_Barrel_signal->Draw();
//   cv->SaveAs("DeltaPhi_Barrel_signal.gif");
//   HOverE_Barrel_signal->Draw();
//   cv->SaveAs("HOverE_Barrel_signal.gif");
//   relIso_Barrel_signal->Draw();
//   cv->SaveAs("relIso_Barrel_signal.gif");

//   sigmaIEtaIEta_Endcap_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Endcap_signal.gif");
   sigmaIPhiIPhi_Endcap->Draw();
   cv->SaveAs("sigmaIPhiIPhi_Endcap.gif");
//   DeltaEta_Endcap_signal->Draw();
//   cv->SaveAs("DeltaEta_Endcap_signal.gif");
//   DeltaPhi_Endcap_signal->Draw();
//   cv->SaveAs("DeltaPhi_Endcap_signal.gif");
//   HOverE_Endcap_signal->Draw();
//   cv->SaveAs("HOverE_Endcap_signal.gif");
//   relIso_Endcap_signal->Draw();
//   cv->SaveAs("relIso_Endcap_signal.gif");


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 


        
  gBenchmark->Show("plotZ");       
} 

void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname) {  
  TFile *file = new TFile("HwwElectronSelectionVariables.root", "READ");

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TH1F *signal = 0;
  TH1F *bkg = 0;
  TLegend *legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
    
 
  signal = (TH1F*)file->Get(signalHistName.c_str());
  bkg = (TH1F*)file->Get(bkgHistName.c_str());
  if (!signal) {cout << "Cannot load histogram " << signalHistName << endl;  assert(signal);}
  if (!bkg) {cout << "Cannot load histogram " << bkgHistName << endl;  assert(bkg);}

  Double_t norm = 0;
  for (int b=0; b<signal->GetXaxis()->GetNbins()+2; ++b) { norm += signal->GetBinContent(b); }
  for (int b=0; b<signal->GetXaxis()->GetNbins()+2; ++b) {
    signal->SetBinContent(b,signal->GetBinContent(b) / norm);
    signal->SetBinError(b,signal->GetBinError(b) / norm);
  }
  norm = 0;
  for (int b=0; b<bkg->GetXaxis()->GetNbins()+2; ++b) { norm += bkg->GetBinContent(b); }
  for (int b=0; b<bkg->GetXaxis()->GetNbins()+2; ++b) {
    bkg->SetBinContent(b,bkg->GetBinContent(b) / norm);
    bkg->SetBinError(b,bkg->GetBinError(b) / norm);
  }

  //Create Efficiency vs Bkg Rejection plots
  const int nPoints = signal->GetXaxis()->GetNbins();
  double sigEff[nPoints];
  double bkgEff[nPoints];
  double x[nPoints];
  double_t minX = signal->GetXaxis()->GetBinLowEdge(1);
  double_t maxX = signal->GetXaxis()->GetBinUpEdge(signal->GetXaxis()->GetNbins());
  double_t binsize = signal->GetXaxis()->GetBinWidth(1);
  for (int i=0; i<nPoints; ++i) {
    x[i] = minX + binsize*(i+1);
    sigEff[i] = 0;
    bkgEff[i] = 0;
    for (int k=0; k<=i; ++k) {
      sigEff[i] += signal->GetBinContent(k);
      bkgEff[i] += bkg->GetBinContent(k);
//       cout << signal->GetBinContent(k) << " ";
    }
//     cout << " -> " << sigEff[i] << endl;
//     sigEff[i] = signal->Integral(0,i);
//     bkgEff[i] = bkg->Integral(0,i);
  }
  TGraph *SignalEffVsBkgEff = new TGraph (nPoints, bkgEff, sigEff);
  SignalEffVsBkgEff->SetMarkerColor(kRed);
  SignalEffVsBkgEff->GetXaxis()->SetTitleOffset(1.02);
  SignalEffVsBkgEff->GetYaxis()->SetTitleOffset(1.05);
  TGraph *SignalEff = new TGraph (nPoints, x, sigEff);
  SignalEff->SetMarkerColor(kRed);
  SignalEff->GetXaxis()->SetTitleOffset(1.02);
  SignalEff->GetYaxis()->SetTitleOffset(1.05);
  SignalEff->GetXaxis()->SetTitle((string(signal->GetXaxis()->GetTitle())+" Cut").c_str());
  SignalEff->GetYaxis()->SetTitle("Efficiency");
  TGraph *BkgEff = new TGraph (nPoints, x, bkgEff);
  BkgEff->SetMarkerColor(kRed);
  BkgEff->GetXaxis()->SetTitleOffset(1.02);
  BkgEff->GetYaxis()->SetTitleOffset(1.05);
  BkgEff->GetXaxis()->SetTitle((string(signal->GetXaxis()->GetTitle())+" Cut").c_str());
  BkgEff->GetYaxis()->SetTitle("Efficiency");

//   signalEff->SetLineStyle(1);
//   signalEff->SetMinimum(XMin);
//   signalEff->SetMaximum(XMax);
//   signalEff->GetXaxis()->SetTitleSize(0.04);
//   signalEff->GetXaxis()->SetLabelSize(0.04);
//   signalEff->GetYaxis()->SetTitleSize(0.04);
//   signalEff->GetYaxis()->SetLabelSize(0.04);
//   signalEff->GetYaxis()->CenterTitle(kTRUE);
  BkgEff->SetMarkerColor(kBlue);
  BkgEff->GetXaxis()->SetTitleOffset(1.02);
  BkgEff->GetYaxis()->SetTitleOffset(1.05);



  //Plot distributions
  Double_t maxY = signal->GetMaximum();
  if (bkg->GetMaximum() > maxY) maxY = bkg->GetMaximum();
  maxY = maxY * 1.2;
  signal->SetMaximum(maxY);
  bkg->SetMaximum(maxY);

  signal->SetLineColor(kRed);
  signal->SetMarkerColor(kRed);
  bkg->SetLineColor(kBlue);
  bkg->SetMarkerColor(kBlue);
  signal->GetXaxis()->SetTitleOffset(1.05);  
  signal->GetYaxis()->SetTitleOffset(1.1);
  bkg->GetXaxis()->SetTitleOffset(1.05);  
  bkg->GetYaxis()->SetTitleOffset(1.1);
  legend->Clear();
  legend->AddEntry(signal, signalLabel.c_str(),  "LP");
  legend->AddEntry(bkg,    bkgLabel.c_str(), "LP");

  signal->Draw();
  bkg->Draw("same");
  legend->Draw();
  cv->SaveAs((plotname+".gif").c_str());

  //Plot Efficiencies
  legend = new TLegend(0.73,0.25,0.93,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->Clear();
  legend->AddEntry(SignalEff, signalLabel.c_str(),  "P");
  legend->AddEntry(BkgEff, bkgLabel.c_str(),  "P");
  SignalEff->Draw("AP");
  BkgEff->Draw("Psame");
  legend->Draw();
  cv->SaveAs((plotname+"_Efficiency.gif").c_str());


  //Plot SignalEff Vs BkgEff
  SignalEffVsBkgEff->Draw("AP");
  SignalEffVsBkgEff->GetXaxis()->SetTitle((bkgLabel+" Efficiency").c_str());
  SignalEffVsBkgEff->GetYaxis()->SetTitle((signalLabel+" Efficiency").c_str());
  cv->SaveAs((plotname+"_SigEffVsBkgEff.gif").c_str());

  file->Close();
}


Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;

//   if (isMuonData) {
//     if ((runNum >= 132440) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
//     } 
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140401)) {
// //       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   } else {
//     //it's electron data
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140041)) {
// //       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( !(triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   }

  return pass;

}



Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.025
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
          )
      ) {
      pass = kFALSE;
    }
  } else {
    pass = kFALSE;
    return pass;
  }

  return pass;
}


Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
    return pass;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.014          
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.5
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.034
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.5
          )
      ) {
      pass = kFALSE;
    }
  }

  return pass;
}

Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

        && (mu->nSeg >= 2)
        && (mu->nValidHits >= 1)
        && (mu->nPixHits >= 1)
        )
    ) pass = kFALSE;

  return pass;
}

//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2)
{
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << mass;
  
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << endl;
  
}
