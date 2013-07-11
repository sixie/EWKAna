//root -l -b -q EWKAna/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Utils/LeptonIDCuts.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

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


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passElectronNumeratorCuts(const mithep::TElectron *ele);
Bool_t passElectronLHNumeratorCuts(const mithep::TElectron *ele, Double_t likelihoodValue, Int_t Option);
Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele, Int_t DenominatorType, string SampleType);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename,  
                               const string smurfOutputFilename, Int_t Option = 0);

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
void ComputeElectronFakeRate_Data() {


//************************************
//Tests With Different WPs
//************************************

//   DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV3.root", 0);
//  DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4L.root", 1);
//   DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4M.root", 2);
//   DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4T.root", 3);
//   DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr-v2_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4.root", 3);

//    DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV3.root", 0);
//  DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4L.root", 1);
//    DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4M.root", 2);
//    DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4T.root", 3);
//   DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4TTruncated.root", 4);
//    DoComputeElectronFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4TTruncated20.root", 5);

//     DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4Cone03L.root", 11);
//      DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4Cone03M.root", 12);
//        DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4Cone03T.root", 13);

//     DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4ConeHybridL.root", 101);
//     DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4ConeHybridM.root", 102);
//     DoComputeElectronFakeRate("WWAnalysisSkimmed_r11a-del-pr_EleFakeRateTriggerSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV4ConeHybridT.root", 103);




//************************************
//Real Stuff
//************************************
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfV6.skim.root", "FakeRates_Electron_V6.root", 0);

//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV6.m10.skim.root", "FakeRates_SmurfV6.m10.root", 0);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV6.prv4.skim.root", "FakeRates_SmurfV6.prv4.root", 0);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV6.a05.skim.root", "FakeRates_SmurfV6.a05.root", 0);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v6_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV6.prv6.skim.root", "FakeRates_SmurfV6.prv6.root", 0);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfV6.r11b.skim.root", "FakeRates_SmurfV6.r11b.root", 0);




//************************************
//Likelihood Stuff
//************************************
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfLH.skim.root", "FakeRates_Electron_LH.root", 10);

//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfLHElectron.m10.skim.root", "FakeRates_SmurfLHElectron.m10.root", 10);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfLHElectron.prv4.skim.root", "FakeRates_SmurfLHElectron.prv4.root", 10);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfLHElectron.a05.skim.root", "FakeRates_SmurfLHElectron.a05.root", 10);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v6_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfLHElectron.prv6.skim.root", "FakeRates_SmurfLHElectron.prv6.root", 10);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfLHElectron.r11b.skim.root", "FakeRates_SmurfLHElectron.r11b.root", 10);


//************************************
//MVA 
//************************************

//------------------
//Fully Inclusive
//------------------
  DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVAIDIsoCombined.skim.root", "FakeRates_Electron_BDTGIDIsoCombined", 25);
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVANoIPInfo.skim.root", "FakeRates_SmurfBDTGNoIPInfoElectron.root", 21);
  DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVAWithIPInfo.skim.root", "FakeRates_Electron_BDTGWithIPInfo", 20);

//------------------
//In Bins of NVtx
//------------------
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVAWithIPInfo.NVtx0To2.root", "FakeRates_Electron_BDTGWithIPInfo_NVtx0To2", 21);
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVAWithIPInfo.NVtx3To5.root", "FakeRates_Electron_BDTGWithIPInfo_NVtx3To5", 22);
//   DoComputeElectronFakeRate("LIST","ElectronFakeRate","ElectronFakeRate.SmurfMVAWithIPInfo.NVtx6To10.root", "FakeRates_Electron_BDTGWithIPInfo_NVtx6To10", 23);
  

//------------------
//2011A Vs 2011B
//------------------
//   DoComputeElectronFakeRate("RUN2011A","ElectronFakeRate","ElectronFakeRate.SmurfMVAWithIPInfo.skim.Run2011A.root", "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.Run2011A", 20);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.r11b.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.NoTriggerObjMatch.r11b", 20);

//------------------
//Different Run Eras
//------------------
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.m10.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.m10.root", 20);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.prv4.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.prv4.root", 20);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.a05.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.a05.root", 20);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v6_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.prv6.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.prv6.root", 20);
//   DoComputeElectronFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root","ElectronFakeRate","ElectronFakeRate.SmurfBDTGWithIPInfoElectron.r11b.skim.root", "FakeRates_SmurfBDTGWithIPInfoElectron.r11b.root", 20);




}



void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               const string smurfOutputFilename,
                               Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 26.5;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *fileLH = TFile::Open("MitPhysics/data/ElectronLikelihoodPdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");

  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useDeltaEta = true;
  defaultSwitches.m_useDeltaPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useSigmaEtaEta = true;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useOneOverEMinusOneOverP = true;
 
  ElectronLikelihood *LH = new ElectronLikelihood(&(*EB0lt15dir),&(*EB1lt15dir), &(*EElt15dir), 
                                                  &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir),
                                                  defaultSwitches,
                                                  std::string("class"),std::string("class"),true,true);

  mithep::ElectronIDMVA *electronIDMVA = 0;

  if (Option >= 20 && Option <= 29) {

    if (Option == 21) {
      electronIDMVA = new mithep::ElectronIDMVA();
      electronIDMVA->Initialize("BDTG method",
                                        "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", 
                                        "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", 
                                        "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", 
                                        "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", 
                                        "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", 
                                        "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml",
                                        mithep::ElectronIDMVA::kNoIPInfo );
    }
    
    if (Option == 20) {
      electronIDMVA = new mithep::ElectronIDMVA();
      electronIDMVA->Initialize("BDTG method",
                                          "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",
                                          mithep::ElectronIDMVA::kWithIPInfo);
    } 

    if (Option == 25) {
      electronIDMVA = new mithep::ElectronIDMVA();
      electronIDMVA->Initialize("BDTG method",
                                          "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml", 
                                          "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml",
                                          mithep::ElectronIDMVA::kIDIsoCombined);

    }

  }


  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
//   ptbins.push_back(10);  
//   ptbins.push_back(80);  


  ptbins.push_back(10);  
  ptbins.push_back(12.5);  
  ptbins.push_back(15);  
  ptbins.push_back(17.5);  
  ptbins.push_back(20);  
  ptbins.push_back(22.5);  
  ptbins.push_back(25);  
  ptbins.push_back(27.5);  
  ptbins.push_back(30);  
  ptbins.push_back(32.5);  
  ptbins.push_back(35);  
  ptbins.push_back(37.5);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
  ptbins.push_back(60);  
  ptbins.push_back(70);  
  ptbins.push_back(80);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(0.0);
  etabins.push_back(0.25);
  etabins.push_back(0.5);
  etabins.push_back(0.75);
  etabins.push_back(1.0);
  etabins.push_back(1.25);
  etabins.push_back(1.479);
  etabins.push_back(1.75);
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.5);
  etabins.push_back(3.0);


  vector<double> phibins;
  phibins.push_back(-3.25);
  phibins.push_back(-2.75);
  phibins.push_back(-2.25);
  phibins.push_back(-1.75);
  phibins.push_back(-1.25);
  phibins.push_back(-0.75);
  phibins.push_back(-0.25);
  phibins.push_back(0.25);
  phibins.push_back(0.75);
  phibins.push_back(1.25);
  phibins.push_back(1.75);
  phibins.push_back(2.25);
  phibins.push_back(2.75);
  phibins.push_back(3.25);

  vector<double> nvtxbins;
  nvtxbins.push_back(0);
  nvtxbins.push_back(2);
  nvtxbins.push_back(4);
  nvtxbins.push_back(6);
  nvtxbins.push_back(8);
  nvtxbins.push_back(10);
  nvtxbins.push_back(12);
  nvtxbins.push_back(14);
  nvtxbins.push_back(16);
  nvtxbins.push_back(18);
  nvtxbins.push_back(20);
  nvtxbins.push_back(22);
  nvtxbins.push_back(24);
  nvtxbins.push_back(26);



  vector<double> rhobins;
  rhobins.push_back(0);
  rhobins.push_back(2);
  rhobins.push_back(4);
  rhobins.push_back(6);
  rhobins.push_back(8);
  rhobins.push_back(10);
  rhobins.push_back(12);
  rhobins.push_back(14);
  rhobins.push_back(16);
  rhobins.push_back(18);
  rhobins.push_back(20);
  rhobins.push_back(22);
  rhobins.push_back(24);
  rhobins.push_back(26);



  vector<double> ptbins2D;
  ptbins2D.push_back(10);  
  ptbins2D.push_back(15);  
  ptbins2D.push_back(20);  
  ptbins2D.push_back(25);  
  ptbins2D.push_back(30);  
  ptbins2D.push_back(35);  
  vector<double> etabins2D;
  etabins2D.push_back(0.0);
  etabins2D.push_back(1.0);
  etabins2D.push_back(1.479);
  etabins2D.push_back(2.0);
  etabins2D.push_back(2.5);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histLeadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH1F*> > > DenominatorVector_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Rho;
  vector<vector<vector<TH2F*> > > DenominatorVector_PtEta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH1F*> > > NumeratorVector_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Rho;
  vector<vector<vector<TH2F*> > > NumeratorVector_PtEta;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_Rho;


  vector<Double_t> denominatorType;
  denominatorType.push_back(4);
  vector<string> sampleLabel;
  sampleLabel.push_back("Ele8Sample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL");
  sampleLabel.push_back("Ele17CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLJet40Sample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLPtCombinedSample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLCombinedSample");
//    sampleLabel.push_back("PhotonJetsSample");
  vector<Double_t> ptThreshold;
//   ptThreshold.push_back(0);
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  ptThreshold.push_back(35);
  ptThreshold.push_back(40);
  ptThreshold.push_back(45);
  ptThreshold.push_back(50);
//   ptThreshold.push_back(70);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH1F*> > tmpDenominatorVector_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Rho;
    vector<vector<TH2F*> > tmpDenominatorVector_PtEta;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH1F*> > tmpNumeratorVector_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Rho;
    vector<vector<TH2F*> > tmpNumeratorVector_PtEta;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_Rho;


    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH1F*>  tmptmpDenominatorVector_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Rho;
      vector<TH2F*>  tmptmpDenominatorVector_PtEta;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH1F*>  tmptmpNumeratorVector_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Rho;
      vector<TH2F*>  tmptmpNumeratorVector_PtEta;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho;
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);

        TH1F *histDenominator_Pt10To20_Barrel_NVtx = new TH1F(("histDenominator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Barrel_Rho = new TH1F(("histDenominator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt10To20_Endcap_NVtx = new TH1F(("histDenominator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Endcap_Rho = new TH1F(("histDenominator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Barrel_NVtx = new TH1F(("histDenominator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Barrel_Rho = new TH1F(("histDenominator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Endcap_NVtx = new TH1F(("histDenominator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Endcap_Rho = new TH1F(("histDenominator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Barrel_NVtx = new TH1F(("histNumerator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Barrel_Rho = new TH1F(("histNumerator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Endcap_NVtx = new TH1F(("histNumerator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Endcap_Rho = new TH1F(("histNumerator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Barrel_NVtx = new TH1F(("histNumerator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Barrel_Rho = new TH1F(("histNumerator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Endcap_NVtx = new TH1F(("histNumerator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Endcap_Rho = new TH1F(("histNumerator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);



        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_NVtx.push_back(histDenominator_NVtx);
        tmptmpDenominatorVector_Rho.push_back(histDenominator_Rho);
        tmptmpDenominatorVector_PtEta.push_back(histDenominator_PtEta);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_NVtx.push_back(histNumerator_NVtx);
        tmptmpNumeratorVector_Rho.push_back(histNumerator_Rho);
        tmptmpNumeratorVector_PtEta.push_back(histNumerator_PtEta);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);
        tmptmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(histDenominator_Pt10To20_Barrel_NVtx);
        tmptmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(histDenominator_Pt10To20_Barrel_Rho);
        tmptmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(histDenominator_Pt10To20_Endcap_NVtx);
        tmptmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(histDenominator_Pt10To20_Endcap_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(histDenominator_Pt20ToInf_Barrel_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(histDenominator_Pt20ToInf_Barrel_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(histDenominator_Pt20ToInf_Endcap_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(histDenominator_Pt20ToInf_Endcap_Rho);
        tmptmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(histNumerator_Pt10To20_Barrel_NVtx);
        tmptmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(histNumerator_Pt10To20_Barrel_Rho);
        tmptmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(histNumerator_Pt10To20_Endcap_NVtx);
        tmptmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(histNumerator_Pt10To20_Endcap_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(histNumerator_Pt20ToInf_Barrel_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(histNumerator_Pt20ToInf_Barrel_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(histNumerator_Pt20ToInf_Endcap_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(histNumerator_Pt20ToInf_Endcap_Rho);

      }
      tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
      tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
      tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
      tmpDenominatorVector_NVtx.push_back(tmptmpDenominatorVector_NVtx);
      tmpDenominatorVector_Rho.push_back(tmptmpDenominatorVector_Rho);
      tmpDenominatorVector_PtEta.push_back(tmptmpDenominatorVector_PtEta);
      tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
      tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
      tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
      tmpNumeratorVector_NVtx.push_back(tmptmpNumeratorVector_NVtx);
      tmpNumeratorVector_Rho.push_back(tmptmpNumeratorVector_Rho);
      tmpNumeratorVector_PtEta.push_back(tmptmpNumeratorVector_PtEta);
      tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
      tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
      tmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_NVtx);
      tmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_Rho);
      tmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_NVtx);
      tmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_Rho);
      tmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
      tmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho);
      tmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
      tmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho);
      tmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_NVtx);
      tmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_Rho);
      tmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_NVtx);
      tmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_Rho);
      tmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
      tmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho);
      tmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
      tmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_NVtx.push_back(tmpDenominatorVector_NVtx);
    DenominatorVector_Rho.push_back(tmpDenominatorVector_Rho);
    DenominatorVector_PtEta.push_back(tmpDenominatorVector_PtEta);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_NVtx.push_back(tmpNumeratorVector_NVtx);
    NumeratorVector_Rho.push_back(tmpNumeratorVector_Rho);
    NumeratorVector_PtEta.push_back(tmpNumeratorVector_PtEta);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
    DenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmpDenominatorVector_Pt10To20_Barrel_NVtx);
    DenominatorVector_Pt10To20_Barrel_Rho.push_back(tmpDenominatorVector_Pt10To20_Barrel_Rho);
    DenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmpDenominatorVector_Pt10To20_Endcap_NVtx);
    DenominatorVector_Pt10To20_Endcap_Rho.push_back(tmpDenominatorVector_Pt10To20_Endcap_Rho);
    DenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
    DenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_Rho);
    DenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
    DenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_Rho);
    NumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmpNumeratorVector_Pt10To20_Barrel_NVtx);
    NumeratorVector_Pt10To20_Barrel_Rho.push_back(tmpNumeratorVector_Pt10To20_Barrel_Rho);
    NumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmpNumeratorVector_Pt10To20_Endcap_NVtx);
    NumeratorVector_Pt10To20_Endcap_Rho.push_back(tmpNumeratorVector_Pt10To20_Endcap_Rho);
    NumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
    NumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_Rho);
    NumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
    NumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_Rho);
  }

  ofstream eventListFile("eventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 

  Int_t NEvents = 0;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-o03-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root");
  } else if (inputFilename == "RUN2011A") {
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-o03-v1_EleFakeRateTriggerAndDenominatorSkim.root");
  }
  else {
    inputfiles.push_back(inputFilename);
  }

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
//       cout << "start event " << ientry << " : " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " \n";
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      
      //For MVA, only use odd event numbers because even event numbers were used for training
      if (Option >= 20 && Option < 30) {
        if (info->evtNum % 2 == 0) continue;

//         if (Option == 21) if (!(info->nPV0 >= 0 && info->nPV0 <= 2)) continue;
//         if (Option == 22) if (!(info->nPV0 >= 3 && info->nPV0 <= 5)) continue;
//         if (Option == 23) if (!(info->nPV0 >= 6 && info->nPV0 <= 10)) continue;

      }

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      Double_t rho = 0;
      if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);


      //********************************************************
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      Int_t NElectrons = electronArr->GetEntries();
      
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

        //likelihood value
        mithep::FourVectorM tmpEleSC;
        tmpEleSC.SetCoordinates(ele->scEt, ele->scEta, ele->scPhi, 0.51099892e-3 );
      
        LikelihoodMeasurements measurements;
        measurements.pt = ele->pt;
        if (ele->isEB && (fabs(ele->eta)<1.0)) measurements.subdet = 0;
        else if (ele->isEB) measurements.subdet = 1;
        else measurements.subdet = 2;
        measurements.deltaPhi = ele->deltaPhiIn;
        measurements.deltaEta = ele->deltaEtaIn;
        measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout;
        measurements.eSuperClusterOverP = ele->EOverP;
        measurements.hadronicOverEm = ele->HoverE;
        measurements.sigmaIEtaIEta = ele->sigiEtaiEta;
        measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->sigiPhiiPhi);
        measurements.fBrem = ele->fBrem;
        measurements.nBremClusters = ele->nBrem;
        measurements.OneOverEMinusOneOverP = (1.0/(ele->EOverP*ele->p)) - 1.0 / ele->p;
        double likelihood = LH->result(measurements);
        Double_t likelihoodValue = 0; 
        if (likelihood <= 0) {
          likelihoodValue = -20.0;
        } else if (likelihood == 1) {
          likelihoodValue = 20.0;
        } else {
          likelihoodValue = TMath::Log(likelihood / (1.0-likelihood));
        }




        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = ele->pt;

//         cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
//         for(Int_t e=0; e<electronArr->GetEntries(); e++) {
//           const mithep::TElectron *el = (mithep::TElectron*)((*electronArr)[e]);
//           cout << "Ele " << e << " : " << el->pt << " " << el->eta << " " << el->phi << " " 
//                << el->ChargedIso04+el->NeutralHadronIso04_10Threshold
//             +el->GammaIso04_10Threshold-el->GammaIsoVetoEtaStrip04_10Threshold
//             -el->ChargedEMIsoVetoEtaStrip04-el->NeutralHadronIso007_10Threshold << " "
//                << endl;
//         }
//         for(Int_t j=0; j<jetArr->GetEntries(); j++) {
//           const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
//           cout << "Jet " << j << " : " << jet->pt << " " << jet->eta << " " << jet->phi << endl;
//         }
        }


        //********************************************************
        // Photons
        //********************************************************
        Int_t NPhotons = 0;      
        Double_t photonEt = -1;
        for(Int_t p=0; p<photonArr->GetEntries(); p++) {
          const mithep::TPhoton *photon = (mithep::TPhoton*)((*photonArr)[p]);
        
          Bool_t isEle = kFALSE;
          for(Int_t e=0; e<electronArr->GetEntries(); e++) {
            const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[e]);   
            if (mithep::MathUtils::DeltaR(photon->phi, photon->eta, tmpEle->phi, tmpEle->eta) < 0.3) isEle = kTRUE;
          }
 


          //photon ID
          if ( photon->et > 20
               && !isEle
               && !photon->hasPixelSeed
               && photon->emIso03 < 2.0 + 0.006*photon->et
               && photon->hadIso03 < 2.0+0.0025*photon->et
               && photon->trkIso03Hollow < 1.5 + 0.001*photon->et
               && ( (fabs(photon->eta) < 1.5 && photon->sigiEtaiEta < 0.01) || (fabs(photon->eta) >= 1.5 && photon->sigiEtaiEta < 0.028) )
            ) {
            continue;
          }


          mithep::FourVectorM phFourVector;
          mithep::FourVectorM eleFourVector;
          eleFourVector.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
          phFourVector.SetCoordinates(photon->et, photon->eta, photon->phi, 0.0 );
          mithep::FourVectorM dilepton = phFourVector+eleFourVector;
        
//         cout << "photon " << p << " : " << photon->et << " " << photon->eta << " " << photon->phi 
//              << "Ele : " << ele->pt << " " << ele->eta << " " << ele->phi << " : " 
//              << dilepton.M() << endl;
        
          if ( fabs(dilepton.M() - 91) > 20 ) {
            photonEt = photon->et;
            NPhotons++;
          }
        }
      

        for ( UInt_t denominatorTypeIndex = 0 ; denominatorTypeIndex < denominatorType.size() ; ++denominatorTypeIndex ) {
          for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
            for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
          
              //********************************************************
              // Event Selection Cuts
              //********************************************************

              if (NElectrons > 1) continue;

              if (sampleLabel[sampleTypeIndex] == "Ele8Sample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL" || sampleLabel[sampleTypeIndex] == "Ele17CaloIdLCaloIsoVLSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLJet40Sample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLCombinedSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLPtCombinedSample") {
                if (met > 20) continue;
                Bool_t passJetSelection = kFALSE;
                if (ptThreshold[ptThresholdIndex] == 0) passJetSelection = kTRUE;
                if (leadingJetPt > ptThreshold[ptThresholdIndex]) {
                  passJetSelection = kTRUE;               
                }
//               cout << leadingJetPt << " " << ptThreshold[ptThresholdIndex] << " " << passJetSelection << endl;
                if (!passJetSelection) continue;
              }
      
              if (sampleLabel[sampleTypeIndex] == "PhotonJetsSample") {
                if (NPhotons != 1) continue;
                if (!(photonEt > ptThreshold[ptThresholdIndex])) continue;  
              }


              if (!(ele->pt > 10 && ele->pt < 35.0)) continue;

              if (passElectronDenominatorCuts(info->triggerBits, ele, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex])) {
              

//               eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << ele->pt << " " << ele->eta << " " << ele->phi << " : " << passNewElectronNumeratorCuts(ele, Option) << endl;
              
//              if (!(ele->pt > 10 && ele->pt < 20)) continue;

                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt,eventweight);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta),eventweight);
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi,eventweight);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta), eventweight);
              
                if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
              

                Double_t iso = ( ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold ) / ele->pt;
                LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(leptonJetPt, 0.01),99.9),eventweight);
              
                DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(iso, 0.000001),0.9999999),eventweight); 
              
                if (
                  (Option == 0 && passElectronNumeratorCuts(ele))
                  ||
                  (Option >= 10 && Option < 20 && passElectronLHNumeratorCuts(ele, likelihoodValue, Option))                
                  ||
                  (Option == 20 && passElectronMVA(ele, 
                                                   electronIDMVA->MVAValue(
                                                     ele->pt,ele->scEta,
                                                     ele->sigiEtaiEta, 
                                                     ele->deltaEtaIn,
                                                     ele->deltaPhiIn, 
                                                     ele->HoverE,
                                                     ele->d0,
                                                     ele->dz, 
                                                     ele->fBrem,
                                                     ele->EOverP,
                                                     ele->ESeedClusterOverPout,
                                                     TMath::Sqrt(ele->sigiPhiiPhi),
                                                     ele->nBrem,
                                                     (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p), 
                                                     ele->ESeedClusterOverPIn,
                                                     ele->ip3d,
                                                     ele->ip3dSig ), 
                                                   2))
                  ||
                  (Option == 25 && passElectronMVAIDIsoCombined(ele, 
                                                                electronIDMVA->MVAValue(
                                                                  ele->pt,ele->scEta,rho,
                                                                  ele->sigiEtaiEta, 
                                                                  ele->deltaEtaIn,
                                                                  ele->deltaPhiIn, 
                                                                  ele->HoverE,
                                                                  ele->d0,
                                                                  ele->dz, 
                                                                  ele->fBrem,
                                                                  ele->EOverP,
                                                                  ele->ESeedClusterOverPout,
                                                                  TMath::Sqrt(ele->sigiPhiiPhi),
                                                                  ele->nBrem,
                                                                  (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p), 
                                                                  ele->ESeedClusterOverPIn,
                                                                  ele->ip3d,
                                                                  ele->ip3dSig,
                                                                  ele->GsfTrackChi2OverNdof,
                                                                  ele->dEtaCalo,
                                                                  ele->dPhiCalo,
                                                                  ele->R9,
                                                                  ele->SCEtaWidth,
                                                                  ele->SCPhiWidth,
                                                                  ele->CovIEtaIPhi,
                                                                  ele->PreShowerOverRaw,
                                                                  ele->ChargedIso03,
                                                                  (ele->NeutralHadronIso03_05Threshold - ele->NeutralHadronIso007_05Threshold),
                                                                  (ele->GammaIso03_05Threshold - ele->GammaIsoVetoEtaStrip03_05Threshold),
                                                                  ele->ChargedIso04 ,
                                                                  (ele->NeutralHadronIso04_05Threshold - ele->NeutralHadronIso007_05Threshold),
                                                                  (ele->GammaIso04_05Threshold - ele->GammaIsoVetoEtaStrip04_05Threshold) 
                                                                  ), 
                                                                0))
                  
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, eventweight);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta), eventweight);
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi, eventweight);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta) , eventweight);   
        
                  if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                  if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                }
              }

            } //loop over denominator types
          } //loop over sample types
        } //loop over ptThresholds

      } //loop over electrons

//      cout << "dddata " << " : " << info->runNum << " " << info->lumiSec << " " << info->evtNum << " ";

    } //end loop over data    
    cout << "done " << inputfiles[f] << endl;

  } //end loop over files

  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        Bool_t printDebug = kFALSE;
        if (denominatorType[denominatorTypeIndex] == 4 
            // && sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" && ptThreshold[ptThresholdIndex] == 30
          ) printDebug = kTRUE;
        cout << label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta" << endl;
        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_PtEta = 
          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", ptbins2D, etabins2D, ErrorType, printDebug);
        

        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);



        TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_nvtx, efficiency_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        
        file->WriteTObject(efficiency_Pt10To20_Barrel_nvtx, efficiency_Pt10To20_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Barrel_rho, efficiency_Pt10To20_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_nvtx, efficiency_Pt10To20_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_rho, efficiency_Pt10To20_Endcap_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_nvtx, efficiency_Pt20ToInf_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_rho, efficiency_Pt20ToInf_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_nvtx, efficiency_Pt20ToInf_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_rho, efficiency_Pt20ToInf_Endcap_rho->GetName(), "WriteDelete");

        file->Close();
        
        //*****************************************************
        // Combined Fake Rate
        //*****************************************************
        if (denominatorType[denominatorTypeIndex] == 4 && 
            (sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLCombinedSample"
              ) 
          ) {
          
          TH2F *eff = 0;
          TH2F *effErrorLow = 0;
          TH2F *effErrorHigh = 0;
          
          file = new TFile((smurfOutputFilename + ".root").c_str(), "UPDATE");

          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          "ElectronFakeRate_V"+IntToString(denominatorType[denominatorTypeIndex])+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D, ErrorType, file);
        
          file->Close();
          delete file;
        }
        


//         //*****************************************************
//         // Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL Fake Rate
//         //*****************************************************
//         if (denominatorType[denominatorTypeIndex] == 4 && 
//             (sampleLabel[sampleTypeIndex] == "Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL"
//               ) 
//           ) {
          
//           TH2F *eff = 0;
//           TH2F *effErrorLow = 0;
//           TH2F *effErrorHigh = 0;
          
//           file = new TFile((smurfOutputFilename + ".CaloIdTTrkIdVLCaloIsoVLTrkIsoVL.root").c_str(), "UPDATE");

//           mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
//                                                           DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
//                                                           "ElectronFakeRate_V"+IntToString(denominatorType[denominatorTypeIndex])+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
//                                                           ptbins2D, etabins2D, ErrorType, file);
        
//           file->Close();
//           delete file;
//         }
        



      }
    }
  }


  cout << "Total Events: " << NEvents << endl;

  gBenchmark->Show("WWTemplate");       
} 

Bool_t passElectronNumeratorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (fabs(ele->eta) >= 2.5) pass = kFALSE;

  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold;
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
    +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_10Threshold;    
  Double_t iso = iso04;
  Double_t isoCutValue = 0;

  if (fabs(ele->scEta) < 1.479) {
    if (ele->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.13;
    }
  } else {
    if (ele->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.09;
    }
  }

  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.10
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }
  } 

  if (ele->pt < 20) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.03
              && ele->HoverE < 0.025
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else  {
      if (! (  (0==0)
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
            )
        ) {
        pass = kFALSE;
      }
    } 

    if (ele->fBrem <= 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 )) pass = kFALSE;
      }
    }
  }




  return pass;
}




Bool_t passElectronLHNumeratorCuts(const mithep::TElectron *ele, Double_t likelihoodValue, Int_t Option) {

  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
      +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
      -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_10Threshold;
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold;
  Double_t iso = iso04;
  
  Double_t isoCutValue = 0;


  if (fabs(ele->scEta) < 1.479) {
    if (ele->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.13;
    }
  } else {
    if (ele->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.09;
    }
  }
  
  Double_t LHCutValue = 0;

  //NEW Tight WP
  if (Option == 10) {

    if(ele->pt > 20){
      if(fabs(ele->scEta) < 1.479){
        if(ele->nBrem == 0)           LHCutValue = 3.5;
	else                          LHCutValue = 4.0;
      }
      else  {                                
        if(ele->nBrem == 0)           LHCutValue = 4.0;
	else                          LHCutValue = 4.0;
      }
    }
    else {
      if(fabs(ele->scEta) < 1.479){
        if(ele->nBrem == 0)           LHCutValue =  4.0;
	else                          LHCutValue =  4.5;
      }
      else  {                                
        if(ele->nBrem == 0)           LHCutValue =  4.0;
	else                          LHCutValue =  4.0;
      }
    }
  }


  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && likelihoodValue > LHCutValue

            && ( ele->trkIso03 ) / ele->pt < 0.2
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
             && ele->HoverE < 0.10
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && likelihoodValue > LHCutValue

             && (ele->trkIso03 ) / ele->pt < 0.2
             && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}



Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele, Int_t DenominatorType, string SampleType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  //match to HLT
  if (SampleType == "Ele8Sample") {
    if (!(triggerBits & kHLT_Ele8 
          && ele->hltMatchBits & kHLTObject_Ele8)) pass = kFALSE;
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLSample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL 
            && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
             )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
             )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele17CaloIdLCaloIsoVLSample") {
    if (!( 
          (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL
            )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLJet40Sample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
          && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
            )          
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLPtCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL)
          ||
          ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL)
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "PhotonJetsSample") {
    if (!(
          (triggerBits & kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL)
          && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL)
      ) pass = kFALSE;
  }

  if (DenominatorType == 1) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->dz) < 0.1
//               && (ele->trkIso03 ) / ele->pt < 0.2
//               && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
//               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->dz) < 0.1
//                && (ele->trkIso03 ) / ele->pt < 0.2
//                && (ele->emIso03 ) / ele->pt < 0.2
//                && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 
  }

  if (DenominatorType == 2) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->dz) < 0.1
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->dz) < 0.1
            )
        ) {
        pass = kFALSE;
      }
    } 
  }


  if (DenominatorType == 3) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->d0) < 0.02
              && fabs(ele->dz) < 0.1
//               && (ele->trkIso03 ) / ele->pt < 0.2
//               && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
//               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.007
               && fabs(ele->deltaPhiIn) < 0.03
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->d0) < 0.02
               && fabs(ele->dz) < 0.1
 //               && (ele->trkIso03 ) / ele->pt < 0.2
//                && (ele->emIso03 ) / ele->pt < 0.2
//                && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 

  }

  if (DenominatorType == 4) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->dz) < 0.1
              && fabs(ele->d0) < 0.02
              && (ele->trkIso03) / ele->pt < 0.2
              && (ele->emIso03) / ele->pt < 0.20
              && (ele->hadIso03) / ele->pt < 0.20
              
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && ele->HoverE < 0.10
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->dz) < 0.1
               && fabs(ele->d0) < 0.02
               && (ele->trkIso03) / ele->pt < 0.2
               && (ele->emIso03) / ele->pt < 0.20
               && (ele->hadIso03) / ele->pt < 0.20
           )
        ) {
        pass = kFALSE;
      }
    } 

  }



  if (DenominatorType == 5) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->dz) < 0.1
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
              && (ele->hadIso03) / ele->pt < 0.2

            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->dz) < 0.1
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
               && (ele->trkIso03 ) / ele->pt < 0.2
               && (ele->emIso03 ) / ele->pt < 0.2
               && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 

  }


  if (DenominatorType == 6) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->nExpHitsInner <= 0
              && passConversionVeto(ele->isConv)
              && fabs(ele->dz) < 0.1
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
              && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else  {
      if (! (  (0==0)
               && ele->nExpHitsInner <= 0
               && passConversionVeto(ele->isConv)
               && fabs(ele->dz) < 0.1
               && (ele->trkIso03 ) / ele->pt < 0.2
               && (ele->emIso03 ) / ele->pt < 0.2
               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }
    } 

  }

//Boris Denominator
  if (DenominatorType == 20) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (ele->emIso03 ) / ele->pt < 0.20
              && (ele->hadIso03) / ele->pt < 0.20
              
            )
        ) {
        pass = kFALSE;
      }      
    }

    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && ele->HoverE < 0.10
               && (ele->trkIso03) / ele->pt < 0.2
               && (ele->emIso03) / ele->pt < 0.20
               && (ele->hadIso03) / ele->pt < 0.20
           )
        ) {
        pass = kFALSE;
      }
    } 

  }




  if (DenominatorType == 100) {
    pass = kTRUE;
  }


  return pass;
}


//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs <<   runNum << " " ;
  ofs <<  lumiSec << " ";
  ofs << evtNum<< " ";
  ofs << mass<< " ";

//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
