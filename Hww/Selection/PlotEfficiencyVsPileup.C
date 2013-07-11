//root -l EWKAna/Hww/Selection/PlotEfficiencyVsPileup.C+\(130,\"\",0\)


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
#include <TGraphAsymmErrors.h>               // 3D vector class

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
Bool_t passFirstElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronCuts(const mithep::TElectron *ele, Int_t ElectronSelectionType, Double_t likelihood);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);

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

Int_t Classify(const mithep::TElectron *ele) {
  
  double eta    = fabs(ele->eta);
  double eOverP = ele->EOverP;
  double fBrem  = ele->fBrem;

  int cat = -1;
  if (ele->isEB == kTRUE) {
    if ((fBrem >= 0.12) and (eOverP > 0.9) and (eOverP < 1.2))
      cat = 0;
    else if (((eta >  .445   and eta <  .45  ) or
  	      (eta >  .79    and eta <  .81  ) or
  	      (eta > 1.137   and eta < 1.157 ) or
  	      (eta > 1.47285 and eta < 1.4744)))
      cat = 6;
    else if (!ele->isEcalDriven)
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) and (eOverP > 0.82) and (eOverP < 1.22))
      cat = 3;
    else if (eta > 1.5 and eta <  1.58)
      cat = 7;
    else if (!ele->isEcalDriven)
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else
      cat = 5;
  }

  return cat;
}


//--------------------------------------------------------------------------------------------------
bool compute_cut(double x, double et, double cut_min, double cut_max, bool gtn = false) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}

//=== MAIN MACRO =================================================================================================

void PlotEfficiencyVsPileup(Double_t mHiggs, const string Label, Int_t ElectronSelectionType) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  string label = Label;
  if (Label != "") label = "_" + Label;
  

  Double_t lumi = 35.5;              // luminosity (pb^-1)
  lumi = 1000;
  Int_t ChargeSelection = 0;
//   ChargeSelection = 1;
    

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



  Double_t fPtMaxLowerCut;
  Double_t fPtMinLowerCut;
  Double_t fDileptonMassUpperCut;
  Double_t fDeltaPhiCut;

  Double_t fHiggsMass = mHiggs;

  //Configure Higgs Mass Dependant Cuts
  if (fHiggsMass == 120) { fPtMaxLowerCut = 20.0;        fPtMinLowerCut = 10.0; 
    fDileptonMassUpperCut = 40.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 130) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 10.0;
    fDileptonMassUpperCut = 45.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 140) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 20.0;
    fDileptonMassUpperCut = 45.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 150) { fPtMaxLowerCut = 27.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 160) { fPtMaxLowerCut = 30.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 170) { fPtMaxLowerCut = 34.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 180) { fPtMaxLowerCut = 36.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 60.0; fDeltaPhiCut = 70.0;
  }
  if (fHiggsMass == 190) { fPtMaxLowerCut = 38.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 80.0; fDeltaPhiCut = 90.0;
  }
  if (fHiggsMass == 200) { fPtMaxLowerCut = 40.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 90.0; fDeltaPhiCut = 100.0;
  }
  if (fHiggsMass == 210) { fPtMaxLowerCut = 44.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 110.0; fDeltaPhiCut = 110.0;
  }
  if (fHiggsMass == 220) { fPtMaxLowerCut = 48.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 120.0; fDeltaPhiCut = 120.0;
  }
  if (fHiggsMass == 230) { fPtMaxLowerCut = 52.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 130.0; fDeltaPhiCut = 130.0;
  }
  if (fHiggsMass == 250) { fPtMaxLowerCut = 55.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 150.0; fDeltaPhiCut = 140.0;
  }
  if (fHiggsMass == 300) { fPtMaxLowerCut = 70.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 200.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 350) { fPtMaxLowerCut = 80.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 250.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 400) { fPtMaxLowerCut = 90.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 300.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 450) { fPtMaxLowerCut = 110.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 350.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 500) { fPtMaxLowerCut = 120.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 400.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 550) { fPtMaxLowerCut = 130.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 450.0; fDeltaPhiCut = 175.0;
  }
  if (fHiggsMass == 600) { fPtMaxLowerCut = 140.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 500.0; fDeltaPhiCut = 175.0;
  }


  vector<vector<string> > inputFiles;
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_full-d22_noskim.root");
//     inputFiles.push_back(vector<string>());
//     inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12_noskim_normalized.root");
     inputFiles.push_back(vector<string>());
     inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-h130ww2l-gf-z2-v8-pu11_noskim_normalized.root");
     inputFiles.push_back(vector<string>());
     inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ggww-z2-v8-pu11_noskim_normalized.root");


//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wz-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zz-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-tt-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-stop-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ttop-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wtop-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm1020-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt1020-powheg-c10-v8-pu11_noskim_normalized.root");
     inputFiles.push_back(vector<string>());
     inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjetsl-z2-v8-pu11-skimmed_TwoRecoLeptonSkim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-wjets-mg-z2-v12-pu-skimmed_TwoRecoLeptonSkim_normalized.root");



//2010 PU
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-z2-v8-pu_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjetsl-z2-v8-pu-skimmed_TwoRecoLeptonSkim_normalized.root");
  


  vector<string> processNames;
//   processNames.push_back("Data");
//    processNames.push_back("Hww130NoPU");
  processNames.push_back("Hww130");
  processNames.push_back("WW"); 
//   processNames.push_back("WZ/ZZ");
//   processNames.push_back("top");
//   processNames.push_back("Z");
  processNames.push_back("WJets");
  
  vector<Int_t> processColors;
  processColors.push_back(kRed);
  processColors.push_back(kBlue);
  processColors.push_back(kBlue);
  

  assert(processNames.size() == inputFiles.size());
  assert(processColors.size() == inputFiles.size());
 //  assert(processNames.size() == 1);
  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector<Double_t> ElectronHOverEEfficiencyNumerator;
  vector<Double_t> ElectronL1CorrectedHOverEEfficiencyNumerator;
  vector<Double_t> ElectronHOverEEfficiencyDenominator;

  vector<Double_t> ElectronIsolationEfficiencyBarrelNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyBarrelNumerator;
  vector<Double_t> ElectronIsolationEfficiencyBarrelDenominator;
  vector<Double_t> ElectronIsolationEfficiencyEndcapNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyEndcapNumerator;
  vector<Double_t> ElectronIsolationEfficiencyEndcapDenominator;
  vector<Double_t> MuonIsolationEfficiencyNumerator;
  vector<Double_t> MuonL1CorrectedIsolationEfficiencyNumerator;
  vector<Double_t> MuonIsolationEfficiencyDenominator;
  for(UInt_t i=0 ; i<20; ++i) {
    ElectronHOverEEfficiencyNumerator.push_back(0.0);
    ElectronL1CorrectedHOverEEfficiencyNumerator.push_back(0.0);
    ElectronHOverEEfficiencyDenominator.push_back(0.0);
    
    ElectronIsolationEfficiencyBarrelNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyBarrelNumerator.push_back(0.0);
    ElectronIsolationEfficiencyBarrelDenominator.push_back(0.0);
    ElectronIsolationEfficiencyEndcapNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyEndcapNumerator.push_back(0.0);
    ElectronIsolationEfficiencyEndcapDenominator.push_back(0.0);
    MuonIsolationEfficiencyNumerator.push_back(0.0);
    MuonL1CorrectedIsolationEfficiencyNumerator.push_back(0.0);
    MuonIsolationEfficiencyDenominator.push_back(0.0);
  }

  vector <TH1F*>  Electron_relIso_Barrel;
  vector <TH1F*>  Electron_relIso_Endcap;
  vector <TH1F*>  Muon_relIso;
  vector <TH1F*>  Electron_relIsoL1Corrected_Barrel;
  vector <TH1F*>  Electron_relIsoL1Corrected_Endcap;
  vector <TH1F*>  Muon_relIsoL1Corrected;
  vector <TH1F*>  Electron_relIso04_Barrel;
  vector <TH1F*>  Electron_relIso04_Endcap;
  vector <TH1F*>  Muon_relIso05;
  vector <TH1F*>  Electron_relIso04L1Corrected_Barrel;
  vector <TH1F*>  Electron_relIso04L1Corrected_Endcap;
  vector <TH1F*>  Muon_relIso05L1Corrected;

  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpElectron_relIso_Barrel = new TH1F((string("Electron_relIso_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIso_Endcap = new TH1F((string("Electron_relIso_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpMuon_relIso = new TH1F((string("Muon_relIso_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIsoL1Corrected_Barrel = new TH1F((string("Electron_relIsoL1Corrected_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIsoL1Corrected_Endcap = new TH1F((string("Electron_relIsoL1Corrected_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpMuon_relIsoL1Corrected = new TH1F((string("Muon_relIsoL1Corrected_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIso04_Barrel = new TH1F((string("Electron_relIso04_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIso04_Endcap = new TH1F((string("Electron_relIso04_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpMuon_relIso05 = new TH1F((string("Muon_relIso05_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIso04L1Corrected_Barrel = new TH1F((string("Electron_relIso04L1Corrected_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpElectron_relIso04L1Corrected_Endcap = new TH1F((string("Electron_relIso04L1Corrected_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    TH1F *tmpMuon_relIso05L1Corrected = new TH1F((string("Muon_relIso05L1Corrected_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
    Electron_relIso_Barrel.push_back(tmpElectron_relIso_Barrel);
    Electron_relIso_Endcap.push_back(tmpElectron_relIso_Endcap);
    Muon_relIso.push_back(tmpMuon_relIso);
    Electron_relIsoL1Corrected_Barrel.push_back(tmpElectron_relIsoL1Corrected_Barrel);
    Electron_relIsoL1Corrected_Endcap.push_back(tmpElectron_relIsoL1Corrected_Endcap);
    Muon_relIsoL1Corrected.push_back(tmpMuon_relIsoL1Corrected);
   Electron_relIso04_Barrel.push_back(tmpElectron_relIso04_Barrel);
    Electron_relIso04_Endcap.push_back(tmpElectron_relIso04_Endcap);
    Muon_relIso05.push_back(tmpMuon_relIso05);
    Electron_relIso04L1Corrected_Barrel.push_back(tmpElectron_relIso04L1Corrected_Barrel);
    Electron_relIso04L1Corrected_Endcap.push_back(tmpElectron_relIso04L1Corrected_Endcap);
    Muon_relIso05L1Corrected.push_back(tmpMuon_relIso05L1Corrected);
  }

  string tmplabel = "";
  TH1F *electron_relIso_Barrel_LowPU = new TH1F((string("electron_relIso_Barrel_LowPU")+tmplabel).c_str(), "; relIso; Fraction of Events ", 200, -1.0, 1.0);
  TH1F *electron_relIso_Barrel_HighPU = new TH1F((string("electron_relIso_Barrel_HighPU")+tmplabel).c_str(), "; relIso; Fraction of Events ", 200, -1.0, 1.0);
  TH1F *electron_relIsoL1Corrected_Barrel_LowPU = new TH1F((string("electron_relIsoL1Corrected_Barrel_LowPU")+tmplabel).c_str(), "; relIso; Fraction of Events ", 200, -1.0, 1.0);
  TH1F *electron_relIsoL1Corrected_Barrel_HighPU = new TH1F((string("electron_relIsoL1Corrected_Barrel_HighPU")+tmplabel).c_str(), "; relIso; Fraction of Events ", 200, -1.0, 1.0);


  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TFile *inputFile=0;
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
  rlrm.AddJSONFile("Cert_132440-144114_7TeV_Sep17ReReco_Collisions10_JSON.txt"); 
  hasJSON = kFALSE;
  
  
  for (int q = 0; q<inputFiles.size() ; ++q) { 
    for (int f = 0; f < inputFiles[q].size() ; ++f) {
      //********************************************************
      // Get Tree
      //********************************************************
      cout << "Reading File " << inputFiles[q][f] << endl;
      eventTree = getTreeFromFile(inputFiles[q][f].c_str(),"Events"); 
      TBranch *infoBr;
      TBranch *electronBr;
      TBranch *muonBr;
      TBranch *jetBr;


      //*****************************************************************************************
      //Loop over muon Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
        if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
		
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
        //for the skimmed input, I already required the HLT bits.
        //    if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;

        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear(); 
        muonArr->Clear(); 
        jetArr->Clear(); 
        electronBr->GetEntry(ientry);
        muonBr->GetEntry(ientry);
        jetBr->GetEntry(ientry);

        //event weight
        Double_t eventweight = info->eventweight * lumi;
//         cout << "eventweight : " << info->eventweight  << " " << eventweight << endl;
//         eventweight = 1;
        //********************************************************
        // TcMet
        //********************************************************
        TVector3 met;        
        if(info->tcMEx!=0 || info->tcMEy!=0) {       
          met.SetXYZ(info->tcMEx, info->tcMEy, 0);
        }

  

        //********************************************************
        // WW Event Selection
        //********************************************************

        Int_t NSoftMuons = 0;
        Int_t NLeptons = 0;
        vector<Int_t> leptonType;
        vector<Int_t> leptonIndex;
        vector<Double_t> leptonPt;
        vector<Double_t> leptonEta;
        vector<Double_t> leptonPhi;
        vector<Int_t> leptonCharge;

        Int_t NJets = 0;
        const mithep::TJet *leadingJet = 0;

        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
          if ( (0==0)
               &&
               passMuonCuts(mu)
               &&
               fabs(mu->eta) < 2.4
               && 
               mu->pt > 0.0
            ) {
            leptonPt.push_back(mu->pt);
            leptonEta.push_back(mu->eta);
            leptonPhi.push_back(mu->phi);
            leptonType.push_back(13);
            leptonIndex.push_back(i);  
            leptonCharge.push_back(mu->q);
          }
        }
        //soft muons
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
          Bool_t isCleanMuon = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isCleanMuon = kTRUE; 
              break;
            }
          }
          if ( mu->pt > 3.0
               && (mu->qualityBits & kTMLastStationAngTight)
               && mu->nTkHits > 10
               && fabs(mu->d0) < 0.2
               && (mu->typeBits & kTracker)
               && !isCleanMuon
            ) {
            NSoftMuons++;
          }
        }
        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
          Bool_t isMuonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isMuonOverlap = kTRUE; 
              break;
            }        
          }

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


          if ( (0==0)
               && 
//                ( (ele->pt > 20 && passFirstElectronCuts(ele)) || 
//                  (ele->pt <= 20 && passElectronCuts(ele, ElectronSelectionType, likelihood))
//                  )
               passFirstElectronCuts(ele)
               &&
               fabs(ele->eta) < 2.5
               && 
               ele->pt > 0.0
               &&
               !isMuonOverlap
            ) {
            leptonPt.push_back(ele->pt);
            leptonEta.push_back(ele->eta);
            leptonPhi.push_back(ele->phi);
            leptonType.push_back(11);
            leptonIndex.push_back(i);
            leptonCharge.push_back(ele->q);
          }
        }


        //sort leptons
        Int_t tempType;
        Int_t tempIndex;
        Double_t tempPt;
        Double_t tempEta;
        Double_t tempPhi;
        Int_t tempCharge;
        for (int l=0; l<leptonIndex.size(); l++) {
          for (int k=0; k < leptonIndex.size() - 1; k++) {
            if (leptonPt[k+1] > leptonPt[k]) {
              tempType = leptonType[k];
              tempIndex = leptonIndex[k];
              tempPt = leptonPt[k];
              tempEta = leptonEta[k];
              tempPhi = leptonPhi[k];
              tempCharge = leptonCharge[k];
          
              leptonType[k] = leptonType[k+1];
              leptonIndex[k] = leptonIndex[k+1];
              leptonPt[k] = leptonPt[k+1];
              leptonEta[k] = leptonEta[k+1];
              leptonPhi[k] = leptonPhi[k+1];
              leptonCharge[k] = leptonCharge[k+1];

              leptonType[k+1] = tempType;
              leptonIndex[k+1] = tempIndex;
              leptonPt[k+1] = tempPt;
              leptonEta[k+1] = tempEta;
              leptonPhi[k+1] = tempPhi;
              leptonCharge[k+1] = tempCharge;
          
            }
          }
        }

        double maxBtag = -99999;
        for(Int_t i=0; i<jetArr->GetEntries(); i++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
              leptonOverlap = kTRUE;
            }
          }

          if (!leptonOverlap) {
            if (jet->pt > 24.5 && fabs(jet->eta) < 5.0 ) {
              if (!leadingJet || jet->pt > leadingJet->pt) {
                leadingJet = jet;
              }
              NJets++;
            } else {
              if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
            }
          }
        }


  


        //******************************************************************************
        //dilepton preselection
        //******************************************************************************
        if (leptonPt.size() < 2) continue;
        if (!(leptonPt[0] > 20.0 && leptonPt[1] > 10.0)) continue;

        for(int i = 0; i < leptonPt.size(); ++i) {
          for(int j = i+1; j < leptonPt.size(); ++j) {

            //require opposite sign
            if ((ChargeSelection == 0 && leptonCharge[i] == leptonCharge[j]) || (ChargeSelection == 1 && leptonCharge[0] != leptonCharge[j])) continue;

            Int_t finalState = -1;
            if (leptonType[i] == 11 && leptonType[j] == 11) {
              finalState = 0;
            } else if (leptonType[i] == 13 && leptonType[j] == 13) {
              finalState = 1;
            } else if (leptonType[i] == 11 && leptonType[j] == 13) {
              finalState = 2;
            } else if (leptonType[i] == 13 && leptonType[j] == 11) {
              finalState = 3;
            }


            //***********************************************************************************************
            //|Z_vert-Z_l| maximum
            //***********************************************************************************************
            double zDiffMax = 0.0;
       
            double dz_i = 0;
            if (leptonType[0] == 11) {
              dz_i = ((mithep::TElectron*)((*electronArr)[leptonIndex[i]]))->dz;
            } else {
              dz_i = ((mithep::TMuon*)((*muonArr)[leptonIndex[i]]))->dz;
            }
            if (dz_i > zDiffMax) zDiffMax = dz_i;
    
            double dz_j;
            if (leptonType[j] == 11) {
              dz_j = ((mithep::TElectron*)((*electronArr)[leptonIndex[j]]))->dz;
            } else {
              dz_j = ((mithep::TMuon*)((*muonArr)[leptonIndex[j]]))->dz;
            }
            if (dz_j > zDiffMax) zDiffMax = dz_j;
            //szDiffMax = fabs(dz_i - dz_j);

            //******************************************************************************
            //construct event variables
            //******************************************************************************
            mithep::FourVectorM lepton1;
            mithep::FourVectorM lepton2;
            if (leptonType[i] == 11) {
              lepton1.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], 0.51099892e-3 );
            } else {
              lepton1.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], 105.658369e-3 );
            }
            if (leptonType[j] == 11) {
              lepton2.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], 0.51099892e-3 );
            } else {
              lepton2.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], 105.658369e-3 );
            }
            mithep::FourVectorM dilepton = lepton1+lepton2;

            double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(lepton1.Phi(), 
                                                                 lepton2.Phi())* 180.0 / TMath::Pi();    
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

//                 if (!(finalState == 0)) continue;
//               if (!(finalState == 0 || finalState == 3)) continue;
//              if (!(finalState == 1 || finalState == 2)) continue;
      //      if (!(lepton2.Pt() < 0 && lepton2.Pt() >= 0) continue;
//                if (!(lepton2.Pt() >= 20)) continue;
//                   if (!(lepton2.Pt() < 20 && lepton2.Pt() >= 10)) continue;
//              if (!(lepton2.Pt() < 15 && lepton2.Pt() >= 10)) continue;
//             if (!(lepton2.Pt() < 10 && lepton2.Pt() >= 5)) continue;
   
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 14;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                                   false, false, false, false};

            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 10.0) passCut[0] = true;

            if(zDiffMax < 1.0)                    passCut[1] = true;
  
            if(met.Pt()    > 20.0)               passCut[2] = true;
  
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
              if(METdeltaPhilEt > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(METdeltaPhilEt > 20) passCut[5] = true;
            }

            if(NJets     < 1)              passCut[6] = true;
        

            if (NSoftMuons == 0 )      passCut[7] = true;

            if (!(leptonPt.size() >= 3 && leptonPt[2] > 10.0)) passCut[8] = true;

            if(maxBtag < 2.1)                     passCut[9] = true;

            if (lepton1.Pt() > fPtMaxLowerCut) passCut[10] = true;
            if (lepton2.Pt() > fPtMinLowerCut) passCut[11] = true;
            if (dilepton.M() < fDileptonMassUpperCut)   passCut[12] = true;
            if (deltaPhiLeptons < fDeltaPhiCut) passCut[13] = true;

//             if (!(passCut[12] && passCut[13])) continue;

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    



            //********************************************************
            // Lepton Efficiencies
            //********************************************************
            if ((0==0) 
//                 && passCut[3]
//                 && info->nPUEvents > 14
              ) {
              Int_t index = info->nPUEvents; if (index > 19) index = 19;
              
              Double_t PUIsolationEnergy = info->PileupEnergyDensity * 3.14159 * pow(0.3,2) * 1.0 ;
              Double_t PUIsolationEnergy04Cone = info->PileupEnergyDensity * 3.14159 * pow(0.4,2) * 1.0 ;
              Double_t PUIsolationEnergy05Cone = info->PileupEnergyDensity * 3.14159 * pow(0.5,2) * 1.0 ;
              Double_t PUHOverE = info->PileupEnergyDensity * 3.14159 * pow(0.15,2) * 1.0 ;

//               if (leptonPhi.size() == 0)  continue;
//               const mithep::TMuon *firstMu = 0;    
//               if (leptonType[0] == 13) firstMu = (mithep::TMuon*)((*muonArr)[leptonIndex[0]]);

//               if (leptonType[0] == 11 || (leptonType[0] == 13 && firstMu->isMCReal )) {
        
                for(Int_t i=0; i<muonArr->GetEntries(); i++) {
                  const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);    

                  Bool_t leptonOverlap = kFALSE;
                  if (mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[0],leptonEta[0]) < 0.3) {
                    leptonOverlap = kTRUE;
                  }

                  if ( (0==0)
                       &&
                       fabs(mu->eta) < 2.4
                       && 
                       mu->pt > 10.0
                       &&
                       ( (processNames[q] != "WJets" && mu->isMCReal) || (processNames[q] == "WJets" && !mu->isMCReal) )
                       &&
                       !leptonOverlap
                       &&
                       passMuonCuts(mu)
                    ) {
                    MuonIsolationEfficiencyDenominator[index]++;  
                    Muon_relIso[q]->Fill( (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt);
                    Muon_relIsoL1Corrected[q]->Fill((mu->trkIso03 + TMath::Max(mu->emIso03 + mu->hadIso03 - PUIsolationEnergy,0.0)) / mu->pt);
                    Muon_relIso05[q]->Fill( (mu->trkIso05 + mu->emIso05 + mu->hadIso05) / mu->pt);
                    Muon_relIso05L1Corrected[q]->Fill((mu->trkIso05 + TMath::Max(mu->emIso05 + mu->hadIso05 - PUIsolationEnergy05Cone,0.0)) / mu->pt);
                    if ((mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15) {
//               if ((mu->trkIso05 + mu->emIso05 + mu->hadIso05) / mu->pt < 0.15) {
                      MuonIsolationEfficiencyNumerator[index]++;
                    }              
                    if ((mu->trkIso03 + TMath::Max(mu->emIso03 + mu->hadIso03 - PUIsolationEnergy,0.0)) / mu->pt < 0.13) {
//               if ((mu->trkIso05 + TMath::Max(mu->emIso05 + mu->hadIso05 - PUIsolationEnergy05Cone,0.0)) / mu->pt < 0.15) {
                      MuonL1CorrectedIsolationEfficiencyNumerator[index]++;
                    }    
            
                  }
                }
 //              }

              const mithep::TElectron *firstEle = 0;
              if (leptonPhi.size() > 0) 
                if (leptonType[0] == 11) firstEle = (mithep::TElectron*)((*electronArr)[leptonIndex[0]]);


              if (leptonPhi.size() == 0 || (leptonPhi.size() > 0 && !(leptonType[0] == 13 || (leptonType[0] == 11 && firstEle->isMCReal)))) continue;
              if (info->VGammaEvent) continue;

              for(Int_t i=0; i<electronArr->GetEntries(); i++) {
                const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
   
                Bool_t leptonOverlap = kFALSE;
                if (mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[0],leptonEta[0]) < 0.3) {
                  leptonOverlap = kTRUE;
                }
                
//                 Double_t relIso = (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03)/ele->pt;
//                 if (fabs(ele->eta) > 1.5) relIso = (ele->trkIso03 + ele->emIso03  + ele->hadIso03)/ele->pt;
//                 Double_t relIsoL1Corrected = (ele->trkIso03 + TMath::Max(TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt;
//                 if (fabs(ele->eta) > 1.5) relIsoL1Corrected = (ele->trkIso03 + TMath::Max(ele->emIso03 + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt;
//                 Double_t relIso04 = (ele->trkIso04 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04)/ele->pt; 
//                 if (fabs(ele->eta) > 1.5) relIso04 = (ele->trkIso04 + ele->emIso04  + ele->hadIso04)/ele->pt;
//                 Double_t relIso04L1Corrected = (ele->trkIso04 + TMath::Max(TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 - PUIsolationEnergy04Cone,0.0)) / ele->pt;
//                 if (fabs(ele->eta) > 1.5) relIso04L1Corrected = (ele->trkIso04 + TMath::Max(ele->emIso04 + ele->hadIso04 - PUIsolationEnergy04Cone,0.0)) / ele->pt;
                Double_t relIso = (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03)/ele->pt;
                if (fabs(ele->eta) > 1.5) relIso = (ele->trkIso03 + ele->emIso03  + ele->hadIso03)/ele->pt;
                Double_t relIsoL1Corrected = (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy) / ele->pt;
                if (fabs(ele->eta) > 1.5) relIsoL1Corrected = (ele->trkIso03 + ele->emIso03 + ele->hadIso03 - PUIsolationEnergy) / ele->pt;
                Double_t relIso04 = (ele->trkIso04 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04)/ele->pt; 
                if (fabs(ele->eta) > 1.5) relIso04 = (ele->trkIso04 + ele->emIso04  + ele->hadIso04)/ele->pt;
                Double_t relIso04L1Corrected = (ele->trkIso04 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 - PUIsolationEnergy04Cone) / ele->pt;
                if (fabs(ele->eta) > 1.5) relIso04L1Corrected = (ele->trkIso04 + ele->emIso04 + ele->hadIso04 - PUIsolationEnergy04Cone) / ele->pt;
           
                if ( (0==0)
                     && 
                     fabs(ele->eta) < 2.5
                     && 
                     ele->pt > 10.0
                     &&
                     ( (processNames[q] != "WJets" && ele->isMCReal) || (processNames[q] == "WJets" && !ele->isMCReal) )
                     &&
                     !leptonOverlap
                     && 
                     passFirstElectronCuts(ele)
                  ) {

                  if (fabs(ele->eta) < 1.5) {
                    ElectronIsolationEfficiencyBarrelDenominator[index]++;
                    Electron_relIso_Barrel[q]->Fill(relIso);
                    Electron_relIsoL1Corrected_Barrel[q]->Fill(relIsoL1Corrected);
                    Electron_relIso04_Barrel[q]->Fill(relIso04);
                    Electron_relIso04L1Corrected_Barrel[q]->Fill(relIso04L1Corrected);
                    if (relIso < 0.1) {
                      ElectronIsolationEfficiencyBarrelNumerator[index]++;              
                    }
                    if (relIsoL1Corrected < 0.065) {
                      ElectronL1CorrectedIsolationEfficiencyBarrelNumerator[index]++;              
                    }   
                  } else {
                    ElectronIsolationEfficiencyEndcapDenominator[index]++;
                    Electron_relIso_Endcap[q]->Fill(relIso);
                    Electron_relIsoL1Corrected_Endcap[q]->Fill(relIsoL1Corrected);
                    Electron_relIso04_Endcap[q]->Fill(relIso04);
                    Electron_relIso04L1Corrected_Endcap[q]->Fill(relIso04L1Corrected);
                    if (relIso < 0.1) {
                      ElectronIsolationEfficiencyEndcapNumerator[index]++;              
                    }
                    if (relIsoL1Corrected < 0.05) {
                      ElectronL1CorrectedIsolationEfficiencyEndcapNumerator[index]++;              
                    }   
                  }

                  Double_t HOverECut = (fabs(ele->eta) < 1.5)?0.04:0.025;
                  ElectronHOverEEfficiencyDenominator[index]++;
                  if (relIso < 0.1) {
                    if (ele->HoverE < HOverECut) {
                      ElectronHOverEEfficiencyNumerator[index]++;                     
                    }
                    if ((ele->HoverE*ele->e - PUHOverE)/ele->e  < HOverECut) {
                      ElectronL1CorrectedHOverEEfficiencyNumerator[index]++;                     
                    }
                  }

                  if (fabs(ele->eta) < 1.5) {
                    if (index < 4) {                
                      electron_relIso_Barrel_LowPU->Fill((ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 ) / ele->pt);
                      electron_relIsoL1Corrected_Barrel_LowPU->Fill((ele->trkIso03 + TMath::Max(TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt);
                    } 
                    if (index > 10 && index < 15) {
                      electron_relIso_Barrel_HighPU->Fill((ele->trkIso03 + ele->emIso03 + ele->hadIso03 ) / ele->pt);
                      electron_relIsoL1Corrected_Barrel_HighPU->Fill((ele->trkIso03 + TMath::Max(ele->emIso03 + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt);   
                    }
                  }

                }
              }

              
            }









 
            if ( (0==0) 
//                  && passCut[6]

                 && (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10] && passCut[11])
//                 && passCut[12]
              ) {
//               if (passCut[12]) {
//               }
            }


            if (passAllCuts 
                //  && lepton2.Pt() < 20.0 
              ) {     
            }

          }
        }
  

      } //end loop over data     
    }
  }


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================


  const int nPoints = ElectronIsolationEfficiencyBarrelNumerator.size();
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.5;     
  }


  double ElectronIsolationBarrelEff[nPoints];
  double ElectronIsolationBarrelEffErrLow[nPoints];
  double ElectronIsolationBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationBarrelEff[i] = ratio;
    ElectronIsolationBarrelEffErrLow[i] = errLow;
    ElectronIsolationBarrelEffErrHigh[i] = errHigh;
    cout << "ElectronIsolationBarrelEff " << i << " : " << ElectronIsolationBarrelEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationBarrelEff, NPileupError, NPileupError, ElectronIsolationBarrelEffErrLow, ElectronIsolationBarrelEffErrHigh);
  ElectronIsolationBarrelEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedIsolationBarrelEff[nPoints];
  double ElectronL1CorrectedIsolationBarrelEffErrLow[nPoints];
  double ElectronL1CorrectedIsolationBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedIsolationBarrelEff[i] = ratio;
    ElectronL1CorrectedIsolationBarrelEffErrLow[i] = errLow;
    ElectronL1CorrectedIsolationBarrelEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedIsolationBarrelEff " << i << " : " << ElectronL1CorrectedIsolationBarrelEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedIsolationBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationBarrelEff,NPileupError, NPileupError, ElectronL1CorrectedIsolationBarrelEffErrLow, ElectronL1CorrectedIsolationBarrelEffErrHigh);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);



  double ElectronIsolationEndcapEff[nPoints];
  double ElectronIsolationEndcapEffErrLow[nPoints];
  double ElectronIsolationEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationEndcapEff[i] = ratio;
    ElectronIsolationEndcapEffErrLow[i] = errLow;
    ElectronIsolationEndcapEffErrHigh[i] = errHigh;
    cout << "ElectronIsolationEndcapEff " << i << " : " << ElectronIsolationEndcapEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationEndcapEff, NPileupError, NPileupError, ElectronIsolationEndcapEffErrLow, ElectronIsolationEndcapEffErrHigh);
  ElectronIsolationEndcapEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedIsolationEndcapEff[nPoints];
  double ElectronL1CorrectedIsolationEndcapEffErrLow[nPoints];
  double ElectronL1CorrectedIsolationEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedIsolationEndcapEff[i] = ratio;
    ElectronL1CorrectedIsolationEndcapEffErrLow[i] = errLow;
    ElectronL1CorrectedIsolationEndcapEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedIsolationEndcapEff " << i << " : " << ElectronL1CorrectedIsolationEndcapEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedIsolationEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationEndcapEff,NPileupError, NPileupError, ElectronL1CorrectedIsolationEndcapEffErrLow, ElectronL1CorrectedIsolationEndcapEffErrHigh);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);




  double ElectronHOverEEff[nPoints];
  double ElectronHOverEEffErrLow[nPoints];
  double ElectronHOverEEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronHOverEEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronHOverEEff[i] = ratio;
    ElectronHOverEEffErrLow[i] = errLow;
    ElectronHOverEEffErrHigh[i] = errHigh;
    cout << "ElectronHOverEEff " << i << " : " << ElectronHOverEEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronHOverEEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronHOverEEff, NPileupError, NPileupError, ElectronHOverEEffErrLow, ElectronHOverEEffErrHigh);
  ElectronHOverEEffVsNPileup->SetMarkerColor(kRed);
  ElectronHOverEEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronHOverEEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedHOverEEff[nPoints];
  double ElectronL1CorrectedHOverEEffErrLow[nPoints];
  double ElectronL1CorrectedHOverEEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedHOverEEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedHOverEEff[i] = ratio;
    ElectronL1CorrectedHOverEEffErrLow[i] = errLow;
    ElectronL1CorrectedHOverEEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedHOverEEff " << i << " : " << ElectronL1CorrectedHOverEEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedHOverEEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedHOverEEff,NPileupError, NPileupError, ElectronL1CorrectedHOverEEffErrLow, ElectronL1CorrectedHOverEEffErrHigh);
  ElectronL1CorrectedHOverEEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedHOverEEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedHOverEEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);



  double MuonIsolationEff[nPoints];
  double MuonIsolationEffErrLow[nPoints];
  double MuonIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonIsolationEff[i] = ratio;
    MuonIsolationEffErrLow[i] = errLow;
    MuonIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonIsolationEff " << i << " : " << MuonIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonIsolationEff, NPileupError, NPileupError, MuonIsolationEffErrLow, MuonIsolationEffErrHigh);
  MuonIsolationEffVsNPileup->SetMarkerColor(kRed);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double MuonL1CorrectedIsolationEff[nPoints];
  double MuonL1CorrectedIsolationEffErrLow[nPoints];
  double MuonL1CorrectedIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonL1CorrectedIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonL1CorrectedIsolationEff[i] = ratio;
    MuonL1CorrectedIsolationEffErrLow[i] = errLow;
    MuonL1CorrectedIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonL1CorrectedIsolationEff " << i << " : " << MuonL1CorrectedIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonL1CorrectedIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonL1CorrectedIsolationEff, NPileupError, NPileupError,MuonL1CorrectedIsolationEffErrLow,MuonL1CorrectedIsolationEffErrHigh  );
  MuonL1CorrectedIsolationEffVsNPileup->SetMarkerColor(kBlue);
  MuonL1CorrectedIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonL1CorrectedIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);




  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  Double_t ymin = 0.0;
  Double_t ymax = 1.1;


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  legend = new TLegend(0.70,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.0;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronIsolationBarrelEffVsNPileup, "NoCorrection", "LP");
    legend->AddEntry(ElectronL1CorrectedIsolationBarrelEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronIsolationBarrelEffVsNPileup->SetTitle("");
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitle("Number of Mixed Pileup Events");
  ElectronIsolationBarrelEffVsNPileup->Draw("AP");
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->Draw("Psame");
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs("ElectronIsolationEfficiency_Barrel_vs_NPU.gif");

  legend = new TLegend(0.70,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.0;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronIsolationEndcapEffVsNPileup, "NoCorrection", "LP");
    legend->AddEntry(ElectronL1CorrectedIsolationEndcapEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronIsolationEndcapEffVsNPileup->SetTitle("");
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitle("Number of Mixed Pileup Events");
  ElectronIsolationEndcapEffVsNPileup->Draw("AP");
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->Draw("Psame");
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs("ElectronIsolationEfficiency_Endcap_vs_NPU.gif");


  legend = new TLegend(0.70,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.0;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronHOverEEffVsNPileup, "NoCorrection", "LP");
    legend->AddEntry(ElectronL1CorrectedHOverEEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronHOverEEffVsNPileup->SetTitle("");
  ElectronHOverEEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronHOverEEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronHOverEEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronHOverEEffVsNPileup->GetXaxis()->SetTitle("Number of Mixed Pileup Events");
  ElectronHOverEEffVsNPileup->Draw("AP");
  ElectronL1CorrectedHOverEEffVsNPileup->Draw("Psame");
  ElectronHOverEEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs("ElectronHOverEEfficiency_PUStudy_vs_NPU.gif");


  legend = new TLegend(0.70,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 0.4;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(MuonIsolationEffVsNPileup, "NoCorrection", "LP");
    legend->AddEntry(MuonL1CorrectedIsolationEffVsNPileup, "FastJetCorrected", "LP");
  }
  MuonIsolationEffVsNPileup->SetTitle("");
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitle("Number of Mixed Pileup Events");
  MuonIsolationEffVsNPileup->Draw("AP");
  MuonL1CorrectedIsolationEffVsNPileup->Draw("Psame");
  MuonIsolationEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);

  legend->Draw();
  cv->SaveAs("MuonIsolationEfficiency_PUStudy_vs_NPU.gif");





  //--------------------------------------------------------------------------------------------------------------
  // Make histograms
  //==============================================================================================================
  Double_t norm = 0;
  for (int b=0; b<electron_relIso_Barrel_LowPU->GetXaxis()->GetNbins()+2; ++b) { norm += electron_relIso_Barrel_LowPU->GetBinContent(b); }
  for (int b=0; b<electron_relIso_Barrel_LowPU->GetXaxis()->GetNbins()+2; ++b) {
    electron_relIso_Barrel_LowPU->SetBinContent(b,electron_relIso_Barrel_LowPU->GetBinContent(b) / norm);
    electron_relIso_Barrel_LowPU->SetBinError(b,electron_relIso_Barrel_LowPU->GetBinError(b) / norm);
  }
  
  norm = 0;
  for (int b=0; b<electron_relIso_Barrel_HighPU->GetXaxis()->GetNbins()+2; ++b) { norm += electron_relIso_Barrel_HighPU->GetBinContent(b); }
  for (int b=0; b<electron_relIso_Barrel_HighPU->GetXaxis()->GetNbins()+2; ++b) {
    electron_relIso_Barrel_HighPU->SetBinContent(b,electron_relIso_Barrel_HighPU->GetBinContent(b) / norm);
    electron_relIso_Barrel_HighPU->SetBinError(b,electron_relIso_Barrel_HighPU->GetBinError(b) / norm);
  } 


  legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(electron_relIso_Barrel_LowPU, "NPU < 4", "LP");
    legend->AddEntry(electron_relIso_Barrel_HighPU, "NPU in [11,14]", "LP");
  }
  electron_relIso_Barrel_LowPU->SetLineColor(kRed);
  electron_relIso_Barrel_LowPU->SetMarkerColor(kRed);
  electron_relIso_Barrel_HighPU->SetLineColor(kBlue);
  electron_relIso_Barrel_HighPU->SetMarkerColor(kBlue);
  electron_relIso_Barrel_LowPU->Draw("hist");
  electron_relIso_Barrel_HighPU->Draw("histsame");
//   electron_relIso_Barrel_LowPU->SetMaximum(1.0);
//   electron_relIso_Barrel_LowPU->SetMinimum(0.0);
  legend->Draw();
  cv->SaveAs("Electron_relIso_Barrel_PUStudy.gif");


norm = 0;
  for (int b=0; b<electron_relIsoL1Corrected_Barrel_LowPU->GetXaxis()->GetNbins()+2; ++b) { norm += electron_relIsoL1Corrected_Barrel_LowPU->GetBinContent(b); }
  for (int b=0; b<electron_relIsoL1Corrected_Barrel_LowPU->GetXaxis()->GetNbins()+2; ++b) {
    electron_relIsoL1Corrected_Barrel_LowPU->SetBinContent(b,electron_relIsoL1Corrected_Barrel_LowPU->GetBinContent(b) / norm);
    electron_relIsoL1Corrected_Barrel_LowPU->SetBinError(b,electron_relIsoL1Corrected_Barrel_LowPU->GetBinError(b) / norm);
  }
  
  norm = 0;
  for (int b=0; b<electron_relIsoL1Corrected_Barrel_HighPU->GetXaxis()->GetNbins()+2; ++b) { norm += electron_relIsoL1Corrected_Barrel_HighPU->GetBinContent(b); }
  for (int b=0; b<electron_relIsoL1Corrected_Barrel_HighPU->GetXaxis()->GetNbins()+2; ++b) {
    electron_relIsoL1Corrected_Barrel_HighPU->SetBinContent(b,electron_relIsoL1Corrected_Barrel_HighPU->GetBinContent(b) / norm);
    electron_relIsoL1Corrected_Barrel_HighPU->SetBinError(b,electron_relIsoL1Corrected_Barrel_HighPU->GetBinError(b) / norm);
  } 


  legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(electron_relIsoL1Corrected_Barrel_LowPU, "NPU < 4", "LP");
    legend->AddEntry(electron_relIsoL1Corrected_Barrel_HighPU, "NPU in [11,14]", "LP");
  }
  electron_relIsoL1Corrected_Barrel_LowPU->SetLineColor(kRed);
  electron_relIsoL1Corrected_Barrel_LowPU->SetMarkerColor(kRed);
  electron_relIsoL1Corrected_Barrel_HighPU->SetLineColor(kBlue);
  electron_relIsoL1Corrected_Barrel_HighPU->SetMarkerColor(kBlue);
  electron_relIsoL1Corrected_Barrel_LowPU->Draw("hist");
  electron_relIsoL1Corrected_Barrel_HighPU->Draw("histsame");
//   electron_relIsoL1Corrected_Barrel_LowPU->SetMaximum(1.0);
//   electron_relIsoL1Corrected_Barrel_LowPU->SetMinimum(0.0);
  legend->Draw();
  cv->SaveAs("Electron_relIsoL1Corrected_Barrel_PUStudy.gif");





  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlots_LeptonEfficiency.root", "UPDATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
    file->WriteTObject(Electron_relIso_Barrel[q],Electron_relIso_Barrel[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIso_Endcap[q],Electron_relIso_Endcap[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Muon_relIso[q],Muon_relIso[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIsoL1Corrected_Barrel[q],Electron_relIsoL1Corrected_Barrel[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIsoL1Corrected_Endcap[q],Electron_relIsoL1Corrected_Endcap[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Muon_relIsoL1Corrected[q],Muon_relIsoL1Corrected[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIso04_Barrel[q],Electron_relIso04_Barrel[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIso04_Endcap[q],Electron_relIso04_Endcap[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Muon_relIso05[q],Muon_relIso05[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIso04L1Corrected_Barrel[q],Electron_relIso04L1Corrected_Barrel[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Electron_relIso04L1Corrected_Endcap[q],Electron_relIso04L1Corrected_Endcap[q]->GetName(), "WriteDelete") ;
    file->WriteTObject(Muon_relIso05L1Corrected[q],Muon_relIso05L1Corrected[q]->GetName(), "WriteDelete") ;
  }

  file->Close();
  delete file;

        
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t isMC = kFALSE;
  Bool_t pass = kFALSE;
  if (isMC) {
    if (triggerBits & kHLT_Mu9) pass = kTRUE;
    if (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) pass = kTRUE;
    if (triggerBits & kHLT_Ele15_LW_L1R) pass = kTRUE;
  } else {
    if (isMuonData) {
      if ((runNum >= 136033) && (runNum <= 147116)) {
        if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
      } 
      if ((runNum >= 136033) && (runNum <= 139980)) {
        if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
      } 
      if ((runNum >= 140058) && (runNum <= 141882)) {
        if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
      } 
      if ((runNum >= 141956) && (runNum <= 144114)) {
        if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
      }
      if ((runNum >= 146428) && (runNum <= 147116)) {
        if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
      }
      if ((runNum >= 147196) && (runNum <= 999999)) {
        if ( (triggerBits & kHLT_Mu15) ) pass = kTRUE;
      }
      if ((runNum >= 147196) && (runNum <= 148058)) {
        if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
      }
      if ((runNum >= 148819) && (runNum <= 149442)) {
        if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
      }
    } else {
      //it's electron data
      if ((runNum >= 136033) && (runNum <= 139980)) {
        if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
      } 
      if ((runNum >= 140058) && (runNum <= 141882)) {
        if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
      } 
      if ((runNum >= 141956) && (runNum <= 144114)) {
        if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
      } 
      if ((runNum >= 146428) && (runNum <= 147116)) {
        if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
      }
      if ((runNum >= 147196) && (runNum <= 148058)) {
        if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
      }
      if ((runNum >= 148819) && (runNum <= 149442)) {
        if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
      }
    }
  }

  return pass;

}


Bool_t passFirstElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
//   if (ele->pt < 20.0) pass = kFALSE;

  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
//             && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
             && fabs(ele->d0) < 0.02
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
//              && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
              && ele->nExpHitsInner <= 0
              && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
              && fabs(ele->d0) < 0.02
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



Bool_t passElectronCuts(const mithep::TElectron *ele, Int_t ElectronSelectionType, Double_t likelihood) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
 if (ElectronSelectionType == -1) {
    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
//               && ele->sigiEtaiEta < 0.01 
//               && fabs(ele->deltaEtaIn) < 0.004
//               && fabs(ele->deltaPhiIn) < 0.06
//               && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
//               && ele->nExpHitsInner <= 0
//               && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//               && fabs(ele->d0) < 0.02
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else if (fabs(ele->eta) > 1.5) {
      if (! (  (0==0)
//                && ele->sigiEtaiEta < 0.03
//                && fabs(ele->deltaEtaIn) < 0.007
//                && fabs(ele->deltaPhiIn) < 0.03
//                && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
//                && ele->nExpHitsInner <= 0
//                && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//                && fabs(ele->d0) < 0.02
            )
        ) {
        pass = kFALSE;
      }
    } else {
      pass = kFALSE;
      return pass;
    }
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
//         && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

        && (mu->nSeg > 1 || mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

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
