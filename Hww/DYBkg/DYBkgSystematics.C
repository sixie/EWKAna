//root -l EWKAna/Hww/DYBkg/DYBkgSystematics.C+\(130,\"\",0\)
//root -l EWKAna/Hww/DYBkg/DYBkgSystematics.C+\(130,\"\",1\)
//root -l EWKAna/Hww/DYBkg/DYBkgSystematics.C+\(130,\"2020Selection\",0\)
//root -l EWKAna/Hww/DYBkg/DYBkgSystematics.C+\(130,\"2020Selection\",1\)


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
#include <TH2F.h>                   // 1D histograms
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
Bool_t passElectronCuts(const mithep::TElectron *ele);
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

void DYBkgSystematics(Double_t mHiggs, const string Label, Int_t FinalState)
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  string label = Label;
  if (Label != "") label = "_" + Label;
  if (FinalState == 0) label += "_ee";
  if (FinalState == 1) label += "_mm";

  
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
  vector<string> processNames;

  if (FinalState == 0) {
    inputFiles.push_back(vector<string>());
    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zee-powheg-c10-v8-pu_noskim_normalized.root");
    processNames.push_back("Zee_2010PU");
  }
  if (FinalState == 1) {
    inputFiles.push_back(vector<string>());
    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm-powheg-c10-v8-pu_noskim_normalized.root");
    processNames.push_back("Zmm_2010PU");
  }

  assert(processNames.size() == inputFiles.size());



//   vector<Int_t> processColors;
//   processColors.push_back(kRed);
//    processColors.push_back(kBlue);
//    processColors.push_back(kRed);
//    processColors.push_back(kBlue);


//   assert(processColors.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  

  vector<vector<Double_t> > CountInZWindow_AfterTCMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterTCMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFTrackMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFTrackMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFNoFwdMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFNoFwdMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFTrackMetPFMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFTrackMetPFMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterTCMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterTCMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFTrackMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFTrackMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFNoFwdMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFNoFwdMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
  for (int q=0; q<processNames.size() ; ++q) {
    vector<Double_t> tmpCountInZWindow_AfterTCMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterTCMetCut;
    vector<Double_t> tmpCountInZWindow_AfterPFMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFMetCut;
    vector<Double_t> tmpCountInZWindow_AfterPFTrackMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFTrackMetCut;
    vector<Double_t> tmpCountInZWindow_AfterPFNoFwdMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFNoFwdMetCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFTrackMetPFMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterTCMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterTCMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFTrackMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFTrackMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFNoFwdMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut;


    for (int b=0; b < 100; ++b) {
      tmpCountInZWindow_AfterTCMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterTCMetCut.push_back(0);
      tmpCountInZWindow_AfterPFMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFMetCut.push_back(0);
      tmpCountInZWindow_AfterPFTrackMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFTrackMetCut.push_back(0);
      tmpCountInZWindow_AfterPFNoFwdMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFNoFwdMetCut.push_back(0);
      tmpCountInZWindow_AfterMinPFTrackMetPFMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetCut.push_back(0);
      tmpCountInZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(0);
      tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterTCMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterTCMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFTrackMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFTrackMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFNoFwdMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(0);
    }

    CountInZWindow_AfterTCMetCut.push_back(tmpCountInZWindow_AfterTCMetCut);
    CountOutOfZWindow_AfterTCMetCut.push_back(tmpCountOutOfZWindow_AfterTCMetCut);
    CountInZWindow_AfterPFMetCut.push_back(tmpCountInZWindow_AfterPFMetCut);
    CountOutOfZWindow_AfterPFMetCut.push_back(tmpCountOutOfZWindow_AfterPFMetCut);
    CountInZWindow_AfterPFTrackMetCut.push_back(tmpCountInZWindow_AfterPFTrackMetCut);
    CountOutOfZWindow_AfterPFTrackMetCut.push_back(tmpCountOutOfZWindow_AfterPFTrackMetCut);
    CountInZWindow_AfterPFNoFwdMetCut.push_back(tmpCountInZWindow_AfterPFNoFwdMetCut);
    CountOutOfZWindow_AfterPFNoFwdMetCut.push_back(tmpCountOutOfZWindow_AfterPFNoFwdMetCut);
    CountInZWindow_AfterMinPFTrackMetPFMetCut.push_back(tmpCountInZWindow_AfterMinPFTrackMetPFMetCut);
    CountOutOfZWindow_AfterMinPFTrackMetPFMetCut.push_back(tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetCut);
    CountInZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(tmpCountInZWindow_AfterMinPFNoFwdMetPFMetCut);
    CountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut);
    CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut);
    CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut);
    CountErrSqrInZWindow_AfterTCMetCut.push_back(tmpCountErrSqrInZWindow_AfterTCMetCut);
    CountErrSqrOutOfZWindow_AfterTCMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterTCMetCut);
    CountErrSqrInZWindow_AfterPFMetCut.push_back(tmpCountErrSqrInZWindow_AfterPFMetCut);
    CountErrSqrOutOfZWindow_AfterPFMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFMetCut);
    CountErrSqrInZWindow_AfterPFTrackMetCut.push_back(tmpCountErrSqrInZWindow_AfterPFTrackMetCut);
    CountErrSqrOutOfZWindow_AfterPFTrackMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFTrackMetCut);
    CountErrSqrInZWindow_AfterPFNoFwdMetCut.push_back(tmpCountErrSqrInZWindow_AfterPFNoFwdMetCut);
    CountErrSqrOutOfZWindow_AfterPFNoFwdMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetCut);
    CountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut);
    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut);
    CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut);
    CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut);
    CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut);
    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut);
    
  }



  vector<vector<Double_t> > CountInZWindow_AfterTCMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterTCMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFTrackMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<vector<Double_t> > CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
  for (int q=0; q<processNames.size() ; ++q) {
    vector<Double_t> tmpCountInZWindow_AfterTCMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterTCMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterPFTrackMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
    vector<Double_t> tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;


    for (int b=0; b < 100; ++b) {
      tmpCountInZWindow_AfterTCMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterTCMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterPFMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
      tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(0);
    }

    CountInZWindow_AfterTCMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterTCMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterTCMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterTCMetdeltaPhilEtCut);
    CountInZWindow_AfterPFMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterPFMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterPFMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterPFMetdeltaPhilEtCut);
    CountInZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterPFTrackMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut);
    CountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut);
    CountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut);
    CountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut);
    CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut);
    CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut);
    CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut);
    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpCountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut);
    
  }


  vector<string> CutLabel;
  CutLabel.push_back("Dilepton20/10");
  CutLabel.push_back("LeptonPt20/20");
  CutLabel.push_back("Lepton dZ");
  CutLabel.push_back("MetPreselection");
  CutLabel.push_back("MassPreselection");
  CutLabel.push_back("ZVeto");
  CutLabel.push_back("ProjectedMet");
  CutLabel.push_back("JetVeto");
  CutLabel.push_back("SoftMuonVeto");
  CutLabel.push_back("ExtraLeptonVeto");
  CutLabel.push_back("BTagVeto");
  CutLabel.push_back("Pt1Cut");
  CutLabel.push_back("Pt2Cut");
  CutLabel.push_back("MassCut");
  CutLabel.push_back("DeltaPhi");
  
  
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
  // Don't need this for MC
  //********************************************************
  
  
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


        //********************************************************
        // Met
        //********************************************************
        TVector3 met;        
        if(info->tcMEx!=0 || info->tcMEy!=0) {       
          met.SetXYZ(info->tcMEx, info->tcMEy, 0);
        }
        TVector3 pfMet;        
        TVector3 pfTrackMet;        
        TVector3 pfNoFwdMet;        
        TVector3 MinPFTrackMetPFMet;
        TVector3 MinPFNoFwdMetPFMet;
        TVector3 MinPFTrackMetPFNoFwdMetPFMet;
        if(info->pfMEx!=0 || info->pfMEy!=0) {       
          pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        }
        if(info->pfTrackMEx!=0 || info->pfTrackMEy!=0) {       
          pfTrackMet.SetXYZ(info->pfTrackMEx, info->pfTrackMEy, 0);
        }
        if(info->pfNeutralNoFwdMEx+info->pfTrackMEx!=0 || info->pfNeutralNoFwdMEy+info->pfTrackMEy!=0) {       
          pfNoFwdMet.SetXYZ(info->pfNeutralNoFwdMEx+info->pfTrackMEx, 
                            info->pfNeutralNoFwdMEy+info->pfTrackMEy, 0);
        }

        if (pfMet.Pt() < pfTrackMet.Pt()) {
          MinPFTrackMetPFMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        } else {
          MinPFTrackMetPFMet.SetXYZ(info->pfTrackMEx, info->pfTrackMEy, 0);     
        }

        if (pfMet.Pt() < pfNoFwdMet.Pt()) {
          MinPFNoFwdMetPFMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        } else {
          MinPFNoFwdMetPFMet.SetXYZ(info->pfNeutralNoFwdMEx+info->pfTrackMEx, 
                                    info->pfNeutralNoFwdMEy+info->pfTrackMEy, 0);
        }
        if (pfMet.Pt() < pfTrackMet.Pt() && pfMet.Pt() < pfNoFwdMet.Pt()) {
          MinPFTrackMetPFNoFwdMetPFMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        } 
        if (pfTrackMet.Pt() < pfMet.Pt() && pfTrackMet.Pt() < pfNoFwdMet.Pt()) {
          MinPFTrackMetPFNoFwdMetPFMet.SetXYZ(info->pfTrackMEx, info->pfTrackMEy, 0);  
        }
        if (pfNoFwdMet.Pt() < pfMet.Pt() && pfNoFwdMet.Pt() < pfTrackMet.Pt() ) {
          MinPFTrackMetPFNoFwdMetPFMet.SetXYZ(info->pfNeutralNoFwdMEx+info->pfTrackMEx, 
                                              info->pfNeutralNoFwdMEy+info->pfTrackMEy, 0);
        }

        //********************************************************
        // TcMet
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
               passElectronCuts(ele)
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

            //*********************************************************************************************
            //Select Final State
            //*********************************************************************************************
            if (!(finalState == FinalState)) continue;


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


            double deltaPhiPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(pfMet.Phi(), lepton1.Phi()),
                                             mithep::MathUtils::DeltaPhi(pfMet.Phi(), lepton2.Phi())};
            double minDeltaPhiPFMetLepton = (deltaPhiPFMetLepton[0] < deltaPhiPFMetLepton[1])?
              deltaPhiPFMetLepton[0]:deltaPhiPFMetLepton[1];
            double PFMETdeltaPhilEt = pfMet.Pt();
            if(minDeltaPhiPFMetLepton < TMath::Pi()/2.)
              PFMETdeltaPhilEt = PFMETdeltaPhilEt * sin(minDeltaPhiPFMetLepton);
            

            double deltaPhiPFTrackMetLepton[2] = {mithep::MathUtils::DeltaPhi(pfTrackMet.Phi(), lepton1.Phi()),
                                                  mithep::MathUtils::DeltaPhi(pfTrackMet.Phi(), lepton2.Phi())};
            double minDeltaPhiPFTrackMetLepton = (deltaPhiPFTrackMetLepton[0] < deltaPhiPFTrackMetLepton[1])?
              deltaPhiPFTrackMetLepton[0]:deltaPhiPFTrackMetLepton[1];
            double PFTrackMETdeltaPhilEt = pfTrackMet.Pt();
            if(minDeltaPhiPFTrackMetLepton < TMath::Pi()/2.)
              PFTrackMETdeltaPhilEt = PFTrackMETdeltaPhilEt * sin(minDeltaPhiPFTrackMetLepton);


            double deltaPhiPFNoFwdMetLepton[2] = {mithep::MathUtils::DeltaPhi(pfNoFwdMet.Phi(), lepton1.Phi()),
                                                  mithep::MathUtils::DeltaPhi(pfNoFwdMet.Phi(), lepton2.Phi())};
            double minDeltaPhiPFNoFwdMetLepton = (deltaPhiPFNoFwdMetLepton[0] < deltaPhiPFNoFwdMetLepton[1])?
              deltaPhiPFNoFwdMetLepton[0]:deltaPhiPFNoFwdMetLepton[1];
            double PFNoFwdMETdeltaPhilEt = pfNoFwdMet.Pt();
            if(minDeltaPhiPFNoFwdMetLepton < TMath::Pi()/2.)
              PFNoFwdMETdeltaPhilEt = PFNoFwdMETdeltaPhilEt * sin(minDeltaPhiPFNoFwdMetLepton);


            double deltaPhiMinPFTrackMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFTrackMetPFMet.Phi(), lepton1.Phi()),
                                                  mithep::MathUtils::DeltaPhi(MinPFTrackMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFTrackMetPFMetLepton = (deltaPhiMinPFTrackMetPFMetLepton[0] < deltaPhiMinPFTrackMetPFMetLepton[1])?
              deltaPhiMinPFTrackMetPFMetLepton[0]:deltaPhiMinPFTrackMetPFMetLepton[1];
            double MinPFTrackMetPFMetdeltaPhilEt = TMath::Min(PFTrackMETdeltaPhilEt,PFMETdeltaPhilEt);

            double deltaPhiMinPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                          mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFNoFwdMetPFMetLepton = (deltaPhiMinPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFNoFwdMetPFMetLepton[1])?
              deltaPhiMinPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFNoFwdMetPFMetLepton[1];
            double MinPFNoFwdMetPFMetdeltaPhilEt = TMath::Min(PFNoFwdMETdeltaPhilEt,PFMETdeltaPhilEt);

            double deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                                    mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton = (deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1])?              
              deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1];
            double MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = TMath::Min(TMath::Min(PFTrackMETdeltaPhilEt, PFNoFwdMETdeltaPhilEt),PFMETdeltaPhilEt);

   
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 14;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                                   false, false, false, false};
            assert(CutLabel.size() == nCuts+1);

            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 20.0) passCut[0] = true;

            if(zDiffMax < 1.0)                    passCut[1] = true;
  
            if(pfMet.Pt()    > 20.0)               passCut[2] = true;
  
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0) passCut[4] = true;   
              if(PFMETdeltaPhilEt > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(PFMETdeltaPhilEt > 20) passCut[5] = true;
            }

            if(NJets     < 1)              passCut[6] = true;
        

            if (NSoftMuons == 0 )      passCut[7] = true;

            if (!(leptonPt.size() >= 3 && leptonPt[2] > 10.0)) passCut[8] = true;

            if(maxBtag < 2.1)                     passCut[9] = true;

            if (lepton1.Pt() > fPtMaxLowerCut) passCut[10] = true;
            if (lepton2.Pt() > fPtMinLowerCut) passCut[11] = true;
            if (dilepton.M() < fDileptonMassUpperCut)   passCut[12] = true;
            if (deltaPhiLeptons < fDeltaPhiCut) passCut[13] = true;


            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
            //*********************************************************************************************
            //R_in/out Plots
            //*********************************************************************************************  
            if (passCut[0] && passCut[1] && passCut[3] && passCut[6] 
                && passCut[7] && passCut[8] && passCut[9]
              ) {

              for (UInt_t b=0; b < 100; ++b) {
                Double_t MetCut = b*1.0;
                //Inside Z Window
                if (fabs(dilepton.M()-91.1876) < 15.0) {
                  if (met.Pt() > MetCut) {
                    CountInZWindow_AfterTCMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterTCMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfMet.Pt() > MetCut) {
                    CountInZWindow_AfterPFMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfTrackMet.Pt() > MetCut) {
                    CountInZWindow_AfterPFTrackMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFTrackMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfNoFwdMet.Pt() > MetCut) {
                    CountInZWindow_AfterPFNoFwdMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFNoFwdMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFMet.Pt() > MetCut) {
                    CountInZWindow_AfterMinPFTrackMetPFMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFNoFwdMetPFMet.Pt() > MetCut) {
                    CountInZWindow_AfterMinPFNoFwdMetPFMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFNoFwdMetPFMet.Pt() > MetCut) {
                    CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b] += eventweight*eventweight;
                  }

                  if (METdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterTCMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFMETdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFTrackMETdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFNoFwdMETdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFMetdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFNoFwdMetPFMetdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt > MetCut) {
                    CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }

                 } else {
                  //outside Z window
                  if (met.Pt() > MetCut) {
                    CountOutOfZWindow_AfterTCMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterTCMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterPFMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfTrackMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterPFTrackMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFTrackMetCut[q][b] += eventweight*eventweight;
                  }
                  if (pfNoFwdMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterPFNoFwdMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFNoFwdMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterMinPFTrackMetPFMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFNoFwdMetPFMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFNoFwdMetPFMet.Pt() > MetCut) {
                    CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b] += eventweight*eventweight;
                  }

                  if (METdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterTCMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFMETdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFTrackMETdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (PFNoFwdMETdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFMetdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFNoFwdMetPFMetdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }
                  if (MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt > MetCut) {
                    CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight;
                    CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b] += eventweight*eventweight;
                  }

                }
              } //loop over met cuts
              
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
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  //--------------------------------------------------------------------------------------------------------------
  //Make R_in/out plots
  //============================================================================================================== 
  vector<TGraphAsymmErrors*> RInOutAfterTCMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFTrackMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFNoFwdMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFTrackMetPFMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFNoFwdMetPFMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut;
  vector<TGraphAsymmErrors*> RInOutAfterTCMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFTrackMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterPFNoFwdMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut;
  vector<TGraphAsymmErrors*> RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut;
  const Int_t nMetPoints = 100;
  for (int q = 0; q<processNames.size() ; ++q) { 
    double cutValue[nMetPoints];
    double cutValueErr[nMetPoints];
    double ratio[nMetPoints];
    double ratioErr[nMetPoints];


    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterTCMetCut[q][b] / CountInZWindow_AfterTCMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterTCMetCut[q][b])/CountOutOfZWindow_AfterTCMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterTCMetCut[q][b])/CountInZWindow_AfterTCMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterTCMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr,cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterTCMetCut->SetName(("RInOutAfterTCMetCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFMetCut[q][b] / CountInZWindow_AfterPFMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFMetCut[q][b])/CountOutOfZWindow_AfterPFMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFMetCut[q][b])/CountInZWindow_AfterPFMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr,cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterPFMetCut->SetName(("RInOutAfterPFMetCut_"+processNames[q]+label).c_str());
   
    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFTrackMetCut[q][b] / CountInZWindow_AfterPFTrackMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFTrackMetCut[q][b])/CountOutOfZWindow_AfterPFTrackMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFTrackMetCut[q][b])/CountInZWindow_AfterPFTrackMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFTrackMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr,ratioErr,ratioErr);
    tmpRInOutAfterPFTrackMetCut->SetName(("RInOutAfterPFTrackMetCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFNoFwdMetCut[q][b] / CountInZWindow_AfterPFNoFwdMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFNoFwdMetCut[q][b])/CountOutOfZWindow_AfterPFNoFwdMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFNoFwdMetCut[q][b])/CountInZWindow_AfterPFNoFwdMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFNoFwdMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterPFNoFwdMetCut->SetName(("RInOutAfterPFNoFwdMetCut_"+processNames[q]+label).c_str());
 
    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFTrackMetPFMetCut[q][b] / CountInZWindow_AfterMinPFTrackMetPFMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetCut[q][b])/CountOutOfZWindow_AfterMinPFTrackMetPFMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFTrackMetPFMetCut[q][b])/CountInZWindow_AfterMinPFTrackMetPFMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFTrackMetPFMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFTrackMetPFMetCut->SetName(("RInOutAfterMinPFTrackMetPFMetCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut[q][b] / CountInZWindow_AfterMinPFNoFwdMetPFMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetCut[q][b])/CountOutOfZWindow_AfterMinPFNoFwdMetPFMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetCut[q][b])/CountInZWindow_AfterMinPFNoFwdMetPFMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFNoFwdMetPFMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFNoFwdMetPFMetCut->SetName(("RInOutAfterMinPFNoFwdMetPFMetCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b] / CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b])/CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b])/CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut->SetName(("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut_"+processNames[q]+label).c_str());


    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterTCMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterTCMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterTCMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterTCMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterTCMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterTCMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterTCMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr,cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterTCMetdeltaPhilEtCut->SetName(("RInOutAfterTCMetdeltaPhilEtCut_"+processNames[q]+label).c_str());


    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterPFMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterPFMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterPFMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr,cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterPFMetdeltaPhilEtCut->SetName(("RInOutAfterPFMetdeltaPhilEtCut_"+processNames[q]+label).c_str());
   
    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterPFTrackMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFTrackMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr,ratioErr,ratioErr);
    tmpRInOutAfterPFTrackMetdeltaPhilEtCut->SetName(("RInOutAfterPFTrackMetdeltaPhilEtCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterPFNoFwdMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterPFNoFwdMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterPFNoFwdMetdeltaPhilEtCut->SetName(("RInOutAfterPFNoFwdMetdeltaPhilEtCut_"+processNames[q]+label).c_str());
 
    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterMinPFTrackMetPFMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut->SetName(("RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut->SetName(("RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut_"+processNames[q]+label).c_str());

    for (UInt_t b = 0 ; b < 100; ++b) {     
      cutValue[b] = b*1.0;
      cutValueErr[b] = 0;        
      ratio[b] = CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b] / CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b];
      ratioErr[b] = ratio[b] * TMath::Sqrt(pow(TMath::Sqrt(CountErrSqrOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b])/CountOutOfZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b],2) + pow(TMath::Sqrt(CountErrSqrInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b])/CountInZWindow_AfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q][b],2));
    }
    TGraphAsymmErrors* tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut = new TGraphAsymmErrors (nMetPoints, cutValue, ratio, cutValueErr, cutValueErr, ratioErr, ratioErr);
    tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut->SetName(("RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut_"+processNames[q]+label).c_str());

    
    RInOutAfterTCMetCut.push_back(tmpRInOutAfterTCMetCut);
    RInOutAfterPFMetCut.push_back(tmpRInOutAfterPFMetCut);
    RInOutAfterPFTrackMetCut.push_back(tmpRInOutAfterPFTrackMetCut);
    RInOutAfterPFNoFwdMetCut.push_back(tmpRInOutAfterPFNoFwdMetCut);
    RInOutAfterMinPFTrackMetPFMetCut.push_back(tmpRInOutAfterMinPFTrackMetPFMetCut);
    RInOutAfterMinPFNoFwdMetPFMetCut.push_back(tmpRInOutAfterMinPFNoFwdMetPFMetCut);
    RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut.push_back(tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut);
    RInOutAfterTCMetdeltaPhilEtCut.push_back(tmpRInOutAfterTCMetdeltaPhilEtCut);
    RInOutAfterPFMetdeltaPhilEtCut.push_back(tmpRInOutAfterPFMetdeltaPhilEtCut);
    RInOutAfterPFTrackMetdeltaPhilEtCut.push_back(tmpRInOutAfterPFTrackMetdeltaPhilEtCut);
    RInOutAfterPFNoFwdMetdeltaPhilEtCut.push_back(tmpRInOutAfterPFNoFwdMetdeltaPhilEtCut);
    RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut.push_back(tmpRInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut);
    RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpRInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut);
    RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut.push_back(tmpRInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut);


  }

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwDYBkgSystematicsPlots.root", "UPDATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 

    file->WriteTObject(RInOutAfterTCMetCut[q],RInOutAfterTCMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFMetCut[q],RInOutAfterPFMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFTrackMetCut[q],RInOutAfterPFTrackMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFNoFwdMetCut[q],RInOutAfterPFNoFwdMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFTrackMetPFMetCut[q],RInOutAfterMinPFTrackMetPFMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFNoFwdMetPFMetCut[q],RInOutAfterMinPFNoFwdMetPFMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut[q],RInOutAfterMinPFTrackMetPFNoFwdMetPFMetCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterTCMetdeltaPhilEtCut[q],RInOutAfterTCMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFMetdeltaPhilEtCut[q],RInOutAfterPFMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFTrackMetdeltaPhilEtCut[q],RInOutAfterPFTrackMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterPFNoFwdMetdeltaPhilEtCut[q],RInOutAfterPFNoFwdMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut[q],RInOutAfterMinPFTrackMetPFMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q],RInOutAfterMinPFNoFwdMetPFMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    file->WriteTObject(RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q],RInOutAfterMinPFTrackMetPFNoFwdMetPFMetdeltaPhilEtCut[q]->GetName(), "WriteDelete");
    
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
             && ele->HoverE < 0.025
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
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
 
  if (ele->pt < 20) {
    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 )) pass = kFALSE;
      }
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

        && (mu->nSeg > 1 || mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  if (mu->pt < 20) {
    if (!
        ( (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10
          && fabs(mu->d0) < 0.01
          )
      ) pass = kFALSE;    
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
