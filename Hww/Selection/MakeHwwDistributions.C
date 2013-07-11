//root -l EWKAna/Hww/Selection/MakeHwwDistributions.C+\(130,\"\",0\)


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

void MakeHwwDistributions(Double_t mHiggs, const string Label, Int_t ElectronSelectionType) 
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
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/HwwAnalysisSkimmed_full-d22_noskim.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-h130ww2l-gf-z2-v8-pu11_noskim_normalized.root");

//       inputFiles.push_back(vector<string>());
//       inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-v8-pu11_noskim_normalized.root");
//     inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ggww-z2-v8-pu11_noskim_normalized.root");
      inputFiles.push_back(vector<string>());
      inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-z2-v8-pu_noskim_normalized.root");

//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wz-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zz-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-tt-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-stop-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ttop-mg-z2-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wtop-mg-z2-v8-pu11_noskim_normalized.root");

//        inputFiles.push_back(vector<string>());
//        inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm-powheg-c10-v8-pu11_noskim_normalized.root");
       inputFiles.push_back(vector<string>());
       inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm-powheg-c10-v8-pu_noskim_normalized.root");

//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm1020-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt1020-powheg-c10-v8-pu11_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjetsl-z2-v8-pu11-skimmed_TwoRecoLeptonSkim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-wjets-mg-z2-v12-pu-skimmed_TwoRecoLeptonSkim_normalized.root");

  vector<string> processNames;
//   processNames.push_back("Data");
//    processNames.push_back("Hww130NoPU");
//   processNames.push_back("Hww130WithPU");
   processNames.push_back("WW"); 
//   processNames.push_back("WZ/ZZ");
//   processNames.push_back("top");
   processNames.push_back("Z");
//   processNames.push_back("WJets");

  vector<Int_t> processColors;
  processColors.push_back(kRed);
  processColors.push_back(kBlue);


  assert(processNames.size() == inputFiles.size());
  assert(processColors.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector <TH1F*>  fHWWSelection; 
  vector <TH1F*>  fHWWToEESelection; 
  vector <TH1F*>  fHWWToMuMuSelection; 
  vector <TH1F*>  fHWWToEMuSelection; 
  vector <TH1F*>  fDeltaPhiLeptons; 
  vector <TH1F*>  fDileptonMass; 
  vector <TH1F*>  fMet; 
  vector <TH1F*>  fPFMet; 
  vector <TH1F*>  fPFTrackMet; 
  vector <TH1F*>  fPFNoFwdMet; 
  vector <TH1F*>  fMinPFTrackMetPFMet; 
  vector <TH1F*>  fMinPFNoFwdMetPFMet; 
  vector <TH1F*>  fMinPFTrackMetPFNoFwdMetPFMet; 
  vector <TH1F*>  fMetDeltaPhilEt; 
  vector <TH1F*>  fPFMetDeltaPhilEt; 
  vector <TH1F*>  fPFTrackMetDeltaPhilEt; 
  vector <TH1F*>  fPFNoFwdMetDeltaPhilEt; 
  vector <TH1F*>  fMinPFTrackMetPFMetDeltaPhilEt; 
  vector <TH1F*>  fMinPFNoFwdMetPFMetDeltaPhilEt; 
  vector <TH1F*>  fMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt; 

  vector <TH2F*>  fPFTrackMetVsPFMet;


  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHWWSelection= new TH1F(("hHWWSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEESelection= new TH1F(("hHWWToEESelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToMuMuSelection= new TH1F(("hHWWToMuMuSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEMuSelection= new TH1F(("hHWWToEMuSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpDeltaPhiLeptons= new TH1F(("hDeltaPhiLeptons"+string("_")+processNames[q]+label).c_str(), ";DeltaPhiLeptons;Number of Events", 25, 0, 180);
    TH1F *tmpDileptonMass = new TH1F(("hDileptonMass"+string("_")+processNames[q]+label).c_str(), ";DileptonMass;Number of Events", 25, 0, 180);
    TH1F *tmpMet = new TH1F(("hMet"+string("_")+processNames[q]+label).c_str(), ";tcMet;Number of Events", 1000, 0, 100);
    TH1F *tmpPFMet = new TH1F(("hPFMet"+string("_")+processNames[q]+label).c_str(), ";PFMet;Number of Events", 1000, 0, 100);
    TH1F *tmpPFTrackMet = new TH1F(("hPFTrackMet"+string("_")+processNames[q]+label).c_str(), ";PFTrackMet;Number of Events", 1000, 0, 100);
    TH1F *tmpPFNoFwdMet = new TH1F(("hPFNoFwdMet"+string("_")+processNames[q]+label).c_str(), ";PFNoFwdMet;Number of Events", 1000, 0, 100);
    TH1F *tmpMinPFTrackMetPFMet = new TH1F(("hMinPFTrackMetPFMet"+string("_")+processNames[q]+label).c_str(), ";Min(PFTrackMet,PFMet);Number of Events", 1000, 0, 100);
     TH1F *tmpMinPFNoFwdMetPFMet = new TH1F(("hMinPFNoFwdMetPFMet"+string("_")+processNames[q]+label).c_str(), ";Min(PFNoFwdMet,PFMet);Number of Events", 1000, 0, 100);
     TH1F *tmpMinPFTrackMetPFNoFwdMetPFMet = new TH1F(("hMinPFTrackMetPFNoFwdMetPFMet"+string("_")+processNames[q]+label).c_str(), ";Min(PFTrackMet,PFNoFwdMet,PFMet);Number of Events", 1000, 0, 100);
    TH1F *tmpMetDeltaPhilEt = new TH1F(("hMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";tcMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpPFMetDeltaPhilEt = new TH1F(("hPFMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";PFMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpPFTrackMetDeltaPhilEt = new TH1F(("hPFTrackMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";PFTrackMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpPFNoFwdMetDeltaPhilEt = new TH1F(("hPFNoFwdMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";PFNoFwdMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpMinPFTrackMetPFMetDeltaPhilEt = new TH1F(("hMinPFTrackMetPFMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";MinPFTrackMetPFMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpMinPFNoFwdMetPFMetDeltaPhilEt = new TH1F(("hMinPFNoFwdMetPFMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";MinPFNoFwdMetPFMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    TH1F *tmpMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt = new TH1F(("hMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";MinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt;Number of Events", 1000, 0, 100);
    fHWWSelection.push_back(tmpHWWSelection);
    fHWWToEESelection.push_back(tmpHWWToEESelection);
    fHWWToMuMuSelection.push_back(tmpHWWToMuMuSelection);
    fHWWToEMuSelection.push_back(tmpHWWToEMuSelection);
    fDeltaPhiLeptons.push_back(tmpDeltaPhiLeptons);
    fDileptonMass.push_back(tmpDileptonMass);
    fMet.push_back(tmpMet);
    fPFMet.push_back(tmpPFMet);
    fPFTrackMet.push_back(tmpPFTrackMet);
    fPFNoFwdMet.push_back(tmpPFNoFwdMet);
    fMinPFTrackMetPFMet.push_back(tmpMinPFTrackMetPFMet);
    fMinPFNoFwdMetPFMet.push_back(tmpMinPFNoFwdMetPFMet);
    fMinPFTrackMetPFNoFwdMetPFMet.push_back(tmpMinPFTrackMetPFNoFwdMetPFMet);
    fMetDeltaPhilEt.push_back(tmpMetDeltaPhilEt);
    fPFMetDeltaPhilEt.push_back(tmpPFMetDeltaPhilEt);
    fPFTrackMetDeltaPhilEt.push_back(tmpPFTrackMetDeltaPhilEt);
    fPFNoFwdMetDeltaPhilEt.push_back(tmpPFNoFwdMetDeltaPhilEt);
    fMinPFTrackMetPFMetDeltaPhilEt.push_back(tmpMinPFTrackMetPFMetDeltaPhilEt);
    fMinPFNoFwdMetPFMetDeltaPhilEt.push_back(tmpMinPFNoFwdMetPFMetDeltaPhilEt);
    fMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt.push_back(tmpMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt);

    TH2F *tmpPFTrackMetVsPFMet = new TH2F(("hPFTrackMetVsPFMet"+string("_")+processNames[q]+label).c_str(), ";PFTrackMet;PFMet;Number of Events", 200, 0, 100, 200, 0, 100);
    fPFTrackMetVsPFMet.push_back(tmpPFTrackMetVsPFMet);

  }


  vector<Double_t> Count_WWTotal;
  vector<Double_t> Count_Total;
  vector<Double_t> Count_EtaBin0;
  vector<Double_t> Count_EtaBin1;
  vector<Double_t> Count_EtaBin2;
  vector<Double_t> Count_Total_statError;
  vector<Double_t> Count_EtaBin0_statError;
  vector<Double_t> Count_EtaBin1_statError;
  vector<Double_t> Count_EtaBin2_statError;
  for (int q=0; q<processNames.size() ; ++q) {
    Count_Total.push_back(0);
    Count_EtaBin0.push_back(0);
    Count_EtaBin1.push_back(0);
    Count_EtaBin2.push_back(0);
    Count_Total_statError.push_back(0);
    Count_EtaBin0_statError.push_back(0);
    Count_EtaBin1_statError.push_back(0);
    Count_EtaBin2_statError.push_back(0);  
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

//           if (!(info->nPUEvents < 2)) continue;
//         if (!(info->nPUEvents >= 5 && info->nPUEvents <= 7)) continue;
//         if (!(info->nPUEvents > 15)) continue;



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
            double MinPFTrackMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
            if(minDeltaPhiMinPFTrackMetPFMetLepton < TMath::Pi()/2.)
              MinPFTrackMetPFMetdeltaPhilEt = MinPFTrackMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFTrackMetPFMetLepton);


            double deltaPhiMinPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                          mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFNoFwdMetPFMetLepton = (deltaPhiMinPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFNoFwdMetPFMetLepton[1])?
              deltaPhiMinPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFNoFwdMetPFMetLepton[1];
            double MinPFNoFwdMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
            if(minDeltaPhiMinPFNoFwdMetPFMetLepton < TMath::Pi()/2.)
              MinPFNoFwdMetPFMetdeltaPhilEt = MinPFNoFwdMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFNoFwdMetPFMetLepton);


            double deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                                    mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton = (deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1])?
              deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1];
            double MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
            if(minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton < TMath::Pi()/2.)
              MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton);


                if (!(finalState == 1)) continue;
//               if (!(finalState == 0 || finalState == 3)) continue;
//              if (!(finalState == 1 || finalState == 2)) continue;
      //      if (!(lepton2.Pt() < 0 && lepton2.Pt() >= 0) continue;
               if (!(lepton2.Pt() >= 20)) continue;
//                   if (!(lepton2.Pt() < 20 && lepton2.Pt() >= 10)) continue;
//              if (!(lepton2.Pt() < 15 && lepton2.Pt() >= 10)) continue;
//             if (!(lepton2.Pt() < 10 && lepton2.Pt() >= 5)) continue;
   
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 14;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                                   false, false, false, false};
            assert(CutLabel.size() == nCuts+1);

            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 10.0) passCut[0] = true;

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

//             if (!(passCut[12] && passCut[13])) continue;

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
            if (  (passCut[0] && passCut[1] && passCut[6])
              ) {
              fMet[q]->Fill(TMath::Min(met.Pt(),99.9), eventweight);  
              fPFMet[q]->Fill(TMath::Min(pfMet.Pt(),99.9), eventweight);  
              fPFTrackMet[q]->Fill(TMath::Min(pfTrackMet.Pt(),99.9), eventweight);  
              fPFNoFwdMet[q]->Fill(TMath::Min(pfNoFwdMet.Pt(),99.9), eventweight);  
              fMinPFTrackMetPFMet[q]->Fill(TMath::Min(TMath::Min(pfTrackMet.Pt(),pfMet.Pt()),99.9), eventweight);  
              fMinPFNoFwdMetPFMet[q]->Fill(TMath::Min(TMath::Min(pfNoFwdMet.Pt(),pfMet.Pt()),99.9), eventweight);  
              fMinPFTrackMetPFNoFwdMetPFMet[q]->Fill(TMath::Min(TMath::Min(TMath::Min(pfTrackMet.Pt(),pfMet.Pt()),pfNoFwdMet.Pt()),99.9), eventweight);  
              fMetDeltaPhilEt[q]->Fill(TMath::Min(METdeltaPhilEt,99.9), eventweight);
              fPFMetDeltaPhilEt[q]->Fill(TMath::Min(PFMETdeltaPhilEt,99.9), eventweight);
              fPFTrackMetDeltaPhilEt[q]->Fill(TMath::Min(PFTrackMETdeltaPhilEt,99.9), eventweight);
              fPFNoFwdMetDeltaPhilEt[q]->Fill(TMath::Min(PFNoFwdMETdeltaPhilEt,99.9), eventweight);              
              fMinPFTrackMetPFMetDeltaPhilEt[q]->Fill(TMath::Min(MinPFTrackMetPFMetdeltaPhilEt,99.9), eventweight);
              fMinPFNoFwdMetPFMetDeltaPhilEt[q]->Fill(TMath::Min(MinPFNoFwdMetPFMetdeltaPhilEt,99.9), eventweight);
              fMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt[q]->Fill(TMath::Min(MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt,99.9), eventweight);

              fPFTrackMetVsPFMet[q]->Fill(TMath::Min(pfTrackMet.Pt(),99.9),TMath::Min(pfMet.Pt(),99.9),eventweight);

             }


 
            if ( (0==0) 
//                  && passCut[6]

                 && (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10] && passCut[11])
//                 && passCut[12]
              ) {
              fDileptonMass[q]->Fill(dilepton.M());              
//               if (passCut[12]) {
                fDeltaPhiLeptons[q]->Fill(deltaPhiLeptons);              
//               }
            }


            if (passAllCuts 
                //  && lepton2.Pt() < 20.0 
              ) {     
               Count_Total[q] += eventweight;
                Count_Total_statError[q] += eventweight*eventweight;
              if (fabs(lepton2.Eta()) >= 0.0 && fabs(lepton2.Eta()) < 1.0) {
                Count_EtaBin0[q] += eventweight;
                Count_EtaBin0_statError[q] += eventweight*eventweight;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_EtaBin1[q] += eventweight;
                Count_EtaBin1_statError[q] += eventweight*eventweight;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_EtaBin2[q] += eventweight;
                 Count_EtaBin2_statError[q] += eventweight*eventweight;
             }
            }

            //Cut Selection Histograms
            fHWWSelection[q]->Fill(-1, eventweight);
            if (finalState == 1 ) {
              fHWWToMuMuSelection[q]->Fill(-1, eventweight);
            } else if(finalState == 0 ) {
              fHWWToEESelection[q]->Fill(-1, eventweight);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelection[q]->Fill(-1, eventweight);
            }
    
            for (int k=0;k<nCuts;k++) {
              bool pass = true;
              bool passPreviousCut = true;
              for (int p=0;p<=k;p++) {
                pass = (pass && passCut[p]);
                if (p<k)
                  passPreviousCut = (passPreviousCut&& passCut[p]);
              }
      
              if (pass) {
                fHWWSelection[q]->Fill(k, eventweight);
                if (finalState == 1 ) {
                  fHWWToMuMuSelection[q]->Fill(k, eventweight);
                } else if(finalState == 0 ) {
                  fHWWToEESelection[q]->Fill(k, eventweight);
                } else if(finalState == 2 || finalState == 3 ) {
                  fHWWToEMuSelection[q]->Fill(k, eventweight);
                }
              }
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






  legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(fDileptonMass[q], processNames[q].c_str(), "LP");
  }

  //normalize distributions
  for (int q = 0; q<processNames.size() ; ++q) { 
    Double_t norm = 0;
    for (int b=0; b<fDileptonMass[q]->GetXaxis()->GetNbins()+2; ++b) { norm += fDileptonMass[q]->GetBinContent(b); }
    for (int b=0; b<fDileptonMass[q]->GetXaxis()->GetNbins()+2; ++b) {
      fDileptonMass[q]->SetBinContent(b,fDileptonMass[q]->GetBinContent(b) / norm);
      fDileptonMass[q]->SetBinError(b,fDileptonMass[q]->GetBinError(b) / norm);
    }
  }

  for (int q = 0; q<processNames.size() ; ++q) { 
    fDileptonMass[q]->SetLineColor(processColors[q]);
    fDileptonMass[q]->SetMarkerColor(processColors[q]);
    if (q==0) {
      fDileptonMass[q]->Draw();
    } else {
      fDileptonMass[q]->Draw("same");
    }
  }
  legend->Draw();
  cv->SaveAs("DileptonMass_PUStudy.gif");


  legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(fDeltaPhiLeptons[q], processNames[q].c_str(), "LP");
  }

  //normalize distributions
  for (int q = 0; q<processNames.size() ; ++q) { 
    Double_t norm = 0;
    for (int b=0; b<fDeltaPhiLeptons[q]->GetXaxis()->GetNbins()+2; ++b) { norm += fDeltaPhiLeptons[q]->GetBinContent(b); }
    for (int b=0; b<fDeltaPhiLeptons[q]->GetXaxis()->GetNbins()+2; ++b) {
      fDeltaPhiLeptons[q]->SetBinContent(b,fDeltaPhiLeptons[q]->GetBinContent(b) / norm);
      fDeltaPhiLeptons[q]->SetBinError(b,fDeltaPhiLeptons[q]->GetBinError(b) / norm);
    }
  }



//Make Met Plots
  const int nPoints = fMetDeltaPhilEt[0]->GetXaxis()->GetNbins();
  double cutValue[nPoints];
  double cutValueErr[nPoints];
  double NSig[nPoints];
  double NSigErr[nPoints];
  double NBkg[nPoints];
  double NBkgErr[nPoints];
  for(UInt_t b=0; b < nPoints; ++b) {
    cutValue[b] = fPFMetDeltaPhilEt[0]->GetXaxis()->GetBinCenter(b);
    cutValueErr[b] = 0;
    Double_t nsig = 0;
    Double_t nsigErrSqr = 0;
    Double_t nbkg = 0;
    Double_t nbkgErrSqr = 0;
    for (UInt_t q=b; q < nPoints; ++q) {
      nsig += fPFMetDeltaPhilEt[0]->GetBinContent(q);
      nsigErrSqr += pow(fPFMetDeltaPhilEt[0]->GetBinError(q),2);
      nbkg += fPFMetDeltaPhilEt[1]->GetBinContent(q);
      nbkgErrSqr += pow(fPFMetDeltaPhilEt[1]->GetBinError(q),2);
    }
    NSig[b] = nsig;
    NSigErr[b] = TMath::Sqrt(nsigErrSqr);
    NBkg[b] = nbkg;
    NBkgErr[b] = TMath::Sqrt(nbkgErrSqr);
  }
  TGraphAsymmErrors *NSigVsCut = new TGraphAsymmErrors (nPoints, cutValue, NSig, cutValueErr, cutValueErr, NSigErr, NSigErr);
  TGraphAsymmErrors *NSigVsNBkg = new TGraphAsymmErrors (nPoints, NBkg, NSig,NBkgErr, NBkgErr, NSigErr, NSigErr );
  
  NSigVsCut->Draw("AP");
  cv->SaveAs("NSigVsCut_AfterProjMetCut_PFMet.gif");

  NSigVsNBkg->Draw("AP");
  cv->SaveAs("NSigVsNBkg_AfterProjMetCut_PFMet.gif");



  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlots.root", "UPDATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
//     file->WriteTObject(fDeltaPhiLeptons[q], fDeltaPhiLeptons[q]->GetName(), "WriteDelete");
//     file->WriteTObject(fDileptonMass[q], fDileptonMass[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMet[q], fMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFMet[q],fPFMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFTrackMet[q],fPFTrackMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFNoFwdMet[q],fPFNoFwdMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFTrackMetPFMet[q],fMinPFTrackMetPFMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFNoFwdMetPFMet[q],fMinPFNoFwdMetPFMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFTrackMetPFNoFwdMetPFMet[q],fMinPFTrackMetPFNoFwdMetPFMet[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMetDeltaPhilEt[q],fMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFMetDeltaPhilEt[q], fPFMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFTrackMetDeltaPhilEt[q],fPFTrackMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPFNoFwdMetDeltaPhilEt[q],fPFNoFwdMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFTrackMetPFMetDeltaPhilEt[q],fMinPFTrackMetPFMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFNoFwdMetPFMetDeltaPhilEt[q],fMinPFNoFwdMetPFMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt[q],fMinPFTrackMetPFNoFwdMetPFMetDeltaPhilEt[q]->GetName(), "WriteDelete");
file->WriteTObject(fPFTrackMetVsPFMet[q],fPFTrackMetVsPFMet[q]->GetName(), "WriteDelete");

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
//             && ele->sigiEtaiEta < 0.01 
//             && fabs(ele->deltaEtaIn) < 0.004
//             && fabs(ele->deltaPhiIn) < 0.06
//             && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
//             && ele->nExpHitsInner <= 0
//             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//             && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
 //             && ele->sigiEtaiEta < 0.03
//              && fabs(ele->deltaEtaIn) < 0.007
//              && fabs(ele->deltaPhiIn) < 0.03
//              && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
//              && ele->nExpHitsInner <= 0
//              && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//              && fabs(ele->d0) < 0.02
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
