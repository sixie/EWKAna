//root -l EWKAna/Hww/FakeLeptonBkg/HwwFakeElectronPrediction.C+\(\"HwwNtuple_r11a-dmu-pr-v1_noskim_0000.root\",130,\"\"\) |& tee debugLog
//root -l EWKAna/Hww/FakeLeptonBkg/HwwFakeElectronPrediction.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full-pr_TightPlusRecoNoTriggerSkim.root\",130,\"\"\)


//root -l EWKAna/Hww/FakeLeptonBkg/HwwFakeElectronPrediction_1JetBin.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full-pr_TightPlusRecoTriggerSkim.root\",130,\"\"\)
//root -l EWKAna/Hww/FakeRate/HwwFakeElectronPrediction_1JetBin.C+\(\"WWAnalysisSkimmed_r10b-mu-d22_TwoRecoLeptonNoTriggerSkim.root\",130,\"\"\)
//root -l EWKAna/Hww/FakeRate/HwwFakeElectronPrediction_1JetBin.C+\(\"WWAnalysisSkimmed_full-d22_TwoRecoLeptonWithTriggerPlusOneTightLeptonSkim.root\",130,\"\"\)
//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
// Debug:
// cat diff | grep ">" | awk '{print "|| (info->runNum == " $2 " && info->evtNum == " $4 ")"}'
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
#include "MitPhysics/FakeMods/interface/FakeRate.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
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


Bool_t passConversionVeto(Int_t isConv) {
 
  Int_t tmp0 = floor(double(isConv) / 2.0);
  Int_t tmp1 = floor(double(tmp0) / 2.0);
  Int_t tmp2 = floor(double(tmp1) / 2.0);
  Int_t tmp3 = floor(double(tmp2) / 2.0);
  Int_t tmp4 = floor(double(tmp3) / 2.0);
  Int_t tmp5 = floor(double(tmp4) / 2.0);
  Int_t tmp6 = floor(double(tmp5) / 2.0);
  Int_t tmp7 = floor(double(tmp6) / 2.0);
  Int_t tmp8 = floor(double(tmp7) / 2.0);
  Int_t tmp9 = floor(double(tmp8) / 2.0);
  Int_t tmp10 = floor(double(tmp9) / 2.0);
  Int_t tmp11 = floor(double(tmp10) / 2.0);
  Int_t tmp12 = floor(double(tmp11) / 2.0);

  Bool_t pass;
  pass =  tmp9 % 2;
  return pass; 
}



//=== MAIN MACRO =================================================================================================

void HwwFakeElectronPrediction_1JetBin(const string inputFilename, Double_t mHiggs, 
                 const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  Double_t lumi;              // luminosity (pb^-1)
  Int_t ChargeSelection = 0;
//   ChargeSelection = 1;

  Double_t fPtMaxLowerCut;
  Double_t fPtMinLowerCut;
  Double_t fDileptonMassUpperCut;
  Double_t fDeltaPhiCut;
  Double_t fMtLowerCut;
  Double_t fMtUpperCut;
 
  Double_t fHiggsMass = mHiggs;

  //Configure Higgs Mass Dependant Cuts
  if (fHiggsMass == 120) { fPtMaxLowerCut = 20.0;        fPtMinLowerCut = 10.0; 
    fDileptonMassUpperCut = 40.0; fDeltaPhiCut = 60.0; fMtLowerCut = 70; fMtUpperCut = 120;
  }
  if (fHiggsMass == 130) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 10.0;
    fDileptonMassUpperCut = 45.0; fDeltaPhiCut = 60.0; fMtLowerCut = 75; fMtUpperCut = 125;
  }
  if (fHiggsMass == 140) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 20.0;
    fDileptonMassUpperCut = 45.0; fDeltaPhiCut = 60.0; fMtLowerCut = 80; fMtUpperCut = 130;
  }
  if (fHiggsMass == 150) { fPtMaxLowerCut = 27.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0; fMtLowerCut = 80; fMtUpperCut = 150;
  }
  if (fHiggsMass == 160) { fPtMaxLowerCut = 30.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0; fMtLowerCut = 90; fMtUpperCut = 160;
  }
  if (fHiggsMass == 170) { fPtMaxLowerCut = 34.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 50.0; fDeltaPhiCut = 60.0; fMtLowerCut = 20; fMtUpperCut = 170;
  }
  if (fHiggsMass == 180) { fPtMaxLowerCut = 36.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 60.0; fDeltaPhiCut = 70.0; fMtLowerCut = 20; fMtUpperCut = 180;
  }
  if (fHiggsMass == 190) { fPtMaxLowerCut = 38.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 80.0; fDeltaPhiCut = 90.0; fMtLowerCut = 20; fMtUpperCut = 190;
  }
  if (fHiggsMass == 200) { fPtMaxLowerCut = 40.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 90.0; fDeltaPhiCut = 100.0; fMtLowerCut = 20; fMtUpperCut = 200;
  }
  if (fHiggsMass == 210) { fPtMaxLowerCut = 44.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 110.0; fDeltaPhiCut = 110.0; fMtLowerCut = 20; fMtUpperCut = 210;
  }
  if (fHiggsMass == 220) { fPtMaxLowerCut = 48.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 120.0; fDeltaPhiCut = 120.0; fMtLowerCut = 20; fMtUpperCut = 220;
  }
  if (fHiggsMass == 230) { fPtMaxLowerCut = 52.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 130.0; fDeltaPhiCut = 130.0; fMtLowerCut = 20; fMtUpperCut = 230; 
  }
  if (fHiggsMass == 250) { fPtMaxLowerCut = 55.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 150.0; fDeltaPhiCut = 140.0; fMtLowerCut = 20; fMtUpperCut = 250; 
  }
  if (fHiggsMass == 300) { fPtMaxLowerCut = 70.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 200.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 300;
  }
  if (fHiggsMass == 350) { fPtMaxLowerCut = 80.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 250.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 350;
  } 
  if (fHiggsMass == 400) { fPtMaxLowerCut = 90.0;         fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 300.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 400;
  }
  if (fHiggsMass == 450) { fPtMaxLowerCut = 110.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 350.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 450;
  }
  if (fHiggsMass == 500) { fPtMaxLowerCut = 120.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 400.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 500;
  }
  if (fHiggsMass == 550) { fPtMaxLowerCut = 130.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 450.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 550;
  }
  if (fHiggsMass == 600) { fPtMaxLowerCut = 140.0;        fPtMinLowerCut = 25.0;
    fDileptonMassUpperCut = 500.0; fDeltaPhiCut = 175.0; fMtLowerCut = 20; fMtUpperCut = 600;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Set up Fake Rate
  //==============================================================================================================
  Bool_t use2DFakeRate = kTRUE;
  Bool_t useFitFunction = kFALSE;
  mithep::FakeRate *fFakeRate = new mithep::FakeRate(
                                                         "ElectronFakeRate.SmurfV6.skim.root",
                                                         "MuonFakeRate.SmurfV6.skim.root",
                                                         "", "",
                                                         "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta",
                                                         "MuonFakeRateDenominatorM2_Mu8PtCombinedSample_ptThreshold15_PtEta",
                                                         use2DFakeRate, useFitFunction );




  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *fLeptonJetPt = new TH1D("histLeptonJetPt_TightPlusFailSample", ";Cut Number;Number of Events", 100, 0,100); 



  TH1D *fHWWSelection= new TH1D("hHWWSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection= new TH1D("hHWWToEESelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection= new TH1D("hHWWToMuMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection= new TH1D("hHWWToEMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWSelection_sysError= new TH1D("hHWWSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection_sysError= new TH1D("hHWWToEESelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection_sysError= new TH1D("hHWWToMuMuSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection_sysError= new TH1D("hHWWToEMuSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);

  TH1D *fLeptonFakePt = new TH1D(       "hLeptonFakePt",";Lepton Fake P_t Min;Number of Events",60,0.,150.);
  TH1D *fLeptonFakeEta = new TH1D(      "hLeptonFakeEta",";Lepton Fake #eta;Number of Events",30,-5.0,5.0);
  TH2D *fLeptonFakePtEta = new TH2D ( "hLeptonFakePtEta", ";Lepton Fake p_{T}; Lepton Fake #eta; Number of Events", 20,0.,100., 30,-5.0,5.0);

  TH1D *fLeptonEta = new TH1D(          "hLeptonEta",";LeptonEta;Number of Events",30,-5.0,5.0);
  TH1D *fLeptonPtMax = new TH1D(        "hLeptonPtMax",";Lepton P_t Max;Number of Events",60,0.,150.);
  TH1D *fLeptonPtMin = new TH1D(        "hLeptonPtMin",";Lepton P_t Min;Number of Events",60,0.,150.);
  TH1D *fMetPtHist = new TH1D(          "hMetPtHist",";Met;Number of Events",30,0.,300.);  
  TH1D *fMetPhiHist = new TH1D(         "hMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  TH1D *fUncorrMetPtHist = new TH1D(    "hUncorrMetPtHist",";Met;Number of Events",30,0.,300.);  
  TH1D *fUncorrMetPhiHist = new TH1D(   "hUncorrMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  TH1D *fDeltaPhiLeptons = new TH1D(    "hDeltaPhiLeptons",";#Delta#phi_{ll};Number of Events",18,0,180);
  TH1D *fDeltaEtaLeptons = new TH1D(    "hDeltaEtaLeptons",";#Delta#eta_{ll};Number of Events",20,-5.0,5.0);

  TH1D *fLeptonEta_sysError = new TH1D( "hLeptonEta_sysError",";LeptonEta;Number of Events",20,-5.0,5.0);
  TH1D *fLeptonPtMax_sysError = new TH1D( "hLeptonPtMax_sysError",";Lepton P_t Max;Number of Events",20,0.,150.);
  TH1D *fLeptonPtMin_sysError = new TH1D( "hLeptonPtMin_sysError",";Lepton P_t Min;Number of Events",20,0.,150.);
  TH1D *fMetPtHist_sysError = new TH1D(  "hMetPtHist_sysError",";Met;Number of Events",30,0.,300.);  
  TH1D *fDeltaPhiLeptons_sysError = new TH1D( "hDeltaPhiLeptons_sysError",";#Delta#phi_{ll};Number of Events",18,0,180);



  TH1D *fMinDeltaPhiLeptonMet_afterCuts = new TH1D(    "hMinDeltaPhiLeptonMet_afterCuts", 
                                            ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  TH1D *fMtLepton1_afterCuts = new TH1D(               "hMtLepton1_afterCuts",
                                             ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  TH1D *fMtLepton2_afterCuts = new TH1D(               "hMtLepton2_afterCuts",
                                             ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  TH1D *fMtHiggs_afterCuts = new TH1D(                 "hMtHiggs_afterCuts",
                                             ";M_t (l1+l2+Met);Number of Events",150,0.,300.);
  TH1D *fLeptonPtPlusMet_afterCuts = new TH1D(         "hLeptonPtPlusMet_afterCuts",
                                             ";LeptonPtPlusMet;Number of Events",150,0., 300.);


  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_ee = new TH1D("dileptonMass_ee", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_emu = new TH1D("dileptonMass_emu", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_mumu = new TH1D("dileptonMass_mumu", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_sysError = new TH1D("dileptonMass_sysError", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_ee_sysError = new TH1D("dileptonMass_ee_sysError", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_emu_sysError = new TH1D("dileptonMass_emu_sysError", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_mumu_sysError = new TH1D("dileptonMass_mumu_sysError", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  Double_t Count_ee = 0;
  Double_t Count_mm = 0;
  Double_t Count_em = 0;
  Double_t Count_me = 0;
  Double_t Count_ee_statError = 0;
  Double_t Count_mm_statError = 0;
  Double_t Count_em_statError = 0;
  Double_t Count_me_statError = 0;
  Double_t Count_ee_sysError = 0;
  Double_t Count_mm_sysError = 0;
  Double_t Count_em_sysError = 0;
  Double_t Count_me_sysError = 0;

  Double_t Count_Pt10To15_EtaBin0 = 0;
  Double_t Count_Pt10To15_EtaBin1 = 0;
  Double_t Count_Pt10To15_EtaBin2 = 0;
  Double_t Count_Pt10To15_EtaBin0_statError = 0;
  Double_t Count_Pt10To15_EtaBin1_statError = 0;
  Double_t Count_Pt10To15_EtaBin2_statError = 0;
  Double_t Count_Pt10To15_EtaBin0_sysError = 0;
  Double_t Count_Pt10To15_EtaBin1_sysError = 0;
  Double_t Count_Pt10To15_EtaBin2_sysError = 0;
  Double_t Count_Pt15To20_EtaBin0 = 0;
  Double_t Count_Pt15To20_EtaBin1 = 0;
  Double_t Count_Pt15To20_EtaBin2 = 0;
  Double_t Count_Pt15To20_EtaBin0_statError = 0;
  Double_t Count_Pt15To20_EtaBin1_statError = 0;
  Double_t Count_Pt15To20_EtaBin2_statError = 0;
  Double_t Count_Pt15To20_EtaBin0_sysError = 0;
  Double_t Count_Pt15To20_EtaBin1_sysError = 0;
  Double_t Count_Pt15To20_EtaBin2_sysError = 0;
  Double_t Count_Pt20To25_EtaBin0 = 0;
  Double_t Count_Pt20To25_EtaBin1 = 0;
  Double_t Count_Pt20To25_EtaBin2 = 0;
  Double_t Count_Pt20To25_EtaBin0_statError = 0;
  Double_t Count_Pt20To25_EtaBin1_statError = 0;
  Double_t Count_Pt20To25_EtaBin2_statError = 0;
  Double_t Count_Pt20To25_EtaBin0_sysError = 0;
  Double_t Count_Pt20To25_EtaBin1_sysError = 0;
  Double_t Count_Pt20To25_EtaBin2_sysError = 0;
  Double_t Count_Pt25To30_EtaBin0 = 0;
  Double_t Count_Pt25To30_EtaBin1 = 0;
  Double_t Count_Pt25To30_EtaBin2 = 0;
  Double_t Count_Pt25To30_EtaBin0_statError = 0;
  Double_t Count_Pt25To30_EtaBin1_statError = 0;
  Double_t Count_Pt25To30_EtaBin2_statError = 0;
  Double_t Count_Pt25To30_EtaBin0_sysError = 0;
  Double_t Count_Pt25To30_EtaBin1_sysError = 0;
  Double_t Count_Pt25To30_EtaBin2_sysError = 0;
  Double_t Count_Pt30ToInf_EtaBin0 = 0;
  Double_t Count_Pt30ToInf_EtaBin1 = 0;
  Double_t Count_Pt30ToInf_EtaBin2 = 0;
  Double_t Count_Pt30ToInf_EtaBin0_statError = 0;
  Double_t Count_Pt30ToInf_EtaBin1_statError = 0;
  Double_t Count_Pt30ToInf_EtaBin2_statError = 0;
  Double_t Count_Pt30ToInf_EtaBin0_sysError = 0;
  Double_t Count_Pt30ToInf_EtaBin1_sysError = 0;
  Double_t Count_Pt30ToInf_EtaBin2_sysError = 0;



  Double_t Count_eF_Pt10To20_Barrel = 0;
  Double_t Count_eF_Pt10To20_Endcap = 0;
  Double_t Count_eF_Pt20ToInf_Barrel = 0;
  Double_t Count_eF_Pt20ToInf_Endcap = 0;
  Double_t Count_mF_Pt10To20_Barrel = 0;
  Double_t Count_mF_Pt10To20_Endcap = 0;
  Double_t Count_mF_Pt20ToInf_Barrel = 0;
  Double_t Count_mF_Pt20ToInf_Endcap = 0;
  Double_t Count_eF_Pt10To20_Barrel_statError = 0;
  Double_t Count_eF_Pt10To20_Endcap_statError = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_statError = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_statError = 0;
  Double_t Count_mF_Pt10To20_Barrel_statError = 0;
  Double_t Count_mF_Pt10To20_Endcap_statError = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_statError = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_statError = 0;
  Double_t Count_eF_Pt10To20_Barrel_sysError = 0;
  Double_t Count_eF_Pt10To20_Endcap_sysError = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_sysError = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_sysError = 0;
  Double_t Count_mF_Pt10To20_Barrel_sysError = 0;
  Double_t Count_mF_Pt10To20_Endcap_sysError = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_sysError = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_sysError = 0;




  ofstream eventListFile("TightPlusFailEventList.sixie.txt");

  Double_t NEventsPreSel = 0;
  Double_t NEventsPreSel_statError = 0;
  Double_t NEventsPreSel_sysError = 0;
  Double_t NEventsPreSelPassHiggsCuts = 0;
  Double_t NEventsPreSelPassHiggsCuts_statError = 0;
  Double_t NEventsPreSelPassHiggsCuts_sysError = 0;


  fHWWSelection->Sumw2();
  fHWWToEESelection->Sumw2();
  fHWWToMuMuSelection->Sumw2();
  fHWWToEMuSelection->Sumw2();
  fLeptonEta->Sumw2();
  fLeptonPtMax->Sumw2();
  fLeptonPtMin->Sumw2();
  fLeptonFakePt->Sumw2();
  fLeptonFakeEta->Sumw2();
  fMetPtHist->Sumw2();
  fMetPhiHist->Sumw2();
  fUncorrMetPtHist->Sumw2();
  fUncorrMetPhiHist->Sumw2();
  fDeltaPhiLeptons->Sumw2();
  fDeltaEtaLeptons->Sumw2();
  fMinDeltaPhiLeptonMet_afterCuts->Sumw2();
  fMtLepton1_afterCuts->Sumw2();
  fMtLepton2_afterCuts->Sumw2();
  fMtHiggs_afterCuts->Sumw2();
  fLeptonPtPlusMet_afterCuts->Sumw2();
  dileptonMass->Sumw2();
  dileptonMass_ee->Sumw2();
  dileptonMass_emu->Sumw2();
  dileptonMass_mumu->Sumw2();


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
  
  inputFile = new TFile(inputFilename.c_str());
  assert(inputFile);

  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//  rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
  rlrm.AddJSONFile("certifiedUCSD.json"); 
//   hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
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
  
  vector<UInt_t> runs;
  vector<UInt_t> events;

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    Bool_t printDebug = kFALSE;
    if ( (0 == 1) 
         || (info->runNum == 163340 && info->evtNum == 78249319)
         || (info->runNum == 163374 && info->evtNum == 340611883)
         || (info->runNum == 163589 && info->evtNum == 46972769)
//         || (info->runNum == 163630 && info->evtNum == 47444749)  //dz cut
//         || (info->runNum == 163664 && info->evtNum == 23060169)  //dz cut
         || (info->runNum == 163758 && info->evtNum == 175891598)   //fires doubleEle but not in doubleEle PD
//         || (info->runNum == 163758 && info->evtNum == 193977241) //dz cut
         || (info->runNum == 163759 && info->evtNum == 324734768)   //not in my ntuples, but should be.
         || (info->runNum == 163796 && info->evtNum == 157371044)  //fires doubleEle but not in doubleEle PD
         || (info->runNum == 163869 && info->evtNum == 67318151)

      ) {
      printDebug = kTRUE;
    }

    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
//     //for the skimmed input, I already required the HLT bits.
//     if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;


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
    // Met
    //********************************************************
    TVector3 pfMet;        
    TVector3 pfTrackMet;        
    TVector3 pfNoFwdMet;        
    if(info->pfMEx!=0 || info->pfMEy!=0) {       
      pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
    }
    if(info->pfTrackMEx!=0 || info->pfTrackMEx!=0) {       
      pfTrackMet.SetXYZ(info->pfTrackMEx, info->pfTrackMEy, 0);
    }
    if(info->pfNeutralNoFwdMEx+info->pfTrackMEx!=0 || info->pfNeutralNoFwdMEy+info->pfTrackMEy!=0) {       
      pfNoFwdMet.SetXYZ(info->pfNeutralNoFwdMEx+info->pfTrackMEx, 
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

    
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
      if ( (0==0)
           &&
           passMuonCuts(mu)
           &&
           fabs(mu->eta) < 2.4
           && 
           mu->pt > 10.0
        ) {
        leptonPt.push_back(mu->pt);
        leptonEta.push_back(mu->eta);
        leptonPhi.push_back(mu->phi);
        leptonType.push_back(13);
        leptonIndex.push_back(i);  
        leptonCharge.push_back(mu->q);
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

      if ( (0==0)
           && 
           passElectronCuts(ele)
           &&
           fabs(ele->eta) < 2.5
           && 
           ele->pt > 10.0
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

    


    
    //******************************************************************************
    //select events with 1lepton only
    //******************************************************************************
    if (leptonPt.size() < 1) continue;
    if (!(leptonPt[0] > 10.0)) continue;
    if (!(leptonPt.size() < 2 || leptonPt[1] < 10.0)) continue;
 

    //******************************************************************************
    //Check duplicate events
    //******************************************************************************
    Bool_t foundRunAndEvent = kFALSE;
    for (UInt_t q=0; q<runs.size(); ++q) {
      if (runs[q] == info->runNum && events[q] == info->evtNum) foundRunAndEvent = kTRUE;
    }
    if (!foundRunAndEvent) {
      runs.push_back(info->runNum);
      events.push_back(info->evtNum);
    } else {
      continue;
    }


    if (printDebug) {
      cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
      cout << "PV: " << info->pvx << " " << info->pvy << " " << info->pvz << endl;

      cout << "Triggers:\n";
      cout << Bool_t(info->triggerBits & kHLT_DoubleMu7) << " "  
           << Bool_t(info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) << " "  
           << Bool_t(info->triggerBits & kHLT_Mu17_Ele8_CaloIdL) << " "  
           << Bool_t(info->triggerBits & kHLT_Mu8_Ele17_CaloIdL) << " "  
           << Bool_t(info->triggerBits & kHLT_Mu15) << " "  
           << Bool_t(info->triggerBits & kHLT_IsoMu17) << " "  
           << Bool_t(info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) << " "  
           << Bool_t(info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) << " "  
           << Bool_t(info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) << " "  
           << endl;
      


      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        Double_t isoCutValue = 0;            
        if (fabs(ele->scEta) < 1.479) {
          if (ele->pt > 20) {
            isoCutValue = 0.15;
          } else {
            isoCutValue = 0.15;
          }
        } else {
          if (ele->pt > 20) {
            isoCutValue = 0.09;
          } else {
            isoCutValue = 0.09;
          }
        }
            
        cout << "Electron" << i << " : " << ele->pt << " " << ele->eta << " " << ele->phi << " "
             << ele->isEB << " " 
             << ele->q << " " 
             << ele->sigiEtaiEta << " " 
             << ele->deltaEtaIn << " " 
             << ele->deltaPhiIn << " " 
             << ele->HoverE << " " 
             << ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold << " "
//                  <<ele->trkIso03 << " " 
//                  << ele->emIso03 << " " 
//                  << ele->hadIso03 << " " 
             << ele->nExpHitsInner << " "
             << passConversionVeto(ele->isConv) << " "
             << ele->d0 << " "
             << ele->dz << " "
             << ele->fBrem << " "
             << ele->EOverP << " "
             << endl;   
          
        if (fabs(ele->scEta) < 1.479) {
          cout << Bool_t(ele->sigiEtaiEta < 0.01) << " " 
               << Bool_t(fabs(ele->deltaEtaIn) < 0.004 ) << " "
               << Bool_t(fabs(ele->deltaPhiIn) < 0.06 ) << " ";
        } else {
          cout << Bool_t(ele->sigiEtaiEta < 0.03) << " " 
               << Bool_t(fabs(ele->deltaEtaIn) < 0.007) << " "
               << Bool_t(fabs(ele->deltaPhiIn) < 0.03)<< " ";
        }        
        cout << Bool_t((ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold) / ele->pt < isoCutValue ) << " "
             << Bool_t(ele->nExpHitsInner <= 0 ) << " "
             << Bool_t(passConversionVeto(ele->isConv) ) << " "
             << Bool_t(fabs(ele->d0) < 0.02 ) << " "
             << Bool_t(fabs(ele->dz) < 0.1 ) << " "
             << Bool_t(ele->isEcalDriven ) << " : "            
             << passElectronCuts(ele) << " "
             << endl;

        
        if (fabs(ele->scEta) < 1.479) {
          cout << Bool_t(ele->sigiEtaiEta < 0.01) << " " 
               << Bool_t(fabs(ele->deltaEtaIn) < 0.007 ) << " "
               << Bool_t(fabs(ele->deltaPhiIn) < 0.15 ) << " "
               << Bool_t(ele->HoverE < 0.12) << " " ;
        } else {
          cout << Bool_t(ele->sigiEtaiEta < 0.03) << " " 
               << Bool_t(fabs(ele->deltaEtaIn) < 0.009) << " "
               << Bool_t(fabs(ele->deltaPhiIn) < 0.10)<< " "
               << Bool_t(ele->HoverE < 0.10) << " " ;
        }        
        cout << Bool_t(ele->nExpHitsInner <= 0 ) << " "
             << Bool_t(passConversionVeto(ele->isConv) ) << " "
             << Bool_t(fabs(ele->dz) < 0.1 ) << " "         
             << Bool_t((ele->trkIso03) / ele->pt < 0.2 ) << " "
             << Bool_t((ele->emIso03) / ele->pt < 0.2 ) << " "
             << Bool_t((ele->hadIso03) / ele->pt < 0.2 ) << " : "
             << passElectronDenominatorCuts(ele) << " "
             << endl;
          
      }
        
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
        Double_t isoCutValue = 0;            
        if (fabs(mu->eta) < 1.479) {
          if (mu->pt > 20) {
            isoCutValue = 0.22;
          } else {
            isoCutValue = 0.11;
          }
        } else {
          if (mu->pt > 20) {
            isoCutValue = 0.20;
          } else {
            isoCutValue = 0.11;
          }
        } 
            
        cout << "Muon" << i << " : " << mu->pt << " " << mu->eta << " " << mu->phi << " "
             << mu->q << " " 
             << Bool_t(mu->typeBits & kGlobal) << " " 
             << Bool_t(mu->typeBits & kTracker) << " " 
             << mu->nTkHits << " " 
             << mu->muNchi2 << " " 
             << Bool_t(mu->qualityBits & kGlobalMuonPromptTight) << " " 
             << Bool_t(mu->qualityBits & kTMLastStationAngTight) << " " 
             << mu->d0 << " " 
             << mu->dz << " " 
             << mu->ChargedIso04 + mu->NeutralIso04_10Threshold << " " 
//                  << mu->trkIso03 << " " 
//                  << mu->emIso03 << " " 
//                  << mu->hadIso03 << " " 
             << mu->nSeg << " " 
             << mu->nMatch << " " 
             << mu->nPixHits << " " 
             << mu->pterr / mu->pt << " " 
             << endl;   

        cout << Bool_t(mu->typeBits & kGlobal) << " " 
             << Bool_t(mu->typeBits & kTracker ) << " "
             << Bool_t(mu->nTkHits > 10 ) << " "
             << Bool_t(mu->muNchi2 < 10.0 ) << " "
             << Bool_t((mu->qualityBits & kGlobalMuonPromptTight)) << " "
             << Bool_t(fabs(mu->d0) < 0.02 ) << " "
             << Bool_t(fabs(mu->dz) < 0.2 ) << " "
             << Bool_t((mu->ChargedIso04 + mu->NeutralIso04_10Threshold ) / mu->pt < isoCutValue ) << " "
//                  << Bool_t((mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15 ) << " "
             << Bool_t(mu->nMatch > 1 ) << " "
             << Bool_t(mu->nPixHits > 0) << " "
             << Bool_t(mu->pterr / mu->pt < 0.1) << " "
             << " : " << passMuonCuts(mu) << " "
             << endl;
           
        Bool_t isCleanMuon = kFALSE;
        for (int k=0; k<leptonPt.size(); ++k) {
          if ( leptonType[k] == 13 
               && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[k],leptonEta[k]) < 0.1
            ) {
            isCleanMuon = kTRUE; 
            break;
          }              
        }
        Bool_t isSoft = kFALSE;
        if ( mu->pt > 3.0
             && (mu->typeBits & kTracker)
             && (mu->qualityBits & kTMLastStationAngTight)
             && mu->nTkHits > 10
             && fabs(mu->d0) < 0.2
             && fabs(mu->dz) < 0.1
             && !isCleanMuon
             && (!(mu->pt > 20 && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10))
          ) {
          isSoft = kTRUE;
        }
        cout  << Bool_t(mu->pt > 3.0) << " "
              << Bool_t(mu->typeBits & kTracker) << " "
              << Bool_t(mu->qualityBits & kTMLastStationAngTight) << " " 
              << Bool_t(mu->nTkHits > 10) << " " 
              << Bool_t(fabs(mu->d0) < 0.2) << " "
              << Bool_t(fabs(mu->dz) < 0.1) << " "
              << !isCleanMuon << " " 
              << Bool_t((!(mu->pt > 20 && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10))) << " : "            
              << isSoft << endl;


      }
    }


    //loop over electron denominators
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
   
      if (!(ele->pt > 10.0 && fabs(ele->eta) < 2.5)) continue;
      if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[0],leptonEta[0]) < 0.1 )
        continue;
      
      //select denominators that fail final selection
      if (!passElectronDenominatorCuts(ele)) continue;
      if (passElectronCuts(ele)) continue;


      //preselection
      if (!(ele->pt > 10.0 && fabs(ele->eta) < 2.5)) continue;
      if ((ChargeSelection == 0 && leptonCharge[0] == ele->q) || (ChargeSelection == 1 && leptonCharge[0] != ele->q)) continue;



      Double_t leptonJetPt = -1;
      for(Int_t jj=0; jj<jetArr->GetEntries(); jj++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[jj]);
        
        if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, jet->phi,jet->eta) < 0.5 ) {
          leptonJetPt = jet->pt;
        }                
      }
      if (leptonJetPt < 0) leptonJetPt = ele->pt;
      if (ele->pt > 10 && ele->pt < 20) {
        fLeptonJetPt->Fill(TMath::Min(TMath::Max(leptonJetPt, 0.01),99.9));
      }


//          if (! (ele->pt > 20)) continue;
//         if (! (ele->pt > 15 && ele->pt < 20)) continue;
//           if (! (ele->pt > 10 && ele->pt < 15)) continue;
//        if (! (fabs(ele->eta) > 0 && fabs(ele->eta) < 1.0)) continue;
//       if (! (fabs(ele->eta) > 1.0 && fabs(ele->eta) < 1.5)) continue;
//       if (! (fabs(ele->eta) > 1.5 && fabs(ele->eta) < 2.5)) continue;
      
     
      //******************************************************************************
      //Calculate Fake Rate
      //******************************************************************************
      double fakeRate = fFakeRate->ElectronFakeRate(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi) / (1-fFakeRate->ElectronFakeRate(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi));
      double fakeRateErrorLow = fFakeRate->ElectronFakeRateStatErrorLow(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi) / pow((1- fFakeRate->ElectronFakeRate(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi)),2);
      double fakeRateErrorHigh = fFakeRate->ElectronFakeRateStatErrorHigh(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi) / pow((1- fFakeRate->ElectronFakeRate(TMath::Min(double(ele->pt), double(34.5)), fabs(ele->eta), ele->phi)),2);

      Double_t eventWeight = fakeRate; 
//       eventWeight = 1.0; 
       Double_t eventWeightError = (  fakeRateErrorLow +   fakeRateErrorHigh) / 2;
       // eventWeightError = 0;
  
       //******************************************************************************
       //Correct the MET for the leptons
       //******************************************************************************
       if  (printDebug) cout << "Track MET Before Correction : " << info->pfTrackMEx  << " , " << info->pfTrackMEy << " : " << pfTrackMet.Pt() << endl;
       
       Double_t MetCorrection_X = 0;
       Double_t MetCorrection_Y = 0;

       if (leptonType[0] == 11) {
         const mithep::TElectron *e = (mithep::TElectron*)((*electronArr)[leptonIndex[0]]);           
         mithep::FourVectorM lepton;
         lepton.SetCoordinates(e->pt, e->eta, e->phi, 0.51099892e-3 );
         mithep::FourVectorM pfLepton;
         pfLepton.SetCoordinates(e->pfPt, e->pfEta, e->pfPhi, 0.51099892e-3 );
         if (e->pfPt >= 0) {
           MetCorrection_X += ( pfLepton.Px() - lepton.Px() );
           MetCorrection_Y += ( pfLepton.Py() - lepton.Py() );       
         } else {
           MetCorrection_X += ( 0.0 - lepton.Px() );
           MetCorrection_Y += ( 0.0 - lepton.Py() );       
         }
         if  (printDebug) cout << "Ele : " << pfLepton.Px() << " , " << pfLepton.Py() << " : " << lepton.Px() << " , " << lepton.Py() << endl;
       } else if (leptonType[0] == 13) {
         const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[0]]);
         mithep::FourVectorM lepton;
         lepton.SetCoordinates(mu->pt, mu->eta, mu->phi, 105.658369e-3 );
         mithep::FourVectorM pfLepton;
         pfLepton.SetCoordinates(mu->pfPt, mu->pfEta, mu->pfPhi, 105.658369e-3 );
         if (mu->pfPt >= 0) {
           MetCorrection_X += ( pfLepton.Px() - lepton.Px() );
           MetCorrection_Y += ( pfLepton.Py() - lepton.Py() );            
         } else {
           MetCorrection_X += ( 0.0 - lepton.Px() );
           MetCorrection_Y += ( 0.0 - lepton.Py() );   
         }
         if  (printDebug) cout << "Mu : " << pfLepton.Px() << " , " << pfLepton.Py() << " : " << lepton.Px() << " , " << lepton.Py() << endl;
       }
       mithep::FourVectorM fakeLepton;
       fakeLepton.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
       mithep::FourVectorM fakepfLepton;
       fakepfLepton.SetCoordinates(ele->pfPt, ele->pfEta, ele->pfPhi, 0.51099892e-3 );
       if (ele->pfPt >= 0) {
         MetCorrection_X += ( fakepfLepton.Px() - fakeLepton.Px() );
         MetCorrection_Y += ( fakepfLepton.Py() - fakeLepton.Py() );       
       } else {
         MetCorrection_X += ( 0.0 - fakeLepton.Px() );
         MetCorrection_Y += ( 0.0 - fakeLepton.Py() );       
       }
       if  (printDebug) cout << "Ele : " << fakepfLepton.Px() << " , " << fakepfLepton.Py() << " : " << fakeLepton.Px() << " , " << fakeLepton.Py() << endl;
       
       if  (printDebug)cout << "Met Correction : " << MetCorrection_X << " , " << MetCorrection_Y << endl;
       
       pfTrackMet.SetXYZ(info->pfTrackMEx + MetCorrection_X, info->pfTrackMEy + MetCorrection_Y , 0);
       
       if  (printDebug)cout << "Corrected Met : " << pfTrackMet.Px() << " , " << pfTrackMet.Py() << " : " << pfTrackMet.Pt() << endl;






      //******************************************************************************
      //Count Jets
      //******************************************************************************
       Int_t NJets = 0;
       const mithep::TJet *leadingJet = 0;       
       double maxBtag = -99999;

       for(Int_t j=0; j<jetArr->GetEntries(); j++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);
        if (jet->pt < 7.0) continue;

        Bool_t leptonOverlap = kFALSE;
        if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[0],leptonEta[0]) < 0.3 ||
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi,ele->eta) < 0.3 ) {
          leptonOverlap = kTRUE;
        }
        
        if (!leptonOverlap) {
          if (jet->pt > 30 && fabs(jet->eta) < 5.0 ) {
            if (!leadingJet || jet->pt > leadingJet->pt) {
              leadingJet = jet;
            }
            NJets++;
          }
          if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
        }
      }
      
     //******************************************************************************
      //soft muons
      //******************************************************************************
      NSoftMuons = 0;
      for(Int_t j=0; j<muonArr->GetEntries(); j++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[j]);
        Bool_t isCleanMuon = kFALSE;        
        if ( leptonType[0] == 13 
             && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[0],leptonEta[0]) < 0.1
          ) {
          isCleanMuon = kTRUE;         
        }
        Bool_t overlapWithElectronDenominator = kFALSE;
        if ( mithep::MathUtils::DeltaR(mu->phi, mu->eta, ele->phi, ele->eta) < 0.3
          ) {
          overlapWithElectronDenominator = kTRUE;         
        }
        if ( mu->pt > 3.0
             && (mu->typeBits & kTracker)
             && (mu->qualityBits & kTMLastStationAngTight)
             && mu->nTkHits > 10
             && fabs(mu->d0) < 0.2
             && fabs(mu->dz) < 0.1
             && !isCleanMuon
             && (!(mu->pt > 20 && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10))
             && !overlapWithElectronDenominator
          ) {
          NSoftMuons++;
        }
      }
      

      //***********************************************************************************************
      //|Z_vert-Z_l| maximum
      //***********************************************************************************************
      double zDiffMax = 0.0;
      
      double dzLepton;
      if (leptonType[0] == 11) {
        dzLepton = ((mithep::TElectron*)((*electronArr)[leptonIndex[0]]))->dz;
      } else {
        dzLepton = ((mithep::TMuon*)((*muonArr)[leptonIndex[0]]))->dz;
      }

      double dz = ele->dz;
//       if (dz > zDiffMax) zDiffMax = dz;
      zDiffMax = fabs(dz - dzLepton);



//       //******************************************************************************
//       //Make electron pt cuts
//       //******************************************************************************
//       //accept electrons only to 15
//       if (leptonType[0] == 11) {
//         if (leptonPt[0] < 15) continue;
//       }
//       if (ele->pt < 15) continue;


      //******************************************************************************
      //construct event variables & Final State
      //******************************************************************************
      Int_t finalState = -1;
      mithep::FourVectorM lepton1;
      mithep::FourVectorM lepton2;
      if (leptonPt[0] > ele->pt) {
        if (leptonType[0] == 11) {
          finalState = 0;
          lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
        } else {
          finalState = 3;
          lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
        }
        lepton2.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
      } else {
        if (leptonType[0] == 11) {
          finalState = 0;
          lepton2.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
        } else {
          finalState = 2;
          lepton2.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
        }
        lepton1.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
      }

      mithep::FourVectorM dilepton = lepton1+lepton2;

      double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(lepton1.Phi(), 
                                                           lepton2.Phi())* 180.0 / TMath::Pi();    




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
      
      double deltaPhiDileptonMet = mithep::MathUtils::DeltaPhi(dilepton.Phi(), pfMet.Phi());
      double mtHiggs = 2.0*dilepton.Pt()*pfMet.Pt()*(1.0 - cos(deltaPhiDileptonMet));
      

      //mtHiggs assuming neutrino pt == dilepton pt
      double enell = TMath::Sqrt(dilepton.Pt()*dilepton.Pt() + dilepton.M()*dilepton.M());
      double enenn = TMath::Sqrt(pfMet.Pt()* pfMet.Pt()  + dilepton.M()*dilepton.M());
      double enex  = dilepton.Px() + pfMet.Px();
      double eney  = dilepton.Py() + pfMet.Py();
      double mll   = dilepton.M();
      double mnu   = mll;
      double mtHiggs_NuPtEqualDileptonPt = mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney);
      


      //*********************************************************************************************
      //Define Cuts
      //*********************************************************************************************
      const int nCuts = 15;
      bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                             false, false, false, false, false};

      Bool_t PreselPtCut = kTRUE;
      if(!(lepton1.Pt() >  20.0 && lepton2.Pt() > 10.0)) PreselPtCut = kFALSE;
//       if (leptonPt[0] > ele->pt) {
//         if (!(ele->pt > 15.0)) PreselPtCut = kFALSE;
//       } else {
//         if (leptonType[0] == 11) {
//           if (!(leptonPt[0] > 15.0)) PreselPtCut = kFALSE;
//         }
//       }


      if(PreselPtCut) passCut[0] = true;
      
      if(zDiffMax < 100000.0) passCut[1] = true;            

      if(pfMet.Pt()    > 20.0)               passCut[2] = true; 
      
      if(dilepton.M() > 12.0)            passCut[3] = true;
   
      if (finalState == 0 || finalState == 1){ // mumu/ee
        if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
        if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 35
           && leadingJet && mithep::MathUtils::DeltaPhi(double(dilepton.Phi()), double(leadingJet->phi)) <= 165./180.*TMath::Pi()
          ) passCut[5] = true;
      }
      else if(finalState == 2 ||finalState == 3 ) { // emu
        passCut[4] = true;
        if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 20) passCut[5] = true;
      }
      
      if(NJets     == 1)              passCut[6] = true;
       
      if (NSoftMuons == 0 )      passCut[7] = true;
      
      passCut[8] = true;
      
      if(maxBtag < 2.1)          passCut[9] = true;
      if (lepton1.Pt() > fPtMaxLowerCut) passCut[10] = true;
      if (lepton2.Pt() > fPtMinLowerCut) passCut[11] = true;
      if (dilepton.M() < fDileptonMassUpperCut)   passCut[12] = true;
      if (deltaPhiLeptons < fDeltaPhiCut) passCut[13] = true;
      if (mtHiggs >= fMtLowerCut && mtHiggs <= fMtUpperCut ) passCut[14] = true;


      //*********************************************************************************************
      //Debug Printout
      //*********************************************************************************************  

      if (printDebug) {
        cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
        cout << "PV: " << info->pvx << " " << info->pvy << " " << info->pvz << endl;
        cout << leptonPt[0] << " " << leptonEta[0] << " " << leptonPhi[0] << " " << leptonType[0] << endl;
        cout << ele->pt << " " << ele->eta << " " << ele->phi << " " << "11" << endl;
        cout << zDiffMax << " " << pfMet.Pt() << " " << dilepton.M() << " " 
             << fabs(dilepton.M()-91.1876) << " " << PFMETdeltaPhilEt << "," <<  PFTrackMETdeltaPhilEt << " " 
             << NJets << " " << NSoftMuons << " " << leptonPt.size() 
             << " " << maxBtag << endl;
        cout << pfMet.Px() << " " << pfMet.Py() << "  : " << pfTrackMet.Px() << " " << pfTrackMet.Py() << endl;
      }

      if (printDebug) {
        for(Int_t q=0; q<jetArr->GetEntries(); q++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[q]);
                
          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (leptonType[k] == 11) {
              const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[leptonIndex[k]]);
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, tmpEle->scEta,tmpEle->scPhi) < 0.3) {
                leptonOverlap = kTRUE;
              }
              cout << "lepton " << k << " " << leptonPt[k] << " " << leptonEta[k] << " " << leptonPhi[k] << " , "
                   << tmpEle->scEta << " " << tmpEle->scPhi << " : DR = " 
                   << mithep::MathUtils::DeltaR(jet->phi, jet->eta, tmpEle->scPhi, tmpEle->scEta)  << endl;
            } else {
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
                leptonOverlap = kTRUE;
              }
              cout << "lepton " << k << " " << leptonPt[k] << " " << leptonEta[k] << " " << leptonPhi[k] << " : DR = "                          
                   << mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k])  << endl;
            }                                  
          }                
          cout << "Jet: " << jet->pt << " " << jet->eta << " " << jet->phi << " " << jet->TrackCountingHighEffBJetTagsDisc << " " << leptonOverlap << endl;
        }
         
        cout << "Event Selection: " << passCut[0] << " " << passCut[1] << " " << passCut[2] << " " << passCut[3] << " " << passCut[4] << " " << passCut[5] << " " << passCut[6] << " " << passCut[7] << " " << passCut[8] << " " << passCut[9] << " : ";
        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
          cout << "PassEvent " << info->runNum << " " << info->lumiSec << " " << info->evtNum  << "\n";
        }
        cout << endl;
      }



      //*********************************************************************************************
      //Make Selection Histograms. Number of events passing each level of cut
      //*********************************************************************************************  
      bool passAllCuts = true;
      for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
      //Cut Selection Histograms
      fHWWSelection->Fill(-1, eventWeight);
      fHWWSelection_sysError->Fill(-1, eventWeightError);
        if (finalState == 1 ) {
          fHWWToMuMuSelection->Fill(-1, eventWeight);
          fHWWToMuMuSelection_sysError->Fill(-1, eventWeightError);
        } else if(finalState == 0 ) {
          fHWWToEESelection->Fill(-1, eventWeight);
          fHWWToEESelection_sysError->Fill(-1, eventWeightError);
        } else { //if(finalState == 2 || finalState == 3 ) {
          fHWWToEMuSelection->Fill(-1, eventWeight);
          fHWWToEMuSelection_sysError->Fill(-1, eventWeightError);
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
            fHWWSelection->Fill(k, eventWeight);
            fHWWSelection_sysError->Fill(k, eventWeightError);
            if (finalState == 1 ) {
              fHWWToMuMuSelection->Fill(k, eventWeight);
              fHWWToMuMuSelection_sysError->Fill(k, eventWeightError);
            } else if(finalState == 0 ) {
              fHWWToEESelection->Fill(k, eventWeight);
              fHWWToEESelection_sysError->Fill(k, eventWeightError);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelection->Fill(k, eventWeight);
              fHWWToEMuSelection_sysError->Fill(k, eventWeightError);
            }
          }
        }

        //*****************************************************************************************
        //Make Preselection Histograms  
        //*****************************************************************************************
//         if (passCut[0] && passCut[1] ) {
          if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
//         if (passAllCuts) {

            NEventsPreSel += eventWeight;
            NEventsPreSel_statError += eventWeight*eventWeight;
            NEventsPreSel_sysError += eventWeightError;
            if (passCut[10] && passCut[11] && passCut[12] && passCut[13]) {
              NEventsPreSelPassHiggsCuts += eventWeight;
              NEventsPreSelPassHiggsCuts_statError += eventWeight*eventWeight;
              NEventsPreSelPassHiggsCuts_sysError += eventWeightError;              
            }

          fLeptonFakePt->Fill(ele->pt, eventWeight);
          fLeptonFakeEta->Fill(ele->eta, eventWeight);
          fLeptonFakePtEta->Fill(ele->pt, ele->eta, eventWeight);


//           fLeptonEta->Fill(lepton1.Eta(), eventWeight); 
          fLeptonEta->Fill(lepton2.Eta(), eventWeight);
          fLeptonPtMax->Fill(lepton1.Pt(), eventWeight);
          fLeptonPtMin->Fill(lepton2.Pt(), eventWeight);
          fMetPtHist->Fill(pfMet.Pt(), eventWeight);                             
          fMetPhiHist->Fill(pfMet.Phi(), eventWeight);                            
          fDeltaPhiLeptons->Fill(deltaPhiLeptons, eventWeight);

//           fLeptonEta_sysError->Fill(lepton2.Eta(), 2*eventWeightError);
          fLeptonEta_sysError->Fill(lepton2.Eta(), eventWeightError);
          fLeptonPtMin_sysError->Fill(lepton2.Pt(), eventWeightError);
          fLeptonPtMax_sysError->Fill(lepton1.Pt(), eventWeightError);
          fMetPtHist_sysError->Fill(pfMet.Pt(), eventWeightError);                             
          fDeltaPhiLeptons_sysError->Fill(deltaPhiLeptons, eventWeightError);
          dileptonMass->Fill(dilepton.M(), eventWeight);
          dileptonMass_sysError->Fill(dilepton.M(), eventWeightError);
          if (finalState == 0) {
            dileptonMass_ee->Fill(dilepton.M(), eventWeight);
            dileptonMass_ee_sysError->Fill(dilepton.M(), eventWeightError);
          } else if (finalState == 1) {
            dileptonMass_mumu->Fill(dilepton.M(), eventWeight);
            dileptonMass_mumu_sysError->Fill(dilepton.M(), eventWeightError);
          } else {
            dileptonMass_emu->Fill(dilepton.M(), eventWeight);
            dileptonMass_emu_sysError->Fill(dilepton.M(), eventWeightError);
          }
        }

        if (
//             passCut[0] && passCut[1] 
          passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//              passAllCuts           
           ) {

          if ((0==0) 
//              && finalState == 0
            ) {
            eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << ele->pt << " " << ele->eta << " : " << fakeRate << " - " << fakeRateErrorLow << " + " << fakeRateErrorHigh << endl;
          }

           //  if (passAllCuts ) {
          if (finalState == 0) {
            Count_ee += eventWeight;
            Count_ee_statError += eventWeight*eventWeight;
            Count_ee_sysError += eventWeightError;
          } else if (finalState == 1) {
            Count_mm += eventWeight;
            Count_mm_statError += eventWeight*eventWeight;
            Count_mm_sysError += eventWeightError;
          } else if (finalState == 2) {
            Count_em += eventWeight;
            Count_em_statError += eventWeight*eventWeight;
            Count_em_sysError += eventWeightError;
          } else {
            Count_me += eventWeight;
            Count_me_statError += eventWeight*eventWeight;
            Count_me_sysError += eventWeightError;
          }

          Double_t ptCategory = lepton2.Pt();
          ptCategory = ele->pt;

          if (ptCategory > 30.0) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt30ToInf_EtaBin0 += eventWeight;
              Count_Pt30ToInf_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt30ToInf_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt30ToInf_EtaBin1 += eventWeight;
              Count_Pt30ToInf_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt30ToInf_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt30ToInf_EtaBin2 += eventWeight;
              Count_Pt30ToInf_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt30ToInf_EtaBin2_sysError += eventWeightError;
            }
          } else if (ptCategory > 25.0 && ptCategory < 30.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt25To30_EtaBin0 += eventWeight;
              Count_Pt25To30_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt25To30_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt25To30_EtaBin1 += eventWeight;
              Count_Pt25To30_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt25To30_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt25To30_EtaBin2 += eventWeight;
              Count_Pt25To30_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt25To30_EtaBin2_sysError += eventWeightError;
            }
          } else if (ptCategory > 20.0 && ptCategory < 25.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt20To25_EtaBin0 += eventWeight;
              Count_Pt20To25_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt20To25_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt20To25_EtaBin1 += eventWeight;
              Count_Pt20To25_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt20To25_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt20To25_EtaBin2 += eventWeight;
              Count_Pt20To25_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt20To25_EtaBin2_sysError += eventWeightError;
            }
          } else if (ptCategory > 15.0 && ptCategory < 20.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt15To20_EtaBin0 += eventWeight;
              Count_Pt15To20_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt15To20_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt15To20_EtaBin1 += eventWeight;
              Count_Pt15To20_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt15To20_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt15To20_EtaBin2 += eventWeight;
              Count_Pt15To20_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt15To20_EtaBin2_sysError += eventWeightError;
            }
          } else if (ptCategory > 10.0 && ptCategory < 15.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt10To15_EtaBin0 += eventWeight;
              Count_Pt10To15_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt10To15_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt10To15_EtaBin1 += eventWeight;
              Count_Pt10To15_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt10To15_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt10To15_EtaBin2 += eventWeight;
              Count_Pt10To15_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt10To15_EtaBin2_sysError += eventWeightError;
            }
          } 


          if (leptonType[0] == 11) {
            if (ele->pt > 20) {
              if (fabs(ele->scEta) < 1.479) {
                Count_eF_Pt20ToInf_Barrel += eventWeight;
                Count_eF_Pt20ToInf_Barrel_statError += eventWeight*eventWeight;
                Count_eF_Pt20ToInf_Barrel_sysError += eventWeightError;
              } else {
                Count_eF_Pt20ToInf_Endcap += eventWeight;
                Count_eF_Pt20ToInf_Endcap_statError += eventWeight*eventWeight;
                Count_eF_Pt20ToInf_Endcap_sysError += eventWeightError;
              }
            } else {
              if (fabs(ele->scEta) < 1.479) {
                Count_eF_Pt10To20_Barrel += eventWeight;
                Count_eF_Pt10To20_Barrel_statError += eventWeight*eventWeight;
                Count_eF_Pt10To20_Barrel_sysError += eventWeightError;
              } else {
                Count_eF_Pt10To20_Endcap += eventWeight;
                Count_eF_Pt10To20_Endcap_statError += eventWeight*eventWeight;
                Count_eF_Pt10To20_Endcap_sysError += eventWeightError;
              }
            }
          } else {
            if (ele->pt > 20) {
              if (fabs(ele->scEta) < 1.479) {
                Count_mF_Pt20ToInf_Barrel += eventWeight;
                Count_mF_Pt20ToInf_Barrel_statError += eventWeight*eventWeight;
                Count_mF_Pt20ToInf_Barrel_sysError += eventWeightError;
              } else {
                Count_mF_Pt20ToInf_Endcap += eventWeight;
                Count_mF_Pt20ToInf_Endcap_statError += eventWeight*eventWeight;
                Count_mF_Pt20ToInf_Endcap_sysError += eventWeightError;
              }
            } else {
              if (fabs(ele->scEta) < 1.479) {
                Count_mF_Pt10To20_Barrel += eventWeight;
                Count_mF_Pt10To20_Barrel_statError += eventWeight*eventWeight;
                Count_mF_Pt10To20_Barrel_sysError += eventWeightError;
              } else {
                Count_mF_Pt10To20_Endcap += eventWeight;
                Count_mF_Pt10To20_Endcap_statError += eventWeight*eventWeight;
                Count_mF_Pt10To20_Endcap_sysError += eventWeightError;
              }
            }
          }

        }

        //*********************************************************************************************
        //Plots after all Cuts
        //*********************************************************************************************
        if (passAllCuts) {
//           fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton, eventWeight);
//           fMtLepton1_afterCuts->Fill(mTW[0], eventWeight);
//           fMtLepton2_afterCuts->Fill(mTW[1], eventWeight);
          fMtHiggs_afterCuts->Fill(mtHiggs, eventWeight);
//           fLeptonPtPlusMet_afterCuts->Fill(lepton1.Pt()+lepton2.Pt()+pfMet.Pt(), eventWeight);
    

        }

    } // end loop over electron denominators


  } //end loop over data     



  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Add systematics Uncertainties To Histgrams
  //==============================================================================================================
//   for (int i=0; i < fHWWSelection->GetXaxis()->GetNbins()+2; ++i) {
//     fHWWSelection->SetBinError(i, TMath::Sqrt(pow(fHWWSelection->GetBinError(i),2) + pow(fHWWSelection_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fHWWToEESelection->GetXaxis()->GetNbins()+2; ++i) {
//     fHWWToEESelection->SetBinError(i, TMath::Sqrt(pow(fHWWToEESelection->GetBinError(i),2) + pow(fHWWToEESelection_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fHWWToMuMuSelection->GetXaxis()->GetNbins()+2; ++i) {
//     fHWWToMuMuSelection->SetBinError(i, TMath::Sqrt(pow(fHWWToMuMuSelection->GetBinError(i),2) + pow(fHWWToMuMuSelection_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fHWWToEMuSelection->GetXaxis()->GetNbins()+2; ++i) {
//     fHWWToEMuSelection->SetBinError(i, TMath::Sqrt(pow(fHWWToEMuSelection->GetBinError(i),2) + pow(fHWWToEMuSelection_sysError->GetBinContent(i),2)));
//   }


//   for (int i=0; i < fLeptonEta->GetXaxis()->GetNbins()+2; ++i) {
//     fLeptonEta->SetBinError(i, TMath::Sqrt(pow(fLeptonEta->GetBinError(i),2) + pow(fLeptonEta_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fLeptonPtMax->GetXaxis()->GetNbins()+2; ++i) {
//     fLeptonPtMax->SetBinError(i, TMath::Sqrt(pow(fLeptonPtMax->GetBinError(i),2) + pow(fLeptonPtMax_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fLeptonPtMin->GetXaxis()->GetNbins()+2; ++i) {
//     fLeptonPtMin->SetBinError(i, TMath::Sqrt(pow(fLeptonPtMin->GetBinError(i),2) + pow(fLeptonPtMin_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fMetPtHist->GetXaxis()->GetNbins()+2; ++i) {
//     fMetPtHist->SetBinError(i, TMath::Sqrt(pow(fMetPtHist->GetBinError(i),2) + pow(fMetPtHist_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < fDeltaPhiLeptons->GetXaxis()->GetNbins()+2; ++i) {
//     fDeltaPhiLeptons->SetBinError(i, TMath::Sqrt(pow(fDeltaPhiLeptons->GetBinError(i),2) + pow(fDeltaPhiLeptons_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < dileptonMass->GetXaxis()->GetNbins()+2; ++i) {
//     dileptonMass->SetBinError(i, TMath::Sqrt(pow(dileptonMass->GetBinError(i),2) + pow(dileptonMass_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < dileptonMass->GetXaxis()->GetNbins()+2; ++i) {
//     dileptonMass_ee->SetBinError(i, TMath::Sqrt(pow(dileptonMass_ee->GetBinError(i),2) + pow(dileptonMass_ee_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < dileptonMass->GetXaxis()->GetNbins()+2; ++i) {
//     dileptonMass_mumu->SetBinError(i, TMath::Sqrt(pow(dileptonMass_mumu->GetBinError(i),2) + pow(dileptonMass_mumu_sysError->GetBinContent(i),2)));
//   }
//   for (int i=0; i < dileptonMass->GetXaxis()->GetNbins()+2; ++i) {
//     dileptonMass_emu->SetBinError(i, TMath::Sqrt(pow(dileptonMass_emu->GetBinError(i),2) + pow(dileptonMass_emu_sysError->GetBinContent(i),2)));
//   }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  dileptonMass->Draw();
  cv->SaveAs("dileptonMass.gif");
  dileptonMass_emu->Draw();
  cv->SaveAs("dileptonMass_emu.gif");
  dileptonMass_ee->Draw();
  cv->SaveAs("dileptonMass_ee.gif");
  dileptonMass_mumu->Draw();
  cv->SaveAs("dileptonMass_mumu.gif");

  fLeptonFakePt->Draw();
  cv->SaveAs("LeptonFakePt.gif");
  fLeptonFakeEta->Draw();
  cv->SaveAs("LeptonFakeEta.gif");
  fLeptonPtMin->Draw();
  cv->SaveAs("LeptonPtMin.gif");
  fLeptonEta->Draw();
  cv->SaveAs("LeptonEtaMin.gif");
  fLeptonPtMax->Draw();
  cv->SaveAs("LeptonPtMax.gif");
  fMetPtHist->Draw();
  cv->SaveAs("MetPtHist.gif");
  fDeltaPhiLeptons->Draw();
  cv->SaveAs("DeltaPhiLeptons.gif");

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  Double_t scaleFactor = 1000 / 35.5;
  scaleFactor=1;
  for (int i=1; i < fHWWToEESelection->GetXaxis()->GetNbins()+1; ++i) {
    cout << fHWWToEESelection->GetBinContent(i) << "+/-" << fHWWToEESelection->GetBinError(i) << "+/-" << fHWWToEESelection_sysError->GetBinContent(i) << " " << fHWWToMuMuSelection->GetBinContent(i) << "+/-" << fHWWToMuMuSelection->GetBinError(i) << "+/-" << fHWWToMuMuSelection_sysError->GetBinContent(i) << " " << fHWWToEMuSelection->GetBinContent(i) << "+/-" << fHWWToEMuSelection->GetBinError(i) << "+/-" << fHWWToEMuSelection_sysError->GetBinContent(i) << " " << fHWWSelection->GetBinContent(i) << "+/-" << fHWWSelection->GetBinError(i) << "+/-" << fHWWSelection_sysError->GetBinContent(i) << endl;
  }

  cout << "NEventsPreSel : " << NEventsPreSel << " +- " << TMath::Sqrt(NEventsPreSel_statError) << " +- " << NEventsPreSel_sysError << endl;
  cout << "NEventsPreSelPassHiggsCuts : " << NEventsPreSelPassHiggsCuts << " +- " << TMath::Sqrt(NEventsPreSelPassHiggsCuts_statError) << " +- " << NEventsPreSelPassHiggsCuts_sysError << endl;
  cout << "Eff = " << NEventsPreSelPassHiggsCuts / NEventsPreSel << endl;
  cout << endl;

  cout << "**************************************************************\n";
  cout << "Event Count : ee final state\n";
  cout << "BB :" << Count_ee << " +/- " << TMath::Sqrt(Count_ee_statError) << " +/- " << Count_ee_sysError <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : mm final state\n";
  cout << "BB :" << Count_mm << " +/- " << TMath::Sqrt(Count_mm_statError) << " +/- " << Count_mm_sysError << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : em final state\n";
  cout << "BB :" << Count_em << " +/- " << TMath::Sqrt(Count_em_statError) << " +/- " << Count_em_sysError<< endl;
  cout << "**************************************************************\n";
  cout << "Event Count : me final state\n";
  cout << "BB :" << Count_me << " +/- " << TMath::Sqrt(Count_me_statError) << " +/- " << Count_me_sysError<< endl;
  cout << "**************************************************************\n";
  cout << endl;

  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin0 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt10To15_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin1 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt10To15_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin2 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt10To15_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin0 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt15To20_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin1 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt15To20_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin2 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt15To20_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin0 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt20To25_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin1 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt20To25_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin2 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt20To25_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin0 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt25To30_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin1 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt25To30_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin2 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt25To30_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin0 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt30ToInf_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin1 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt30ToInf_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin2 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt30ToInf_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";


  cout << "**************************************************************\n";
  cout << "Event Count : Ele + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_eF_Pt10To20_Barrel << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Barrel_statError) << " +/- " << Count_eF_Pt10To20_Barrel_sysError <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_eF_Pt10To20_Endcap << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Endcap_statError) << " +/- " << Count_eF_Pt10To20_Endcap_sysError <<  endl;
  cout << "Barrel Pt [20+] : " << Count_eF_Pt20ToInf_Barrel << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Barrel_statError) << " +/- " << Count_eF_Pt20ToInf_Barrel_sysError <<  endl;
  cout << "Endcap Pt [20+] : " << Count_eF_Pt20ToInf_Endcap << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Endcap_statError) << " +/- " << Count_eF_Pt20ToInf_Endcap_sysError <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Mu + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_mF_Pt10To20_Barrel << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Barrel_statError) << " +/- " << Count_mF_Pt10To20_Barrel_sysError <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_mF_Pt10To20_Endcap << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Endcap_statError) << " +/- " << Count_mF_Pt10To20_Endcap_sysError <<  endl;
  cout << "Barrel Pt [20+] : " << Count_mF_Pt20ToInf_Barrel << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Barrel_statError) << " +/- " << Count_mF_Pt20ToInf_Barrel_sysError <<  endl;
  cout << "Endcap Pt [20+] : " << Count_mF_Pt20ToInf_Endcap << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Endcap_statError) << " +/- " << Count_mF_Pt20ToInf_Endcap_sysError <<  endl;
  cout << "**************************************************************\n";
  cout << endl;





  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("WWSelectionPlotsFakePrediction.root", "RECREATE");
//   TFile *file = new TFile("WWSelectionPlots_IsoDenominator.root", "RECREATE");

  file->WriteTObject(fHWWSelection, fHWWSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEESelection, fHWWToEESelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToMuMuSelection , fHWWToMuMuSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEMuSelection,fHWWToEMuSelection->GetName(), "WriteDelete");
  
  file->WriteTObject(fLeptonFakeEta ,fLeptonFakeEta->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePt ,fLeptonFakePt->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePtEta ,fLeptonFakePtEta->GetName(), "WriteDelete");


  file->WriteTObject(fLeptonEta ,fLeptonEta->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtMax ,fLeptonPtMax->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtMin ,fLeptonPtMin->GetName(), "WriteDelete");
  file->WriteTObject(fMetPtHist ,fMetPtHist->GetName(), "WriteDelete");
  file->WriteTObject(fMetPhiHist ,fMetPhiHist->GetName(), "WriteDelete");
  file->WriteTObject(fUncorrMetPtHist ,fUncorrMetPtHist->GetName(), "WriteDelete");
  file->WriteTObject(fUncorrMetPhiHist ,fUncorrMetPhiHist->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaPhiLeptons ,fDeltaPhiLeptons->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaEtaLeptons ,fDeltaEtaLeptons->GetName(), "WriteDelete");
                     
  file->WriteTObject(fMinDeltaPhiLeptonMet_afterCuts ,fMinDeltaPhiLeptonMet_afterCuts->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton1_afterCuts ,fMtLepton1_afterCuts->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton2_afterCuts ,fMtLepton2_afterCuts->GetName(), "WriteDelete");
  file->WriteTObject(fMtHiggs_afterCuts ,fMtHiggs_afterCuts->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtPlusMet_afterCuts ,fLeptonPtPlusMet_afterCuts->GetName(), "WriteDelete");
   
  file->WriteTObject(dileptonMass ,dileptonMass->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_ee ,dileptonMass_ee ->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_emu ,dileptonMass_emu->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_mumu ,dileptonMass_mumu->GetName(), "WriteDelete");
  file->Close();
  delete file;

  file = new TFile("ElectronFakeRate.SmurfV5.root", "UPDATE");
  file->WriteTObject(fLeptonJetPt, fLeptonJetPt->GetName(), "WriteDelete");
  file->Close();
  delete file;

     
  gBenchmark->Show("WWTemplate");       
} 




Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
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
            && fabs(ele->dz) < 0.2
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
             && fabs(ele->dz) < 0.2
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




Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) > 2.4) pass = kFALSE;

  Double_t iso04 = mu->ChargedIso04 + mu->NeutralIso04_10Threshold;
  Double_t iso03 = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
  Double_t iso = iso03;
  Double_t isoCutValue = 0;

  if (fabs(mu->eta) < 1.479) {
    if (mu->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (mu->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  } 

  if (! 
      ( (
          (Bool_t(mu->typeBits & kGlobal) 
           && mu->muNchi2 < 10.0
           && (mu->nValidHits > 0)
           && (mu->nMatch > 1 )
         )
           || 
          ( mu->typeBits & kTracker            
          && Bool_t(mu->qualityBits & kTMLastStationTight) 
            )
        )
                
        && mu->nTkHits > 10
        && (mu->nPixHits > 0)
        && fabs(mu->d0) < 0.02
        && fabs(mu->dz) < 0.1
        && iso / mu->pt < isoCutValue
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;


  if (mu->pt < 20) {
    if (!
        ( fabs(mu->d0) < 0.01
          )
      ) {
      pass = kFALSE;
    }    
  }

  return pass;
}







Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;


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
            && (ele->trkIso03 ) / ele->pt < 0.2
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
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->dz) < 0.1
             && (ele->trkIso03) / ele->pt < 0.2
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20
          )
      ) {
      pass = kFALSE;
    }
  } 

  return pass;
}



Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

//   //M1 Denominator 
//   if (! 
//       ( (
//           (Bool_t(mu->typeBits & kGlobal) 
//            && mu->muNchi2 < 10.0
//            && (mu->nValidHits > 0)
//            && (mu->nMatch > 1 )
//             )
//           || 
//           ( mu->typeBits & kTracker            
//             && Bool_t(mu->qualityBits & kTMLastStationTight) 
//             )
//         )
//         && mu->nTkHits > 10
//         && (mu->nPixHits > 0)
//         && fabs(mu->d0) < 0.2
//         && fabs(mu->dz) < 0.1
//         && (mu->ChargedIso03 + mu->NeutralIso03_10Threshold) / mu->pt < 1.0
//         && (mu->pterr / mu->pt < 0.1)
//         )
//     ) pass = kFALSE;    

  //M2 Denominator 
  if (! 
      ( (
          (Bool_t(mu->typeBits & kGlobal) 
           && mu->muNchi2 < 10.0
           && (mu->nValidHits > 0)
           && (mu->nMatch > 1 )
            )
          || 
          ( mu->typeBits & kTracker            
            && Bool_t(mu->qualityBits & kTMLastStationTight) 
            )
        )        
        && mu->nTkHits > 10        
        && ( mu->nPixHits > 0)
        && fabs(mu->d0) < 0.2
        && fabs(mu->dz) < 0.1
        && (mu->ChargedIso03 + mu->NeutralIso03_10Threshold) / mu->pt < 0.4
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;    


  return pass;
}



//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
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
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
