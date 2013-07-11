

// cat PassPassEventList.smurfV6.txt | grep "Pass" | sort -n --key=1 > ! k ; mv k PassPassEventList.smurfV6.txt
// cat EventList.passpass.sixie.txt | sort -n --key=1 > ! k ; mv k EventList.passpass.sixie.txt
// cat PassPassEventList.smurfV6.txt | awk '{print $1}' > ! pp.smurf.txt
// cat EventList.passpass.sixie.txt | awk '{print $1}' > ! pp.sixie.txt
// sdiff pp.smurf.txt pp.sixie.txt > ! diff
// cat diff | grep "<" | awk '{print "|| info->evtNum == " $1}'
// cat diff | grep ">" | awk '{print "|| info->evtNum == " $2}'

// cat PassFailEventList.smurfv6.txt | grep "Pass" | sort -n --key=1 > ! k ; mv k PassFailEventList.smurfv6.txt
// cat EventList.passfail.sixie.txt | sort -n --key=1 > ! k ; mv k EventList.passfail.sixie.txt
// cat PassFailEventList.smurfv6.txt | awk '{print $1}' > ! pf.smurf.txt
// cat EventList.passfail.sixie.txt | awk '{print $1}' > ! pf.sixie.txt
// sdiff pf.smurf.txt pf.sixie.txt > ! diff
// cat diff | grep "<" | awk '{print "|| info->evtNum == " $1}'
// cat diff | grep ">" | awk '{print "|| info->evtNum == " $2}'

//root -l EWKAna/Hww/FakeLeptonBkg/HwwFakeElectronPrediction_MC.C+\(\"/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjets-z2-v8-pu11-2l_noskim_normalized.root\",130,\"\"\)
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

Bool_t isRealLepton(Int_t value) {
  
  Bool_t answer = kFALSE;

  Int_t isEle;
  Int_t isMu;
  Int_t isTau;
  Int_t isPhoton;
  
  isEle = value % 2;
  Int_t tmp0 = floor(double(value) / 2.0);
  isMu = tmp0 % 2;
  Int_t tmp1 = floor(double(tmp0) / 2.0);
  isTau = tmp1 % 2;
  Int_t tmp2 = floor(double(tmp1) / 2.0);
  isPhoton = tmp2 % 2;
  
  if (isEle || isMu) answer = kTRUE;
  return answer;

}


//=== MAIN MACRO =================================================================================================

void HwwFakeElectronPrediction_MC(const string inputFilename, Double_t mHiggs, 
                                  const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  ofstream passpassFile("EventList.passpass.sixie.txt");
  ofstream passfailFile("EventList.passfail.sixie.txt");


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
//                                                           "ElectronFakeRate.SmurfV3.root",
//                                                       "ElectronFakeRate.SmurfV4M.root",
//                                                       "ElectronFakeRate.SmurfV4T.root",
//                                                       "ElectronFakeRate.SmurfV4Cone03L.root",
//                                                       "ElectronFakeRate.SmurfV4Cone03M.root",
//                                                        "ElectronFakeRate.SmurfV4Cone03T.root",
//                                                        "ElectronFakeRate.SmurfV4ConeHybridL.root",
//                                                        "ElectronFakeRate.SmurfV4ConeHybridM.root",
//                                                         "ElectronFakeRate.SmurfV4ConeHybridT.root",
//                                                           "ElectronFakeRate.SmurfV4TTruncated.root",
//                                                           "ElectronFakeRate.SmurfV4TTruncated20.root",
    "ElectronFakeRate.SmurfV6.skim.root",
    "MuonFakeRate.SmurfV6.skim.root",
    "", "",
    "ElectronFakeRateDenominatorV4_QCDMCCombined_ptThreshold30_PtEta",
    "MuonFakeRateDenominatorM3_QCDMCCombined_ptThreshold15_PtEta",
    use2DFakeRate, useFitFunction );



  //--------------------------------------------------------------------------------------------------------------
  // Histograms from Prediction
  //==============================================================================================================  
  TH1D *fLeptonFakePt_preSelection_predicted = new TH1D(       "hLeptonFakePt_preSelection_predicted",";Lepton Fake P_t Min;Number of Events",30,0.,150.);
  TH1D *fLeptonFakeEta_preSelection_predicted = new TH1D(      "hLeptonFakeEta_preSelection_predicted",";Lepton Fake #eta;Number of Events",20,-5.0,5.0);




  TH1D *fHWWSelection_predicted = new TH1D("hHWWSelection_predicted", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection_predicted = new TH1D("hHWWToEESelection_predicted", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection_predicted = new TH1D("hHWWToMuMuSelection_predicted", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection_predicted = new TH1D("hHWWToEMuSelection_predicted", ";Cut Number;Number of Events", 15, -1.5, 13.5);

  TH1D *fLeptonFakePt_WWSelection_predicted = new TH1D(       "hLeptonFakePt_WWSelection_predicted",";Lepton Fake P_t Min;Number of Events",15,0.,150.);
  TH1D *fLeptonFakeEta_WWSelection_predicted = new TH1D(      "hLeptonFakeEta_WWSelection_predicted",";Lepton Fake #eta;Number of Events",11,-5.5,5.5);
  TH2D *fLeptonFakePtEta_WWSelection_predicted = new TH2D ( "hLeptonFakePtEta_WWSelection_predicted", ";Lepton Fake p_{T}; Lepton Fake #eta; Number of Events", 20,0.,100., 30,-5.0,5.0);

  TH1D *fLeptonPtMax_WWSelection_predicted = new TH1D(        "hLeptonPtMax_WWSelection_predicted",";Lepton P_t Max;Number of Events",15,0.,150.);
  TH1D *fLeptonPtMin_WWSelection_predicted = new TH1D(        "hLeptonPtMin_WWSelection_predicted",";Lepton P_t Min;Number of Events",15,0.,150.);
  TH1D *fMetPtHist_WWSelection_predicted = new TH1D(          "hMetPtHist_WWSelection_predicted",";Met;Number of Events",30,0.,300.);  
  TH1D *fDeltaPhiLeptons_WWSelection_predicted = new TH1D(    "hDeltaPhiLeptons_WWSelection_predicted",";#Delta#phi_{ll};Number of Events",18,0,180);
  TH1D *fDeltaEtaLeptons_WWSelection_predicted = new TH1D(    "hDeltaEtaLeptons_WWSelection_predicted",";#Delta#eta_{ll};Number of Events",20,-5.0,5.0);


  TH1D *fMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted = new TH1D(    "hMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted", 
                                                                             ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  TH1D *fMtLepton1_afterCuts_WWSelection_predicted = new TH1D(               "hMtLepton1_afterCuts_WWSelection_predicted",
                                                                             ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  TH1D *fMtLepton2_afterCuts_WWSelection_predicted = new TH1D(               "hMtLepton2_afterCuts_WWSelection_predicted",
                                                                             ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  TH1D *fMtHiggs_afterCuts_WWSelection_predicted = new TH1D(                 "hMtHiggs_afterCuts_WWSelection_predicted",
                                                                             ";M_t (l1+l2+Met);Number of Events",30,0.,300.);
  TH1D *fLeptonPtPlusMet_afterCuts_WWSelection_predicted = new TH1D(         "hLeptonPtPlusMet_afterCuts_WWSelection_predicted",
                                                                             ";LeptonPtPlusMet;Number of Events",150,0., 300.);

  TH1D *dileptonMass_WWSelection_predicted = new TH1D("dileptonMass_WWSelection_predicted", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_ee_WWSelection_predicted = new TH1D("dileptonMass_ee_WWSelection_predicted", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_emu_WWSelection_predicted = new TH1D("dileptonMass_emu_WWSelection_predicted", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_mumu_WWSelection_predicted = new TH1D("dileptonMass_mumu_WWSelection_predicted", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);

  Double_t Count_ee_WWSelection_predicted = 0;
  Double_t Count_mm_WWSelection_predicted = 0;
  Double_t Count_em_WWSelection_predicted = 0;
  Double_t Count_me_WWSelection_predicted = 0;
  Double_t Count_ee_statError_WWSelection_predicted = 0;
  Double_t Count_mm_statError_WWSelection_predicted = 0;
  Double_t Count_em_statError_WWSelection_predicted = 0;
  Double_t Count_me_statError_WWSelection_predicted = 0;

  Double_t Count_Pt10To15_EtaBin0_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin1_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin2_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin0_statError_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin1_statError_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin2_statError_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin0_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin1_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt10To15_EtaBin2_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin0_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin1_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin2_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin0_statError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin1_statError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin2_statError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin0_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin1_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt15To20_EtaBin2_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin0_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin1_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin2_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin0_statError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin1_statError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin2_statError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin0_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin1_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt20To25_EtaBin2_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin0_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin1_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin2_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin0_statError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin1_statError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin2_statError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin0_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin1_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt25To30_EtaBin2_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin0_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin1_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin2_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin0_statError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin1_statError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin2_statError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin0_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin1_sysError_WWSelection_predicted = 0;
  Double_t Count_Pt30ToInf_EtaBin2_sysError_WWSelection_predicted = 0;

  Double_t Count_eF_Pt10To20_Barrel_WWSelection_predicted = 0;
  Double_t Count_eF_Pt10To20_Endcap_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Barrel_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Endcap_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_WWSelection_predicted = 0;
  Double_t Count_eF_Pt10To20_Barrel_statError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt10To20_Endcap_statError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_statError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_statError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Barrel_statError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Endcap_statError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_statError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_statError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt10To20_Barrel_sysError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt10To20_Endcap_sysError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_sysError_WWSelection_predicted = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_sysError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Barrel_sysError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt10To20_Endcap_sysError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_sysError_WWSelection_predicted = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_sysError_WWSelection_predicted = 0;


  fHWWSelection_predicted->Sumw2();
  fHWWToEESelection_predicted->Sumw2();
  fHWWToMuMuSelection_predicted->Sumw2();
  fHWWToEMuSelection_predicted->Sumw2();
  fLeptonPtMax_WWSelection_predicted->Sumw2();
  fLeptonPtMin_WWSelection_predicted->Sumw2();
  fLeptonFakePt_WWSelection_predicted->Sumw2();
  fLeptonFakeEta_WWSelection_predicted->Sumw2();
  fMetPtHist_WWSelection_predicted->Sumw2();
  fDeltaPhiLeptons_WWSelection_predicted->Sumw2();
  fDeltaEtaLeptons_WWSelection_predicted->Sumw2();
  fMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted->Sumw2();
  fMtLepton1_afterCuts_WWSelection_predicted->Sumw2();
  fMtLepton2_afterCuts_WWSelection_predicted->Sumw2();
  fMtHiggs_afterCuts_WWSelection_predicted->Sumw2();
  fLeptonPtPlusMet_afterCuts_WWSelection_predicted->Sumw2();
  dileptonMass_WWSelection_predicted->Sumw2();
  dileptonMass_ee_WWSelection_predicted->Sumw2();
  dileptonMass_emu_WWSelection_predicted->Sumw2();
  dileptonMass_mumu_WWSelection_predicted->Sumw2();

  //--------------------------------------------------------------------------------------------------------------
  // Histograms from FullSim
  //==============================================================================================================  
  TH1D *fLeptonFakePt_preSelection_fullSim = new TH1D(       "hLeptonFakePt_preSelection_fullSim",";Lepton Fake P_t Min;Number of Events",30,0.,150.);
  TH1D *fLeptonFakeEta_preSelection_fullSim = new TH1D(      "hLeptonFakeEta_preSelection_fullSim",";Lepton Fake #eta;Number of Events",20,-5.0,5.0);

  TH1D *fHWWSelection_fullSim = new TH1D("hHWWSelection_fullSim", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection_fullSim = new TH1D("hHWWToEESelection_fullSim", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection_fullSim = new TH1D("hHWWToMuMuSelection_fullSim", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection_fullSim = new TH1D("hHWWToEMuSelection_fullSim", ";Cut Number;Number of Events", 15, -1.5, 13.5);

  TH1D *fLeptonFakePt_WWSelection_fullSim = new TH1D(       "hLeptonFakePt_WWSelection_fullSim",";Lepton Fake P_t Min;Number of Events",15,0.,150.);
  TH1D *fLeptonFakeEta_WWSelection_fullSim = new TH1D(      "hLeptonFakeEta_WWSelection_fullSim",";Lepton Fake #eta;Number of Events",11,-5.5,5.5);
  TH2D *fLeptonFakePtEta_WWSelection_fullSim = new TH2D ( "hLeptonFakePtEta_WWSelection_fullSim", ";Lepton Fake p_{T}; Lepton Fake #eta; Number of Events", 20,0.,100., 30,-5.0,5.0);

  TH1D *fLeptonPtMax_WWSelection_fullSim = new TH1D(        "hLeptonPtMax_WWSelection_fullSim",";Lepton P_t Max;Number of Events",15,0.,150.);
  TH1D *fLeptonPtMin_WWSelection_fullSim = new TH1D(        "hLeptonPtMin_WWSelection_fullSim",";Lepton P_t Min;Number of Events",15,0.,150.);
  TH1D *fMetPtHist_WWSelection_fullSim = new TH1D(          "hMetPtHist_WWSelection_fullSim",";Met;Number of Events",30,0.,300.);  
  TH1D *fDeltaPhiLeptons_WWSelection_fullSim = new TH1D(    "hDeltaPhiLeptons_WWSelection_fullSim",";#Delta#phi_{ll};Number of Events",18,0,180);
  TH1D *fDeltaEtaLeptons_WWSelection_fullSim = new TH1D(    "hDeltaEtaLeptons_WWSelection_fullSim",";#Delta#eta_{ll};Number of Events",20,-5.0,5.0);


  TH1D *fMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim = new TH1D(    "hMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim", 
                                                                           ";Min #Delta#phi_{l,Met};Number of Events",90,0.,180);
  TH1D *fMtLepton1_afterCuts_WWSelection_fullSim = new TH1D(               "hMtLepton1_afterCuts_WWSelection_fullSim",
                                                                           ";M_t (Lepton1,Met);Number of Events",100,0.,200.);
  TH1D *fMtLepton2_afterCuts_WWSelection_fullSim = new TH1D(               "hMtLepton2_afterCuts_WWSelection_fullSim",
                                                                           ";M_t (Lepton2,Met);Number of Events",100,0.,200.);
  TH1D *fMtHiggs_afterCuts_WWSelection_fullSim = new TH1D(                 "hMtHiggs_afterCuts_WWSelection_fullSim",
                                                                           ";M_t (l1+l2+Met);Number of Events",30,0.,300.);
  TH1D *fLeptonPtPlusMet_afterCuts_WWSelection_fullSim = new TH1D(         "hLeptonPtPlusMet_afterCuts_WWSelection_fullSim",
                                                                           ";LeptonPtPlusMet;Number of Events",150,0., 300.);

  TH1D *dileptonMass_WWSelection_fullSim = new TH1D("dileptonMass_WWSelection_fullSim", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_ee_WWSelection_fullSim = new TH1D("dileptonMass_ee_WWSelection_fullSim", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_emu_WWSelection_fullSim = new TH1D("dileptonMass_emu_WWSelection_fullSim", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);
  TH1D *dileptonMass_mumu_WWSelection_fullSim = new TH1D("dileptonMass_mumu_WWSelection_fullSim", "; Mass [GeV/c^{2}]; Number of Events", 20, 0, 200);

  Double_t Count_ee_WWSelection_fullSim = 0;
  Double_t Count_mm_WWSelection_fullSim = 0;
  Double_t Count_em_WWSelection_fullSim = 0;
  Double_t Count_me_WWSelection_fullSim = 0;
  Double_t Count_ee_statError_WWSelection_fullSim = 0;
  Double_t Count_mm_statError_WWSelection_fullSim = 0;
  Double_t Count_em_statError_WWSelection_fullSim = 0;
  Double_t Count_me_statError_WWSelection_fullSim = 0;

  Double_t Count_Pt10To15_EtaBin0_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin1_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin2_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin0_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin1_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin2_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin0_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin1_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt10To15_EtaBin2_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin0_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin1_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin2_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin0_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin1_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin2_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin0_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin1_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt15To20_EtaBin2_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin0_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin1_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin2_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin0_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin1_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin2_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin0_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin1_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt20To25_EtaBin2_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin0_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin1_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin2_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin0_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin1_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin2_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin0_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin1_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt25To30_EtaBin2_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin0_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin1_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin2_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin0_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin1_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin2_statError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin0_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin1_sysError_WWSelection_fullSim = 0;
  Double_t Count_Pt30ToInf_EtaBin2_sysError_WWSelection_fullSim = 0;

  Double_t Count_eF_Pt10To20_Barrel_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt10To20_Endcap_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Barrel_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Endcap_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt10To20_Barrel_statError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt10To20_Endcap_statError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_statError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_statError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Barrel_statError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Endcap_statError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_statError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_statError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt10To20_Barrel_sysError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt10To20_Endcap_sysError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Barrel_sysError_WWSelection_fullSim = 0;
  Double_t Count_eF_Pt20ToInf_Endcap_sysError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Barrel_sysError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt10To20_Endcap_sysError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Barrel_sysError_WWSelection_fullSim = 0;
  Double_t Count_mF_Pt20ToInf_Endcap_sysError_WWSelection_fullSim = 0;


  fHWWSelection_fullSim->Sumw2();
  fHWWToEESelection_fullSim->Sumw2();
  fHWWToMuMuSelection_fullSim->Sumw2();
  fHWWToEMuSelection_fullSim->Sumw2();
  fLeptonPtMax_WWSelection_fullSim->Sumw2();
  fLeptonPtMin_WWSelection_fullSim->Sumw2();
  fLeptonFakePt_WWSelection_fullSim->Sumw2();
  fLeptonFakeEta_WWSelection_fullSim->Sumw2();
  fMetPtHist_WWSelection_fullSim->Sumw2();
  fDeltaPhiLeptons_WWSelection_fullSim->Sumw2();
  fDeltaEtaLeptons_WWSelection_fullSim->Sumw2();
  fMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim->Sumw2();
  fMtLepton1_afterCuts_WWSelection_fullSim->Sumw2();
  fMtLepton2_afterCuts_WWSelection_fullSim->Sumw2();
  fMtHiggs_afterCuts_WWSelection_fullSim->Sumw2();
  fLeptonPtPlusMet_afterCuts_WWSelection_fullSim->Sumw2();
  dileptonMass_WWSelection_fullSim->Sumw2();
  dileptonMass_ee_WWSelection_fullSim->Sumw2();
  dileptonMass_emu_WWSelection_fullSim->Sumw2();
  dileptonMass_mumu_WWSelection_fullSim->Sumw2();

  vector<string> CutLabel;
  CutLabel.push_back("Dilepton20/10");
  CutLabel.push_back("PreselectionPtCuts");
  CutLabel.push_back("zDiffMaxCut");
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
  CutLabel.push_back("mtHiggsCut");
  
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
    if ((0==1)

|| info->evtNum == 7755124
|| info->evtNum == 19668349
|| info->evtNum == 215345879
|| info->evtNum == 288792650


|| info->evtNum == 6086587
|| info->evtNum == 71683939
|| info->evtNum == 251947246
|| info->evtNum == 295317759

      ) printDebug = kTRUE;


    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;


    if (printDebug) {
      cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
      cout << "PV: " << info->pvx << " " << info->pvy << " " << info->pvz << endl;
      cout << "VGamma : " << info->VGammaEvent << endl;
    }

//     //do not consider VGamma
//     if (info->VGammaEvent) continue;

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
           fabs(ele->eta) <= 2.5
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
    //Debug
    //******************************************************************************
    if (printDebug) {

      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
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
            



        cout << "Electron" << i << " : " << ele->pt << " " << ele->eta << " " << ele->phi << " "
             << ele->isEB << " " 
             << ele->q << " : " 
             << ele->sigiEtaiEta << " " 
             << ele->deltaEtaIn << " " 
             << ele->deltaPhiIn << " " 
             << ele->HoverE << " " 
             << ele->trkIso03 << " " 
             << (TMath::Max(ele->emIso03 - 1.0, 0.0)) << " " 
             << ele->hadIso03 << " " 

             << ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold << " "
             << ele->ChargedIso04 << " " 
             << ele->NeutralHadronIso04_10Threshold << " " 
             << ele->GammaIso04_10Threshold << " " 
             << ele->GammaIsoVetoEtaStrip04_10Threshold << " " 
             << ele->ChargedEMIsoVetoEtaStrip04 << " " 
             << ele->NeutralHadronIso007_10Threshold << " " 
             << " : " 
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
          
        cout << Bool_t(ele->sigiEtaiEta < 0.03) << " " 
             << Bool_t(fabs(ele->deltaEtaIn) < 0.007 ) << " "
             << Bool_t(fabs(ele->deltaPhiIn) < 0.03 ) << " "
//                  << Bool_t((ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1 ) << " "
             << Bool_t((ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->NeutralHadronIso007_10Threshold) / ele->pt < isoCutValue ) << " "
             << Bool_t(ele->nExpHitsInner <= 0 ) << " "
             << Bool_t(passConversionVeto(ele->isConv) ) << " "
             << Bool_t(fabs(ele->d0) < 0.02 ) << " "
             << Bool_t(fabs(ele->dz) < 0.1 ) << " "
             << Bool_t(ele->isEcalDriven ) << " : "
             
             << passElectronCuts(ele) << " "
             << passElectronDenominatorCuts(ele) << " "
             << isRealLepton(ele->isMCReal) << " " 
             << endl;


        if (0==0) {
          Bool_t pass = kTRUE;


          cout << pass << " ";
          
          if (fabs(ele->eta) > 2.5) pass = kFALSE;
        
          cout << pass << " ";

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


          cout << pass << " ";

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

          cout << pass << " ";

            if (ele->fBrem <= 0.15) {
              if (fabs(ele->eta) > 1.0) {
                pass = kFALSE;
              } else {
                if (!( ele->EOverP > 0.95 )) pass = kFALSE;
              }
            }
          cout << pass << " ";

          }
          cout << pass << " ";
          cout << endl;
        }

      }
          
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
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
             << mu->ChargedIso03 + mu->NeutralIso03_10Threshold << " " 
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
             << Bool_t(fabs(mu->dz) < 0.1 ) << " "
             << Bool_t((mu->ChargedIso03 + mu->NeutralIso03_10Threshold ) / mu->pt < isoCutValue ) << " "
//                  << Bool_t((mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15 ) << " "
             << Bool_t(mu->nMatch > 1 ) << " "
             << Bool_t(mu->nPixHits > 0) << " "
             << Bool_t(mu->pterr / mu->pt < 0.1) << " "
             << " : " << passMuonCuts(mu) << " "
             << isRealLepton(mu->isMCReal) << " " 
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

      cout << "Leptons : " << leptonPt.size() << endl;
    }












    //******************************************************************************
    //For Prediction: select events with 1lepton only 
    //******************************************************************************
    if (!(leptonPt.size() < 1) && (leptonPt[0] > 10.0) 
        && (leptonPt.size() < 2 || leptonPt[1] < 10.0)
      ) {

      if  (printDebug) cout << "PASS FAIL\n";

      //loop over electron denominators
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
   
        if (!(ele->pt > 10.0 && fabs(ele->eta) < 2.5)) continue;
        if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[0],leptonEta[0]) < 0.1 ) {
          continue;
        }


        //select denominators that fail final selection
        if (!passElectronDenominatorCuts(ele)) continue;
        if (passElectronCuts(ele)) continue;


       //preselection
        if (!(ele->pt > 10.0 && fabs(ele->eta) < 2.5)) continue;
        if ((ChargeSelection == 0 && leptonCharge[0] == ele->q) || (ChargeSelection == 1 && leptonCharge[0] != ele->q)) continue;


//          if (! (ele->pt > 20)) continue;
//         if (! (ele->pt > 15 && ele->pt < 20)) continue;
//           if (! (ele->pt > 10 && ele->pt < 15)) continue;
//        if (! (fabs(ele->eta) > 0 && fabs(ele->eta) < 1.0)) continue;
//       if (! (fabs(ele->eta) > 1.0 && fabs(ele->eta) < 1.5)) continue;
//       if (! (fabs(ele->eta) > 1.5 && fabs(ele->eta) < 2.5)) continue;
    

        //******************************************************************************
        //Use only Fake Electron Events
        //******************************************************************************
        Bool_t isFakeEle = kFALSE;
        if (!isRealLepton(ele->isMCReal)) isFakeEle = kTRUE;
        if (!isFakeEle) continue;


  
        //******************************************************************************
        //Calculate Fake Rate
        //******************************************************************************
        double fakeRate = fFakeRate->ElectronFakeRate(ele->pt, fabs(ele->eta), ele->phi) / (1-fFakeRate->ElectronFakeRate(ele->pt, fabs(ele->eta), ele->phi));
        double fakeRateErrorLow = fFakeRate->ElectronFakeRateStatErrorLow(ele->pt, fabs(ele->eta), ele->phi) / pow((1- fFakeRate->ElectronFakeRate(ele->pt, fabs(ele->eta), ele->phi)),2);
        double fakeRateErrorHigh = fFakeRate->ElectronFakeRateStatErrorHigh(ele->pt, fabs(ele->eta), ele->phi) / pow((1- fFakeRate->ElectronFakeRate(ele->pt, fabs(ele->eta), ele->phi)),2);

        Double_t eventWeight = fakeRate; 
        //eventWeight = 1.0; 
        Double_t eventWeightError = (  fakeRateErrorLow +   fakeRateErrorHigh) / 2;
        // eventWeightError = 0;

        //******************************************************************************
        //Correct the MET for the leptons
        //******************************************************************************
        if  (printDebug) cout << "PASS FAIL\n";
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
        double maxBtag = -99999;
        NJets = 0;
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);
//           if (jet->pt < 7.0) continue;

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
            } else {
              if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
            }
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
        zDiffMax = fabs(dz - dzLepton);


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
        double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt()*pfMet.Pt()*(1.0 - cos(deltaPhiDileptonMet)));
      

        //mtHiggs assuming neutrino pt == dilepton pt
        double enell = TMath::Sqrt(dilepton.Pt()*dilepton.Pt() + dilepton.M()*dilepton.M());
        double enenn = TMath::Sqrt(pfMet.Pt()* pfMet.Pt()  + dilepton.M()*dilepton.M());
        double enex  = dilepton.Px() + pfMet.Px();
        double eney  = dilepton.Py() + pfMet.Py();
        double mll   = dilepton.M();
        double mnu   = mll;
        double mtHiggs_NuPtEqualDileptonPt = TMath::Sqrt(mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney));
      


        //*********************************************************************************************
        //Define Cuts
        //*********************************************************************************************
        const int nCuts = 15;
        bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                               false, false, false, false, false};

        Bool_t PreselPtCut = kTRUE;
        if(!(lepton1.Pt() >  20.0 && lepton2.Pt() > 10.0)) PreselPtCut = kFALSE;

        if(PreselPtCut) passCut[0] = true;
      
        if(zDiffMax < 100000.0) passCut[1] = true;            

        if(pfMet.Pt()    > 20.0)               passCut[2] = true; 
      
        if(dilepton.M() > 12.0)            passCut[3] = true;
   
        if (finalState == 0 || finalState == 1) { // mumu/ee
          if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
          if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 35) passCut[5] = true;
        }
        else if(finalState == 2 ||finalState == 3 ) { // emu
          passCut[4] = true;
          if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 20) passCut[5] = true;
        }
      
        if(NJets     == 0)              passCut[6] = true;
       
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
          cout << "Final state: " << finalState << endl;
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
        fHWWSelection_predicted->Fill(-1, eventWeight);
        if (finalState == 1 ) {
          fHWWToMuMuSelection_predicted->Fill(-1, eventWeight);
        } else if(finalState == 0 ) {
          fHWWToEESelection_predicted->Fill(-1, eventWeight);
        } else { //if(finalState == 2 || finalState == 3 ) {
          fHWWToEMuSelection_predicted->Fill(-1, eventWeight);
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
            fHWWSelection_predicted->Fill(k, eventWeight);
            if (finalState == 1 ) {
              fHWWToMuMuSelection_predicted->Fill(k, eventWeight);
            } else if(finalState == 0 ) {
              fHWWToEESelection_predicted->Fill(k, eventWeight);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelection_predicted->Fill(k, eventWeight);
            }
          }
        }

        //*****************************************************************************************
        //Print EventList
        //*****************************************************************************************
        if (
          passCut[0] && passCut[1] 
          && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//           && passCut[10] && passCut[11] && passCut[12] && passCut[13] && passCut[14]
          ) {

          if ((0==0) 
//              && finalState == 0
            ) {
            passfailFile << info->evtNum << " : " << ele->pt << " " << ele->eta << " : " << fakeRate  << endl;
          }
        }

        //*****************************************************************************************
        //Make Preselection Histograms  
        //*****************************************************************************************
        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] ) {
          fLeptonFakePt_preSelection_predicted->Fill(ele->pt, eventWeight);
          fLeptonFakeEta_preSelection_predicted->Fill(ele->eta, eventWeight);
        }        


        //*****************************************************************************************
        //Make WW Selection Histograms  
        //*****************************************************************************************
//         if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[5]) {
        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
//         if (passAllCuts) {


          fLeptonFakePt_WWSelection_predicted->Fill(ele->pt, eventWeight);
          fLeptonFakeEta_WWSelection_predicted->Fill(ele->eta, eventWeight);
          fLeptonFakePtEta_WWSelection_predicted->Fill(ele->pt, ele->eta, eventWeight);

          fLeptonPtMax_WWSelection_predicted->Fill(lepton1.Pt(), eventWeight);
          fLeptonPtMin_WWSelection_predicted->Fill(lepton2.Pt(), eventWeight);
          fMetPtHist_WWSelection_predicted->Fill(pfMet.Pt(), eventWeight);                             
          fDeltaPhiLeptons_WWSelection_predicted->Fill(deltaPhiLeptons, eventWeight);

          dileptonMass_WWSelection_predicted->Fill(dilepton.M(), eventWeight);
          if (finalState == 0) {
            dileptonMass_ee_WWSelection_predicted->Fill(dilepton.M(), eventWeight);
          } else if (finalState == 1) {
            dileptonMass_mumu_WWSelection_predicted->Fill(dilepton.M(), eventWeight);
          } else {
            dileptonMass_emu_WWSelection_predicted->Fill(dilepton.M(), eventWeight);
          }

          fMtHiggs_afterCuts_WWSelection_predicted->Fill(mtHiggs, eventWeight);

        }

        if (
          passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//              passAllCuts           
          ) {


          //  if (passAllCuts ) {
          if (finalState == 0) {
            Count_ee_WWSelection_predicted += eventWeight;
            Count_ee_statError_WWSelection_predicted += eventWeight*eventWeight;
          } else if (finalState == 1) {
            Count_mm_WWSelection_predicted += eventWeight;
            Count_mm_statError_WWSelection_predicted += eventWeight*eventWeight;
          } else if (finalState == 2) {
            Count_em_WWSelection_predicted += eventWeight;
            Count_em_statError_WWSelection_predicted += eventWeight*eventWeight;
          } else {
            Count_me_WWSelection_predicted += eventWeight;
            Count_me_statError_WWSelection_predicted += eventWeight*eventWeight;
          }

          Double_t ptCategory = lepton2.Pt();
          ptCategory = ele->pt;

          if (ptCategory > 30.0) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt30ToInf_EtaBin0_WWSelection_predicted += eventWeight;
              Count_Pt30ToInf_EtaBin0_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt30ToInf_EtaBin1_WWSelection_predicted += eventWeight;
              Count_Pt30ToInf_EtaBin1_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt30ToInf_EtaBin2_WWSelection_predicted += eventWeight;
              Count_Pt30ToInf_EtaBin2_statError_WWSelection_predicted += eventWeight*eventWeight;
            }
          } else if (ptCategory > 25.0 && ptCategory < 30.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt25To30_EtaBin0_WWSelection_predicted += eventWeight;
              Count_Pt25To30_EtaBin0_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt25To30_EtaBin1_WWSelection_predicted += eventWeight;
              Count_Pt25To30_EtaBin1_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt25To30_EtaBin2_WWSelection_predicted += eventWeight;
              Count_Pt25To30_EtaBin2_statError_WWSelection_predicted += eventWeight*eventWeight;
            }
          } else if (ptCategory > 20.0 && ptCategory < 25.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt20To25_EtaBin0_WWSelection_predicted += eventWeight;
              Count_Pt20To25_EtaBin0_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt20To25_EtaBin1_WWSelection_predicted += eventWeight;
              Count_Pt20To25_EtaBin1_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt20To25_EtaBin2_WWSelection_predicted += eventWeight;
              Count_Pt20To25_EtaBin2_statError_WWSelection_predicted += eventWeight*eventWeight;
            }
          } else if (ptCategory > 15.0 && ptCategory < 20.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt15To20_EtaBin0_WWSelection_predicted += eventWeight;
              Count_Pt15To20_EtaBin0_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt15To20_EtaBin1_WWSelection_predicted += eventWeight;
              Count_Pt15To20_EtaBin1_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt15To20_EtaBin2_WWSelection_predicted += eventWeight;
              Count_Pt15To20_EtaBin2_statError_WWSelection_predicted += eventWeight*eventWeight;
            }
          } else if (ptCategory > 10.0 && ptCategory < 15.0 ) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt10To15_EtaBin0_WWSelection_predicted += eventWeight;
              Count_Pt10To15_EtaBin0_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt10To15_EtaBin1_WWSelection_predicted += eventWeight;
              Count_Pt10To15_EtaBin1_statError_WWSelection_predicted += eventWeight*eventWeight;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt10To15_EtaBin2_WWSelection_predicted += eventWeight;
              Count_Pt10To15_EtaBin2_statError_WWSelection_predicted += eventWeight*eventWeight;
            }
          } 


          if (leptonType[0] == 11) {
            if (ele->pt > 20) {
              if (fabs(ele->scEta) < 1.479) {
                Count_eF_Pt20ToInf_Barrel_WWSelection_predicted += eventWeight;
                Count_eF_Pt20ToInf_Barrel_statError_WWSelection_predicted += eventWeight*eventWeight;
              } else {
                Count_eF_Pt20ToInf_Endcap_WWSelection_predicted += eventWeight;
                Count_eF_Pt20ToInf_Endcap_statError_WWSelection_predicted += eventWeight*eventWeight;
              }
            } else {
              if (fabs(ele->scEta) < 1.479) {
                Count_eF_Pt10To20_Barrel_WWSelection_predicted += eventWeight;
                Count_eF_Pt10To20_Barrel_statError_WWSelection_predicted += eventWeight*eventWeight;
              } else {
                Count_eF_Pt10To20_Endcap_WWSelection_predicted += eventWeight;
                Count_eF_Pt10To20_Endcap_statError_WWSelection_predicted += eventWeight*eventWeight;
              }
            }
          } else {
            if (ele->pt > 20) {
              if (fabs(ele->scEta) < 1.479) {
                Count_mF_Pt20ToInf_Barrel_WWSelection_predicted += eventWeight;
                Count_mF_Pt20ToInf_Barrel_statError_WWSelection_predicted += eventWeight*eventWeight;

              } else {
                Count_mF_Pt20ToInf_Endcap_WWSelection_predicted += eventWeight;
                Count_mF_Pt20ToInf_Endcap_statError_WWSelection_predicted += eventWeight*eventWeight;
              }
            } else {
              if (fabs(ele->scEta) < 1.479) {
                Count_mF_Pt10To20_Barrel_WWSelection_predicted += eventWeight;
                Count_mF_Pt10To20_Barrel_statError_WWSelection_predicted += eventWeight*eventWeight;
              } else {
                Count_mF_Pt10To20_Endcap_WWSelection_predicted += eventWeight;
                Count_mF_Pt10To20_Endcap_statError_WWSelection_predicted += eventWeight*eventWeight;
              }
            }
          }

        }

        //*********************************************************************************************
        //Plots after all Cuts
        //*********************************************************************************************
        if (passAllCuts) {
//           fMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted->Fill(minDeltaPhiMetLepton, eventWeight);
//           fMtLepton1_afterCuts_WWSelection_predicted->Fill(mTW[0], eventWeight);
//           fMtLepton2_afterCuts_WWSelection_predicted->Fill(mTW[1], eventWeight);
//           fLeptonPtPlusMet_afterCuts_WWSelection_predicted->Fill(lepton1.Pt()+lepton2.Pt()+pfMet.Pt(), eventWeight);
    

        }

      } // end loop over electron denominators
    }


    //******************************************************************************
    //For FullSim: MC selection
    //******************************************************************************
    if (
      !(leptonPt.size() < 2) 
      ) {


      //******************************************************************************
      //Correct the MET for the leptons
      //******************************************************************************
      if  (printDebug) cout << "PASS PASS\n";
      if  (printDebug) cout << "Track MET Before Correction : " << info->pfTrackMEx  << " , " << info->pfTrackMEy << " : " << pfTrackMet.Pt() << endl;
        

      Double_t MetCorrection_X = 0;
      Double_t MetCorrection_Y = 0;
      for (Int_t l = 0; l < leptonPt.size(); ++l) {
        if (leptonType[l] == 11) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[l]]);           
          mithep::FourVectorM lepton;
          lepton.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
          mithep::FourVectorM pfLepton;
          pfLepton.SetCoordinates(ele->pfPt, ele->pfEta, ele->pfPhi, 0.51099892e-3 );
          if (ele->pfPt >= 0) {
            MetCorrection_X += ( pfLepton.Px() - lepton.Px() );
            MetCorrection_Y += ( pfLepton.Py() - lepton.Py() );       
          } else {
            MetCorrection_X += ( 0.0 - lepton.Px() );
            MetCorrection_Y += ( 0.0 - lepton.Py() );       
          }
          if  (printDebug) cout << "Ele : " << pfLepton.Px() << " , " << pfLepton.Py() << " : " << lepton.Px() << " , " << lepton.Py() << endl;
        } else if (leptonType[l] == 13) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[l]]);
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
      }
      if  (printDebug)cout << "Met Correction : " << MetCorrection_X << " , " << MetCorrection_Y << endl;

      pfTrackMet.SetXYZ(info->pfTrackMEx + MetCorrection_X, info->pfTrackMEy + MetCorrection_Y , 0);

      if  (printDebug)cout << "Corrected Met : " << pfTrackMet.Px() << " , " << pfTrackMet.Py() << " : " << pfTrackMet.Pt() << endl;



      //******************************************************************************
      //Count Jets
      //******************************************************************************
      double maxBtag = -99999;
      NJets = 0;
      for(Int_t j=0; j<jetArr->GetEntries(); j++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);
//         if (jet->pt < 7.0) continue;

        Bool_t leptonOverlap = kFALSE;
        for (int k=0; k<leptonPt.size(); ++k) {            
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.5) {
            leptonOverlap = kTRUE;
          }
        }
                     
        if (!leptonOverlap) {
          if (jet->pt > 30 && fabs(jet->eta) < 5.0 ) {
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
      //soft muons
      //******************************************************************************
      NSoftMuons = 0;
      for(Int_t j=0; j<muonArr->GetEntries(); j++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[j]);
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
             && (mu->typeBits & kTracker)
             && (mu->qualityBits & kTMLastStationAngTight)
             && mu->nTkHits > 10
             && fabs(mu->d0) < 0.2
             && fabs(mu->dz) < 0.1
             && !isCleanMuon
             && (!(mu->pt > 20 && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10))           
          ) {
          NSoftMuons++;
        }
      }

      //******************************************************************************
      //Pt Preselection
      //******************************************************************************
      if (!(leptonPt[0] > 20.0 && leptonPt[1] > 10.0)) continue;

      if (printDebug) cout << "pass preselection \n";

      for(int i = 0; i < leptonPt.size(); ++i) {
        for(int j = i+1; j < leptonPt.size(); ++j) {

          //require opposite sign
          if ((ChargeSelection == 0 && leptonCharge[i] == leptonCharge[j]) || (ChargeSelection == 1 && leptonCharge[0] != leptonCharge[j])) continue;

          if (printDebug) cout << "pass charge selection \n";

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

          Int_t RealLeptonIndex = -1;
          if (leptonType[i] == 11) {
            const mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);
            if ( isRealLepton(electron->isMCReal) ) RealLeptonIndex = i;
          } else {
            const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);
            if ( isRealLepton(muon->isMCReal) ) RealLeptonIndex = i;
          }
          if (leptonType[j] == 11) {
            const mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[leptonIndex[j]]);
            if ( isRealLepton(electron->isMCReal) ) RealLeptonIndex = j;
          } else {
            const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[leptonIndex[j]]);
            if ( isRealLepton(muon->isMCReal) ) RealLeptonIndex = j;
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
          double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt()*pfMet.Pt()*(1.0 - cos(deltaPhiDileptonMet)));
            
          //mtHiggs assuming neutrino pt == dilepton pt
          double enell = TMath::Sqrt(dilepton.Pt()*dilepton.Pt() + dilepton.M()*dilepton.M());
          double enenn = TMath::Sqrt(pfMet.Pt()* pfMet.Pt()  + dilepton.M()*dilepton.M());
          double enex  = dilepton.Px() + pfMet.Px();
          double eney  = dilepton.Py() + pfMet.Py();
          double mll   = dilepton.M();
          double mnu   = mll;
          double mtHiggs_NuPtEqualDileptonPt = TMath::Sqrt(mll*mll + mnu*mnu + 2.0*(enell*enenn - enex*enex - eney*eney));




//             if (!(finalState == 1)) continue;
//               if (!(finalState == 0 || finalState == 3)) continue;
//              if (!(finalState == 1 || finalState == 2)) continue;
          //      if (!(lepton2.Pt() < 0 && lepton2.Pt() >= 0) continue;
//             if (!(lepton2.Pt() >= 20)) continue;
//                           if (!(lepton2.Pt() < 20 && lepton2.Pt() >= 10)) continue;
//              if (!(lepton2.Pt() < 15 && lepton2.Pt() >= 10)) continue;
//             if (!(lepton2.Pt() < 10 && lepton2.Pt() >= 5)) continue;
   
          //*********************************************************************************************
          //Define Cuts
          //*********************************************************************************************
          const int nCuts = 15;
          bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                                 false, false, false, false, false};
          assert(CutLabel.size() == nCuts+1);

          Bool_t PreselPtCut = kTRUE;
          if(!(lepton1.Pt() >  20.0 && lepton2.Pt() > 10.0)) PreselPtCut = kFALSE;

          if(PreselPtCut) passCut[0] = true;
          if(zDiffMax < 100000.0) passCut[1] = true;            
          if(pfMet.Pt()    > 20.0)               passCut[2] = true; 
  
          if(dilepton.M() > 12.0)            passCut[3] = true;
   
          if (finalState == 0 || finalState == 1){ // mumu/ee
            if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
            if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 35) passCut[5] = true;
          }
          else if(finalState == 2 ||finalState == 3 ) { // emu
            passCut[4] = true;
            if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 20) passCut[5] = true;
          }

          if(NJets     < 1)              passCut[6] = true; 
        

          if (NSoftMuons == 0 )      passCut[7] = true;

          if (!(leptonPt.size() >= 3 && leptonPt[2] > 10.0)) passCut[8] = true;

          if(maxBtag < 2.1)                     passCut[9] = true;

          if (lepton1.Pt() > fPtMaxLowerCut) passCut[10] = true;
          if (lepton2.Pt() > fPtMinLowerCut) passCut[11] = true;
          if (dilepton.M() < fDileptonMassUpperCut)   passCut[12] = true;
          if (deltaPhiLeptons < fDeltaPhiCut) passCut[13] = true;
          if (mtHiggs >= fMtLowerCut && mtHiggs <= fMtUpperCut ) passCut[14] = true;

          if (printDebug) {
            cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
            cout << "PV: " << info->pvx << " " << info->pvy << " " << info->pvz << endl;
            cout << lepton1.Pt() << " " << lepton1.Eta() << " " << lepton1.Phi() << " " << leptonType[i] << endl;
            cout << lepton2.Pt() << " " << lepton2.Eta() << " " << lepton2.Phi() << " " << leptonType[j] << endl;
            cout << zDiffMax << " " << pfMet.Pt() << " " << dilepton.M() << " " 
                 << fabs(dilepton.M()-91.1876) << " " << PFMETdeltaPhilEt << "," <<  PFTrackMETdeltaPhilEt << " " 
                 << NJets << " " << NSoftMuons << " " << leptonPt.size() 
                 << " " << maxBtag << endl;
            cout << pfMet.Pt() << " " << pfMet.Px() << " " << pfMet.Py() << "  : " << pfTrackMet.Pt() << " "  << pfTrackMet.Px() << " " << pfTrackMet.Py() << endl;
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
              
            cout << "Event Selection: " << passCut[0] << " " << passCut[1] << " " << passCut[2] << " " << passCut[3] << " " << passCut[4] << " " << passCut[5] << " " << passCut[6] << " " << passCut[7] << " " << passCut[8] << " " << passCut[9] << " : " << passCut[10] << " " << passCut[11] << " "<< passCut[12] << " "<< passCut[13] << " "<< passCut[14] << " : "  ;

            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
              cout << "PassWWEvent " << info->runNum << " " << info->lumiSec << " " << info->evtNum  ;
            }
            cout << endl;
          }
            

          //*********************************************************************************************
          //Classify Fake Electrons
          //*********************************************************************************************  
          Double_t ptCategory = lepton2.Pt();
          Double_t etaCategory = lepton2.Eta();
          Bool_t isFakeEle = kFALSE;
          if (RealLeptonIndex == i) {
            ptCategory = leptonPt[j];
            etaCategory = leptonEta[j];
            if (leptonType[j] == 11) isFakeEle = kTRUE;
          } else {
            ptCategory = leptonPt[i];  
            etaCategory = leptonEta[i];
            if (leptonType[i] == 11) isFakeEle = kTRUE;
          }

          if (printDebug) cout << "Fake Ele: " << isFakeEle << " : " << RealLeptonIndex << " " << i << " " << j << " : " << leptonType[i] << " " << leptonType[j] << " " << endl;

          if (!isFakeEle) continue;
          if (printDebug) cout << "is Fake ele event\n";

          //*********************************************************************************************
          //Make Selection Histograms. Number of events passing each level of cut
          //*********************************************************************************************  
          bool passAllCuts = true;
          for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
          //Cut Selection Histograms
          fHWWSelection_fullSim->Fill(-1);
          if (finalState == 1 ) {
            fHWWToMuMuSelection_fullSim->Fill(-1);
          } else if(finalState == 0 ) {
            fHWWToEESelection_fullSim->Fill(-1);
          } else { //if(finalState == 2 || finalState == 3 ) {
            fHWWToEMuSelection_fullSim->Fill(-1);
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
              fHWWSelection_fullSim->Fill(k);
              if (finalState == 1 ) {
                fHWWToMuMuSelection_fullSim->Fill(k);
              } else if(finalState == 0 ) {
                fHWWToEESelection_fullSim->Fill(k);
              } else if(finalState == 2 || finalState == 3 ) {
                fHWWToEMuSelection_fullSim->Fill(k);
              }
            }
          }

          //*****************************************************************************************
          //Print EventList
          //*****************************************************************************************
          if (
            passCut[0] && passCut[1] 
            && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//           && passCut[10] && passCut[11] && passCut[12] && passCut[13] && passCut[14]
            ) {

            if ( (0==0) 
//              && finalState == 0
              ) {
//             passpassFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << ptCategory << " " << etaCategory << endl;
              passpassFile << info->evtNum << " : " << ptCategory << " " << etaCategory << endl;
            }
          }


          //*****************************************************************************************
          //Make Preselection Histograms  
          //*****************************************************************************************
          if (passCut[0] && passCut[1] && passCut[2] && passCut[3]) {
            if (RealLeptonIndex == j) {
              fLeptonFakePt_preSelection_fullSim->Fill(leptonPt[i]);
              fLeptonFakeEta_preSelection_fullSim->Fill(leptonEta[i]);
            } else {
              fLeptonFakePt_preSelection_fullSim->Fill(leptonPt[j]);
              fLeptonFakeEta_preSelection_fullSim->Fill(leptonEta[j]);
            }
          }

          //*****************************************************************************************
          //Make Preselection Histograms  
          //*****************************************************************************************
//         if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[5]) {
          if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
//         if (passAllCuts) {

            if (RealLeptonIndex == j) {
              fLeptonFakePt_WWSelection_fullSim->Fill(leptonPt[i]);
              fLeptonFakeEta_WWSelection_fullSim->Fill(leptonEta[i]);
              fLeptonFakePtEta_WWSelection_fullSim->Fill(leptonPt[i], leptonEta[i]);
            } else {
              fLeptonFakePt_WWSelection_fullSim->Fill(leptonPt[j]);
              fLeptonFakeEta_WWSelection_fullSim->Fill(leptonEta[j]);
              fLeptonFakePtEta_WWSelection_fullSim->Fill(leptonPt[j], leptonEta[j]);  
            }

            fLeptonPtMax_WWSelection_fullSim->Fill(lepton1.Pt());
            fLeptonPtMin_WWSelection_fullSim->Fill(lepton2.Pt());
            fMetPtHist_WWSelection_fullSim->Fill(pfMet.Pt());                             
            fDeltaPhiLeptons_WWSelection_fullSim->Fill(deltaPhiLeptons);

            dileptonMass_WWSelection_fullSim->Fill(dilepton.M());
            if (finalState == 0) {
              dileptonMass_ee_WWSelection_fullSim->Fill(dilepton.M());
            } else if (finalState == 1) {
              dileptonMass_mumu_WWSelection_fullSim->Fill(dilepton.M());
            } else {
              dileptonMass_emu_WWSelection_fullSim->Fill(dilepton.M());
            }

            fMtHiggs_afterCuts_WWSelection_fullSim->Fill(mtHiggs);

          }

          if (
            passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//              passAllCuts           
            ) {


            //  if (passAllCuts ) {
            if (finalState == 0) {
              Count_ee_WWSelection_fullSim += 1;
              Count_ee_statError_WWSelection_fullSim += 1*1;
            } else if (finalState == 1) {
              Count_mm_WWSelection_fullSim += 1;
              Count_mm_statError_WWSelection_fullSim += 1*1;
            } else if (finalState == 2) {
              Count_em_WWSelection_fullSim += 1;
              Count_em_statError_WWSelection_fullSim += 1*1;
            } else {
              Count_me_WWSelection_fullSim += 1;
              Count_me_statError_WWSelection_fullSim += 1*1;
            }


            if (ptCategory > 30.0) {
              if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
                Count_Pt30ToInf_EtaBin0_WWSelection_fullSim += 1;
                Count_Pt30ToInf_EtaBin0_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_Pt30ToInf_EtaBin1_WWSelection_fullSim += 1;
                Count_Pt30ToInf_EtaBin1_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_Pt30ToInf_EtaBin2_WWSelection_fullSim += 1;
                Count_Pt30ToInf_EtaBin2_statError_WWSelection_fullSim += 1*1;
              }
            } else if (ptCategory > 25.0 && ptCategory < 30.0 ) {
              if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
                Count_Pt25To30_EtaBin0_WWSelection_fullSim += 1;
                Count_Pt25To30_EtaBin0_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_Pt25To30_EtaBin1_WWSelection_fullSim += 1;
                Count_Pt25To30_EtaBin1_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_Pt25To30_EtaBin2_WWSelection_fullSim += 1;
                Count_Pt25To30_EtaBin2_statError_WWSelection_fullSim += 1*1;
              }
            } else if (ptCategory > 20.0 && ptCategory < 25.0 ) {
              if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
                Count_Pt20To25_EtaBin0_WWSelection_fullSim += 1;
                Count_Pt20To25_EtaBin0_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_Pt20To25_EtaBin1_WWSelection_fullSim += 1;
                Count_Pt20To25_EtaBin1_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_Pt20To25_EtaBin2_WWSelection_fullSim += 1;
                Count_Pt20To25_EtaBin2_statError_WWSelection_fullSim += 1*1;
              }
            } else if (ptCategory > 15.0 && ptCategory < 20.0 ) {
              if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
                Count_Pt15To20_EtaBin0_WWSelection_fullSim += 1;
                Count_Pt15To20_EtaBin0_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_Pt15To20_EtaBin1_WWSelection_fullSim += 1;
                Count_Pt15To20_EtaBin1_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_Pt15To20_EtaBin2_WWSelection_fullSim += 1;
                Count_Pt15To20_EtaBin2_statError_WWSelection_fullSim += 1*1;
              }
            } else if (ptCategory > 10.0 && ptCategory < 15.0 ) {
              if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
                Count_Pt10To15_EtaBin0_WWSelection_fullSim += 1;
                Count_Pt10To15_EtaBin0_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
                Count_Pt10To15_EtaBin1_WWSelection_fullSim += 1;
                Count_Pt10To15_EtaBin1_statError_WWSelection_fullSim += 1*1;
              } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
                Count_Pt10To15_EtaBin2_WWSelection_fullSim += 1;
                Count_Pt10To15_EtaBin2_statError_WWSelection_fullSim += 1*1;
              }
            } 


            if (leptonType[0] == 11) {
              if (ptCategory > 20) {
                if (fabs(etaCategory) < 1.479) {
                  Count_eF_Pt20ToInf_Barrel_WWSelection_fullSim += 1;
                  Count_eF_Pt20ToInf_Barrel_statError_WWSelection_fullSim += 1*1;
                } else {
                  Count_eF_Pt20ToInf_Endcap_WWSelection_fullSim += 1;
                  Count_eF_Pt20ToInf_Endcap_statError_WWSelection_fullSim += 1*1;
                }
              } else {
                if (fabs(etaCategory) < 1.479) {
                  Count_eF_Pt10To20_Barrel_WWSelection_fullSim += 1;
                  Count_eF_Pt10To20_Barrel_statError_WWSelection_fullSim += 1*1;
                } else {
                  Count_eF_Pt10To20_Endcap_WWSelection_fullSim += 1;
                  Count_eF_Pt10To20_Endcap_statError_WWSelection_fullSim += 1*1;
                }
              }
            } else {
              if (ptCategory > 20) {
                if (fabs(etaCategory) < 1.479) {
                  Count_mF_Pt20ToInf_Barrel_WWSelection_fullSim += 1;
                  Count_mF_Pt20ToInf_Barrel_statError_WWSelection_fullSim += 1*1;
                } else {
                  Count_mF_Pt20ToInf_Endcap_WWSelection_fullSim += 1;
                  Count_mF_Pt20ToInf_Endcap_statError_WWSelection_fullSim += 1*1;
                }
              } else {
                if (fabs(etaCategory) < 1.479) {
                  Count_mF_Pt10To20_Barrel_WWSelection_fullSim += 1;
                  Count_mF_Pt10To20_Barrel_statError_WWSelection_fullSim += 1*1;
                } else {
                  Count_mF_Pt10To20_Endcap_WWSelection_fullSim += 1;
                  Count_mF_Pt10To20_Endcap_statError_WWSelection_fullSim += 1*1;
                }
              }
            }

          }

          //*********************************************************************************************
          //Plots after all Cuts
          //*********************************************************************************************
          if (passAllCuts) {
//           fMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim->Fill(minDeltaPhiMetLepton);
//           fMtLepton1_afterCuts_WWSelection_fullSim->Fill(mTW[0]);
//           fMtLepton2_afterCuts_WWSelection_fullSim->Fill(mTW[1]);
//           fLeptonPtPlusMet_afterCuts_WWSelection_fullSim->Fill(lepton1.Pt()+lepton2.Pt()+pfMet.Pt());
    

          }



          
        } //loop over leptons
      }//loop over leptons
        
  


    }




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
  // Summary print out
  //============================================================================================================== 
  cout << "*******************************************************************************************************\n";
  cout << "\nPredicted Yields\n";
  cout << "*******************************************************************************************************\n";


  Double_t scaleFactor = 1000 / 35.5;
  scaleFactor=1;
  for (int i=1; i < fHWWToEESelection_predicted->GetXaxis()->GetNbins()+1; ++i) {
    cout << CutLabel[i] << " : " << fHWWToEESelection_predicted->GetBinContent(i) << "+/-" << fHWWToEESelection_predicted->GetBinError(i)  << " " << fHWWToMuMuSelection_predicted->GetBinContent(i) << "+/-" << fHWWToMuMuSelection_predicted->GetBinError(i) << " " << fHWWToEMuSelection_predicted->GetBinContent(i) << "+/-" << fHWWToEMuSelection_predicted->GetBinError(i)  << " " << fHWWSelection_predicted->GetBinContent(i) << "+/-" << fHWWSelection_predicted->GetBinError(i)  << endl;
  }

 
  cout << "**************************************************************\n";
  cout << "Event Count : ee final state\n";
  cout << "BB :" << Count_ee_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_ee_statError_WWSelection_predicted) <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : mm final state\n";
  cout << "BB :" << Count_mm_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_mm_statError_WWSelection_predicted) << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : em final state\n";
  cout << "BB :" << Count_em_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_em_statError_WWSelection_predicted) << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : me final state\n";
  cout << "BB :" << Count_me_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_me_statError_WWSelection_predicted) << endl;
  cout << "**************************************************************\n";
  cout << endl;

  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin0 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin0_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin0_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin1 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin1_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin1_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin2 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin2_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin2_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin0 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin0_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin0_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin1 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin1_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin1_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin2 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin2_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin2_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin0 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin0_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin0_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin1 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin1_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin1_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin2 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin2_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin2_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin0 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin0_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin0_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin1 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin1_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin1_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin2 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin2_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin2_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin0 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin0_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin0_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin1 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin1_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin1_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin2 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin2_WWSelection_predicted*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin2_statError_WWSelection_predicted)*scaleFactor <<  endl;
  cout << "**************************************************************\n";



  cout << "**************************************************************\n";
  cout << "Event Count : Ele + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_eF_Pt10To20_Barrel_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Barrel_statError_WWSelection_predicted) <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_eF_Pt10To20_Endcap_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Endcap_statError_WWSelection_predicted) <<  endl;
  cout << "Barrel Pt [20+] : " << Count_eF_Pt20ToInf_Barrel_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Barrel_statError_WWSelection_predicted) <<  endl;
  cout << "Endcap Pt [20+] : " << Count_eF_Pt20ToInf_Endcap_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Endcap_statError_WWSelection_predicted) <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Mu + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_mF_Pt10To20_Barrel_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Barrel_statError_WWSelection_predicted) <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_mF_Pt10To20_Endcap_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Endcap_statError_WWSelection_predicted) <<  endl;
  cout << "Barrel Pt [20+] : " << Count_mF_Pt20ToInf_Barrel_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Barrel_statError_WWSelection_predicted) <<  endl;
  cout << "Endcap Pt [20+] : " << Count_mF_Pt20ToInf_Endcap_WWSelection_predicted << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Endcap_statError_WWSelection_predicted) <<  endl;
  cout << "**************************************************************\n";
  cout << endl;


  cout << "*******************************************************************************************************\n";
  cout << "\nFullSim Yields\n";
  cout << "*******************************************************************************************************\n";

  for (int i=1; i < fHWWToEESelection_fullSim->GetXaxis()->GetNbins()+1; ++i) {
    cout << CutLabel[i] << " : " << fHWWToEESelection_fullSim->GetBinContent(i) << "+/-" << fHWWToEESelection_fullSim->GetBinError(i) << " " << fHWWToMuMuSelection_fullSim->GetBinContent(i) << "+/-" << fHWWToMuMuSelection_fullSim->GetBinError(i) << " " << fHWWToEMuSelection_fullSim->GetBinContent(i) << "+/-" << fHWWToEMuSelection_fullSim->GetBinError(i) << " " << fHWWSelection_fullSim->GetBinContent(i) << "+/-" << fHWWSelection_fullSim->GetBinError(i) << endl;
  }

 
  cout << "**************************************************************\n";
  cout << "Event Count : ee final state\n";
  cout << "BB :" << Count_ee_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_ee_statError_WWSelection_fullSim) <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : mm final state\n";
  cout << "BB :" << Count_mm_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_mm_statError_WWSelection_fullSim) << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : em final state\n";
  cout << "BB :" << Count_em_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_em_statError_WWSelection_fullSim) << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : me final state\n";
  cout << "BB :" << Count_me_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_me_statError_WWSelection_fullSim) << endl;
  cout << "**************************************************************\n";
  cout << endl;

  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin0 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin0_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin0_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin1 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin1_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin1_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt10To15 EtaBin2 \n";
  cout << "BB :" << Count_Pt10To15_EtaBin2_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt10To15_EtaBin2_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin0 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin0_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin0_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin1 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin1_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin1_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt15To20 EtaBin2 \n";
  cout << "BB :" << Count_Pt15To20_EtaBin2_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt15To20_EtaBin2_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin0 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin0_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin0_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin1 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin1_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin1_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20To25 EtaBin2 \n";
  cout << "BB :" << Count_Pt20To25_EtaBin2_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20To25_EtaBin2_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin0 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin0_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin0_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin1 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin1_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin1_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt25To30 EtaBin2 \n";
  cout << "BB :" << Count_Pt25To30_EtaBin2_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt25To30_EtaBin2_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";

  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin0 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin0_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin0_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin1 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin1_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin1_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt30ToInf EtaBin2 \n";
  cout << "BB :" << Count_Pt30ToInf_EtaBin2_WWSelection_fullSim*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt30ToInf_EtaBin2_statError_WWSelection_fullSim)*scaleFactor <<  endl;
  cout << "**************************************************************\n";



  cout << "**************************************************************\n";
  cout << "Event Count : Ele + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_eF_Pt10To20_Barrel_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Barrel_statError_WWSelection_fullSim) <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_eF_Pt10To20_Endcap_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_eF_Pt10To20_Endcap_statError_WWSelection_fullSim) <<  endl;
  cout << "Barrel Pt [20+] : " << Count_eF_Pt20ToInf_Barrel_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Barrel_statError_WWSelection_fullSim) <<  endl;
  cout << "Endcap Pt [20+] : " << Count_eF_Pt20ToInf_Endcap_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_eF_Pt20ToInf_Endcap_statError_WWSelection_fullSim) <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Mu + Fake\n";
  cout << "Barrel Pt [10,20] : " << Count_mF_Pt10To20_Barrel_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Barrel_statError_WWSelection_fullSim) <<  endl;
  cout << "Endcap Pt [10,20] : " << Count_mF_Pt10To20_Endcap_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_mF_Pt10To20_Endcap_statError_WWSelection_fullSim) <<  endl;
  cout << "Barrel Pt [20+] : " << Count_mF_Pt20ToInf_Barrel_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Barrel_statError_WWSelection_fullSim) <<  endl;
  cout << "Endcap Pt [20+] : " << Count_mF_Pt20ToInf_Endcap_WWSelection_fullSim << " +/- " << TMath::Sqrt(Count_mF_Pt20ToInf_Endcap_statError_WWSelection_fullSim) <<  endl;
  cout << "**************************************************************\n";
  cout << endl;


  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("ElectronFakeBkgClosureTest.root", "RECREATE");
//   TFile *file = new TFile("WWSelectionPlots_IsoDenominator.root", "RECREATE");

  file->WriteTObject(fHWWSelection_predicted, fHWWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEESelection_predicted, fHWWToEESelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToMuMuSelection_predicted , fHWWToMuMuSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEMuSelection_predicted,fHWWToEMuSelection_predicted->GetName(), "WriteDelete");
  
  file->WriteTObject(fLeptonFakeEta_WWSelection_predicted ,fLeptonFakeEta_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePt_WWSelection_predicted ,fLeptonFakePt_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePtEta_WWSelection_predicted ,fLeptonFakePtEta_WWSelection_predicted->GetName(), "WriteDelete");


  file->WriteTObject(fLeptonPtMax_WWSelection_predicted ,fLeptonPtMax_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtMin_WWSelection_predicted ,fLeptonPtMin_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fMetPtHist_WWSelection_predicted ,fMetPtHist_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaPhiLeptons_WWSelection_predicted ,fDeltaPhiLeptons_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaEtaLeptons_WWSelection_predicted ,fDeltaEtaLeptons_WWSelection_predicted->GetName(), "WriteDelete");
                     
  file->WriteTObject(fMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted ,fMinDeltaPhiLeptonMet_afterCuts_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton1_afterCuts_WWSelection_predicted ,fMtLepton1_afterCuts_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton2_afterCuts_WWSelection_predicted ,fMtLepton2_afterCuts_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fMtHiggs_afterCuts_WWSelection_predicted ,fMtHiggs_afterCuts_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtPlusMet_afterCuts_WWSelection_predicted ,fLeptonPtPlusMet_afterCuts_WWSelection_predicted->GetName(), "WriteDelete");
   
  file->WriteTObject(dileptonMass_WWSelection_predicted ,dileptonMass_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_ee_WWSelection_predicted ,dileptonMass_ee_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_emu_WWSelection_predicted ,dileptonMass_emu_WWSelection_predicted->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_mumu_WWSelection_predicted ,dileptonMass_mumu_WWSelection_predicted->GetName(), "WriteDelete");



  file->WriteTObject(fHWWSelection_fullSim, fHWWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEESelection_fullSim, fHWWToEESelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToMuMuSelection_fullSim , fHWWToMuMuSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEMuSelection_fullSim,fHWWToEMuSelection_fullSim->GetName(), "WriteDelete");
  
  file->WriteTObject(fLeptonFakeEta_WWSelection_fullSim ,fLeptonFakeEta_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePt_WWSelection_fullSim ,fLeptonFakePt_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonFakePtEta_WWSelection_fullSim ,fLeptonFakePtEta_WWSelection_fullSim->GetName(), "WriteDelete");


  file->WriteTObject(fLeptonPtMax_WWSelection_fullSim ,fLeptonPtMax_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtMin_WWSelection_fullSim ,fLeptonPtMin_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fMetPtHist_WWSelection_fullSim ,fMetPtHist_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaPhiLeptons_WWSelection_fullSim ,fDeltaPhiLeptons_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fDeltaEtaLeptons_WWSelection_fullSim ,fDeltaEtaLeptons_WWSelection_fullSim->GetName(), "WriteDelete");
                     
  file->WriteTObject(fMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim ,fMinDeltaPhiLeptonMet_afterCuts_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton1_afterCuts_WWSelection_fullSim ,fMtLepton1_afterCuts_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fMtLepton2_afterCuts_WWSelection_fullSim ,fMtLepton2_afterCuts_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fMtHiggs_afterCuts_WWSelection_fullSim ,fMtHiggs_afterCuts_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(fLeptonPtPlusMet_afterCuts_WWSelection_fullSim ,fLeptonPtPlusMet_afterCuts_WWSelection_fullSim->GetName(), "WriteDelete");
   
  file->WriteTObject(dileptonMass_WWSelection_fullSim ,dileptonMass_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_ee_WWSelection_fullSim ,dileptonMass_ee_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_emu_WWSelection_fullSim ,dileptonMass_emu_WWSelection_fullSim->GetName(), "WriteDelete");
  file->WriteTObject(dileptonMass_mumu_WWSelection_fullSim ,dileptonMass_mumu_WWSelection_fullSim->GetName(), "WriteDelete");




  file->Close();
  delete file;



  //*****************************************************************************************
  //Draw ClosureTest Plots
  //*****************************************************************************************
  file = new TFile("ElectronFakeBkgClosureTestPlots.root", "RECREATE");

  TCanvas *cv = new TCanvas("cv", "cv", 800,600);
  TLegend *legend = new TLegend(0.20,0.75,0.60,0.90);



  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonFakePt_preSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonFakePt_preSelection_predicted,"Predicted", "LP");
  fLeptonFakePt_preSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonFakePt_preSelection_fullSim->SetLineColor(kBlue);
  fLeptonFakePt_preSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonFakePt_preSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonFakePt_preSelection_fullSim->GetYaxis()->SetRangeUser(0,80);
  fLeptonFakePt_preSelection_fullSim->Draw("E1");
  fLeptonFakePt_preSelection_predicted->SetMarkerColor(kRed);
  fLeptonFakePt_preSelection_predicted->SetLineColor(kRed);
  fLeptonFakePt_preSelection_predicted->SetFillColor(kRed);
  fLeptonFakePt_preSelection_predicted->SetFillStyle(3001);
  fLeptonFakePt_preSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_FakeElePt_PreSelection.gif");
  file->WriteTObject(cv , "LeptonFakePt", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonFakeEta_preSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonFakeEta_preSelection_predicted,"Predicted", "LP");
  fLeptonFakeEta_preSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonFakeEta_preSelection_fullSim->SetLineColor(kBlue);
  fLeptonFakeEta_preSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonFakeEta_preSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonFakeEta_preSelection_fullSim->GetYaxis()->SetRangeUser(0,80);
  fLeptonFakeEta_preSelection_fullSim->Draw("E1");
  fLeptonFakeEta_preSelection_predicted->SetMarkerColor(kRed);
  fLeptonFakeEta_preSelection_predicted->SetLineColor(kRed);
  fLeptonFakeEta_preSelection_predicted->SetFillColor(kRed);
  fLeptonFakeEta_preSelection_predicted->SetFillStyle(3001);
  fLeptonFakeEta_preSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_FakeEleEta_PreSelection.gif");
  file->WriteTObject(cv , "LeptonFakeEta", "WriteDelete");








  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonFakePt_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonFakePt_WWSelection_predicted,"Predicted", "LP");
  fLeptonFakePt_WWSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonFakePt_WWSelection_fullSim->SetLineColor(kBlue);
  fLeptonFakePt_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonFakePt_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonFakePt_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,60);
  fLeptonFakePt_WWSelection_fullSim->Draw("E1");
  fLeptonFakePt_WWSelection_predicted->SetMarkerColor(kRed);
  fLeptonFakePt_WWSelection_predicted->SetLineColor(kRed);
  fLeptonFakePt_WWSelection_predicted->SetFillColor(kRed);
  fLeptonFakePt_WWSelection_predicted->SetFillStyle(3001);
  fLeptonFakePt_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_FakeElePt.gif");
  file->WriteTObject(cv , "LeptonFakePt", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonFakeEta_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonFakeEta_WWSelection_predicted,"Predicted", "LP");
  fLeptonFakeEta_WWSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonFakeEta_WWSelection_fullSim->SetLineColor(kBlue);
  fLeptonFakeEta_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonFakeEta_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonFakeEta_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,60);
  fLeptonFakeEta_WWSelection_fullSim->Draw("E1");
  fLeptonFakeEta_WWSelection_predicted->SetMarkerColor(kRed);
  fLeptonFakeEta_WWSelection_predicted->SetLineColor(kRed);
  fLeptonFakeEta_WWSelection_predicted->SetFillColor(kRed);
  fLeptonFakeEta_WWSelection_predicted->SetFillStyle(3001);
  fLeptonFakeEta_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_FakeEleEta.gif");
  file->WriteTObject(cv , "LeptonFakeEta", "WriteDelete");



  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonPtMax_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonPtMax_WWSelection_predicted,"Predicted", "LP");
  fLeptonPtMax_WWSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonPtMax_WWSelection_fullSim->SetLineColor(kBlue);
  fLeptonPtMax_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonPtMax_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonPtMax_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,80);
  fLeptonPtMax_WWSelection_fullSim->Draw("E1");
  fLeptonPtMax_WWSelection_predicted->SetMarkerColor(kRed);
  fLeptonPtMax_WWSelection_predicted->SetLineColor(kRed);
  fLeptonPtMax_WWSelection_predicted->SetFillColor(kRed);
  fLeptonPtMax_WWSelection_predicted->SetFillStyle(3001);
  fLeptonPtMax_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_PtMax.gif");
  file->WriteTObject(cv , "LeptonPtMax", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fLeptonPtMin_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fLeptonPtMin_WWSelection_predicted,"Predicted", "LP");
  fLeptonPtMin_WWSelection_fullSim->SetMarkerColor(kBlue);
  fLeptonPtMin_WWSelection_fullSim->SetLineColor(kBlue);
  fLeptonPtMin_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fLeptonPtMin_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fLeptonPtMin_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,80);
  fLeptonPtMin_WWSelection_fullSim->Draw("E1");
  fLeptonPtMin_WWSelection_predicted->SetMarkerColor(kRed);
  fLeptonPtMin_WWSelection_predicted->SetLineColor(kRed);
  fLeptonPtMin_WWSelection_predicted->SetFillColor(kRed);
  fLeptonPtMin_WWSelection_predicted->SetFillStyle(3001);
  fLeptonPtMin_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_PtMin.gif");
  file->WriteTObject(cv , "LeptonPtMin", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fMetPtHist_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fMetPtHist_WWSelection_predicted,"Predicted", "LP");
  fMetPtHist_WWSelection_fullSim->SetMarkerColor(kBlue);
  fMetPtHist_WWSelection_fullSim->SetLineColor(kBlue);
  fMetPtHist_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fMetPtHist_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fMetPtHist_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,60);
  fMetPtHist_WWSelection_fullSim->Draw("E1");
  fMetPtHist_WWSelection_predicted->SetMarkerColor(kRed);
  fMetPtHist_WWSelection_predicted->SetLineColor(kRed);
  fMetPtHist_WWSelection_predicted->SetFillColor(kRed);
  fMetPtHist_WWSelection_predicted->SetFillStyle(3001);
  fMetPtHist_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_Met.gif");
  file->WriteTObject(cv , "MetPtHist", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fDeltaPhiLeptons_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fDeltaPhiLeptons_WWSelection_predicted,"Predicted", "LP");
  fDeltaPhiLeptons_WWSelection_fullSim->SetMarkerColor(kBlue);
  fDeltaPhiLeptons_WWSelection_fullSim->SetLineColor(kBlue);
  fDeltaPhiLeptons_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fDeltaPhiLeptons_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fDeltaPhiLeptons_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,35);
  fDeltaPhiLeptons_WWSelection_fullSim->Draw("E1");
  fDeltaPhiLeptons_WWSelection_predicted->SetMarkerColor(kRed);
  fDeltaPhiLeptons_WWSelection_predicted->SetLineColor(kRed);
  fDeltaPhiLeptons_WWSelection_predicted->SetFillColor(kRed);
  fDeltaPhiLeptons_WWSelection_predicted->SetFillStyle(3001);
  fDeltaPhiLeptons_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_DeltaPhi.gif");
  file->WriteTObject(cv , "DeltaPhiLeptons", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(dileptonMass_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(dileptonMass_WWSelection_predicted,"Predicted", "LP");
  dileptonMass_WWSelection_fullSim->SetMarkerColor(kBlue);
  dileptonMass_WWSelection_fullSim->SetLineColor(kBlue);
  dileptonMass_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  dileptonMass_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  dileptonMass_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,30);
  dileptonMass_WWSelection_fullSim->Draw("E1");
  dileptonMass_WWSelection_predicted->SetMarkerColor(kRed);
  dileptonMass_WWSelection_predicted->SetLineColor(kRed);
  dileptonMass_WWSelection_predicted->SetFillColor(kRed);
  dileptonMass_WWSelection_predicted->SetFillStyle(3001);
  dileptonMass_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_DileptonMass.gif");
  file->WriteTObject(cv , "dileptonMass", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fMtHiggs_afterCuts_WWSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fMtHiggs_afterCuts_WWSelection_predicted,"Predicted", "LP");
  fMtHiggs_afterCuts_WWSelection_fullSim->SetMarkerColor(kBlue);
  fMtHiggs_afterCuts_WWSelection_fullSim->SetLineColor(kBlue);
  fMtHiggs_afterCuts_WWSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fMtHiggs_afterCuts_WWSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fMtHiggs_afterCuts_WWSelection_fullSim->GetYaxis()->SetRangeUser(0,30);
  fMtHiggs_afterCuts_WWSelection_fullSim->Draw("E1");
  fMtHiggs_afterCuts_WWSelection_predicted->SetMarkerColor(kRed);
  fMtHiggs_afterCuts_WWSelection_predicted->SetLineColor(kRed);
  fMtHiggs_afterCuts_WWSelection_predicted->SetFillColor(kRed);
  fMtHiggs_afterCuts_WWSelection_predicted->SetFillStyle(3001);
  fMtHiggs_afterCuts_WWSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_MtHiggs.gif");
  file->WriteTObject(cv , "mtHiggs", "WriteDelete");



  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fHWWToEESelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fHWWToEESelection_predicted,"Predicted", "LP");
  fHWWToEESelection_fullSim->SetMarkerColor(kBlue);
  fHWWToEESelection_fullSim->SetLineColor(kBlue);
  fHWWToEESelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fHWWToEESelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fHWWToEESelection_fullSim->GetYaxis()->SetRangeUser(0,250);
  fHWWToEESelection_fullSim->Draw("E1");
  fHWWToEESelection_predicted->SetMarkerColor(kRed);
  fHWWToEESelection_predicted->SetLineColor(kRed);
  fHWWToEESelection_predicted->SetFillColor(kRed);
  fHWWToEESelection_predicted->SetFillStyle(3001);
  fHWWToEESelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_CutFlowEE.gif");
  file->WriteTObject(cv , "fHWWToEESelection", "WriteDelete");


  legend = new TLegend(0.50,0.75,0.85,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(fHWWToEMuSelection_fullSim,"Simulation", "LP");
  legend->AddEntry(fHWWToEMuSelection_predicted,"Predicted", "LP");
  fHWWToEMuSelection_fullSim->SetMarkerColor(kBlue);
  fHWWToEMuSelection_fullSim->SetLineColor(kBlue);
  fHWWToEMuSelection_fullSim->GetYaxis()->SetTitleOffset(1.2);
  fHWWToEMuSelection_fullSim->GetXaxis()->SetTitleOffset(1.05);
  fHWWToEMuSelection_fullSim->GetYaxis()->SetRangeUser(0,300);
  fHWWToEMuSelection_fullSim->Draw("E1");
  fHWWToEMuSelection_predicted->SetMarkerColor(kRed);
  fHWWToEMuSelection_predicted->SetLineColor(kRed);
  fHWWToEMuSelection_predicted->SetFillColor(kRed);
  fHWWToEMuSelection_predicted->SetFillStyle(3001);
  fHWWToEMuSelection_predicted->Draw("E3,same");
  legend->Draw();
  cv->SaveAs("FakeElectronClosureTest_CutFlowEMu.gif");
  file->WriteTObject(cv , "fHWWToEMuSelection", "WriteDelete");

  
  file->Close();
  delete file;
     
  gBenchmark->Show("WWTemplate");       
} 





Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (!(fabs(ele->eta) <= 2.5)) pass = kFALSE;

  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->NeutralHadronIso007_10Threshold;
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
    +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
    -ele->NeutralHadronIso007_10Threshold;
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
      if (fabs(ele->scEta) > 1.0) {
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
  if (!(ele->pt > 10 && fabs(ele->eta) <= 2.5)) pass = kFALSE;

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
            && (ele->trkIso03) / ele->pt < 0.20
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

  return pass;
}



Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

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
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && ( mu->nPixHits > 0)
        && fabs(mu->d0) < 0.2
        && fabs(mu->dz) < 0.1
        && (mu->ChargedIso03 + mu->NeutralIso03_10Threshold) / mu->pt < 0.4
        && mu->trkIso03 / mu->pt < 0.3
        && mu->emIso03 / mu->pt < 0.3
        && mu->hadIso03 / mu->pt < 0.3
        && ( mu->pterr / mu->pt < 0.1)
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
