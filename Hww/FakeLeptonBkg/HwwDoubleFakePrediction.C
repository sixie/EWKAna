//root -l EWKAna/Hww/FakeRate/HwwDoubleFakePrediction.C+\(\"/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_full-d22_TwoV3DenominatorWithTriggerSkim.root\",130,\"\"\)
//root -l EWKAna/Hww/FakeRate/HwwDoubleFakePrediction.C+\(\"/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_full-d22_TwoRecoLeptonWithTriggerSkim.root\",130,\"\"\)
//root -l EWKAna/Hww/FakeRate/HwwDoubleFakePrediction.C+\(\"/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-v8-pu11_noskim_normalized.root\",130,\"\"\)
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
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
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



//=== MAIN MACRO =================================================================================================

void HwwDoubleFakePrediction(const string inputFilename, Double_t mHiggs, 
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

  Double_t fHiggsMass = mHiggs;

  //Configure Higgs Mass Dependant Cuts
  if (fHiggsMass == 120) { fPtMaxLowerCut = 20.0;        fPtMinLowerCut = 20.0; 
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

  //--------------------------------------------------------------------------------------------------------------
  // Set up Fake Rate
  //==============================================================================================================
  Bool_t use2DFakeRate = kTRUE;
  Bool_t useFitFunction = kFALSE;
  mithep::FakeRate *fFakeRate = new mithep::FakeRate(
//                                                         "ElectronFakeRate.SmurfV3.root",
                                                            "ElectronFakeRate.SmurfV4M.root",
//                                                         "ElectronFakeRate.SmurfV4VT.root",
                                                     "MuonFakeRate_V4.root",
                                                     "", "",
                                                     "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLSample_ptThreshold15_PtEta",
                                                     "MuonFakeRate_v4_Mu9Jet30_PtEta",
                                                     use2DFakeRate, useFitFunction );

  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *fHWWSelection= new TH1D("hHWWSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection= new TH1D("hHWWToEESelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection= new TH1D("hHWWToMuMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection= new TH1D("hHWWToEMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWSelection_sysError= new TH1D("hHWWSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5); 
  TH1D *fHWWToEESelection_sysError= new TH1D("hHWWToEESelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection_sysError= new TH1D("hHWWToMuMuSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection_sysError= new TH1D("hHWWToEMuSelection_sysError", ";Cut Number;Number of Events", 15, -1.5, 13.5);

  TH1D *fLeptonEta = new TH1D(          "hLeptonEta",";LeptonEta;Number of Events",20,-5.0,5.0);
  TH1D *fLeptonPtMax = new TH1D(        "hLeptonPtMax",";Lepton P_t Max;Number of Events",20,0.,150.);
  TH1D *fLeptonPtMin = new TH1D(        "hLeptonPtMin",";Lepton P_t Min;Number of Events",20,0.,150.);
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
  Double_t Count_Pt20ToInf_EtaBin0 = 0;
  Double_t Count_Pt20ToInf_EtaBin1 = 0;
  Double_t Count_Pt20ToInf_EtaBin2 = 0;
  Double_t Count_Pt20ToInf_EtaBin0_statError = 0;
  Double_t Count_Pt20ToInf_EtaBin1_statError = 0;
  Double_t Count_Pt20ToInf_EtaBin2_statError = 0;
  Double_t Count_Pt20ToInf_EtaBin0_sysError = 0;
  Double_t Count_Pt20ToInf_EtaBin1_sysError = 0;
  Double_t Count_Pt20ToInf_EtaBin2_sysError = 0;





  ofstream eventListFile("eventList.txt");

  Double_t NEventsPreSel = 0;
  Double_t NEventsPreSelPassHiggsCuts = 0;


  fHWWSelection->Sumw2();
  fHWWToEESelection->Sumw2();
  fHWWToMuMuSelection->Sumw2();
  fHWWToEMuSelection->Sumw2();
  fLeptonEta->Sumw2();
  fLeptonPtMax->Sumw2();
  fLeptonPtMin->Sumw2();
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
  rlrm.AddJSONFile("Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON_v3.txt"); 
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

    Int_t NSoftMuons = 0;
    Int_t NLeptons = 0;
    vector<Int_t> leptonType;
    vector<Int_t> leptonIndex;
    vector<Double_t> leptonPt;
    vector<Double_t> leptonEta;
    vector<Double_t> leptonPhi;
    vector<Int_t> leptonCharge;

    const mithep::TJet *leadingJet = 0;
    
    Int_t NTightLeptons = 0;

    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
      if ( (0==0)
           &&
           passMuonDenominatorCuts(mu)
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
        if (passMuonCuts(mu)) NTightLeptons++;
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
           passElectronDenominatorCuts(ele)
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
        if (passElectronCuts(ele)) NTightLeptons++;
      }
      
      if (info->runNum == 148822 && info->evtNum ==  269137037) {
        cout << "Ele " << ele->pt << " " << ele->eta << " " << ele->phi << " : " << passElectronCuts(ele) << endl;
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
    //select events with 0 tight leptons
    //******************************************************************************
    if (NTightLeptons > 0) continue;

    for (Int_t i=0; i<leptonPt.size(); ++i) {

      for(Int_t j=i; j<leptonPt.size(); ++j) {
        
        if (!(leptonPt[i] > 10.0)) continue;
        
        //preselection
        if ((ChargeSelection == 0 && leptonCharge[i] == leptonCharge[j]) || (ChargeSelection == 1 && leptonCharge[i] != leptonCharge[j])) continue;

        if (leptonType[i] == 11) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[i]]);
          if (passElectronCuts(ele)) continue;
        } else {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[i]]);      
          if (passMuonCuts(mu)) continue;
        }
        if (leptonType[j] == 11) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[j]]);
          if (passElectronCuts(ele)) continue;
        } else {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[j]]);      
          if (passMuonCuts(mu)) continue;
        }


//             if (! (leptonPt[j] > 20)) continue;
//            if (! (leptonPt[j] > 10 && leptonPt[j] < 20)) continue;
//        if (! (fabs(leptonEta[j]) > 0 && fabs(leptonEta[j]) < 1.0)) continue;
//       if (! (fabs(leptonEta[j]) > 1.0 && fabs(leptonEta[j]) < 1.5)) continue;
//       if (! (fabs(leptonEta[j]) > 1.5 && fabs(leptonEta[j]) < 2.5)) continue;

        if (leptonType[i] == 11) {
          if (leptonPt[i] < 15) continue;
        } 
        if (leptonType[j] == 11) {
          if (leptonPt[j] < 15) continue;
        } 
      
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
        

      //******************************************************************************
      //Calculate Fake Rate
      //******************************************************************************
      double fakeRate1;
      double fakeRate1ErrorLow ;
      double fakeRate1ErrorHigh ;
      double fakeRate2;
      double fakeRate2ErrorLow ;
      double fakeRate2ErrorHigh ;

      if (leptonType[i] == 11) {
        fakeRate1 = fFakeRate->ElectronFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i]) / (1-fFakeRate->ElectronFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i]));
        fakeRate1ErrorLow = fFakeRate->ElectronFakeRateStatErrorLow(leptonPt[i], leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->ElectronFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i])),2);
        fakeRate1ErrorHigh = fFakeRate->ElectronFakeRateStatErrorHigh(leptonPt[i], leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->ElectronFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i])),2);        
      } else {
        if (leptonPt[i] < 35) {
          fakeRate1 = fFakeRate->MuonFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i]) / (1-fFakeRate->MuonFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i]));
          fakeRate1ErrorLow = fFakeRate->MuonFakeRateStatErrorLow(leptonPt[i], leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->MuonFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i])),2);
          fakeRate1ErrorHigh = fFakeRate->MuonFakeRateStatErrorHigh(leptonPt[i], leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->MuonFakeRate(leptonPt[i], leptonEta[i], leptonPhi[i])),2);
        } else {
          fakeRate1 = fFakeRate->MuonFakeRate(30.0, leptonEta[i], leptonPhi[i]) / (1-fFakeRate->MuonFakeRate(30.0, leptonEta[i], leptonPhi[i]));
          fakeRate1ErrorLow = fFakeRate->MuonFakeRateStatErrorLow(30.0, leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->MuonFakeRate(30.0, leptonEta[i], leptonPhi[i])),2);
          fakeRate1ErrorHigh = fFakeRate->MuonFakeRateStatErrorHigh(30.0, leptonEta[i], leptonPhi[i]) / pow((1- fFakeRate->MuonFakeRate(30.0, leptonEta[i], leptonPhi[i])),2);
        }
      }
      
      if (leptonType[j] == 11) {
        fakeRate2 = fFakeRate->ElectronFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j]) / (1-fFakeRate->ElectronFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j]));
        fakeRate2ErrorLow = fFakeRate->ElectronFakeRateStatErrorLow(leptonPt[j], leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->ElectronFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j])),2);
        fakeRate2ErrorHigh = fFakeRate->ElectronFakeRateStatErrorHigh(leptonPt[j], leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->ElectronFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j])),2);        
      } else {
        if (leptonPt[j] < 35) {
          fakeRate2 = fFakeRate->MuonFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j]) / (1-fFakeRate->MuonFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j]));
          fakeRate2ErrorLow = fFakeRate->MuonFakeRateStatErrorLow(leptonPt[j], leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->MuonFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j])),2);
          fakeRate2ErrorHigh = fFakeRate->MuonFakeRateStatErrorHigh(leptonPt[j], leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->MuonFakeRate(leptonPt[j], leptonEta[j], leptonPhi[j])),2);
        } else {
          fakeRate2 = fFakeRate->MuonFakeRate(30.0, leptonEta[j], leptonPhi[j]) / (1-fFakeRate->MuonFakeRate(30.0, leptonEta[j], leptonPhi[j]));
          fakeRate2ErrorLow = fFakeRate->MuonFakeRateStatErrorLow(30.0, leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->MuonFakeRate(30.0, leptonEta[j], leptonPhi[j])),2);
          fakeRate2ErrorHigh = fFakeRate->MuonFakeRateStatErrorHigh(30.0, leptonEta[j], leptonPhi[j]) / pow((1- fFakeRate->MuonFakeRate(30.0, leptonEta[j], leptonPhi[j])),2);
        }
      }

       
       Double_t eventWeight = fakeRate1 * fakeRate2; 
       //For MC
       //eventWeight = 35.5 * info->eventweight * fakeRate1 * fakeRate2; 
       //eventWeight = 1.0; 
      Double_t eventWeightErrorLow = eventWeight * TMath::Sqrt(pow(fakeRate1ErrorLow/fakeRate1,2) + pow(fakeRate2ErrorLow/fakeRate2,2));
      Double_t eventWeightErrorHigh = eventWeight * TMath::Sqrt(pow(fakeRate1ErrorHigh/fakeRate1,2) + pow(fakeRate2ErrorHigh/fakeRate2,2));
       Double_t eventWeightError = ( eventWeightErrorLow +   eventWeightErrorHigh) / 2;
       // eventWeightError = 0;
  
      //******************************************************************************
      //Count Jets
      //******************************************************************************
      double maxBtag = -99999;
      Int_t NJets = 0;
      for(Int_t k=0; k<jetArr->GetEntries(); k++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[k]);
        
        Bool_t leptonOverlap = kFALSE;
        if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[i],leptonEta[i]) < 0.3 ||
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[j],leptonEta[j]) < 0.3 ) {
          leptonOverlap = kTRUE;
        }
        
        if (!leptonOverlap) {
          if (jet->pt > 25 && fabs(jet->eta) < 5.0 ) {
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
      for(Int_t m=0; m<muonArr->GetEntries(); m++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[m]);
        Bool_t isOverlap = kFALSE;        
        if ( leptonType[i] == 13 
             && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[i],leptonEta[i]) < 0.1
          ) {
          isOverlap = kTRUE;         
        }
        if ( leptonType[j] == 13 
             && mithep::MathUtils::DeltaR(mu->phi, mu->eta, leptonPhi[j],leptonEta[j]) < 0.1
          ) {
          isOverlap = kTRUE;         
        }
       
        if ( fabs(mu->eta) < 2.4                   
             && mu->pt > 3.0
             && (mu->qualityBits & kTMLastStationAngTight)
             && mu->nTkHits > 10
             && fabs(mu->d0) < 0.2
             && (mu->typeBits & kTracker)
             && !isOverlap
          ) {
          NSoftMuons++;
        }
      }
      

      //***********************************************************************************************
      //|Z_vert-Z_l| maximum
      //***********************************************************************************************
      double zDiffMax = 0.0;
      
        double dz_i = 0;
        if (leptonType[i] == 11) {
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
        
        zDiffMax = fabs(dz_i - dz_j);



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

      //*********************************************************************************************
      //Define Cuts
      //*********************************************************************************************
      const int nCuts = 14;
      bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false};
      
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
      
      passCut[8] = true;
      
      if(maxBtag < 2.1)          passCut[9] = true;
      if (lepton1.Pt() > fPtMaxLowerCut) passCut[10] = true;
      if (lepton2.Pt() > fPtMinLowerCut) passCut[11] = true;
      if (dilepton.M() < fDileptonMassUpperCut)   passCut[12] = true;
      if (deltaPhiLeptons < fDeltaPhiCut) passCut[13] = true;

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
//         if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[5]) {
         if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
//         if (passAllCuts) {

          NEventsPreSel += eventWeight;
          if (passCut[10] && passCut[11] && passCut[12] && passCut[13]) NEventsPreSelPassHiggsCuts += eventWeight;

//           fLeptonEta->Fill(lepton1.Eta(), eventWeight); 
          fLeptonEta->Fill(lepton2.Eta(), eventWeight);
          fLeptonPtMax->Fill(lepton1.Pt(), eventWeight);
          fLeptonPtMin->Fill(lepton2.Pt(), eventWeight);
          fMetPtHist->Fill(met.Pt(), eventWeight);                             
          fMetPhiHist->Fill(met.Phi(), eventWeight);                            
          fDeltaPhiLeptons->Fill(deltaPhiLeptons, eventWeight);

//           fLeptonEta_sysError->Fill(lepton2.Eta(), 2*eventWeightError);
          fLeptonEta_sysError->Fill(lepton2.Eta(), eventWeightError);
          fLeptonPtMin_sysError->Fill(lepton2.Pt(), eventWeightError);
          fLeptonPtMax_sysError->Fill(lepton1.Pt(), eventWeightError);
          fMetPtHist_sysError->Fill(met.Pt(), eventWeightError);                             
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
          passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]
//           passAllCuts 
          ) {
          cout << info->runNum << " " << info->evtNum << " " << info->lumiSec << " : Fake rate: " << fakeRate1 << " * " << fakeRate2 << " = " << eventWeight << " " << finalState << endl;
          cout << leptonType[i] << " " << leptonPt[i] << " " << leptonEta[i] << " : " << fakeRate1 << " - " << fakeRate1ErrorLow << " + " << fakeRate1ErrorHigh << "   :    " 
               << leptonType[j] << " " << leptonPt[j] << " " << leptonEta[j] << " : " << fakeRate2 << " - " << fakeRate2ErrorLow << " + " << fakeRate2ErrorHigh << "   :    " 
               << eventWeight << " - " <<  eventWeightErrorLow << " + " << eventWeightErrorHigh << " " 
               << endl;
          for(Int_t k=0; k<jetArr->GetEntries(); k++) {
            const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[k]);
            cout << "jet " << k << " " << jet->pt << " " << jet->eta << " " << jet->phi << endl;
          }


          if (finalState == 0) {
            Count_ee += eventWeight;
            Count_ee_statError += eventWeight*eventWeight;
            Count_ee_sysError += eventWeightError;
          } else if (finalState == 1) {
            Count_mm += eventWeight;
            Count_mm_statError += eventWeight*eventWeight;
            Count_mm_sysError += eventWeightError;
          } else {
            Count_em += eventWeight;
            Count_em_statError += eventWeight*eventWeight;
            Count_em_sysError += eventWeightError;
          }

          if (lepton2.Pt() > 20.0) {
            if (fabs(lepton2.Eta()) >= 0 && fabs(lepton2.Eta()) < 1.0) {
              Count_Pt20ToInf_EtaBin0 += eventWeight;
              Count_Pt20ToInf_EtaBin0_statError += eventWeight*eventWeight;
              Count_Pt20ToInf_EtaBin0_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.0 && fabs(lepton2.Eta()) < 1.5) {
              Count_Pt20ToInf_EtaBin1 += eventWeight;
              Count_Pt20ToInf_EtaBin1_statError += eventWeight*eventWeight;
              Count_Pt20ToInf_EtaBin1_sysError += eventWeightError;
            } else if (fabs(lepton2.Eta()) >= 1.5 && fabs(lepton2.Eta()) < 2.5) {
              Count_Pt20ToInf_EtaBin2 += eventWeight;
              Count_Pt20ToInf_EtaBin2_statError += eventWeight*eventWeight;
              Count_Pt20ToInf_EtaBin2_sysError += eventWeightError;
            }
          } else if (lepton2.Pt() > 15.0 && lepton2.Pt() < 20.0 ) {
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
          } else if (lepton2.Pt() > 10.0 && lepton2.Pt() < 15.0 ) {
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

        }

        //*********************************************************************************************
        //Plots after all Cuts
        //*********************************************************************************************
        if (passAllCuts) {
          fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton, eventWeight);
          fMtLepton1_afterCuts->Fill(mTW[0], eventWeight);
          fMtLepton2_afterCuts->Fill(mTW[1], eventWeight);
          fMtHiggs_afterCuts->Fill(mtHiggs, eventWeight);
          fLeptonPtPlusMet_afterCuts->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt(), eventWeight);
    

        }

    } // end loop over electron denominators


  } //end loop over data     
  }


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
  scaleFactor = 1;
  for (int i=1; i < fHWWToEESelection->GetXaxis()->GetNbins()+1; ++i) {
    cout << fHWWToEESelection->GetBinContent(i) << "+/-" << fHWWToEESelection->GetBinError(i) << "+/-" << fHWWToEESelection_sysError->GetBinContent(i) << " " << fHWWToMuMuSelection->GetBinContent(i) << "+/-" << fHWWToMuMuSelection->GetBinError(i) << "+/-" << fHWWToMuMuSelection_sysError->GetBinContent(i) << " " << fHWWToEMuSelection->GetBinContent(i) << "+/-" << fHWWToEMuSelection->GetBinError(i) << "+/-" << fHWWToEMuSelection_sysError->GetBinContent(i) << " " << fHWWSelection->GetBinContent(i) << "+/-" << fHWWSelection->GetBinError(i) << "+/-" << fHWWSelection_sysError->GetBinContent(i) << endl;
  }

  cout << "NEventsPreSel : " << NEventsPreSel << endl;
  cout << "NEventsPreSelPassHiggsCuts : " << NEventsPreSelPassHiggsCuts << endl;
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
  cout << "Event Count : Pt20ToInf EtaBin0 \n";
  cout << "BB :" << Count_Pt20ToInf_EtaBin0*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20ToInf_EtaBin0_statError)*scaleFactor << " +/- " << Count_Pt20ToInf_EtaBin0_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20ToInf EtaBin1 \n";
  cout << "BB :" << Count_Pt20ToInf_EtaBin1*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20ToInf_EtaBin1_statError)*scaleFactor << " +/- " << Count_Pt20ToInf_EtaBin1_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";
  cout << "Event Count : Pt20ToInf EtaBin2 \n";
  cout << "BB :" << Count_Pt20ToInf_EtaBin2*scaleFactor << " +/- " << TMath::Sqrt(Count_Pt20ToInf_EtaBin2_statError)*scaleFactor << " +/- " << Count_Pt20ToInf_EtaBin2_sysError*scaleFactor <<  endl;
  cout << "**************************************************************\n";









  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("WWSelectionPlotsFakePrediction.root", "RECREATE");
//   TFile *file = new TFile("WWSelectionPlots_IsoDenominator.root", "RECREATE");

  file->WriteTObject(fHWWSelection, fHWWSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEESelection, fHWWToEESelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToMuMuSelection , fHWWToMuMuSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEMuSelection,fHWWToEMuSelection->GetName(), "WriteDelete");
  
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

        
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;
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
    pass = kTRUE;
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


Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }

  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
             && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
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

  if (mu->pt < 20) {
    if (!
        ( (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10
          && fabs(mu->d0) < 0.01
          )
      ) pass = kFALSE;    
  }

  return pass;
}


Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.2
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 1.0
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
