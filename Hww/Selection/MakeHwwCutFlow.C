//root -l EWKAna/Hww/Selection/MakeHwwCutFlowSync.C+\(140,\"\",-1\)


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
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passNewElectronCuts(const mithep::TElectron *ele, Double_t lilelihood, Int_t Option);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);

Bool_t passConversionVeto(Int_t isConv);


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

void MakeHwwCutFlow(Double_t mHiggs, const string Label, Int_t ElectronSelectionType) 
{  
  gBenchmark->Start("WWTemplate");

  ofstream eventListFile("eventList.txt");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

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

  






  Bool_t printDebug = kTRUE;
  printDebug = kFALSE;

  Double_t lumi = 35.5;              // luminosity (pb^-1)
  lumi = 1000;
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


  vector<vector<string> > inputFiles;
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full-m10-v1_TightPlusDenominatorTriggerSkim.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full-pr-v4_TightPlusDenominatorNoTriggerSkim.root");
;

//   inputFiles.push_back(vector<string>());

//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h115ww2l-gf-v1g1-pu_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h130ww2l-gf-v1g1-pu_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h160ww2l-gf-v1g1-pu_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h200ww2l-gf-v1g1-pu_noskim_normalized.root");
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjets-z2-v8-pu11-2l_noskim_normalized.root");

//   inputFiles.back().push_back("HwwNtuple_p11-h130ww2l-gf-v1g1-pu_noskim_0000.root");


//  inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("HwwAnalysis_p11-ww2l-v1g1-pu_noskim_0000.root");
// //   inputFiles.push_back(vector<string>());
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-ggww-v1g1-pu_noskim_normalized.root");

// //   inputFiles.push_back(vector<string>());
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wz-v1g1-pu_noskim_normalized.root");
// //   inputFiles.push_back(vector<string>());
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zz-v1g1-pu_noskim_normalized.root");

// //   inputFiles.push_back(vector<string>());
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-ttj-v1g1-pu_noskim_normalized.root");

// //   inputFiles.push_back(vector<string>());
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-stop-v1g1-pu_noskim_normalized.root");
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-ttop-v1g1-pu_noskim_normalized.root");
// //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wtop-v1g1-pu_noskim_normalized.root");

//    inputFiles.push_back(vector<string>());
//    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root");
// //    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zee1020-powheg-c10-v8-pu11_noskim_normalized.root");
//    inputFiles.push_back(vector<string>());
//    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root");
// //    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm1020-powheg-c10-v8-pu11_noskim_normalized.root");
// //    inputFiles.push_back(vector<string>());
// //    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zttm20-v1g1-pu_noskim_normalized.root");
// //    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt1020-powheg-c10-v8-pu11_noskim_normalized.root");

//   inputFiles.push_back(vector<string>());
//    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjets-z2-v8-pu11-2l_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-wjetsl-v1g1-pu-skimmed_TightPlusRecoNoTriggerSkim_normalized.root");








  vector<string> processNames;
  processNames.push_back("Data");
//   processNames.push_back("HWW130");
//      processNames.push_back("HWWLowPU");
//    processNames.push_back("HWWHighPU");
//    processNames.push_back("qqWW"); 
// //   processNames.push_back("ggWW"); 
// //   processNames.push_back("WZ");
// //   processNames.push_back("ZZ");
// //   processNames.push_back("ttbar");
//   processNames.push_back("DYee");
//   processNames.push_back("DYmm");
// //   processNames.push_back("DYtt");
//   processNames.push_back("WJetsEnriched");
//   processNames.push_back("WJetsMadgraph");

  assert(processNames.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  Double_t Count_ee = 0;
  Double_t Count_em = 0;
  Double_t Count_me = 0;
  Double_t Count_mm = 0;
  Double_t Count_ee_hh = 0;
  Double_t Count_ee_hl = 0;
  Double_t Count_em_hh = 0;
  Double_t Count_em_hl = 0;
  Double_t Count_me_hh = 0;
  Double_t Count_me_hl = 0;
  Double_t Count_mm_hh = 0;
  Double_t Count_mm_hl = 0;
  Double_t Count_total_hh = 0;
  Double_t Count_total_hl = 0;
  Double_t Count_total_all = 0;

  vector <TH1F*>  fHWWSelection; 
  vector <TH1F*>  fHWWToEESelection; 
  vector <TH1F*>  fHWWToMuMuSelection; 
  vector <TH1F*>  fHWWToEMuSelection; 
  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHWWSelection= new TH1F(("hHWWSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEESelection= new TH1F(("hHWWToEESelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToMuMuSelection= new TH1F(("hHWWToMuMuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEMuSelection= new TH1F(("hHWWToEMuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    tmpHWWSelection->Sumw2();
    tmpHWWToEESelection->Sumw2();
    tmpHWWToMuMuSelection->Sumw2();
    tmpHWWToEMuSelection->Sumw2();
    
    fHWWSelection.push_back(tmpHWWSelection);
    fHWWToEESelection.push_back(tmpHWWToEESelection);
    fHWWToMuMuSelection.push_back(tmpHWWToMuMuSelection);
    fHWWToEMuSelection.push_back(tmpHWWToEMuSelection);
  }

  vector<Double_t> Count_Total_WWSelection;
  vector<Double_t> Count_Total_HWWSelection;
  vector<Double_t> Count_Total;
  vector<Double_t> Count_EtaBin0;
  vector<Double_t> Count_EtaBin1;
  vector<Double_t> Count_EtaBin2;
  vector<Double_t> Count_Total_WWSelection_statError;
  vector<Double_t> Count_Total_HWWSelection_statError;
  vector<Double_t> Count_Total_statError;
  vector<Double_t> Count_EtaBin0_statError;
  vector<Double_t> Count_EtaBin1_statError;
  vector<Double_t> Count_EtaBin2_statError;
  for (int q=0; q<processNames.size() ; ++q) {
    Count_Total_WWSelection.push_back(0);
    Count_Total_HWWSelection.push_back(0);
    Count_Total.push_back(0);
    Count_EtaBin0.push_back(0);
    Count_EtaBin1.push_back(0);
    Count_EtaBin2.push_back(0);
    Count_Total_WWSelection_statError.push_back(0);
    Count_Total_HWWSelection_statError.push_back(0);
    Count_Total_statError.push_back(0);
    Count_EtaBin0_statError.push_back(0);
    Count_EtaBin1_statError.push_back(0);
    Count_EtaBin2_statError.push_back(0);  
  }



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
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
  rlrm.AddJSONFile("Cert_160404-166861_7TeV_PromptReco_Collisions11_JSON.txt"); 
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
  
      cout << "Reading File : " << inputFiles[q][f] << " : " << eventTree->GetEntries() << " Events\n";
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
        if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
		
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
        //for the skimmed input, I already required the HLT bits.
        //    if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;



        //********************************************************
        // Printdebug
        //********************************************************
        printDebug = kFALSE;
        if ((0 == 1) 
            || (info->evtNum == 9316)
             || (info->evtNum == 9318)
             || (info->evtNum == 9319)
            || (info->evtNum == 9320)
            || (info->evtNum == 9321)
            || (info->evtNum == 9322)
            || (info->evtNum == 9323)
            || (info->evtNum == 9324)
            || (info->evtNum == 9325)
            || (info->evtNum == 9326)
            || (info->evtNum == 9327)
            || (info->evtNum == 9328)
            || (info->evtNum == 9329)
            || (info->evtNum == 9330)
            || (info->evtNum == 9331)
            || (info->evtNum == 9332)
          ) printDebug = kTRUE;
        



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
//                  &&
//                   passNewMuonCuts(mu)
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
          measurements.OneOverEMinusOneOverP = (1.0/tmpEleSC.E())*(1 - ele->EOverP);
          double likelihood = LH->result(measurements);
          Double_t likelihoodValue = 0; 
          if (likelihood <= 0) {
            likelihoodValue = -20.0;
          } else if (likelihood == 1) {
            likelihoodValue = 20.0;
          } else {
            likelihoodValue = TMath::Log(likelihood / (1.0-likelihood));
          }



          if ( (0==0)
               &&
               passElectronCuts(ele)
//                   &&
//                  passNewElectronCuts(ele, likelihoodValue, 13 )
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
        //Count Jets
        //******************************************************************************
        double maxBtag = -99999;
        for(Int_t i=0; i<jetArr->GetEntries(); i++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
          if (jet->pt < 7.0) continue;

          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (leptonType[k] == 11) {
              const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[leptonIndex[k]]);
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta,tmpEle->phi, tmpEle->eta) < 0.5) {
                leptonOverlap = kTRUE;
              }
            } else {
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.5) {
                leptonOverlap = kTRUE;
              }
            }
          }

          if (!leptonOverlap) {
            if (jet->pt > 30.0 && fabs(jet->eta) < 5.0 ) {
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

//         eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, 0,
//                   leptonPt[0], leptonEta[0], leptonPhi[0], leptonCharge[0],  leptonPt[1], leptonEta[1], leptonPhi[1], leptonCharge[1]);
        
        
        //******************************************************************************
        //Correct the MET for the leptons
        //******************************************************************************

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
          }
        }

        pfTrackMet.SetXYZ(info->pfTrackMEx + MetCorrection_X, info->pfTrackMEy + MetCorrection_Y , 0);



        //******************************************************************************
        //Pt Preselection
        //******************************************************************************
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
//             if (leptonType[j] == 11) {
//               if (!(lepton2.Pt() > 15.0)) PreselPtCut = kFALSE;
//             }
//             if (!(mithep::MathUtils::DeltaR(lepton1.Phi(), lepton1.Eta() , lepton2.Phi(), lepton2.Eta()) > 0.5)) PreselPtCut = kFALSE;

            if(PreselPtCut) passCut[0] = true;
            if(zDiffMax < 100000.0) passCut[1] = true;            
            if(pfMet.Pt()    > 20.0)               passCut[2] = true; 
  
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(!(fabs(dilepton.M()-91.1876)   > 15.0))   passCut[4] = true;
              if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              if(!(fabs(dilepton.M()-91.1876)   > 15.0))   passCut[4] = true;
              //passCut[4] = true;
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

            //For event list matching
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
              eventListFile << info->evtNum << endl;
            }

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
 
            if ( (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9])
              ) {
              Count_Total_WWSelection[q]++;
              Count_Total_WWSelection_statError[q]++;
            }


            if ( (passCut[0] && passCut[1] 
                  && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] 
//                   && passCut[10] && passCut[11] && passCut[12] && passCut[13] && passCut[14]
                   )) {

              if (0==0 
//                   && 
//                   fabs(lepton2.Eta()) > 1.5
                ) {
                if (finalState == 0) {
                  Count_ee++;
                  if (lepton2.Pt() > 20) Count_ee_hh++;
                  else Count_ee_hl++;
                }
                if (finalState == 1) {
                  Count_mm++;
                  if (lepton2.Pt() > 20) Count_mm_hh++;
                  else Count_mm_hl++;
                }
                if (finalState == 2) {
                  Count_em++;
                  if (lepton2.Pt() > 20) Count_em_hh++;
                  else Count_em_hl++;
                }
                if (finalState == 3) {
                  Count_me++;
                  if (lepton2.Pt() > 20) Count_me_hh++;
                  else Count_me_hl++;
                }

                if (lepton2.Pt() > 20) Count_total_hh++;
                else Count_total_hl++;
                
                Count_total_all++;
                

              }
            }
            
            if (
              (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9])
                //passAllCuts               
              ) {    


               Count_Total_HWWSelection[q]++;
              Count_Total_HWWSelection_statError[q]++;
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

  //--------------------------------------------------------------------------------------------------------------
  // make tex cutflow table
  //============================================================================================================== 
  std::ofstream texFile("HwwCutFlowTableSync.tex");
  
  texFile << "\\begin{table}[!ht]" << endl;
  texFile << "  \\begin{center}" << endl;
  texFile << " {\\small" << endl;
  texFile << "  \\begin{tabular} {|c|c|c|c|c|c|c|c|}" << endl;
  texFile << "\\hline" << endl;
//   texFile << "  Cut & " << processNames[0] << " & " << processNames[1] << " & " << processNames[2] << " & " << processNames[3] << " & " << processNames[4] << " & " << processNames[5] << " & " << processNames[5] << " & " << processNames[6] << " \\" << endl;
  texFile << "  \\hline" << endl;
  texFile << "  \\hline" << endl;

  for(int c=0; c<CutLabel.size() ; ++c) {
    texFile << CutLabel[c] << "   & "  ;
    for (int q = 0; q<processNames.size() ; ++q) { 
      char tmp[20];
      sprintf(tmp, "%.4f +/- %.4f", fHWWSelection[q]->GetBinContent(1+c),fHWWSelection[q]->GetBinError(1+c));
//        sprintf(tmp, "%.4f", fHWWSelection[q]->GetBinContent(1+c) / 106350 );
      texFile << tmp << " & " << endl;
    }
    texFile << endl;
  }

  texFile << " \\hline" << endl;
  texFile << "  \\end{tabular}" << endl;
  texFile << "  }" << endl;
  texFile << "  \\caption{}" << endl;
  texFile << "   \\label{tab:HwwCutFlow}" << endl;
  texFile << "  \\end{center}" << endl;
  texFile << "\\end{table}" << endl;

  double scaleFactor = 1.0;
  for (int q = 0; q<processNames.size() ; ++q) { 
    cout << "Process : " << processNames[q] << endl;
    cout << "**************************************************************\n";
    cout << "Event Count : Total \n";
    cout << "BB :" << Count_Total[q]*scaleFactor << " +/- " << TMath::Sqrt(Count_Total_statError[q])*scaleFactor <<  endl;
    cout << "**************************************************************\n";
    cout << "Event Count : EtaBin0 \n";
    cout << "BB :" << Count_EtaBin0[q]*scaleFactor << " +/- " << TMath::Sqrt(Count_EtaBin0_statError[q])*scaleFactor <<  endl;
    cout << "**************************************************************\n";
    cout << "Event Count : EtaBin1 \n";
    cout << "BB :" << Count_EtaBin1[q]*scaleFactor << " +/- " << TMath::Sqrt(Count_EtaBin1_statError[q])*scaleFactor  <<  endl;
    cout << "**************************************************************\n";
    cout << "Event Count : EtaBin2 \n";
    cout << "BB :" << Count_EtaBin2[q]*scaleFactor << " +/- " << TMath::Sqrt(Count_EtaBin2_statError[q])*scaleFactor <<  endl;
    cout << "**************************************************************\n";


    Double_t num = Count_Total_HWWSelection[q];
    Double_t den = Count_Total_WWSelection[q];
    Double_t ratio = 0.0;
    Double_t errLow = 0.0;
    Double_t errHigh = 0.0;
    mithep::MathUtils::CalcRatio(num , den, ratio, errLow, errHigh, 2);      


    cout << "\nEfficiency of HWW Cuts\n";
    cout << num << " / " << den  << " = " << ratio << " + " << errHigh << " - " << errLow << endl;
  }

  cout << "COUNTS:\n\n";
  cout << "EE: " << Count_ee << " : " << Count_ee_hh << " , " << Count_ee_hl << endl;
  cout << "EM: " << Count_em << " : " << Count_em_hh << " , " << Count_em_hl << endl;
  cout << "ME: " << Count_me << " : " << Count_me_hh << " , " << Count_me_hl << endl;
  cout << "MM: " << Count_mm << " : " << Count_mm_hh << " , " << Count_mm_hl << endl;
  cout << "All: " << Count_total_all << " : " << Count_total_hh << " , " << Count_total_hl << endl;
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlots.root", "RECREATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
    file->WriteTObject(fHWWSelection[q], fHWWSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToEESelection[q], fHWWToEESelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToMuMuSelection[q], fHWWToMuMuSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToEMuSelection[q],fHWWToEMuSelection[q]->GetName(), "WriteDelete");
  }

  file->Close();
  delete file;

        
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t isMC = kFALSE;
  Bool_t pass = kTRUE;

  return pass;

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

Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (!(fabs(ele->eta) <= 2.5)) pass = kFALSE;

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
      if (fabs(ele->scEta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 )) pass = kFALSE;
      }
    }
  }

  return pass;
}



Bool_t passNewElectronCuts(const mithep::TElectron *ele, Double_t likelihoodValue, Int_t Option) {
  
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

  Double_t LHCutValue = 0;
  

  if (Option == 0) {
    //Baseline 90%
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -4.274;
        } else {
          LHCutValue = -4.274;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -3.773;
        } else {
          LHCutValue = -3.773;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -5.092;
        } else {
          LHCutValue = -5.092;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -2.796;
        } else {
          LHCutValue = -2.796;
        }
      }
    }
  }


  if (Option == 1) {
    //Baseline 90%
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -1.497;
        } else {
          LHCutValue = -1.497;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -1.521;
        } else {
          LHCutValue = -1.521;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -2.571;
        } else {
          LHCutValue = -2.571;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -0.657;
        } else {
          LHCutValue = -0.657;
        }
      }
    }
  }

  if (Option == 2) {
    //Baseline 85%
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 0.163;
        } else {
          LHCutValue = 0.163;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 0.065;
        } else {
          LHCutValue = 0.065;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -0.683;
        } else {
          LHCutValue = -0.683;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 1.564;
        } else {
          LHCutValue = 1.564;
        }
      }
    }
  }

  if (Option == 3) {
    //Baseline 80%
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 1.193;
        } else {
          LHCutValue = 1.193;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 1.345;
        } else {
          LHCutValue = 1.345;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 0.810;
        } else {
          LHCutValue = 0.810;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 3.021;
        } else {
          LHCutValue = 3.021;
        }
      }
    }
  }

  if (Option == 4) {
    //Baseline 70%  
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 1.781;
        } else {
          LHCutValue = 1.781;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 2.397;
        } else {
          LHCutValue = 2.397;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 2.361;
        } else {
          LHCutValue = 2.361;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 4.052;
        } else {
          LHCutValue = 4.052;
        }
      }
    }
  }

  //90% for pt >20, 80% for pt < 20
  if (Option == 13) {
    //Baseline 70%  
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -1.497;
        } else {
          LHCutValue = 1.193;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -1.521;
        } else {
          LHCutValue = 1.345;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -2.571;
        } else {
          LHCutValue = 0.810;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -0.657;
        } else {
          LHCutValue = 3.021;
        }
      }
    }
  }


  //90% for pt >20, 70% for pt < 20
  if (Option == 14) {
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -1.497;
        } else {
          LHCutValue = 1.781;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -1.521;
        } else {
          LHCutValue = 2.397;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = -2.571;
        } else {
          LHCutValue = 2.361;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = -0.657;
        } else {
          LHCutValue = 4.052;
        }
      }
    }
  }

  //85% for pt >20, 70% for pt < 20
  if (Option == 24) {
    if (fabs(ele->scEta) < 1.479) {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue = 0.163;
        } else {
          LHCutValue = 1.781;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 0.065;
        } else {
          LHCutValue = 2.397;
        }
      }
    } else {
      if (ele->nBrem == 0) {
        if (ele->pt > 20) {
          LHCutValue =  -0.683;
        } else {
          LHCutValue = 2.361;
        }
      } else {
        if (ele->pt > 20) {
          LHCutValue = 1.564;
        } else {
          LHCutValue = 4.052;
        }
      }
    }
  }






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


Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) > 2.4) pass = kFALSE;

  Double_t iso04 = mu->ChargedIso04 + mu->NeutralIso04_10Threshold;
  Double_t iso03 = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
  Double_t iso = iso03;
  Double_t isoCutValue = 0;


//   TightCone03  WP
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
