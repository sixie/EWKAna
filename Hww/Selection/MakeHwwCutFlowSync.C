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
#include "EWKAna/Utils/LeptonIDCuts.hh"
#include "EWKAna/Utils/HiggsCuts.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
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

void MakeHwwCutFlowSync(Double_t mHiggs = 120.0, Int_t NJetBin = 0, const string Label = "") 
{  
  gBenchmark->Start("WWTemplate");

  ofstream eventListFile("eventList.HWWNtuple.txt");
  
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

  mithep::ElectronIDMVA *electronIDMVANoIPInfo = new mithep::ElectronIDMVA();
  electronIDMVANoIPInfo->Initialize("BDTG method",
                              "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml",
                              mithep::ElectronIDMVA::kNoIPInfo );

  mithep::ElectronIDMVA *electronIDMVAWithIPInfo = new mithep::ElectronIDMVA();
  electronIDMVAWithIPInfo->Initialize("BDTG method",
                              "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml", 
                              "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",
                              mithep::ElectronIDMVA::kWithIPInfo);


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
 
  Double_t mH = mHiggs;

  Int_t channel = -1;
  if     (mH == 115) channel = 0;
  else if(mH == 120) channel = 1;
  else if(mH == 130) channel = 2;
  else if(mH == 140) channel = 3;
  else if(mH == 150) channel = 4;
  else if(mH == 160) channel = 5;
  else if(mH == 170) channel = 6;
  else if(mH == 180) channel = 7;
  else if(mH == 190) channel = 8;
  else if(mH == 200) channel = 9;
  else if(mH == 210) channel = 10;
  else if(mH == 220) channel = 11;
  else if(mH == 230) channel = 12;
  else if(mH == 250) channel = 13;
  else if(mH == 300) channel = 14;
  else if(mH == 350) channel = 15;
  else if(mH == 400) channel = 16;
  else if(mH == 450) channel = 17;
  else if(mH == 500) channel = 18;
  else if(mH == 550) channel = 19;
  else if(mH == 600) channel = 20;

  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};

  
  float dilmass_cut = 10000;   
  if     ( mH == 115 ) dilmass_cut =  70.0;
  else if( mH == 120 ) dilmass_cut =  70.0;
  else if( mH == 130 ) dilmass_cut =  80.0;
  else if( mH == 140 ) dilmass_cut =  90.0;
  else if( mH == 150 ) dilmass_cut = 100.0;
  else if( mH == 160 ) dilmass_cut = 100.0;
  else if( mH == 165 ) dilmass_cut = 100.0;
  else if( mH == 170 ) dilmass_cut = 100.0;
  else if( mH == 180 ) dilmass_cut = 110.0;
  else if( mH == 190 ) dilmass_cut = 120.0;
  else if( mH == 200 ) dilmass_cut = 130.0;
  else if( mH == 210 ) dilmass_cut = 140.0;
  else if( mH == 220 ) dilmass_cut = 150.0;
  else                 dilmass_cut = mH;


  vector<vector<string> > inputFiles;
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full-m10-v1_TightPlusDenominatorTriggerSkim.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/021/SmurfV6/WWAnalysisSkimmed_r11a-full-pr-v4_TightPlusDenominatorNoTriggerSkim.root");
;

   inputFiles.push_back(vector<string>());

//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-h115ww2l-gf-v1g1-pu_noskim_normalized.root");
//    inputFiles.back().push_back("/data/blue/sixie/ntuples/HWW/mc/HwwAnalysis_s11-h120ww2l-gf-v11-pu_noskim_normalized.root");
   inputFiles.back().push_back("/data/blue/sixie/ntuples/HWW/synchronization/AllNtuplerTest_HWWNtuple_s11-h120ww2l-gf-v11-pu_noskim_0000.root");


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
//   processNames.push_back("Data");
   processNames.push_back("HWW120");
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
  CutLabel.push_back("dPhiDiLepJet1");
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
  rlrm.AddJSONFile("/data/smurf/sixie/data/auxiliar/HWWLP11Cert.json"); 
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
        // Printdebug
        //********************************************************
        Bool_t printDebug = kFALSE;
        if ((0 == 1) ||
            (info->evtNum == 43336) ||
            (info->evtNum == 99049) 
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
        const mithep::TJet *jet2 = 0;
    
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
          if ( (0==0)
               &&
               passMuonID(mu)
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
               && fabs(mu->dz) < 0.2
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

          Double_t BDTGNoIPInfo = electronIDMVANoIPInfo->MVAValue(
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
            ele->ip3dSig
            );
          Double_t BDTGWithIPInfo = electronIDMVAWithIPInfo->MVAValue(
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
            ele->ip3dSig
            );

          if ( (0==0)
               &&
               passEleID(ele)
//                passElectronMVA(ele, BDTGNoIPInfo, 1)
//                passElectronMVA(ele, BDTGWithIPInfo, 2)
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
          if (jet->rawPt < 7.0) continue;

          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (leptonType[k] == 11) {
              const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[leptonIndex[k]]);
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta,tmpEle->phi, tmpEle->eta) < 0.3) {
                leptonOverlap = kTRUE;
              }
            } else {
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
                leptonOverlap = kTRUE;
              }
            }
          }

          if (!leptonOverlap) {
            if (fabs(jet->eta) < 5.0) {
              if (!leadingJet || jet->pt > leadingJet->pt) {
                if (leadingJet) jet2 = leadingJet;
                leadingJet = jet;
              } else if (!jet2 || jet->pt > jet2->pt) {
                jet2 = jet;
              }
            }
          
            if (fabs(jet->eta) < 5.0 && jet->pt > 30.0) {
              NJets++;
            }
            if (jet->TrackCountingHighEffBJetTagsDisc > maxBtag ) maxBtag = jet->TrackCountingHighEffBJetTagsDisc;
          }
        }


        if (printDebug) {
          cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
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
              if(ele->pt > 20) {
                cout << Bool_t(ele->sigiEtaiEta < 0.01) << " " 
                     << Bool_t(fabs(ele->deltaEtaIn) < 0.004 ) << " "
                     << Bool_t(fabs(ele->deltaPhiIn) < 0.06 ) << " "
                     << Bool_t(ele->HoverE < 0.04) << " ";

              } else {
                cout << Bool_t(ele->sigiEtaiEta < 0.01) << " " 
                     << Bool_t(fabs(ele->deltaEtaIn) < 0.004 ) << " "
                     << Bool_t(fabs(ele->deltaPhiIn) < 0.03 ) << " "
                     << Bool_t(ele->HoverE < 0.025) << " ";
              }
            } else {
              if (ele->pt > 20) {
                cout << Bool_t(ele->sigiEtaiEta < 0.03) << " " 
                     << Bool_t(fabs(ele->deltaEtaIn) < 0.007 ) << " "
                     << Bool_t(fabs(ele->deltaPhiIn) < 0.03 ) << " "
                     << Bool_t(ele->HoverE < 0.10) << " ";
              } else {
                cout << Bool_t(ele->sigiEtaiEta < 0.03) << " " 
                     << Bool_t(fabs(ele->deltaEtaIn) < 0.005 ) << " "
                     << Bool_t(fabs(ele->deltaPhiIn) < 0.02 ) << " "
                     << Bool_t(ele->HoverE < 0.10) << " ";
              }
              cout << " - " << Bool_t(ele->fBrem > 0.15 || (fabs(ele->eta) < 1.0 && ele->EOverP > 0.95)) << " - " ;
            }

            cout << Bool_t((ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold+ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold-ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold) / ele->pt < isoCutValue ) << " "
                 << Bool_t(ele->nExpHitsInner <= 0 ) << " "
                 << Bool_t(passConversionVeto(ele->isConv) ) << " "
                 << Bool_t(fabs(ele->d0) < 0.02 ) << " "
                 << Bool_t(fabs(ele->dz) < 0.1 ) << " "
                 << Bool_t(ele->isEcalDriven ) << " : "
                 << passEleID(ele) << " "
                 << endl;


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
            
            cout << "Likelihood : " <<  info->evtNum << " " << ele->pt << " " << ele->eta << " : " << measurements.subdet << " " << measurements.deltaPhi << " " << measurements.deltaEta << " "
                 << measurements.eSeedClusterOverPout << " " 
                 << measurements.eSuperClusterOverP << " " 
                 << measurements.hadronicOverEm << " " 
                 << measurements.sigmaIEtaIEta << " "  
                 << measurements.sigmaIPhiIPhi << " "  
                 << measurements.fBrem << " " 
                 << measurements.nBremClusters << " " 
                 << measurements.OneOverEMinusOneOverP << " " 
                 << " : " << likelihood << " " << likelihoodValue << "  === " << tmpEleSC.E() << " " << ele->scEt << " " << ele->scEta << " " << ele->scPhi << " " << ele->EOverP << endl;
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
                 << Bool_t(fabs(mu->dz) < 0.2 ) << " "
                 << Bool_t((mu->ChargedIso03 + mu->NeutralIso03_10Threshold ) / mu->pt < isoCutValue ) << " "
                 << Bool_t(mu->nMatch > 1 ) << " "
                 << Bool_t(mu->nPixHits > 0) << " "
                 << Bool_t(mu->pterr / mu->pt < 0.1) << " "
                  << " : " << passMuonID(mu) << " "
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
                 && fabs(mu->dz) < 0.2
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
                  << Bool_t(fabs(mu->dz) < 0.2) << " "
                  << !isCleanMuon << " " 
                  << Bool_t((!(mu->pt > 20 && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.10))) << " : "            
                  << isSoft << endl;


          }


          for(Int_t i=0; i<jetArr->GetEntries(); i++) {
            const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
          
            cout << "Jet" << i << " : " << jet->pt << " " << jet->eta << " " << jet->phi << endl;
            if (jet->rawPt < 7.0) continue;

            Bool_t leptonOverlap = kFALSE;
            for (int k=0; k<leptonPt.size(); ++k) {
              if (leptonType[k] == 11) {
                const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[leptonIndex[k]]);
                if (mithep::MathUtils::DeltaR(jet->phi, jet->eta,tmpEle->phi, tmpEle->eta) < 0.3) {
                  leptonOverlap = kTRUE;
                }
              } else {
                if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.3) {
                  leptonOverlap = kTRUE;
                }
              }
            }
            cout << "leptonOverlap: " << leptonOverlap << " : " << jet->TrackCountingHighEffBJetTagsDisc << endl;
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


            //leadingJet angle with dilepton system
            bool passdPhiDiLepJet1Cut = (!leadingJet) || leadingJet->pt <= 15.0 ||
              (mithep::MathUtils::DeltaPhi( Double_t(leadingJet->phi), Double_t(dilepton.Phi()) )*180.0/TMath::Pi() < 165. || 
               (finalState == 2 || finalState == 3));
               

            Bool_t passVBFHWWSelection = kFALSE;
            if(NJetBin == 2){
              int centrality = 0;             
              if( leadingJet && jet2 ) {
                
                if(((leadingJet->eta - lepton1.Eta() > 0 && jet2->eta - lepton1.Eta() < 0) ||
                    (jet2->eta - lepton1.Eta() > 0 && leadingJet->eta - lepton1.Eta() < 0)) &&
                   ((leadingJet->eta - lepton2.Eta() > 0 && jet2->eta - lepton2.Eta() < 0) ||
                    (jet2->eta - lepton2.Eta() > 0 && leadingJet->eta - lepton2.Eta() < 0))) centrality = 1; 
                
              mithep::FourVectorM j1;
              mithep::FourVectorM j2;
              j1.SetCoordinates(leadingJet->pt, leadingJet->eta, leadingJet->phi, leadingJet->mass );
              j2.SetCoordinates(jet2->pt, jet2->eta, jet2->phi, jet2->mass );
              mithep::FourVectorM dijet = j1+j2;

              passVBFHWWSelection = dijet.M() > 450. &&
                TMath::Abs(leadingJet->eta - jet2->eta) > 3.5 &&
                (mH > 200 || dilepton.M() < 100.) &&
                centrality == 1;
              }            
            }


   
            //*********************************************************************************************
            //Define Cuts
            //*********************************************************************************************
            const int nCuts = 16;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false, false,
                                   false, false, false, false, false};
            assert(CutLabel.size() == nCuts+1);

            Bool_t PreselPtCut = kTRUE;
            if(!(lepton1.Pt() >  20.0 && lepton2.Pt() > 10.0)) 
              PreselPtCut = kFALSE;
            if(PreselPtCut)                                                        passCut[0] = true;
            if(zDiffMax < 100000.0)                                                passCut[1] = true;            
            if(pfMet.Pt()    > 20.0)                                               passCut[2] = true; 
  
            if(dilepton.M() > 12.0 && dilepton.M() <= dilmass_cut)                 passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)                              passCut[4] = true;
              if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 40)          passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(TMath::Min(PFMETdeltaPhilEt,PFTrackMETdeltaPhilEt) > 20)          passCut[5] = true;
            }

            if(NJets == NJetBin)                                                   passCut[6] = true; 
        

            if (NSoftMuons == 0 )                                                  passCut[7] = true;

            if (!(leptonPt.size() >= 3 && leptonPt[2] > 10.0))                     passCut[8] = true;

            if(maxBtag < 2.1)                                                      passCut[9] = true;
            if (passdPhiDiLepJet1Cut)                                              passCut[10] = true;

            if (NJetBin == 0 || NJetBin == 1) {
              if (lepton1.Pt() > cutPtMaxLow[channel])                             passCut[11] = true;
              if (lepton2.Pt() > cutPtMinLow[channel])                             passCut[12] = true;
              if (dilepton.M() < cutMassHigh[channel])                             passCut[13] = true;
              if (deltaPhiLeptons < cutDeltaphilHigh[channel])                     passCut[14] = true;
              if (mtHiggs >= cutMTLow[channel] && mtHiggs <= cutMTHigh[channel] )  passCut[15] = true;
            } 
            if (NJetBin == 2 ) {
              if (passVBFHWWSelection) {
                passCut[11] = true;
                passCut[12] = true;
                passCut[13] = true;
                passCut[14] = true;
                passCut[15] = true;
              }              
            }
            
            //*********************************************************************************************************
            //Final State Selection
            //*********************************************************************************************************
//             if (!(finalState == 2 || finalState == 3)) continue;
//             if (!(finalState == 0 || finalState == 1)) continue;


            if (printDebug) {
              cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
              cout << "PV: " << info->pvx << " " << info->pvy << " " << info->pvz << endl;
              cout << lepton1.Pt() << " " << lepton1.Eta() << " " << lepton1.Phi() << " " << leptonType[i] << endl;
              cout << lepton2.Pt() << " " << lepton2.Eta() << " " << lepton2.Phi() << " " << leptonType[j] << endl;
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
                    if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, tmpEle->scPhi, tmpEle->scEta) < 0.3) {
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
              
              if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10]) {
                cout << "PassEvent " << info->runNum << " " << info->lumiSec << " " << info->evtNum  << "\n";
              }
            }
            
            //For event list matching
            if (passCut[0] && passCut[1]) {
              eventListFile << info->evtNum << " : "
                            << passCut[2] << " "
                            << passCut[3] << " "
                            << passCut[4] << " "
                            << passCut[5] << " "
                            << passCut[6] << " "
                            << passCut[7] << " "
                            << passCut[8] << " "
                            << passCut[9] << " : "
                            << passCut[10] << " "
                            << passCut[11] << " "
                            << passCut[12] << " "
                            << passCut[13] << " "
                            << passCut[14] << " : "
                            << passCut[15] << " : ";
              if (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10]
                  && passCut[11] && passCut[12] && passCut[13] && passCut[14] && passCut[15]
                ) eventListFile << " pass ";
              eventListFile << endl;
            }
            

            //*********************************************************************************************
            //Make Selection Histograms. Number of events passing each level of cut
            //*********************************************************************************************  
            bool passAllCuts = true;
            for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
 
            if ( (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10])
              ) {
              Count_Total_WWSelection[q]++;
              Count_Total_WWSelection_statError[q]++;
            }


            if ( (passCut[0] && passCut[1] 
                  && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10] 
//                   && passCut[11] && passCut[12] && passCut[13] && passCut[14] && passCut[15]
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
              (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9] && passCut[10])
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
