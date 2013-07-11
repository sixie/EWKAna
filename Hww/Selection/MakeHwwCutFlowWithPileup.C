//root -l EWKAna/Hww/Selection/MakeHwwCutFlowWithPileup.C+\(130,\"\",-1\)


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
Bool_t passFirstElectronCuts(const mithep::TElectron *ele, Double_t PileupEnergyDensity);
Bool_t passElectronCuts(const mithep::TElectron *ele, Int_t ElectronSelectionType, Double_t likelihood);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu, Double_t PileupEnergyDensity);
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

void MakeHwwCutFlowWithPileup(Double_t mHiggs, const string Label, Int_t ElectronSelectionType) 
{  
  gBenchmark->Start("WWTemplate");

  string label = Label; if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi = 35.5;              // luminosity (pb^-1)
//   lumi = 1000;
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
  if (fHiggsMass == 130) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 15.0;
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
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12_noskim_normalized.root");
    inputFiles.push_back(vector<string>());
    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-h130ww2l-gf-z2-v8-pu11_noskim_normalized.root");



  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-vv-mg-v8-pu11_noskim_normalized.root");
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ggww-z2-v8-pu11_noskim_normalized.root");

  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wz-z2-v8-pu11_noskim_normalized.root");
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zz-z2-v8-pu11_noskim_normalized.root");

  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-tt-mg-z2-v8-pu11_noskim_normalized.root");

  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-stop-mg-z2-v8-pu11_noskim_normalized.root");
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ttop-mg-z2-v8-pu11_noskim_normalized.root");
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wtop-mg-z2-v8-pu11_noskim_normalized.root");
  
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zee-powheg-c10-v8-pu11_noskim_normalized.root");
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zee1020-powheg-c10-v8-pu11_noskim_normalized.root");
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm-powheg-c10-v8-pu11_noskim_normalized.root");
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-zmm1020-powheg-c10-v8-pu11_noskim_normalized.root");
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt-powheg-c10-v8-pu11_noskim_normalized.root");
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-ztt1020-powheg-c10-v8-pu11_noskim_normalized.root");
  
  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_w10-wjetsl-z2-v8-pu11-skimmed_TwoRecoLeptonSkim_normalized.root");








  vector<string> processNames;
//   processNames.push_back("Data");
  processNames.push_back("Hww130");
  processNames.push_back("qqWW"); 
  processNames.push_back("ggWW"); 
  processNames.push_back("WZ");
  processNames.push_back("ZZ");
  processNames.push_back("ttbar");
  processNames.push_back("SingleTop");
  processNames.push_back("DYee");
  processNames.push_back("DYmm");
  processNames.push_back("DYtt");
  processNames.push_back("WJets");

  assert(processNames.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector <TH1F*>  fHWWSelection; 
  vector <TH1F*>  fHWWToEESelection; 
  vector <TH1F*>  fHWWToMuMuSelection; 
  vector <TH1F*>  fHWWToEMuSelection; 
  vector <TH1F*>  fPtMax; 
  vector <TH1F*>  fPtMin; 
  vector <TH1F*>  fMinPFNoFwdMetPFMetDeltaPhilEt; 
  vector <TH1F*>  fNJets; 
  vector <TH1F*>  fNSoftMuons; 
  vector <TH1F*>  fBTagDiscriminator; 
  vector <TH1F*>  fDileptonMass; 
  vector <TH1F*>  fDeltaPhiLeptons; 
  vector <TH1F*>  fMtHiggs;

  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHWWSelection= new TH1F(("hHWWSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEESelection= new TH1F(("hHWWToEESelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToMuMuSelection= new TH1F(("hHWWToMuMuSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHWWToEMuSelection= new TH1F(("hHWWToEMuSelection"+string("_")+processNames[q]+label).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);

    TH1F *tmpPtMax = new TH1F(("hPtMax"+string("_")+processNames[q]+label).c_str(), ";PtMax [GeV/c];Number of Events", 100, 0, 100);
    TH1F *tmpPtMin = new TH1F(("hPtMin"+string("_")+processNames[q]+label).c_str(), ";PtMin [GeV/c];Number of Events", 100, 0, 100);
    TH1F *tmpMinPFNoFwdMetPFMetDeltaPhilEt = new TH1F(("hMinPFNoFwdMetPFMetDeltaPhilEt"+string("_")+processNames[q]+label).c_str(), ";MinPFNoFwdMetPFMetDeltaPhilEt;Number of Events", 100, 0, 100);
    TH1F *tmpNJets = new TH1F(("hNJets"+string("_")+processNames[q]+label).c_str(), ";NJets;Number of Events", 6,-0.5,5.5);
    TH1F *tmpNSoftMuons = new TH1F(("hNSoftMuons"+string("_")+processNames[q]+label).c_str(), ";NSoftMuons;Number of Events", 6,-0.5,5.5);
    TH1F *tmpBTagDiscriminator = new TH1F(("hBTagDiscriminator"+string("_")+processNames[q]+label).c_str(), ";BTagDiscriminator;Number of Events", 100, -5.0,5.0);
    TH1F *tmpDileptonMass = new TH1F(("hDileptonMass"+string("_")+processNames[q]+label).c_str(), ";DileptonMass;Number of Events", 100, 0, 200);
    TH1F *tmpDeltaPhiLeptons= new TH1F(("hDeltaPhiLeptons"+string("_")+processNames[q]+label).c_str(), ";DeltaPhiLeptons;Number of Events", 100, 0, 180);
    TH1F *tmpMtHiggs = new TH1F(("hMtHiggs"+string("_")+processNames[q]+label).c_str(), ";MtHiggs [GeV/c];Number of Events", 200, 0, 400);


    tmpHWWSelection->Sumw2();
    tmpHWWToEESelection->Sumw2();
    tmpHWWToMuMuSelection->Sumw2();
    tmpHWWToEMuSelection->Sumw2();
    tmpPtMax->Sumw2();
    tmpPtMin->Sumw2();
    tmpMinPFNoFwdMetPFMetDeltaPhilEt->Sumw2();
    tmpNJets->Sumw2();
    tmpNSoftMuons->Sumw2();
    tmpBTagDiscriminator->Sumw2();
    tmpDileptonMass->Sumw2();
    tmpDeltaPhiLeptons->Sumw2();
    tmpMtHiggs->Sumw2();
    
    fHWWSelection.push_back(tmpHWWSelection);
    fHWWToEESelection.push_back(tmpHWWToEESelection);
    fHWWToMuMuSelection.push_back(tmpHWWToMuMuSelection);
    fHWWToEMuSelection.push_back(tmpHWWToEMuSelection);
    fPtMax.push_back(tmpPtMax); 
    fPtMin.push_back(tmpPtMin); 
    fMinPFNoFwdMetPFMetDeltaPhilEt.push_back(tmpMinPFNoFwdMetPFMetDeltaPhilEt); 
    fNJets.push_back(tmpNJets); 
    fNSoftMuons.push_back(tmpNSoftMuons); 
    fBTagDiscriminator.push_back(tmpBTagDiscriminator); 
    fDileptonMass.push_back(tmpDileptonMass); 
    fDeltaPhiLeptons.push_back(tmpDeltaPhiLeptons); 
    fMtHiggs.push_back(tmpMtHiggs);


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
        //********************************************************
        // TcMet
        //********************************************************
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
               passMuonCuts(mu, info->PileupEnergyDensity)
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
//                ( (ele->pt > 20 && passFirstElectronCuts(ele)) 
//                  || 
//                  (ele->pt <= 20 && passElectronCuts(ele, ElectronSelectionType, likelihood))
//                  )
               passFirstElectronCuts(ele, info->PileupEnergyDensity)
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
            double deltaPhiDileptonMet = mithep::MathUtils::DeltaPhi(pfMet.Phi(), 
                                                                     dilepton.Phi())*180.0 / TMath::Pi();    
            double mtHiggs = TMath::Sqrt(2.0*dilepton.Pt() * pfMet.Phi()*
                                         (1.0 - cos(deltaPhiDileptonMet * TMath::Pi() / 180.0)));

            //angle between MET and closest lepton
            double deltaPhiMetLepton[2] = {mithep::MathUtils::DeltaPhi(pfMet.Phi(), lepton1.Phi()),
                                           mithep::MathUtils::DeltaPhi(pfMet.Phi(), lepton2.Phi())};
  
            double mTW[2] = {TMath::Sqrt(2.0*lepton1.Pt()*pfMet.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[0]))),
                             TMath::Sqrt(2.0*lepton2.Pt()*pfMet.Pt()*
                                         (1.0 - cos(deltaPhiMetLepton[1])))};

            double minDeltaPhiMetLepton = (deltaPhiMetLepton[0] < deltaPhiMetLepton[1])?
              deltaPhiMetLepton[0]:deltaPhiMetLepton[1];

            double METdeltaPhilEt = pfMet.Pt();
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
  //           double MinPFTrackMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
//             if(minDeltaPhiMinPFTrackMetPFMetLepton < TMath::Pi()/2.)
//               MinPFTrackMetPFMetdeltaPhilEt = MinPFTrackMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFTrackMetPFMetLepton);
            double MinPFTrackMetPFMetdeltaPhilEt = TMath::Min(PFTrackMETdeltaPhilEt,PFMETdeltaPhilEt);

            double deltaPhiMinPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                          mithep::MathUtils::DeltaPhi(MinPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFNoFwdMetPFMetLepton = (deltaPhiMinPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFNoFwdMetPFMetLepton[1])?
              deltaPhiMinPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFNoFwdMetPFMetLepton[1];
//             double MinPFNoFwdMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
//             if(minDeltaPhiMinPFNoFwdMetPFMetLepton < TMath::Pi()/2.)
//               MinPFNoFwdMetPFMetdeltaPhilEt = MinPFNoFwdMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFNoFwdMetPFMetLepton);
            double MinPFNoFwdMetPFMetdeltaPhilEt = TMath::Min(PFNoFwdMETdeltaPhilEt,PFMETdeltaPhilEt);

            double deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[2] = {mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton1.Phi()),
                                                                    mithep::MathUtils::DeltaPhi(MinPFTrackMetPFNoFwdMetPFMet.Phi(), lepton2.Phi())};
            double minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton = (deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0] < deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1])?
              deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[0]:deltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton[1];
//             double MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = pfNoFwdMet.Pt();
//             if(minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton < TMath::Pi()/2.)
//               MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt * sin(minDeltaPhiMinPFTrackMetPFNoFwdMetPFMetLepton);
            double MinPFTrackMetPFNoFwdMetPFMetdeltaPhilEt = TMath::Min(TMath::Min(PFTrackMETdeltaPhilEt, PFNoFwdMETdeltaPhilEt),PFMETdeltaPhilEt);





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
            const int nCuts = 14;
            bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                                   false, false, false, false};
            assert(CutLabel.size() == nCuts+1);

            if(lepton1.Pt() >  20.0 &&
               lepton2.Pt() >= 10.0) passCut[0] = true;

            if(zDiffMax < 1.0)                    passCut[1] = true;
  
            if(MinPFNoFwdMetPFMet.Pt()    > 20.0)               passCut[2] = true; 
  
            if(dilepton.M() > 12.0)            passCut[3] = true;
   
            if (finalState == 0 || finalState == 1){ // mumu/ee
              if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[4] = true;
              if(MinPFTrackMetPFMetdeltaPhilEt > 35) passCut[5] = true;
            }
            else if(finalState == 2 ||finalState == 3 ) { // emu
              passCut[4] = true;
              if(MinPFTrackMetPFMetdeltaPhilEt > 20) passCut[5] = true;
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
    
            //after Z veto
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] ) {
              fPtMax[q]->Fill(lepton1.Pt(), eventweight);
              fPtMin[q]->Fill(lepton2.Pt(), eventweight);
              fMinPFNoFwdMetPFMetDeltaPhilEt[q]->Fill(MinPFTrackMetPFMetdeltaPhilEt, eventweight);
            }

            //after Met cut
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5]  ) {
              fNJets[q]->Fill(NJets, eventweight);
            }

            //after jet veto
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] ) {
              fNSoftMuons[q]->Fill(NSoftMuons, eventweight);
              fBTagDiscriminator[q]->Fill(maxBtag, eventweight);
            }
           
            //after all WW cuts
            if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]  ) {
              fDileptonMass[q]->Fill(dilepton.M(), eventweight);
              fDeltaPhiLeptons[q]->Fill(deltaPhiLeptons, eventweight);
            }
           
            //after all HWW cuts
          if (passAllCuts  ) {
            fMtHiggs[q]->Fill(mtHiggs,eventweight);
          }


            if ( (passCut[0] && passCut[1] && passCut[2] && passCut[3] &&  passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9])
              ) {
              Count_Total_WWSelection[q]++;
              Count_Total_WWSelection_statError[q]++;
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
  std::ofstream texFile("HwwCutFlowTable_2011PU.tex");
  
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

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionWithPileupPlots.root", "RECREATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
    file->WriteTObject(fHWWSelection[q], fHWWSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToEESelection[q], fHWWToEESelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToMuMuSelection[q], fHWWToMuMuSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHWWToEMuSelection[q],fHWWToEMuSelection[q]->GetName(), "WriteDelete");

    file->WriteTObject(fPtMax[q],fPtMax[q]->GetName(), "WriteDelete");
    file->WriteTObject(fPtMin[q],fPtMin[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMinPFNoFwdMetPFMetDeltaPhilEt[q],fMinPFNoFwdMetPFMetDeltaPhilEt[q]->GetName(), "WriteDelete");
    file->WriteTObject(fNJets[q],fNJets[q]->GetName(), "WriteDelete");
    file->WriteTObject(fNSoftMuons[q],fNSoftMuons[q]->GetName(), "WriteDelete");
    file->WriteTObject(fBTagDiscriminator[q],fBTagDiscriminator[q]->GetName(), "WriteDelete");
    file->WriteTObject(fDileptonMass[q],fDileptonMass[q]->GetName(), "WriteDelete");
    file->WriteTObject(fDeltaPhiLeptons[q],fDeltaPhiLeptons[q]->GetName(), "WriteDelete");
    file->WriteTObject(fMtHiggs[q],fMtHiggs[q]->GetName(), "WriteDelete");

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


Bool_t passFirstElectronCuts(const mithep::TElectron *ele, Double_t PileupEnergyDensity) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  Double_t PileupIsoEnergy = PileupEnergyDensity * 3.1416 * pow(0.3,2);
  Double_t PileupHOverE = PileupEnergyDensity * 3.1416 * pow(0.15,2) / ele->e;

  Double_t isoCut = 0.07;
  if (ele->pt > 20) {
    if (fabs(ele->eta) > 1.5) isoCut = 0.05;
  } else {
    if (fabs(ele->eta) > 1.5) {
      isoCut = 0.05;
    } else {
      isoCut = 0.035;
    }
  }

  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && (ele->trkIso03 + TMath::Max(TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03  - PileupIsoEnergy,0.0) )  / ele->pt < isoCut
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
             && (ele->trkIso03 + TMath::Max(ele->emIso03  + ele->hadIso03 - PileupIsoEnergy,0.0) ) / ele->pt < isoCut
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

 if (ElectronSelectionType == 0) {
    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < 0.1
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < 0.1
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

  if (ElectronSelectionType == 1) {
    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / TMath::Max(ele->pt,Float_t(20)) < 0.1
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / TMath::Max(ele->pt,Float_t(20)) < 0.1
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

  if (ElectronSelectionType == 2) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

  if (ElectronSelectionType == 2000) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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


  if (ElectronSelectionType == -2) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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


  if (ElectronSelectionType == 20) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + ele->emIso04 + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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


  if (ElectronSelectionType == 3) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.02;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

  if (ElectronSelectionType == 4) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.01;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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


  //TighterIso2 + VBTF70
  if (ElectronSelectionType == 21) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.03
              && ele->HoverE < 0.025
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
               && ele->HoverE < 0.025
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

  //TighterIso2 + VBTF60
  if (ElectronSelectionType == 22) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.025
              && ele->HoverE < 0.025
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
               && ele->HoverE < 0.025
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

 if (ElectronSelectionType == 3) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.02;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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


 if (ElectronSelectionType == 4) {

    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.01;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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



//   //WW +  ConvVeto WrongHits0
//   if (ElectronSelectionType == 100) {
//     //Barrel 
//     if (fabs(ele->eta) < 1.5) {
//       if (! ( (0==0)
//               && ele->sigiEtaiEta < 0.01 
//               && fabs(ele->deltaEtaIn) < 0.004
//               && fabs(ele->deltaPhiIn) < 0.06
//               && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
//               && ele->nExpHitsInner <= 0
// //               && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//               && fabs(ele->d0) < 0.02
//               && ele->isConv == kTRUE
//             )
//         ) {
//         pass = kFALSE;
//       }      
//     }
//     //Endcap
//     else if (fabs(ele->eta) > 1.5) {
//       if (! (  (0==0)
//                && ele->sigiEtaiEta < 0.03
//                && fabs(ele->deltaEtaIn) < 0.007
//                && fabs(ele->deltaPhiIn) < 0.03
//                && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
//                && ele->nExpHitsInner <= 0
// //                && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//                && fabs(ele->d0) < 0.02
//                && ele->isConv == kTRUE
//             )
//         ) {
//         pass = kFALSE;
//       }
//     } else {
//       pass = kFALSE;
//       return pass;
//     }
//   }


//   //WW +  ConvVeto WrongHits1
//   if (ElectronSelectionType == 101) {

//     //Barrel 
//     if (fabs(ele->eta) < 1.5) {
//       if (! ( (0==0)
//               && ele->sigiEtaiEta < 0.01 
//               && fabs(ele->deltaEtaIn) < 0.004
//               && fabs(ele->deltaPhiIn) < 0.06
//               && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
//               && ele->nExpHitsInner <= 0
// //               && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//               && fabs(ele->d0) < 0.02
//               && ele->isConvOneWrongHit == kTRUE
//             )
//         ) {
//         pass = kFALSE;
//       }      
//     }
//     //Endcap
//     else if (fabs(ele->eta) > 1.5) {
//       if (! (  (0==0)
//                && ele->sigiEtaiEta < 0.03
//                && fabs(ele->deltaEtaIn) < 0.007
//                && fabs(ele->deltaPhiIn) < 0.03
//                && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
//                && ele->nExpHitsInner <= 0
// //                && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//                && fabs(ele->d0) < 0.02
//                && ele->isConvOneWrongHit == kTRUE
//             )
//         ) {
//         pass = kFALSE;
//       }
//     } else {
//       pass = kFALSE;
//       return pass;
//     }
//   }


  //WW + Guiseppe 1a : fbrem/eta
  if (ElectronSelectionType == 201) {
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) pass = kFALSE;
    }
  }

 //WW + Guiseppe 1b : fbrem/eta
  if (ElectronSelectionType == 202) {
    Double_t relIsoCut = 0.1;
//     if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
//     if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }
  }

 //WW + Guiseppe 4 : mva
  if (ElectronSelectionType == 204) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

    if (ele->mva < 0.4) pass = kFALSE;

  }


//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == 210) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }


//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == 211) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.02;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }

//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == 212) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.01;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }



//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == -210) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }


//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == -211) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.02;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }

//WW + Guiseppe 1ab + slidingIsolation (0.4)
  if (ElectronSelectionType == -212) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.01;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.03;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < relIsoCut
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

  }




 //WW + Guiseppe 1ab + CiC Supertight Isolation
  if (ElectronSelectionType == 220) {
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

    if (ele->fBrem < 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 && fabs(ele->deltaPhiIn*ele->q) < 0.006 )) pass = kFALSE;
      }
    }

    Int_t originalvalue = ele->passSuperTightId;
    Int_t bit0;
    Int_t bit1;
    Int_t bit2;
    Int_t bit3;
    bit0 = originalvalue % 2;
    Int_t tmp0 = floor(double(originalvalue) / 2.0);
    bit1 = tmp0 % 2;
    Int_t tmp1 = floor(double(tmp0) / 2.0);
    bit2 = tmp1 % 2;
    Int_t tmp2 = floor(double(tmp1) / 2.0);
    bit3 = tmp2 % 2;

    if (!(bit1 == 1)) {
        pass = kFALSE;      
    }

  }




  //CiC v6 Isolation only
  if (ElectronSelectionType == 300 || ElectronSelectionType == 310 || ElectronSelectionType == 320 || ElectronSelectionType == 330 || ElectronSelectionType == 340 ||
      ElectronSelectionType == 301 || ElectronSelectionType == 311 || ElectronSelectionType == 321 || ElectronSelectionType == 331 || ElectronSelectionType == 341 ||
      ElectronSelectionType == 302 || ElectronSelectionType == 312 || ElectronSelectionType == 322 || ElectronSelectionType == 332 || ElectronSelectionType == 342) {

    Double_t relIsoCut = 0.1;
    if (ElectronSelectionType == 301 || ElectronSelectionType == 311 || ElectronSelectionType == 321 || ElectronSelectionType == 331 || ElectronSelectionType == 341) {
      if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
      if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;
    }

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    Double_t fBrem = ele->fBrem;
    Double_t hOverE = ele->HoverE;
    Double_t sigmaee = ele->sigiEtaiEta;
    Double_t deltaPhiIn = TMath::Abs(ele->deltaPhiIn);
    Double_t deltaEtaIn = TMath::Abs(ele->deltaEtaIn);
    Double_t eSeedOverPin = ele->ESeedClusterOverPout; 
    
    Int_t mishits = ele->nExpHitsInner;
    Double_t tkIso   = ele->trkIso03;
    Double_t ecalIso = ele->emIso04;
    Double_t hcalIso  = ele->hadIso04;
    
    int cat = Classify(ele);

    Double_t cutdcotdistSuperTight[9] = {
      2.11e-02, 1.86e-02, 1.55e-02, 3.40e-02, 2.85e-02, 3.32e-02, 1.64e-02, 3.75e-02, 1.30e-04};
    Double_t cutdetainSuperTight[9] = {
      7.84e-03, 3.67e-03, 7.00e-03, 1.28e-02, 5.65e-03, 9.53e-03, 1.08e-02, 2.97e-02, 7.24e-03};
    Double_t cutdetainlSuperTight[9] = {
      7.61e-03, 3.28e-03, 6.57e-03, 1.03e-02, 5.05e-03, 8.55e-03, 1.07e-02, 2.94e-02, 4.10e-03};
    Double_t cutdphiinSuperTight[9] = {
      4.83e-02, 7.39e-02, 2.38e-01, 5.74e-02, 1.29e-01, 2.13e-01, 3.31e-01, 3.93e-01, 2.84e-01};
    Double_t cutdphiinlSuperTight[9] = {
      5.79e-02, 7.21e-02, 2.18e-01, 7.70e-02, 1.41e-01, 2.11e-01, 2.43e-01, 3.53e-01, 2.89e-01};
    Double_t cuteseedopcorSuperTight[9] = {
      7.32e-01, 9.77e-01, 9.83e-01, 8.55e-01, 4.31e-01, 7.35e-01, 4.18e-01, 9.99e-01, 5.89e-01};
    Double_t cutfmishitsSuperTight[9] = {
      3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
    Double_t cuthoeSuperTight[9] = {
      9.19e-02, 4.11e-02, 1.42e-01, 3.35e-01, 3.82e-02, 1.41e-01, 4.29e-01, 4.01e-01, 3.99e-01};
    Double_t cuthoelSuperTight[9] = {
      7.51e-02, 3.81e-02, 1.41e-01, 3.32e-01, 3.10e-02, 1.43e-01, 2.35e-01, 3.80e-01, 1.32e-01};
    Double_t cutip_gsfSuperTight[9] = {
      1.42e-02, 2.66e-02, 1.06e-01, 3.38e-02, 3.23e-01, 1.07e-01, 7.74e-02, 2.32e-01, 7.80e-02};
    Double_t cutip_gsflSuperTight[9] = {
      1.15e-02, 2.72e-02, 8.41e-02, 2.49e-02, 4.17e-01, 1.02e-01, 7.90e-02, 1.69e-01, 4.79e-02};
    Double_t cutiso_sumSuperTight[9] = {
      8.95e+00, 8.18e+00, 8.75e+00, 7.47e+00, 5.43e+00, 5.87e+00, 8.16e+00, 1.02e+01, 1.78e+00};
    Double_t cutiso_sumoetSuperTight[9] = {
      6.45e+00, 5.14e+00, 4.99e+00, 5.21e+00, 2.65e+00, 3.12e+00, 4.52e+00, 4.72e+00, 3.68e+00};
    Double_t cutiso_sumoetlSuperTight[9] = {
      6.02e+00, 3.96e+00, 4.23e+00, 4.73e+00, 1.99e+00, 2.64e+00, 3.72e+00, 3.81e+00, 1.44e+00};
    Double_t cutseeSuperTight[9] = {
      1.09e-02, 1.05e-02, 1.05e-02, 3.24e-02, 2.81e-02, 2.95e-02, 9.77e-03, 2.75e-02, 2.95e-02};
    Double_t cutseelSuperTight[9] = {
      1.12e-02, 1.05e-02, 1.07e-02, 3.51e-02, 2.75e-02, 2.87e-02, 9.59e-03, 2.67e-02, 2.98e-02};
  // HyperTight1 cuts
  Double_t cutdcotdistHyperTight1[9] = {
  1.48e-02, 1.50e-02, 8.25e-03, 3.16e-02, 2.85e-02, 3.15e-02, 6.62e-03, 3.48e-02, 3.63e-06};
  Double_t cutdetainHyperTight1[9] = {
  6.51e-03, 3.51e-03, 5.53e-03, 9.16e-03, 5.30e-03, 8.28e-03, 1.08e-02, 2.97e-02, 7.24e-03};
  Double_t cutdetainlHyperTight1[9] = {
  6.05e-03, 3.23e-03, 4.93e-03, 8.01e-03, 4.93e-03, 7.91e-03, 1.03e-02, 2.94e-02, 4.10e-03};
  Double_t cutdphiinHyperTight1[9] = {
  4.83e-02, 4.91e-02, 2.30e-01, 3.48e-02, 7.44e-02, 2.04e-01, 9.95e-02, 3.93e-01, 2.84e-01};
  Double_t cutdphiinlHyperTight1[9] = {
  4.74e-02, 4.51e-02, 2.18e-01, 2.99e-02, 7.37e-02, 2.11e-01, 9.99e-02, 3.53e-01, 2.89e-01};
  Double_t cuteseedopcorHyperTight1[9] = {
  7.72e-01, 9.90e-01, 1.01e+00, 8.55e-01, 9.11e-01, 7.72e-01, 9.17e-01, 1.06e+00, 7.63e-01};
  Double_t cutfmishitsHyperTight1[9] = {
  3.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01};
  Double_t cuthoeHyperTight1[9] = {
  6.17e-02, 3.70e-02, 1.41e-01, 2.91e-01, 3.82e-02, 1.34e-01, 4.19e-01, 3.87e-01, 3.93e-01};
  Double_t cuthoelHyperTight1[9] = {
  4.43e-02, 3.57e-02, 1.41e-01, 2.81e-01, 3.07e-02, 1.28e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight1[9] = {
  1.21e-02, 1.76e-02, 6.01e-02, 2.96e-02, 1.74e-01, 9.70e-02, 7.74e-02, 1.33e-01, 7.80e-02};
  Double_t cutip_gsflHyperTight1[9] = {
  1.01e-02, 1.56e-02, 6.87e-02, 2.13e-02, 1.25e-01, 8.16e-02, 7.90e-02, 1.30e-01, 4.79e-02};
  Double_t cutiso_sumHyperTight1[9] = {
  7.92e+00, 6.85e+00, 7.87e+00, 6.77e+00, 4.47e+00, 5.28e+00, 6.57e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight1[9] = {
  5.20e+00, 3.93e+00, 3.88e+00, 4.10e+00, 2.40e+00, 2.43e+00, 3.49e+00, 3.94e+00, 3.01e+00};
  Double_t cutiso_sumoetlHyperTight1[9] = {
  4.18e+00, 3.12e+00, 3.44e+00, 3.25e+00, 1.77e+00, 2.06e+00, 2.83e+00, 3.12e+00, 1.43e+00};
  Double_t cutseeHyperTight1[9] = {
  1.05e-02, 1.04e-02, 1.01e-02, 3.24e-02, 2.80e-02, 2.85e-02, 9.67e-03, 2.61e-02, 2.95e-02};
  Double_t cutseelHyperTight1[9] = {
  1.04e-02, 1.03e-02, 1.01e-02, 3.04e-02, 2.74e-02, 2.78e-02, 9.58e-03, 2.54e-02, 2.83e-02};

  // HyperTight2 cuts
  Double_t cutdcotdistHyperTight2[9] = {
  1.15e-02, 1.07e-02, 4.01e-03, 2.97e-02, 2.85e-02, 3.10e-02, 9.34e-04, 3.40e-02, 2.82e-07};
  Double_t cutdetainHyperTight2[9] = {
  5.29e-03, 2.56e-03, 4.89e-03, 7.89e-03, 5.30e-03, 7.37e-03, 8.91e-03, 9.36e-03, 5.94e-03};
  Double_t cutdetainlHyperTight2[9] = {
  4.48e-03, 2.59e-03, 4.42e-03, 6.54e-03, 4.93e-03, 6.98e-03, 8.49e-03, 9.06e-03, -4.81e-03};
  Double_t cutdphiinHyperTight2[9] = {
  2.41e-02, 3.83e-02, 1.48e-01, 2.91e-02, 3.15e-02, 1.57e-01, 8.90e-02, 1.02e-01, 2.81e-01};
  Double_t cutdphiinlHyperTight2[9] = {
  2.13e-02, 3.79e-02, 1.25e-01, 2.24e-02, 3.69e-02, 1.64e-01, 9.99e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight2[9] = {
  1.03e+00, 9.95e-01, 1.03e+00, 1.01e+00, 9.46e-01, 9.03e-01, 9.97e-01, 1.14e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight2[9] = {
  1.50e+00, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight2[9] = {
  4.94e-02, 3.45e-02, 1.40e-01, 2.02e-01, 3.82e-02, 1.19e-01, 1.23e-01, 3.82e-01, 2.50e-01};
  Double_t cuthoelHyperTight2[9] = {
  4.04e-02, 3.42e-02, 1.31e-01, 1.85e-01, 3.01e-02, 1.27e-01, 2.27e-01, 3.80e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight2[9] = {
  1.14e-02, 1.38e-02, 5.29e-02, 1.87e-02, 1.31e-01, 8.63e-02, 7.74e-02, 1.04e-01, 2.42e-02};
  Double_t cutip_gsflHyperTight2[9] = {
  9.83e-03, 1.35e-02, 4.27e-02, 1.72e-02, 1.25e-01, 7.92e-02, 7.90e-02, 1.30e-01, 3.40e-02};
  Double_t cutiso_sumHyperTight2[9] = {
  6.40e+00, 5.77e+00, 6.54e+00, 5.22e+00, 3.86e+00, 4.63e+00, 6.31e+00, 1.02e+01, 1.78e+00};
  Double_t cutiso_sumoetHyperTight2[9] = {
  4.03e+00, 3.03e+00, 3.24e+00, 3.13e+00, 2.05e+00, 2.01e+00, 2.99e+00, 3.44e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight2[9] = {
  3.08e+00, 2.31e+00, 2.84e+00, 2.53e+00, 1.65e+00, 1.72e+00, 2.34e+00, 3.11e+00, 1.35e+00};
  Double_t cutseeHyperTight2[9] = {
  1.03e-02, 1.03e-02, 9.88e-03, 3.03e-02, 2.79e-02, 2.79e-02, 9.67e-03, 2.52e-02, 2.58e-02};
  Double_t cutseelHyperTight2[9] = {
  1.02e-02, 1.02e-02, 9.80e-03, 2.90e-02, 2.74e-02, 2.75e-02, 9.58e-03, 2.49e-02, 2.50e-02};

  // HyperTight3 cuts
  Double_t cutdcotdistHyperTight3[9] = {
  9.63e-03, 5.11e-03, 1.95e-04, 2.97e-02, 2.85e-02, 2.18e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight3[9] = {
  4.86e-03, 2.29e-03, 4.40e-03, 7.79e-03, 4.07e-03, 6.33e-03, 7.70e-03, 7.93e-03, 5.94e-03};
  Double_t cutdetainlHyperTight3[9] = {
  4.48e-03, 2.30e-03, 4.14e-03, 6.04e-03, 3.87e-03, 6.09e-03, 7.97e-03, 8.04e-03, -4.81e-03};
  Double_t cutdphiinHyperTight3[9] = {
  2.41e-02, 2.88e-02, 7.39e-02, 2.91e-02, 1.91e-02, 1.14e-01, 3.61e-02, 8.92e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight3[9] = {
  1.95e-02, 3.42e-02, 8.06e-02, 2.22e-02, 2.26e-02, 9.73e-02, 4.51e-02, 9.23e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight3[9] = {
  1.07e+00, 1.01e+00, 1.08e+00, 1.01e+00, 9.69e-01, 9.10e-01, 1.04e+00, 1.20e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight3[9] = {
  5.00e-01, 1.50e+00, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight3[9] = {
  3.52e-02, 3.45e-02, 1.33e-01, 1.88e-01, 2.72e-02, 1.19e-01, 9.28e-02, 2.46e-01, 2.50e-01};
  Double_t cuthoelHyperTight3[9] = {
  4.04e-02, 3.40e-02, 1.31e-01, 1.84e-01, 2.64e-02, 1.18e-01, 9.76e-02, 2.53e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight3[9] = {
  1.14e-02, 1.26e-02, 3.79e-02, 1.68e-02, 1.21e-01, 5.29e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight3[9] = {
  9.83e-03, 1.18e-02, 3.59e-02, 1.56e-02, 1.20e-01, 5.36e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight3[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.87e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight3[9] = {
  3.03e+00, 2.50e+00, 2.58e+00, 2.44e+00, 1.91e+00, 1.76e+00, 2.92e+00, 3.13e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight3[9] = {
  2.36e+00, 2.02e+00, 2.29e+00, 1.89e+00, 1.65e+00, 1.69e+00, 2.03e+00, 2.79e+00, 1.35e+00};
  Double_t cutseeHyperTight3[9] = {
  1.03e-02, 1.01e-02, 9.84e-03, 2.89e-02, 2.74e-02, 2.73e-02, 9.47e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight3[9] = {
  1.02e-02, 1.00e-02, 9.73e-03, 2.79e-02, 2.73e-02, 2.69e-02, 9.40e-03, 2.46e-02, 2.50e-02};

  // HyperTight4 cuts
  Double_t cutdcotdistHyperTight4[9] = {
  2.70e-04, 1.43e-04, 1.95e-04, 2.64e-03, 2.82e-02, 1.64e-02, 2.61e-05, 2.57e-02, 2.82e-07};
  Double_t cutdetainHyperTight4[9] = {
  2.44e-03, 1.67e-03, 2.26e-03, 3.43e-03, 3.51e-03, 3.52e-03, 2.98e-03, 4.79e-03, 5.94e-03};
  Double_t cutdetainlHyperTight4[9] = {
  2.34e-03, 1.29e-03, 2.30e-03, 3.30e-03, 3.61e-03, 3.84e-03, 2.53e-03, 3.66e-03, -4.81e-03};
  Double_t cutdphiinHyperTight4[9] = {
  8.44e-03, 5.21e-03, 2.18e-02, 1.39e-02, 7.82e-03, 1.52e-02, 2.59e-02, 3.87e-02, 2.81e-01};
  Double_t cutdphiinlHyperTight4[9] = {
  5.77e-03, 3.20e-03, 2.85e-02, 2.22e-02, 7.00e-03, 1.84e-02, 2.91e-02, 4.40e-02, 2.37e-01};
  Double_t cuteseedopcorHyperTight4[9] = {
  1.15e+00, 1.01e+00, 1.21e+00, 1.07e+00, 9.69e-01, 9.10e-01, 1.08e+00, 1.36e+00, 8.00e-01};
  Double_t cutfmishitsHyperTight4[9] = {
  5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, 5.00e-01, -5.00e-01};
  Double_t cuthoeHyperTight4[9] = {
  2.39e-02, 2.68e-02, 2.12e-02, 1.03e-01, 9.92e-03, 7.07e-02, 7.12e-02, 1.48e-01, 2.50e-01};
  Double_t cuthoelHyperTight4[9] = {
  2.87e-02, 1.94e-02, 2.16e-02, 5.68e-02, 1.35e-02, 4.04e-02, 7.98e-02, 1.50e-01, 1.32e-01};
  Double_t cutip_gsfHyperTight4[9] = {
  7.61e-03, 5.22e-03, 3.79e-02, 1.02e-02, 4.62e-02, 1.82e-02, 7.74e-02, 3.35e-02, 2.42e-02};
  Double_t cutip_gsflHyperTight4[9] = {
  7.81e-03, 4.25e-03, 3.08e-02, 1.04e-02, 2.35e-02, 2.45e-02, 7.90e-02, 2.88e-02, 3.40e-02};
  Double_t cutiso_sumHyperTight4[9] = {
  5.40e+00, 5.41e+00, 5.88e+00, 4.32e+00, 3.86e+00, 4.33e+00, 5.86e+00, 9.05e+00, 1.78e+00};
  Double_t cutiso_sumoetHyperTight4[9] = {
  2.53e+00, 2.10e+00, 1.87e+00, 1.84e+00, 1.79e+00, 1.61e+00, 2.53e+00, 1.98e+00, 2.76e+00};
  Double_t cutiso_sumoetlHyperTight4[9] = {
  2.28e+00, 2.02e+00, 2.04e+00, 1.69e+00, 1.65e+00, 1.61e+00, 2.03e+00, 1.82e+00, 1.35e+00};
  Double_t cutseeHyperTight4[9] = {
  9.99e-03, 9.61e-03, 9.65e-03, 2.75e-02, 2.61e-02, 2.64e-02, 9.18e-03, 2.44e-02, 2.58e-02};
  Double_t cutseelHyperTight4[9] = {
  9.66e-03, 9.69e-03, 9.58e-03, 2.73e-02, 2.66e-02, 2.66e-02, 8.64e-03, 2.46e-02, 2.50e-02};

    Double_t cutdcotdist[9];
    Double_t cutdetain[9];
    Double_t cutdetainl[9];
    Double_t cutdphiin[9];
    Double_t cutdphiinl[9];
    Double_t cuteseedopcor[9];
    Double_t cutfmishits[9];
    Double_t cuthoe[9];
    Double_t cuthoel[9];
    Double_t cutip_gsf[9];
    Double_t cutip_gsfl[9];
    Double_t cutiso_sum[9];
    Double_t cutiso_sumoet[9];
    Double_t cutiso_sumoetl[9];
    Double_t cutsee[9];
    Double_t cutseel[9];
 
    if(ElectronSelectionType == 300 || ElectronSelectionType == 301 || ElectronSelectionType == 302 ) {
      memcpy(cutdcotdist   ,cutdcotdistSuperTight   ,sizeof(cutdcotdistSuperTight));
      memcpy(cutdetain     ,cutdetainSuperTight     ,sizeof(cutdetainSuperTight));
      memcpy(cutdetainl    ,cutdetainlSuperTight    ,sizeof(cutdetainlSuperTight));
      memcpy(cutdphiin     ,cutdphiinSuperTight     ,sizeof(cutdphiinSuperTight));
      memcpy(cutdphiinl    ,cutdphiinlSuperTight    ,sizeof(cutdphiinlSuperTight));
      memcpy(cuteseedopcor ,cuteseedopcorSuperTight ,sizeof(cuteseedopcorSuperTight));
      memcpy(cutfmishits   ,cutfmishitsSuperTight   ,sizeof(cutfmishitsSuperTight));
      memcpy(cuthoe        ,cuthoeSuperTight	  ,sizeof(cuthoeSuperTight));
      memcpy(cuthoel       ,cuthoelSuperTight	  ,sizeof(cuthoelSuperTight));
      memcpy(cutip_gsf     ,cutip_gsfSuperTight     ,sizeof(cutip_gsfSuperTight));
      memcpy(cutip_gsfl    ,cutip_gsflSuperTight    ,sizeof(cutip_gsflSuperTight));
      memcpy(cutiso_sum    ,cutiso_sumSuperTight    ,sizeof(cutiso_sumSuperTight));
      memcpy(cutiso_sumoet ,cutiso_sumoetSuperTight ,sizeof(cutiso_sumoetSuperTight));
      memcpy(cutiso_sumoetl,cutiso_sumoetlSuperTight,sizeof(cutiso_sumoetlSuperTight));
      memcpy(cutsee        ,cutseeSuperTight	  ,sizeof(cutseeSuperTight));
      memcpy(cutseel       ,cutseelSuperTight	  ,sizeof(cutseelSuperTight));
    }
    else if(ElectronSelectionType == 310 || ElectronSelectionType == 311 || ElectronSelectionType == 312) {
      memcpy(cutdcotdist   ,cutdcotdistHyperTight1   ,sizeof(cutdcotdistHyperTight1));
      memcpy(cutdetain     ,cutdetainHyperTight1     ,sizeof(cutdetainHyperTight1));
      memcpy(cutdetainl    ,cutdetainlHyperTight1    ,sizeof(cutdetainlHyperTight1));
      memcpy(cutdphiin     ,cutdphiinHyperTight1     ,sizeof(cutdphiinHyperTight1));
      memcpy(cutdphiinl    ,cutdphiinlHyperTight1    ,sizeof(cutdphiinlHyperTight1));
      memcpy(cuteseedopcor ,cuteseedopcorHyperTight1 ,sizeof(cuteseedopcorHyperTight1));
      memcpy(cutfmishits   ,cutfmishitsHyperTight1   ,sizeof(cutfmishitsHyperTight1));
      memcpy(cuthoe        ,cuthoeHyperTight1	  ,sizeof(cuthoeHyperTight1));
      memcpy(cuthoel       ,cuthoelHyperTight1	  ,sizeof(cuthoelHyperTight1));
      memcpy(cutip_gsf     ,cutip_gsfHyperTight1     ,sizeof(cutip_gsfHyperTight1));
      memcpy(cutip_gsfl    ,cutip_gsflHyperTight1    ,sizeof(cutip_gsflHyperTight1));
      memcpy(cutiso_sum    ,cutiso_sumHyperTight1    ,sizeof(cutiso_sumHyperTight1));
      memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight1 ,sizeof(cutiso_sumoetHyperTight1));
      memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight1,sizeof(cutiso_sumoetlHyperTight1));
      memcpy(cutsee        ,cutseeHyperTight1	  ,sizeof(cutseeHyperTight1));
      memcpy(cutseel       ,cutseelHyperTight1	  ,sizeof(cutseelHyperTight1));
    }
    else if(ElectronSelectionType == 320 || ElectronSelectionType == 321 || ElectronSelectionType == 322) {
      memcpy(cutdcotdist   ,cutdcotdistHyperTight2   ,sizeof(cutdcotdistHyperTight2));
      memcpy(cutdetain     ,cutdetainHyperTight2     ,sizeof(cutdetainHyperTight2));
      memcpy(cutdetainl    ,cutdetainlHyperTight2    ,sizeof(cutdetainlHyperTight2));
      memcpy(cutdphiin     ,cutdphiinHyperTight2     ,sizeof(cutdphiinHyperTight2));
      memcpy(cutdphiinl    ,cutdphiinlHyperTight2    ,sizeof(cutdphiinlHyperTight2));
      memcpy(cuteseedopcor ,cuteseedopcorHyperTight2 ,sizeof(cuteseedopcorHyperTight2));
      memcpy(cutfmishits   ,cutfmishitsHyperTight2   ,sizeof(cutfmishitsHyperTight2));
      memcpy(cuthoe        ,cuthoeHyperTight2	  ,sizeof(cuthoeHyperTight2));
      memcpy(cuthoel       ,cuthoelHyperTight2	  ,sizeof(cuthoelHyperTight2));
      memcpy(cutip_gsf     ,cutip_gsfHyperTight2     ,sizeof(cutip_gsfHyperTight2));
      memcpy(cutip_gsfl    ,cutip_gsflHyperTight2    ,sizeof(cutip_gsflHyperTight2));
      memcpy(cutiso_sum    ,cutiso_sumHyperTight2    ,sizeof(cutiso_sumHyperTight2));
      memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight2 ,sizeof(cutiso_sumoetHyperTight2));
      memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight2,sizeof(cutiso_sumoetlHyperTight2));
      memcpy(cutsee        ,cutseeHyperTight2	  ,sizeof(cutseeHyperTight2));
      memcpy(cutseel       ,cutseelHyperTight2	  ,sizeof(cutseelHyperTight2));
    }
    else if(ElectronSelectionType == 330 || ElectronSelectionType == 331 || ElectronSelectionType == 332) {
      memcpy(cutdcotdist   ,cutdcotdistHyperTight3   ,sizeof(cutdcotdistHyperTight3));
      memcpy(cutdetain     ,cutdetainHyperTight3     ,sizeof(cutdetainHyperTight3));
      memcpy(cutdetainl    ,cutdetainlHyperTight3    ,sizeof(cutdetainlHyperTight3));
      memcpy(cutdphiin     ,cutdphiinHyperTight3     ,sizeof(cutdphiinHyperTight3));
      memcpy(cutdphiinl    ,cutdphiinlHyperTight3    ,sizeof(cutdphiinlHyperTight3));
      memcpy(cuteseedopcor ,cuteseedopcorHyperTight3 ,sizeof(cuteseedopcorHyperTight3));
      memcpy(cutfmishits   ,cutfmishitsHyperTight3   ,sizeof(cutfmishitsHyperTight3));
      memcpy(cuthoe        ,cuthoeHyperTight3	  ,sizeof(cuthoeHyperTight3));
      memcpy(cuthoel       ,cuthoelHyperTight3	  ,sizeof(cuthoelHyperTight3));
      memcpy(cutip_gsf     ,cutip_gsfHyperTight3     ,sizeof(cutip_gsfHyperTight3));
      memcpy(cutip_gsfl    ,cutip_gsflHyperTight3    ,sizeof(cutip_gsflHyperTight3));
      memcpy(cutiso_sum    ,cutiso_sumHyperTight3    ,sizeof(cutiso_sumHyperTight3));
      memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight3 ,sizeof(cutiso_sumoetHyperTight3));
      memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight3,sizeof(cutiso_sumoetlHyperTight3));
      memcpy(cutsee        ,cutseeHyperTight3	  ,sizeof(cutseeHyperTight3));
      memcpy(cutseel       ,cutseelHyperTight3	  ,sizeof(cutseelHyperTight3));
    }
    else if(ElectronSelectionType == 340 || ElectronSelectionType == 341 || ElectronSelectionType == 342) {
      memcpy(cutdcotdist   ,cutdcotdistHyperTight4   ,sizeof(cutdcotdistHyperTight4));
      memcpy(cutdetain     ,cutdetainHyperTight4     ,sizeof(cutdetainHyperTight4));
      memcpy(cutdetainl    ,cutdetainlHyperTight4    ,sizeof(cutdetainlHyperTight4));
      memcpy(cutdphiin     ,cutdphiinHyperTight4     ,sizeof(cutdphiinHyperTight4));
      memcpy(cutdphiinl    ,cutdphiinlHyperTight4    ,sizeof(cutdphiinlHyperTight4));
      memcpy(cuteseedopcor ,cuteseedopcorHyperTight4 ,sizeof(cuteseedopcorHyperTight4));
      memcpy(cutfmishits   ,cutfmishitsHyperTight4   ,sizeof(cutfmishitsHyperTight4));
      memcpy(cuthoe        ,cuthoeHyperTight4	  ,sizeof(cuthoeHyperTight4));
      memcpy(cuthoel       ,cuthoelHyperTight4	  ,sizeof(cuthoelHyperTight4));
      memcpy(cutip_gsf     ,cutip_gsfHyperTight4     ,sizeof(cutip_gsfHyperTight4));
      memcpy(cutip_gsfl    ,cutip_gsflHyperTight4    ,sizeof(cutip_gsflHyperTight4));
      memcpy(cutiso_sum    ,cutiso_sumHyperTight4    ,sizeof(cutiso_sumHyperTight4));
      memcpy(cutiso_sumoet ,cutiso_sumoetHyperTight4 ,sizeof(cutiso_sumoetHyperTight4));
      memcpy(cutiso_sumoetl,cutiso_sumoetlHyperTight4,sizeof(cutiso_sumoetlHyperTight4));
      memcpy(cutsee        ,cutseeHyperTight4	  ,sizeof(cutseeHyperTight4));
      memcpy(cutseel       ,cutseelHyperTight4	  ,sizeof(cutseelHyperTight4));
    }   

    int result = 0;
    
    const int ncuts = 10;
    std::vector<bool> cut_results(ncuts, false);
    
    float iso_sum = tkIso + ecalIso + hcalIso;
    if(fabs(ele->scEta)>1.5)
      iso_sum += (fabs(ele->scEta)-1.5)*1.09;
    float iso_sumoet = iso_sum*(40./ele->scEt);
    
    float eseedopincor = eSeedOverPin + fBrem;
    if(fBrem < 0)
      eseedopincor = eSeedOverPin;
    float dist = (TMath::Abs(ele->partnerDist)      == -9999.? 9999:TMath::Abs(ele->partnerDist));
    float dcot = (TMath::Abs(ele->partnerDeltaCot) == -9999.? 9999:TMath::Abs(ele->partnerDeltaCot));

    float dcotdistcomb = ((0.04 - std::max(dist, dcot)) > 0?(0.04 - std::max(dist, dcot)):0);

   for (int cut=0; cut<ncuts; cut++) {
    switch (cut) {
    case 0:
      cut_results[cut] = compute_cut(fabs(deltaEtaIn), ele->scEt, cutdetainl[cat], cutdetain[cat]);
      break;
    case 1:
      cut_results[cut] = compute_cut(fabs(deltaPhiIn), ele->scEt, cutdphiinl[cat], cutdphiin[cat]);
      break;
    case 2:
      cut_results[cut] = (eseedopincor > cuteseedopcor[cat]);
      break;
    case 3:
      cut_results[cut] = compute_cut(hOverE, ele->scEt, cuthoel[cat], cuthoe[cat]);
      break;
    case 4:
      cut_results[cut] = compute_cut(sigmaee, ele->scEt, cutseel[cat], cutsee[cat]);
      break;
    case 5:
      cut_results[cut] = compute_cut(iso_sumoet, ele->scEt, cutiso_sumoetl[cat], cutiso_sumoet[cat]);
      break;
    case 6:
      cut_results[cut] = (iso_sum < cutiso_sum[cat]);
      break;
    case 7:
      cut_results[cut] = compute_cut(fabs(ele->d0), ele->scEt, cutip_gsfl[cat], cutip_gsf[cat]);
      break;
    case 8:
      cut_results[cut] = (mishits < cutfmishits[cat]);
      break;
    case 9:
      cut_results[cut] = (dcotdistcomb < cutdcotdist[cat]);
      break;
    }
  }
    
   Bool_t passID = kFALSE;
   Bool_t passIso = kFALSE;
   Bool_t passIP = kFALSE;
   Bool_t passConversion = kFALSE;

   // ID part
   if (cut_results[0] & cut_results[1] & cut_results[2] & cut_results[3] & cut_results[4]) {
     passID = kTRUE;
     result = result + 1;
   }
   
   // ISO part
   if (cut_results[5] & cut_results[6]) {
     result = result + 2;
     passIso = kTRUE;
   }
   
   // IP part
   if (cut_results[7]) {
     result = result + 8;
     passIP = kTRUE;
   }

   // Conversion part
   if (cut_results[8] and cut_results[9]) {
     result = result + 4;
     passConversion = kTRUE;
   }
  


//     if (result != ele->passSuperTightId) { 
//       cout << a << "," << b << "," << c << " : " << ele->pt << " " << ele->eta << " " << ele->phi << "  : " << result << " =? " << ele->passSuperTightId << " < [cat=" << cat << "] " << endl;
//         }

   if (ElectronSelectionType == 300 || ElectronSelectionType == 310 || ElectronSelectionType == 320 || ElectronSelectionType == 330 || ElectronSelectionType == 340) {
     if (!(passID && passIso && passIP && passConversion)) pass = kFALSE;
   }
   if (ElectronSelectionType == 301 || ElectronSelectionType == 311 || ElectronSelectionType == 321 || ElectronSelectionType == 331 || ElectronSelectionType == 341) {
     if (!(passID)) pass = kFALSE;
   }
   if (ElectronSelectionType == 302 || ElectronSelectionType == 312 || ElectronSelectionType == 322 || ElectronSelectionType == 332 || ElectronSelectionType == 342) {
     if (!(passIso)) pass = kFALSE;
   }   

  }

//  //WW + CicTight
//   if (ElectronSelectionType == 310) {
//     //Barrel 
//     if (fabs(ele->eta) < 1.5) {
//       if (! ( (0==0)
//               && ele->sigiEtaiEta < 0.01 
//               && fabs(ele->deltaEtaIn) < 0.004
//               && fabs(ele->deltaPhiIn) < 0.06
//               && ele->HoverE < 0.04
//               && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
//               && ele->nExpHitsInner <= 0
//               && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//               && fabs(ele->d0) < 0.02
//             )
//         ) {
//         pass = kFALSE;
//       }      
//     }
//     //Endcap
//     else if (fabs(ele->eta) > 1.5) {
//       if (! (  (0==0)
//                && ele->sigiEtaiEta < 0.03
//                && fabs(ele->deltaEtaIn) < 0.007
//                && fabs(ele->deltaPhiIn) < 0.03
//                && ele->HoverE < 0.025
//                && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
//                && ele->nExpHitsInner <= 0
//                && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
//                && fabs(ele->d0) < 0.02
//             )
//         ) {
//         pass = kFALSE;
//       }
//     } else {
//       pass = kFALSE;
//       return pass;
//     }

//     if (ele->passCustomTightId == kFALSE) {
//         pass = kFALSE;      
//     }
//   }






 //likelihood VeryLoose
  if (ElectronSelectionType == 401) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.221) pass = kFALSE;
      } else {
        if (likelihood <= 0.180) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.103) pass = kFALSE;
      } else {
        if (likelihood <= 0.355) pass = kFALSE;
      }    
    }

  }

 //likelihood Loose
  if (ElectronSelectionType == 402) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.370) pass = kFALSE;
      } else {
        if (likelihood <= 0.390) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.433) pass = kFALSE;
      } else {
        if (likelihood <= 0.745) pass = kFALSE;
      }    
    }

  }

 //likelihood Medium
  if (ElectronSelectionType == 403) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.592) pass = kFALSE;
      } else {
        if (likelihood <= 0.732) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.741) pass = kFALSE;
      } else {
        if (likelihood <= 0.919) pass = kFALSE;
      }    
    }

  }



 //likelihood Tight
  if (ElectronSelectionType == 404) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.820) pass = kFALSE;
      } else {
        if (likelihood <= 0.925) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.930) pass = kFALSE;
      } else {
        if (likelihood <= 0.979) pass = kFALSE;
      }    
    }

  }


 //likelihood SiTight1 
  if (ElectronSelectionType == 410) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.950) pass = kFALSE;
      } else {
        if (likelihood <= 0.950) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.950) pass = kFALSE;
      } else {
        if (likelihood <= 0.950) pass = kFALSE;
      }    
    }

  }


//likelihood SiTight2
  if (ElectronSelectionType == 411) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.975) pass = kFALSE;
      } else {
        if (likelihood <= 0.975) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.975) pass = kFALSE;
      } else {
        if (likelihood <= 0.975) pass = kFALSE;
      }    
    }

  }


//likelihood SiTight3
  if (ElectronSelectionType == 412) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.990) pass = kFALSE;
      } else {
        if (likelihood <= 0.990) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.990) pass = kFALSE;
      } else {
        if (likelihood <= 0.990) pass = kFALSE;
      }    
    }

  }


//likelihood SiTight3
  if (ElectronSelectionType == 413) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.995) pass = kFALSE;
      } else {
        if (likelihood <= 0.995) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.995) pass = kFALSE;
      } else {
        if (likelihood <= 0.995) pass = kFALSE;
      }    
    }

  }


//likelihood SiTight3
  if (ElectronSelectionType == 414) {
    Double_t relIsoCut = 0.1;
    if (ele->pt >= 10 && ele->pt < 15) relIsoCut = 0.05;
    if (ele->pt >= 15 && ele->pt < 20) relIsoCut = 0.07;

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < relIsoCut
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < relIsoCut
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.999) pass = kFALSE;
      } else {
        if (likelihood <= 0.999) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.999) pass = kFALSE;
      } else {
        if (likelihood <= 0.999) pass = kFALSE;
      }    
    }

  }


  if (ElectronSelectionType == 423) {
     //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < 0.1
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < 0.1
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.995) pass = kFALSE;
      } else {
        if (likelihood <= 0.995) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.995) pass = kFALSE;
      } else {
        if (likelihood <= 0.995) pass = kFALSE;
      }    
    }


    Int_t originalvalue = ele->passSuperTightId;
    Int_t bit0;
    Int_t bit1;
    Int_t bit2;
    Int_t bit3;
    bit0 = originalvalue % 2;
    Int_t tmp0 = floor(double(originalvalue) / 2.0);
    bit1 = tmp0 % 2;
    Int_t tmp1 = floor(double(tmp0) / 2.0);
    bit2 = tmp1 % 2;
    Int_t tmp2 = floor(double(tmp1) / 2.0);
    bit3 = tmp2 % 2;

    if (!(bit1 == 1 )) {
        pass = kFALSE;      
    }

  }


//likelihood SiTight4
  if (ElectronSelectionType == 424) {

    //Barrel 
    if (fabs(ele->eta) < 1.5) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && (ele->trkIso03 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt < 0.1
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
               && (ele->trkIso03 + ele->emIso04  + ele->hadIso04) / ele->pt < 0.1
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

    if (fabs(ele->eta) < 1.5) {
      if (ele->nBrem == 0) {
        if (likelihood <= 0.999) pass = kFALSE;
      } else {
        if (likelihood <= 0.999) pass = kFALSE;
      }
    } else if (fabs(ele->eta) > 1.5) {
        if (ele->nBrem == 0) {
        if (likelihood <= 0.999) pass = kFALSE;
      } else {
        if (likelihood <= 0.999) pass = kFALSE;
      }    
    }


    Int_t originalvalue = ele->passSuperTightId;
    Int_t bit0;
    Int_t bit1;
    Int_t bit2;
    Int_t bit3;
    bit0 = originalvalue % 2;
    Int_t tmp0 = floor(double(originalvalue) / 2.0);
    bit1 = tmp0 % 2;
    Int_t tmp1 = floor(double(tmp0) / 2.0);
    bit2 = tmp1 % 2;
    Int_t tmp2 = floor(double(tmp1) / 2.0);
    bit3 = tmp2 % 2;

    if (!(bit1 == 1 )) {
        pass = kFALSE;      
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


Bool_t passMuonCuts(const mithep::TMuon *mu, Double_t PileupEnergyDensity) {
  
  Bool_t pass = kTRUE;

  Double_t PileupIsoEnergy = PileupEnergyDensity * 3.1416 * pow(0.3,2);

  Double_t isoCut = 0.13;
  if (mu->pt < 20) isoCut = 0.08;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + TMath::Max(mu->emIso03 + mu->hadIso03 - PileupIsoEnergy,0.0)) / mu->pt < isoCut

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
