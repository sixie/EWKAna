//root -l EWKAna/Hww/LeptonSelection/MakeBkgElectronDistributions.C+\(\)
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
#include <TLegend.h>     
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
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void MakeBkgElectron(const string inputFilename,
                     const string label);

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

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


//=== MAIN MACRO =================================================================================================
void MakeBkgElectronDistributions() {

//   MakeBkgElectron("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr-v1_EleFakeRateTriggerSkim.root","Ele10Jet30");
  MakeBkgElectron("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-ww2l-v1g1-pu_noskim_normalized.root","WW");


}



void MakeBkgElectron(const string inputFilename,
                     const string Label) 
{  

  string label = Label; 
  
  if (Label != "") label = "_"+Label;
  
  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Double_t jetPtThreshold = 30;

  //*****************************************************************************************
  //Histogram
  //*****************************************************************************************
  vector<Double_t> ElectronIsolationEfficiencyNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyNumerator;
  vector<Double_t> ElectronIsolationEfficiencyDenominator;
  for(UInt_t i=0 ; i<20; ++i) {
    ElectronIsolationEfficiencyNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyNumerator.push_back(0.0);
    ElectronIsolationEfficiencyDenominator.push_back(0.0);
  }

  TH1F *Electron_relIso_Barrel = new TH1F((string("Electron_relIso_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIsoL1Corrected_Barrel = new TH1F((string("Electron_relIsoL1Corrected_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIso04_Barrel = new TH1F((string("Electron_relIso04_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIso04L1Corrected_Barrel = new TH1F((string("Electron_relIso04L1Corrected_Barrel_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIso_Endcap = new TH1F((string("Electron_relIso_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIsoL1Corrected_Endcap = new TH1F((string("Electron_relIsoL1Corrected_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIso04_Endcap = new TH1F((string("Electron_relIso04_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
  TH1F *Electron_relIso04L1Corrected_Endcap = new TH1F((string("Electron_relIso04L1Corrected_Endcap_BkgData")+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);


  vector<TH1F*>   Rho;
  vector<TH1F*>   Electron_caloIso_Barrel;
  vector<TH1F*>   Electron_caloIso_Endcap;
  vector<TH1F*>   Electron_trkIso_Barrel;
  vector<TH1F*>   Electron_trkIso_Endcap;
  vector<TH1F*>   Electron_caloIso04_Barrel;
  vector<TH1F*>   Electron_caloIso04_Endcap;
  vector<TH1F*>   Electron_trkIso04_Barrel;
  vector<TH1F*>   Electron_trkIso04_Endcap;

  for (int n=0; n < 20; ++n) {     
    TH1F *tmpRho = new TH1F((string("RhoElectron_") + IntToString(n) + label).c_str(), "; Rho; Number of Events ", 2000, -20, 20);
    TH1F *tmpElectron_caloIso_Barrel = new TH1F((string("Electron_caloIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso_Endcap = new TH1F((string("Electron_caloIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso_Barrel = new TH1F((string("Electron_trkIso_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso_Endcap = new TH1F((string("Electron_trkIso_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso04_Barrel = new TH1F((string("Electron_caloIso04_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_caloIso04_Endcap = new TH1F((string("Electron_caloIso04_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso04_Barrel = new TH1F((string("Electron_trkIso04_Barrel_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);
    TH1F *tmpElectron_trkIso04_Endcap = new TH1F((string("Electron_trkIso04_Endcap_")+ IntToString(n) +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200, 200.0);

    tmpRho->StatOverflows(kTRUE);
    tmpElectron_caloIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_caloIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_trkIso_Barrel->StatOverflows(kTRUE);
    tmpElectron_trkIso_Endcap->StatOverflows(kTRUE);
    tmpElectron_caloIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_caloIso04_Endcap->StatOverflows(kTRUE);
    tmpElectron_trkIso04_Barrel->StatOverflows(kTRUE);
    tmpElectron_trkIso04_Endcap->StatOverflows(kTRUE);

    Rho.push_back(tmpRho);
    Electron_caloIso_Barrel.push_back(tmpElectron_caloIso_Barrel);
    Electron_trkIso_Barrel.push_back(tmpElectron_trkIso_Barrel);
    Electron_caloIso04_Barrel.push_back(tmpElectron_caloIso04_Barrel);       
    Electron_trkIso04_Barrel.push_back(tmpElectron_trkIso04_Barrel);       
    Electron_caloIso_Endcap.push_back(tmpElectron_caloIso_Endcap);
    Electron_trkIso_Endcap.push_back(tmpElectron_trkIso_Endcap);
    Electron_caloIso04_Endcap.push_back(tmpElectron_caloIso04_Endcap);       
    Electron_trkIso04_Endcap.push_back(tmpElectron_trkIso04_Endcap);       
  }

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
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
   rlrm.AddJSONFile("260311_tmp_JSON.txt"); 
   hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

  cout << "Total Events: " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
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
	
    Double_t PUIsolationEnergy = info->PileupEnergyDensity * 3.14159 * pow(0.3,2) * 1.0 ;
    Double_t PUIsolationEnergy04Cone = info->PileupEnergyDensity * 3.14159 * pow(0.4,2) * 1.0 ;
    Double_t PUIsolationEnergy05Cone = info->PileupEnergyDensity * 3.14159 * pow(0.5,2) * 1.0 ;
    Int_t NVertex = info->nPV0; if (NVertex > 19) NVertex = 19;

    //********************************************************
    // Event Selection Cuts
    //********************************************************
//     if (met.Pt() > 10) continue;
//     if (muonArr->GetEntries() > 1) continue;

    Double_t tempLeadingJetPt = 0;
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
      if( jet->pt > tempLeadingJetPt) tempLeadingJetPt = jet->pt;
    }
    
 



    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);

 //      //pass HLT selection
//       if (!passHLT(info->triggerBits, info->runNum, SelectionType)) continue;
      
      //pass event selection
      Bool_t passJetSelection = kFALSE;
      for(Int_t i=0; i<jetArr->GetEntries(); i++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
        
        if (jet->pt > jetPtThreshold &&
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 0.5) {
          passJetSelection = kTRUE;
          break;
        }
      }
//       if (!passJetSelection) continue;

      if (!passElectronDenominatorCuts(ele)) continue;


      Double_t relIso = (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt;
      if (fabs(ele->eta) < 1.5) relIso = (ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt;
      Double_t relIso04 = (ele->trkIso04 + TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04) / ele->pt;
      if (fabs(ele->eta) < 1.5) relIso04 = (ele->trkIso04 + ele->emIso04 + ele->hadIso04) / ele->pt;
      Double_t relIsoL1Corrected = (ele->trkIso03 + TMath::Max(TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03 - PUIsolationEnergy, 0.0)) / ele->pt;
      if (fabs(ele->eta) < 1.5) relIso = (ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt;
      Double_t relIso04L1Corrected = (ele->trkIso04 + TMath::Max(TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04 - PUIsolationEnergy04Cone,0.0)) / ele->pt;
      if (fabs(ele->eta) < 1.5) relIso04 = (ele->trkIso04 + ele->emIso04 + ele->hadIso04) / ele->pt;

      Rho[NVertex]->Fill(info->PileupEnergyDensity);
      if (fabs(ele->eta) < 1.5) {
        Electron_caloIso_Barrel[NVertex]->Fill(TMath::Min(double( TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03), double(199.99)));
        Electron_trkIso_Barrel[NVertex]->Fill(TMath::Min(double( ele->trkIso03 ), double(199.99)));
        Electron_caloIso04_Barrel[NVertex]->Fill(TMath::Min( double(TMath::Max(ele->emIso04 - 1.0, 0.0) + ele->hadIso04), double(199.99)));
        Electron_trkIso04_Barrel[NVertex]->Fill(TMath::Min( double(ele->trkIso04 ) , double(199.99)));
        Electron_relIso_Barrel->Fill( TMath::Min(double(relIso), double(1.99)) );
        Electron_relIsoL1Corrected_Barrel->Fill(TMath::Min(double(relIsoL1Corrected), double(1.99)));
        Electron_relIso04_Barrel->Fill(TMath::Min(double(relIso04), double(1.99)) );
        Electron_relIso04L1Corrected_Barrel->Fill(TMath::Min(double(relIso04L1Corrected), double(1.99)));
      } else {
        Electron_caloIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->emIso03  + ele->hadIso03), double(199.99)));
        Electron_trkIso_Endcap[NVertex]->Fill(TMath::Min(double( ele->trkIso03 ), double(199.99)));
        Electron_caloIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->emIso04  + ele->hadIso04), double(199.99)));
        Electron_trkIso04_Endcap[NVertex]->Fill(TMath::Min( double(ele->trkIso04 ) , double(199.99)));
        Electron_relIso_Endcap->Fill( TMath::Min(double(relIso), double(1.99)) );
        Electron_relIsoL1Corrected_Endcap->Fill(TMath::Min(double(relIsoL1Corrected), double(1.99)));
        Electron_relIso04_Endcap->Fill(TMath::Min(double(relIso04), double(1.99)) );
        Electron_relIso04L1Corrected_Endcap->Fill(TMath::Min(double(relIso04L1Corrected), double(1.99)));
     }


      ElectronIsolationEfficiencyDenominator[NVertex]++;  
      if ((ele->trkIso03 + ele->emIso03 + ele->hadIso03) / ele->pt < 0.10) {
//                 if ((ele->trkIso04 + ele->emIso04 + ele->hadIso04) / ele->pt < 0.15) {
        ElectronIsolationEfficiencyNumerator[NVertex]++;
      }              
      if ((ele->trkIso03 + TMath::Max(ele->emIso03 + ele->hadIso03 - PUIsolationEnergy,0.0)) / ele->pt < 0.10) {
//                 if ((ele->trkIso04 + TMath::Max(ele->emIso04 + ele->hadIso04 - PUIsolationEnergy04Cone,0.0)) / ele->pt < 0.15) {
        ElectronL1CorrectedIsolationEfficiencyNumerator[NVertex]++;
      }
      
    }

  } //end loop over data     


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  
  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================
  const int nPoints = 20;
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.0;     
  }

  double ElectronIsolationEff[nPoints];
  double ElectronIsolationEffErrLow[nPoints];
  double ElectronIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationEff[i] = ratio;
    ElectronIsolationEffErrLow[i] = errLow;
    ElectronIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronIsolationEff " << i << " : " << ElectronIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationEff, NPileupError, NPileupError, ElectronIsolationEffErrLow, ElectronIsolationEffErrHigh);
  ElectronIsolationEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedIsolationEff[nPoints];
  double ElectronL1CorrectedIsolationEffErrLow[nPoints];
  double ElectronL1CorrectedIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedIsolationEff[i] = ratio;
    ElectronL1CorrectedIsolationEffErrLow[i] = errLow;
    ElectronL1CorrectedIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedIsolationEff " << i << " : " << ElectronL1CorrectedIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronL1CorrectedIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationEff, NPileupError, NPileupError,ElectronL1CorrectedIsolationEffErrLow,ElectronL1CorrectedIsolationEffErrHigh  );
  ElectronL1CorrectedIsolationEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  Double_t ymin = 0.0;
  Double_t ymax = 1.1;

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  legend = new TLegend(0.20,0.55,0.43,0.70);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 0.4;
  legend->AddEntry(ElectronIsolationEffVsNPileup, "NoCorrection", "LP");
  legend->AddEntry(ElectronL1CorrectedIsolationEffVsNPileup, "FastJetCorrected", "LP");
  ElectronIsolationEffVsNPileup->SetTitle("");
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronIsolationEffVsNPileup->Draw("AP");
  ElectronL1CorrectedIsolationEffVsNPileup->Draw("Psame");
  ElectronIsolationEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);

  legend->Draw();
  cv->SaveAs(("ElectronIsolationEfficiency_BkgData_vs_NVertices"+label+".gif").c_str());




  //*****************************************************************************************
  //Save Efficiency Plots
  //*****************************************************************************************
  TFile *file = new TFile("HwwSelectionPlots_LeptonEfficiency.root", "UPDATE");
  file->cd();

  file->WriteTObject(Electron_relIso_Barrel, Electron_relIso_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIsoL1Corrected_Barrel, Electron_relIsoL1Corrected_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIso04_Barrel, Electron_relIso04_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIso04L1Corrected_Barrel, Electron_relIso04L1Corrected_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIso_Endcap, Electron_relIso_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIsoL1Corrected_Endcap, Electron_relIsoL1Corrected_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIso04_Endcap, Electron_relIso04_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(Electron_relIso04L1Corrected_Endcap, Electron_relIso04L1Corrected_Endcap->GetName(), "WriteDelete");

   for (int n=0; n < 20 ; ++n) {
      file->WriteTObject(Rho[n],Rho[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Barrel[n],Electron_caloIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Barrel[n],Electron_trkIso_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Barrel[n],Electron_caloIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Barrel[n],Electron_trkIso04_Barrel[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Endcap[n],Electron_caloIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Endcap[n],Electron_trkIso_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Endcap[n],Electron_caloIso04_Endcap[n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Endcap[n],Electron_trkIso04_Endcap[n]->GetName(), "WriteDelete") ;
    }


  file->Close();
  delete file;

    
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType) {


  Bool_t pass = kFALSE;


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
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
//             && ele->HoverE < 0.04
//             && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
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
//              && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
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
