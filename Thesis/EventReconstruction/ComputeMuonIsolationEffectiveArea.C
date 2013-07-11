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
#include "TLegend.h"
#include "TProfile.h"
#include "TPaveLabel.h"
#include "TF1.h"

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"

#endif

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


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
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


//*************************************************************************************************
//Make Projection Graph
//*************************************************************************************************
TGraphAsymmErrors *MakeMeanVsNPUGraph(TH2F *hist2D, string name, string axisLabel = "") {

  const UInt_t nPoints = hist2D->GetYaxis()->GetNbins();
  double NPU[nPoints];
  double NPUErr[nPoints];
  double MeanValue[nPoints];
  double MeanErr[nPoints];

//   cout << "Npoints: " << nPoints << endl;

  for(UInt_t b=0; b < nPoints; ++b) {
    TH1F *Projection = (TH1F*)hist2D->ProjectionX("_px", b,b+1);
    NPU[b] = (hist2D->GetYaxis()->GetBinUpEdge(b) + hist2D->GetYaxis()->GetBinLowEdge(b)) / 2;
    NPUErr[b] = 0.1;
    MeanValue[b] = Projection->GetMean();
    MeanErr[b] = Projection->GetMeanError();

//     cout << b << " : " << NPU[b] << " : " << MeanValue[b] << " " << MeanErr[b] << endl;

  }

  TGraphAsymmErrors *Graph = new TGraphAsymmErrors (nPoints,  NPU, MeanValue, NPUErr, NPUErr, MeanErr, MeanErr);
  Graph->SetName(name.c_str());
  Graph->SetTitle("");
  Graph->SetMarkerColor(kBlack);
  Graph->GetXaxis()->SetTitleOffset(1.02);
  Graph->GetXaxis()->SetTitle("Number of Pileup Events");
  Graph->GetYaxis()->SetTitleOffset(1.05);
  Graph->GetYaxis()->SetTitle(axisLabel.c_str());

  return Graph;

}


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void FillIsolationHistograms(string MuonFile, string Label, Int_t EtaBin = -1)
{  

  string label = "";
  if (Label != "") label = "_" + Label;

   if (EtaBin == 0) {
    label = label + "_Eta0To1";
  }
  if (EtaBin == 1) {
    label = label + "_Eta1To1p479";
  }
  if (EtaBin == 2) {
    label = label + "_Eta1p479To2";
  }
  if (EtaBin == 3) {
    label = label + "_Eta2To2p25";
  }
  if (EtaBin == 4) {
    label = label + "_Eta2p25To2p5";
  }

 //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH2F *Mu_ChargedIso03_VS_NPU = new TH2F(("Mu_ChargedIso03_VS_NPU"+label).c_str(), "; ChargedIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_NeutralIso03_VS_NPU = new TH2F(("Mu_ChargedIso03_VS_NPU"+label).c_str(), "; ChargedIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_ChargedIso04_VS_NPU = new TH2F(("Mu_ChargedIso04_VS_NPU"+label).c_str(), "; ChargedIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_NeutralIso04_VS_NPU = new TH2F(("Mu_ChargedIso04_VS_NPU"+label).c_str(), "; ChargedIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadEnergy_VS_NPU = new TH2F(("Mu_HadEnergy_VS_NPU"+label).c_str(), "; HadEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HoEnergy_VS_NPU = new TH2F(("Mu_HoEnergy_VS_NPU"+label).c_str(), "; HoEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EmEnergy_VS_NPU = new TH2F(("Mu_EmEnergy_VS_NPU"+label).c_str(), "; EmEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadS9Energy_VS_NPU = new TH2F(("Mu_HadS9Energy_VS_NPU"+label).c_str(), "; HadS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HoS9Energy_VS_NPU = new TH2F(("Mu_HoS9Energy_VS_NPU"+label).c_str(), "; HoS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EmS9Energy_VS_NPU = new TH2F(("Mu_EmS9Energy_VS_NPU"+label).c_str(), "; EmS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_TrkIso03_VS_NPU = new TH2F(("Mu_TrkIso03_VS_NPU"+label).c_str(), "; TrkIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EMIso03_VS_NPU = new TH2F(("Mu_EMIso03_VS_NPU"+label).c_str(), "; EMIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadIso03_VS_NPU = new TH2F(("Mu_HadIso03_VS_NPU"+label).c_str(), "; HadIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_TrkIso05_VS_NPU = new TH2F(("Mu_TrkIso05_VS_NPU"+label).c_str(), "; TrkIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EMIso05_VS_NPU = new TH2F(("Mu_EMIso05_VS_NPU"+label).c_str(), "; EMIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadIso05_VS_NPU = new TH2F(("Mu_HadIso05_VS_NPU"+label).c_str(), "; HadIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_Rho_VS_NPU = new TH2F(("Mu_Rho_VS_NPU"+label).c_str(), "; Rho ; Number of Pileup Events ; Number of Events ",  1000, 0 , 200, 50, -0.5, 49.5);

  TH2F *Mu_ChargedIso03_VS_NVtx = new TH2F(("Mu_ChargedIso03_VS_NVtx"+label).c_str(), "; ChargedIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_NeutralIso03_VS_NVtx = new TH2F(("Mu_ChargedIso03_VS_NVtx"+label).c_str(), "; ChargedIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_ChargedIso04_VS_NVtx = new TH2F(("Mu_ChargedIso04_VS_NVtx"+label).c_str(), "; ChargedIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_NeutralIso04_VS_NVtx = new TH2F(("Mu_ChargedIso04_VS_NVtx"+label).c_str(), "; ChargedIso04 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadEnergy_VS_NVtx = new TH2F(("Mu_HadEnergy_VS_NVtx"+label).c_str(), "; HadEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HoEnergy_VS_NVtx = new TH2F(("Mu_HoEnergy_VS_NVtx"+label).c_str(), "; HoEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EmEnergy_VS_NVtx = new TH2F(("Mu_EmEnergy_VS_NVtx"+label).c_str(), "; EmEnergy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadS9Energy_VS_NVtx = new TH2F(("Mu_HadS9Energy_VS_NVtx"+label).c_str(), "; HadS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HoS9Energy_VS_NVtx = new TH2F(("Mu_HoS9Energy_VS_NVtx"+label).c_str(), "; HoS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EmS9Energy_VS_NVtx = new TH2F(("Mu_EmS9Energy_VS_NVtx"+label).c_str(), "; EmS9Energy ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_TrkIso03_VS_NVtx = new TH2F(("Mu_TrkIso03_VS_NVtx"+label).c_str(), "; TrkIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EMIso03_VS_NVtx = new TH2F(("Mu_EMIso03_VS_NVtx"+label).c_str(), "; EMIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadIso03_VS_NVtx = new TH2F(("Mu_HadIso03_VS_NVtx"+label).c_str(), "; HadIso03 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_TrkIso05_VS_NVtx = new TH2F(("Mu_TrkIso05_VS_NVtx"+label).c_str(), "; TrkIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_EMIso05_VS_NVtx = new TH2F(("Mu_EMIso05_VS_NVtx"+label).c_str(), "; EMIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_HadIso05_VS_NVtx = new TH2F(("Mu_HadIso05_VS_NVtx"+label).c_str(), "; HadIso05 ; Number of Pileup Events ; Number of Events ",  1000, 0 , 1000, 50, -0.5, 49.5);
  TH2F *Mu_Rho_VS_NVtx = new TH2F(("Mu_Rho_VS_NVtx"+label).c_str(), "; Rho ; Number of Pileup Events ; Number of Events ",  1000, 0 , 200, 50, -0.5, 49.5);

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

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
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  Int_t NEvents = 0;
  Bool_t isMC = kTRUE;

  vector<string> inputfiles;
  if (MuonFile == "Fall11ZmmMC") {
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0000.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0004.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0012.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0015.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0016.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/mc025/AllNtuple_HWWNtuple_f11-zmmm20-powheg-v14b-pu_noskim_0018.root");
  } else if (MuonFile == "Data2011") {
    isMC = kFALSE;
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-m10-v1_TightPlusRecoSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-pr-v4_TightPlusRecoSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-a05-v1_TightPlusRecoSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-o03-v1_TightPlusRecoSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-smu-pr-v1_TightPlusRecoSkim.root");    
  } else {
    inputfiles.push_back(MuonFile);
  }

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);


      //********************************************************
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      if (isMC) {
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
          
          //********************************************************
          // Select MC Truth Muons
          //********************************************************
          
          if (!((UInt_t(abs(max(0,mu->isMCReal))) & 2) == 2)) continue;

          if (EtaBin == 0) {
            if (!(fabs(mu->eta) < 1.0)) continue;
          }
          if (EtaBin == 1) {
            if (!(fabs(mu->eta) >= 1.0 && fabs(mu->eta) < 1.479)) continue;
          }
          if (EtaBin == 2) {
             if (!(fabs(mu->eta) >= 1.479 && fabs(mu->eta) < 2.0)) continue;
          }
          if (EtaBin == 3) {
            if (!(fabs(mu->eta) >= 2.0 && fabs(mu->eta) < 2.25)) continue;
          }
          if (EtaBin == 4) {
            if (!(fabs(mu->eta) >= 2.25 && fabs(mu->eta) < 2.5)) continue;
          }

          Mu_ChargedIso03_VS_NPU->Fill(mu->ChargedIso03,info->nPUEvents);
          Mu_NeutralIso03_VS_NPU->Fill(mu->NeutralIso03_01Threshold,info->nPUEvents); 
          Mu_ChargedIso04_VS_NPU->Fill(mu->ChargedIso04,info->nPUEvents); 
          Mu_NeutralIso04_VS_NPU->Fill(mu->NeutralIso04_01Threshold,info->nPUEvents); 
          Mu_HadEnergy_VS_NPU->Fill(mu->HadEnergy,info->nPUEvents); 
          Mu_EmEnergy_VS_NPU->Fill(mu->EmEnergy,info->nPUEvents); 
          Mu_HoEnergy_VS_NPU->Fill(mu->HoEnergy,info->nPUEvents); 
          Mu_HadS9Energy_VS_NPU->Fill(mu->HadS9Energy,info->nPUEvents); 
          Mu_EmS9Energy_VS_NPU->Fill(mu->EmS9Energy,info->nPUEvents); 
          Mu_HoS9Energy_VS_NPU->Fill(mu->HoS9Energy,info->nPUEvents); 
          Mu_TrkIso03_VS_NPU->Fill(mu->trkIso03,info->nPUEvents); 
          Mu_EMIso03_VS_NPU->Fill(mu->emIso03,info->nPUEvents); 
          Mu_HadIso03_VS_NPU->Fill(mu->hadIso03,info->nPUEvents); 
          Mu_TrkIso05_VS_NPU->Fill(mu->trkIso05,info->nPUEvents); 
          Mu_EMIso05_VS_NPU->Fill(mu->emIso05,info->nPUEvents); 
          Mu_HadIso05_VS_NPU->Fill(mu->hadIso05,info->nPUEvents); 
          Mu_Rho_VS_NPU->Fill(info->PileupEnergyDensity,info->nPUEvents);
  
          Mu_ChargedIso03_VS_NVtx->Fill(mu->ChargedIso03,info->nPV0);
          Mu_NeutralIso03_VS_NVtx->Fill(mu->NeutralIso03_01Threshold,info->nPV0); 
          Mu_ChargedIso04_VS_NVtx->Fill(mu->ChargedIso04,info->nPV0); 
          Mu_NeutralIso04_VS_NVtx->Fill(mu->NeutralIso04_01Threshold,info->nPV0); 
          Mu_HadEnergy_VS_NVtx->Fill(mu->HadEnergy,info->nPV0); 
          Mu_EmEnergy_VS_NVtx->Fill(mu->EmEnergy,info->nPV0); 
          Mu_HoEnergy_VS_NVtx->Fill(mu->HoEnergy,info->nPV0); 
          Mu_HadS9Energy_VS_NVtx->Fill(mu->HadS9Energy,info->nPV0); 
          Mu_EmS9Energy_VS_NVtx->Fill(mu->EmS9Energy,info->nPV0); 
          Mu_HoS9Energy_VS_NVtx->Fill(mu->HoS9Energy,info->nPV0); 
          Mu_TrkIso03_VS_NVtx->Fill(mu->trkIso03,info->nPV0); 
          Mu_EMIso03_VS_NVtx->Fill(mu->emIso03,info->nPV0); 
          Mu_HadIso03_VS_NVtx->Fill(mu->hadIso03,info->nPV0); 
          Mu_TrkIso05_VS_NVtx->Fill(mu->trkIso05,info->nPV0); 
          Mu_EMIso05_VS_NVtx->Fill(mu->emIso05,info->nPV0); 
          Mu_HadIso05_VS_NVtx->Fill(mu->hadIso05,info->nPV0); 
          Mu_Rho_VS_NVtx->Fill(info->PileupEnergyDensity,info->nPV0);

        } //loop over muons
      } else {

        Bool_t FilledRho = kFALSE;

        //Require single mu triggers
        ULong_t trigger  = kHLT_IsoMu17 | kHLT_IsoMu24 | kHLT_IsoMu30 | kHLT_Mu15 | kHLT_Mu30;
        if(!(info->triggerBits & trigger)) continue;      
        
        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *tag = (mithep::TMuon*)((*muonArr)[i]);
          
          if(tag->pt        < 20)  continue;
          if(fabs(tag->eta) > 2.4) continue;
          if(!passMuonID(tag))     continue;
          
          if(!((info->triggerBits & kHLT_IsoMu17) && (tag->hltMatchBits & kHLTObject_IsoMu17)) &&
             !((info->triggerBits & kHLT_IsoMu24) && (tag->hltMatchBits & kHLTObject_IsoMu24)) &&
             !((info->triggerBits & kHLT_IsoMu30) && (tag->hltMatchBits & kHLTObject_IsoMu30)) &&
             !((info->triggerBits & kHLT_Mu15)    && (tag->hltMatchBits & kHLTObject_Mu15)) &&
             !((info->triggerBits & kHLT_Mu30)    && (tag->hltMatchBits & kHLTObject_Mu30)) ) 
            continue;
          
          const Double_t m = 0.105658369;
          TLorentzVector vtag;
          vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
          
          for(Int_t j=0; j<muonArr->GetEntriesFast(); j++) {
            if(i==j) continue;
            
            const mithep::TMuon *probe = (mithep::TMuon*)((*muonArr)[j]);
            if(probe->q == tag->q) continue;
            
            if(!(probe->typeBits & kGlobal) && !(probe->typeBits & kTracker)) continue;
            
//             if (!(probe->pt > 20 && probe->pt <= 30)) continue;
              if (!(probe->pt > 30 && probe->pt <= 40)) continue;
//             if (!(probe->pt > 40 && probe->pt <= 50)) continue;


            if (EtaBin == 0) {
              if (!(fabs(probe->eta) < 1.0)) continue;
            }
            if (EtaBin == 1) {
              if (!(fabs(probe->eta) >= 1.0 && fabs(probe->eta) < 1.479)) continue;
            }
            if (EtaBin == 2) {
              if (!(fabs(probe->eta) >= 1.479 && fabs(probe->eta) < 2.0)) continue;
            }
            if (EtaBin == 3) {
              if (!(fabs(probe->eta) >= 2.0 && fabs(probe->eta) < 2.25)) continue;
            }
            if (EtaBin == 4) {
              if (!(fabs(probe->eta) >= 2.25 && fabs(probe->eta) < 2.5)) continue;
            }
            

            TLorentzVector vprobe;
            vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
            
            TLorentzVector vdimuon = vtag + vprobe;
            if((vdimuon.M()<75) || (vdimuon.M()>110)) continue;
            
            //probe is ok
            Mu_ChargedIso03_VS_NVtx->Fill(probe->ChargedIso03,info->nPV0);
            Mu_NeutralIso03_VS_NVtx->Fill(probe->NeutralIso03_01Threshold,info->nPV0); 
            Mu_ChargedIso04_VS_NVtx->Fill(probe->ChargedIso04,info->nPV0); 
            Mu_NeutralIso04_VS_NVtx->Fill(probe->NeutralIso04_01Threshold,info->nPV0); 
            Mu_HadEnergy_VS_NVtx->Fill(probe->HadEnergy,info->nPV0); 
            Mu_EmEnergy_VS_NVtx->Fill(probe->EmEnergy,info->nPV0); 
            Mu_HoEnergy_VS_NVtx->Fill(probe->HoEnergy,info->nPV0); 
            Mu_HadS9Energy_VS_NVtx->Fill(probe->HadS9Energy,info->nPV0); 
            Mu_EmS9Energy_VS_NVtx->Fill(probe->EmS9Energy,info->nPV0); 
            Mu_HoS9Energy_VS_NVtx->Fill(probe->HoS9Energy,info->nPV0); 
            Mu_TrkIso03_VS_NVtx->Fill(probe->trkIso03,info->nPV0); 
            Mu_EMIso03_VS_NVtx->Fill(probe->emIso03,info->nPV0); 
            Mu_HadIso03_VS_NVtx->Fill(probe->hadIso03,info->nPV0); 
            Mu_TrkIso05_VS_NVtx->Fill(probe->trkIso05,info->nPV0); 
            Mu_EMIso05_VS_NVtx->Fill(probe->emIso05,info->nPV0); 
            Mu_HadIso05_VS_NVtx->Fill(probe->hadIso05,info->nPV0); 

            if (!FilledRho) {
              Mu_Rho_VS_NVtx->Fill(info->PileupEnergyDensity,info->nPV0);
              FilledRho = kTRUE;
            }
            
          }//loop over probes
        } //loop over tags        
      } // if data

    } //loop over events

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  gBenchmark->Show("WWTemplate");       



  TGraphAsymmErrors *Mu_ChargedIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_NeutralIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_ChargedIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_NeutralIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HadEnergy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HoEnergy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_EmEnergy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HadS9Energy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HoS9Energy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_EmS9Energy_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_TrkIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_EMIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HadIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_TrkIso05_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_EMIso05_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_HadIso05_Vs_NPU_Graph;
  TGraphAsymmErrors *Mu_Rho_Vs_NPU_Graph;
  if(isMC) {
    Mu_ChargedIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_ChargedIso03_VS_NPU,("Mu_ChargedIso03_VS_NPU"+label+"_Graph").c_str(),"ChargedIso03 [GeV]");
    Mu_NeutralIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_NeutralIso03_VS_NPU,("Mu_NeutralIso03_VS_NPU"+label+"_Graph").c_str(),"NeutralIso03 [GeV]");
    Mu_ChargedIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_ChargedIso04_VS_NPU,("Mu_ChargedIso04_VS_NPU"+label+"_Graph").c_str(),"ChargedIso04 [GeV]");
    Mu_NeutralIso04_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_NeutralIso04_VS_NPU,("Mu_NeutralIso04_VS_NPU"+label+"_Graph").c_str(),"NeutralIso04 [GeV]");
    Mu_HadEnergy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HadEnergy_VS_NPU,("Mu_HadEnergy_VS_NPU"+label+"_Graph").c_str(),"HadEnergy [GeV]");
    Mu_HoEnergy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HoEnergy_VS_NPU,("Mu_HoEnergy_VS_NPU"+label+"_Graph").c_str(),"HoEnergy [GeV]");
    Mu_EmEnergy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_EmEnergy_VS_NPU,("Mu_EmEnergy_VS_NPU"+label+"_Graph").c_str(),"EmEnergy [GeV]");
    Mu_HadS9Energy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HadS9Energy_VS_NPU,("Mu_HadS9Energy_VS_NPU"+label+"_Graph").c_str(),"HadS9Energy [GeV]");
    Mu_HoS9Energy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HoS9Energy_VS_NPU,("Mu_HoS9Energy_VS_NPU"+label+"_Graph").c_str(),"HoS9Energy [GeV]");
    Mu_EmS9Energy_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_EmS9Energy_VS_NPU,("Mu_EmS9Energy_VS_NPU"+label+"_Graph").c_str(),"EmS9Energy [GeV]");
    Mu_TrkIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_TrkIso03_VS_NPU,("Mu_TrkIso03_VS_NPU"+label+"_Graph").c_str(),"TrkIso03 [GeV]");
    Mu_EMIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_EMIso03_VS_NPU,("Mu_EMIso03_VS_NPU"+label+"_Graph").c_str(),"EMIso03 [GeV]");
    Mu_HadIso03_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HadIso03_VS_NPU,("Mu_HadIso03_VS_NPU"+label+"_Graph").c_str(),"HadIso03 [GeV]");
    Mu_TrkIso05_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_TrkIso05_VS_NPU,("Mu_TrkIso05_VS_NPU"+label+"_Graph").c_str(),"TrkIso05 [GeV]");
    Mu_EMIso05_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_EMIso05_VS_NPU,("Mu_EMIso05_VS_NPU"+label+"_Graph").c_str(),"EMIso05 [GeV]");
    Mu_HadIso05_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_HadIso05_VS_NPU,("Mu_HadIso05_VS_NPU"+label+"_Graph").c_str(),"HadIso05 [GeV]");
    Mu_Rho_Vs_NPU_Graph = MakeMeanVsNPUGraph(Mu_Rho_VS_NPU,("Mu_Rho_VS_NPU"+label+"_Graph").c_str(),"#rho [GeV]");
  }

  TGraphAsymmErrors *Mu_ChargedIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_ChargedIso03_VS_NVtx,("Mu_ChargedIso03_VS_NVtx"+label+"_Graph").c_str(),"ChargedIso03 [GeV]");
  TGraphAsymmErrors *Mu_NeutralIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_NeutralIso03_VS_NVtx,("Mu_NeutralIso03_VS_NVtx"+label+"_Graph").c_str(),"NeutralIso03 [GeV]");
  TGraphAsymmErrors *Mu_ChargedIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_ChargedIso04_VS_NVtx,("Mu_ChargedIso04_VS_NVtx"+label+"_Graph").c_str(),"ChargedIso04 [GeV]");
  TGraphAsymmErrors *Mu_NeutralIso04_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_NeutralIso04_VS_NVtx,("Mu_NeutralIso04_VS_NVtx"+label+"_Graph").c_str(),"NeutralIso04 [GeV]");
  TGraphAsymmErrors *Mu_HadEnergy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HadEnergy_VS_NVtx,("Mu_HadEnergy_VS_NVtx"+label+"_Graph").c_str(),"HadEnergy [GeV]");
  TGraphAsymmErrors *Mu_HoEnergy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HoEnergy_VS_NVtx,("Mu_HoEnergy_VS_NVtx"+label+"_Graph").c_str(),"HoEnergy [GeV]");
  TGraphAsymmErrors *Mu_EmEnergy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_EmEnergy_VS_NVtx,("Mu_EmEnergy_VS_NVtx"+label+"_Graph").c_str(),"EmEnergy [GeV]");
  TGraphAsymmErrors *Mu_HadS9Energy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HadS9Energy_VS_NVtx,("Mu_HadS9Energy_VS_NVtx"+label+"_Graph").c_str(),"HadS9Energy [GeV]");
  TGraphAsymmErrors *Mu_HoS9Energy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HoS9Energy_VS_NVtx,("Mu_HoS9Energy_VS_NVtx"+label+"_Graph").c_str(),"HoS9Energy [GeV]");
  TGraphAsymmErrors *Mu_EmS9Energy_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_EmS9Energy_VS_NVtx,("Mu_EmS9Energy_VS_NVtx"+label+"_Graph").c_str(),"EmS9Energy [GeV]");
  TGraphAsymmErrors *Mu_TrkIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_TrkIso03_VS_NVtx,("Mu_TrkIso03_VS_NVtx"+label+"_Graph").c_str(),"TrkIso03 [GeV]");  
  TGraphAsymmErrors *Mu_EMIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_EMIso03_VS_NVtx,("Mu_EMIso03_VS_NVtx"+label+"_Graph").c_str(),"EMIso03 [GeV]");  
  TGraphAsymmErrors *Mu_HadIso03_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HadIso03_VS_NVtx,("Mu_HadIso03_VS_NVtx"+label+"_Graph").c_str(),"HadIso03 [GeV]");  
  TGraphAsymmErrors *Mu_TrkIso05_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_TrkIso05_VS_NVtx,("Mu_TrkIso05_VS_NVtx"+label+"_Graph").c_str(),"TrkIso05 [GeV]");  
  TGraphAsymmErrors *Mu_EMIso05_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_EMIso05_VS_NVtx,("Mu_EMIso05_VS_NVtx"+label+"_Graph").c_str(),"EMIso05 [GeV]");  
  TGraphAsymmErrors *Mu_HadIso05_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_HadIso05_VS_NVtx,("Mu_HadIso05_VS_NVtx"+label+"_Graph").c_str(),"HadIso05 [GeV]");  
  TGraphAsymmErrors *Mu_Rho_Vs_NVtx_Graph = MakeMeanVsNPUGraph(Mu_Rho_VS_NVtx,("Mu_Rho_VS_NVtx"+label+"_Graph").c_str(),"#rho [GeV]");
 


  //*****************************************************************************************
  //Save Histograms in file
  //*****************************************************************************************
  TFile *file = new TFile("EffectiveArea.root", "UPDATE");
  file->cd();
  if(isMC) {
    file->WriteTObject(Mu_ChargedIso03_Vs_NPU_Graph , Mu_ChargedIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_NeutralIso03_Vs_NPU_Graph , Mu_NeutralIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_ChargedIso04_Vs_NPU_Graph , Mu_ChargedIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_NeutralIso04_Vs_NPU_Graph , Mu_NeutralIso04_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HadEnergy_Vs_NPU_Graph , Mu_HadEnergy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HoEnergy_Vs_NPU_Graph , Mu_HoEnergy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_EmEnergy_Vs_NPU_Graph , Mu_EmEnergy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HadS9Energy_Vs_NPU_Graph , Mu_HadS9Energy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HoS9Energy_Vs_NPU_Graph , Mu_HoS9Energy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_EmS9Energy_Vs_NPU_Graph , Mu_EmS9Energy_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_TrkIso03_Vs_NPU_Graph , Mu_TrkIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_EMIso03_Vs_NPU_Graph , Mu_EMIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HadIso03_Vs_NPU_Graph , Mu_HadIso03_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_TrkIso05_Vs_NPU_Graph , Mu_TrkIso05_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_EMIso05_Vs_NPU_Graph , Mu_EMIso05_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_HadIso05_Vs_NPU_Graph , Mu_HadIso05_Vs_NPU_Graph->GetName(), "WriteDelete");  
    file->WriteTObject(Mu_Rho_Vs_NPU_Graph , Mu_Rho_Vs_NPU_Graph->GetName(), "WriteDelete");  
  }

  file->WriteTObject(Mu_ChargedIso03_Vs_NVtx_Graph , Mu_ChargedIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_NeutralIso03_Vs_NVtx_Graph , Mu_NeutralIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_ChargedIso04_Vs_NVtx_Graph , Mu_ChargedIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_NeutralIso04_Vs_NVtx_Graph , Mu_NeutralIso04_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HadEnergy_Vs_NVtx_Graph , Mu_HadEnergy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HoEnergy_Vs_NVtx_Graph , Mu_HoEnergy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_EmEnergy_Vs_NVtx_Graph , Mu_EmEnergy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HadS9Energy_Vs_NVtx_Graph , Mu_HadS9Energy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HoS9Energy_Vs_NVtx_Graph , Mu_HoS9Energy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_EmS9Energy_Vs_NVtx_Graph , Mu_EmS9Energy_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_TrkIso03_Vs_NVtx_Graph , Mu_TrkIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_EMIso03_Vs_NVtx_Graph , Mu_EMIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HadIso03_Vs_NVtx_Graph , Mu_HadIso03_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_TrkIso05_Vs_NVtx_Graph , Mu_TrkIso05_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_EMIso05_Vs_NVtx_Graph , Mu_EMIso05_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_HadIso05_Vs_NVtx_Graph , Mu_HadIso05_Vs_NVtx_Graph->GetName(), "WriteDelete");  
  file->WriteTObject(Mu_Rho_Vs_NVtx_Graph , Mu_Rho_Vs_NVtx_Graph->GetName(), "WriteDelete");  

  file->WriteTObject(Mu_Rho_VS_NVtx , Mu_Rho_VS_NVtx->GetName(), "WriteDelete");  

  file->Close();
  delete file;

  gBenchmark->Show("WWTemplate");       
} 



void ComputeEffectiveArea(string Label, Bool_t isMC = kTRUE, Int_t EtaBin = -1) {

  string label = "";
  if (Label != "") label = "_" + Label;

  if (EtaBin == 0) {
    label = label + "_Eta0To1";
  }
  if (EtaBin == 1) {
    label = label + "_Eta1To1p479";
  }
  if (EtaBin == 2) {
    label = label + "_Eta1p479To2";
  }
  if (EtaBin == 3) {
    label = label + "_Eta2To2p25";
  }
  if (EtaBin == 4) {
    label = label + "_Eta2p25To2p5";
  }

   TFile *file = new TFile("EffectiveArea.root", "READ");

   TGraphAsymmErrors *Mu_ChargedIso03_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_NeutralIso03_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_ChargedIso04_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_NeutralIso04_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HadEnergy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HoEnergy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_EmEnergy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HadS9Energy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HoS9Energy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_EmS9Energy_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_TrkIso03_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_EMIso03_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HadIso03_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_TrkIso05_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_EMIso05_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_HadIso05_Vs_NPU_Graph = 0;
   TGraphAsymmErrors *Mu_Rho_Vs_NPU_Graph = 0;
   if (isMC) {
     Mu_ChargedIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_ChargedIso03_VS_NPU"+label+"_Graph").c_str());
     Mu_NeutralIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_NeutralIso03_VS_NPU"+label+"_Graph").c_str());
     Mu_ChargedIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_ChargedIso04_VS_NPU"+label+"_Graph").c_str());
     Mu_NeutralIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_NeutralIso04_VS_NPU"+label+"_Graph").c_str());
     Mu_HadEnergy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadEnergy_VS_NPU"+label+"_Graph").c_str());
     Mu_HoEnergy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HoEnergy_VS_NPU"+label+"_Graph").c_str());
     Mu_EmEnergy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EmEnergy_VS_NPU"+label+"_Graph").c_str());
     Mu_HadS9Energy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadS9Energy_VS_NPU"+label+"_Graph").c_str());
     Mu_HoS9Energy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HoS9Energy_VS_NPU"+label+"_Graph").c_str());
     Mu_EmS9Energy_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EmS9Energy_VS_NPU"+label+"_Graph").c_str());
     Mu_TrkIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_TrkIso03_VS_NPU"+label+"_Graph").c_str());
     Mu_EMIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EMIso03_VS_NPU"+label+"_Graph").c_str());
     Mu_HadIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadIso03_VS_NPU"+label+"_Graph").c_str());
     Mu_TrkIso05_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_TrkIso05_VS_NPU"+label+"_Graph").c_str());
     Mu_EMIso05_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EMIso05_VS_NPU"+label+"_Graph").c_str());
     Mu_HadIso05_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadIso05_VS_NPU"+label+"_Graph").c_str());
     Mu_Rho_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Mu_Rho_VS_NPU"+label+"_Graph").c_str());
   }

   TGraphAsymmErrors *Mu_ChargedIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_ChargedIso03_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_NeutralIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_NeutralIso03_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_ChargedIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_ChargedIso04_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_NeutralIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_NeutralIso04_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HadEnergy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadEnergy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HoEnergy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HoEnergy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_EmEnergy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EmEnergy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HadS9Energy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadS9Energy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HoS9Energy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HoS9Energy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_EmS9Energy_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EmS9Energy_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_TrkIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_TrkIso03_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_EMIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EMIso03_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HadIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadIso03_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_TrkIso05_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_TrkIso05_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_EMIso05_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_EMIso05_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_HadIso05_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_HadIso05_VS_NVtx"+label+"_Graph").c_str());
   TGraphAsymmErrors *Mu_Rho_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Mu_Rho_VS_NVtx"+label+"_Graph").c_str());
   
   Double_t Mu_ChargedIso03_Vs_NPU_EffArea = 0;
   Double_t Mu_ChargedIso03_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_NeutralIso03_Vs_NPU_EffArea = 0;
   Double_t Mu_NeutralIso03_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_ChargedIso04_Vs_NPU_EffArea = 0;
   Double_t Mu_ChargedIso04_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_NeutralIso04_Vs_NPU_EffArea = 0;
   Double_t Mu_NeutralIso04_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_HadEnergy_Vs_NPU_EffArea = 0;
   Double_t Mu_HadEnergy_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_HoEnergy_Vs_NPU_EffArea = 0;
   Double_t Mu_HoEnergy_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_EmEnergy_Vs_NPU_EffArea = 0;
   Double_t Mu_EmEnergy_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_HadS9Energy_Vs_NPU_EffArea = 0;
   Double_t Mu_HadS9Energy_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_HoS9Energy_Vs_NPU_EffArea = 0;
   Double_t Mu_HoS9Energy_Vs_NPU_EffAreaErr = 0;
   Double_t Mu_EmS9Energy_Vs_NPU_EffArea = 0;
   Double_t Mu_EmS9Energy_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_TrkIso03_Vs_NPU_EffArea = 0;
   Double_t Mu_TrkIso03_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_EMIso03_Vs_NPU_EffArea = 0;
   Double_t Mu_EMIso03_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_HadIso03_Vs_NPU_EffArea = 0;
   Double_t Mu_HadIso03_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_TrkIso05_Vs_NPU_EffArea = 0;
   Double_t Mu_TrkIso05_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_EMIso05_Vs_NPU_EffArea = 0;
   Double_t Mu_EMIso05_Vs_NPU_EffAreaErr = 0;  
   Double_t Mu_HadIso05_Vs_NPU_EffArea = 0;
   Double_t Mu_HadIso05_Vs_NPU_EffAreaErr = 0;  

   Double_t Mu_Rho_Vs_NPU_Slope = 0;
   Double_t Mu_Rho_Vs_NPU_SlopeErr = 0;
   Double_t Mu_ChargedIso03_Vs_NPU_Slope = 0;
   Double_t Mu_ChargedIso03_Vs_NPU_SlopeErr = 0;
   Double_t Mu_NeutralIso03_Vs_NPU_Slope = 0;
   Double_t Mu_NeutralIso03_Vs_NPU_SlopeErr = 0;
   Double_t Mu_ChargedIso04_Vs_NPU_Slope = 0;
   Double_t Mu_ChargedIso04_Vs_NPU_SlopeErr = 0;
   Double_t Mu_NeutralIso04_Vs_NPU_Slope = 0;
   Double_t Mu_NeutralIso04_Vs_NPU_SlopeErr = 0;
   Double_t Mu_HadEnergy_Vs_NPU_Slope = 0;
   Double_t Mu_HadEnergy_Vs_NPU_SlopeErr = 0;
   Double_t Mu_HoEnergy_Vs_NPU_Slope = 0;
   Double_t Mu_HoEnergy_Vs_NPU_SlopeErr = 0;
   Double_t Mu_EmEnergy_Vs_NPU_Slope = 0;
   Double_t Mu_EmEnergy_Vs_NPU_SlopeErr = 0;
   Double_t Mu_HadS9Energy_Vs_NPU_Slope = 0;
   Double_t Mu_HadS9Energy_Vs_NPU_SlopeErr = 0;
   Double_t Mu_HoS9Energy_Vs_NPU_Slope = 0;
   Double_t Mu_HoS9Energy_Vs_NPU_SlopeErr = 0;
   Double_t Mu_EmS9Energy_Vs_NPU_Slope = 0;
   Double_t Mu_EmS9Energy_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_TrkIso03_Vs_NPU_Slope = 0;
   Double_t Mu_TrkIso03_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_EMIso03_Vs_NPU_Slope = 0;
   Double_t Mu_EMIso03_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_HadIso03_Vs_NPU_Slope = 0;
   Double_t Mu_HadIso03_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_TrkIso05_Vs_NPU_Slope = 0;
   Double_t Mu_TrkIso05_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_EMIso05_Vs_NPU_Slope = 0;
   Double_t Mu_EMIso05_Vs_NPU_SlopeErr = 0;  
   Double_t Mu_HadIso05_Vs_NPU_Slope = 0;
   Double_t Mu_HadIso05_Vs_NPU_SlopeErr = 0;  



   TF1 *f1= 0;

   if (isMC) {
     //*****************************************************************************************************
     //*****************************************************************************************************
     //Vs NPU
     //*****************************************************************************************************
     //*****************************************************************************************************

     //----------------------
     // Rho
     //----------------------
     f1 = new TF1(("Mu_Rho_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_Rho_Vs_NPU_Graph->Fit(("Mu_Rho_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_Rho_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_Rho_Vs_NPU_SlopeErr = f1->GetParError(1);

     //----------------------
     // Energy Observables
     //----------------------
     f1 = new TF1(("Mu_ChargedIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_ChargedIso03_Vs_NPU_Graph->Fit(("Mu_ChargedIso03_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_ChargedIso03_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_ChargedIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_ChargedIso03_Vs_NPU_EffArea = Mu_ChargedIso03_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_ChargedIso03_Vs_NPU_EffAreaErr = Mu_ChargedIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_ChargedIso03_Vs_NPU_SlopeErr/Mu_ChargedIso03_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_NeutralIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_NeutralIso03_Vs_NPU_Graph->Fit(("Mu_NeutralIso03_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_NeutralIso03_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_NeutralIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_NeutralIso03_Vs_NPU_EffArea = Mu_NeutralIso03_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_NeutralIso03_Vs_NPU_EffAreaErr = Mu_NeutralIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_NeutralIso03_Vs_NPU_SlopeErr/Mu_NeutralIso03_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_ChargedIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_ChargedIso04_Vs_NPU_Graph->Fit(("Mu_ChargedIso04_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_ChargedIso04_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_ChargedIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_ChargedIso04_Vs_NPU_EffArea = Mu_ChargedIso04_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_ChargedIso04_Vs_NPU_EffAreaErr = Mu_ChargedIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_ChargedIso04_Vs_NPU_SlopeErr/Mu_ChargedIso04_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_NeutralIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_NeutralIso04_Vs_NPU_Graph->Fit(("Mu_NeutralIso04_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_NeutralIso04_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_NeutralIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_NeutralIso04_Vs_NPU_EffArea = Mu_NeutralIso04_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_NeutralIso04_Vs_NPU_EffAreaErr = Mu_NeutralIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_NeutralIso04_Vs_NPU_SlopeErr/Mu_NeutralIso04_Vs_NPU_Slope,2));


     f1 = new TF1(("Mu_HadEnergy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HadEnergy_Vs_NPU_Graph->Fit(("Mu_HadEnergy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HadEnergy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HadEnergy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HadEnergy_Vs_NPU_EffArea = Mu_HadEnergy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HadEnergy_Vs_NPU_EffAreaErr = Mu_HadEnergy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HadEnergy_Vs_NPU_SlopeErr/Mu_HadEnergy_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_HoEnergy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HoEnergy_Vs_NPU_Graph->Fit(("Mu_HoEnergy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HoEnergy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HoEnergy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HoEnergy_Vs_NPU_EffArea = Mu_HoEnergy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HoEnergy_Vs_NPU_EffAreaErr = Mu_HoEnergy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HoEnergy_Vs_NPU_SlopeErr/Mu_HoEnergy_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_EmEnergy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_EmEnergy_Vs_NPU_Graph->Fit(("Mu_EmEnergy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_EmEnergy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_EmEnergy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_EmEnergy_Vs_NPU_EffArea = Mu_EmEnergy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_EmEnergy_Vs_NPU_EffAreaErr = Mu_EmEnergy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_EmEnergy_Vs_NPU_SlopeErr/Mu_EmEnergy_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_HadS9Energy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HadS9Energy_Vs_NPU_Graph->Fit(("Mu_HadS9Energy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HadS9Energy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HadS9Energy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HadS9Energy_Vs_NPU_EffArea = Mu_HadS9Energy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HadS9Energy_Vs_NPU_EffAreaErr = Mu_HadS9Energy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HadS9Energy_Vs_NPU_SlopeErr/Mu_HadS9Energy_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_HoS9Energy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HoS9Energy_Vs_NPU_Graph->Fit(("Mu_HoS9Energy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HoS9Energy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HoS9Energy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HoS9Energy_Vs_NPU_EffArea = Mu_HoS9Energy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HoS9Energy_Vs_NPU_EffAreaErr = Mu_HoS9Energy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HoS9Energy_Vs_NPU_SlopeErr/Mu_HoS9Energy_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_EmS9Energy_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_EmS9Energy_Vs_NPU_Graph->Fit(("Mu_EmS9Energy_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_EmS9Energy_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_EmS9Energy_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_EmS9Energy_Vs_NPU_EffArea = Mu_EmS9Energy_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_EmS9Energy_Vs_NPU_EffAreaErr = Mu_EmS9Energy_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_EmS9Energy_Vs_NPU_SlopeErr/Mu_EmS9Energy_Vs_NPU_Slope,2));


     f1 = new TF1(("Mu_TrkIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_TrkIso03_Vs_NPU_Graph->Fit(("Mu_TrkIso03_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_TrkIso03_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_TrkIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_TrkIso03_Vs_NPU_EffArea = Mu_TrkIso03_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_TrkIso03_Vs_NPU_EffAreaErr = Mu_TrkIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_TrkIso03_Vs_NPU_SlopeErr/Mu_TrkIso03_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_EMIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_EMIso03_Vs_NPU_Graph->Fit(("Mu_EMIso03_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_EMIso03_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_EMIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_EMIso03_Vs_NPU_EffArea = Mu_EMIso03_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_EMIso03_Vs_NPU_EffAreaErr = Mu_EMIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_EMIso03_Vs_NPU_SlopeErr/Mu_EMIso03_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_HadIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HadIso03_Vs_NPU_Graph->Fit(("Mu_HadIso03_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HadIso03_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HadIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HadIso03_Vs_NPU_EffArea = Mu_HadIso03_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HadIso03_Vs_NPU_EffAreaErr = Mu_HadIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HadIso03_Vs_NPU_SlopeErr/Mu_HadIso03_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_TrkIso05_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_TrkIso05_Vs_NPU_Graph->Fit(("Mu_TrkIso05_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_TrkIso05_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_TrkIso05_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_TrkIso05_Vs_NPU_EffArea = Mu_TrkIso05_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_TrkIso05_Vs_NPU_EffAreaErr = Mu_TrkIso05_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_TrkIso05_Vs_NPU_SlopeErr/Mu_TrkIso05_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_EMIso05_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_EMIso05_Vs_NPU_Graph->Fit(("Mu_EMIso05_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_EMIso05_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_EMIso05_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_EMIso05_Vs_NPU_EffArea = Mu_EMIso05_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_EMIso05_Vs_NPU_EffAreaErr = Mu_EMIso05_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_EMIso05_Vs_NPU_SlopeErr/Mu_EMIso05_Vs_NPU_Slope,2));

     f1 = new TF1(("Mu_HadIso05_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
     f1->SetLineColor(kRed);
     Mu_HadIso05_Vs_NPU_Graph->Fit(("Mu_HadIso05_VS_NPU"+label+"_Fit").c_str(),"R");
     Mu_HadIso05_Vs_NPU_Slope = f1->GetParameter(1);
     Mu_HadIso05_Vs_NPU_SlopeErr = f1->GetParError(1);
     Mu_HadIso05_Vs_NPU_EffArea = Mu_HadIso05_Vs_NPU_Slope / Mu_Rho_Vs_NPU_Slope;
     Mu_HadIso05_Vs_NPU_EffAreaErr = Mu_HadIso05_Vs_NPU_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NPU_SlopeErr/Mu_Rho_Vs_NPU_Slope,2) + pow(Mu_HadIso05_Vs_NPU_SlopeErr/Mu_HadIso05_Vs_NPU_Slope,2));


   }

   //*****************************************************************************************************
   //*****************************************************************************************************
   //Vs NVtx
   //*****************************************************************************************************
   //*****************************************************************************************************

   //----------------------
   // Rho
   //----------------------
   f1 = new TF1(("Mu_Rho_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_Rho_Vs_NVtx_Graph->Fit(("Mu_Rho_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_Rho_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_Rho_Vs_NVtx_SlopeErr = f1->GetParError(1);

   //----------------------
   // Energy Observables
   //----------------------
   f1 = new TF1(("Mu_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_ChargedIso03_Vs_NVtx_Graph->Fit(("Mu_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_ChargedIso03_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_ChargedIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_ChargedIso03_Vs_NVtx_EffArea = Mu_ChargedIso03_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_ChargedIso03_Vs_NVtx_EffAreaErr = Mu_ChargedIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_ChargedIso03_Vs_NVtx_SlopeErr/Mu_ChargedIso03_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_NeutralIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_NeutralIso03_Vs_NVtx_Graph->Fit(("Mu_NeutralIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_NeutralIso03_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_NeutralIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_NeutralIso03_Vs_NVtx_EffArea = Mu_NeutralIso03_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_NeutralIso03_Vs_NVtx_EffAreaErr = Mu_NeutralIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_NeutralIso03_Vs_NVtx_SlopeErr/Mu_NeutralIso03_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_ChargedIso04_Vs_NVtx_Graph->Fit(("Mu_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_ChargedIso04_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_ChargedIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_ChargedIso04_Vs_NVtx_EffArea = Mu_ChargedIso04_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_ChargedIso04_Vs_NVtx_EffAreaErr = Mu_ChargedIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_ChargedIso04_Vs_NVtx_SlopeErr/Mu_ChargedIso04_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_NeutralIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_NeutralIso04_Vs_NVtx_Graph->Fit(("Mu_NeutralIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_NeutralIso04_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_NeutralIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_NeutralIso04_Vs_NVtx_EffArea = Mu_NeutralIso04_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_NeutralIso04_Vs_NVtx_EffAreaErr = Mu_NeutralIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_NeutralIso04_Vs_NVtx_SlopeErr/Mu_NeutralIso04_Vs_NVtx_Slope,2));


   f1 = new TF1(("Mu_HadEnergy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HadEnergy_Vs_NVtx_Graph->Fit(("Mu_HadEnergy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HadEnergy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HadEnergy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HadEnergy_Vs_NVtx_EffArea = Mu_HadEnergy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HadEnergy_Vs_NVtx_EffAreaErr = Mu_HadEnergy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HadEnergy_Vs_NVtx_SlopeErr/Mu_HadEnergy_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_HoEnergy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HoEnergy_Vs_NVtx_Graph->Fit(("Mu_HoEnergy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HoEnergy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HoEnergy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HoEnergy_Vs_NVtx_EffArea = Mu_HoEnergy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HoEnergy_Vs_NVtx_EffAreaErr = Mu_HoEnergy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HoEnergy_Vs_NVtx_SlopeErr/Mu_HoEnergy_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_EmEnergy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_EmEnergy_Vs_NVtx_Graph->Fit(("Mu_EmEnergy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_EmEnergy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_EmEnergy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_EmEnergy_Vs_NVtx_EffArea = Mu_EmEnergy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_EmEnergy_Vs_NVtx_EffAreaErr = Mu_EmEnergy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_EmEnergy_Vs_NVtx_SlopeErr/Mu_EmEnergy_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_HadS9Energy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HadS9Energy_Vs_NVtx_Graph->Fit(("Mu_HadS9Energy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HadS9Energy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HadS9Energy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HadS9Energy_Vs_NVtx_EffArea = Mu_HadS9Energy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HadS9Energy_Vs_NVtx_EffAreaErr = Mu_HadS9Energy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HadS9Energy_Vs_NVtx_SlopeErr/Mu_HadS9Energy_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_HoS9Energy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HoS9Energy_Vs_NVtx_Graph->Fit(("Mu_HoS9Energy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HoS9Energy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HoS9Energy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HoS9Energy_Vs_NVtx_EffArea = Mu_HoS9Energy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HoS9Energy_Vs_NVtx_EffAreaErr = Mu_HoS9Energy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HoS9Energy_Vs_NVtx_SlopeErr/Mu_HoS9Energy_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_EmS9Energy_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_EmS9Energy_Vs_NVtx_Graph->Fit(("Mu_EmS9Energy_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_EmS9Energy_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_EmS9Energy_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_EmS9Energy_Vs_NVtx_EffArea = Mu_EmS9Energy_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_EmS9Energy_Vs_NVtx_EffAreaErr = Mu_EmS9Energy_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_EmS9Energy_Vs_NVtx_SlopeErr/Mu_EmS9Energy_Vs_NVtx_Slope,2));


   f1 = new TF1(("Mu_TrkIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_TrkIso03_Vs_NVtx_Graph->Fit(("Mu_TrkIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_TrkIso03_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_TrkIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_TrkIso03_Vs_NVtx_EffArea = Mu_TrkIso03_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_TrkIso03_Vs_NVtx_EffAreaErr = Mu_TrkIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_TrkIso03_Vs_NVtx_SlopeErr/Mu_TrkIso03_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_EMIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_EMIso03_Vs_NVtx_Graph->Fit(("Mu_EMIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_EMIso03_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_EMIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_EMIso03_Vs_NVtx_EffArea = Mu_EMIso03_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_EMIso03_Vs_NVtx_EffAreaErr = Mu_EMIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_EMIso03_Vs_NVtx_SlopeErr/Mu_EMIso03_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_HadIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HadIso03_Vs_NVtx_Graph->Fit(("Mu_HadIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HadIso03_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HadIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HadIso03_Vs_NVtx_EffArea = Mu_HadIso03_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HadIso03_Vs_NVtx_EffAreaErr = Mu_HadIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HadIso03_Vs_NVtx_SlopeErr/Mu_HadIso03_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_TrkIso05_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_TrkIso05_Vs_NVtx_Graph->Fit(("Mu_TrkIso05_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_TrkIso05_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_TrkIso05_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_TrkIso05_Vs_NVtx_EffArea = Mu_TrkIso05_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_TrkIso05_Vs_NVtx_EffAreaErr = Mu_TrkIso05_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_TrkIso05_Vs_NVtx_SlopeErr/Mu_TrkIso05_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_EMIso05_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_EMIso05_Vs_NVtx_Graph->Fit(("Mu_EMIso05_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_EMIso05_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_EMIso05_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_EMIso05_Vs_NVtx_EffArea = Mu_EMIso05_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_EMIso05_Vs_NVtx_EffAreaErr = Mu_EMIso05_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_EMIso05_Vs_NVtx_SlopeErr/Mu_EMIso05_Vs_NVtx_Slope,2));

   f1 = new TF1(("Mu_HadIso05_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
   f1->SetLineColor(kRed);
   Mu_HadIso05_Vs_NVtx_Graph->Fit(("Mu_HadIso05_VS_NVtx"+label+"_Fit").c_str(),"R");
   Double_t Mu_HadIso05_Vs_NVtx_Slope = f1->GetParameter(1);
   Double_t Mu_HadIso05_Vs_NVtx_SlopeErr = f1->GetParError(1);
   Double_t Mu_HadIso05_Vs_NVtx_EffArea = Mu_HadIso05_Vs_NVtx_Slope / Mu_Rho_Vs_NVtx_Slope;
   Double_t Mu_HadIso05_Vs_NVtx_EffAreaErr = Mu_HadIso05_Vs_NVtx_EffArea * TMath::Sqrt(pow( Mu_Rho_Vs_NVtx_SlopeErr/Mu_Rho_Vs_NVtx_Slope,2) + pow(Mu_HadIso05_Vs_NVtx_SlopeErr/Mu_HadIso05_Vs_NVtx_Slope,2));

   cout << endl << endl 
        << "***********************************************************************************"
        << endl << endl ;

   char buffer[200];
   char buffer2[200];

   if (isMC) {
     sprintf(buffer,"%.3f +/- %.3f",Mu_ChargedIso03_Vs_NPU_EffArea,Mu_ChargedIso03_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_ChargedIso03_Vs_NVtx_EffArea,Mu_ChargedIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_ChargedIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_NeutralIso03_Vs_NPU_EffArea,Mu_NeutralIso03_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_NeutralIso03_Vs_NVtx_EffArea,Mu_NeutralIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_NeutralIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_ChargedIso04_Vs_NPU_EffArea,Mu_ChargedIso04_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_ChargedIso04_Vs_NVtx_EffArea,Mu_ChargedIso04_Vs_NVtx_EffAreaErr);
     cout << "Mu_ChargedIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_NeutralIso04_Vs_NPU_EffArea,Mu_NeutralIso04_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_NeutralIso04_Vs_NVtx_EffArea,Mu_NeutralIso04_Vs_NVtx_EffAreaErr);
     cout << "Mu_NeutralIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HadEnergy_Vs_NPU_EffArea,Mu_HadEnergy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadEnergy_Vs_NVtx_EffArea,Mu_HadEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadEnergy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HoEnergy_Vs_NPU_EffArea,Mu_HoEnergy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HoEnergy_Vs_NVtx_EffArea,Mu_HoEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HoEnergy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_EmEnergy_Vs_NPU_EffArea,Mu_EmEnergy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EmEnergy_Vs_NVtx_EffArea,Mu_EmEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_EmEnergy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HadS9Energy_Vs_NPU_EffArea,Mu_HadS9Energy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadS9Energy_Vs_NVtx_EffArea,Mu_HadS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadS9Energy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HoS9Energy_Vs_NPU_EffArea,Mu_HoS9Energy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HoS9Energy_Vs_NVtx_EffArea,Mu_HoS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HoS9Energy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_EmS9Energy_Vs_NPU_EffArea,Mu_EmS9Energy_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EmS9Energy_Vs_NVtx_EffArea,Mu_EmS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_EmS9Energy_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_TrkIso03_Vs_NPU_EffArea,Mu_TrkIso03_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_TrkIso03_Vs_NVtx_EffArea,Mu_TrkIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_TrkIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_EMIso03_Vs_NPU_EffArea,Mu_EMIso03_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EMIso03_Vs_NVtx_EffArea,Mu_EMIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_EMIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HadIso03_Vs_NPU_EffArea,Mu_HadIso03_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadIso03_Vs_NVtx_EffArea,Mu_HadIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_TrkIso05_Vs_NPU_EffArea,Mu_TrkIso05_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_TrkIso05_Vs_NVtx_EffArea,Mu_TrkIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_TrkIso05_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_EMIso05_Vs_NPU_EffArea,Mu_EMIso05_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EMIso05_Vs_NVtx_EffArea,Mu_EMIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_EMIso05_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

     sprintf(buffer,"%.3f +/- %.3f",Mu_HadIso05_Vs_NPU_EffArea,Mu_HadIso05_Vs_NPU_EffAreaErr);
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadIso05_Vs_NVtx_EffArea,Mu_HadIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadIso05_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;


     //*************************************************
     //Make Some Plots
     //*************************************************
     TCanvas *cv = 0;
     TPaveLabel *label = 0;

     cv = new TCanvas("cv","cv",800,600);
     Mu_NeutralIso04_Vs_NPU_Graph->Draw("AP");
     Mu_NeutralIso04_Vs_NPU_Graph->GetXaxis()->SetRangeUser(0,35);
     label = new TPaveLabel(5.0,5.5,20.0,6.0, Form("Slope = %.3f +/- %.3f", Mu_NeutralIso04_Vs_NPU_Slope, Mu_NeutralIso04_Vs_NPU_SlopeErr));
     label->SetFillStyle(0);
     label->SetBorderSize(0);
     label->SetTextSize(0.6);
     label->Draw();
     cv->SaveAs("Mu_NeutralIso04_Vs_NPU_Graph.gif");

     cv = new TCanvas("cv","cv",800,600);
     Mu_Rho_Vs_NPU_Graph->Draw("AP");
     Mu_Rho_Vs_NPU_Graph->GetXaxis()->SetRangeUser(0,35);
     label = new TPaveLabel(5.0,20,20.0,22.5, Form("Slope = %.3f +/- %.3f", Mu_Rho_Vs_NPU_Slope, Mu_Rho_Vs_NPU_SlopeErr));
     label->SetFillStyle(0);
     label->SetBorderSize(0);
     label->SetTextSize(0.6);
     label->Draw();
     cv->SaveAs("Mu_Rho_Vs_NPU_Graph.gif");
     
     

   } else {

     sprintf(buffer2,"%.3f +/- %.3f",Mu_ChargedIso03_Vs_NVtx_EffArea,Mu_ChargedIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_ChargedIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_NeutralIso03_Vs_NVtx_EffArea,Mu_NeutralIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_NeutralIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_ChargedIso04_Vs_NVtx_EffArea,Mu_ChargedIso04_Vs_NVtx_EffAreaErr);
     cout << "Mu_ChargedIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_NeutralIso04_Vs_NVtx_EffArea,Mu_NeutralIso04_Vs_NVtx_EffAreaErr);
     cout << "Mu_NeutralIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadEnergy_Vs_NVtx_EffArea,Mu_HadEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadEnergy_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HoEnergy_Vs_NVtx_EffArea,Mu_HoEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HoEnergy_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EmEnergy_Vs_NVtx_EffArea,Mu_EmEnergy_Vs_NVtx_EffAreaErr);
     cout << "Mu_EmEnergy_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadS9Energy_Vs_NVtx_EffArea,Mu_HadS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadS9Energy_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_HoS9Energy_Vs_NVtx_EffArea,Mu_HoS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_HoS9Energy_Vs_NVtx_EffArea : " << buffer2 << endl;
     
     sprintf(buffer2,"%.3f +/- %.3f",Mu_EmS9Energy_Vs_NVtx_EffArea,Mu_EmS9Energy_Vs_NVtx_EffAreaErr);
     cout << "Mu_EmS9Energy_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_TrkIso03_Vs_NVtx_EffArea,Mu_TrkIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_TrkIso03_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_EMIso03_Vs_NVtx_EffArea,Mu_EMIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_EMIso03_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadIso03_Vs_NVtx_EffArea,Mu_HadIso03_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadIso03_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_TrkIso05_Vs_NVtx_EffArea,Mu_TrkIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_TrkIso05_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_EMIso05_Vs_NVtx_EffArea,Mu_EMIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_EMIso05_Vs_NVtx_EffArea : " << buffer2 << endl;

     sprintf(buffer2,"%.3f +/- %.3f",Mu_HadIso05_Vs_NVtx_EffArea,Mu_HadIso05_Vs_NVtx_EffAreaErr);
     cout << "Mu_HadIso05_Vs_NVtx_EffArea : " << buffer2 << endl;

   }

}

void ComputeMuonIsolationEffectiveArea() {
  
//   FillIsolationHistograms("/data/blue/sixie/ntuples/HWW/mc/HwwAnalysis_s11-h115ww2l-gf-v11-pu_noskim_normalized.root", "HWW115");


  //Fall11 Z->mumu
//   FillIsolationHistograms("Fall11ZmmMC", "Fall11ZmmMC",  0);
//   FillIsolationHistograms("Fall11ZmmMC", "Fall11ZmmMC", 1);
//   FillIsolationHistograms("Fall11ZmmMC", "Fall11ZmmMC", 2);
//   FillIsolationHistograms("Fall11ZmmMC", "Fall11ZmmMC", 3);
//   FillIsolationHistograms("Fall11ZmmMC", "Fall11ZmmMC", 4);
//   ComputeEffectiveArea("Fall11ZmmMC",kTRUE,0);
//   ComputeEffectiveArea("Fall11ZmmMC",kTRUE,1);
//   ComputeEffectiveArea("Fall11ZmmMC",kTRUE,2);
//   ComputeEffectiveArea("Fall11ZmmMC",kTRUE,3);
//   ComputeEffectiveArea("Fall11ZmmMC",kTRUE,4);

  //Data
  FillIsolationHistograms("Data2011", "Data2011", 0);
  FillIsolationHistograms("Data2011", "Data2011", 1);
  FillIsolationHistograms("Data2011", "Data2011", 2);
  FillIsolationHistograms("Data2011", "Data2011", 3);
  FillIsolationHistograms("Data2011", "Data2011", 4);
   ComputeEffectiveArea("Data2011",kFALSE, 0);
   ComputeEffectiveArea("Data2011",kFALSE, 1);
   ComputeEffectiveArea("Data2011",kFALSE, 2);
   ComputeEffectiveArea("Data2011",kFALSE, 3);
   ComputeEffectiveArea("Data2011",kFALSE, 4);

//   ComputeEffectiveArea("Zmm",kTRUE);

}
