//root -l EWKAna/Hww/LeptonSelection/MakeSignalElectronDistributions.C+\(\)
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

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2);
void runMakeElectronPlots(const string dataInputFilename,       
                          const string Label, Bool_t isMuonData);
void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname);

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

void MakeSignalElectronDistributions() 
{  
  
  runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_f10-h130ww2l-gf-z2-v12-pu_noskim_normalized.root","Hww130_lowPt",kFALSE);


}



void runMakeElectronPlots(const string dataInputFilename,       
                          const string Label, Bool_t isMuonData) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
  
  Double_t fPtMaxLowerCut;
  Double_t fPtMinLowerCut;
  Double_t fDileptonMassUpperCut;
  Double_t fDeltaPhiCut;
  Double_t fHiggsMass = 130;
  
  //Configure Higgs Mass Dependant Cuts
  if (fHiggsMass == 120) { fPtMaxLowerCut = 20.0;        fPtMinLowerCut = 20.0; 
    fDileptonMassUpperCut = 40.0; fDeltaPhiCut = 60.0;
  }
  if (fHiggsMass == 130) { fPtMaxLowerCut = 25.0;        fPtMinLowerCut = 20.0;
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
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_ee = new TH1D("dileptonMass_ee", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_emu = new TH1D("dileptonMass_emu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_mumu = new TH1D("dileptonMass_mumu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  Int_t Count_ee_BB = 0;
  Int_t Count_ee_BE = 0;
  Int_t Count_ee_EE = 0;
  Int_t Count_mm_BB = 0;
  Int_t Count_mm_BE = 0;
  Int_t Count_mm_EE = 0;
  Int_t Count_em_BB = 0;
  Int_t Count_em_BE = 0;
  Int_t Count_em_EE = 0;
  ofstream eventListFile("eventList.txt");

  TH1D *ptMin = new TH1D("ptMin", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  TH1D *ptMax = new TH1D("ptMax", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  string tmplabel = Label;
  if (tmplabel != "") tmplabel = "_" + Label;
  //electron variables
  TH1F *sigmaIEtaIEta_Barrel = new TH1F((string("sigmaIEtaIEta_Barrel")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.05);
  TH1F *DeltaEta_Barrel = new TH1F((string("DeltaEta_Barrel")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Barrel = new TH1F((string("DeltaPhi_Barrel")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.03, 0.03);
  TH1F *HOverE_Barrel = new TH1F((string("HOverE_Barrel")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *relIso_Barrel = new TH1F((string("relIso_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Barrel = new TH1F((string("fBrem_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);

  TH1F *sigmaIEtaIEta_Endcap = new TH1F((string("sigmaIEtaIEta_Endcap")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.1);
  TH1F *DeltaEta_Endcap = new TH1F((string("DeltaEta_Endcap")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Endcap = new TH1F((string("DeltaPhi_Endcap")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.03, 0.03);
  TH1F *HOverE_Endcap = new TH1F((string("HOverE_Endcap")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *relIso_Endcap = new TH1F((string("relIso_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Endcap = new TH1F((string("fBrem_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  

  infile = new TFile(dataInputFilename.c_str());
  assert(infile);

    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kFALSE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("merged_JsonReRecoSep17_JsonStreamExpressV2.txt"); 
  

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(), "Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
//     mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//     if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
    
//     if (!passHLT(info->triggerBits, info->runNum, isMuonData)) continue;
    
    
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
           mu->pt > 20.0
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

      if ( (0==0)
           && 
           passElectronCuts(ele)
           &&
           fabs(ele->eta) < 2.5
           && 
           ele->pt > 20.0
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
    //dilepton preselection
    //******************************************************************************
    if (leptonPt.size() < 2) continue;
    if (!(leptonPt[0] > 20.0 && leptonPt[1] > 20.0)) continue;

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
        if (leptonType[0] == 11) {
          lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 0.51099892e-3 );
        } else {
          lepton1.SetCoordinates(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
        }
        if (leptonType[1] == 11) {
          lepton2.SetCoordinates(leptonPt[1], leptonEta[1], leptonPhi[1], 0.51099892e-3 );
        } else {
          lepton2.SetCoordinates(leptonPt[1], leptonEta[1], leptonPhi[1], 105.658369e-3 );
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
        bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                               false, false, false, false};
  
        if(lepton1.Pt() >  20.0 &&
           lepton2.Pt() >= 20.0) passCut[0] = true;

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

        if (!(leptonPt.size() >= 3 && leptonPt[2] > 20.0)) passCut[8] = true;

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
        //Fill Histograms
        //*********************************************************************************************
        if (
//       passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] 
//         && 
          passAllCuts
          &&
      
          (finalState == 0 || finalState == 3 ) 
          ) {
          ptMin->Fill(lepton2.Pt());
          ptMax->Fill(lepton1.Pt());

          //fill electron variables for 2nd electron
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[leptonIndex[1]]);
          if (ele->pt >= 10 && ele->pt < 15 
              && ele->isUsed == kTRUE
            ) {
            if (fabs(ele->eta) < 1.5) {
              sigmaIEtaIEta_Barrel->Fill( ele->sigiEtaiEta );
              DeltaEta_Barrel->Fill( ele->deltaEtaIn );
              DeltaPhi_Barrel->Fill( ele->deltaPhiIn );
              HOverE_Barrel->Fill( ele->HoverE );
              relIso_Barrel->Fill( (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt );      
            } else {          
              sigmaIEtaIEta_Endcap->Fill( ele->sigiEtaiEta );
              DeltaEta_Endcap->Fill( ele->deltaEtaIn );
              DeltaPhi_Endcap->Fill( ele->deltaPhiIn );
              HOverE_Endcap->Fill( ele->HoverE );
              relIso_Endcap->Fill( (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt );
            }
          }
        }
      }
    }
  } //end loop over data     

  infile->Close(); delete infile;
  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlotsFakePrediction.root", "UPDATE");

  file->WriteTObject(ptMin,ptMin->GetName(), "WriteDelete");
  file->WriteTObject(ptMax,ptMax->GetName(), "WriteDelete");


  file->WriteTObject(sigmaIEtaIEta_Barrel,sigmaIEtaIEta_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(DeltaEta_Barrel,DeltaEta_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(DeltaPhi_Barrel,DeltaPhi_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(HOverE_Barrel,HOverE_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(relIso_Barrel,relIso_Barrel->GetName(), "WriteDelete");

  file->WriteTObject(sigmaIEtaIEta_Endcap,sigmaIEtaIEta_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(DeltaEta_Endcap,DeltaEta_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(DeltaPhi_Endcap,DeltaPhi_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(HOverE_Endcap,HOverE_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(relIso_Endcap,relIso_Endcap->GetName(), "WriteDelete");
  file->Close();
//   //--------------------------------------------------------------------------------------------------------------
//   // Make plots
//   //==============================================================================================================
//    TCanvas *cv = new TCanvas("cv","cv", 800,600);
//    ptMin->Draw();
//   cv->SaveAs("ptMin.gif");
//   ptMax->Draw();
//   cv->SaveAs("ptMax.gif");

//   sigmaIEtaIEta_Barrel_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Barrel_signal.gif");
//   DeltaEta_Barrel_signal->Draw();
//   cv->SaveAs("DeltaEta_Barrel_signal.gif");
//   DeltaPhi_Barrel_signal->Draw();
//   cv->SaveAs("DeltaPhi_Barrel_signal.gif");
//   HOverE_Barrel_signal->Draw();
//   cv->SaveAs("HOverE_Barrel_signal.gif");
//   relIso_Barrel_signal->Draw();
//   cv->SaveAs("relIso_Barrel_signal.gif");

//   sigmaIEtaIEta_Endcap_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Endcap_signal.gif");
//   DeltaEta_Endcap_signal->Draw();
//   cv->SaveAs("DeltaEta_Endcap_signal.gif");
//   DeltaPhi_Endcap_signal->Draw();
//   cv->SaveAs("DeltaPhi_Endcap_signal.gif");
//   HOverE_Endcap_signal->Draw();
//   cv->SaveAs("HOverE_Endcap_signal.gif");
//   relIso_Endcap_signal->Draw();
//   cv->SaveAs("relIso_Endcap_signal.gif");


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 


        
  gBenchmark->Show("plotZ");       
} 

void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname) {  
  TFile *file = new TFile("WWSelectionPlotsFakePrediction.root", "READ");

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TH1F *signal = 0;
  TH1F *bkg = 0;
  TLegend *legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
    
 
  signal = (TH1F*)file->Get(signalHistName.c_str());
  bkg = (TH1F*)file->Get(bkgHistName.c_str());
  if (!signal) {cout << "Cannot load histogram " << signalHistName << endl;  assert(signal);}
  if (!bkg) {cout << "Cannot load histogram " << bkgHistName << endl;  assert(bkg);}

  Double_t norm = 0;
  norm = signal->Integral();
  for (int b=0; b<signal->GetXaxis()->GetNbins()+2; ++b) {
    signal->SetBinContent(b,signal->GetBinContent(b) / norm);
    signal->SetBinError(b,signal->GetBinError(b) / norm);
  }
  norm = bkg->Integral();
  for (int b=0; b<bkg->GetXaxis()->GetNbins()+2; ++b) {
    bkg->SetBinContent(b,bkg->GetBinContent(b) / norm);
    bkg->SetBinError(b,bkg->GetBinError(b) / norm);
  }

  Double_t maxY = signal->GetMaximum();
  if (bkg->GetMaximum() > maxY) maxY = bkg->GetMaximum();
  maxY = maxY * 1.2;
  signal->SetMaximum(maxY);
  bkg->SetMaximum(maxY);

  signal->SetLineColor(kRed);
  signal->SetMarkerColor(kRed);
  bkg->SetLineColor(kBlue);
  bkg->SetMarkerColor(kBlue);
  legend->Clear();
  legend->AddEntry(signal, signalLabel.c_str(),  "LP");
  legend->AddEntry(bkg,    bkgLabel.c_str(), "LP");

  signal->Draw();
  bkg->Draw("same");
  legend->Draw();
  cv->SaveAs((plotname+".gif").c_str());

//   sigmaIEtaIEta_Barrel_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Barrel_signal.gif");
//   DeltaEta_Barrel_signal->Draw();
//   cv->SaveAs("DeltaEta_Barrel_signal.gif");
//   DeltaPhi_Barrel_signal->Draw();
//   cv->SaveAs("DeltaPhi_Barrel_signal.gif");
//   HOverE_Barrel_signal->Draw();
//   cv->SaveAs("HOverE_Barrel_signal.gif");
//   relIso_Barrel_signal->Draw();
//   cv->SaveAs("relIso_Barrel_signal.gif");

//   sigmaIEtaIEta_Endcap_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Endcap_signal.gif");
//   DeltaEta_Endcap_signal->Draw();
//   cv->SaveAs("DeltaEta_Endcap_signal.gif");
//   DeltaPhi_Endcap_signal->Draw();
//   cv->SaveAs("DeltaPhi_Endcap_signal.gif");
//   HOverE_Endcap_signal->Draw();
//   cv->SaveAs("HOverE_Endcap_signal.gif");
//   relIso_Endcap_signal->Draw();
//   cv->SaveAs("relIso_Endcap_signal.gif");

  file->Close();
}


Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;

//   if (isMuonData) {
//     if ((runNum >= 132440) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
//     } 
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140401)) {
// //       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   } else {
//     //it's electron data
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140041)) {
// //       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( !(triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   }

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
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / TMath::Max(ele->pt,Float_t(20)) < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
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

        && (mu->nSeg >= 2)
        && (mu->nValidHits >= 1)
        && (mu->nPixHits >= 1)
        )
    ) pass = kFALSE;

  return pass;
}

//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2)
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
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << endl;
  
}
