//root -l EWKAna/Hww/Selection/HwwSelection.C+\(\"/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_full-d22_triggerv4skim.root\",130,\"\"\)



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

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passElectronCuts(const mithep::TElectron *ele);
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


//=== MAIN MACRO =================================================================================================

void HwwSelection(const string inputFilename, Double_t mHiggs,
                 const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi = 35.5;              // luminosity (pb^-1)
  Int_t ChargeSelection = 0;
//   ChargeSelection = 1;
    
  Double_t fPtMaxLowerCut;
  Double_t fPtMinLowerCut;
  Double_t fDileptonMassUpperCut;
  Double_t fDeltaPhiCut;

  Double_t fHiggsMass = mHiggs;

  //Configure Higgs Mass Dependant Cuts
  if (fHiggsMass == 120) { fPtMaxLowerCut = 20.0;        fPtMinLowerCut = 15.0; 
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


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
 TH1D *fHWWSelection= new TH1D("hHWWSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEESelection= new TH1D("hHWWToEESelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToMuMuSelection= new TH1D("hHWWToMuMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);
  TH1D *fHWWToEMuSelection= new TH1D("hHWWToEMuSelection", ";Cut Number;Number of Events", 15, -1.5, 13.5);

  TH1D *fLeptonEta = new TH1D(          "hLeptonEta",";LeptonEta;Number of Events",100,-5.,5.0);
  TH1D *fLeptonPtMax = new TH1D(        "hLeptonPtMax",";Lepton P_t Max;Number of Events",150,0.,150.);
  TH1D *fLeptonPtMin = new TH1D(        "hLeptonPtMin",";Lepton P_t Min;Number of Events",150,0.,150.);
  TH1D *fMetPtHist = new TH1D(          "hMetPtHist",";Met;Number of Events",150,0.,300.);  
  TH1D *fMetPhiHist = new TH1D(         "hMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  TH1D *fUncorrMetPtHist = new TH1D(    "hUncorrMetPtHist",";Met;Number of Events",150,0.,300.);  
  TH1D *fUncorrMetPhiHist = new TH1D(   "hUncorrMetPhiHist",";#phi;Number of Events",28,-3.5,3.5);
  TH1D *fDeltaPhiLeptons = new TH1D(    "hDeltaPhiLeptons",";#Delta#phi_{ll};Number of Events",90,0,180);
  TH1D *fDeltaEtaLeptons = new TH1D(    "hDeltaEtaLeptons",";#Delta#eta_{ll};Number of Events",100,-50.,5.0);

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

    //event weight
    Double_t eventweight = info->eventweight * lumi;

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

    
        //*********************************************************************************************
        //Define Cuts
        //*********************************************************************************************
        const int nCuts = 14;
        bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false,
                               false, false, false, false};
  
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
    



        //Cut Selection Histograms
        fHWWSelection->Fill(-1, eventweight);
        if (finalState == 1 ) {
          fHWWToMuMuSelection->Fill(-1, eventweight);
        } else if(finalState == 0 ) {
          fHWWToEESelection->Fill(-1, eventweight);
        } else if(finalState == 2 || finalState == 3 ) {
          fHWWToEMuSelection->Fill(-1, eventweight);
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
            fHWWSelection->Fill(k, eventweight);
            if (finalState == 1 ) {
              fHWWToMuMuSelection->Fill(k, eventweight);
            } else if(finalState == 0 ) {
              fHWWToEESelection->Fill(k, eventweight);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelection->Fill(k, eventweight);
            }
          }
        }

        //*****************************************************************************************
        //Make Preselection Histograms  
        //*****************************************************************************************
        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5]) {
          fLeptonEta->Fill(lepton1.Eta(), eventweight); 
          fLeptonEta->Fill(lepton2.Eta(), eventweight);
          fLeptonPtMax->Fill(lepton1.Pt(), eventweight);
          fLeptonPtMin->Fill(lepton2.Pt(), eventweight);
          fMetPtHist->Fill(met.Pt(), eventweight);                             
          fMetPhiHist->Fill(met.Phi(), eventweight);                            
          fDeltaPhiLeptons->Fill(deltaPhiLeptons, eventweight);
        
        }

    //     if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7] && passCut[8] && passCut[9]) {
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);
//         }



        //*********************************************************************************************
        //Plots after all Cuts
        //*********************************************************************************************
        if (passAllCuts) {
          fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton, eventweight);
          fMtLepton1_afterCuts->Fill(mTW[0], eventweight);
          fMtLepton2_afterCuts->Fill(mTW[1], eventweight);
          fMtHiggs_afterCuts->Fill(mtHiggs, eventweight);
          fLeptonPtPlusMet_afterCuts->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt(), eventweight);
    
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);

        }
      }
    }
  

  } //end loop over data     



  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  for (int i=1; i < fHWWToEESelection->GetXaxis()->GetNbins()+1; ++i) {
    cout << fHWWToEESelection->GetBinContent(i) << " " << fHWWToMuMuSelection->GetBinContent(i) << " " << fHWWToEMuSelection->GetBinContent(i) << " " << fHWWSelection->GetBinContent(i) << endl;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlots.root", "RECREATE");
  
  file->WriteTObject(fHWWSelection, fHWWSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEESelection, fHWWToEESelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToMuMuSelection , fHWWToMuMuSelection->GetName(), "WriteDelete");
  file->WriteTObject(fHWWToEMuSelection,fHWWToEMuSelection->GetName(), "WriteDelete");
  
  file->Close();
  delete file;

        
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t isMC = kFALSE;
  Bool_t pass = kTRUE;

  return pass;

}



Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  Double_t isoCut = 0.1;
  if (ele->pt > 20) {
    if (fabs(ele->eta) > 1.5) isoCut = 0.1;
  } else {
    if (fabs(ele->eta) > 1.5) {
      isoCut = 0.07;
    } else {
      isoCut = 0.07;
    }
  }


  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < isoCut
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
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < isoCut
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
