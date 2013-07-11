//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-eg-s17_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-el-pr-v2_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-mu-pr-v2_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-mu-s17_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_full_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/WWSelection.C+\(\"HwwNtuple_f10-ww2l-z2-v12-pu_noskim_0000.root\",\"\"\)
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
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
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

void WWSelection(const string inputFilename, 
                 const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
  Int_t ChargeSelection = 0;
//   ChargeSelection = 1;
    

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *fHWWSelection= new TH1D("hHWWSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  TH1D *fHWWToEESelection= new TH1D("hHWWToEESelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  TH1D *fHWWToMuMuSelection= new TH1D("hHWWToMuMuSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);
  TH1D *fHWWToEMuSelection= new TH1D("hHWWToEMuSelection", ";Cut Number;Number of Events", 11, -1.5, 9.5);

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
  ofstream eventListFile("eventList.txt");

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
   rlrm.AddJSONFile("Cert_TopNov5_Merged_135821-149442_allPVT.txt"); 
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
        const int nCuts = 10;
        bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false};
  
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

        //*********************************************************************************************
        //Make Selection Histograms. Number of events passing each level of cut
        //*********************************************************************************************  
        bool passAllCuts = true;
        for(int c=0; c<nCuts; c++) passAllCuts = passAllCuts & passCut[c];
    
        if (info->evtNum == 25372 ) cout << "debug : " << passCut[9] << " " << passAllCuts << endl;



        //Cut Selection Histograms
        fHWWSelection->Fill(-1);
        if (finalState == 1 ) {
          fHWWToMuMuSelection->Fill(-1);
        } else if(finalState == 0 ) {
          fHWWToEESelection->Fill(-1);
        } else if(finalState == 2 || finalState == 3 ) {
          fHWWToEMuSelection->Fill(-1);
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
            fHWWSelection->Fill(k);
            if (finalState == 1 ) {
              fHWWToMuMuSelection->Fill(k);
            } else if(finalState == 0 ) {
              fHWWToEESelection->Fill(k);
            } else if(finalState == 2 || finalState == 3 ) {
              fHWWToEMuSelection->Fill(k);
            }
          }
        }

        //*****************************************************************************************
        //Make Preselection Histograms  
        //*****************************************************************************************
        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5]) {
          fLeptonEta->Fill(lepton1.Eta()); 
          fLeptonEta->Fill(lepton2.Eta());
          fLeptonPtMax->Fill(lepton1.Pt());
          fLeptonPtMin->Fill(lepton2.Pt());
          fMetPtHist->Fill(met.Pt());                             
          fMetPhiHist->Fill(met.Phi());                            
          fDeltaPhiLeptons->Fill(deltaPhiLeptons);


          dileptonMass->Fill(dilepton.M());
          if (finalState == 0) {
            dileptonMass_ee->Fill(dilepton.M());
            if (dilepton.M() > 60 && dilepton.M() < 120) {
              if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                Count_ee_BB++;
              } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                Count_ee_EE++;
              } else {
                Count_ee_BE++;
              }
            }
          } else if (finalState == 1) {
            dileptonMass_mumu->Fill(dilepton.M());
            if (dilepton.M() > 60 && dilepton.M() < 120) {
              if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                Count_mm_BB++;
              } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                Count_mm_EE++;
              } else {
                Count_mm_BE++;
              }
            }
          } else {
            dileptonMass_emu->Fill(dilepton.M());
            if (dilepton.M() > 60 && dilepton.M() < 120) {
              if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
                Count_em_BB++;
              } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
                Count_em_EE++;
              } else {
                Count_em_BE++;
              }
            }
          }
        }

    //     if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[6] && passCut[7]) {
//           eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
//                     lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);
//         }



        //*********************************************************************************************
        //Plots after all Cuts
        //*********************************************************************************************
        if (passAllCuts) {
          fMinDeltaPhiLeptonMet_afterCuts->Fill(minDeltaPhiMetLepton);
          fMtLepton1_afterCuts->Fill(mTW[0]);
          fMtLepton2_afterCuts->Fill(mTW[1]);
          fMtHiggs_afterCuts->Fill(mtHiggs);
          fLeptonPtPlusMet_afterCuts->Fill(lepton1.Pt()+lepton2.Pt()+met.Pt());
    
          eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
                    lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), leptonCharge[i], lepton2.Pt(), lepton2.Eta(), lepton2.Phi(), leptonCharge[j]);

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
  dileptonMass->Draw();
  cv->SaveAs("dileptonMass.gif");
  dileptonMass_emu->Draw();
  cv->SaveAs("dileptonMass_emu.gif");
  dileptonMass_ee->Draw();
  cv->SaveAs("dileptonMass_ee.gif");
  dileptonMass_mumu->Draw();
  cv->SaveAs("dileptonMass_mumu.gif");

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  for (int i=1; i < fHWWToEESelection->GetXaxis()->GetNbins()+1; ++i) {
    cout << fHWWToEESelection->GetBinContent(i) << " " << fHWWToMuMuSelection->GetBinContent(i) << " " << fHWWToEMuSelection->GetBinContent(i) << " " << fHWWSelection->GetBinContent(i) << endl;
  }


//   cout << "**************************************************************\n";
//   cout << "Event Count : ee final state\n";
//   cout << "BB :" << Count_ee_BB << endl;
//   cout << "BE :" << Count_ee_BE << endl;
//   cout << "EE :" << Count_ee_EE << endl;
//   cout << "Total :" << Count_ee_BB + Count_ee_BE + Count_ee_EE << endl;
//   cout << "**************************************************************\n";
//   cout << "Event Count : mm final state\n";
//   cout << "BB :" << Count_mm_BB << endl;
//   cout << "BE :" << Count_mm_BE << endl;
//   cout << "EE :" << Count_mm_EE << endl;
//   cout << "Total :" << Count_mm_BB + Count_mm_BE + Count_mm_EE << endl;
//   cout << "**************************************************************\n";
//   cout << "Event Count : em final state\n";
//   cout << "BB :" << Count_em_BB << endl;
//   cout << "BE :" << Count_em_BE << endl;
//   cout << "EE :" << Count_em_EE << endl;
//   cout << "Total :" << Count_em_BB + Count_em_BE + Count_em_EE << endl;









  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("WWSelectionPlots_SS.root", "RECREATE");

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
