//root -l EWKAna/Hww/Acceptance/HwwJetVetoEfficiencyPDFSystematics.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_full_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Acceptance/HwwJetVetoEfficiencyPDFSystematics.C+\(\"HwwNtuple_f10-ww2l-z2-v12-pu_noskim_0000.root\",\"WW\"\)
//root -l EWKAna/Hww/Acceptance/HwwJetVetoEfficiencyPDFSystematics.C+\(\"/home/sixie/hist/HwwAcceptance/HwwNtuple_f10-h160ww2l-gf-z2-v12_noskim_0000.root\",\"/home/sixie/hist/HwwAcceptance/HwwNtuple_f10-zmm-powheg-c10-v12_noskim.root\",\"Hww160\"\)
//root -l EWKAna/Hww/Acceptance/HwwJetVetoEfficiencyPDFSystematics.C+\(\"/home/sixie/hist/HwwAcceptance/HwwNtuple_f10-h400ww2l-gf-z2-v12_noskim_0000.root\",\"/home/sixie/hist/HwwAcceptance/HwwNtuple_f10-zmm-powheg-c10-v12_noskim.root\",\"Hww400\"\)
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
#include "TGraphAsymmErrors.h"

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
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


//=== MAIN MACRO =================================================================================================

void HwwJetVetoEfficiencyPDFSystematics(const string HwwFilename, const string ZllFilename,
                 const string Label) 
{  
  gBenchmark->Start("WWTemplate");
  string label = "";
  if (Label !=  "") label = "_" + Label;

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================



  Double_t lumi;              // luminosity (pb^-1)
  Int_t ChargeSelection = 0;
//   ChargeSelection = 1;
    

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector<TH1D*> HwwLeadingJetPtHists;
  vector<TH1D*> ZllLeadingJetPtHists;
  


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
  mithep::TGenInfo *genInfo    = new mithep::TGenInfo();
  Int_t                   fNPDFMembers;     
  Float_t                 fPDFWeights[100]; 


  //********************************************************
  // Get Tree
  //********************************************************
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;
  TBranch *NPDFMembersBr;
  TBranch *PDFWeightsBr;
  TBranch *genInfoBr;

  //*****************************************************************************************
  //Loop over Hww Data Tree
  //*****************************************************************************************
  eventTree = getTreeFromFile(HwwFilename.c_str(),"Events"); 
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("Gen", &genInfo);          genInfoBr = eventTree->GetBranch("Gen");
  eventTree->SetBranchAddress("NPDFMembers", &fNPDFMembers); NPDFMembersBr = eventTree->GetBranch("NPDFMembers");
  eventTree->SetBranchAddress("PDFWeights", &fPDFWeights);   PDFWeightsBr = eventTree->GetBranch("PDFWeights");
 
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);
		
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
    NPDFMembersBr->GetEntry(ientry);
    PDFWeightsBr->GetEntry(ientry);


    //********************************************************
    // PDF weights
    //********************************************************
//     cout << "NPDFMembers: " << fNPDFMembers << endl;
    for (UInt_t i=0; i<fNPDFMembers; ++i) {
//       cout << "PDF : " << i << " " << fPDFWeights[i] << endl;
      if (ientry == 0) {
        TH1D *tmpHist = new TH1D((string("HwwLeadingJetPtHist") + label + "_PDF" + IntToString(i)).c_str(), "; Leading Jet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
        HwwLeadingJetPtHists.push_back(tmpHist);
      }
    }
    assert(HwwLeadingJetPtHists.size() == fNPDFMembers);
 



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
    double leadingJetPt = 0;
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

      if (jet->pt > 10 && fabs(jet->eta) < 5.0 ) {
        if (jet->pt > leadingJetPt) {
          leadingJetPt = jet->pt;
        }
      }
      

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

    Bool_t passEventSelection = kFALSE;

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
        bool passCut[nCuts] = {false, false, false, false, false, false, false, false, false, false, false, false, false, false};
  
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

        if (lepton1.Pt() > 30)     passCut[10] = true;
        if (lepton2.Pt() > 25)     passCut[11] = true;
        if (dilepton.M() < 50)     passCut[12] = true;
        if (deltaPhiLeptons < 60)  passCut[13] = true;
//         if (lepton1.Pt() > 90)     passCut[10] = true;
//         if (lepton2.Pt() > 25)     passCut[11] = true;
//         if (dilepton.M() < 300)     passCut[12] = true;
//         if (deltaPhiLeptons < 175)  passCut[13] = true;

        //*********************************************************************************************
        //Make Selection Histograms. Number of events passing each level of cut
        //*********************************************************************************************  
         if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[7] && passCut[8] && passCut[10] && passCut[11] && passCut[12] && passCut[13] ) {
 //        if (passCut[0] && passCut[1] && passCut[2] && passCut[3] && passCut[4] && passCut[5] && passCut[7] && passCut[8]  ) {
          passEventSelection = kTRUE;
        }
      } //end loop over leptons
    } //end loop over leptons

    if (passEventSelection) {
      for (UInt_t i=0; i<fNPDFMembers; ++i) {
        if (i != 0) {
          HwwLeadingJetPtHists[i]->Fill(leadingJetPt, fPDFWeights[i]);      
        } else {
//           HwwLeadingJetPtHists[i]->Fill(leadingJetPt, 1.0);      
          HwwLeadingJetPtHists[i]->Fill(leadingJetPt, fPDFWeights[i]);      
        }
      }
    }    
    
  } //end loop over data     

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;





  //*****************************************************************************************
  //Loop over Zee Data Tree
  //*****************************************************************************************
  // Access samples and fill histograms
  TTree *ZllEventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info_Zll    = new mithep::TEventInfo();
  TClonesArray *electronArr_Zll   = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr_Zll       = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr_Zll        = new TClonesArray("mithep::TJet");
  mithep::TGenInfo *genInfo_Zll   = new mithep::TGenInfo();
  Int_t                   fNPDFMembers_Zll;     
  Float_t                 fPDFWeights_Zll[100]; 


  //********************************************************
  // Get Tree
  //********************************************************
  TBranch *infoBr_Zll;
  TBranch *electronBr_Zll;
  TBranch *muonBr_Zll;
  TBranch *jetBr_Zll;
  TBranch *NPDFMembersBr_Zll;
  TBranch *PDFWeightsBr_Zll;
  TBranch *genInfoBr_Zll;

  ZllEventTree = getTreeFromFile(ZllFilename.c_str(),"Events"); 
  // Set branch address to structures that will store the info  
  ZllEventTree->SetBranchAddress("Info",       &info_Zll);      infoBr_Zll       = ZllEventTree->GetBranch("Info");
  ZllEventTree->SetBranchAddress("Electron", &electronArr_Zll); electronBr_Zll = ZllEventTree->GetBranch("Electron");
  ZllEventTree->SetBranchAddress("Muon", &muonArr_Zll);         muonBr_Zll = ZllEventTree->GetBranch("Muon");
  ZllEventTree->SetBranchAddress("PFJet", &jetArr_Zll);         jetBr_Zll = ZllEventTree->GetBranch("PFJet");
  ZllEventTree->SetBranchAddress("Gen", &genInfo_Zll);          genInfoBr_Zll = ZllEventTree->GetBranch("Gen");
  ZllEventTree->SetBranchAddress("NPDFMembers", &fNPDFMembers_Zll); NPDFMembersBr_Zll = ZllEventTree->GetBranch("NPDFMembers");
  ZllEventTree->SetBranchAddress("PDFWeights", &fPDFWeights_Zll);   PDFWeightsBr_Zll = ZllEventTree->GetBranch("PDFWeights");
 
  for(UInt_t ientry=0; ientry<ZllEventTree->GetEntries(); ientry++) {       	
    infoBr_Zll->GetEntry(ientry);
    genInfoBr_Zll->GetEntry(ientry);
		
    //for the skimmed input, I already required the HLT bits.
    //    if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;


    //********************************************************
    // Load the branches
    //********************************************************
    electronArr_Zll->Clear(); 
    muonArr_Zll->Clear(); 
    jetArr_Zll->Clear(); 
    electronBr_Zll->GetEntry(ientry);
    muonBr_Zll->GetEntry(ientry);
    jetBr_Zll->GetEntry(ientry);
    NPDFMembersBr_Zll->GetEntry(ientry);
    PDFWeightsBr_Zll->GetEntry(ientry);


    //********************************************************
    // PDF weights
    //********************************************************
//     cout << "NPDFMembers: " << fNPDFMembers << endl;
    for (UInt_t i=0; i<fNPDFMembers; ++i) {
//       cout << "PDF : " << i << " " << fPDFWeights[i] << endl;
      if (ientry == 0) {
        TH1D *tmpHist = new TH1D((string("ZllLeadingJetPtHist") + label + "_PDF" + IntToString(i)).c_str(), "; Leading Jet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
        ZllLeadingJetPtHists.push_back(tmpHist);
      }
    }
    assert(ZllLeadingJetPtHists.size() == fNPDFMembers);
 



    //********************************************************
    // TcMet
    //********************************************************
    TVector3 met;        
    if(info_Zll->tcMEx!=0 || info_Zll->tcMEy!=0) {       
      met.SetXYZ(info_Zll->tcMEx, info_Zll->tcMEy, 0);
    }
	
     //********************************************************
    // Jets
    //********************************************************
    Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    Double_t leadingJetPt = 0;

    for(Int_t i=0; i<jetArr_Zll->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr_Zll)[i]);

      if (!(jet->pt > 10 && fabs(jet->eta) < 5.0)) continue;
      if (!leadingJet || jet->pt > leadingJet->pt) {
        leadingJet = jet;
        leadingJetPt = jet->pt;
      }
      NJets++;
    }

    //******************************************************************************
    //dilepton preselection
    //******************************************************************************
    if (muonArr_Zll->GetEntries() < 2) continue;

    //******************************************************************************
    //loop over muon pairs
    //******************************************************************************
    Bool_t passZSelection = kFALSE;

   for(Int_t i=0; i<muonArr_Zll->GetEntries(); i++) {
      const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr_Zll)[i]);
      if ( !(
             mu1->pt > 20.0
             &&
             passMuonCuts(mu1)
             )
        ) continue;
      

      for(Int_t j=i+1; j<muonArr_Zll->GetEntries(); j++) {
        const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr_Zll)[j]);
        if ( !(
               mu2->pt > 20.0
               && passMuonCuts(mu2)
               )
          ) continue;
        
        mithep::FourVectorM lepton1;
        mithep::FourVectorM lepton2;
        lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 0.51099892e-3 );
        lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 0.51099892e-3 );
        mithep::FourVectorM dilepton = lepton1+lepton2;

         if (dilepton.M() > 76 
              && dilepton.M() < 106
          ) {
          passZSelection = kTRUE;
        }        
      }
    }

    if (passZSelection) {
      for (UInt_t i=0; i<fNPDFMembers; ++i) {
        if (i!=0) {
          ZllLeadingJetPtHists[i]->Fill(leadingJetPt, fPDFWeights[i]);      
        } else {
//           ZllLeadingJetPtHists[i]->Fill(leadingJetPt, 1.0);      
          ZllLeadingJetPtHists[i]->Fill(leadingJetPt, fPDFWeights[i]);      
        }
      }
    }    
    
  } //end loop over data     

  delete info_Zll;
  delete electronArr_Zll;
  delete muonArr_Zll;
  delete jetArr_Zll;


  //--------------------------------------------------------------------------------------------------------------
  // Compute Jet Veto Efficiency
  //==============================================================================================================
  vector<TH1D*> HwwJetVetoEfficiency;
  vector<TH1D*> ZllJetVetoEfficiency;

  for(int i=0 ; i<HwwLeadingJetPtHists.size() ; ++i) {
    TH1D *tmpJetVetoEfficiency = (TH1D*)HwwLeadingJetPtHists[0]->Clone((string("HwwJetVetoEfficiency") + label + "_PDF" + IntToString(i)).c_str());
    tmpJetVetoEfficiency->GetYaxis()->SetTitle("Jet Veto Efficiency");
    tmpJetVetoEfficiency->GetYaxis()->SetTitleOffset(1.5);
    for (int j=1; j<tmpJetVetoEfficiency->GetXaxis()->GetNbins()+1; ++j) {
      Int_t n = HwwLeadingJetPtHists[i]->Integral(1,j);
      Int_t d = HwwLeadingJetPtHists[i]->Integral(1,HwwLeadingJetPtHists[0]->GetXaxis()->GetNbins()+1);
      Double_t ratio;
      Double_t errLow;
      Double_t errHigh;     
      Int_t errorType = 2;
      mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);
      
      tmpJetVetoEfficiency->SetBinContent(j,ratio);
      tmpJetVetoEfficiency->SetBinError(j,(errLow+errHigh)/2); 
//       cout << "Hww: " << tmpJetVetoEfficiency->GetXaxis()->GetBinLowEdge(j) << " : " << n << " / " << d << " = " 
//            << ratio << " + " << errHigh << " - " << errLow << endl;
    }    
    HwwJetVetoEfficiency.push_back(tmpJetVetoEfficiency);
  }

  for(int i=0 ; i<ZllLeadingJetPtHists.size() ; ++i) {
    TH1D *tmpJetVetoEfficiency = (TH1D*)ZllLeadingJetPtHists[0]->Clone((string("ZllJetVetoEfficiency") + label + "_PDF" + IntToString(i)).c_str());
    tmpJetVetoEfficiency->GetYaxis()->SetTitle("Jet Veto Efficiency");
    tmpJetVetoEfficiency->GetYaxis()->SetTitleOffset(1.5);
    for (int j=1; j<tmpJetVetoEfficiency->GetXaxis()->GetNbins()+1; ++j) {
      Int_t n = ZllLeadingJetPtHists[i]->Integral(1,j);
      Int_t d = ZllLeadingJetPtHists[i]->Integral(1,ZllLeadingJetPtHists[0]->GetXaxis()->GetNbins()+1);
      Double_t ratio;
      Double_t errLow;
      Double_t errHigh;     
      Int_t errorType = 2;
      mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);
      
      tmpJetVetoEfficiency->SetBinContent(j,ratio);
      tmpJetVetoEfficiency->SetBinError(j,(errLow+errHigh)/2); 
//       cout << "Zll : " << tmpJetVetoEfficiency->GetXaxis()->GetBinLowEdge(j) << " : " << n << " / " << d << " = " 
//            << ratio << " + " << errHigh << " - " << errLow << endl;
    }    
    ZllJetVetoEfficiency.push_back(tmpJetVetoEfficiency);    
  }

  //--------------------------------------------------------------------------------------------------------------
  // Compute Jet Veto Efficiency Scale Factor
  //==============================================================================================================
  Double_t x[201];
  Double_t y[201];
  Double_t xErr[201];
  Double_t yErrLow[201];
  Double_t yErrHigh[201];
  for (int j=0; j<201; ++j) {
    x[j] = 0;
    y[j] = 0;
    xErr[j] = 0;
    yErrLow[j] = 0;
    yErrHigh[j] = 0;
  }


  TH1D *jetVetoEfficiencySFHist = (TH1D*)HwwJetVetoEfficiency[0]->Clone("jetVetoEfficiencySFHist");

//   for (int j=1; j<HwwJetVetoEfficiency[0]->GetXaxis()->GetNbins(); ++j) {
   for (int j=1; j<HwwJetVetoEfficiency[0]->GetXaxis()->GetNbins()+1; ++j) {
    Int_t NComponents = (HwwJetVetoEfficiency.size() - 1) / 2;
    Double_t TotalSystematicHighSqr = 0;
    Double_t TotalSystematicLowSqr = 0;
    Double_t SFDefault = HwwJetVetoEfficiency[0]->GetBinContent(j) / ZllJetVetoEfficiency[0]->GetBinContent(j);
    for(int i=0 ; i<NComponents ; ++i) {
      Double_t SFUp = HwwJetVetoEfficiency[2*i + 1]->GetBinContent(j) / ZllJetVetoEfficiency[2*i + 1]->GetBinContent(j);
      Double_t SFDown = HwwJetVetoEfficiency[2*i + 2]->GetBinContent(j) / ZllJetVetoEfficiency[2*i + 2]->GetBinContent(j);
      Double_t systematicsHigh = (SFUp - SFDefault) / SFDefault;
      Double_t systematicsLow = (SFDown - SFDefault) / SFDefault;
      if (systematicsHigh > systematicsLow) {
        if (systematicsHigh < 0) systematicsHigh = 0;
        if (systematicsLow > 0) systematicsLow = 0;
        TotalSystematicHighSqr += pow(systematicsHigh,2);
        TotalSystematicLowSqr += pow(systematicsLow,2);
      } else {
        if (systematicsHigh > 0) systematicsHigh = 0;
        if (systematicsLow < 0) systematicsLow = 0;
        TotalSystematicHighSqr += pow(systematicsLow,2);
        TotalSystematicLowSqr += pow(systematicsHigh,2);
      }
    }
    x[j] = HwwJetVetoEfficiency[0]->GetXaxis()->GetBinLowEdge(j) ;
    xErr[j] = 0;
    y[j] = SFDefault;
    yErrLow[j] = TMath::Sqrt(TotalSystematicLowSqr) * SFDefault;
    yErrHigh[j] = TMath::Sqrt(TotalSystematicHighSqr) * SFDefault;
//     cout << "Bin " << j << " " << HwwJetVetoEfficiency[0]->GetXaxis()->GetBinLowEdge(j) << " : " 
//          << SFDefault << " + " << TMath::Sqrt(TotalSystematicHighSqr)*SFDefault << " " << TMath::Sqrt(TotalSystematicLowSqr)*SFDefault << endl;
    cout << "Bin " << j << " " << x[j] << " " << y[j] << " " << yErrLow[j] << " " << yErrHigh[j] << endl;

    jetVetoEfficiencySFHist->SetBinContent(j,SFDefault);
    jetVetoEfficiencySFHist->SetBinError(j,TMath::Max(yErrLow[j],yErrHigh[j]));

  }

  TGraphAsymmErrors *jetVetoEfficiencySF = new TGraphAsymmErrors(201, x, y, xErr, xErr, yErrLow,yErrHigh );
  jetVetoEfficiencySF->SetName("jetVetoEfficiencySF_PDFSystematics");
  jetVetoEfficiencySF->SetTitle("");
  jetVetoEfficiencySF->GetXaxis()->SetTitle("Jet Threshold");
  jetVetoEfficiencySF->GetYaxis()->SetTitle("Jet Veto Efficiency Scale Factor");


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
//   jetVetoEfficiencySF->SetFillColor(kBlue);
//   jetVetoEfficiencySF->SetFillStyle(3001);
//   jetVetoEfficiencySF->GetXaxis()->SetRangeUser(20,50);  
//   jetVetoEfficiencySF->GetYaxis()->SetTitleOffset(1.3);
//   jetVetoEfficiencySF->GetYaxis()->SetRangeUser(0.5, 1.0);
//   jetVetoEfficiencySF->Draw("A3");
//   cv->SaveAs("jetVetoEfficiencySF_PDFSystematics.gif");


  jetVetoEfficiencySFHist->SetMarkerSize(0);
  jetVetoEfficiencySFHist->SetFillColor(kBlue);
  jetVetoEfficiencySFHist->SetFillStyle(3001);
  jetVetoEfficiencySFHist->GetXaxis()->SetRangeUser(20,50);  
  jetVetoEfficiencySFHist->GetYaxis()->SetTitle("Jet Veto Efficiency Ratio");
  jetVetoEfficiencySFHist->GetYaxis()->SetTitleOffset(1.3);
  jetVetoEfficiencySFHist->GetYaxis()->SetRangeUser(0.5, 0.8);
  jetVetoEfficiencySFHist->Draw("E3");
  cv->SaveAs("jetVetoEfficiencySF_PDFSystematics.gif");
  cv->SaveAs("jetVetoEfficiencySF_PDFSystematics.pdf");

 
//   //--------------------------------------------------------------------------------------------------------------
//   // Save Histograms;
//   //============================================================================================================== 
//   TFile *file = new TFile("WWSelectionPlots_SS.root", "RECREATE");
//   file->Close();
//   delete file;

        
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
