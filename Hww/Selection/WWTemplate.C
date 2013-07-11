//root -l EWKAna/Hww/Selection/WWTemplate.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-eg-s17_noskim.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/Selection/WWTemplate.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-el-pr-v2_noskim.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/Selection/WWTemplate.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-mu-pr-v2_noskim.root\",\"\",kTRUE\)
//root -l EWKAna/Hww/Selection/WWTemplate.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-mu-s17_noskim.root\",\"\",kTRUE\)
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
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2);


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

void WWTemplate(const string dataInputFilename,       
                const string Label, Bool_t isMuonData) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    

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
  rlrm.AddJSONFile("Cert_132440-147454_7TeV_StreamExpress_Collisions10_JSON.txt"); 
  

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
		
    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
    if (!passHLT(info->triggerBits, info->runNum, isMuonData)) continue;


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

    Int_t NLeptons = 0;
    vector<Int_t> leptonType;
    vector<Int_t> leptonIndex;
    vector<Double_t> leptonPt;
    vector<Double_t> leptonEta;
    vector<Double_t> leptonPhi;

    Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      if ( (0==0)
           && 
           passElectronCuts(ele)
           &&
           fabs(ele->scEta) < 2.5
           &&
           !(fabs(ele->scEta) > 1.4442 && fabs(ele->scEta) < 1.566)
           && 
           ele->scEt > 20.0
        ) {
        leptonPt.push_back(ele->pt);
        leptonEta.push_back(ele->eta);
        leptonPhi.push_back(ele->phi);
        leptonType.push_back(11);
        leptonIndex.push_back(i);
      }
    }
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
      if ( (0==0)
           &&
           passMuonCuts(mu)
           &&
           fabs(mu->eta) < 2.1
           && 
           mu->pt > 20.0
        ) {
        leptonPt.push_back(mu->pt);
        leptonEta.push_back(mu->eta);
        leptonPhi.push_back(mu->phi);
        leptonType.push_back(13);
        leptonIndex.push_back(i);  
      }
    }
    //sort leptons
    Int_t tempType;
    Int_t tempIndex;
    Double_t tempPt;
    Double_t tempEta;
    Double_t tempPhi;
    for (int l=0; l<leptonIndex.size(); l++) {
      for (int k=0; k < leptonIndex.size() - 1; k++) {
        if (leptonPt[k+1] > leptonPt[k]) {
          tempType = leptonType[k];
          tempIndex = leptonIndex[k];
          tempPt = leptonPt[k];
          tempEta = leptonEta[k];
          tempPhi = leptonPhi[k];
          
          leptonType[k] = leptonType[k+1];
          leptonIndex[k] = leptonIndex[k+1];
          leptonPt[k] = leptonPt[k+1];
          leptonEta[k] = leptonEta[k+1];
          leptonPhi[k] = leptonPhi[k+1];

          leptonType[k+1] = tempType;
          leptonIndex[k+1] = tempIndex;
          leptonPt[k+1] = tempPt;
          leptonEta[k+1] = tempEta;
          leptonPhi[k+1] = tempPhi;
          
        }
      }
    }
    
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

      if (!(jet->et > 25 && fabs(jet->eta) < 5.0)) continue;
      if (!leadingJet || jet->et > leadingJet->et) {
        leadingJet = jet;
      }
      NJets++;
    }


    //******************************************************************************
    //dilepton preselection
    //******************************************************************************
    if (leptonPt.size() < 2) continue;
    if (!(leptonPt[0] > 20.0 && leptonPt[1] > 10.0)) continue;
    if (leptonPt.size() >= 3 || leptonPt[2]>10.0) continue;



    Int_t finalState = -1;
    if (leptonType[0] == 11 && leptonType[1] == 11) {
      finalState = 0;
    } else if (leptonType[0] == 13 && leptonType[1] == 13) {
      finalState = 1;
    } else if (leptonType[0] == 11 && leptonType[1] == 13) {
      finalState = 2;
    } else if (leptonType[0] == 13 && leptonType[1] == 11) {
      finalState = 3;
    }

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

    double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(leptonPhi[0], 
                                                         leptonPhi[1])* 180.0 / TMath::Pi();    
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
    const int nCuts = 7;
    bool passCut[nCuts] = {false, false, false, false, false, false, false};
  
    if(lepton1.Pt() >  20.0 &&
       lepton2.Pt() >= 20.0) passCut[0] = true;
    
    if(met.Pt()    > 20.0)               passCut[1] = true;
  
    if(dilepton.M() > 12.0)            passCut[2] = true;
   
    if(NJets     < 1)              passCut[5] = true;

//     if(CleanLeptons->GetEntries() == 2 &&
//        SoftMuons->GetEntries() == 0)      passCut[6] = true;

    if (finalState == 0 || finalState == 1){ // mumu/ee
      if(fabs(dilepton.M()-91.1876)   > 15.0)   passCut[3] = true;
      if(METdeltaPhilEt > 35) passCut[4] = true;
    }
    else if(finalState == 2 ||finalState == 3 ) { // emu
      passCut[3] = true;
      if(METdeltaPhilEt > 20) passCut[4] = true;
    }
 
    //*********************************************************************************************
    //Fill Histograms
    //*********************************************************************************************
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

      eventDump(eventListFile, info->runNum, info->lumiSec, info->evtNum, dilepton.M(),
                lepton1.Pt(), lepton1.Eta(), lepton1.Phi(), lepton2.Pt(), lepton2.Eta(), lepton2.Phi());


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
  cout << "**************************************************************\n";
  cout << "Event Count : ee final state\n";
  cout << "BB :" << Count_ee_BB << endl;
  cout << "BE :" << Count_ee_BE << endl;
  cout << "EE :" << Count_ee_EE << endl;
  cout << "Total :" << Count_ee_BB + Count_ee_BE + Count_ee_EE << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : mm final state\n";
  cout << "BB :" << Count_mm_BB << endl;
  cout << "BE :" << Count_mm_BE << endl;
  cout << "EE :" << Count_mm_EE << endl;
  cout << "Total :" << Count_mm_BB + Count_mm_BE + Count_mm_EE << endl;
  cout << "**************************************************************\n";
  cout << "Event Count : em final state\n";
  cout << "BB :" << Count_em_BB << endl;
  cout << "BE :" << Count_em_BE << endl;
  cout << "EE :" << Count_em_EE << endl;
  cout << "Total :" << Count_em_BB + Count_em_BE + Count_em_EE << endl;



        
  gBenchmark->Show("plotZ");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;

  if (isMuonData) {
    if ((runNum >= 132440) && (runNum <= 147119)) {
      if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
    } 
    if ((runNum >= 132440) && (runNum <= 135058)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 135059) && (runNum <= 140401)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140042) && (runNum <= 141900)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 141901) && (runNum <= 146427)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 146428) && (runNum <= 147119)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15_v1) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
    }
  } else {
    //it's electron data
    if ((runNum >= 132440) && (runNum <= 135058)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 135059) && (runNum <= 140041)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140042) && (runNum <= 141900)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 141901) && (runNum <= 146427)) {
      if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 146428) && (runNum <= 147119)) {
      if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147120) && (runNum <= 999999)) {
      if ( !(triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
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
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.06
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
