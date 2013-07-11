//root -l EWKAna/Hww/Selection/ZeeSelection.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full_TightPlusRecoTriggerSkim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/ZeeSelection.C+\(\"/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root\",\"\"\)



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
      cout << "Cannot get Directory ZeeAnalysisMod from file " << infname << endl;
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

void ZeeSelection(const string dataInputFilename, const string Label) {   

  gBenchmark->Start("ZeeSelection");

  string label = Label;
  if (Label != "") label = "_" + Label;

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_ee = new TH1D("dileptonMass_ee", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_emu = new TH1D("dileptonMass_emu", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_mumu = new TH1D("dileptonMass_mumu", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_BB = new TH1D("dileptonMass_BB", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_BE = new TH1D("dileptonMass_BE", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  TH1D *dileptonMass_EE = new TH1D("dileptonMass_EE", "; Mass [GeV/c^{2}]; Number of Events", 80, 0, 200);
  Double_t Count_ee_BB = 0;
  Double_t Count_ee_BE = 0;
  Double_t Count_ee_EE = 0;
  Double_t Count_mm_BB = 0;
  Double_t Count_mm_BE = 0;
  Double_t Count_mm_EE = 0;
  Double_t Count_em_BB = 0;
  Double_t Count_em_BE = 0;
  Double_t Count_em_EE = 0;
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
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  
    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("json_DCSONLY.txt_160404-163233"); 
  hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(),"Events"); 
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  

  vector<double> RunNumber;
  vector<double> EventsInRun;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
	
    if (ientry % 100000 == 0) cout << ientry << endl;
	
    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
//     if (!passHLT(info->triggerBits, info->runNum, kFALSE)) continue;

    //Find the Run
    Int_t runIndex = -1;
    for (int r=0;r < RunNumber.size(); ++r) {
      if (RunNumber[r] == info->runNum) {
        runIndex = r;
        EventsInRun[r]++;
        break;
      }
    }
    
    if (runIndex == -1) {
      RunNumber.push_back(info->runNum);
      EventsInRun.push_back(1);
    }
    

    Double_t eventweight = info->eventweight * (26.5);
//     eventweight = 1;
    
    

    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear();    
    jetArr->Clear(); 
    electronBr->GetEntry(ientry);
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

     Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

      if (!(jet->pt > 25 && fabs(jet->eta) < 5.0)) continue;
      if (!leadingJet || jet->pt > leadingJet->pt) {
        leadingJet = jet;
      }
      NJets++;
    }


    //******************************************************************************
    //dilepton preselection
    //******************************************************************************
    if (electronArr->GetEntries() < 2) continue;


    //******************************************************************************
    //loop over electron pairs
    //******************************************************************************
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[i]);
      if ( !(
             fabs(ele1->scEta) < 2.5
             &&
             !(fabs(ele1->scEta) > 1.4442 && fabs(ele1->scEta) < 1.566)
             && 
             ele1->scEt > 20.0
             &&
             passElectronCuts(ele1)
             )
        ) continue;
      

      for(Int_t j=i+1; j<electronArr->GetEntries(); j++) {
        const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[j]);
        if ( !(
               fabs(ele2->scEta) < 2.5
               && !(fabs(ele2->scEta) > 1.4442 && fabs(ele2->scEta) < 1.566)
               && ele2->scEt > 20.0
               && passElectronCuts(ele2)
               )
          ) continue;
        
        mithep::FourVectorM lepton1;
        mithep::FourVectorM lepton2;
        lepton1.SetCoordinates(ele1->pt, ele1->eta, ele1->phi, 0.51099892e-3 );
        lepton2.SetCoordinates(ele2->pt, ele2->eta, ele2->phi, 0.51099892e-3 );
        mithep::FourVectorM dilepton = lepton1+lepton2;
 

        dileptonMass->Fill(dilepton.M(),eventweight);
        dileptonMass_ee->Fill(dilepton.M(),eventweight);       
        if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
          dileptonMass_BB->Fill(dilepton.M(),eventweight);
        } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
          dileptonMass_EE->Fill(dilepton.M(),eventweight);
        } else {          
          dileptonMass_BE->Fill(dilepton.M(),eventweight);
        }
        

//         if (ele1->pt > 40 && ele2->pt > 30) {
        if (dilepton.M() > 60 && dilepton.M() < 120) {
          if (fabs(ele1->eta) < 1.5 && fabs(ele2->eta) < 1.5) {
            Count_ee_BB += eventweight;
          } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
            Count_ee_EE += eventweight;
          } else {
            Count_ee_BE += eventweight;
          }
        }
        
      }
    }


  } //end loop over data     

  delete info;
  delete electronArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  dileptonMass->Draw("E");
  cv->SaveAs("dileptonMass_ee.gif");
  dileptonMass_BB->Draw("E");
  cv->SaveAs("dileptonMass_ee_BB.gif");
  dileptonMass_BE->Draw("E");
  cv->SaveAs("dileptonMass_ee_BE.gif");
  dileptonMass_EE->Draw("E");
  cv->SaveAs("dileptonMass_ee_EE.gif");

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
 cout << "Runs Analyzed\n";
  //sort by run number first
  int i, j, flag = 1;    // set flag to 1 to start first pass
  int tempRun, tempNEvents;             // holding variable
  for(i = 0; i < RunNumber.size() && flag; i++) {
    flag = 0;
    for (j=0; j < (RunNumber.size() -1); j++)
    {
      if (RunNumber[j+1] < RunNumber[j])      // ascending order simply changes to <
      { 
        tempRun = RunNumber[j];             // swap elements
        RunNumber[j] = RunNumber[j+1];
        RunNumber[j+1] = tempRun;
        tempNEvents = EventsInRun[j];             // swap elements
        EventsInRun[j] = EventsInRun[j+1];
        EventsInRun[j+1] = tempNEvents;
        flag = 1;               // indicates that a swap occurred.
      }
    }
  }

  double totalTriggers = 0;
  for (int r=0;r < RunNumber.size(); ++r) {
    cout << "Run : " << RunNumber[r] << " : " << EventsInRun[r] << endl;    
    totalTriggers += EventsInRun[r];
  }
  cout << "Total: " << totalTriggers << endl;


  cout << "**************************************************************\n";
  cout << "Event Count : ee final state\n";
  cout << "BB :" << Count_ee_BB << endl;
  cout << "BE :" << Count_ee_BE << endl;
  cout << "EE :" << Count_ee_EE << endl;
  cout << "Total :" << Count_ee_BB + Count_ee_BE + Count_ee_EE << endl;


  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("ZXSelectionPlots.root", "UPDATE");

  file->WriteTObject(dileptonMass, (string(dileptonMass->GetName()) + label).c_str(), "WriteDelete");
  file->WriteTObject(dileptonMass_BB, (string(dileptonMass_BB->GetName()) + label).c_str(), "WriteDelete");
  file->WriteTObject(dileptonMass_BE, (string(dileptonMass_BE->GetName()) + label).c_str(), "WriteDelete");
  file->WriteTObject(dileptonMass_EE, (string(dileptonMass_EE->GetName()) + label).c_str(), "WriteDelete");


  file->Close();
  delete file;

        
        
  gBenchmark->Show("ZeeSelection");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {

  Bool_t isMC = kFALSE;
  Bool_t pass = kFALSE;

  if ( (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) ) pass = kTRUE;



//   if (isMC) {
//     if (triggerBits & kHLT_Mu9) pass = kTRUE;
//     if (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) pass = kTRUE;
//     if (triggerBits & kHLT_Ele15_LW_L1R) pass = kTRUE;
//   } else {
//     if (isMuonData) {
//       if ((runNum >= 136033) && (runNum <= 147116)) {
//         if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
//       } 
//       if ((runNum >= 136033) && (runNum <= 139980)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 140058) && (runNum <= 141882)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 141956) && (runNum <= 144114)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 146428) && (runNum <= 147116)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 999999)) {
//         if ( (triggerBits & kHLT_Mu15) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 148058)) {
//         if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 148819) && (runNum <= 149442)) {
//         if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
//       }
//     } else {
//       //it's electron data
//       if ((runNum >= 136033) && (runNum <= 139980)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 140058) && (runNum <= 141882)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 141956) && (runNum <= 144114)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 146428) && (runNum <= 147116)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 148058)) {
//         if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 148819) && (runNum <= 149442)) {
//         if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
//       }
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
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
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
