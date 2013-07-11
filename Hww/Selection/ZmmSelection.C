//root -l EWKAna/Hww/Selection/ZmmSelection.C+\(\"/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_f10-zmm-powheg-c10-v12_noskim_normalized.root\",\"\"\)
//root -l EWKAna/Hww/Selection/ZmmSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-mu-pr-v2_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/ZmmSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-mu-s17_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/ZmmSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_mu_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Selection/ZmmSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_full_noskim.root\",\"\"\)



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

void ZmmSelection(const string dataInputFilename, const string Label) {   

  gBenchmark->Start("ZeeSelection");

  
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
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  

  infile = new TFile(dataInputFilename.c_str());
  assert(infile);

    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("Cert_TopNov5_Merged_135821-149442_allPVT.txt"); 

  

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(),"Events"); 
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  
  cout << "Total: " << eventTree->GetEntries() << endl;
  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
	
//     Double_t eventweight = info->eventweight*3.1;
    double eventweight = 1;
//    cout << eventweight << endl;
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
    if (!passHLT(info->triggerBits, info->runNum, kTRUE)) continue;


    //********************************************************
    // Load the branches
    //********************************************************
    muonArr->Clear();    
    jetArr->Clear(); 
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

     Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

      if (!(jet->pt > 25 && fabs(jet->pt) < 5.0)) continue;
      if (!leadingJet || jet->pt > leadingJet->pt) {
        leadingJet = jet;
      }
      NJets++;
    }


    //******************************************************************************
    //dilepton preselection
    //******************************************************************************
    if (muonArr->GetEntries() < 2) continue;


    //******************************************************************************
    //loop over muon pairs
    //******************************************************************************
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
      if ( !(
//              fabs(mu1->eta) < 2.1
//              && 
             mu1->pt > 20.0
             &&
             passMuonCuts(mu1)
             )
        ) continue;
      

      for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
        const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
        if ( !(
//                fabs(mu2->eta) < 2.1
//                && 
               mu2->pt > 20.0
               && passMuonCuts(mu2)
               )
          ) continue;
        
        mithep::FourVectorM lepton1;
        mithep::FourVectorM lepton2;
        lepton1.SetCoordinates(mu1->pt, mu1->eta, mu1->phi, 0.51099892e-3 );
        lepton2.SetCoordinates(mu2->pt, mu2->eta, mu2->phi, 0.51099892e-3 );
        mithep::FourVectorM dilepton = lepton1+lepton2;
 

        dileptonMass->Fill(dilepton.M());
        dileptonMass_mumu->Fill(dilepton.M());
        if (dilepton.M() > 60 
            //  && dilepton.M() < 120
          ) {
          if (fabs(lepton1.Eta()) < 1.5 && fabs(lepton2.Eta()) < 1.5) {
            Count_mm_BB += eventweight;
          } else if (fabs(lepton1.Eta()) > 1.5 && fabs(lepton2.Eta()) > 1.5) {
            Count_mm_EE += eventweight;
          } else {
            Count_mm_BE += eventweight;
          }
        }
        
      }
    }


  } //end loop over data     

  delete info;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  dileptonMass->Draw("E");
  cv->SaveAs("dileptonMass_mm.gif");

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  cout << "**************************************************************\n";
  cout << "Event Count : mm final state\n";
  cout << "BB :" << Count_mm_BB << endl;
  cout << "BE :" << Count_mm_BE << endl;
  cout << "EE :" << Count_mm_EE << endl;
  cout << "Total :" << Count_mm_BB + Count_mm_BE + Count_mm_EE << endl;



        
  gBenchmark->Show("ZeeSelection");       
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
      if ((runNum >= 147196) && (runNum <= 149442)) {
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



Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->typeBits & kTracker
         && mu->nTkHits > 10
         && mu->muNchi2 < 10.0
//         && (mu->qualityBits & kGlobalMuonPromptTight)
//         && fabs(mu->d0) < 0.02
//         && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

//         && (mu->nSeg > 1 || mu->nMatch > 1 )
//         && (mu->nPixHits > 0)


//          && (mu->nPixHits > 0)
//          && (mu->nValidHits > 0 
// //               && mu->nMatch > 1
//            )
//         && mu->trkIso03 < 10
        
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
