//root -l EWKAna/Hww/Acceptance/ZmmJetVetoEfficiency.C+\(\"HwwNtuple_f10-zmm-powheg-c10-v12_noskim_0000.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/Acceptance/ZmmJetVetoEfficiency.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-mu-pr-v2_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Acceptance/ZmmJetVetoEfficiency.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-mu-s17_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Acceptance/ZmmJetVetoEfficiency.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_mu_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Acceptance/ZmmJetVetoEfficiency.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_full_noskim.root\",\"\"\,kTRUE)



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

void ZmmJetVetoEfficiency(const string dataInputFilename, const string Label, Bool_t isData) {   

  gBenchmark->Start("ZeeSelection");

  string label = "";
  if (Label !=  "") label = "_"+Label;
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *leadingJetPtHist = new TH1D((string("leadingJetPtHist") + "_" + label).c_str(), "; Leading Jet p_{T} [GeV/c]; Number of Events", 200, 0, 200);

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
  eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  
  cout << "Total: " << eventTree->GetEntries() << endl;

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
	
    double eventweight = 1;
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(isData && hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
    if (isData && !passHLT(info->triggerBits, info->runNum, kTRUE)) continue;


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
    // Jets
    //********************************************************
    Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    Double_t leadingJetPt = 0;

    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

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
    if (muonArr->GetEntries() < 2) continue;


    //******************************************************************************
    //loop over muon pairs
    //******************************************************************************
    Bool_t passZSelection = kFALSE;

    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
      if ( !(
             mu1->pt > 20.0
             &&
             passMuonCuts(mu1)
             )
        ) continue;
      

      for(Int_t j=i+1; j<muonArr->GetEntries(); j++) {
        const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
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
        leadingJetPtHist->Fill(leadingJetPt);
     }

  } //end loop over data     

  

  delete info;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Compute Jet Veto Efficiency
  //==============================================================================================================
  TH1D *jetVetoEfficiency = (TH1D*)leadingJetPtHist->Clone("jetVetoEfficiency"); 
  jetVetoEfficiency->GetYaxis()->SetTitle("Jet Veto Efficiency");
  jetVetoEfficiency->GetYaxis()->SetTitleOffset(1.5);
  for (int i=1; i<jetVetoEfficiency->GetXaxis()->GetNbins()+1; ++i) {
    Int_t n = leadingJetPtHist->Integral(1,i);
    Int_t d = leadingJetPtHist->Integral(1,leadingJetPtHist->GetXaxis()->GetNbins()+1);
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Int_t errorType = 2;
    mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);

    jetVetoEfficiency->SetBinContent(i,ratio);
    jetVetoEfficiency->SetBinError(i,(errLow+errHigh)/2); 
    cout << jetVetoEfficiency->GetXaxis()->GetBinLowEdge(i) << " : " << n << " / " << d << " = " 
         << ratio << " + " << errHigh << " - " << errLow << endl;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  leadingJetPtHist->Draw();
  cv->SaveAs((string("leadingJetPt") + "_" + label).c_str());

  jetVetoEfficiency->Draw();
  cv->SaveAs((string("JetVetoEfficiency") + "_" + label+".gif").c_str());
 
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  cout << "**************************************************************\n";


        
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
