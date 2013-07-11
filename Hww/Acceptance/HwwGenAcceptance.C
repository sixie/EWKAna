//root -l EWKAna/Hww/Acceptance/HwwGenAcceptance.C+\(\"/home/sixie/hist/HwwAcceptance/ggHWW160GenOnly_default.root\",\"\"\)
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
#include "EWKAna/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================
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

void HwwGenAcceptance(const string inputFilename,  const string Label) 
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
  TH1D *genMet = new TH1D("genMet", ";  Met [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonPtMin = new TH1D("leptonPtMin", ";  Lepton p_{T} [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonEta = new TH1D("leptonEta", "; Lepton #eta; Number of Events", 50, -5, 5);
  TH1D *bosonSystemPt = new TH1D("bosonSystemPt", ";  Higgs p_{T} [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leadingGenJetPt = new TH1D("leadingGenJetPt", "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *nGenJets = new TH1D("nGenJets", "; Mass [GeV/c^{2}]; Number of Events", 10, -0.5, 9.5);
  TH1D *leptonDeltaPhi = new TH1D("leptonDeltaPhi", "; #Delta #phi (l,l); Number of Events", 50, 0, 180);



  Double_t NEventsGenerated = 0;
  Double_t NEventsAccepted = 0;

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
  mithep::TGenInfo *genInfo    = new mithep::TGenInfo();


  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &genInfo); TBranch *genInfoBr = eventTree->GetBranch("Gen");
   

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);
		
    NEventsGenerated++;

    mithep::FourVectorM lepton1;
    mithep::FourVectorM lepton2;
    lepton1.SetCoordinates(genInfo->pt_1, genInfo->eta_1, genInfo->phi_1, 0.51099892e-3 );
    lepton2.SetCoordinates(genInfo->pt_2, genInfo->eta_2, genInfo->phi_2, 0.51099892e-3 );
    mithep::FourVectorM dilepton = lepton1+lepton2;
    Double_t ptMax = genInfo->pt_1; 
    Double_t ptMin = genInfo->pt_2; 
    if (genInfo->pt_1 < genInfo->pt_2) {
      ptMax = genInfo->pt_2; 
      ptMin = genInfo->pt_1; 
    }
    double deltaPhiLeptons = mithep::MathUtils::DeltaPhi(Double_t(genInfo->phi_1), Double_t(genInfo->phi_2));

    if (genInfo->pt_1 > 20.0 && genInfo->pt_2 > 20.0
        && fabs(genInfo->eta_1) < 2.5 && fabs(genInfo->eta_2) < 2.5
        && genInfo->met > 30.0
        && fabs(dilepton.M() - 91.2) > 15.0
        && genInfo->jetpt_1 < 25.0    
        && ptMax > 30
        && ptMin > 25
        && dilepton.M() < 50
        && deltaPhiLeptons < 60
 //        && ptMax > 90
//         && ptMin > 25
//         && dilepton.M() < 300
//         && deltaPhiLeptons < 175
      ) {
      NEventsAccepted++;
    }
  } //end loop over data     

  delete info;
  delete genInfo;

  //--------------------------------------------------------------------------------------------------------------
  // Compute Jet Veto Efficiency
  //==============================================================================================================
    Int_t n = NEventsAccepted;
    Int_t d = NEventsGenerated;
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Int_t errorType = 2;
    mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);

    cout << "Gen Acceptance : " << n << " / " << d << " = " 
         << ratio << " + " << errHigh << " - " << errLow << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
//   TCanvas *cv = new TCanvas("cv","cv", 800,600);
//   leadingGenJetPt->Draw();
//   cv->SaveAs("leadingGenJetPt.gif");

//   jetVetoEfficiency->Draw();
//   cv->SaveAs("jetVetoEfficiency.gif");
 

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
 


        
  gBenchmark->Show("plotZ");       
} 

