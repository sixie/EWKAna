//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-dmu-pr-v4_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-dmu-pr-v4_TightPlusRecoTriggerSkim.root\",kFALSE,0\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr-v4_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-pr-v4_TightPlusRecoTriggerSkim.root\",kFALSE,1\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-mueg-pr-v4_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-mueg-pr-v4_TightPlusRecoTriggerSkim.root\",kFALSE,2\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v4_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v4_TightPlusRecoTriggerSkim.root\",kFALSE,3\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-sel-pr-v4_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-sel-pr-v4_TightPlusRecoTriggerSkim.root\",kFALSE,4\)

//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-dmu-m10-v1_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-dmu-m10-v1_TightPlusRecoTriggerSkim.root\",kFALSE,0\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-m10-v1_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-del-m10-v1_TightPlusRecoTriggerSkim.root\",kFALSE,1\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-mueg-m10-v1_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-mueg-m10-v1_TightPlusRecoTriggerSkim.root\",kFALSE,2\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-m10-v1_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-m10-v1_TightPlusRecoTriggerSkim.root\",kFALSE,3\)
//root -l -b -q EWKAna/Hww/Skim/SkimTrigger.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-sel-m10-v1_TightPlusRecoNoTriggerSkim.root\",\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-sel-m10-v1_TightPlusRecoTriggerSkim.root\",kFALSE,4\)

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#endif

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"

//=== FUNCTION DECLARATIONS ======================================================================================
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMC, Int_t PD);

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


// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimTrigger(string inputFilename, string outputFilename, Bool_t isMC, Int_t PD ) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();

  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr     = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr    = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr     = new TClonesArray("mithep::TPhoton");

  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  
  TFile* outfile = new TFile(outputFilename.c_str(), "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",     &info);
  outEventTree->Branch("Muon",     &muonArr);
  outEventTree->Branch("Electron", &electronArr);
  outEventTree->Branch("PFJet",    &jetArr);
  outEventTree->Branch("Photon",   &photonArr);

  TTree *eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
   assert(eventTree);
    
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr    = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("Photon",    &photonArr);  TBranch *photonBr   = eventTree->GetBranch("Photon");

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    muonArr->Clear(); muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);      
    jetArr->Clear(); jetBr->GetEntry(ientry);
    photonArr->Clear(); photonBr->GetEntry(ientry);

    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
      
    Bool_t keep = kTRUE;
    if (!passHLT(info->triggerBits, info->runNum, isMC, PD)) keep = kFALSE;
    
    if(keep) {
      outEventTree->Fill();
      nPassEvts++;
    }
  }
  
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
    
  std::cout << outputFilename << " created!" << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;
  
  gBenchmark->Show("SkimNtuples");
}  

Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMC , Int_t PD) {


  Bool_t pass = kFALSE;
  if (isMC) {
  } else {
    if (PD == 0) {
      //doubleMu
      if ((runNum >= 160329) && (runNum <= 999999)) {
        if ( (triggerBits & kHLT_DoubleMu7) ) pass = kTRUE;
      } 
      
    } else if (PD == 1) {
      //doubleEle
      if ((runNum >= 160329) && (runNum <= 999999)) {
        if ( !(triggerBits & kHLT_DoubleMu7 ) 
             && (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) 
          ) pass = kTRUE;
      } 
    } else if (PD == 2) {
      //MuEG
      if ((runNum >= 160329) && (runNum <= 999999)) {
        if ( !(triggerBits & kHLT_DoubleMu7)
             && !(triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL)
             && ( (triggerBits & kHLT_Mu17_Ele8_CaloIdL)  
                  || (triggerBits & kHLT_Mu8_Ele17_CaloIdL)
               )
          ) pass = kTRUE;
      }
    } else if (PD == 3) {
      //SingleMu
      if ((runNum >= 160329) && (runNum <= 999999)) {
        if ( !(triggerBits & kHLT_DoubleMu7)
             && !(triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL)
             && !( (triggerBits & kHLT_Mu17_Ele8_CaloIdL)  
                  || (triggerBits & kHLT_Mu8_Ele17_CaloIdL)
               )
             && ( triggerBits & kHLT_Mu15 || triggerBits & kHLT_IsoMu17 )
          ) pass = kTRUE;
      }
    } else if (PD == 4) {
        //SingleEle
      if ((runNum >= 160329) && (runNum <= 999999)) {
        if ( !(triggerBits & kHLT_DoubleMu7)
             && !(triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL)
             && !( (triggerBits & kHLT_Mu17_Ele8_CaloIdL)  
                  || (triggerBits & kHLT_Mu8_Ele17_CaloIdL)
               )
             && !( triggerBits & kHLT_Mu15 || triggerBits & kHLT_IsoMu17 )
             && ( triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT)
          ) pass = kTRUE;
      }
    } else {
      cout << "PD == " << PD << " is not recognized.\n";
    }
  }
  return pass;

}

