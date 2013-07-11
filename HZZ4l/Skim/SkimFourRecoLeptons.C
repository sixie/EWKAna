//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-dmu-pr.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-dmu-pr_FourRecoLeptonSkim.root\"\)
//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-del-pr.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-del-pr_FourRecoLeptonSkim.root\"\)
//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-mueg-pr.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-mueg-pr_FourRecoLeptonSkim.root\"\)
//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-sel-pr.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-sel-pr_FourRecoLeptonSkim.root\"\)
//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-smu-pr-v1.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-smu-pr-v1_FourRecoLeptonSkim.root\"\)
//root -l -b -q EWKAna/HZZ4l/Skim/SkimFourRecoLeptons.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-smu-pr-v2.txt\",\"/home/sixie/hist/HwwAnalysis/2011Data/HZZ4lAnalysis_r11a-smu-pr-v2_FourRecoLeptonSkim.root\"\)



#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
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

Double_t getNEvents(string filename) {

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  if (!dir) {
    cout << "Could not find directory AnaFwkMod"
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }  
  return hist->Integral();
}

// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimFourRecoLeptons(string inputFilename, string outputFilename) 
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
  UInt_t nEventsTotal = 0;

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


  // list input ntuple files to be skimmed
  vector<string> infilenames;  
  ifstream ifs;
  ifs.open(inputFilename.c_str()); 
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();

  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Skimming " << infilenames[ifile] << "..." << endl;
    TTree *eventTree = getTreeFromFile(infilenames[ifile].c_str(),"Events"); 
    nEventsTotal += getNEvents(infilenames[ifile].c_str()); 
    assert(eventTree);
    
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("Photon",    &photonArr);  TBranch *photonBr   = eventTree->GetBranch("Photon");
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
      infoBr->GetEntry(ientry);
      muonArr->Clear();     muonBr->GetEntry(ientry);
      electronArr->Clear(); electronBr->GetEntry(ientry);      
      jetArr->Clear(); jetBr->GetEntry(ientry);
      photonArr->Clear(); photonBr->GetEntry(ientry);
     
      if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

      nInputEvts++;
      
      Bool_t keep = kTRUE;
    
      Int_t NRecoLeptons = 0;
      Int_t NMuons = 0;
      Int_t NDenominatorMuons = 0;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
        if ( fabs(mu->eta) < 2.4 && mu->pt > 10.0) {
          NRecoLeptons++;
        }
      }
      Int_t NElectrons = 0;
      Int_t NDenominatorElectrons = 0;
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        if ( fabs(ele->eta) < 2.5 && ele->pt > 10.0 ) {
          NRecoLeptons++;
        }
      }

     //One Tight AND Two Loose
     if (!(NRecoLeptons >= 4 )) {
       keep = kFALSE;
     }

      if(keep) {
        outEventTree->Fill();
        nPassEvts++;
      }
    }
  }
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
    
  std::cout << outputFilename << " created!" << std::endl;
  std::cout << " >>> Total Number of Events: " << nEventsTotal << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;

  gBenchmark->Show("SkimNtuples");
}  


