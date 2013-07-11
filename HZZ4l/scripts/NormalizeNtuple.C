//root -l EWKAna/Hww/scripts/NormalizeNtuple.C+\(\"/home/sixie/hist/WWAnalysis/merged/WWAnalysis_p10-stop-mg-v26_noskim.root\",\"p10-stop-mg-v26\",\"/home/sixie/hist/WWAnalysis/normalized/WWAnalysis_p10-stop-mg-v26_noskim.root\"\)

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
#include "TH1F.h"
#include "MitAna/Utils/interface/SimpleTable.h"
#endif

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"

//=== FUNCTION DECLARATIONS ======================================================================================
Double_t getNormalizationWeight(string filename, string datasetName);

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
void NormalizeNtuple(string inputFilename, string datasetName, string outputFilename) 
{
  gBenchmark->Start("NormalizeNtuple");
    
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

  Double_t normalizationWeight = getNormalizationWeight(inputFilename, datasetName);

  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    muonArr->Clear();     muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);
    photonArr->Clear(); photonBr->GetEntry(ientry);
    jetArr->Clear();      jetBr->GetEntry(ientry);
      
    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
     
    info->eventweight = info->eventweight * normalizationWeight;
   outEventTree->Fill();
    nPassEvts++;
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
  
  gBenchmark->Show("NormalizeNtuple");
}  

 
//--------------------------------------------------------------------------------------------------
// Get Total Number of Events in the sample
//--------------------------------------------------------------------------------------------------
Double_t getNormalizationWeight(string filename, string datasetName) {
  // Get Normalization Weight Factor

  //Get Number of Events in the Sample
  Double_t NEvents = 1;
  TFile *file = new TFile(filename.c_str(),"READ");
  TDirectory *dir = 0;
  TH1F *hist = 0;
  if (!file) {
    cout << "Could not open file " << filename << endl;
  } else {
    
    hist = (TH1F*)file->Get("hDAllEvents");
    if (!hist) {
      dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
      if (!dir) {
        cout << "Could not find directory AnaFwkMod"
             << " in file " << filename << endl;
        delete file;
      } else {
        
        hist = (TH1F*)dir->Get("hDAllEvents");
        if (!hist) {
          cout << "Could not find histogram hDEvents in directory AnaFwkMod"
               << " in file " << filename << endl;
          delete dir;
          delete file;
        }
      }
    }
  }

  if (hist) {
    NEvents = hist->Integral();
  }

  //Get CrossSection
  Double_t Weight = 1.0;
  if (datasetName.find("data") == string::npos) {
    mithep::SimpleTable xstab("/home/sixie/CMSSW_analysis/src/MitPhysics/data/xs.dat");
    cerr << "Use xs table: " << "$CMSSW_BASE/src/MitPhysics/data/xs.dat" << endl;
    Double_t CrossSection = xstab.Get(datasetName.c_str());
    Weight = CrossSection / NEvents;
    cerr << datasetName << " : " << CrossSection << " " << NEvents << " " << Weight << endl;
  }

//   if (dir) delete dir;
//   if (file) delete file;
  return Weight;

}
