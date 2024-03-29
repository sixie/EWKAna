//================================================================================================
//
// Select probes for single muon trigger efficiency with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "EWKAna/Utils/LeptonIDCuts.hh"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

// structure for output ntuple
#include "EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectMuEGMuonLeadingLegTrigEffTP(const TString conf,              // input file
                                    const TString outputDir,         // output directory
                                    Int_t RunRange = 0,              // Run Range
                                    const Bool_t  matchGen = kFALSE  // match to generator muons
  ) {
  gBenchmark->Start("selectSingleMuEffTP");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //*****************************************************************************************
  //Setup MVA
  //*****************************************************************************************
  mithep::ElectronIDMVA *electronIDMVAWithIPInfo = new mithep::ElectronIDMVA();
  electronIDMVAWithIPInfo->Initialize("BDTG method",
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_WithIPInfo_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_WithIPInfo_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_WithIPInfo_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_WithIPInfo_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_WithIPInfo_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_WithIPInfo_BDTG.weights.xml",
                              mithep::ElectronIDMVA::kWithIPInfo);

  
  // mass region
  Double_t massLo;
  Double_t massHi;

  Double_t lumi;              // luminosity (pb^-1)
  
  vector<TString>  fnamev;    // sample files 
  vector<Int_t>    typev;     // dataset type 
  vector<Double_t> xsecv;     // per file cross section
  vector<TString>  jsonv;     // per file JSON file

  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> massLo >> massHi; 
      
    } else if(state==1) {  // define data sample
      string fname;
      Int_t type;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> type >> xsec >> json;
      fnamev.push_back(fname);
      typev.push_back(type);
      xsecv.push_back(xsec);
      jsonv.push_back(json);        
    }
  }
  ifs.close();
  
  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F");

  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr= new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr     = new TClonesArray("mithep::TJet");
  
  // loop over files  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    Bool_t hasJSON = kFALSE;
    mithep::RunLumiRangeMap rlrm;
    if(jsonv[ifile].CompareTo("NONE")!=0) { 
      hasJSON = kTRUE;
      rlrm.AddJSONFile(jsonv[ifile].Data()); 
    }
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info", &info);            TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet", &jetArr);         TBranch *jetBr      = eventTree->GetBranch("PFJet");
    TBranch *genBr = 0;
    if(matchGen) {
      eventTree->SetBranchAddress("Gen", &gen);
      genBr = eventTree->GetBranch("Gen");
    }
    
    // Determine maximum number of events to consider
    // *** CASES ***
    // <> lumi < 0 => use all events in the sample
    // <> xsec = 0 => for data (use all events)
    const Double_t xsec = xsecv[ifile];
    Double_t weight = 1;
    if(lumi>0) { 
      if(xsec>0) { weight = lumi*xsec/(Double_t)eventTree->GetEntries(); }      
    }

    // loop over events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      infoBr->GetEntry(ientry);
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      //********************************************************************************
      // check for certified runs
      //********************************************************************************
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

      //********************************************************************************
      //RunRange Requirement
      //********************************************************************************
      Bool_t passRunRange = kTRUE;
      if (RunRange == 1) {
        if (!(info->runNum <= 170053)) passRunRange = kFALSE;
      } else if (RunRange == 2) {
        if (!(info->runNum > 173198)) passRunRange = kFALSE;
      }
      if (!passRunRange) continue;

      //********************************************************************************
      // Tag & Probe Trigger Requirement
      //********************************************************************************
      ULong_t  tnpTrigger = kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT ;
      if(!(info->triggerBits & tnpTrigger)) continue;      
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      if(matchGen) genBr->GetEntry(ientry);
      
      muonArr->Clear();
      electronArr->Clear();
      jetArr->Clear();
      muonBr->GetEntry(ientry);
      electronBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
        
	if(matchGen) {
	  Bool_t match1 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5);
	  Bool_t match2 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5);
	  if(!match1 && !match2)
	    continue;
	}
	
        //********************************************************************************
        // Tag Requirement
        //********************************************************************************
	if(tag->pt        < 20)  continue;
	if(fabs(tag->scEta) > 2.5) continue;
        if (!passElectronMVA(tag, electronIDMVAWithIPInfo->MVAValue(
                               tag->pt,tag->scEta,
                               tag->sigiEtaiEta, 
                               tag->deltaEtaIn,
                               tag->deltaPhiIn, 
                               tag->HoverE,
                               tag->d0,
                               tag->dz, 
                               tag->fBrem,
                               tag->EOverP,
                               tag->ESeedClusterOverPout,
                               TMath::Sqrt(tag->sigiPhiiPhi),
                               tag->nBrem,
                               (1.0/(tag->scEt * TMath::CosH(tag->scEta)) - 1/tag->p), 
                               tag->ESeedClusterOverPIn,
                               tag->ip3d,
                               tag->ip3dSig ), 
                             2)) continue;
   
        if(
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT))
          &&
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT))
          &&
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele52_CaloIdVT_TrkIdT))

          &&
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele65_CaloIdVT_TrkIdT))
          &&
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele80_CaloIdVT_TrkIdT))
          &&
          !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            (tag->hltMatchBits & kHLTObject_Ele32_WP70))
          )
          continue;   
        
	
	for(Int_t j=0; j<muonArr->GetEntriesFast(); j++) {
	  
	  const mithep::TMuon *probe = (mithep::TMuon*)((*muonArr)[j]);
	  if(probe->q == tag->q) continue;
          if (toolbox::deltaR(tag->eta, tag->phi, probe->eta, probe->phi) < 0.5) continue;

	  if(matchGen) {
	    Bool_t match1 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5);
	    Bool_t match2 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5);
	    if(!match1 && !match2)
	      continue;
	  }
          if(!passMuonID(probe))     continue;


         //ttbar selection : 2 jets, at least 1 b-tagged
          Int_t NJets = 0;
          Int_t NBTags = 0;
          for(Int_t jetIndex=0; jetIndex<jetArr->GetEntriesFast(); ++jetIndex) {
            const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[jetIndex]);

//             cout << "Jet " << jetIndex << " : " << jet->pt << " " << jet->eta << " " << jet->phi << " : " 
//                  << toolbox::deltaR(jet->eta, jet->phi, tag->eta, tag->phi) << " " 
//                  << toolbox::deltaR(jet->eta, jet->phi, probe->eta, probe->phi) << " " 
//                  << jet->TrackCountingHighEffBJetTagsDisc << " "
//                  << " ----------> : " << ( (toolbox::deltaR(jet->eta, jet->phi, tag->eta, tag->phi) > 0.5) && 
//                                            (toolbox::deltaR(jet->eta, jet->phi, probe->eta, probe->phi) > 0.5) &&
//                                            (jet->pt > 30.0) && (jet->TrackCountingHighEffBJetTagsDisc > 3.5)) << " "
//                  << endl;
            if (toolbox::deltaR(jet->eta, jet->phi, tag->eta, tag->phi) < 0.5) continue;
            if (toolbox::deltaR(jet->eta, jet->phi, probe->eta, probe->phi) < 0.5) continue;
            if (jet->pt > 30.0) {
              NJets++;
              if (jet->TrackCountingHighEffBJetTagsDisc > 3.5) { NBTags++; }
            }
          }
          if (NJets != 2 && NBTags >= 1) continue;
	  
	  nProbes += weight;
	  
	  Bool_t pass = kFALSE;

	  if(info->runNum >= 150000 && info->runNum <= 170053) pass = ( (info->triggerBits & kHLT_Mu17_Ele8_CaloIdL)         
                                                                               // && (probe->hltMatchBits & kHLTObject_Ele17_CaloIdL)
            );	  
	  if(info->runNum >  170053 && info->runNum <= 173199) pass = ( (info->triggerBits & kHLT_Mu17_Ele8_CaloIdL)
                                                                               //&& (probe->hltMatchBits & kHLTObject_Ele17_CaloIdT_CaloIsoVL)
            );
	  if(info->runNum >  173199                          ) pass = ( (info->triggerBits & kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL)
                                                                               //&& (probe->hltMatchBits & kHLTObject_Ele17_CaloIdT_CaloIsoVL)
            );
	  if(typev[ifile]==eMC) 
	    pass = ( ((info->triggerBits & kHLT_Mu17_Ele8_CaloIdL) && (probe->hltMatchBits & kHLTObject_Ele8_CaloIdL)) 
                     || 
                     ((info->triggerBits & kHLT_Mu8_Ele17_CaloIdL) && (probe->hltMatchBits & kHLTObject_Ele17_CaloIdL)) 
              ); 
	  
          // Fill tree
	  data.mass    = 90;
	  data.pt      = probe->pt;
          data.eta     = probe->eta;
          data.phi     = probe->phi;
          data.weight  = weight;
	  data.q       = probe->q;
          data.npv     = info->nPV0;
	  data.npu     = info->nPUEvents;
	  data.pass    = (pass) ? 1 : 0;
	  data.runNum  = info->runNum;
	  data.lumiSec = info->lumiSec;
	  data.evtNum  = info->evtNum;
          data.rho     = info->PileupEnergyDensity;
	  outTree->Fill();
	}
      }
    }
    delete infile;
    infile=0, eventTree=0;    
  }
  delete info;
  delete gen;
  delete muonArr;
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  if(lumi>0) {
    cout << " L_int = " << lumi << "/pb" << endl;
    cout << endl;
  }
  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectSingleMuEffTP"); 
}
