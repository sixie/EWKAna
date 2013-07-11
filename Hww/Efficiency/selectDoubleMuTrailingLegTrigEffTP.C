//================================================================================================
//
// Select probes for double muon trigger efficiency with Tag&Probe method
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
#include <TRandom.h>                // random number generator
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

// structure for output ntuple
#include "EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectDoubleMuTrailingLegTrigEffTP(const TString conf,              // input file
                         const TString outputDir,         // output directory
                         Int_t RunRange = 0,              // Run Range
		         const Bool_t  matchGen = kFALSE  // match to generator muons
) {
  gBenchmark->Start("selectDoubleMuEffTP");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
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

  TRandom rg;

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
    eventTree->SetBranchAddress("Info", &info);    TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
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
        if (!(info->runNum <= 164237)) passRunRange = kFALSE;
      } else if (RunRange == 2) {
        if (!(info->runNum > 164237 && info->runNum <= 178380)) passRunRange = kFALSE;
      } else if (RunRange == 3) {
        if (!(info->runNum > 178380 && info->runNum <= 999999)) passRunRange = kFALSE;
      }

      if (!passRunRange) continue;

      //********************************************************************************
      // Tag & Probe Trigger Requirement
      //********************************************************************************
      ULong_t  tnpTrigger = kHLT_Mu15 | kHLT_Mu30 | kHLT_IsoMu17 | kHLT_IsoMu24 | kHLT_IsoMu30 ;
      if(!(info->triggerBits & tnpTrigger)) continue;                  
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      if(matchGen) genBr->GetEntry(ientry);
      
      muonArr->Clear();
      muonBr->GetEntry(ientry);
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
        const mithep::TMuon *mu1 = (mithep::TMuon*)((*muonArr)[i]);
        
	if(matchGen) {
	  Bool_t match1 = (toolbox::deltaR(mu1->eta, mu1->phi, gen->eta_1, gen->phi_1) < 0.5);
	  Bool_t match2 = (toolbox::deltaR(mu1->eta, mu1->phi, gen->eta_2, gen->phi_2) < 0.5);
	  if(!match1 && !match2)
	    continue;
	}
		
	if(fabs(mu1->eta) > 2.4) continue;
	if(!passMuonID(mu1))     continue;
    
        const Double_t m = 0.105658369;
	TLorentzVector vmu1;
	vmu1.SetPtEtaPhiM(mu1->pt, mu1->eta, mu1->phi, m);
	
	for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) { 
	  const mithep::TMuon *mu2 = (mithep::TMuon*)((*muonArr)[j]);
	  if(mu2->q == mu1->q)     continue;
	  	  
	  if(matchGen) {
	    Bool_t match1 = (toolbox::deltaR(mu2->eta, mu2->phi, gen->eta_1, gen->phi_1) < 0.5);
	    Bool_t match2 = (toolbox::deltaR(mu2->eta, mu2->phi, gen->eta_2, gen->phi_2) < 0.5);
	    if(!match1 && !match2)
	      continue;
	  }
	  
	  if(fabs(mu2->eta) > 2.4) continue;
	  if(!passMuonID(mu2)) continue;
	 	  
	  TLorentzVector vmu2;
	  vmu2.SetPtEtaPhiM(mu2->pt, mu2->eta, mu2->phi, m);
	  
	  TLorentzVector vdimuon = vmu1 + vmu2;
	  if((vdimuon.M()<massLo) || (vdimuon.M()>massHi)) continue;
	  
          //********************************************************************************
          // Randomly choose the probe
          //********************************************************************************
	  const mithep::TMuon *tag=0, *probe=0;
	  if(rg.Binomial(1,0.5)<1) { tag = mu1; probe = mu2; }
	  else                     { tag = mu2; probe = mu1; }
	  if(tag->pt < 20)  continue;
	  
	  if(!((info->triggerBits & kHLT_IsoMu17) && (tag->hltMatchBits & kHLTObject_IsoMu17)) &&
	     !((info->triggerBits & kHLT_IsoMu24) && (tag->hltMatchBits & kHLTObject_IsoMu24)) &&
	     !((info->triggerBits & kHLT_IsoMu30) && (tag->hltMatchBits & kHLTObject_IsoMu30)) &&
	     !((info->triggerBits & kHLT_Mu15)    && (tag->hltMatchBits & kHLTObject_Mu15)) &&
	     !((info->triggerBits & kHLT_Mu30)    && (tag->hltMatchBits & kHLTObject_Mu30)) &&
	     !((info->triggerBits & kHLT_Mu30)    && (tag->hltMatchBits & kHLTObject_Mu40))
            )
	   continue;   
	  
	  nProbes += weight;
	  
	  Bool_t pass = kFALSE;
	  if(info->runNum >= 150000 && info->runNum <= 164237) 
	    pass = (info->triggerBits & kHLT_DoubleMu7) && (probe->hltMatchBits & kHLTObject_Mu7);
	  if(info->runNum >= 165085 && info->runNum <= 178380)
	    pass = (info->triggerBits & kHLT_Mu17_Mu8) && (probe->hltMatchBits & (kHLTObject_Mu8));	  
	  if(info->runNum >= 178381 )
	    pass = (info->triggerBits & kHLT_Mu17_Mu8) && (probe->hltMatchBits & (kHLTObject_Mu8));	  
          if(typev[ifile]==eMC)
	    pass = (info->triggerBits & kHLT_DoubleMu7) && (probe->hltMatchBits & kHLTObject_Mu7);
	  
	  // Fill tree
	  data.mass    = vdimuon.M();
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
      
  gBenchmark->Show("selectDoubleMuEffTP"); 
}
