//================================================================================================
//
// Select probes for double electron trigger efficiency with Tag&Probe method
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
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

// structure for output ntuple
#include "EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectDoubleElectronTrailingLegTrigEffTP(const TString conf,            // input file
                               const TString outputDir,         // output directory
                               Int_t RunRange = 0,              // Run Range
		               const Bool_t  matchGen = kFALSE  // match to generator electrons
) {
  gBenchmark->Start("selectDoubleElectronTrigEffTP");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  //*****************************************************************************************
  //Setup MVA
  //*****************************************************************************************
  mithep::ElectronIDMVA *electronIDMVAIDIsoCombined = new mithep::ElectronIDMVA();
  electronIDMVAIDIsoCombined->Initialize("BDTG method",
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml", 
                              "/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml",
                              mithep::ElectronIDMVA::kIDIsoCombined);

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
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo   *gen   = new mithep::TGenInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  
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
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
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

      Double_t rho = 0;
      if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;

      //********************************************************************************
      //RunRange Requirement
      //********************************************************************************
      Bool_t passRunRange = kTRUE;
      if (RunRange == 1) {
        if (!(info->runNum <= 170053)) passRunRange = kFALSE;
      } else if (RunRange == 2) {
        if (!(info->runNum > 170053)) passRunRange = kFALSE;
      }
      if (!passRunRange) continue;


      //********************************************************************************
      // Tag & Probe Trigger Requirement
      //********************************************************************************
      ULong_t tnpTrigger = 0;
      if(typev[ifile]==eDiEl) {
        tnpTrigger = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 |
		     kHLT_Ele32_CaloIdL_CaloIsoVL_SC17 |
	             kHLT_Ele17_CaloIdL_CaloIsoVL;
      } else if(typev[ifile]==eMC) {
        tnpTrigger = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30;
      }
      if(!(info->triggerBits & tnpTrigger)) continue;        
      

      if(matchGen) genBr->GetEntry(ientry);
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);
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
	if(tag->pt          < 20)  continue;
	if(fabs(tag->scEta) > 2.5) continue;

        if (!passElectronMVAIDIsoCombined(tag, electronIDMVAIDIsoCombined->MVAValue(
                                            tag->pt,tag->scEta,rho,
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
                                            tag->ip3dSig,
                                            tag->GsfTrackChi2OverNdof,
                                            tag->dEtaCalo,
                                            tag->dPhiCalo,
                                            tag->R9,
                                            tag->SCEtaWidth,
                                            tag->SCPhiWidth,
                                            tag->CovIEtaIPhi,
                                            tag->PreShowerOverRaw,
                                            tag->ChargedIso03,
                                            (tag->NeutralHadronIso03_05Threshold - tag->NeutralHadronIso007_05Threshold),
                                            (tag->GammaIso03_05Threshold - tag->GammaIsoVetoEtaStrip03_05Threshold),
                                            tag->ChargedIso04 ,
                                            (tag->NeutralHadronIso04_05Threshold - tag->NeutralHadronIso007_05Threshold),
                                            (tag->GammaIso04_05Threshold - tag->GammaIsoVetoEtaStrip04_05Threshold) 
                                            ), 
                                          0)) continue;
    
        if(!((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT)) &&
	   !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdL_CaloIsoVL)) &&
	   !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT)) &&
	   !((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) && (tag->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL)) ) 
	  continue;
    
        const Double_t m = 0.000511;
	TLorentzVector vtag;
	vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
	
	for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
	  if(i==j) continue;
	  
	  const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
	  if(probe->q == tag->q) continue;
	  
	  if(matchGen) {
	    Bool_t match1 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5);
	    Bool_t match2 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5);
	    if(!match1 && !match2)
	      continue;
	  }
	  
          //********************************************************************************
          // Probe Requirements
          //********************************************************************************
          if (!passElectronMVAIDIsoCombined(probe, electronIDMVAIDIsoCombined->MVAValue(
                                              probe->pt,probe->scEta,rho,
                                              probe->sigiEtaiEta, 
                                              probe->deltaEtaIn,
                                              probe->deltaPhiIn, 
                                              probe->HoverE,
                                              probe->d0,
                                              probe->dz, 
                                              probe->fBrem,
                                              probe->EOverP,
                                              probe->ESeedClusterOverPout,
                                              TMath::Sqrt(probe->sigiPhiiPhi),
                                              probe->nBrem,
                                              (1.0/(probe->scEt * TMath::CosH(probe->scEta)) - 1/probe->p), 
                                              probe->ESeedClusterOverPIn,
                                              probe->ip3d,
                                              probe->ip3dSig,
                                              probe->GsfTrackChi2OverNdof,
                                              probe->dEtaCalo,
                                              probe->dPhiCalo,
                                              probe->R9,
                                              probe->SCEtaWidth,
                                              probe->SCPhiWidth,
                                              probe->CovIEtaIPhi,
                                              probe->PreShowerOverRaw,
                                              probe->ChargedIso03,
                                              (probe->NeutralHadronIso03_05Threshold - probe->NeutralHadronIso007_05Threshold),
                                              (probe->GammaIso03_05Threshold - probe->GammaIsoVetoEtaStrip03_05Threshold),
                                              probe->ChargedIso04 ,
                                              (probe->NeutralHadronIso04_05Threshold - probe->NeutralHadronIso007_05Threshold),
                                              (probe->GammaIso04_05Threshold - probe->GammaIsoVetoEtaStrip04_05Threshold) 
                                              ), 
                                            0)) continue;
          
          
          if(!((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_SC8)) &&
             !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLTObject_Ele8)) &&
             !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLTObject_SC17)) &&
             !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLTObject_Ele17)) 
            )
	    continue;
 	  
	  TLorentzVector vprobe;
	  vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
	  
	  TLorentzVector vdielectron = vtag + vprobe;
	  if((vdielectron.M()<massLo) || (vdielectron.M()>massHi)) continue;
	  
	  nProbes += weight;
	  
          
          //********************************************************************************
          // Pass Trigger Definition
          //********************************************************************************
          ULong_t myTrigBit = 0;
          ULong_t myTrigObj = 0;
          if (info->runNum <= 170053) {
            myTrigObj = kHLTObject_Ele8_CaloIdL_CaloIsoVL;          
            myTrigBit = kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL;
          } else if (info->runNum <= 999999) {
            myTrigObj = kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;       
            myTrigBit = kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL;
          }
	  Bool_t pass = ( (probe->hltMatchBits & myTrigObj) && (info->triggerBits & myTrigBit));
	  
          // Fill tree
	  data.mass    = vdielectron.M();
	  data.pt      = probe->pt;
          data.eta     = probe->scEta;
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
  delete electronArr;
  
    
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
      
  gBenchmark->Show("selectDoubleElectronEffTP"); 
}
