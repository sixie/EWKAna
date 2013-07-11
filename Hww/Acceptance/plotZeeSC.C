#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TMath.h>                  // ROOT math utilities
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O

#include "MitCommon/MathTools/interface/MathUtils.h"

// define classes and constants to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TDielectron.hh"   
#endif


//=== MAIN MACRO =================================================================================================

void plotZeeSC(const TString fname = "/server/02a/ksung/EWKAna/p10-zee-powhep-cteq66-v26_M0_ntuple.root") 
{  
  gBenchmark->Start("plotZeeSC");  
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  // Generated Z boson mass window
  const Double_t massLow  = 60;
  const Double_t massHigh = 120;

  // ECAL gap
  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;
  
  // matching cone
  const Double_t matchR = 0.3;

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Set up histograms
  //
  Double_t nGenEvents = 0;  // number of events generated in mass window
  
  // overall
  Double_t nSelEvents = 0;
  Double_t acc;
  Double_t accErr;

  // barrel-barrel
  Double_t nSelEventsBB = 0;
  Double_t accBB;
  Double_t accErrBB;
  
  // barrel-endcap
  Double_t nSelEventsBE = 0;
  Double_t accBE;
  Double_t accErrBE;

  // endcap-endcap
  Double_t nSelEventsEE = 0;
  Double_t accEE;
  Double_t accErrEE;
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo   *gen     = new mithep::TGenInfo();
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");

  // Read input file
  cout << "Processing " << fname << "..." << endl;
  infile = new TFile(fname); 
  assert(infile);
        
  // Get the TTrees
  eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &gen);                 TBranch *genBr        = eventTree->GetBranch("Gen");
  eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");

  // loop through generator events
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      
    genBr->GetEntry(ientry);  
    if((gen->vmass < massLow) || (gen->vmass > massHigh)) continue;
      
    nGenEvents += gen->weight;

    infoBr->GetEntry(ientry);
    const UInt_t trigger = kHLT_Ele15_LW_L1R;
//    if(!(info->triggerBits & trigger)) continue;  // no trigger accept? Skip to next event...
	    
    // loop through dielectrons
    dielectronArr->Clear();
    dielectronBr->GetEntry(ientry);    
    for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
      const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
    
      // Exclude ECAL gap region (should already be done for ntuple, but just to make sure...)
      if((fabs(dielectron->scEta_1)>kGAP_LOW) && (fabs(dielectron->scEta_1)<kGAP_HIGH)) continue;
      if((fabs(dielectron->scEta_2)>kGAP_LOW) && (fabs(dielectron->scEta_2)<kGAP_HIGH)) continue;

      Bool_t isB1 = (fabs(dielectron->scEta_1)<kGAP_LOW);
      Bool_t isB2 = (fabs(dielectron->scEta_2)<kGAP_LOW);         
      
      // requirements on BOTH electrons
//      if((dielectron->mass < 60) || (dielectron->mass > 120))		                     continue;  // outside mass region? Skip to next event...
//      if(!(dielectron->hltMatchBits_1 & trigger) && !(dielectron->hltMatchBits_2 & trigger)) continue;  // no match to HLT object? Skip to next event... 
      if((dielectron->scEt_1 < 20)               || (dielectron->scEt_2 < 20))	             continue;  // below supercluster ET cut? Skip to next event...
      if((fabs(dielectron->scEta_1) > 2.5)       || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...

      // SC-Gen matching
      Double_t match11 = mithep::MathUtils::DeltaR(gen->phi_1, gen->eta_1, dielectron->scPhi_1, dielectron->scEta_1) < matchR;
      Double_t match12 = mithep::MathUtils::DeltaR(gen->phi_1, gen->eta_1, dielectron->scPhi_2, dielectron->scEta_2) < matchR;
      Double_t match21 = mithep::MathUtils::DeltaR(gen->phi_2, gen->eta_2, dielectron->scPhi_1, dielectron->scEta_1) < matchR;
      Double_t match22 = mithep::MathUtils::DeltaR(gen->phi_2, gen->eta_2, dielectron->scPhi_2, dielectron->scEta_2) < matchR;
      if((!match11 || !match22) && (!match12 || !match21)) continue;
      

      /******** We have a Z candidate! HURRAY! ********/
      
      nSelEvents += gen->weight;        
            
      if(isB1 && isB2)        { nSelEventsBB += gen->weight; } // barrel-barrel
      else if(!isB1 && !isB2) { nSelEventsEE += gen->weight; } // endcap-endcap
      else                    { nSelEventsBE += gen->weight; } // barrel-endcap            
    }
  }
  delete infile;
  infile=0, eventTree=0;
  delete info;
  delete gen;
  delete dielectronArr;        
  
  //
  // Compute acceptance for this sample
  //
  acc    = nSelEvents/nGenEvents;
  accErr = sqrt(acc*(1-acc)/nGenEvents);
    
  // barrel-barrel
  accBB    = nSelEventsBB/nGenEvents;
  accErrBB = sqrt(accBB*(1-accBB)/nGenEvents);
    
  // endcap-endcap
  accEE    = nSelEventsEE/nGenEvents;
  accErrEE = sqrt(accEE*(1-accEE)/nGenEvents);
    
  // barrel-endcap
  accBE    = nSelEventsBE/nGenEvents;
  accErrBE = sqrt(accBE*(1-accBE)/nGenEvents);     

 
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
      
  cout << "   Overall acceptance = " << setw(8) << nSelEvents   << " / " << setw(8) << nGenEvents << " = " << setw(10) << acc   << " +/- " << accErr   << endl;
  cout << "        BB acceptance = " << setw(8) << nSelEventsBB << " / " << setw(8) << nGenEvents << " = " << setw(10) << accBB << " +/- " << accErrBB << endl;
  cout << "        EE acceptance = " << setw(8) << nSelEventsEE << " / " << setw(8) << nGenEvents << " = " << setw(10) << accEE << " +/- " << accErrEE << endl;
  cout << "        BE acceptance = " << setw(8) << nSelEventsBE << " / " << setw(8) << nGenEvents << " = " << setw(10) << accBE << " +/- " << accErrBE << endl;
  cout << endl;

  gBenchmark->Show("plotZeeSC");
}
