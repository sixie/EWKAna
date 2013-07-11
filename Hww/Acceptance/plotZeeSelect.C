#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TMath.h>                  // ROOT math utilities
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/CPlot.hh"          // helper class for plots
#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/MyTools.hh"        // miscellaneous helper functions

// define classes and constants to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TDielectron.hh"   
#endif

//#define FILTER 40


//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);


//=== MAIN MACRO =================================================================================================

void plotZeeSelect(const TString input) 
{  
  gBenchmark->Start("plotZeeSelect");  
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Bool_t doSave  = true;    // save plots?
  TString format = "png";   // output file format
  
  TString outputDir;        // output directory
  vector<TString> fnamev;   // file names   
  vector<TString> labelv;   // legend label
  vector<Int_t>   colorv;   // color in plots
  vector<Int_t>   linev;    // line style

  Double_t massLow, massHigh;
  
  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    
    if(state==0) {
      outputDir = line;
      getline(ifs,line); stringstream ss1(line); ss1 >> massLow;
      getline(ifs,line); stringstream ss2(line); ss2 >> massHigh;
      state++;
    } else if(state==1) {
      string fname;
      Int_t color, linesty;
      stringstream ss(line);
      ss >> fname >> color >> linesty;
      string label = line.substr(line.find('@')+1);
      fnamev.push_back(fname);
      labelv.push_back(label);
      colorv.push_back(color);
      linev.push_back(linesty);
    }
  }
  ifs.close();
  
  CPlot::sOutDir = outputDir + TString("/select");

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  vector<Double_t> nGenEventsv;
  
  //
  // Set up histograms and counters
  //
  vector<TH1F*>    hPtv;
  vector<TH1F*>    hEtav;
  vector<TH1F*>    hPhiv;
  vector<TH1F*>    hDiEleMassv;
  vector<TH1F*>    hDiElePtv;
  vector<TH1F*>    hDiEleyv;
  vector<TH1F*>    hDiElePhiv;  
  vector<Double_t> nSelEventsv;
  vector<Double_t> accv;
  vector<Double_t> accErrv;

  vector<TH1F*>    hPtBBv;
  vector<TH1F*>    hEtaBBv;
  vector<TH1F*>    hPhiBBv;
  vector<TH1F*>    hDiEleMassBBv;
  vector<TH1F*>    hDiElePtBBv;
  vector<TH1F*>    hDiEleyBBv;
  vector<TH1F*>    hDiElePhiBBv;
  vector<Double_t> nSelEventsBBv;
  vector<Double_t> accBBv;
  vector<Double_t> accErrBBv;
  
  vector<TH1F*>    hPtBEv;
  vector<TH1F*>    hEtaBEv;
  vector<TH1F*>    hPhiBEv;
  vector<TH1F*>    hDiEleMassBEv;
  vector<TH1F*>    hDiElePtBEv;
  vector<TH1F*>    hDiEleyBEv;
  vector<TH1F*>    hDiElePhiBEv;
  vector<Double_t> nSelEventsBEv;
  vector<Double_t> accBEv;
  vector<Double_t> accErrBEv;

  vector<TH1F*>    hPtEEv;
  vector<TH1F*>    hEtaEEv;
  vector<TH1F*>    hPhiEEv;  
  vector<TH1F*>    hDiEleMassEEv;
  vector<TH1F*>    hDiElePtEEv;
  vector<TH1F*>    hDiEleyEEv;
  vector<TH1F*>    hDiElePhiEEv;
  vector<Double_t> nSelEventsEEv;
  vector<Double_t> accEEv;
  vector<Double_t> accErrEEv;
  
  char hname[100];
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    nGenEventsv.push_back(0);
    
    // overall
    sprintf(hname,"hPt_%i",ifile);        hPtv.push_back(new TH1F(hname,"",150,0,150));         hPtv[ifile]->Sumw2();
    sprintf(hname,"hEta_%i",ifile);       hEtav.push_back(new TH1F(hname,"",150,-3,3));         hEtav[ifile]->Sumw2();
    sprintf(hname,"hPhi_%i",ifile);       hPhiv.push_back(new TH1F(hname,"",40,-3.2,3.2));      hPhiv[ifile]->Sumw2();    
    sprintf(hname,"hDiEleMass_%i",ifile); hDiEleMassv.push_back(new TH1F(hname,"",180,20,200)); hDiEleMassv[ifile]->Sumw2();
    sprintf(hname,"hDiElePt_%i",ifile);   hDiElePtv.push_back(new TH1F(hname,"",100,0,100));    hDiElePtv[ifile]->Sumw2();
    sprintf(hname,"hDiEley_%i",ifile);    hDiEleyv.push_back(new TH1F(hname,"",50,-5,5));       hDiEleyv[ifile]->Sumw2();
    sprintf(hname,"hDiElePhi_%i",ifile);  hDiElePhiv.push_back(new TH1F(hname,"",50,-3.2,3.2)); hDiElePhiv[ifile]->Sumw2();        
    nSelEventsv.push_back(0); 
    
    // barrel-barrel
    sprintf(hname,"hPtBB_%i",ifile);        hPtBBv.push_back(new TH1F(hname,"",150,0,150));         hPtBBv[ifile]->Sumw2();
    sprintf(hname,"hEtaBB_%i",ifile);       hEtaBBv.push_back(new TH1F(hname,"",150,-3,3));         hEtaBBv[ifile]->Sumw2();
    sprintf(hname,"hPhiBB_%i",ifile);       hPhiBBv.push_back(new TH1F(hname,"",40,-3.2,3.2));      hPhiBBv[ifile]->Sumw2();    
    sprintf(hname,"hDiEleMassBB_%i",ifile); hDiEleMassBBv.push_back(new TH1F(hname,"",180,20,200)); hDiEleMassBBv[ifile]->Sumw2();
    sprintf(hname,"hDiElePtBB_%i",ifile);   hDiElePtBBv.push_back(new TH1F(hname,"",100,0,100));    hDiElePtBBv[ifile]->Sumw2();
    sprintf(hname,"hDiEleyBB_%i",ifile);    hDiEleyBBv.push_back(new TH1F(hname,"",50,-5,5));       hDiEleyBBv[ifile]->Sumw2();
    sprintf(hname,"hDiElePhiBB_%i",ifile);  hDiElePhiBBv.push_back(new TH1F(hname,"",50,-3.2,3.2)); hDiElePhiBBv[ifile]->Sumw2();    
    nSelEventsBBv.push_back(0); 

    // barrel-endcap
    sprintf(hname,"hPtBE_%i",ifile);        hPtBEv.push_back(new TH1F(hname,"",150,0,150));         hPtBEv[ifile]->Sumw2();
    sprintf(hname,"hEtaBE_%i",ifile);       hEtaBEv.push_back(new TH1F(hname,"",150,-3,3));         hEtaBEv[ifile]->Sumw2();
    sprintf(hname,"hPhiBE_%i",ifile);       hPhiBEv.push_back(new TH1F(hname,"",40,-3.2,3.2));      hPhiBEv[ifile]->Sumw2();    
    sprintf(hname,"hDiEleMassBE_%i",ifile); hDiEleMassBEv.push_back(new TH1F(hname,"",180,20,200)); hDiEleMassBEv[ifile]->Sumw2();
    sprintf(hname,"hDiElePtBE_%i",ifile);   hDiElePtBEv.push_back(new TH1F(hname,"",100,0,100));    hDiElePtBEv[ifile]->Sumw2();
    sprintf(hname,"hDiEleyBE_%i",ifile);    hDiEleyBEv.push_back(new TH1F(hname,"",50,-5,5));       hDiEleyBEv[ifile]->Sumw2();
    sprintf(hname,"hDiElePhiBE_%i",ifile);  hDiElePhiBEv.push_back(new TH1F(hname,"",50,-3.2,3.2)); hDiElePhiBEv[ifile]->Sumw2();    
    nSelEventsBEv.push_back(0); 
    
    // endcap-endcap
    sprintf(hname,"hPtEE_%i",ifile);        hPtEEv.push_back(new TH1F(hname,"",150,0,150));         hPtEEv[ifile]->Sumw2();
    sprintf(hname,"hEtaEE_%i",ifile);       hEtaEEv.push_back(new TH1F(hname,"",150,-3,3));         hEtaEEv[ifile]->Sumw2();
    sprintf(hname,"hPhiEE_%i",ifile);       hPhiEEv.push_back(new TH1F(hname,"",40,-3.2,3.2));      hPhiEEv[ifile]->Sumw2();    
    sprintf(hname,"hDiEleMassEE_%i",ifile); hDiEleMassEEv.push_back(new TH1F(hname,"",180,20,200)); hDiEleMassEEv[ifile]->Sumw2();
    sprintf(hname,"hDiElePtEE_%i",ifile);   hDiElePtEEv.push_back(new TH1F(hname,"",100,0,100));    hDiElePtEEv[ifile]->Sumw2();
    sprintf(hname,"hDiEleyEE_%i",ifile);    hDiEleyEEv.push_back(new TH1F(hname,"",50,-5,5));       hDiEleyEEv[ifile]->Sumw2();
    sprintf(hname,"hDiElePhiEE_%i",ifile);  hDiElePhiEEv.push_back(new TH1F(hname,"",50,-3.2,3.2)); hDiElePhiEEv[ifile]->Sumw2();    
    nSelEventsEEv.push_back(0);    
  }
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  
  // loop over samples
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
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
#ifdef FILTER
      if(gen->mass < FILTER) continue;
#endif   
      if((gen->vmass < massLow) || (gen->vmass > massHigh)) continue;
      
      nGenEventsv[ifile] += gen->weight;

      infoBr->GetEntry(ientry);
      const UInt_t trigger = kHLT_Ele15_LW_L1R;
      if(!(info->triggerBits & trigger)) continue;  // no trigger accept? Skip to next event...
	    
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
     
        // outside mass region? Skip to next event...
        if((dielectron->mass < 60) || (dielectron->mass > 120))	continue;
    
//********** WP95 **********	  
	// requirements on BOTH electrons
        if(!(dielectron->hltMatchBits_1 & trigger) && !(dielectron->hltMatchBits_2 & trigger)) continue;  // no match to HLT object? Skip to next event... 
        if((dielectron->scEt_1 < 20)		   || (dielectron->scEt_2 < 20))	       continue;  // below supercluster ET cut? Skip to next event...
        if((fabs(dielectron->scEta_1) > 2.5)	   || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
        if((dielectron->nExpHitsInner_1 > 1)	   || (dielectron->nExpHitsInner_2 > 1))       continue;  // more hits expected inside first hit? Skip to next event...
        if(!(dielectron->isEcalDriven_1)	   || !(dielectron->isEcalDriven_2))	       continue;  // not ECAL seeded electrons? Skip to next event...
        
	// barrel/endcap dependent requirments      
        if(isB1) {  // barrel
//          Double_t iso = (dielectron->trkIso03_1 + TMath::Max(dielectron->emIso03_1-1,Float_t(0)) + dielectron->hadIso03_1)/dielectron->pt_1;
//          if(iso>0.15) continue;
          if(dielectron->trkIso03_1	    > 0.15*(dielectron->pt_1)) continue;
	  if(dielectron->emIso03_1	    > 2.00*(dielectron->pt_1)) continue;
          if(dielectron->hadIso03_1	    > 0.12*(dielectron->pt_1)) continue;
	  if(dielectron->sigiEtaiEta_1      > 0.01)		       continue;
//	    if(fabs(dielectron->deltaPhiIn_1) > 0.8)			 continue;
	  if(fabs(dielectron->deltaEtaIn_1) > 0.007)		       continue;
	  if(dielectron->HoverE_1	    > 0.15)		       continue;
        } else {  // endcap
//          Double_t iso = (dielectron->trkIso03_1 + dielectron->emIso03_1 + dielectron->hadIso03_1)/dielectron->pt_1;
//          if(iso>0.1) continue;
          if(dielectron->trkIso03_1	    > 0.08*(dielectron->pt_1)) continue;
	  if(dielectron->emIso03_1	    > 0.06*(dielectron->pt_1)) continue;
	  if(dielectron->hadIso03_1	    > 0.05*(dielectron->pt_1)) continue;
	  if(dielectron->sigiEtaiEta_1      > 0.03)		       continue;
//	    if(fabs(dielectron->deltaPhiIn_1) > 0.7)			 continue;
//	    if(fabs(dielectron->deltaEtaIn_1) > 0.01)			 continue;
	  if(dielectron->HoverE_1	    > 0.07)		       continue;
        }
      
        if(isB2) {  // barrel
//          Double_t iso = (dielectron->trkIso03_2 + TMath::Max(dielectron->emIso03_2-1,Float_t(0)) + dielectron->hadIso03_2)/dielectron->pt_2;
//          if(iso>0.15) continue;
          if(dielectron->trkIso03_2	    > 0.15*(dielectron->pt_2)) continue;
	  if(dielectron->emIso03_2	    > 2.00*(dielectron->pt_2)) continue;
	  if(dielectron->hadIso03_2	    > 0.12*(dielectron->pt_2)) continue;
	  if(dielectron->sigiEtaiEta_2      > 0.01)		       continue;
//	    if(fabs(dielectron->deltaPhiIn_2) > 0.8)			 continue;
	  if(fabs(dielectron->deltaEtaIn_2) > 0.007)		       continue;
	  if(dielectron->HoverE_2	    > 0.15)		       continue;
        } else {  // endcap
//          Double_t iso = (dielectron->trkIso03_2 + dielectron->emIso03_2 + dielectron->hadIso03_2)/dielectron->pt_2;
//          if(iso>0.1) continue;
          if(dielectron->trkIso03_2	    > 0.08*(dielectron->pt_2)) continue;
	  if(dielectron->emIso03_2	    > 0.06*(dielectron->pt_2)) continue;
	  if(dielectron->hadIso03_2	    > 0.05*(dielectron->pt_2)) continue;
	  if(dielectron->sigiEtaiEta_2      > 0.03)		       continue;
//	    if(fabs(dielectron->deltaPhiIn_2) > 0.7)			 continue;
//	    if(fabs(dielectron->deltaEtaIn_2) > 0.01)			 continue;
	  if(dielectron->HoverE_2	    > 0.07)		       continue;
        }	   
//*/
/******** WP80 **********	
	// requirements on BOTH electrons
        if(!(dielectron->hltMatchBits_1 & trigger) && !(dielectron->hltMatchBits_2 & trigger)) continue;  // no match to HLT object? Skip to next event... 
        if((dielectron->scEt_1 < 20)		   || (dielectron->scEt_2 < 20))	       continue;  // below supercluster ET cut? Skip to next event...
	if((fabs(dielectron->scEta_1) > 2.5)	   || (fabs(dielectron->scEta_2) > 2.5))       continue;  // outside eta range? Skip to next event...
        if(!(dielectron->isEcalDriven_1)	   || !(dielectron->isEcalDriven_2))	       continue;  // not ECAL seeded electrons? Skip to next event...
	
	// conversion rejection
	if((dielectron->nExpHitsInner_1 > 0)	   || (dielectron->nExpHitsInner_2 > 0))	 continue;
        if((fabs(dielectron->partnerDist_1)<0.02) && (fabs(dielectron->partnerDeltaCot_1)<0.02)) continue;
	if((fabs(dielectron->partnerDist_2)<0.02) && (fabs(dielectron->partnerDeltaCot_2)<0.02)) continue;	
		
	// barrel/endcap dependent requirments      
        if(isB1) {  // barrel
          if(dielectron->trkIso03_1	    > 0.09*(dielectron->pt_1)) continue;
	  if(dielectron->emIso03_1	    > 0.07*(dielectron->pt_1)) continue;
          if(dielectron->hadIso03_1	    > 0.10*(dielectron->pt_1)) continue;
	  if(dielectron->sigiEtaiEta_1      > 0.01)		       continue;
          if(fabs(dielectron->deltaPhiIn_1) > 0.06)		       continue;
	  if(fabs(dielectron->deltaEtaIn_1) > 0.004)		       continue;
	  if(dielectron->HoverE_1	    > 0.04)		       continue;
        } else {  // endcap
          if(dielectron->trkIso03_1	    > 0.04*(dielectron->pt_1))  continue;
	  if(dielectron->emIso03_1	    > 0.05*(dielectron->pt_1))  continue;
	  if(dielectron->hadIso03_1	    > 0.025*(dielectron->pt_1)) continue;
	  if(dielectron->sigiEtaiEta_1      > 0.03)			continue;
	  if(fabs(dielectron->deltaPhiIn_1) > 0.03)			continue;
	  if(dielectron->HoverE_1	    > 0.025)			continue;
        }
      
        if(isB2) {  // barrel
          if(dielectron->trkIso03_2	    > 0.09*(dielectron->pt_2)) continue;
	  if(dielectron->emIso03_2	    > 0.07*(dielectron->pt_2)) continue;
	  if(dielectron->hadIso03_2	    > 0.10*(dielectron->pt_2)) continue;
	  if(dielectron->sigiEtaiEta_2      > 0.01)		       continue;
	  if(fabs(dielectron->deltaPhiIn_2) > 0.06)		       continue;
	  if(fabs(dielectron->deltaEtaIn_2) > 0.004)		       continue;
	  if(dielectron->HoverE_2	    > 0.04)		       continue;
        } else {  // endcap
          if(dielectron->trkIso03_2	    > 0.04*(dielectron->pt_2))  continue;
	  if(dielectron->emIso03_2	    > 0.05*(dielectron->pt_2))  continue;
	  if(dielectron->hadIso03_2	    > 0.025*(dielectron->pt_2)) continue;
	  if(dielectron->sigiEtaiEta_2      > 0.03)			continue;
	  if(fabs(dielectron->deltaPhiIn_2) > 0.03)			continue;
	  if(dielectron->HoverE_2	    > 0.025)			continue;
        }	 
//*/


        /******** We have a Z candidate! HURRAY! ********/
      
	nSelEventsv[ifile] += gen->weight;      
      
        hDiEleMassv[ifile]->Fill(dielectron->mass,gen->weight);
        hDiElePtv[ifile]  ->Fill(dielectron->pt,  gen->weight);
        hDiEleyv[ifile]   ->Fill(dielectron->y,   gen->weight);
        hDiElePhiv[ifile] ->Fill(dielectron->phi, gen->weight);      
      
        hPtv[ifile] ->Fill(dielectron->pt_1, gen->weight); hPtv[ifile]->Fill(dielectron->pt_2,  gen->weight);
        hEtav[ifile]->Fill(dielectron->eta_1,gen->weight); hEtav[ifile]->Fill(dielectron->eta_2,gen->weight);
        hPhiv[ifile]->Fill(dielectron->phi_1,gen->weight); hPhiv[ifile]->Fill(dielectron->phi_2,gen->weight);      
            
        if(isB1 && isB2) {  // barrel-barrel        
 	  nSelEventsBBv[ifile] += gen->weight;      
          hDiEleMassBBv[ifile]->Fill(dielectron->mass,gen->weight);
          hDiElePtBBv[ifile]  ->Fill(dielectron->pt,  gen->weight);
          hDiEleyBBv[ifile]   ->Fill(dielectron->y,   gen->weight);
          hDiElePhiBBv[ifile] ->Fill(dielectron->phi, gen->weight);      
        
	  hPtBBv[ifile] ->Fill(dielectron->pt_1, gen->weight); hPtBBv[ifile] ->Fill(dielectron->pt_2, gen->weight);
          hEtaBBv[ifile]->Fill(dielectron->eta_1,gen->weight); hEtaBBv[ifile]->Fill(dielectron->eta_2,gen->weight);
          hPhiBBv[ifile]->Fill(dielectron->phi_1,gen->weight); hPhiBBv[ifile]->Fill(dielectron->phi_2,gen->weight);	
      
        } else if(!isB1 && !isB2) {  // endcap-endcap
          nSelEventsEEv[ifile] += gen->weight;      
          hDiEleMassEEv[ifile]->Fill(dielectron->mass,gen->weight);
          hDiElePtEEv[ifile]  ->Fill(dielectron->pt,  gen->weight);
          hDiEleyEEv[ifile]   ->Fill(dielectron->y,   gen->weight);
          hDiElePhiEEv[ifile] ->Fill(dielectron->phi, gen->weight);      
        
	  hPtEEv[ifile] ->Fill(dielectron->pt_1, gen->weight); hPtEEv[ifile] ->Fill(dielectron->pt_2, gen->weight);
          hEtaEEv[ifile]->Fill(dielectron->eta_1,gen->weight); hEtaEEv[ifile]->Fill(dielectron->eta_2,gen->weight);
          hPhiEEv[ifile]->Fill(dielectron->phi_1,gen->weight); hPhiEEv[ifile]->Fill(dielectron->phi_2,gen->weight);
      
        } else {  // barrel-endcap
          nSelEventsBEv[ifile] += gen->weight;      
          hDiEleMassBEv[ifile]->Fill(dielectron->mass,gen->weight);
          hDiElePtBEv[ifile]  ->Fill(dielectron->pt,  gen->weight);
          hDiEleyBEv[ifile]   ->Fill(dielectron->y,   gen->weight);
          hDiElePhiBEv[ifile] ->Fill(dielectron->phi, gen->weight);      
        
	  hPtBEv[ifile] ->Fill(dielectron->pt_1, gen->weight); hPtBEv[ifile] ->Fill(dielectron->pt_2, gen->weight);
          hEtaBEv[ifile]->Fill(dielectron->eta_1,gen->weight); hEtaBEv[ifile]->Fill(dielectron->eta_2,gen->weight);
          hPhiBEv[ifile]->Fill(dielectron->phi_1,gen->weight); hPhiBEv[ifile]->Fill(dielectron->phi_2,gen->weight);
        }            
      }
    }
    
    delete infile;
    infile=0, eventTree=0;
        
    //
    // Compute acceptance for this sample
    //
    accv.push_back(nSelEventsv[ifile]/nGenEventsv[ifile]);
    accErrv.push_back(sqrt(accv[ifile]*(1-accv[ifile])/nGenEventsv[ifile]));
    
    // barrel-barrel
    accBBv.push_back(nSelEventsBBv[ifile]/nGenEventsv[ifile]);
    accErrBBv.push_back(sqrt(accBBv[ifile]*(1-accBBv[ifile])/nGenEventsv[ifile]));
    
    // endcap-endcap
    accEEv.push_back(nSelEventsEEv[ifile]/nGenEventsv[ifile]);
    accErrEEv.push_back(sqrt(accEEv[ifile]*(1-accEEv[ifile])/nGenEventsv[ifile]));
    
    // barrel-endcap
    accBEv.push_back(nSelEventsBEv[ifile]/nGenEventsv[ifile]);
    accErrBEv.push_back(sqrt(accBEv[ifile]*(1-accBEv[ifile])/nGenEventsv[ifile]));     
  }
  delete info;
  delete gen;
  delete dielectronArr;

  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  // *** all histograms normalized to 1 unless stated otherwise
  //==============================================================================================================
  TCanvas *c = MakeCanvas("c","c",800,600);
  
  // string buffers
  char ylabel[100];   // y-axis label
  
  // electron pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","p_{T}(e) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPtv[i]->Scale(1.0/hPtv[i]->Integral());
    plotPt.AddHist1D(hPtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPt.Draw(c,doSave,format);
  
  // electron eta
  sprintf(ylabel,"a.u. / %.2f",hEtav[0]->GetBinWidth(1)); 
  CPlot plotEta("eta","","#eta(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hEtav[i]->Scale(1.0/hEtav[i]->Integral());
    plotEta.AddHist1D(hEtav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotEta.TransLegend(-0.15,-0.5);
  plotEta.Draw(c,doSave,format);
  
  // electron phi
  sprintf(ylabel,"a.u. / %.1f",hPhiv[0]->GetBinWidth(1));
  CPlot plotPhi("phi","","#phi(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhiv[i]->Scale(1.0/hPhiv[i]->Integral());
    plotPhi.AddHist1D(hPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhi.SetYRange(0.8*(hPhiv[0]->GetMaximum()),1.2*(hPhiv[0]->GetMaximum()));
  plotPhi.Draw(c,doSave,format);

  // dielectron mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDiEleMassv[0]->GetBinWidth(1));
  CPlot plotDiEleMass("dielectronmass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleMassv[i]->Scale(1.0/hDiEleMassv[i]->Integral());
    plotDiEleMass.AddHist1D(hDiEleMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleMass.SetLogy();
  plotDiEleMass.Draw(c,doSave,format);
      
  // dielectron pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDiElePtv[0]->GetBinWidth(1));
  CPlot plotDiElePt("dielectronpt","","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePtv[i]->Scale(1.0/hDiElePtv[i]->Integral());
    plotDiElePt.AddHist1D(hDiElePtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePt.SetLogy();
  plotDiElePt.Draw(c,doSave,format);
  
  // dielectron rapidity
  sprintf(ylabel,"a.u. / %.2f",hDiEleyv[0]->GetBinWidth(1));
  CPlot plotDiEley("dielectrony","","y(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleyv[i]->Scale(1.0/hDiEleyv[i]->Integral());
    plotDiEley.AddHist1D(hDiEleyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEley.TransLegend(0.1,0);
  plotDiEley.Draw(c,doSave,format);
  
  // dielectron phi
  sprintf(ylabel,"a.u. / %.2f",hDiElePhiv[0]->GetBinWidth(1));
  CPlot plotDiElePhi("dielectronphi","","#phi(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePhiv[i]->Scale(1.0/hDiElePhiv[i]->Integral());
    plotDiElePhi.AddHist1D(hDiElePhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePhi.SetYRange(0.8*(hDiElePhiv[0]->GetMaximum()),1.2*(hDiElePhiv[0]->GetMaximum()));
  plotDiElePhi.Draw(c,doSave,format);

  // electron pT (barrel-barrel)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPtBBv[0]->GetBinWidth(1));
  CPlot plotPtBB("ptBB","barrel-barrel","p_{T}(e) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPtBBv[i]->Scale(1.0/hPtBBv[i]->Integral());
    plotPtBB.AddHist1D(hPtBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPtBB.Draw(c,doSave,format);
  
  // electron eta (barrel-barrel)
  sprintf(ylabel,"a.u. / %.2f",hEtaBBv[0]->GetBinWidth(1)); 
  CPlot plotEtaBB("etaBB","barrel-barrel","#eta(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hEtaBBv[i]->Scale(1.0/hEtaBBv[i]->Integral());
    plotEtaBB.AddHist1D(hEtaBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotEtaBB.TransLegend(-0.15,-0.5);
  plotEtaBB.Draw(c,doSave,format);
  
  // electron phi (barrel-barrel)
  sprintf(ylabel,"a.u. / %.1f",hPhiBBv[0]->GetBinWidth(1));
  CPlot plotPhiBB("phiBB","barrel-barrel","#phi(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhiBBv[i]->Scale(1.0/hPhiBBv[i]->Integral());
    plotPhiBB.AddHist1D(hPhiBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhiBB.SetYRange(0.8*(hPhiBBv[0]->GetMaximum()),1.2*(hPhiBBv[0]->GetMaximum()));
  plotPhiBB.Draw(c,doSave,format);

  // dielectron mass (barrel-barrel)
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDiEleMassBBv[0]->GetBinWidth(1));
  CPlot plotDiEleMassBB("dielectronmassBB","barrel-barrel","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleMassBBv[i]->Scale(1.0/hDiEleMassBBv[i]->Integral());
    plotDiEleMassBB.AddHist1D(hDiEleMassBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleMassBB.SetLogy();
  plotDiEleMassBB.Draw(c,doSave,format);
      
  // dielectron pT (barrel-barrel)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDiElePtBBv[0]->GetBinWidth(1));
  CPlot plotDiElePtBB("dielectronptBB","barrel-barrel","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePtBBv[i]->Scale(1.0/hDiElePtBBv[i]->Integral());
    plotDiElePtBB.AddHist1D(hDiElePtBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePtBB.SetLogy();
  plotDiElePtBB.Draw(c,doSave,format);
  
  // dielectron rapidity (barrel-barrel)
  sprintf(ylabel,"a.u. / %.2f",hDiEleyBBv[0]->GetBinWidth(1));
  CPlot plotDiEleyBB("dielectronyBB","barrel-barrel","y(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleyBBv[i]->Scale(1.0/hDiEleyBBv[i]->Integral());
    plotDiEleyBB.AddHist1D(hDiEleyBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleyBB.TransLegend(0.1,0);
  plotDiEleyBB.Draw(c,doSave,format);
  
  // dielectron phi (barrel-barrel)
  sprintf(ylabel,"a.u. / %.2f",hDiElePhiBBv[0]->GetBinWidth(1));
  CPlot plotDiElePhiBB("dielectronphiBB","barrel-barrel","#phi(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePhiBBv[i]->Scale(1.0/hDiElePhiBBv[i]->Integral());
    plotDiElePhiBB.AddHist1D(hDiElePhiBBv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePhiBB.SetYRange(0.8*(hDiElePhiBBv[0]->GetMaximum()),1.2*(hDiElePhiBBv[0]->GetMaximum()));
  plotDiElePhiBB.Draw(c,doSave,format);
  
  // electron pT (endcap-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPtEEv[0]->GetBinWidth(1));
  CPlot plotPtEE("ptEE","endcap-endcap","p_{T}(e) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPtEEv[i]->Scale(1.0/hPtEEv[i]->Integral());
    plotPtEE.AddHist1D(hPtEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPtEE.Draw(c,doSave,format);
  
  // electron eta (endcap-endcap)
  sprintf(ylabel,"a.u. / %.2f",hEtaEEv[0]->GetBinWidth(1)); 
  CPlot plotEtaEE("etaEE","endcap-endcap","#eta(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hEtaEEv[i]->Scale(1.0/hEtaEEv[i]->Integral());
    plotEtaEE.AddHist1D(hEtaEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotEtaEE.TransLegend(-0.15,-0.5);
  plotEtaEE.Draw(c,doSave,format);
  
  // electron phi (endcap-endcap)
  sprintf(ylabel,"a.u. / %.1f",hPhiEEv[0]->GetBinWidth(1));
  CPlot plotPhiEE("phiEE","endcap-endcap","#phi(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhiEEv[i]->Scale(1.0/hPhiEEv[i]->Integral());
    plotPhiEE.AddHist1D(hPhiEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhiEE.SetYRange(0.8*(hPhiEEv[0]->GetMaximum()),1.2*(hPhiEEv[0]->GetMaximum()));
  plotPhiEE.Draw(c,doSave,format);

  // dielectron mass (endcap-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDiEleMassEEv[0]->GetBinWidth(1));
  CPlot plotDiEleMassEE("dielectronmassEE","endcap-endcap","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleMassEEv[i]->Scale(1.0/hDiEleMassEEv[i]->Integral());
    plotDiEleMassEE.AddHist1D(hDiEleMassEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleMassEE.SetLogy();
  plotDiEleMassEE.Draw(c,doSave,format);
      
  // dielectron pT (endcap-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDiElePtEEv[0]->GetBinWidth(1));
  CPlot plotDiElePtEE("dielectronptEE","endcap-endcap","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePtEEv[i]->Scale(1.0/hDiElePtEEv[i]->Integral());
    plotDiElePtEE.AddHist1D(hDiElePtEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePtEE.SetLogy();
  plotDiElePtEE.Draw(c,doSave,format);
  
  // dielectron rapidity (endcap-endcap)
  sprintf(ylabel,"a.u. / %.2f",hDiEleyEEv[0]->GetBinWidth(1));
  CPlot plotDiEleyEE("dielectronyEE","endcap-endcap","y(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleyEEv[i]->Scale(1.0/hDiEleyEEv[i]->Integral());
    plotDiEleyEE.AddHist1D(hDiEleyEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleyEE.TransLegend(0.1,0);
  plotDiEleyEE.Draw(c,doSave,format);
  
  // dielectron phi (endcap-endcap)
  sprintf(ylabel,"a.u. / %.2f",hDiElePhiEEv[0]->GetBinWidth(1));
  CPlot plotDiElePhiEE("dielectronphiEE","endcap-endcap","#phi(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePhiEEv[i]->Scale(1.0/hDiElePhiEEv[i]->Integral());
    plotDiElePhiEE.AddHist1D(hDiElePhiEEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePhiEE.SetYRange(0.8*(hDiElePhiEEv[0]->GetMaximum()),1.2*(hDiElePhiEEv[0]->GetMaximum()));
  plotDiElePhiEE.Draw(c,doSave,format);

  // electron pT (barrel-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPtBEv[0]->GetBinWidth(1));
  CPlot plotPtBE("ptBE","barrel-endcap","p_{T}(e) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPtBEv[i]->Scale(1.0/hPtBEv[i]->Integral());
    plotPtBE.AddHist1D(hPtBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPtBE.Draw(c,doSave,format);
  
  // electron eta (barrel-endcap)
  sprintf(ylabel,"a.u. / %.2f",hEtaBEv[0]->GetBinWidth(1)); 
  CPlot plotEtaBE("etaBE","barrel-endcap","#eta(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hEtaBEv[i]->Scale(1.0/hEtaBEv[i]->Integral());
    plotEtaBE.AddHist1D(hEtaBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotEtaBE.TransLegend(-0.15,-0.5);
  plotEtaBE.Draw(c,doSave,format);
  
  // electron phi (barrel-endcap)
  sprintf(ylabel,"a.u. / %.1f",hPhiBEv[0]->GetBinWidth(1));
  CPlot plotPhiBE("phiBE","barrel-endcap","#phi(e)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhiBEv[i]->Scale(1.0/hPhiBEv[i]->Integral());
    plotPhiBE.AddHist1D(hPhiBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhiBE.SetYRange(0.8*(hPhiBEv[0]->GetMaximum()),1.2*(hPhiBEv[0]->GetMaximum()));
  plotPhiBE.Draw(c,doSave,format);

  // dielectron mass (barrel-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDiEleMassBEv[0]->GetBinWidth(1));
  CPlot plotDiEleMassBE("dielectronmassBE","barrel-endcap","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleMassBEv[i]->Scale(1.0/hDiEleMassBEv[i]->Integral());
    plotDiEleMassBE.AddHist1D(hDiEleMassBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleMassBE.SetLogy();
  plotDiEleMassBE.Draw(c,doSave,format);
      
  // dielectron pT (barrel-endcap)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDiElePtBEv[0]->GetBinWidth(1));
  CPlot plotDiElePtBE("dielectronptBE","barrel-endcap","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePtBEv[i]->Scale(1.0/hDiElePtBEv[i]->Integral());
    plotDiElePtBE.AddHist1D(hDiElePtBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePtBE.SetLogy();
  plotDiElePtBE.Draw(c,doSave,format);
  
  // dielectron rapidity (barrel-endcap)
  sprintf(ylabel,"a.u. / %.2f",hDiEleyBEv[0]->GetBinWidth(1));
  CPlot plotDiEleyBE("dielectronyBE","barrel-endcap","y(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiEleyBEv[i]->Scale(1.0/hDiEleyBEv[i]->Integral());
    plotDiEleyBE.AddHist1D(hDiEleyBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiEleyBE.TransLegend(0.1,0);
  plotDiEleyBE.Draw(c,doSave,format);
  
  // dielectron phi (barrel-endcap)
  sprintf(ylabel,"a.u. / %.2f",hDiElePhiBEv[0]->GetBinWidth(1));
  CPlot plotDiElePhiBE("dielectronphiBE","barrel-endcap","#phi(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDiElePhiBEv[i]->Scale(1.0/hDiElePhiBEv[i]->Integral());
    plotDiElePhiBE.AddHist1D(hDiElePhiBEv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDiElePhiBE.SetYRange(0.8*(hDiElePhiBEv[0]->GetMaximum()),1.2*(hDiElePhiBEv[0]->GetMaximum()));
  plotDiElePhiBE.Draw(c,doSave,format);

  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  for(UInt_t i=0; i<fnamev.size(); i++) {    
    cout << labelv[i] << " file: " << fnamev[i] << endl;
    cout << "   Overall acceptance = " << setw(8) << nSelEventsv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accv[i] << " +/- " << accErrv[i];
    if(i>0) {
      Double_t da    = 100.*(accv[i]-accv[0])/accv[0];
      Double_t daErr = 100.*sqrt((accErrv[i])*(accErrv[i])+(accErrv[0])*(accErrv[0]))/accv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
    } else {
      cout << endl;
    }
    cout << "        BB acceptance = " << setw(8) << nSelEventsBBv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accBBv[i] << " +/- " << accErrBBv[i];
    if(i>0) {
      Double_t da    = 100.*(accBBv[i]-accBBv[0])/accBBv[0];
      Double_t daErr = 100.*sqrt((accErrBBv[i])*(accErrBBv[i])+(accErrBBv[0])*(accErrBBv[0]))/accBBv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
    } else {
      cout << endl;
    }
    cout << "        EE acceptance = " << setw(8) << nSelEventsEEv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accEEv[i] << " +/- " << accErrEEv[i];
    if(i>0) {
      Double_t da    = 100.*(accEEv[i]-accEEv[0])/accEEv[0];
      Double_t daErr = 100.*sqrt((accErrEEv[i])*(accErrEEv[i])+(accErrEEv[0])*(accErrEEv[0]))/accEEv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
    } else {
      cout << endl;
    }
    cout << "        BE acceptance = " << setw(8) << nSelEventsBEv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accBEv[i] << " +/- " << accErrBEv[i];
    if(i>0) {
      Double_t da    = 100.*(accBEv[i]-accBEv[0])/accBEv[0];
      Double_t daErr = 100.*sqrt((accErrBEv[i])*(accErrBEv[i])+(accErrBEv[0])*(accErrBEv[0]))/accBEv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
    } else {
      cout << endl;
    }
    cout << endl;
  }
  cout << endl;
   
  if(doSave) {
    makeHTML(outputDir);
    
    ofstream txtfile;
    char txtfname[100];    
    sprintf(txtfname,"%s/select.txt",outputDir.Data());
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    txtfile << endl;
  
    for(UInt_t i=0; i<fnamev.size(); i++) {    
      txtfile << labelv[i] << " file: " << fnamev[i] << endl;
      txtfile << "   Overall acceptance = " << setw(8) << nSelEventsv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accv[i] << " +/- " << accErrv[i];
      if(i>0) {
        Double_t da    = 100.*(accv[i]-accv[0])/accv[0];
        Double_t daErr = 100.*sqrt((accErrv[i])*(accErrv[i])+(accErrv[0])*(accErrv[0]))/accv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
      } else {
        txtfile << endl;
      }
      txtfile << "        BB acceptance = " << setw(8) << nSelEventsBBv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accBBv[i] << " +/- " << accErrBBv[i];
      if(i>0) {
        Double_t da    = 100.*(accBBv[i]-accBBv[0])/accBBv[0];
        Double_t daErr = 100.*sqrt((accErrBBv[i])*(accErrBBv[i])+(accErrBBv[0])*(accErrBBv[0]))/accBBv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
      } else {
        txtfile << endl;
      }
      txtfile << "        EE acceptance = " << setw(8) << nSelEventsEEv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accEEv[i] << " +/- " << accErrEEv[i];
      if(i>0) {
        Double_t da    = 100.*(accEEv[i]-accEEv[0])/accEEv[0];
        Double_t daErr = 100.*sqrt((accErrEEv[i])*(accErrEEv[i])+(accErrEEv[0])*(accErrEEv[0]))/accEEv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
      } else {
       txtfile << endl;
      } 
      txtfile << "        BE acceptance = " << setw(8) << nSelEventsBEv[i] << " / " << setw(8) << nGenEventsv[i] << " = " << setw(10) << accBEv[i] << " +/- " << accErrBEv[i];
      if(i>0) {
        Double_t da    = 100.*(accBEv[i]-accBEv[0])/accBEv[0];
        Double_t daErr = 100.*sqrt((accErrBEv[i])*(accErrBEv[i])+(accErrBEv[0])*(accErrBEv[0]))/accBEv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;
      } else {
        txtfile << endl;
      }
      txtfile << endl;
    }
    txtfile << endl;
    
    cout << " <> Output saved in " << outputDir << "/" << endl;
    cout << endl;
  }
  
  gBenchmark->Show("plotZeeSelect");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/select.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>SELECTION</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">SELECTION</h3>" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/pt.png\"><img src=\"select/pt.png\" alt=\"select/pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/eta.png\"><img src=\"select/eta.png\" alt=\"select/eta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/phi.png\"><img src=\"select/phi.png\" alt=\"select/phi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronmass.png\"><img src=\"select/dielectronmass.png\" alt=\"select/dielectronmass.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronpt.png\"><img src=\"select/dielectronpt.png\" alt=\"select/dielectronpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectrony.png\"><img src=\"select/dielectrony.png\" alt=\"select/dielectrony.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronphi.png\"><img src=\"select/dielectronphi.png\" alt=\"select/dielectronphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/ptBB.png\"><img src=\"select/ptBB.png\" alt=\"select/ptBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/etaBB.png\"><img src=\"select/etaBB.png\" alt=\"select/etaBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/phiBB.png\"><img src=\"select/phiBB.png\" alt=\"select/phiBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronmassBB.png\"><img src=\"select/dielectronmassBB.png\" alt=\"select/dielectronmassBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronptBB.png\"><img src=\"select/dielectronptBB.png\" alt=\"select/dielectronptBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronyBB.png\"><img src=\"select/dielectronyBB.png\" alt=\"select/dielectronyBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronphiBB.png\"><img src=\"select/dielectronphiBB.png\" alt=\"select/dielectronphiBB.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/ptEE.png\"><img src=\"select/ptEE.png\" alt=\"select/ptEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/etaEE.png\"><img src=\"select/etaEE.png\" alt=\"select/etaEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/phiEE.png\"><img src=\"select/phiEE.png\" alt=\"select/phiEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronmassEE.png\"><img src=\"select/dielectronmassEE.png\" alt=\"select/dielectronmassEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronptEE.png\"><img src=\"select/dielectronptEE.png\" alt=\"select/dielectronptEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronyEE.png\"><img src=\"select/dielectronyEE.png\" alt=\"select/dielectronyEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronphiEE.png\"><img src=\"select/dielectronphiEE.png\" alt=\"select/dielectronphiEE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/ptBE.png\"><img src=\"select/ptBE.png\" alt=\"select/ptBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/etaBE.png\"><img src=\"select/etaBE.png\" alt=\"select/etaBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/phiBE.png\"><img src=\"select/phiBE.png\" alt=\"select/phiBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronmassBE.png\"><img src=\"select/dielectronmassBE.png\" alt=\"select/dielectronmassBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronptBE.png\"><img src=\"select/dielectronptBE.png\" alt=\"select/dielectronptBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronyBE.png\"><img src=\"select/dielectronyBE.png\" alt=\"select/dielectronyBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"select/dielectronphiBE.png\"><img src=\"select/dielectronphiBE.png\" alt=\"select/dielectronphiBE.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
        
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
