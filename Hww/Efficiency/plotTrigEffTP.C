//================================================================================================
//
// Analyze single electron trigger efficiency with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes, 1D efficiency projections, efficiency in pT-eta
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/CPlot.hh"          // helper class for plots
#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"  

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// structure for output ntuple
#include "EffData.hh" 
#endif

//#define __WP95__
#define __WP80__

//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);

// Make 2D efficiency map
void makeEffHist2D(TH2F *hEff, TH2F *hErrl, TH2F *hErrh,
                   const TH2F *hNumer, const TH2F *hDenom, const Int_t method);


//=== MAIN MACRO ================================================================================================= 

void plotTrigEffTP(const TString conf      = "zee.conf",      // input file
                   const TString outputDir = "test/WP80_HLT", // output directory
	           const TString format    = "png"            // plot format
) {
  gBenchmark->Start("plotTrigEffTP");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // trigger requirement
  const UInt_t trigger  = 0x8000; //kHLT_Ele15_LW_L1R;
 
  // histogram ranges
  const Double_t muPtMin  = 20;
  const Double_t muPtMax  = 200;
  const Double_t muEtaMin = -2.5;
  const Double_t muEtaMax = 2.5;
  
  // mass region
  const Double_t massLo = 60;
  const Double_t massHi = 120;
  
  // cuts for computing overall electron efficiency in smaller phase space
  const Double_t ptCut  = 20;
  const Double_t etaCut = 2.5;
  
  // efficiency error calculation method
  // method: 0 -> Bayes Divide
  //         1 -> Feldman-Cousins 
  //         2 -> Clopper-Pearson  
  const Int_t method=0;

  Double_t lumi;              // luminosity (pb^-1)
  Bool_t   doWeight;          // weight events?  
  
  vector<TString>  fnamev;    // sample files 
  vector<Double_t> xsecv;     // per file cross section
  vector<TString>  jsonv;     // per file JSON file
  vector<Double_t> weightv;   // per file event weight

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
      stringstream ss2(line); ss2 >> doWeight;
      
    } else if(state==1) {  // define data sample
      string fname;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> xsec >> json;
      fnamev.push_back(fname);
      xsecv.push_back(xsec);
      jsonv.push_back(json);    
    }
  }
  ifs.close();

  CPlot::sOutDir = outputDir + TString("/plots");

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
 
  //  
  // Set up histograms
  // 
  TH1F *hNumerPt    = new TH1F("hNumerPt","",50,muPtMin,muPtMax);
  TH1F *hDenomPt    = new TH1F("hDenomPt","",50,muPtMin,muPtMax);
  TH1F *hNumerEta   = new TH1F("hNumerEta","",50,muEtaMin,muEtaMax);
  TH1F *hDenomEta   = new TH1F("hDenomEta","",50,muEtaMin,muEtaMax);
  TH1F *hNumerPhi   = new TH1F("hNumerPhi","",40,-3.2,3.2);
  TH1F *hDenomPhi   = new TH1F("hDenomPhi","",40,-3.2,3.2);
  TH2F *hNumerPtEta = new TH2F("hNumerPtEta","",50,muPtMin,muPtMax,50,muEtaMin,muEtaMax);
  TH2F *hDenomPtEta = new TH2F("hDenomPtEta","",50,muPtMin,muPtMax,50,muEtaMin,muEtaMax);   
  
  //
  // Set up output ntuple
  //
  TTree *outTree = new TTree("Events","Events");
  outTree->SetDirectory(0);
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:pass/i:tag");

  //
  // Counters
  //
  Double_t numer=0, denom=0;
  
  
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
    eventTree->SetBranchAddress("Info",    &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",     &gen);         TBranch *genBr      = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Electron",&electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    
    // Determine maximum number of events to consider
    // *** CASES ***
    // <> lumi < 0                             => use all events in the sample
    // <> xsec = 0                             => for data (use all events)
    // <> lumi > 0, xsec > 0, doWeight = true  => use all events and scale to lumi
    // <> lumi > 0, xsec > 0, doWeight = false => compute expected number of events (n)
    UInt_t maxEvents = eventTree->GetEntries();
    const Double_t xsec = xsecv[ifile];
    Double_t weight = 1;
    if(lumi>0) { 
      if(xsec>0) {
        if(doWeight) {
          weight = lumi*xsec/(Double_t)eventTree->GetEntries();
        } else {
          maxEvents = (UInt_t)(lumi*xsec);
        }
      }      
    }
    if(maxEvents > eventTree->GetEntries()) {
      cout << "Not enough events for " << lumi << " pb^-1 in file: " << fnamev[ifile] << endl;
      return;
    }
    weightv.push_back(weight);

    // loop over events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if(ientry >= maxEvents) break;

      infoBr->GetEntry(ientry);
      
      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

      // trigger requirement
      if(!(info->triggerBits & trigger)) continue;      
  
      genBr->GetEntry(ientry);
      const Double_t dRMax = 0.3;
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
        
	if((fabs(tag->scEta)>kGAP_LOW) && (fabs(tag->scEta)<kGAP_HIGH)) continue;
	Bool_t tagEB = (fabs(tag->scEta)<kGAP_LOW);
	
	if(tag->scEt        < 20)  continue;
	if(fabs(tag->scEta) > 2.5) continue;
	if(!(tag->isEcalDriven))   continue;
	
//	if((toolbox::deltaR(gen->eta_1, gen->phi_1, tag->eta, tag->phi)>dRMax) 
//	    && (toolbox::deltaR(gen->eta_2, gen->phi_2, tag->eta, tag->phi)>dRMax)) continue;

#ifdef __WP95__		 
        if(tag->nExpHitsInner > 1) continue;	               
	if(tagEB) {
	  if(tag->trkIso03         > 0.15*(tag->pt)) continue;
	  if(tag->emIso03          > 2.00*(tag->pt)) continue;
	  if(tag->hadIso03         > 0.12*(tag->pt)) continue;
	  if(tag->sigiEtaiEta      > 0.01)           continue;
//	  if(fabs(tag->deltaPhiIn) > 0.8)            continue;
	  if(fabs(tag->deltaEtaIn) > 0.007)          continue;
	  if(tag->HoverE           > 0.15)           continue;
	} else {
	  if(tag->trkIso03         > 0.08*(tag->pt)) continue;
	  if(tag->emIso03          > 0.06*(tag->pt)) continue;
	  if(tag->hadIso03         > 0.05*(tag->pt)) continue;
	  if(tag->sigiEtaiEta      > 0.03)           continue;
//	  if(fabs(tag->deltaPhiIn) > 0.7)            continue;
//	  if(fabs(tag->deltaEtaIn) > 0.01)           continue;
	  if(tag->HoverE           > 0.07)           continue;	  
	}            
#endif
#ifdef __WP80__ 
        if(tag->nExpHitsInner > 0) continue;
	if((fabs(tag->partnerDist)<0.02) && (fabs(tag->partnerDeltaCot)<0.02)) continue; 	        
      	if(tagEB) {
	  if(tag->trkIso03         > 0.09*(tag->pt)) continue;
	  if(tag->emIso03          > 0.07*(tag->pt)) continue;
	  if(tag->hadIso03         > 0.10*(tag->pt)) continue;
	  if(tag->sigiEtaiEta      > 0.01)           continue;
	  if(fabs(tag->deltaPhiIn) > 0.06)           continue;
	  if(fabs(tag->deltaEtaIn) > 0.004)          continue;
	  if(tag->HoverE           > 0.04)           continue;
	} else {
	  if(tag->trkIso03         > 0.04*(tag->pt))  continue;
	  if(tag->emIso03          > 0.05*(tag->pt))  continue;
	  if(tag->hadIso03         > 0.025*(tag->pt)) continue;
	  if(tag->sigiEtaiEta      > 0.03)            continue;
	  if(fabs(tag->deltaPhiIn) > 0.03)            continue;
//	  if(fabs(tag->deltaEtaIn) > 0.007)           continue;
	  if(tag->HoverE           > 0.025)           continue;	  
	}  
#endif    
       
	if(!(tag->hltMatchBits & trigger)) continue;
	
	const Double_t m = 0.105658369;
	TLorentzVector vtag;
	vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
	
	for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
	  if(i==j) continue;
	  
	  const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);

	  if((fabs(probe->scEta)>kGAP_LOW) && (fabs(probe->scEta)<kGAP_HIGH)) continue;
	  Bool_t probeEB = (fabs(probe->scEta)<kGAP_LOW); 
	  
	  if(probe->scEt        < 20)  continue;
	  if(fabs(probe->scEta) > 2.5) continue;
	  if(!(probe->isEcalDriven))   continue;  

//          if((toolbox::deltaR(gen->eta_1, gen->phi_1, probe->eta, probe->phi)>dRMax) 
//	    && (toolbox::deltaR(gen->eta_2, gen->phi_2, probe->eta, probe->phi)>dRMax)) continue;
	    	
#ifdef __WP95__		  
          if(probe->nExpHitsInner > 1) continue;	                  
	  if(probeEB) {
	    if(probe->trkIso03         > 0.15*(probe->pt)) continue;
	    if(probe->emIso03          > 2.00*(probe->pt)) continue;
	    if(probe->hadIso03         > 0.12*(probe->pt)) continue;
	    if(probe->sigiEtaiEta      > 0.01)             continue;
//	    if(fabs(probe->deltaPhiIn) > 0.8)              continue;
	    if(fabs(probe->deltaEtaIn) > 0.007)            continue;
	    if(probe->HoverE           > 0.15)             continue;
	  } else {
	    if(probe->trkIso03         > 0.08*(probe->pt)) continue;
	    if(probe->emIso03          > 0.06*(probe->pt)) continue;
	    if(probe->hadIso03         > 0.05*(probe->pt)) continue;
	    if(probe->sigiEtaiEta      > 0.03)             continue;
//	    if(fabs(probe->deltaPhiIn) > 0.7)              continue;
//	    if(fabs(probe->deltaEtaIn) > 0.001)            continue;
	    if(probe->HoverE           > 0.07)             continue;	  
	  }
#endif
#ifdef __WP80__
          if(probe->nExpHitsInner > 0) continue;
	  if((fabs(probe->partnerDist)<0.02) && (fabs(probe->partnerDeltaCot)<0.02)) continue; 	        
      	  if(probeEB) {
	    if(probe->trkIso03         > 0.09*(probe->pt)) continue;
	    if(probe->emIso03          > 0.07*(probe->pt)) continue;
	    if(probe->hadIso03         > 0.10*(probe->pt)) continue;
	    if(probe->sigiEtaiEta      > 0.01)             continue;
	    if(fabs(probe->deltaPhiIn) > 0.06)             continue;
	    if(fabs(probe->deltaEtaIn) > 0.004)            continue;
	    if(probe->HoverE           > 0.04)             continue;
	  } else {
	    if(probe->trkIso03         > 0.04*(probe->pt))  continue;
	    if(probe->emIso03          > 0.05*(probe->pt))  continue;
	    if(probe->hadIso03         > 0.025*(probe->pt)) continue;
	    if(probe->sigiEtaiEta      > 0.03)              continue;
	    if(fabs(probe->deltaPhiIn) > 0.03)              continue;
//	    if(fabs(probe->deltaEtaIn) > 0.007)             continue;
	    if(probe->HoverE           > 0.025)             continue;	  
	  }  
#endif   
	  
	  TLorentzVector vprobe;
	  vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
	  
	  TLorentzVector vDiEle = vtag + vprobe;
	  if((vDiEle.M()<massLo) || (vDiEle.M()>massHi)) continue;
	  
	  Bool_t pass = (probe->hltMatchBits & trigger);
	    
          // Fill histograms
          hDenomPt->Fill(probe->pt); 
          hDenomEta->Fill(probe->eta);
          hDenomPhi->Fill(probe->phi);
          hDenomPtEta->Fill(probe->pt,probe->eta);
          if(pass) {
            hNumerPt->Fill(probe->pt); 
            hNumerEta->Fill(probe->eta);
            hNumerPhi->Fill(probe->phi);
            hNumerPtEta->Fill(probe->pt,probe->eta);
          }
         
          // Fill tree
	  data.mass   = vDiEle.M();
	  data.pt     = probe->pt;
          data.eta    = probe->eta;
          data.phi    = probe->phi;
          data.weight = weightv[ifile];
          data.pass   = (pass) ? 1 : 0;
          data.tag    = (pass) ? 1 : 0;
	  outTree->Fill();
    
          if(probe->pt>ptCut && fabs(probe->eta)<etaCut) { 
            denom+=weight;
            if(pass) numer+=weight; 
          }
	}
      }
    }
    delete infile;
    infile=0, eventTree=0;    
  }
  delete info;
  delete electronArr;


  // Compute efficiencies
  cout << "Computing efficiencies..."; 
  cout.flush(); 

  TH2F *hEffPtEta = new TH2F("hEffPtEta","",
                             hNumerPtEta->GetNbinsX(),hNumerPtEta->GetXaxis()->GetXmin(),hNumerPtEta->GetXaxis()->GetXmax(),
                             hNumerPtEta->GetNbinsY(),hNumerPtEta->GetYaxis()->GetXmin(),hNumerPtEta->GetYaxis()->GetXmax());
  TH2F *hErrlPtEta = (TH2F*)hEffPtEta->Clone("hErrlPtEta");
  TH2F *hErrhPtEta = (TH2F*)hEffPtEta->Clone("hErrhPtEta");                                                  
  makeEffHist2D(hEffPtEta,hErrlPtEta,hErrhPtEta,hNumerPtEta,hDenomPtEta,method);
  
  TGraphAsymmErrors *grEffPt  = toolbox::makeEffGraph(hNumerPt, hDenomPt, method); grEffPt->SetName("effPt");
  TGraphAsymmErrors *grEffEta = toolbox::makeEffGraph(hNumerEta,hDenomEta,method); grEffEta->SetName("effEta");
  TGraphAsymmErrors *grEffPhi = toolbox::makeEffGraph(hNumerPhi,hDenomPhi,method); grEffPhi->SetName("effPhi");
      
  Double_t eff, errl, errh;  
  eff = toolbox::calcEff(numer,denom,&errl,&errh,method);

  cout << "DONE!" << endl;
  cout << endl; 
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);

  // Electron efficiency in pT
  CPlot plotEffPt("effpt","","electron p_{T} [GeV/c]","#varepsilon");
  plotEffPt.AddGraph(grEffPt,"",kBlack,kFullDotMedium);
  plotEffPt.Draw(c,kTRUE,format);
  
  // Electron efficiency in eta
  CPlot plotEffEta("effeta","","electron #eta","#varepsilon");
  plotEffEta.AddGraph(grEffEta,"",kBlack,kFullDotMedium);
  plotEffEta.Draw(c,kTRUE,format);
  
  // Electron efficiency in phi
  CPlot plotEffPhi("effphi","","electron #phi","#varepsilon");
  plotEffPhi.AddGraph(grEffPhi,"",kBlack,kFullDotMedium);
  plotEffPhi.Draw(c,kTRUE,format);
  
  //
  // pT-eta electron efficiency maps
  //     
  gStyle->SetPalette(1);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);

  hEffPtEta->SetTitleOffset(1.2,"Y");
  CPlot plotEffPtEta("effpteta","","p_{T} [GeV/c]","#eta");
  plotEffPtEta.AddHist2D(hEffPtEta,"COLZ");
  plotEffPtEta.Draw(c,kTRUE,format);    

  hErrlPtEta->SetTitleOffset(1.2,"Y");
  CPlot plotErrlPtEta("errlpteta","","p_{T} [GeV/c]","#eta");
  plotErrlPtEta.AddHist2D(hErrlPtEta,"COLZ");
  plotErrlPtEta.Draw(c,kTRUE,format);
  
  hErrhPtEta->SetTitleOffset(1.2,"Y");
  CPlot plotErrhPtEta("errhpteta","","p_{T} [GeV/c]","#eta");
  plotErrhPtEta.AddHist2D(hErrhPtEta,"COLZ");
  plotErrhPtEta.Draw(c,kTRUE,format);
  
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  cout << "   electron efficiency (pT>" << ptCut << ", |eta|<" << etaCut << "):" << endl;
  cout << "      " << numer << " / " << denom << " = " << eff << " (+" << errh << ", -" << errl << ")" << endl;
  cout << endl;
  if(lumi>0) {
    cout << " L_int = " << lumi << "/pb" << endl;
    cout << endl;
  }
        
  // make webpage
  makeHTML(outputDir);

  // print summary text file
  ofstream txtfile;
  char txtfname[100];
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl; 
  if(lumi>0) {
    txtfile << " L_int = " << lumi << "/pb" << endl;
    txtfile << endl;
  }
  txtfile << "   electron efficiency (pT>" << ptCut << ", |eta|<" << etaCut << "):" << endl;
  txtfile << "      " << numer << " / " << denom << " = " << eff << " (+" << errh << ", -" << errl << ")" << endl;
  txtfile << endl; 
  txtfile.close();
        
  // write output ROOT file  
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/eff.root"),"RECREATE"); 
  grEffPt->Write();
  grEffEta->Write();
  grEffPhi->Write();
  hEffPtEta->Write();
  hErrlPtEta->Write();
  hErrhPtEta->Write();
  outTree->Write();
  outFile->Close();
  delete outFile; 
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("plotTrigEffTP"); 
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeEffHist2D(TH2F *hEff, TH2F *hErrl, TH2F *hErrh,
                   const TH2F *hNumer, const TH2F *hDenom, const Int_t method)
{
  const Int_t npt  = hNumer->GetNbinsX();
  const Int_t neta = hNumer->GetNbinsY();  
                                                                          
  for(Int_t ipt=1; ipt<=npt; ipt++) {
    for(Int_t ieta=1; ieta<=neta; ieta++) {
      const Int_t ibin = hNumer->GetBin(ipt,ieta);
      Double_t numer = hNumer->GetBinContent(ibin);
      Double_t denom = hDenom->GetBinContent(ibin);
      Double_t errl, errh;
      Double_t eff = toolbox::calcEff(numer,denom,&errl,&errh,method);
      hEff->SetBinContent(ibin,eff);
      hErrl->SetBinContent(ibin,errl);
      hErrh->SetBinContent(ibin,errh);
    }
  }
}

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effpt.png\"><img src=\"plots/effpt.png\" alt=\"plots/effpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effeta.png\"><img src=\"plots/effeta.png\" alt=\"plots/effeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effphi.png\"><img src=\"plots/effphi.png\" alt=\"plots/effphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;    
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/effpteta.png\"><img src=\"plots/effpteta.png\" alt=\"plots/effpteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errlpteta.png\"><img src=\"plots/errlpteta.png\" alt=\"plots/errlpteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/errhpteta.png\"><img src=\"plots/errhpteta.png\" alt=\"plots/errhpteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
}
