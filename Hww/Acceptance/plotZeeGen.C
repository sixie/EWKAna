#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/CPlot.hh"          // helper class for plots
#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/MyTools.hh"        // miscellaneous helper functions

// define classes and constants to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#endif

//#define FILTER 40


//=== FUNCTION DECLARATIONS ======================================================================================

// Make graphs from ratio of distributions 
TGraphErrors* makeRatioGraph(TH1F* href, TH1F* h2);

// generate web page
void makeHTML(const TString outDir);


//=== MAIN MACRO =================================================================================================

void plotZeeGen(const TString input) 
{
  gBenchmark->Start("plotZeeGen");

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
  
  CPlot::sOutDir = outputDir + TString("/gen");

  const Double_t kGAP_LOW  = 1.4442;
  const Double_t kGAP_HIGH = 1.566;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  vector<TH1F*> hScale1v;
  vector<TH1F*> hScale2v;
  
  vector<TH1F*> hZMass1v;
  vector<TH1F*> hZMass2v;
  vector<TH1F*> hZPtv;   
  vector<TH1F*> hZyv;    
  vector<TH1F*> hZPhiv;  
  
  vector<TH1F*> hDielectronMass1v;
  vector<TH1F*> hDielectronMass2v;
  vector<TH1F*> hDielectronPtv;   
  vector<TH1F*> hDielectronyv;    
  vector<TH1F*> hDielectronPhiv;  
  
  vector<TH1F*> hElectron1Ptv; 
  vector<TH1F*> hElectron1Etav; 
  vector<TH1F*> hElectron1Phiv; 
  
  vector<TH1F*> hElectron2Ptv; 
  vector<TH1F*> hElectron2Etav; 
  vector<TH1F*> hElectron2Phiv; 
  
  vector<TH1F*> hNPhov;
  vector<TH1F*> hPhoPtv;
  vector<TH1F*> hPhoEtav;
  vector<TH1F*> hPhoPhiv;
  vector<TH1F*> hPhoDRv;
  vector<TH1F*> hPhoE1v;
  vector<TH1F*> hPhoE2v;
  vector<TH1F*> hPhoE3v;
  
  vector<TH1F*> hDecxv;
  vector<TH1F*> hDecyv;
  vector<TH1F*> hDeczv;
  
  vector<TH1F*> hNPartv;
  vector<TH1F*> hNChv;
    
  vector<TProfile*> pr_nGenPart_ptv;  
    
  vector<Double_t> nEventsv;  
  vector<Double_t> nPassv;
  vector<UInt_t>   nNegEvtsv;
  vector<UInt_t>   nZv;
  vector<Double_t> accv;
  vector<Double_t> accErrv;

  vector<Double_t> nPassBBv;
  vector<Double_t> accBBv;
  vector<Double_t> accErrBBv;
 
  vector<Double_t> nPassBEv;
  vector<Double_t> accBEv;
  vector<Double_t> accErrBEv;
 
  vector<Double_t> nPassEEv;
  vector<Double_t> accEEv;
  vector<Double_t> accErrEEv;
    
  const Int_t nbinsPhoE3 = 100;
  Double_t binxPhoE3[nbinsPhoE3+1];
  for(Int_t i=0; i<=nbinsPhoE3; i++) { binxPhoE3[i] = pow(10,-8+i*7/(Double_t)nbinsPhoE3); }

  char hname[100];
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    sprintf(hname,"hScale1_%i",ifile); hScale1v.push_back(new TH1F(hname,"",100,0,1000)); hScale1v[ifile]->Sumw2();
    sprintf(hname,"hScale2_%i",ifile); hScale2v.push_back(new TH1F(hname,"",50,0,200));   hScale2v[ifile]->Sumw2();
    
    sprintf(hname,"hZMass1_%i",ifile); hZMass1v.push_back(new TH1F(hname,"",100,0,1000)); hZMass1v[ifile]->Sumw2();
    sprintf(hname,"hZMass2_%i",ifile); hZMass2v.push_back(new TH1F(hname,"",100,0,200));  hZMass2v[ifile]->Sumw2();
    sprintf(hname,"hZPt_%i",ifile);    hZPtv.push_back(new TH1F(hname,"",50,0,100));      hZPtv[ifile]->Sumw2();
    sprintf(hname,"hZy_%i",ifile);     hZyv.push_back(new TH1F(hname,"",80,-8,8));        hZyv[ifile]->Sumw2();
    sprintf(hname,"hZPhi_%i",ifile);   hZPhiv.push_back(new TH1F(hname,"",80,-3.2,3.2));  hZPhiv[ifile]->Sumw2();

    sprintf(hname,"hDielectronMass1_%i",ifile); hDielectronMass1v.push_back(new TH1F(hname,"",100,0,1000)); hDielectronMass1v[ifile]->Sumw2();
    sprintf(hname,"hDielectronMass2_%i",ifile); hDielectronMass2v.push_back(new TH1F(hname,"",100,0,200));  hDielectronMass2v[ifile]->Sumw2();
    sprintf(hname,"hDielectronPt_%i",ifile);    hDielectronPtv.push_back(new TH1F(hname,"",50,0,100));      hDielectronPtv[ifile]->Sumw2();
    sprintf(hname,"hDielectrony_%i",ifile);     hDielectronyv.push_back(new TH1F(hname,"",80,-8,8));        hDielectronyv[ifile]->Sumw2();
    sprintf(hname,"hDielectronPhi_%i",ifile);   hDielectronPhiv.push_back(new TH1F(hname,"",80,-3.2,3.2));  hDielectronPhiv[ifile]->Sumw2();

    sprintf(hname,"hElectron1Pt_%i",ifile);  hElectron1Ptv.push_back(new TH1F(hname,"",50,0,200));     hElectron1Ptv[ifile]->Sumw2();
    sprintf(hname,"hElectron1Eta_%i",ifile); hElectron1Etav.push_back(new TH1F(hname,"",80,-8,8));     hElectron1Etav[ifile]->Sumw2();
    sprintf(hname,"hElectron1Phi_%i",ifile); hElectron1Phiv.push_back(new TH1F(hname,"",80,-3.2,3.2)); hElectron1Phiv[ifile]->Sumw2();
    
    sprintf(hname,"hElectron2Pt_%i",ifile);  hElectron2Ptv.push_back(new TH1F(hname,"",50,0,200));     hElectron2Ptv[ifile]->Sumw2();
    sprintf(hname,"hElectron2Eta_%i",ifile); hElectron2Etav.push_back(new TH1F(hname,"",80,-8,8));     hElectron2Etav[ifile]->Sumw2();
    sprintf(hname,"hElectron2Phi_%i",ifile); hElectron2Phiv.push_back(new TH1F(hname,"",80,-3.2,3.2)); hElectron2Phiv[ifile]->Sumw2();

    sprintf(hname,"hNPho_%i",ifile);   hNPhov.push_back(new TH1F(hname,"",5,-0.5,4.5));            hNPhov[ifile]->Sumw2();
    sprintf(hname,"hPhoPt_%i",ifile);  hPhoPtv.push_back(new TH1F(hname,"",100,0,20));             hPhoPtv[ifile]->Sumw2();
    sprintf(hname,"hPhoEta_%i",ifile); hPhoEtav.push_back(new TH1F(hname,"",80,-8,8));             hPhoEtav[ifile]->Sumw2();
    sprintf(hname,"hPhoPhi_%i",ifile); hPhoPhiv.push_back(new TH1F(hname,"",80,-3.2,3.2));         hPhoPhiv[ifile]->Sumw2();
    sprintf(hname,"hPhoDR_%i",ifile);  hPhoDRv.push_back(new TH1F(hname,"",100,0,1));              hPhoDRv[ifile]->Sumw2();
    sprintf(hname,"hPhoE1_%i",ifile);  hPhoE1v.push_back(new TH1F(hname,"",100,0,20));             hPhoE1v[ifile]->Sumw2();
    sprintf(hname,"hPhoE2_%i",ifile);  hPhoE2v.push_back(new TH1F(hname,"",100,0,0.2));            hPhoE2v[ifile]->Sumw2();
    sprintf(hname,"hPhoE3_%i",ifile);  hPhoE3v.push_back(new TH1F(hname,"",nbinsPhoE3,binxPhoE3)); hPhoE3v[ifile]->Sumw2();
    
    sprintf(hname,"hDecx_%i",ifile); hDecxv.push_back(new TH1F(hname,"",100,0,0.06));     hDecxv[ifile]->Sumw2();
    sprintf(hname,"hDecy_%i",ifile); hDecyv.push_back(new TH1F(hname,"",100,-0.02,0.02)); hDecyv[ifile]->Sumw2();
    sprintf(hname,"hDecz_%i",ifile); hDeczv.push_back(new TH1F(hname,"",100,-20,20));     hDeczv[ifile]->Sumw2();
    
    sprintf(hname,"hNPart_%i",ifile); hNPartv.push_back(new TH1F(hname,"",50,0,500)); hNPartv[ifile]->Sumw2();
    sprintf(hname,"hNCh_%i",ifile);   hNChv.push_back(new TH1F(hname,"",50,0,300));   hNChv[ifile]->Sumw2(); 
    
    sprintf(hname,"pr_nGenPart_pt_%i",ifile); pr_nGenPart_ptv.push_back(new TProfile(hname,"",40,0,200,0,500)); 
    
    nEventsv.push_back(0);
    nNegEvtsv.push_back(0);
    nPassv.push_back(0); 
    nPassBBv.push_back(0);
    nPassBEv.push_back(0); 
    nPassEEv.push_back(0);     
  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
    
  // Data structures to store info from TTrees
  mithep::TGenInfo *gen  = new mithep::TGenInfo();
  
  // loop over samples  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);
        
    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Gen",&gen);
    TBranch *genBr = eventTree->GetBranch("Gen");
  
    // loop over events    
    nZv.push_back(eventTree->GetEntries());
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      genBr->GetEntry(ientry);
#ifdef FILTER
      if(gen->mass < FILTER) continue;
#endif           
      if((gen->vmass < massLow) || (gen->vmass > massHigh)) continue;
      
      nEventsv[ifile] += gen->weight;
      if(gen->weight<0)
        nNegEvtsv[ifile]++;

      Bool_t isB1 = (fabs(gen->eta_1)<kGAP_LOW);
      Bool_t isB2 = (fabs(gen->eta_2)<kGAP_LOW);
      Bool_t isE1 = (fabs(gen->eta_1)>kGAP_HIGH);
      Bool_t isE2 = (fabs(gen->eta_2)>kGAP_HIGH);

      if((gen->pt_1>20) && (gen->pt_2>20)
         && ((fabs(gen->eta_1)<kGAP_LOW) || (fabs(gen->eta_1)>kGAP_HIGH))
	 && ((fabs(gen->eta_2)<kGAP_LOW) || (fabs(gen->eta_2)>kGAP_HIGH))   
	 && (fabs(gen->eta_1)<2.5) && (fabs(gen->eta_2)<2.5)) {
        
	nPassv[ifile] += gen->weight;
        if(isB1 && isB2)                          { nPassBBv[ifile] += gen->weight; } 
        else if(isE1 && isE2)                     { nPassEEv[ifile] += gen->weight; } 
        else if((isB1 && isE2) || (isE1 && isB2)) { nPassBEv[ifile] += gen->weight; }
      }
      
      hScale1v[ifile]->Fill(gen->scale,gen->weight);
      hScale2v[ifile]->Fill(gen->scale,gen->weight);
      
      hZMass1v[ifile]->Fill(gen->vmass,gen->weight);
      hZMass2v[ifile]->Fill(gen->vmass,gen->weight);
      hZPtv[ifile]   ->Fill(gen->vpt,  gen->weight);
      hZyv[ifile]    ->Fill(gen->vy,   gen->weight);
      hZPhiv[ifile]  ->Fill(gen->vphi, gen->weight);
    
      hDielectronMass1v[ifile]->Fill(gen->mass,gen->weight);
      hDielectronMass2v[ifile]->Fill(gen->mass,gen->weight);
      hDielectronPtv[ifile]   ->Fill(gen->pt,  gen->weight);
      hDielectronyv[ifile]    ->Fill(gen->y,   gen->weight);
      hDielectronPhiv[ifile]  ->Fill(gen->phi, gen->weight);

      hDecxv[ifile]->Fill(gen->decx,gen->weight);
      hDecyv[ifile]->Fill(gen->decy,gen->weight);
      hDeczv[ifile]->Fill(gen->decz,gen->weight);
    
      hElectron1Ptv[ifile] ->Fill(gen->pt_1, gen->weight); hElectron2Ptv[ifile] ->Fill(gen->pt_2, gen->weight);
      hElectron1Etav[ifile]->Fill(gen->eta_1,gen->weight); hElectron2Etav[ifile]->Fill(gen->eta_2,gen->weight);
      hElectron1Phiv[ifile]->Fill(gen->phi_1,gen->weight); hElectron2Phiv[ifile]->Fill(gen->phi_2,gen->weight);
      
      hNPhov[ifile]->Fill(gen->npho,gen->weight);
      if(gen->npho>0) {
        hPhoPtv[ifile] ->Fill(gen->phopt, gen->weight);
	hPhoEtav[ifile]->Fill(gen->phoeta,gen->weight);
	hPhoPhiv[ifile]->Fill(gen->phophi,gen->weight);
      
        TLorentzVector vecpho;
	vecpho.SetPtEtaPhiM(gen->phopt,gen->phoeta,gen->phophi,0);
	hPhoE1v[ifile]->Fill(vecpho.Energy(),gen->weight);
	hPhoE2v[ifile]->Fill(vecpho.Energy(),gen->weight);
	hPhoE3v[ifile]->Fill(vecpho.Energy(),gen->weight);
	
	Double_t dR1 = toolbox::deltaR(gen->phoeta,gen->phophi,gen->eta_1,gen->phi_1);
	Double_t dR2 = toolbox::deltaR(gen->phoeta,gen->phophi,gen->eta_2,gen->phi_2);
	hPhoDRv[ifile]->Fill((dR1<dR2) ? dR1 : dR2,gen->weight);
      }      
      
      hNPartv[ifile]->Fill(gen->nGenPart,gen->weight);
      hNChv[ifile]  ->Fill(gen->nGenCh,  gen->weight);   
      
      pr_nGenPart_ptv[ifile]->Fill(gen->vpt,gen->nGenPart,gen->weight);   
    }
   
    delete infile;
    infile=0, eventTree=0;
    
    accv.push_back(nPassv[ifile]/nEventsv[ifile]);
    accErrv.push_back(sqrt(accv[ifile]*(1-accv[ifile])/nEventsv[ifile])); 
    
    accBBv.push_back(nPassBBv[ifile]/nEventsv[ifile]);
    accErrBBv.push_back(sqrt(accBBv[ifile]*(1-accBBv[ifile])/nEventsv[ifile]));
    
    accBEv.push_back(nPassBEv[ifile]/nEventsv[ifile]);
    accErrBEv.push_back(sqrt(accBEv[ifile]*(1-accBEv[ifile])/nEventsv[ifile]));
    
    accEEv.push_back(nPassEEv[ifile]/nEventsv[ifile]);
    accErrEEv.push_back(sqrt(accEEv[ifile]*(1-accEEv[ifile])/nEventsv[ifile]));
  }
  delete gen;
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Event scale (Q)
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hScale1v[0]->GetBinWidth(1));
  CPlot plotScale1("scale1","","Q [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hScale1v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotScale1.AddHist1D(hScale1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotScale1.SetLogy();
  plotScale1.Draw(c,doSave,format);

  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hScale2v[0]->GetBinWidth(1));
  CPlot plotScale2("scale2","","Q [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hScale2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotScale2.AddHist1D(hScale2v[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotScale2.Draw(c,doSave,format);
     
  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMass1v[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZMass1v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZMass1.AddHist1D(hZMass1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c,doSave,format);
    
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMass2v[0]->GetBinWidth(1));
  CPlot plotZMass2("zmass2","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZMass2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZMass2.AddHist1D(hZMass2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass2.SetLogy();
  plotZMass2.Draw(c,doSave,format);
        
  // Z pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hZPtv[0]->GetBinWidth(1));
  CPlot plotZPt("zpt","","p_{T}(Z) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZPtv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZPt.AddHist1D(hZPtv[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotZPt.Draw(c,doSave,format);
  
  // Z rapidity
  sprintf(ylabel,"a.u. / %.2f",hZyv[0]->GetBinWidth(1));
  CPlot plotZy("zy","","y(Z)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZyv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZy.AddHist1D(hZyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZy.TransLegend(0.1,0);
  plotZy.Draw(c,doSave,format);
  
  // Z phi
  sprintf(ylabel,"a.u. / %.2f",hZPhiv[0]->GetBinWidth(1));
  CPlot plotZPhi("zphi","","#phi(Z)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZPhiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZPhi.AddHist1D(hZPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZPhi.SetYRange(0.8*(hZPhiv[0]->GetMaximum()),1.2*(hZPhiv[0]->GetMaximum()));
  plotZPhi.Draw(c,doSave,format); 

  
  // dielectron mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDielectronMass1v[0]->GetBinWidth(1));
  CPlot plotDielectronMass1("dielectronmass1","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDielectronMass1v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDielectronMass1.AddHist1D(hDielectronMass1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDielectronMass1.SetLogy();
  plotDielectronMass1.Draw(c,doSave,format);
  
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDielectronMass2v[0]->GetBinWidth(1));
  CPlot plotDielectronMass2("dielectronmass2","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDielectronMass2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDielectronMass2.AddHist1D(hDielectronMass2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDielectronMass2.SetLogy();
  plotDielectronMass2.Draw(c,doSave,format);
      
  // dielectron pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDielectronPtv[0]->GetBinWidth(1));
  CPlot plotDielectronPt("dielectronpt","","p_{T}(e^{+}e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDielectronPtv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDielectronPt.AddHist1D(hDielectronPtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDielectronPt.Draw(c,doSave,format);
  
  // dielectron rapidity
  sprintf(ylabel,"a.u. / %.2f",hDielectronyv[0]->GetBinWidth(1));
  CPlot plotDielectrony("dielectrony","","y(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDielectronyv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDielectrony.AddHist1D(hDielectronyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDielectrony.TransLegend(0.1,0);
  plotDielectrony.Draw(c,doSave,format);
  
  // dielectron phi
  sprintf(ylabel,"a.u. / %.2f",hDielectronPhiv[0]->GetBinWidth(1));
  CPlot plotDielectronPhi("dielectronphi","","#phi(e^{+}e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDielectronPhiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDielectronPhi.AddHist1D(hDielectronPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDielectronPhi.SetYRange(0.8*(hDielectronPhiv[0]->GetMaximum()),1.2*(hDielectronPhiv[0]->GetMaximum()));
  plotDielectronPhi.Draw(c,doSave,format);
  
  
  // electron pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hElectron1Ptv[0]->GetBinWidth(1));
  CPlot plotElectron1Pt("mu1pt","","p_{T}(e^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron1Ptv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron1Pt.AddHist1D(hElectron1Ptv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron1Pt.Draw(c,doSave,format);
  
  // electron eta
  sprintf(ylabel,"a.u. / %.2f",hElectron1Etav[0]->GetBinWidth(1)); 
  CPlot plotElectron1Eta("mu1eta","","#eta(e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron1Etav[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron1Eta.AddHist1D(hElectron1Etav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron1Eta.TransLegend(0.1,0);
  plotElectron1Eta.Draw(c,doSave,format);
  
  // electron phi
  sprintf(ylabel,"a.u. / %.2f",hElectron1Phiv[0]->GetBinWidth(1));
  CPlot plotElectron1Phi("mu1phi","","#phi(e^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron1Phiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron1Phi.AddHist1D(hElectron1Phiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron1Phi.SetYRange(0.8*(hElectron1Phiv[0]->GetMaximum()),1.2*(hElectron1Phiv[0]->GetMaximum()));
  plotElectron1Phi.Draw(c,doSave,format);  


  // positron pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hElectron2Ptv[0]->GetBinWidth(1));
  CPlot plotElectron2Pt("mu2pt","","p_{T}(e^{+}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron2Ptv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron2Pt.AddHist1D(hElectron2Ptv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron2Pt.Draw(c,doSave,format);
  
  // positron eta
  sprintf(ylabel,"a.u. / %.2f",hElectron2Etav[0]->GetBinWidth(1)); 
  CPlot plotElectron2Eta("mu2eta","","#eta(e^{+})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron2Etav[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron2Eta.AddHist1D(hElectron2Etav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron2Eta.TransLegend(0.1,0);
  plotElectron2Eta.Draw(c,doSave,format);
  
  // positron phi
  sprintf(ylabel,"a.u. / %.2f",hElectron2Phiv[0]->GetBinWidth(1));
  CPlot plotElectron2Phi("mu2phi","","#phi(e^{+})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hElectron2Phiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotElectron2Phi.AddHist1D(hElectron2Phiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotElectron2Phi.SetYRange(0.8*(hElectron2Phiv[0]->GetMaximum()),1.2*(hElectron2Phiv[0]->GetMaximum()));
  plotElectron2Phi.Draw(c,doSave,format); 
  

  // number of FSR photons
  CPlot plotNPho("npho","","N_{#gamma}","a.u");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNPhov[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotNPho.AddHist1D(hNPhov[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotNPho.Draw(c,doSave,format);
  
  // leading photon pT (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPhoPtv[0]->GetBinWidth(1));
  CPlot plotPhoPt("phopt","","leading photon p_{T} [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoPtv[i]->Scale(1.0/(Double_t)hPhoPtv[i]->GetEntries());
    plotPhoPt.AddHist1D(hPhoPtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoPt.SetLogy();
  plotPhoPt.Draw(c,doSave,format);
  
  // leading photon eta (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoEtav[0]->GetBinWidth(1)); 
  CPlot plotPhoEta("phoeta","","leading photon #eta",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoEtav[i]->Scale(1.0/(Double_t)hPhoEtav[i]->GetEntries());
    plotPhoEta.AddHist1D(hPhoEtav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoEta.TransLegend(0.1,0);
  plotPhoEta.Draw(c,doSave,format);
  
  // leading photon phi (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoPhiv[0]->GetBinWidth(1));
  CPlot plotPhoPhi("phophi","","leading photon #phi",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoPhiv[i]->Scale(1.0/(Double_t)hPhoPhiv[i]->GetEntries());
    plotPhoPhi.AddHist1D(hPhoPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoPhi.SetYRange(0.8*(hPhoPhiv[0]->GetMaximum()),1.2*(hPhoPhiv[0]->GetMaximum()));
  plotPhoPhi.Draw(c,doSave,format); 
  
  // leading photon-lepton dR (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoDRv[0]->GetBinWidth(1));
  CPlot plotPhoDR("phodr","","#DeltaR(#gamma,l)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoDRv[i]->Scale(1.0/(Double_t)hPhoDRv[i]->GetEntries());
    plotPhoDR.AddHist1D(hPhoDRv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoDR.SetLogy();
  plotPhoDR.Draw(c,doSave,format);
  
  // leading photon energy (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f GeV",hPhoE1v[0]->GetBinWidth(1)); 
  CPlot plotPhoE1("phoe1","","leading photon E [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE1v[i]->Scale(1.0/(Double_t)hPhoE1v[i]->GetEntries());
    plotPhoE1.AddHist1D(hPhoE1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE1.SetLogy();
  plotPhoE1.Draw(c,doSave,format);
  
  sprintf(ylabel,"a.u. / %.2f MeV",1000.*(hPhoE2v[0]->GetBinWidth(1)));
  CPlot plotPhoE2("phoe2","","leading photon E [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE2v[i]->Scale(1.0/(Double_t)hPhoE2v[i]->GetEntries());
    plotPhoE2.AddHist1D(hPhoE2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE2.SetLogy();
  plotPhoE2.Draw(c,doSave,format); 
  
  CPlot plotPhoE3("phoe3","","leading photon E [GeV]","a.u.");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE3v[i]->Scale(1.0/(Double_t)hPhoE3v[i]->GetEntries());
    plotPhoE3.AddHist1D(hPhoE3v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE3.SetLogy();
  plotPhoE3.SetLogx();
  plotPhoE3.Draw(c,doSave,format); 
  
    
  // Decay vertex
  sprintf(ylabel,"a.u. / %.5f cm",hDecxv[0]->GetBinWidth(1));
  CPlot plotDecx("decx","","decay vertex x [cm]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDecxv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDecx.AddHist1D(hDecxv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDecx.TransLegend(0.1,0);
  plotDecx.Draw(c,doSave,format);

  sprintf(ylabel,"a.u. / %.5f cm",hDecyv[0]->GetBinWidth(1));
  CPlot plotDecy("decy","","decay vertex y [cm]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDecyv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDecy.AddHist1D(hDecyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDecy.TransLegend(0.1,0);
  plotDecy.Draw(c,doSave,format);
  
  sprintf(ylabel,"a.u. / %.2f cm",hDeczv[0]->GetBinWidth(1));
  CPlot plotDecz("decz","","decay vertex z [cm]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDeczv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDecz.AddHist1D(hDeczv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDecz.TransLegend(0.1,0);
  plotDecz.Draw(c,doSave,format);  
  

  // Particle multiplicity
  sprintf(ylabel,"a.u. / %.1f",hNPartv[0]->GetBinWidth(1));
  CPlot plotNPart("npart","","N_{p} (status 1)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNPartv[i]->Scale(1.0/nEventsv[i]);
    plotNPart.AddHist1D(hNPartv[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotNPart.Draw(c,doSave,format);
  
  sprintf(ylabel,"a.u. / %.1f",hNChv[0]->GetBinWidth(1));
  CPlot plotNCh("nch","","N_{ch}",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNChv[i]->Scale(1.0/nEventsv[i]);
    plotNCh.AddHist1D(hNChv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotNCh.Draw(c,doSave,format);    

  // particle multiplicity vs. boson pT
  CPlot plotNGenPartVsScale("npart_v_pt","","p_{T}(Z) [GeV]","N_{p}");
  for(UInt_t i=0; i<fnamev.size(); i++) {
    plotNGenPartVsScale.AddProfile(pr_nGenPart_ptv[i],labelv[i],"",colorv[i],kFullDotSmall);
  }
  plotNGenPartVsScale.SetYRange(50,350);
  plotNGenPartVsScale.Draw(c,doSave,format);

  //
  // Ratio plots
  //
  // event scale (Q)
  CPlot plotScaleRatio1("scaleRatio1","","Q [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hScale1v[0],hScale1v[i]);
    plotScaleRatio1.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotScaleRatio1.SetYRange(0.5,1.5);
  plotScaleRatio1.AddLine(0,1,1100,1,kBlack,2);
  plotScaleRatio1.Draw(c,doSave,format);
    
  CPlot plotScaleRatio2("scaleRatio2","","Q [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hScale2v[0],hScale2v[i]);
    plotScaleRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotScaleRatio2.SetYRange(0.5,1.5);
  plotScaleRatio2.AddLine(0,1,220,1,kBlack,2);
  plotScaleRatio2.Draw(c,doSave,format);

  // Z mass
  CPlot plotZMassRatio1("zmassRatio1","","m(Z) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZMass1v[0],hZMass1v[i]);
    plotZMassRatio1.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZMassRatio1.SetYRange(0.5,1.5);
  plotZMassRatio1.AddLine(0,1,1100,1,kBlack,2);
  plotZMassRatio1.Draw(c,doSave,format);
    
  CPlot plotZMassRatio2("zmassRatio2","","m(Z) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZMass2v[0],hZMass2v[i]);
    plotZMassRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZMassRatio2.SetYRange(0.5,1.5);
  plotZMassRatio2.AddLine(0,1,220,1,kBlack,2);
  plotZMassRatio2.Draw(c,doSave,format);
        
  // Z pT
  CPlot plotZPtRatio("zptRatio","","p_{T}(Z) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZPtv[0],hZPtv[i]);
    plotZPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZPtRatio.SetYRange(0.5,1.5);
  plotZPtRatio.AddLine(0,1,110,1,kBlack,2);
  plotZPtRatio.Draw(c,doSave,format);
  
  // Z rapidity
  CPlot plotZyRatio("zyRatio","","y(Z)","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZyv[0],hZyv[i]);
    plotZyRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZyRatio.SetYRange(0.5,1.5);
  plotZyRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotZyRatio.Draw(c,doSave,format);

  // Dielectron mass
  CPlot plotDielectronMassRatio1("dielectronmassRatio1","","m(e^{+}e^{-}) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDielectronMass1v[0],hDielectronMass1v[i]);
    plotDielectronMassRatio1.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDielectronMassRatio1.SetYRange(0.5,1.5);
  plotDielectronMassRatio1.AddLine(0,1,1100,1,kBlack,2);
  plotDielectronMassRatio1.Draw(c,doSave,format);
    
  CPlot plotDielectronMassRatio2("dielectronmassRatio2","","m(e^{+}e^{-}) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDielectronMass2v[0],hDielectronMass2v[i]);
    plotDielectronMassRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDielectronMassRatio2.SetYRange(0.5,1.5);
  plotDielectronMassRatio2.AddLine(0,1,220,1,kBlack,2);
  plotDielectronMassRatio2.Draw(c,doSave,format);
        
  // Dielectron pT
  CPlot plotDielectronPtRatio("dielectronptRatio","","p_{T}(e^{+}e^{-}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDielectronPtv[0],hDielectronPtv[i]);
    plotDielectronPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDielectronPtRatio.SetYRange(0.5,1.5);
  plotDielectronPtRatio.AddLine(0,1,110,1,kBlack,2);
  plotDielectronPtRatio.Draw(c,doSave,format);
  
  // Dielectron rapidity
  CPlot plotDielectronyRatio("dielectronyRatio","","y(e^{+}e^{-})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDielectronyv[0],hDielectronyv[i]);
    plotDielectronyRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDielectronyRatio.SetYRange(0.5,1.5);
  plotDielectronyRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotDielectronyRatio.Draw(c,doSave,format);

  // electron pT
  CPlot plotElectron1PtRatio("mu1ptRatio","","p_{T}(e^{-}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hElectron1Ptv[0],hElectron1Ptv[i]);
    plotElectron1PtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotElectron1PtRatio.SetYRange(0.5,1.5);
  plotElectron1PtRatio.AddLine(0,1,220,1,kBlack,2);
  plotElectron1PtRatio.Draw(c,doSave,format);
  
  // electron eta 
  CPlot plotElectron1EtaRatio("mu1etaRatio","","#eta(e^{-})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hElectron1Etav[0],hElectron1Etav[i]);
    plotElectron1EtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotElectron1EtaRatio.SetYRange(0.5,1.5);
  plotElectron1EtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotElectron1EtaRatio.Draw(c,doSave,format);
  
  // positron pT
  CPlot plotElectron2PtRatio("mu2ptRatio","","p_{T}(e^{+}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hElectron2Ptv[0],hElectron2Ptv[i]);
    plotElectron2PtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotElectron2PtRatio.SetYRange(0.5,1.5);
  plotElectron2PtRatio.AddLine(0,1,220,1,kBlack,2);
  plotElectron2PtRatio.Draw(c,doSave,format);
  
  // positron eta 
  CPlot plotElectron2EtaRatio("mu2etaRatio","","#eta(e^{+})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hElectron2Etav[0],hElectron2Etav[i]);
    plotElectron2EtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotElectron2EtaRatio.SetYRange(0.5,1.5);
  plotElectron2EtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotElectron2EtaRatio.Draw(c,doSave,format);
  
  
  // photon pT
  CPlot plotPhoPtRatio("phoptRatio","","leading photon p_{T} [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoPtv[0],hPhoPtv[i]);
    plotPhoPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoPtRatio.SetYRange(0.5,1.5);
  plotPhoPtRatio.AddLine(0,1,22,1,kBlack,2);
  plotPhoPtRatio.Draw(c,doSave,format);
  
  // photon eta 
  CPlot plotPhoEtaRatio("phoetaRatio","","leading photon #eta","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoEtav[0],hPhoEtav[i]);
    plotPhoEtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoEtaRatio.SetYRange(0.5,1.5);
  plotPhoEtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotPhoEtaRatio.Draw(c,doSave,format);
  
  // photon-lepton dR
  CPlot plotPhoDRRatio("phodrRatio","","#DeltaR(#gamma,l)","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoDRv[0],hPhoDRv[i]);
    plotPhoDRRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoDRRatio.SetYRange(0.5,1.5);
  plotPhoDRRatio.AddLine(0,1,1.1,1,kBlack,2);
  plotPhoDRRatio.Draw(c,doSave,format);
  
  // photon energy 
  CPlot plotPhoE1Ratio("phoe1Ratio","","leading photon E [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoE1v[0],hPhoE1v[i]);
    plotPhoE1Ratio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoE1Ratio.SetYRange(0.5,1.5);
  plotPhoE1Ratio.AddLine(0,1,22,1,kBlack,2);
  plotPhoE1Ratio.Draw(c,doSave,format);
   
  CPlot plotPhoE2Ratio("phoe2Ratio","","leading photon E [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoE2v[0],hPhoE2v[i]);
    plotPhoE2Ratio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoE2Ratio.SetYRange(0.5,1.5);
  plotPhoE2Ratio.AddLine(0,1,0.22,1,kBlack,2);
  plotPhoE2Ratio.Draw(c,doSave,format);
            
      
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
    cout << labelv[ifile] << " file: " << fnamev[ifile] << endl;
    cout << "     Number of generated events: " << nZv[ifile] << endl;
    cout << "     Number of neg. wgt. events: " << nNegEvtsv[ifile] << endl;
    cout << "   Number of events preselected: " << nEventsv[ifile] << endl;
    cout << "       Number of events passing: " << nPassv[ifile] << endl; 
    cout << "             Overall Acceptance: " << setw(10) << accv[ifile] << " +/- " << accErrv[ifile]; 
    if(ifile>0) {
      Double_t da    = 100.*(accv[ifile]-accv[0])/accv[0];
      Double_t daErr = 100.*sqrt((accErrv[ifile])*(accErrv[ifile])+(accErrv[0])*(accErrv[0]))/accv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" ;        
    }
    cout << endl;
    cout << "                  BB Acceptance: " << setw(10) << accBBv[ifile] << " +/- " << accErrBBv[ifile]; 
    if(ifile>0) {
      Double_t da    = 100.*(accBBv[ifile]-accBBv[0])/accBBv[0];
      Double_t daErr = 100.*sqrt((accErrBBv[ifile])*(accErrBBv[ifile])+(accErrBBv[0])*(accErrBBv[0]))/accBBv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
    }
    cout << endl;
    cout << "                  EE Acceptance: " << setw(10) << accEEv[ifile] << " +/- " << accErrEEv[ifile]; 
    if(ifile>0) {
      Double_t da    = 100.*(accEEv[ifile]-accEEv[0])/accEEv[0];
      Double_t daErr = 100.*sqrt((accErrEEv[ifile])*(accErrEEv[ifile])+(accErrEEv[0])*(accErrEEv[0]))/accEEv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
    }
    cout << endl;
    cout << "                  BE Acceptance: " << setw(10) << accBEv[ifile] << " +/- " << accErrBEv[ifile]; 
    if(ifile>0) {
      Double_t da    = 100.*(accBEv[ifile]-accBEv[0])/accBEv[0];
      Double_t daErr = 100.*sqrt((accErrBEv[ifile])*(accErrBEv[ifile])+(accErrBEv[0])*(accErrBEv[0]))/accBEv[0];
      cout << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]" << endl;       
    }
    cout << endl;
  }
  cout << endl;
  
  if(doSave) {
    makeHTML(outputDir);
    
    ofstream txtfile;
    char txtfname[100];
    sprintf(txtfname,"%s/gen.txt",outputDir.Data());
    txtfile.open(txtfname);
    txtfile << "*" << endl;
    txtfile << "* SUMMARY" << endl;
    txtfile << "*--------------------------------------------------" << endl;
    txtfile << endl; 
    for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
      txtfile << labelv[ifile] << " file: " << fnamev[ifile] << endl;
      txtfile << "     Number of generated events: " << nZv[ifile] << endl;
      txtfile << "     Number of neg. wgt. events: " << nNegEvtsv[ifile] << endl;
      txtfile << "   Number of events preselected: " << nEventsv[ifile] << endl;
      txtfile << "       Number of events passing: " << nPassv[ifile] << endl; 
      txtfile << "             Overall Acceptance: " << setw(10) << accv[ifile] << " +/- " << accErrv[ifile]; 
      if(ifile>0) {
        Double_t da    = 100.*(accv[ifile]-accv[0])/accv[0];
        Double_t daErr = 100.*sqrt((accErrv[ifile])*(accErrv[ifile])+(accErrv[0])*(accErrv[0]))/accv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
      }
      txtfile << endl;
      txtfile << "                  BB Acceptance: " << setw(10) << accBBv[ifile] << " +/- " << accErrBBv[ifile]; 
      if(ifile>0) {
        Double_t da    = 100.*(accBBv[ifile]-accBBv[0])/accBBv[0];
        Double_t daErr = 100.*sqrt((accErrBBv[ifile])*(accErrBBv[ifile])+(accErrBBv[0])*(accErrBBv[0]))/accBBv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
      }
      txtfile << endl;
      txtfile << "                  EE Acceptance: " << setw(10) << accEEv[ifile] << " +/- " << accErrEEv[ifile]; 
      if(ifile>0) {
        Double_t da    = 100.*(accEEv[ifile]-accEEv[0])/accEEv[0];
        Double_t daErr = 100.*sqrt((accErrEEv[ifile])*(accErrEEv[ifile])+(accErrEEv[0])*(accErrEEv[0]))/accEEv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
      }
      txtfile << endl;
      txtfile << "                  BE Acceptance: " << setw(10) << accBEv[ifile] << " +/- " << accErrBEv[ifile]; 
      if(ifile>0) {
        Double_t da    = 100.*(accBEv[ifile]-accBEv[0])/accBEv[0];
        Double_t daErr = 100.*sqrt((accErrBEv[ifile])*(accErrBEv[ifile])+(accErrBEv[0])*(accErrBEv[0]))/accBEv[0];
        txtfile << setw(15) << "[ dA/A = (" << da << " +/- " << daErr << ") % ]";        
      }
      txtfile << endl;
    }      
    txtfile << endl;
    
    cout << " <> Output saved in " << outputDir << "/" << endl;
    cout << endl;
  }
  
  gBenchmark->Show("plotZeeGen");
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TGraphErrors* makeRatioGraph(TH1F* href, TH1F* h2) 
{
  char gname[50];
  sprintf(gname,"ratio_%s",h2->GetName());
  Double_t xval[href->GetNbinsX()];
  Double_t xerr[href->GetNbinsX()];
  Double_t yval[href->GetNbinsX()];
  Double_t yerr[href->GetNbinsX()];
  
  for(Int_t ibin=1; ibin<=href->GetNbinsX(); ibin++) {
    Double_t numer    = h2->GetBinContent(ibin);
    Double_t numerErr = h2->GetBinError(ibin);
    Double_t denom    = href->GetBinContent(ibin);
    Double_t denomErr = href->GetBinError(ibin);
    if(denom == 0) {
      if(numer == 0) {
        yval[ibin-1] = 1;
	yerr[ibin-1] = 0;
      } else {
        yval[ibin-1] = 0;
	yerr[ibin-1] = 0;
      }
    } else {
      yval[ibin-1] = numer/denom;
      yerr[ibin-1] =(numer/denom)*sqrt(numerErr*numerErr/numer/numer + denomErr*denomErr/denom/denom);
    }
    xval[ibin-1] = href->GetBinCenter(ibin);
    xerr[ibin-1] = 0.5*(href->GetBinWidth(ibin));    
  }
  
  return new TGraphErrors(href->GetNbinsX(),xval,yval,xerr,yerr);
}

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/gen.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>GEN</title></head>" << endl;
  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">GEN</h3>" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scale1.png\"><img src=\"gen/scale1.png\" alt=\"gen/scale1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scale2.png\"><img src=\"gen/scale2.png\" alt=\"gen/scale2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/npart.png\"><img src=\"gen/npart.png\" alt=\"gen/npart.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/nch.png\"><img src=\"gen/nch.png\" alt=\"gen/nch.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scaleRatio1.png\"><img src=\"gen/scaleRatio1.png\" alt=\"gen/scaleRatio1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scaleRatio2.png\"><img src=\"gen/scaleRatio2.png\" alt=\"gen/scaleRatio2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/npart_v_pt.png\"><img src=\"gen/npart_v_pt.png\" alt=\"gen/npart_v_pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zmass1.png\"><img src=\"gen/zmass1.png\" alt=\"gen/zmass1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zmass2.png\"><img src=\"gen/zmass2.png\" alt=\"gen/zmass2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zmassRatio1.png\"><img src=\"gen/zmassRatio1.png\" alt=\"gen/zmassRatio1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zmassRatio2.png\"><img src=\"gen/zmassRatio2.png\" alt=\"gen/zmassRatio2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zpt.png\"><img src=\"gen/zpt.png\" alt=\"gen/zpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zy.png\"><img src=\"gen/zy.png\" alt=\"gen/zy.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zphi.png\"><img src=\"gen/zphi.png\" alt=\"gen/zphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zptRatio.png\"><img src=\"gen/zptRatio.png\" alt=\"gen/zptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/zyRatio.png\"><img src=\"gen/zyRatio.png\" alt=\"gen/zyRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/decx.png\"><img src=\"gen/decx.png\" alt=\"gen/decx.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/decy.png\"><img src=\"gen/decy.png\" alt=\"gen/decy.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/decz.png\"><img src=\"gen/decz.png\" alt=\"gen/decz.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronmass1.png\"><img src=\"gen/dielectronmass1.png\" alt=\"gen/dielectronmass1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronmass2.png\"><img src=\"gen/dielectronmass2.png\" alt=\"gen/dielectronmass2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronmassRatio1.png\"><img src=\"gen/dielectronmassRatio1.png\" alt=\"gen/dielectronmassRatio1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronmassRatio2.png\"><img src=\"gen/dielectronmassRatio2.png\" alt=\"gen/dielectronmassRatio2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronpt.png\"><img src=\"gen/dielectronpt.png\" alt=\"gen/dielectronpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectrony.png\"><img src=\"gen/dielectrony.png\" alt=\"gen/dielectrony.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronphi.png\"><img src=\"gen/dielectronphi.png\" alt=\"gen/dielectronphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronptRatio.png\"><img src=\"gen/dielectronptRatio.png\" alt=\"gen/dielectronptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dielectronyRatio.png\"><img src=\"gen/dielectronyRatio.png\" alt=\"gen/dielectronyRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu1pt.png\"><img src=\"gen/mu1pt.png\" alt=\"gen/mu1pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu1eta.png\"><img src=\"gen/mu1eta.png\" alt=\"gen/mu1eta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu1phi.png\"><img src=\"gen/mu1phi.png\" alt=\"gen/mu1phi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu1ptRatio.png\"><img src=\"gen/mu1ptRatio.png\" alt=\"gen/mu1ptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu1etaRatio.png\"><img src=\"gen/mu1etaRatio.png\" alt=\"gen/mu1etaRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu2pt.png\"><img src=\"gen/mu2pt.png\" alt=\"gen/mu2pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu2eta.png\"><img src=\"gen/mu2eta.png\" alt=\"gen/mu2eta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu2phi.png\"><img src=\"gen/mu2phi.png\" alt=\"gen/mu2phi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu2ptRatio.png\"><img src=\"gen/mu2ptRatio.png\" alt=\"gen/mu2ptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/mu2etaRatio.png\"><img src=\"gen/mu2etaRatio.png\" alt=\"gen/mu2etaRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;  
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phopt.png\"><img src=\"gen/phopt.png\" alt=\"gen/phopt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoeta.png\"><img src=\"gen/phoeta.png\" alt=\"gen/phoeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phophi.png\"><img src=\"gen/phophi.png\" alt=\"gen/phophi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/npho.png\"><img src=\"gen/npho.png\" alt=\"gen/npho.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
    htmlfile << "<tr>" << endl;  
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoptRatio.png\"><img src=\"gen/phoptRatio.png\" alt=\"gen/phoptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoetaRatio.png\"><img src=\"gen/phoetaRatio.png\" alt=\"gen/phoetaRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phodr.png\"><img src=\"gen/phodr.png\" alt=\"gen/phodr.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoe1.png\"><img src=\"gen/phoe1.png\" alt=\"gen/phoe1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoe2.png\"><img src=\"gen/phoe2.png\" alt=\"gen/phoe2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoe3.png\"><img src=\"gen/phoe3.png\" alt=\"gen/phoe3.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phodrRatio.png\"><img src=\"gen/phodrRatio.png\" alt=\"gen/phodrRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoe1Ratio.png\"><img src=\"gen/phoe1Ratio.png\" alt=\"gen/phoe1Ratio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/phoe2Ratio.png\"><img src=\"gen/phoe2Ratio.png\" alt=\"gen/phoe2Ratio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
