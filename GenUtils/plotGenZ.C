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

#include "GenStructDefs.hh"
#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// Make graphs from ratio of distributions   
TGraphErrors* makeRatioGraph(TH1F* href, TH1F* h2);

// generate web page
void makeHTML(const TString outDir);


//=== MAIN MACRO =================================================================================================

void plotGenZ(const TString input) 
{
  gBenchmark->Start("plotGenZ");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  TString format = "png";  // output file format
  
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
  
  vector<TH1F*> hDileptonMass1v;
  vector<TH1F*> hDileptonMass2v;
  vector<TH1F*> hDileptonPtv;   
  vector<TH1F*> hDileptonyv;    
  vector<TH1F*> hDileptonPhiv;  
  
  vector<TH1F*> hLepton1Ptv; 
  vector<TH1F*> hLepton1Etav; 
  vector<TH1F*> hLepton1Phiv; 
  
  vector<TH1F*> hLepton2Ptv; 
  vector<TH1F*> hLepton2Etav; 
  vector<TH1F*> hLepton2Phiv; 
  
  vector<TH1F*> hNPhov;
  vector<TH1F*> hPhoPtv;
  vector<TH1F*> hPhoEtav;
  vector<TH1F*> hPhoPhiv;
  vector<TH1F*> hPhoDRv;
  vector<TH1F*> hPhoE1v;
  vector<TH1F*> hPhoE2v;
  vector<TH1F*> hPhoE3v;
    
  vector<UInt_t> nEventsv;
  vector<UInt_t> nNegEvtsv;
  vector<UInt_t> nPassv;
  vector<UInt_t> nZv;
  vector<Double_t> accv;
  vector<Double_t> accErrv;

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

    sprintf(hname,"hDileptonMass1_%i",ifile); hDileptonMass1v.push_back(new TH1F(hname,"",100,0,1000)); hDileptonMass1v[ifile]->Sumw2();
    sprintf(hname,"hDileptonMass2_%i",ifile); hDileptonMass2v.push_back(new TH1F(hname,"",100,0,200));  hDileptonMass2v[ifile]->Sumw2();
    sprintf(hname,"hDileptonPt_%i",ifile);    hDileptonPtv.push_back(new TH1F(hname,"",50,0,100));      hDileptonPtv[ifile]->Sumw2();
    sprintf(hname,"hDileptony_%i",ifile);     hDileptonyv.push_back(new TH1F(hname,"",80,-8,8));        hDileptonyv[ifile]->Sumw2();
    sprintf(hname,"hDileptonPhi_%i",ifile);   hDileptonPhiv.push_back(new TH1F(hname,"",80,-3.2,3.2));  hDileptonPhiv[ifile]->Sumw2();

    sprintf(hname,"hLepton1Pt_%i",ifile);  hLepton1Ptv.push_back(new TH1F(hname,"",50,0,200));     hLepton1Ptv[ifile]->Sumw2();
    sprintf(hname,"hLepton1Eta_%i",ifile); hLepton1Etav.push_back(new TH1F(hname,"",80,-8,8));     hLepton1Etav[ifile]->Sumw2();
    sprintf(hname,"hLepton1Phi_%i",ifile); hLepton1Phiv.push_back(new TH1F(hname,"",80,-3.2,3.2)); hLepton1Phiv[ifile]->Sumw2();
    
    sprintf(hname,"hLepton2Pt_%i",ifile);  hLepton2Ptv.push_back(new TH1F(hname,"",50,0,200));     hLepton2Ptv[ifile]->Sumw2();
    sprintf(hname,"hLepton2Eta_%i",ifile); hLepton2Etav.push_back(new TH1F(hname,"",80,-8,8));     hLepton2Etav[ifile]->Sumw2();
    sprintf(hname,"hLepton2Phi_%i",ifile); hLepton2Phiv.push_back(new TH1F(hname,"",80,-3.2,3.2)); hLepton2Phiv[ifile]->Sumw2();

    sprintf(hname,"hNPho_%i",ifile);   hNPhov.push_back(new TH1F(hname,"",5,-0.5,4.5));            hNPhov[ifile]->Sumw2();
    sprintf(hname,"hPhoPt_%i",ifile);  hPhoPtv.push_back(new TH1F(hname,"",100,0,20));             hPhoPtv[ifile]->Sumw2();
    sprintf(hname,"hPhoEta_%i",ifile); hPhoEtav.push_back(new TH1F(hname,"",80,-8,8));             hPhoEtav[ifile]->Sumw2();
    sprintf(hname,"hPhoPhi_%i",ifile); hPhoPhiv.push_back(new TH1F(hname,"",80,-3.2,3.2));         hPhoPhiv[ifile]->Sumw2();
    sprintf(hname,"hPhoDR_%i",ifile);  hPhoDRv.push_back(new TH1F(hname,"",100,0,1));              hPhoDRv[ifile]->Sumw2();
    sprintf(hname,"hPhoE1_%i",ifile);  hPhoE1v.push_back(new TH1F(hname,"",100,0,20));             hPhoE1v[ifile]->Sumw2();
    sprintf(hname,"hPhoE2_%i",ifile);  hPhoE2v.push_back(new TH1F(hname,"",100,0,0.2));            hPhoE2v[ifile]->Sumw2();
    sprintf(hname,"hPhoE3_%i",ifile);  hPhoE3v.push_back(new TH1F(hname,"",nbinsPhoE3,binxPhoE3)); hPhoE3v[ifile]->Sumw2();
        
    nEventsv.push_back(0);
    nNegEvtsv.push_back(0);
    nPassv.push_back(0);    
  }

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
    
  // Data structures to store info
  TGenData data;
  
  // loop over samples  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);
        
    // Get the TTree
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info 
    eventTree->SetBranchAddress("Events", &data.npho); 
  
    // loop over events    
    nZv.push_back(eventTree->GetEntries());
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);
     
      if((data.vmass < massLow) || (data.vmass > massHigh)) continue;
      
      nEventsv[ifile] += data.weight;
      if(data.weight<0)
        nNegEvtsv[ifile]++;

      if((data.pt_1>20) && (data.pt_2>20) 
         && (fabs(data.eta_1)<2.4) && (fabs(data.eta_2)<2.4) 
//	 && (data.mass>60) && (data.mass<120)
        )
        nPassv[ifile] += data.weight;
      
      hScale1v[ifile]->Fill(data.scalePdf,data.weight);
      hScale2v[ifile]->Fill(data.scalePdf,data.weight);
      
      hZMass1v[ifile]->Fill(data.vmass,data.weight);
      hZMass2v[ifile]->Fill(data.vmass,data.weight);
      hZPtv[ifile]   ->Fill(data.vpt,  data.weight);
      hZyv[ifile]    ->Fill(data.vy,   data.weight);
      hZPhiv[ifile]  ->Fill(data.vphi, data.weight);
    
      hDileptonMass1v[ifile]->Fill(data.mass,data.weight);
      hDileptonMass2v[ifile]->Fill(data.mass,data.weight);
      hDileptonPtv[ifile]   ->Fill(data.pt,  data.weight);
      hDileptonyv[ifile]    ->Fill(data.y,   data.weight);
      hDileptonPhiv[ifile]  ->Fill(data.phi, data.weight);
    
      hLepton1Ptv[ifile] ->Fill(data.pt_1,  data.weight); hLepton2Ptv[ifile] ->Fill(data.pt_2, data.weight);
      hLepton1Etav[ifile]->Fill(data.eta_1,data.weight);  hLepton2Etav[ifile]->Fill(data.eta_2,data.weight);
      hLepton1Phiv[ifile]->Fill(data.phi_1,data.weight);  hLepton2Phiv[ifile]->Fill(data.phi_2,data.weight);
      
      hNPhov[ifile]->Fill(data.npho,data.weight);
      if(data.npho>0) {
        hPhoPtv[ifile] ->Fill(data.phopt, data.weight);
	hPhoEtav[ifile]->Fill(data.phoeta,data.weight);
	hPhoPhiv[ifile]->Fill(data.phophi,data.weight);
      
        TLorentzVector vecpho;
	vecpho.SetPtEtaPhiM(data.phopt,data.phoeta,data.phophi,0);
	hPhoE1v[ifile]->Fill(vecpho.Energy(),data.weight);
	hPhoE2v[ifile]->Fill(vecpho.Energy(),data.weight);
	hPhoE3v[ifile]->Fill(vecpho.Energy(),data.weight);
	
	Double_t dR1 = toolbox::deltaR(data.phoeta,data.phophi,data.eta_1,data.phi_1);
	Double_t dR2 = toolbox::deltaR(data.phoeta,data.phophi,data.eta_2,data.phi_2);
	hPhoDRv[ifile]->Fill((dR1<dR2) ? dR1 : dR2,data.weight);
      }      
    }
   
    delete infile;
    infile=0, eventTree=0;
    
    accv.push_back((Double_t)nPassv[ifile]/(Double_t)nEventsv[ifile]);
    accErrv.push_back(sqrt(accv[ifile]*(1-accv[ifile])/(Double_t)nEventsv[ifile])); 
  }
  
  
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
  plotScale1.Draw(c,kTRUE,format);

  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hScale2v[0]->GetBinWidth(1));
  CPlot plotScale2("scale2","","Q [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hScale2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotScale2.AddHist1D(hScale2v[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotScale2.Draw(c,kTRUE,format);
     
  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMass1v[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZMass1v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZMass1.AddHist1D(hZMass1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMass2v[0]->GetBinWidth(1));
  CPlot plotZMass2("zmass2","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZMass2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZMass2.AddHist1D(hZMass2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass2.SetLogy();
  plotZMass2.Draw(c,kTRUE,format);
        
  // Z pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hZPtv[0]->GetBinWidth(1));
  CPlot plotZPt("zpt","","p_{T}(Z) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZPtv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZPt.AddHist1D(hZPtv[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotZPt.Draw(c,kTRUE,format);
  
  // Z rapidity
  sprintf(ylabel,"a.u. / %.2f",hZyv[0]->GetBinWidth(1));
  CPlot plotZy("zy","","y(Z)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZyv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZy.AddHist1D(hZyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZy.TransLegend(0.1,0);
  plotZy.Draw(c,kTRUE,format);
  
  // Z phi
  sprintf(ylabel,"a.u. / %.2f",hZPhiv[0]->GetBinWidth(1));
  CPlot plotZPhi("zphi","","#phi(Z)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hZPhiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotZPhi.AddHist1D(hZPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
//  plotZPhi.SetYRange(0.8*(hZPhiv[0]->GetMaximum()),1.2*(hZPhiv[0]->GetMaximum()));
  plotZPhi.Draw(c,kTRUE,format); 

  
  // dilepton mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDileptonMass1v[0]->GetBinWidth(1));
  CPlot plotDileptonMass1("dileptonmass1","","m(l^{+}l^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDileptonMass1v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDileptonMass1.AddHist1D(hDileptonMass1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDileptonMass1.SetLogy();
  plotDileptonMass1.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hDileptonMass2v[0]->GetBinWidth(1));
  CPlot plotDileptonMass2("dileptonmass2","","m(l^{+}l^{-}) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDileptonMass2v[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDileptonMass2.AddHist1D(hDileptonMass2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDileptonMass2.SetLogy();
  plotDileptonMass2.Draw(c,kTRUE,format);
      
  // dilepton pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hDileptonPtv[0]->GetBinWidth(1));
  CPlot plotDileptonPt("dileptonpt","","p_{T}(l^{+}l^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDileptonPtv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDileptonPt.AddHist1D(hDileptonPtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDileptonPt.Draw(c,kTRUE,format);
  
  // dilepton rapidity
  sprintf(ylabel,"a.u. / %.2f",hDileptonyv[0]->GetBinWidth(1));
  CPlot plotDileptony("dileptony","","y(l^{+}l^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDileptonyv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDileptony.AddHist1D(hDileptonyv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotDileptony.TransLegend(0.1,0);
  plotDileptony.Draw(c,kTRUE,format);
  
  // dilepton phi
  sprintf(ylabel,"a.u. / %.2f",hDileptonPhiv[0]->GetBinWidth(1));
  CPlot plotDileptonPhi("dileptonphi","","#phi(l^{+}l^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hDileptonPhiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotDileptonPhi.AddHist1D(hDileptonPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
//  plotDileptonPhi.SetYRange(0.8*(hDileptonPhiv[0]->GetMaximum()),1.2*(hDileptonPhiv[0]->GetMaximum()));
  plotDileptonPhi.Draw(c,kTRUE,format);
  
  
  // lepton pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hLepton1Ptv[0]->GetBinWidth(1));
  CPlot plotLepton1Pt("mu1pt","","p_{T}(l^{-}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton1Ptv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton1Pt.AddHist1D(hLepton1Ptv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton1Pt.Draw(c,kTRUE,format);
  
  // lepton eta
  sprintf(ylabel,"a.u. / %.2f",hLepton1Etav[0]->GetBinWidth(1)); 
  CPlot plotLepton1Eta("mu1eta","","#eta(l^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton1Etav[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton1Eta.AddHist1D(hLepton1Etav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton1Eta.TransLegend(0.1,0);
  plotLepton1Eta.Draw(c,kTRUE,format);
  
  // lepton phi
  sprintf(ylabel,"a.u. / %.2f",hLepton1Phiv[0]->GetBinWidth(1));
  CPlot plotLepton1Phi("mu1phi","","#phi(l^{-})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton1Phiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton1Phi.AddHist1D(hLepton1Phiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton1Phi.SetYRange(0.8*(hLepton1Phiv[0]->GetMaximum()),1.2*(hLepton1Phiv[0]->GetMaximum()));
  plotLepton1Phi.Draw(c,kTRUE,format);  


  // anti-lepton pT
  sprintf(ylabel,"a.u. / %.1f GeV/c",hLepton2Ptv[0]->GetBinWidth(1));
  CPlot plotLepton2Pt("mu2pt","","p_{T}(l^{+}) [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton2Ptv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton2Pt.AddHist1D(hLepton2Ptv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton2Pt.Draw(c,kTRUE,format);
  
  // anti-lepton eta
  sprintf(ylabel,"a.u. / %.2f",hLepton2Etav[0]->GetBinWidth(1)); 
  CPlot plotLepton2Eta("mu2eta","","#eta(l^{+})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton2Etav[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton2Eta.AddHist1D(hLepton2Etav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton2Eta.TransLegend(0.1,0);
  plotLepton2Eta.Draw(c,kTRUE,format);
  
  // anti-lepton phi
  sprintf(ylabel,"a.u. / %.2f",hLepton2Phiv[0]->GetBinWidth(1));
  CPlot plotLepton2Phi("mu2phi","","#phi(l^{+})",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hLepton2Phiv[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotLepton2Phi.AddHist1D(hLepton2Phiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotLepton2Phi.SetYRange(0.8*(hLepton2Phiv[0]->GetMaximum()),1.2*(hLepton2Phiv[0]->GetMaximum()));
  plotLepton2Phi.Draw(c,kTRUE,format); 
  

  // number of FSR photons
  CPlot plotNPho("npho","","N_{#gamma}","a.u");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hNPhov[i]->Scale(1.0/(Double_t)nEventsv[i]);
    plotNPho.AddHist1D(hNPhov[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotNPho.Draw(c,kTRUE,format);
  
  // leading photon pT (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.1f GeV/c",hPhoPtv[0]->GetBinWidth(1));
  CPlot plotPhoPt("phopt","","leading photon p_{T} [GeV/c]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoPtv[i]->Scale(1.0/(Double_t)hPhoPtv[i]->GetEntries());
    plotPhoPt.AddHist1D(hPhoPtv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoPt.SetLogy();
  plotPhoPt.Draw(c,kTRUE,format);
  
  // leading photon eta (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoEtav[0]->GetBinWidth(1)); 
  CPlot plotPhoEta("phoeta","","leading photon #eta",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoEtav[i]->Scale(1.0/(Double_t)hPhoEtav[i]->GetEntries());
    plotPhoEta.AddHist1D(hPhoEtav[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoEta.TransLegend(0.1,0);
  plotPhoEta.Draw(c,kTRUE,format);
  
  // leading photon phi (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoPhiv[0]->GetBinWidth(1));
  CPlot plotPhoPhi("phophi","","leading photon #phi",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoPhiv[i]->Scale(1.0/(Double_t)hPhoPhiv[i]->GetEntries());
    plotPhoPhi.AddHist1D(hPhoPhiv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoPhi.SetYRange(0.8*(hPhoPhiv[0]->GetMaximum()),1.2*(hPhoPhiv[0]->GetMaximum()));
  plotPhoPhi.Draw(c,kTRUE,format); 
  
  // leading photon-lepton dR (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f",hPhoDRv[0]->GetBinWidth(1));
  CPlot plotPhoDR("phodr","","#DeltaR(#gamma,l)",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoDRv[i]->Scale(1.0/(Double_t)hPhoDRv[i]->GetEntries());
    plotPhoDR.AddHist1D(hPhoDRv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoDR.SetLogy();
  plotPhoDR.Draw(c,kTRUE,format);
  
  // leading photon energy (normalized to number of photons)
  sprintf(ylabel,"a.u. / %.2f GeV",hPhoE1v[0]->GetBinWidth(1)); 
  CPlot plotPhoE1("phoe1","","leading photon E [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE1v[i]->Scale(1.0/(Double_t)hPhoE1v[i]->GetEntries());
    plotPhoE1.AddHist1D(hPhoE1v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE1.SetLogy();
  plotPhoE1.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"a.u. / %.2f MeV",1000.*(hPhoE2v[0]->GetBinWidth(1)));
  CPlot plotPhoE2("phoe2","","leading photon E [GeV]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE2v[i]->Scale(1.0/(Double_t)hPhoE2v[i]->GetEntries());
    plotPhoE2.AddHist1D(hPhoE2v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE2.SetLogy();
  plotPhoE2.Draw(c,kTRUE,format); 

  CPlot plotPhoE3("phoe3","","leading photon E [GeV]","a.u.");
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    hPhoE3v[i]->Scale(1.0/(Double_t)hPhoE3v[i]->GetEntries());
    plotPhoE3.AddHist1D(hPhoE3v[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotPhoE3.SetLogy();
  plotPhoE3.SetLogx();  
  plotPhoE3.Draw(c,kTRUE,format); 


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
  plotScaleRatio1.Draw(c,kTRUE,format);
    
  CPlot plotScaleRatio2("scaleRatio2","","Q [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hScale2v[0],hScale2v[i]);
    plotScaleRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotScaleRatio2.SetYRange(0.5,1.5);
  plotScaleRatio2.AddLine(0,1,220,1,kBlack,2);
  plotScaleRatio2.Draw(c,kTRUE,format);

  // Z mass
  CPlot plotZMassRatio1("zmassRatio1","","m(Z) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZMass1v[0],hZMass1v[i]);
    plotZMassRatio1.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZMassRatio1.SetYRange(0.5,1.5);
  plotZMassRatio1.AddLine(0,1,1100,1,kBlack,2);
  plotZMassRatio1.Draw(c,kTRUE,format);
    
  CPlot plotZMassRatio2("zmassRatio2","","m(Z) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZMass2v[0],hZMass2v[i]);
    plotZMassRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZMassRatio2.SetYRange(0.5,1.5);
  plotZMassRatio2.AddLine(0,1,220,1,kBlack,2);
  plotZMassRatio2.Draw(c,kTRUE,format);
        
  // Z pT
  CPlot plotZPtRatio("zptRatio","","p_{T}(Z) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZPtv[0],hZPtv[i]);
    plotZPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZPtRatio.SetYRange(0.5,1.5);
  plotZPtRatio.AddLine(0,1,110,1,kBlack,2);
  plotZPtRatio.Draw(c,kTRUE,format);
  
  // Z rapidity
  CPlot plotZyRatio("zyRatio","","y(Z)","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hZyv[0],hZyv[i]);
    plotZyRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotZyRatio.SetYRange(0.5,1.5);
  plotZyRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotZyRatio.Draw(c,kTRUE,format);

  // Dilepton mass
  CPlot plotDileptonMassRatio1("dileptonmassRatio1","","m(l^{+}l^{-}) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDileptonMass1v[0],hDileptonMass1v[i]);
    plotDileptonMassRatio1.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDileptonMassRatio1.SetYRange(0.5,1.5);
  plotDileptonMassRatio1.AddLine(0,1,1100,1,kBlack,2);
  plotDileptonMassRatio1.Draw(c,kTRUE,format);
    
  CPlot plotDileptonMassRatio2("dileptonmassRatio2","","m(l^{+}l^{-}) [GeV/c^{2}]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDileptonMass2v[0],hDileptonMass2v[i]);
    plotDileptonMassRatio2.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDileptonMassRatio2.SetYRange(0.5,1.5);
  plotDileptonMassRatio2.AddLine(0,1,220,1,kBlack,2);
  plotDileptonMassRatio2.Draw(c,kTRUE,format);
        
  // Dilepton pT
  CPlot plotDileptonPtRatio("dileptonptRatio","","p_{T}(l^{+}l^{-}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDileptonPtv[0],hDileptonPtv[i]);
    plotDileptonPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDileptonPtRatio.SetYRange(0.5,1.5);
  plotDileptonPtRatio.AddLine(0,1,110,1,kBlack,2);
  plotDileptonPtRatio.Draw(c,kTRUE,format);
  
  // Dilepton rapidity
  CPlot plotDileptonyRatio("dileptonyRatio","","y(l^{+}l^{-})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hDileptonyv[0],hDileptonyv[i]);
    plotDileptonyRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotDileptonyRatio.SetYRange(0.5,1.5);
  plotDileptonyRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotDileptonyRatio.Draw(c,kTRUE,format);
  
  // lepton pT
  CPlot plotLepton1PtRatio("mu1ptRatio","","p_{T}(l^{-}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hLepton1Ptv[0],hLepton1Ptv[i]);
    plotLepton1PtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotLepton1PtRatio.SetYRange(0.5,1.5);
  plotLepton1PtRatio.AddLine(0,1,220,1,kBlack,2);
  plotLepton1PtRatio.Draw(c,kTRUE,format);
  
  // lepton eta 
  CPlot plotLepton1EtaRatio("mu1etaRatio","","#eta(l^{-})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hLepton1Etav[0],hLepton1Etav[i]);
    plotLepton1EtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotLepton1EtaRatio.SetYRange(0.5,1.5);
  plotLepton1EtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotLepton1EtaRatio.Draw(c,kTRUE,format);
  
  // anti-lepton pT
  CPlot plotLepton2PtRatio("mu2ptRatio","","p_{T}(l^{+}) [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hLepton2Ptv[0],hLepton2Ptv[i]);
    plotLepton2PtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotLepton2PtRatio.SetYRange(0.5,1.5);
  plotLepton2PtRatio.AddLine(0,1,220,1,kBlack,2);
  plotLepton2PtRatio.Draw(c,kTRUE,format);
  
  // anti-lepton eta 
  CPlot plotLepton2EtaRatio("mu2etaRatio","","#eta(l^{+})","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hLepton2Etav[0],hLepton2Etav[i]);
    plotLepton2EtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotLepton2EtaRatio.SetYRange(0.5,1.5);
  plotLepton2EtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotLepton2EtaRatio.Draw(c,kTRUE,format);
  
  
  // photon pT
  CPlot plotPhoPtRatio("phoptRatio","","leading photon p_{T} [GeV/c]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoPtv[0],hPhoPtv[i]);
    plotPhoPtRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoPtRatio.SetYRange(0.5,1.5);
  plotPhoPtRatio.AddLine(0,1,22,1,kBlack,2);
  plotPhoPtRatio.Draw(c,kTRUE,format);
  
  // photon eta 
  CPlot plotPhoEtaRatio("phoetaRatio","","leading photon #eta","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoEtav[0],hPhoEtav[i]);
    plotPhoEtaRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoEtaRatio.SetYRange(0.5,1.5);
  plotPhoEtaRatio.AddLine(-9.5,1,9.5,1,kBlack,2);
  plotPhoEtaRatio.Draw(c,kTRUE,format);
  
  // photon-lepton dR
  CPlot plotPhoDRRatio("phodrRatio","","#DeltaR(#gamma,l)","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoDRv[0],hPhoDRv[i]);
    plotPhoDRRatio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoDRRatio.SetYRange(0.5,1.5);
  plotPhoDRRatio.AddLine(0,1,1.1,1,kBlack,2);
  plotPhoDRRatio.Draw(c,kTRUE,format);
  
  // photon energy 
  CPlot plotPhoE1Ratio("phoe1Ratio","","leading photon E [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoE1v[0],hPhoE1v[i]);
    plotPhoE1Ratio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoE1Ratio.SetYRange(0.5,1.5);
  plotPhoE1Ratio.AddLine(0,1,22,1,kBlack,2);
  plotPhoE1Ratio.Draw(c,kTRUE,format);
   
  CPlot plotPhoE2Ratio("phoe2Ratio","","leading photon E [GeV]","ratio");
  for(UInt_t i=1; i<fnamev.size(); i++) { 
    TGraphErrors* ratio = makeRatioGraph(hPhoE2v[0],hPhoE2v[i]);
    plotPhoE2Ratio.AddGraph(ratio,labelv[i],"XL",colorv[i],kFullDotMedium,1); 
  }
  plotPhoE2Ratio.SetYRange(0.5,1.5);
  plotPhoE2Ratio.AddLine(0,1,0.22,1,kBlack,2);
  plotPhoE2Ratio.Draw(c,kTRUE,format);
            
      
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
    cout << "                     Acceptance: " << accv[ifile] << " +/- " << accErrv[ifile] << endl;
    if(ifile>0) {
      Double_t da    = 100.*(accv[ifile]-accv[0])/accv[0];
      Double_t daErr = 100.*sqrt((accErrv[ifile])*(accErrv[ifile])+(accErrv[0])*(accErrv[0]))/accv[0];
      cout << "                          dA/A : (" << da << " +/- " << daErr << ") %" << endl;
    }
    cout << endl;
  }
  cout << endl;
  
  if(kTRUE) {
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
      txtfile << "                     Acceptance: " << accv[ifile] << " +/- " << accErrv[ifile] << endl; 
      if(ifile>0) {
        Double_t da    = 100.*(accv[ifile]-accv[0])/accv[0];
        Double_t daErr = 100.*sqrt((accErrv[ifile])*(accErrv[ifile])+(accErrv[0])*(accErrv[0]))/accv[0];
        txtfile << "                          dA/A : (" << da << " +/- " << daErr << ") %" << endl;
      }
      txtfile << endl;
    }
    txtfile << endl;
    
    cout << " <> Output saved in " << outputDir << "/" << endl;
    cout << endl;
  }
  
  gBenchmark->Show("plotGenZ");
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
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scaleRatio1.png\"><img src=\"gen/scaleRatio1.png\" alt=\"gen/scaleRatio1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/scaleRatio2.png\"><img src=\"gen/scaleRatio2.png\" alt=\"gen/scaleRatio2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
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
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonmass1.png\"><img src=\"gen/dileptonmass1.png\" alt=\"gen/dileptonmass1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonmass2.png\"><img src=\"gen/dileptonmass2.png\" alt=\"gen/dileptonmass2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonmassRatio1.png\"><img src=\"gen/dileptonmassRatio1.png\" alt=\"gen/dileptonmassRatio1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonmassRatio2.png\"><img src=\"gen/dileptonmassRatio2.png\" alt=\"gen/dileptonmassRatio2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonpt.png\"><img src=\"gen/dileptonpt.png\" alt=\"gen/dileptonpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptony.png\"><img src=\"gen/dileptony.png\" alt=\"gen/dileptony.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonphi.png\"><img src=\"gen/dileptonphi.png\" alt=\"gen/dileptonphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonptRatio.png\"><img src=\"gen/dileptonptRatio.png\" alt=\"gen/dileptonptRatio.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"gen/dileptonyRatio.png\"><img src=\"gen/dileptonyRatio.png\" alt=\"gen/dileptonyRatio.png\" width=\"100%\"></a></td>" << endl;
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
