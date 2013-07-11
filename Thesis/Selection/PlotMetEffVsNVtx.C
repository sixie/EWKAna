#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TGraphErrors.h>           // Graph class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <TH1F.h>                   // ROOT math functions
#include <TGraphAsymmErrors.h>      // ROOT math functions
#include <TLegend.h>      // ROOT math functions

// #include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "EWKAna/Thesis/Background/DY/SmurfTree.h"
#endif

static const unsigned int basic_selection = SmurfTree::BaseLine | 
                                            SmurfTree::ChargeMatch | 
                                            SmurfTree::Lep1FullSelection |
                                            SmurfTree::Lep2FullSelection |
                                            SmurfTree::ExtraLeptonVeto |
                                            SmurfTree::TopVeto;


// compute projected MET
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi) 
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = acos(cos(lepPhi-metPhi));
  if(dphi > 0.5*pi)
    return met;
    
  return met*sin(dphi);
}


//=== MAIN MACRO =================================================================================================
// Options:
//
// --- ZWindowSubtractionMethod == 0 : Opposite Flavor Background Subtraction
// --- ZWindowSubtractionMethod == 1 : Data Corrected MC Background Subtraction
// 
// --- MassCutLow  : lower  mass cut with respect to z pole mass
// --- MassCutHigh : higher mass cut with respect to z pole mass
//

void SaveMetEffVsNVtx(Int_t period = 13, Int_t ZWindowSubtractionMethod = 0, 
                      Double_t MassCutLow = 15, Double_t MassCutHigh = 15)
{

  //*******************************************************
  // Settings 
  //*******************************************************

  Double_t lumi = 1;
  TString filesPath   = "dummy";
  unsigned int minRun = 0;
  unsigned int maxRun = 999999;

  //*******************************************************
  //Input Files
  //*******************************************************
  vector<TString> infilenamev;  

  infilenamev.push_back("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/dymm.root");
  infilenamev.push_back("/data/smurf/sixie/data/Run2011_Fall11_MVAIDIsoCombinedDetIsoSameSigWP/mitf-alljets/dyee.root");

  //*******************************************************
  //mH dependent Cuts
  //*******************************************************
  const Double_t mZ = 91.1876;
  
  const Int_t nbins = 4;  
  const Float_t bins[nbins+1] = {20, 25, 30, 37, 50};  

  vector<Double_t> binEdges;
  for (UInt_t k=0; k<nbins+1; ++k) binEdges.push_back(bins[k]);
  
  const Int_t nmass = 19;
  const Double_t mH[nmass] = {0,115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};  
    
  //*******************************************************
  //Yields and  histograms
  //*******************************************************
  TH1F *NVtx = new TH1F("NVtx","; Number of Primary Vertices; Number of Events", 30,-0.5,29.5);
  TH1F *NVtx_MinProjMetCut40 = new TH1F("NVtx_ProjMetCut40","; Number of Primary Vertices; Number of Events", 30,-0.5,29.5);
  TH1F *NVtx_MinProjMetCutNVtxDependent = new TH1F("NVtx_ProjMetCutNVtxDependent","; Number of Primary Vertices; Number of Events", 30,-0.5,29.5);


  //*******************************************************
  //Systematic Error on VZ Normalization
  //*******************************************************
  const Double_t vzNormSystematic = 0.10; //10% systematic on the VZ cross section inside Z window

  //*******************************************************
  //Event Loop
  //*******************************************************
  SmurfTree tree;
  for(UInt_t ifile=0; ifile<infilenamev.size(); ifile++) {
    cout << "Processing " << infilenamev[ifile] << "..." << endl;
    tree.LoadTree(infilenamev[ifile]);
    tree.InitTree(0);
    
    for(UInt_t ientry = 0; ientry <tree.tree_->GetEntries(); ientry++){
      tree.tree_->GetEntry(ientry);
        
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      if(tree.dstype_ == SmurfTree::data && tree.run_ <  minRun) continue;
      if(tree.dstype_ == SmurfTree::data && tree.run_ >  maxRun) continue;

      //apply trigger requirement
      if(tree.dstype_ == SmurfTree::data)
        if ( (tree.cuts_ & SmurfTree::Trigger) != SmurfTree::Trigger ) continue; 
        
      if(!((tree.cuts_ & SmurfTree::BaseLine) == SmurfTree::BaseLine)            ) continue;
      if(!((tree.cuts_ & SmurfTree::ChargeMatch) == SmurfTree::ChargeMatch)      ) continue;
      if( (tree.cuts_ & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
      if( !((tree.cuts_ & SmurfTree::TopVeto) == SmurfTree::TopVeto)             ) continue; // cut on btagging
      if(tree.dilep_.M() < 20) continue;
        
      Int_t ijet = tree.njets_;
      if(ijet >= 2){
        if(tree.jet3_.Pt() <= 30)					                       ijet = 2;
        else if(tree.jet3_.Pt() > 30 && (
    	  (tree.jet1_.Eta()-tree.jet3_.Eta() > 0 && tree.jet2_.Eta()-tree.jet3_.Eta() < 0) ||
    	  (tree.jet2_.Eta()-tree.jet3_.Eta() > 0 && tree.jet1_.Eta()-tree.jet3_.Eta() < 0)))   ijet = 3;
        else							                               ijet = 2;
        if(tree.njets_ < 2 || tree.njets_ > 3)                                                 ijet = 3;

	if(TMath::Abs(tree.jet1_.Eta()) >= 4.5 ||TMath::Abs(tree.jet2_.Eta()) >= 4.5)          ijet = 3;
      }
      if(ijet>2) continue;

      if(tree.lep1_.Pt() < 20) continue;
      if(tree.lep2_.Pt() < 20) continue;
 
      if (!(tree.type_==SmurfTree::mm || tree.type_==SmurfTree::ee)) continue;
      if(mZ - tree.dilep_.M() < MassCutLow && tree.dilep_.M() - mZ < MassCutHigh) continue;

      Double_t pfmet     = tree.met_;
      Double_t pfmetphi  = tree.metPhi_;
      Double_t trkmet    = tree.trackMet_;
      Double_t trkmetphi = tree.trackMetPhi_;
      Double_t minpfmet  = TMath::Min(projectedMET(pfmet,pfmetphi,tree.lep1_.Phi()),  projectedMET(pfmet,pfmetphi,tree.lep2_.Phi()));
      Double_t mintrkmet = TMath::Min(projectedMET(trkmet,trkmetphi,tree.lep1_.Phi()),projectedMET(trkmet,trkmetphi,tree.lep2_.Phi()));
      Double_t minmet    = TMath::Min(minpfmet,mintrkmet);

      NVtx->Fill(tree.nvtx_);
      if (minmet > 40) NVtx_MinProjMetCut40->Fill(tree.nvtx_);
      if(minmet > 37 + 0.5 * tree.nvtx_ ) NVtx_MinProjMetCutNVtxDependent->Fill(tree.nvtx_);
     
 
    } //loop over events
  } // loop over input files 

   //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  vector<double> nvtxbins;
  for (UInt_t k = 0; k < 15; ++k) nvtxbins.push_back(k*2);

  TGraphAsymmErrors *Eff_MinProjMetCut40 = mithep::EfficiencyUtils::createEfficiencyGraph(NVtx_MinProjMetCut40, NVtx, "MinProjMetCut40Efficiency_NVtx", nvtxbins, 2, -99, -99, 0, 1);
  TGraphAsymmErrors *Eff_MinProjMetCutNVtxDependent = mithep::EfficiencyUtils::createEfficiencyGraph(NVtx_MinProjMetCutNVtxDependent, NVtx, "MinProjMetCutNVtxDependentEfficiency_NVtx", nvtxbins, 2, -99, -99, 0, 1);


  TFile *file = new TFile("MetEfficiency.root", "UPDATE");
  file->cd();
  file->WriteTObject(NVtx, NVtx->GetName(), "WriteDelete");
  file->WriteTObject(NVtx_MinProjMetCut40, NVtx_MinProjMetCut40->GetName(), "WriteDelete");
  file->WriteTObject(NVtx_MinProjMetCutNVtxDependent, NVtx_MinProjMetCutNVtxDependent->GetName(), "WriteDelete");
  file->WriteTObject(Eff_MinProjMetCut40, Eff_MinProjMetCut40->GetName(), "WriteDelete");
  file->WriteTObject(Eff_MinProjMetCutNVtxDependent, Eff_MinProjMetCutNVtxDependent->GetName(), "WriteDelete");
  file->Close();
  delete file;

}


void MakePlot() {

  TFile *file = new TFile("MetEfficiency.root", "READ");
  TGraphAsymmErrors *Eff_MinProjMetCut40 = (TGraphAsymmErrors*)file->Get("MinProjMetCut40Efficiency_NVtx");
  TGraphAsymmErrors *Eff_MinProjMetCutNVtxDependent = (TGraphAsymmErrors*)file->Get("MinProjMetCutNVtxDependentEfficiency_NVtx");
  assert(Eff_MinProjMetCut40);
  assert(Eff_MinProjMetCutNVtxDependent);

  Eff_MinProjMetCut40->SetLineColor(kRed);
  Eff_MinProjMetCut40->SetMarkerColor(kRed);

  TLegend *legend = 0;
  legend = new TLegend(0.2,0.75,0.95,0.90);   
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(Eff_MinProjMetCut40, "Min ProjMET > 40 GeV/c" , "LP");
  legend->AddEntry(Eff_MinProjMetCutNVtxDependent, "Min ProjMET > (37 + 0.5 #times NPV) GeV/c" , "LP");
 

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  Eff_MinProjMetCut40->SetTitle("");
  Eff_MinProjMetCut40->Draw("AP");
  Eff_MinProjMetCutNVtxDependent->Draw("Psame");
  Eff_MinProjMetCut40->GetXaxis()->SetTitle("Number of Primary Vertices (NPV)");
  Eff_MinProjMetCut40->GetXaxis()->SetRangeUser(0.5,20.5);
  Eff_MinProjMetCut40->GetYaxis()->SetRangeUser(0.0,0.004);
  Eff_MinProjMetCut40->GetXaxis()->SetTitleOffset(1.05);
  Eff_MinProjMetCut40->GetYaxis()->SetLabelSize(0.04);
  legend->Draw();
  cv->SaveAs("MinMetEfficiencyVsNVtx.gif");
  cv->SaveAs("MinMetEfficiencyVsNVtx.eps");

}


void PlotMetEffVsNVtx() {

//   SaveMetEffVsNVtx();
  MakePlot();

}

