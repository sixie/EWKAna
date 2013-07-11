#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TLegend.h"
#include "THStack.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TFile.h"
#include "TTree.h"
#include "TArrow.h"
#include "TStyle.h"
#include <vector>
#include <math.h>
#include <iostream>
#include <fstream>
#include "TH2F.h"
#include "TF1.h"
#include "Math/LorentzVector.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "Enum.h"
#include "cuts.h"

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void createPlot(std::vector<int> samples, std::vector<TString> files, std::vector<TString> legend);


// WW process for the 6 final states
// 61 '  f(p1)+f(p2) --> W^+(-->nu(p3)+e^+(p4)) +W^-(-->e^-(p5)+nu~(p6))' 'N'


TFile *output_file = new TFile("WW_mcfm_mc.root", "RECREATE");
ofstream text; 

void getmcfmmcvar()
{  
  // load macros  
  gROOT->ProcessLine(".L ~/tdrstyle.C");
  gROOT->ProcessLine("setTDRStyle()"); 
  gROOT->ProcessLine(".L Enum.h++");
  gROOT->ProcessLine(".L cuts.h++");

  // load macros  
  std::vector<int> samples;
  std::vector<TString> files;
  std::vector<TString> legend;

  /*
  samples.push_back(mcfm);
  files.push_back("WWqqbr_tota_mstw8lo_80__80__test.root");
  legend.push_back("MCFM Default Scale");

  samples.push_back(mcfmup);
  files.push_back("WWqqbr_tota_mstw8lo_160_160_test.root");
  legend.push_back("MCFM Scaled Up");

  samples.push_back(mcfmdown);
  files.push_back("WWqqbr_tota_mstw8lo_40__40__test.root");
  legend.push_back("MCFM Scaled Down");
  
  samples.push_back(madgraph);
  files.push_back("/smurf/data/Run2011_Summer11_SmurfV7_42X/mitf-alljets/qqww.root");
  legend.push_back("Madgraph");
  
  samples.push_back(mcnlo);
  files.push_back("/smurf/data/Run2011_Summer11_SmurfV7_42X/mitf-alljets/ww_mcnlo.root");
  legend.push_back("MC@NLO");
  */

  samples.push_back(mcfm);
  files.push_back("WWqqbr_tota_cteq66m_80__80__test.root");
  legend.push_back("MCFM Default Scale");

  samples.push_back(mcfmup);
  files.push_back("WWqqbr_tota_cteq66m_160_160_test.root");
  legend.push_back("MCFM Scaled Up");

  samples.push_back(mcfmdown);
  files.push_back("WWqqbr_tota_cteq66m_40__40__test.root");
  legend.push_back("MCFM Scaled Down");


  samples.push_back(madgraph);
  files.push_back("/smurf/cerati/genww/qqww.root");
  legend.push_back("Madgraph");


  samples.push_back(mcnlo);
  files.push_back("/smurf/cerati/genww/qqww_mcnlo.root");
  legend.push_back("MC@NLO");

  samples.push_back(mcnloup);
  files.push_back("/smurf/cerati/genww/qqww_mcnlo_up_new.root");
  legend.push_back("MC@NLO Scaled Up");


  samples.push_back(mcnlodown);
  files.push_back("/smurf/cerati/genww/qqww_mcnlo_down_new.root");
  legend.push_back("MC@NLO Scaled Down");
  
  
  createPlot(samples, files, legend);
  

}

void createPlot(std::vector<int> samples, std::vector<TString> files, std::vector<TString> legend)
{

  TString y_title = "Number of Entries";
  const int nHist = files.size(); // number of files

  // Declare the histograms to be saved 
  TH1F *hdilmass[nHist];
  TH1F *hdilpt[nHist];
  TH1F *hdphi[nHist];
  TH1F *hleadleppt[nHist];
  TH1F *htrailleppt[nHist];
  TH1F *hmet[nHist];
  TH1F *hmt[nHist];
  TH1F *hdr[nHist];
  
  TH1F *hdilmass_0j[nHist];
  TH1F *hdilpt_0j[nHist];
  TH1F *hdphi_0j[nHist];
  TH1F *hleadleppt_0j[nHist];
  TH1F *htrailleppt_0j[nHist];
  TH1F *hmet_0j[nHist];
  TH1F *hmt_0j[nHist];
  TH1F *hdr_0j[nHist];
  
  TH1F *hmt_sig[nHist];
  TH1F *hmt_bkg[nHist];

  TH1F *hmt_sig_0j[nHist];
  TH1F *hmt_bkg_0j[nHist];

  TH1F *hmt_sig_1j[nHist];
  TH1F *hmt_bkg_1j[nHist];

  
  // Get the histograms from the ntuples
  for (int i=0;i<nHist;i++) {
    TString treeName = "h10";
    if ( samples[i] & MC) treeName = "tree";

    TChain *chain = new TChain(treeName);
    chain->Add(files[i]);
    assert(chain);
  
    // declare histograms  to fill
    Color_t color = kBlack;
    TString sampleName = "mcfm";
    setSample(int(samples[i]), color, sampleName);

    // define the histograms to plot
    
    // dilmass 
    hdilmass[i] = new TH1F(TString("WW_"+sampleName+"_hdilmass"), TString("WW_"+sampleName+"_hdilmass"), 40, 0, 200);
    if ( samples[i] & MC )
      hdilmass[i]->Sumw2();
    hdilmass[i]->SetLineColor(color);
    hdilmass[i]->SetMarkerColor(color);
    
    hdilmass_0j[i] = new TH1F(TString("WW_"+sampleName+"_hdilmass_0j"), TString("WW_"+sampleName+"_hdilmass_0j"), 40, 0, 200);
    if ( samples[i] & MC )
      hdilmass_0j[i]->Sumw2();
    hdilmass_0j[i]->SetLineColor(color);
    hdilmass_0j[i]->SetMarkerColor(color);

    // leading lepton pT
    hleadleppt[i] = new TH1F(TString("WW_"+sampleName+"_hleadleppt"), TString("WW_"+sampleName+"_hleadleppt"), 20, 0, 100);
    if ( samples[i] & MC )
      hleadleppt[i]->Sumw2();
    hleadleppt[i]->SetLineColor(color);
    hleadleppt[i]->SetMarkerColor(color);
    
    hleadleppt_0j[i] = new TH1F(TString("WW_"+sampleName+"_hleadleppt_0j"), TString("WW_"+sampleName+"_hleadleppt_0j"), 20, 0, 100);
    if ( samples[i] & MC )
      hleadleppt_0j[i]->Sumw2();
    hleadleppt_0j[i]->SetLineColor(color);
    hleadleppt_0j[i]->SetMarkerColor(color);

    // trailing lepton pT
    htrailleppt[i] = new TH1F(TString("WW_"+sampleName+"_htrailleppt"), TString("WW_"+sampleName+"_htrailleppt"), 20, 0, 100);
    if ( samples[i] & MC )
      htrailleppt[i]->Sumw2();
    htrailleppt[i]->SetLineColor(color);
    htrailleppt[i]->SetMarkerColor(color);
    
    htrailleppt_0j[i] = new TH1F(TString("WW_"+sampleName+"_htrailleppt_0j"), TString("WW_"+sampleName+"_htrailleppt_0j"), 20, 0, 100);
    if ( samples[i] & MC )
      htrailleppt_0j[i]->Sumw2();
    htrailleppt_0j[i]->SetLineColor(color);
    htrailleppt_0j[i]->SetMarkerColor(color);

    // MET
    hmet[i] = new TH1F(TString("WW_"+sampleName+"_hmet"), TString("WW_"+sampleName+"_hmet"), 40, 0, 100);
    if ( samples[i] & MC )
      hmet[i]->Sumw2();
    hmet[i]->SetLineColor(color);
    hmet[i]->SetMarkerColor(color);
    
    hmet_0j[i] = new TH1F(TString("WW_"+sampleName+"_hmet_0j"), TString("WW_"+sampleName+"_hmet_0j"), 40, 0, 100);
    if ( samples[i] & MC )
      hmet_0j[i]->Sumw2();
    hmet_0j[i]->SetLineColor(color);
    hmet_0j[i]->SetMarkerColor(color);

    // dilepton pT
    hdilpt[i] = new TH1F(TString("WW_"+sampleName+"_hdilpt"), TString("WW_"+sampleName+"_hdilpt"), 20, 0, 100);
    if ( samples[i] & MC )
      hdilpt[i]->Sumw2();
    hdilpt[i]->SetLineColor(color);
    hdilpt[i]->SetMarkerColor(color);
    
    hdilpt_0j[i] = new TH1F(TString("WW_"+sampleName+"_hdilpt_0j"), TString("WW_"+sampleName+"_hdilpt_0j"), 20, 0, 100);
    if ( samples[i] & MC )
      hdilpt_0j[i]->Sumw2();
    hdilpt_0j[i]->SetLineColor(color);
    hdilpt_0j[i]->SetMarkerColor(color);

    // deltaphi (ll)
    hdphi[i] = new TH1F(TString("WW_"+sampleName+"_hdphi"), TString("WW_"+sampleName+"_hdphi"), 18, 0, 180.0);
    if ( samples[i] & MC )
      hdphi[i]->Sumw2();
    hdphi[i]->SetLineColor(color);
    hdphi[i]->SetMarkerColor(color);
    
    hdphi_0j[i] = new TH1F(TString("WW_"+sampleName+"_hdphi_0j"), TString("WW_"+sampleName+"_hdphi_0j"), 18, 0, 180);
    if ( samples[i] & MC )
      hdphi_0j[i]->Sumw2();
    hdphi_0j[i]->SetLineColor(color);
    hdphi_0j[i]->SetMarkerColor(color);

    // deltaR
    hdr[i] = new TH1F(TString("WW_"+sampleName+"_hdr"), TString("WW_"+sampleName+"_hdr"), 25, 0, 5);
    if ( samples[i] & MC )
      hdr[i]->Sumw2();
    hdr[i]->SetLineColor(color);
    hdr[i]->SetMarkerColor(color);
    
    hdr_0j[i] = new TH1F(TString("WW_"+sampleName+"_hdr_0j"), TString("WW_"+sampleName+"_hdr_0j"), 25, 0, 5);
    if ( samples[i] & MC )
      hdr_0j[i]->Sumw2();
    hdr_0j[i]->SetLineColor(color);
    hdr_0j[i]->SetMarkerColor(color);
    
    // transverse mass
    hmt[i] = new TH1F(TString("WW_"+sampleName+"_hmt"), TString("WW_"+sampleName+"_hmt"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt[i]->Sumw2();
    hmt[i]->SetLineColor(color);
    hmt[i]->SetMarkerColor(color);
    
    hmt_0j[i] = new TH1F(TString("WW_"+sampleName+"_hmt_0j"), TString("WW_"+sampleName+"_hmt_0j"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_0j[i]->Sumw2();
    hmt_0j[i]->SetLineColor(color);
    hmt_0j[i]->SetMarkerColor(color);
    
    // transverse mass
    hmt_sig[i] = new TH1F(TString("WW_"+sampleName+"_hmt_sig"), TString("WW_"+sampleName+"_hmt_sig"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_sig[i]->Sumw2();
    hmt_sig[i]->SetLineColor(color);
    hmt_sig[i]->SetMarkerColor(color);

    // transverse mass
    hmt_sig_0j[i] = new TH1F(TString("WW_"+sampleName+"_hmt_sig_0j"), TString("WW_"+sampleName+"_hmt_sig_0j"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_sig_0j[i]->Sumw2();
    hmt_sig_0j[i]->SetLineColor(color);
    hmt_sig_0j[i]->SetMarkerColor(color);

    // transverse mass
    hmt_sig_1j[i] = new TH1F(TString("WW_"+sampleName+"_hmt_sig_1j"), TString("WW_"+sampleName+"_hmt_sig_1j"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_sig_1j[i]->Sumw2();
    hmt_sig_1j[i]->SetLineColor(color);
    hmt_sig_1j[i]->SetMarkerColor(color);
    
    // transverse mass
    hmt_bkg[i] = new TH1F(TString("WW_"+sampleName+"_hmt_bkg"), TString("WW_"+sampleName+"_hmt_bkg"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_bkg[i]->Sumw2();
    hmt_bkg[i]->SetLineColor(color);
    hmt_bkg[i]->SetMarkerColor(color);

    // transverse mass
    hmt_bkg_0j[i] = new TH1F(TString("WW_"+sampleName+"_hmt_bkg_0j"), TString("WW_"+sampleName+"_hmt_bkg_0j"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_bkg_0j[i]->Sumw2();
    hmt_bkg_0j[i]->SetLineColor(color);
    hmt_bkg_0j[i]->SetMarkerColor(color);

    // transverse mass
    hmt_bkg_1j[i] = new TH1F(TString("WW_"+sampleName+"_hmt_bkg_1j"), TString("WW_"+sampleName+"_hmt_bkg_1j"), 50, 0, 250);
    if ( samples[i] & MC )
      hmt_bkg_1j[i]->Sumw2();
    hmt_bkg_1j[i]->SetLineColor(color);
    hmt_bkg_1j[i]->SetMarkerColor(color);
    
    std::cout  << "Processing " << chain->GetEntries() << " entries. \n";
    
    // mcfm variables to be used
    float px3_ = 1.0;
    float py3_ = 1.0;
    float pz3_ = 1.0;
    float E3_ = 1.0;
    float px4_ = 1.0;
    float py4_ = 1.0;
    float pz4_ = 1.0;
    float E4_ = 1.0;
    float px5_ = 1.0;
    float py5_ = 1.0;
    float pz5_ = 1.0;
    float E5_ = 1.0;
    float px6_ = 1.0;
    float py6_ = 1.0;
    float pz6_ = 1.0;
    float E6_ = 1.0;
    float px7_ = 1.0;
    float py7_ = 1.0;
    float pz7_ = 1.0;
    float E7_ = 1.0;
    float wt_ALL_ = 1.0;

    if (chain->GetBranchStatus("px3"))
      chain->SetBranchAddress("px3", &px3_);
    if (chain->GetBranchStatus("py3"))
      chain->SetBranchAddress("py3", &py3_);
    if (chain->GetBranchStatus("pz3"))
      chain->SetBranchAddress("pz3", &pz3_);
    if (chain->GetBranchStatus("E_3"))
      chain->SetBranchAddress("E_3", &E3_);

    if (chain->GetBranchStatus("px4"))
      chain->SetBranchAddress("px4", &px4_);
    if (chain->GetBranchStatus("py4"))
      chain->SetBranchAddress("py4", &py4_);
    if (chain->GetBranchStatus("pz4"))
      chain->SetBranchAddress("pz4", &pz4_);
    if (chain->GetBranchStatus("E_4"))
      chain->SetBranchAddress("E_4", &E4_);

    if (chain->GetBranchStatus("px5"))
      chain->SetBranchAddress("px5", &px5_);
    if (chain->GetBranchStatus("py5"))
      chain->SetBranchAddress("py5", &py5_);
    if (chain->GetBranchStatus("pz5"))
      chain->SetBranchAddress("pz5", &pz5_);
    if (chain->GetBranchStatus("E_5"))
      chain->SetBranchAddress("E_5", &E5_);

    if (chain->GetBranchStatus("px6"))
      chain->SetBranchAddress("px6", &px6_);
    if (chain->GetBranchStatus("py6"))
      chain->SetBranchAddress("py6", &py6_);
    if (chain->GetBranchStatus("pz6"))
      chain->SetBranchAddress("pz6", &pz6_);
    if (chain->GetBranchStatus("E_6"))
      chain->SetBranchAddress("E_6", &E6_);

    if (chain->GetBranchStatus("px7"))
      chain->SetBranchAddress("px7", &px7_);
    if (chain->GetBranchStatus("py7"))
      chain->SetBranchAddress("py7", &py7_);
    if (chain->GetBranchStatus("pz7"))
      chain->SetBranchAddress("pz7", &pz7_);
    if (chain->GetBranchStatus("E_7"))
      chain->SetBranchAddress("E_7", &E7_);

    if (chain->GetBranchStatus("wt_ALL"))
      chain->SetBranchAddress("wt_ALL", &wt_ALL_);
      
    // smurf tree variables to be used
    LorentzVector*  jet1_ = 0;  
    float scale1fb_ = 1.;
    int type_ = 0;
    LorentzVector*  dilep_ = 0;
    LorentzVector*  lep1_ = 0;
    LorentzVector*  lep2_ = 0;
    unsigned int cuts_ = 0;
    float pmet_ = 0;
    float pTrackMet_ = 0;
    float mt_ = 0;
    float met_ = 0;
    int nvtx_ = 0;
    int njets_ = 0;
    float dPhi_ = 0;
    float dR_ = 0;
    int lep1MotherMcId_ = 0;
    int lep2MotherMcId_ = 0;
    
    
    if ( chain->GetBranchStatus("jet1") )
      chain->SetBranchAddress( "jet1", &jet1_);
    if ( chain->GetBranchStatus("scale1fb") )
      chain->SetBranchAddress( "scale1fb",  &scale1fb_);
    if ( chain->GetBranchStatus("type"))
      chain->SetBranchAddress( "type",  &type_);
    if ( chain->GetBranchStatus("dilep"))
      chain->SetBranchAddress( "dilep", &dilep_);
    if ( chain->GetBranchStatus("lep1"))
      chain->SetBranchAddress( "lep1", &lep1_);
    if ( chain->GetBranchStatus("lep2"))
      chain->SetBranchAddress( "lep2", &lep2_);
    if ( chain->GetBranchStatus("cuts"))
      chain->SetBranchAddress( "cuts"      , &cuts_     );     
    if ( chain->GetBranchStatus("met"))
      chain->SetBranchAddress( "met"      , &met_     );     
    if ( chain->GetBranchStatus("pmet"))
      chain->SetBranchAddress( "pmet"      , &pmet_     );     
    if ( chain->GetBranchStatus("pTrackMet"))
      chain->SetBranchAddress( "pTrackMet"      , &pTrackMet_     );     
    if ( chain->GetBranchStatus("nvtx"))
      chain->SetBranchAddress( "nvtx"      , &nvtx_     );     
    if ( chain->GetBranchStatus("njets"))
      chain->SetBranchAddress( "njets"      , &njets_     );     
    if ( chain->GetBranchStatus("dPhi"))
      chain->SetBranchAddress( "dPhi"      , &dPhi_     );    
    if ( chain->GetBranchStatus("dR"))
      chain->SetBranchAddress( "dR"      , &dR_     );    
    if ( chain->GetBranchStatus("mt"))
      chain->SetBranchAddress( "mt"      , &mt_     );     
    if ( chain->GetBranchStatus("lep1MotherMcId"))
      chain->SetBranchAddress( "lep1MotherMcId"      , &lep1MotherMcId_     );      
    if ( chain->GetBranchStatus("lep2MotherMcId"))
      chain->SetBranchAddress( "lep2MotherMcId"      , &lep2MotherMcId_     );      
    
      
    for ( int ievt = 0; ievt < chain->GetEntries(); ievt++) {
      chain->GetEntry(ievt);
      
      // variables to fill
      float dilmsq(0.0); 
      float dilmass(0.0);
      float leadleppt(0.);
      float trailleppt(0.);
      float met(0.);
      float dilpt(0.);
      float dphi(0.);
      float mt(0.);
      float dr(0.);
      
      // variables to use
      int njets(999);
      float weight(1.0);
      
      
      // Fill the variables respectively for MCFM and MC trees
      // assess smurf based kinematics
      if ( samples[i] & MC) {      
	// use mctruth
	if ( abs(lep1MotherMcId_) != 24 || abs(lep2MotherMcId_) != 24 ) continue;
	leadleppt = lep1_->Pt();
	trailleppt = lep2_->Pt();
	// acceptance cuts
	if ( TMath::Abs(lep1_->Eta()) > 2.5 ) continue;
	if ( TMath::Abs(lep2_->Eta()) > 2.5 ) continue;
	dilmass = dilep_->mass();
	met = met_;
	dilpt = dilep_->Pt();
	dphi = dPhi_;
	dr = dR_;
	mt = mt_;
	weight = scale1fb_;
	njets = njets_;

      }
      // assess MCFM based kinematics
      else {
	// initialize the lorentz vectors
	LorentzVector p3(px3_, py3_, pz3_, E3_);
	LorentzVector p4(px4_, py4_, pz4_, E4_);
	LorentzVector p5(px5_, py5_, pz5_, E5_);
	LorentzVector p6(px6_, py6_, pz6_, E6_);
	// std::cout << Form("p = (%.1f, %.1f, %.1f, %.1f): pT is %.1f ", px3_, py3_, pz3_, E3_, p1.Pt()) ;
	
	leadleppt = p4.Pt();
	trailleppt = p5.Pt();
	if ( p4.Pt() < p5.Pt() ) {
	  leadleppt = p5.Pt();
	  trailleppt = p4.Pt();
	}
	if ( TMath::Abs(p4.Eta()) > 2.5 ) continue;
	if ( TMath::Abs(p5.Eta()) > 2.5 ) continue;
	dilpt = (p4+p5).Pt();
	met = sqrt(pow(px3_+px6_,2) + pow(py3_+py6_,2));
	dilmsq = (p4+p5).mass2();
	dilmass = dilmsq > 0 ? (p4+p5).mass() : 0;

	// deltaphi
	dphi =TMath::Abs(ROOT::Math::VectorUtil::DeltaPhi(p4, p5));
	dr =  ROOT::Math::VectorUtil::DeltaR(p4, p5); 
	float dphidilmet = ROOT::Math::VectorUtil::DeltaPhi(p4+p5, p3+p6); 
	// mt
	mt = 2*sqrt(dilpt * met ) * fabs(sin(dphidilmet / 2.));
	weight = wt_ALL_;
	if ( (p3+p4+p5+p6).Pt() < 30 ) njets = 0;
      }

      if ( dr < 0.3) continue;

      // if ( pass_sig_atlas(leadleppt, trailleppt, met, dilmass, dphi, dilpt) ) {
      if ( pass_sig_cms(leadleppt, trailleppt, met, dilmass, dphi, mt, dilpt) ) {
	hmt_sig[i]->Fill(mt, weight);
	if ( njets == 0 )
	  hmt_sig_0j[i]->Fill(mt, weight);
	if ( samples[i] & MC && njets == 1)
	  hmt_sig_1j[i]->Fill(mt, weight); 
	  
      }
      
      //if ( pass_bkg_atlas(leadleppt, trailleppt, met, dilmass, dilpt) ) {
      if ( pass_bkg_cms(leadleppt, trailleppt, met, dilmass, dilpt) ) {
	hmt_bkg[i]->Fill(mt, weight);	
	if ( njets == 0 )
	  hmt_bkg_0j[i]->Fill(mt, weight);
	if ( samples[i] & MC && njets == 1)
	  hmt_bkg_1j[i]->Fill(mt, weight); 
      }
      
      // fill plots at the baseline value
      // if ( pass_baseline_atlas(leadleppt, trailleppt, met, dilmass) ) {
      //if ( pass_baseline_cms_ww(leadleppt, trailleppt, met, dilmass) ) {
      if ( pass_baseline_cms_hww(leadleppt, trailleppt, met, dilmass) ) {
	hdilmass[i]->Fill(dilmass, weight);
	hdilpt[i]->Fill(dilpt, weight);
	hdphi[i]->Fill(dphi*180./TMath::Pi(), weight);
	hleadleppt[i]->Fill(leadleppt, weight);
	htrailleppt[i]->Fill(trailleppt, weight);
	hmt[i]->Fill(mt, weight);
	hmet[i]->Fill(met, weight);
	hdr[i]->Fill(dr, weight);
	
	if ( njets == 0 ) {
	  hdilmass_0j[i]->Fill(dilmass, weight);
	  hdilpt_0j[i]->Fill(dilpt, weight);
	  hdphi_0j[i]->Fill(dphi*180./TMath::Pi(), weight);
	  hleadleppt_0j[i]->Fill(leadleppt, weight);
	  htrailleppt_0j[i]->Fill(trailleppt, weight);
	  hmt_0j[i]->Fill(mt, weight);
	  hmet_0j[i]->Fill(met, weight);
	  hdr_0j[i]->Fill(dr, weight);
	}
      }
      
    }
  }
  
  output_file->cd();
  
  for(int i=0;i<nHist;i++) {
    hdilmass[i]->Write();
    hdilpt[i]->Write();
    hdphi[i]->Write();
    hleadleppt[i]->Write();
    htrailleppt[i]->Write();
    hmet[i]->Write();
    hmt[i]->Write();
    hdphi[i]->Write();
    hdr[i]->Write();
    hmt_sig[i]->Write();
    hmt_bkg[i]->Write();
    hmt_sig_0j[i]->Write();
    hmt_bkg_0j[i]->Write();
    hmt_sig_1j[i]->Write();
    hmt_bkg_1j[i]->Write();

    hdilmass_0j[i]->Write();
    hdilpt_0j[i]->Write();
    hdphi_0j[i]->Write();
    hleadleppt_0j[i]->Write();
    htrailleppt_0j[i]->Write();
    hmet_0j[i]->Write();
    hmt_0j[i]->Write();
    hdr_0j[i]->Write();
  }
  
  // tidy up
  
  for ( int i = 0; i<nHist;i++) {
    delete hdilmass[i];
    delete hdilpt[i];
    delete hdphi[i];
    delete hleadleppt[i];
    delete htrailleppt[i];
    delete hmet[i];
    delete hmt[i];
    delete hdr[i];

    delete hdilmass_0j[i];
    delete hdilpt_0j[i];
    delete hdphi_0j[i];
    delete hleadleppt_0j[i];
    delete htrailleppt_0j[i];
    delete hmet_0j[i];
    delete hmt_0j[i];
    delete hdr_0j[i];
  }

}

