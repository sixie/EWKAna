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

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void plotMCFM() {

  // Declare the histograms to be saved 
  TH1F *hdilmass = new TH1F( "WW_hdilmass", "WW_hdilmass", 40, 0, 200);
  TH1F *hdilpt = new TH1F("WW_hdilpt", "WW_hdilpt", 20, 0, 100);
  TH1F *hdphi  = new TH1F("WW_hdphi", "WW_hdphi", 18, 0, 180.0);
  TH1F *hleadleppt = new TH1F("WW_hleadleppt", "WW_hleadleppt", 20, 0, 100);
  TH1F *htrailleppt= new TH1F("WW_htrailleppt", "WW_htrailleppt", 20, 0, 100);
  TH1F *hmet = new TH1F("WW_hmet", "WW_hmet", 40, 0, 100);
  TH1F *hmt = new TH1F("WW_hmt", "WW_hmt", 50, 0, 250);
  TH1F *hdr = new TH1F("WW_hdr", "WW_hdr", 25, 0, 5);


  TChain *chain = new TChain("h10");
  chain->Add("/data/smurf/yygao/data/MCFM/WWqqbr_tota_mstw8lo_80__80__test.root");
  assert(chain);

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
      

  for ( int ievt = 0; ievt < chain->GetEntries(); ievt++) {
    chain->GetEntry(ievt);
 
    // initialize the lorentz vectors
    LorentzVector p3(px3_, py3_, pz3_, E3_);
    LorentzVector p4(px4_, py4_, pz4_, E4_);
    LorentzVector p5(px5_, py5_, pz5_, E5_);
    LorentzVector p6(px6_, py6_, pz6_, E6_);

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

    if ( dr < 0.3) continue;
      
    hdilmass->Fill(dilmass, weight);
    hdilpt->Fill(dilpt, weight);
    hdphi->Fill(dphi*180./TMath::Pi(), weight);
    hleadleppt->Fill(leadleppt, weight);
    htrailleppt->Fill(trailleppt, weight);
    hmt->Fill(mt, weight);
    hmet->Fill(met, weight);
    hdr->Fill(dr, weight);
  }

    
  TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
  file->cd();
  file->WriteTObject(hdilmass, hdilmass->GetName(), "WriteDelete");
  file->WriteTObject(hdilpt, hdilpt->GetName(), "WriteDelete");
  file->WriteTObject(hdphi, hdphi->GetName(), "WriteDelete");
  file->WriteTObject(hleadleppt, hleadleppt->GetName(), "WriteDelete");
  file->WriteTObject(htrailleppt, htrailleppt->GetName(), "WriteDelete");
  file->WriteTObject(hmet, hmet->GetName(), "WriteDelete");
  file->WriteTObject(hmt, hmt->GetName(), "WriteDelete");
  file->WriteTObject(hdr, hdr->GetName(), "WriteDelete");
  file->Close();    


}
