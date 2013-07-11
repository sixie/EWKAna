//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <MitStyle.h>
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"
#include "TPaveLabel.h"

// // lumi section selection with JSON files
// #include "MitCommon/DataFormats/interface/Types.h"
// #include "MitAna/DataCont/interface/RunLumiRangeMap.h"
// #include "MitCommon/MathTools/interface/MathUtils.h"
// #include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
// #include "MitHiggs/Utils/interface/EfficiencyUtils.h"
// #include "MitHiggs/Utils/interface/PlotUtils.h"

// helper functions for lepton ID selection
// #include "EWKAna/Utils/LeptonIDCuts.hh"
// #include "MitPhysics/Utils/interface/ElectronIDMVA.h"
// #include "TMVA/Tools.h"
// #include "TMVA/Reader.h"

#endif

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}



//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


//*************************************************************************************************
//Make Projection Graph
//*************************************************************************************************
TGraphAsymmErrors *MakeMeanVsNPUGraph(TH2F *hist2D, string name, string axisLabel = "") {

  const UInt_t nPoints = hist2D->GetYaxis()->GetNbins();
  double NPU[nPoints];
  double NPUErr[nPoints];
  double MeanValue[nPoints];
  double MeanErr[nPoints];

//   cout << "Npoints: " << nPoints << endl;

  TFile *file = new TFile("Test.root", "UPDATE");
  
  for(UInt_t b=0; b < nPoints; ++b) {
    TH1F *Projection = (TH1F*)hist2D->ProjectionX("_px", b,b+1);
    NPU[b] = (hist2D->GetYaxis()->GetBinUpEdge(b) + hist2D->GetYaxis()->GetBinLowEdge(b)) / 2;
    NPUErr[b] = 0.1;
    MeanValue[b] = Projection->GetMean();
    MeanErr[b] = Projection->GetMeanError();

//     cout << b << " : " << NPU[b] << " : " << MeanValue[b] << " " << MeanErr[b] << endl;
    char buffer[200];
    sprintf(buffer,"%d",b);
    if (name == "Ele_HoverE_VS_NPU_Fall11ZeeMC_Graph" || name == "Ele_HcalDepth1OverEcal_VS_NPU_Fall11ZeeMC_Graph" ||name == "Ele_HcalDepth2OverEcal_VS_NPU_Fall11ZeeMC_Graph") {
      file->WriteTObject(Projection , (name + "_Proj" + string(buffer)).c_str(), "WriteDelete");  
    }
  }
  file->Close();  

  TGraphAsymmErrors *Graph = new TGraphAsymmErrors (nPoints,  NPU, MeanValue, NPUErr, NPUErr, MeanErr, MeanErr);
  Graph->SetName(name.c_str());
  Graph->SetTitle("");
  Graph->SetMarkerColor(kBlack);
  Graph->GetXaxis()->SetTitleOffset(1.02);
  Graph->GetXaxis()->SetTitle("Number of Pileup Events");
  Graph->GetYaxis()->SetTitleOffset(1.05);
  Graph->GetYaxis()->SetTitle(axisLabel.c_str());

  return Graph;

}



void ComputeElectronEffectiveArea(string Label, Bool_t isMC = kTRUE, Int_t EtaBin = -1) {

  string label = "";
  if (Label != "") label = "_" + Label;
  
  if (EtaBin == 0) {
    label = label + "_Eta0To1";
  }
  if (EtaBin == 1) {
    label = label + "_Eta1To1p479";
  }
  if (EtaBin == 2) {
    label = label + "_Eta1p479To2";
  }
  if (EtaBin == 3) {
    label = label + "_Eta2To2p25";
  }
  if (EtaBin == 4) {
    label = label + "_Eta2p25To2p5";
  }

  TFile *file = new TFile("EffectiveArea.root", "READ");
    
  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HoverE_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_EMIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HadIso03_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_EMIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_HadIso04_Vs_NPU_Graph;
  TGraphAsymmErrors *Ele_Rho_Vs_NPU_Graph;
  if (isMC) {
    Ele_ChargedIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Graph").c_str());
    Ele_ChargedIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Graph").c_str());
    Ele_NeutralHadronIso007_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Graph").c_str());
    Ele_HoverE_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HoverE_VS_NPU"+label+"_Graph").c_str());
    Ele_HcalDepth1OverEcal_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Graph").c_str());
    Ele_HcalDepth2OverEcal_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Graph").c_str());
    Ele_TrkIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_EMIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_HadIso03_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso03_VS_NPU"+label+"_Graph").c_str());
    Ele_TrkIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_EMIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_HadIso04_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso04_VS_NPU"+label+"_Graph").c_str());
    Ele_Rho_Vs_NPU_Graph = (TGraphAsymmErrors*)file->Get(("Ele_Rho_VS_NPU"+label+"_Graph").c_str());
  }

  TGraphAsymmErrors *Ele_ChargedIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_ChargedIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_ChargedIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_NeutralHadronIso007_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HoverE_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HoverE_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HcalDepth1OverEcal_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HcalDepth2OverEcal_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_TrkIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_EMIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HadIso03_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso03_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_TrkIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_TrkIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_EMIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_EMIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_HadIso04_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_HadIso04_VS_NVtx"+label+"_Graph").c_str());
  TGraphAsymmErrors *Ele_Rho_Vs_NVtx_Graph = (TGraphAsymmErrors*)file->Get(("Ele_Rho_VS_NVtx"+label+"_Graph").c_str());
   
  Double_t Ele_ChargedIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_ChargedIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_ChargedIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_ChargedIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea = 0;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_NeutralHadronIso007_Vs_NPU_EffArea = 0;
  Double_t Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HoverE_Vs_NPU_EffArea = 0;
  Double_t Ele_HoverE_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HcalDepth1OverEcal_Vs_NPU_EffArea = 0;
  Double_t Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HcalDepth2OverEcal_Vs_NPU_EffArea = 0;
  Double_t Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_TrkIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_TrkIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_EMIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_EMIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HadIso03_Vs_NPU_EffArea = 0;
  Double_t Ele_HadIso03_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_TrkIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_TrkIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_EMIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_EMIso04_Vs_NPU_EffAreaErr = 0;
  Double_t Ele_HadIso04_Vs_NPU_EffArea = 0;
  Double_t Ele_HadIso04_Vs_NPU_EffAreaErr = 0;

  TF1 *f1= 0;

  if (isMC) {
    //*****************************************************************************************************
    //*****************************************************************************************************
    //Vs NPU
    //*****************************************************************************************************
    //*****************************************************************************************************

    //----------------------
    // Rho
    //----------------------
    f1 = new TF1(("Ele_Rho_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_Rho_Vs_NPU_Graph->Fit(("Ele_Rho_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_Rho_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_Rho_Vs_NPU_SlopeErr = f1->GetParError(1);

    //----------------------
    // Energy Observables
    //----------------------
    f1 = new TF1(("Ele_ChargedIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_ChargedIso03_Vs_NPU_Graph->Fit(("Ele_ChargedIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_ChargedIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_ChargedIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_ChargedIso03_Vs_NPU_EffArea = Ele_ChargedIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_ChargedIso03_Vs_NPU_EffAreaErr = Ele_ChargedIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_ChargedIso03_Vs_NPU_SlopeErr/Ele_ChargedIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso03_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso03_Vs_NPU_EffArea = Ele_NeutralHadronIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso03_Vs_NPU_SlopeErr/Ele_NeutralHadronIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIso03_Vs_NPU_Graph->Fit(("Ele_GammaIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIso03_Vs_NPU_EffArea = Ele_GammaIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIso03_Vs_NPU_EffAreaErr = Ele_GammaIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIso03_Vs_NPU_SlopeErr/Ele_GammaIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Graph->Fit(("Ele_GammaIsoVetoEtaStrip03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea = Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr = Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip03_Vs_NPU_SlopeErr/Ele_GammaIsoVetoEtaStrip03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_ChargedIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_ChargedIso04_Vs_NPU_Graph->Fit(("Ele_ChargedIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_ChargedIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_ChargedIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_ChargedIso04_Vs_NPU_EffArea = Ele_ChargedIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_ChargedIso04_Vs_NPU_EffAreaErr = Ele_ChargedIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_ChargedIso04_Vs_NPU_SlopeErr/Ele_ChargedIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso04_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso04_Vs_NPU_EffArea = Ele_NeutralHadronIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso04_Vs_NPU_SlopeErr/Ele_NeutralHadronIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIso04_Vs_NPU_Graph->Fit(("Ele_GammaIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIso04_Vs_NPU_EffArea = Ele_GammaIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIso04_Vs_NPU_EffAreaErr = Ele_GammaIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIso04_Vs_NPU_SlopeErr/Ele_GammaIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Graph->Fit(("Ele_GammaIsoVetoEtaStrip04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea = Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr = Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip04_Vs_NPU_SlopeErr/Ele_GammaIsoVetoEtaStrip04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_NeutralHadronIso007_Vs_NPU_Graph->Fit(("Ele_NeutralHadronIso007_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_NeutralHadronIso007_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_NeutralHadronIso007_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_NeutralHadronIso007_Vs_NPU_EffArea = Ele_NeutralHadronIso007_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr = Ele_NeutralHadronIso007_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_NeutralHadronIso007_Vs_NPU_SlopeErr/Ele_NeutralHadronIso007_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HoverE_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HoverE_Vs_NPU_Graph->Fit(("Ele_HoverE_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HoverE_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HoverE_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HoverE_Vs_NPU_EffArea = Ele_HoverE_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HoverE_Vs_NPU_EffAreaErr = Ele_HoverE_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HoverE_Vs_NPU_SlopeErr/Ele_HoverE_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HcalDepth1OverEcal_Vs_NPU_Graph->Fit(("Ele_HcalDepth1OverEcal_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HcalDepth1OverEcal_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HcalDepth1OverEcal_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HcalDepth1OverEcal_Vs_NPU_EffArea = Ele_HcalDepth1OverEcal_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr = Ele_HcalDepth1OverEcal_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HcalDepth1OverEcal_Vs_NPU_SlopeErr/Ele_HcalDepth1OverEcal_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HcalDepth2OverEcal_Vs_NPU_Graph->Fit(("Ele_HcalDepth2OverEcal_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HcalDepth2OverEcal_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HcalDepth2OverEcal_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HcalDepth2OverEcal_Vs_NPU_EffArea = Ele_HcalDepth2OverEcal_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr = Ele_HcalDepth2OverEcal_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HcalDepth2OverEcal_Vs_NPU_SlopeErr/Ele_HcalDepth2OverEcal_Vs_NPU_Slope,2));


    f1 = new TF1(("Ele_TrkIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_TrkIso03_Vs_NPU_Graph->Fit(("Ele_TrkIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_TrkIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_TrkIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_TrkIso03_Vs_NPU_EffArea = Ele_TrkIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_TrkIso03_Vs_NPU_EffAreaErr = Ele_TrkIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_TrkIso03_Vs_NPU_SlopeErr/Ele_TrkIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_EMIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_EMIso03_Vs_NPU_Graph->Fit(("Ele_EMIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_EMIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_EMIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_EMIso03_Vs_NPU_EffArea = Ele_EMIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_EMIso03_Vs_NPU_EffAreaErr = Ele_EMIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_EMIso03_Vs_NPU_SlopeErr/Ele_EMIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HadIso03_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HadIso03_Vs_NPU_Graph->Fit(("Ele_HadIso03_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HadIso03_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HadIso03_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HadIso03_Vs_NPU_EffArea = Ele_HadIso03_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HadIso03_Vs_NPU_EffAreaErr = Ele_HadIso03_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HadIso03_Vs_NPU_SlopeErr/Ele_HadIso03_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_TrkIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_TrkIso04_Vs_NPU_Graph->Fit(("Ele_TrkIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_TrkIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_TrkIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_TrkIso04_Vs_NPU_EffArea = Ele_TrkIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_TrkIso04_Vs_NPU_EffAreaErr = Ele_TrkIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_TrkIso04_Vs_NPU_SlopeErr/Ele_TrkIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_EMIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_EMIso04_Vs_NPU_Graph->Fit(("Ele_EMIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_EMIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_EMIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_EMIso04_Vs_NPU_EffArea = Ele_EMIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_EMIso04_Vs_NPU_EffAreaErr = Ele_EMIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_EMIso04_Vs_NPU_SlopeErr/Ele_EMIso04_Vs_NPU_Slope,2));

    f1 = new TF1(("Ele_HadIso04_VS_NPU"+label+"_Fit").c_str(), "pol1", 0.5, 30.5);
    Ele_HadIso04_Vs_NPU_Graph->Fit(("Ele_HadIso04_VS_NPU"+label+"_Fit").c_str(),"R");
    Double_t Ele_HadIso04_Vs_NPU_Slope = f1->GetParameter(1);
    Double_t Ele_HadIso04_Vs_NPU_SlopeErr = f1->GetParError(1);
    Ele_HadIso04_Vs_NPU_EffArea = Ele_HadIso04_Vs_NPU_Slope / Ele_Rho_Vs_NPU_Slope;
    Ele_HadIso04_Vs_NPU_EffAreaErr = Ele_HadIso04_Vs_NPU_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NPU_SlopeErr/Ele_Rho_Vs_NPU_Slope,2) + pow(Ele_HadIso04_Vs_NPU_SlopeErr/Ele_HadIso04_Vs_NPU_Slope,2));


  }

  //*****************************************************************************************************
  //*****************************************************************************************************
  //Vs NVtx
  //*****************************************************************************************************
  //*****************************************************************************************************

  //----------------------
  // Rho
  //----------------------
  f1 = new TF1(("Ele_Rho_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_Rho_Vs_NVtx_Graph->Fit(("Ele_Rho_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_Rho_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_Rho_Vs_NVtx_SlopeErr = f1->GetParError(1);

  //----------------------
  // Energy Observables
  //----------------------
  f1 = new TF1(("Ele_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_ChargedIso03_Vs_NVtx_Graph->Fit(("Ele_ChargedIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso03_Vs_NVtx_EffArea = Ele_ChargedIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso03_Vs_NVtx_EffAreaErr = Ele_ChargedIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso03_Vs_NVtx_SlopeErr/Ele_ChargedIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso03_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_EffArea = Ele_NeutralHadronIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso03_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIso03_Vs_NVtx_Graph->Fit(("Ele_GammaIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso03_Vs_NVtx_EffArea = Ele_GammaIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso03_Vs_NVtx_EffAreaErr = Ele_GammaIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso03_Vs_NVtx_SlopeErr/Ele_GammaIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Graph->Fit(("Ele_GammaIsoVetoEtaStrip03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea = Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr = Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_SlopeErr/Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_ChargedIso04_Vs_NVtx_Graph->Fit(("Ele_ChargedIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_ChargedIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_ChargedIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_ChargedIso04_Vs_NVtx_EffArea = Ele_ChargedIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_ChargedIso04_Vs_NVtx_EffAreaErr = Ele_ChargedIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_ChargedIso04_Vs_NVtx_SlopeErr/Ele_ChargedIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso04_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_EffArea = Ele_NeutralHadronIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso04_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIso04_Vs_NVtx_Graph->Fit(("Ele_GammaIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIso04_Vs_NVtx_EffArea = Ele_GammaIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIso04_Vs_NVtx_EffAreaErr = Ele_GammaIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIso04_Vs_NVtx_SlopeErr/Ele_GammaIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Graph->Fit(("Ele_GammaIsoVetoEtaStrip04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea = Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr = Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_SlopeErr/Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_NeutralHadronIso007_Vs_NVtx_Graph->Fit(("Ele_NeutralHadronIso007_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_EffArea = Ele_NeutralHadronIso007_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr = Ele_NeutralHadronIso007_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_NeutralHadronIso007_Vs_NVtx_SlopeErr/Ele_NeutralHadronIso007_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HoverE_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HoverE_Vs_NVtx_Graph->Fit(("Ele_HoverE_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HoverE_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HoverE_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HoverE_Vs_NVtx_EffArea = Ele_HoverE_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HoverE_Vs_NVtx_EffAreaErr = Ele_HoverE_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HoverE_Vs_NVtx_SlopeErr/Ele_HoverE_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HcalDepth1OverEcal_Vs_NVtx_Graph->Fit(("Ele_HcalDepth1OverEcal_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea = Ele_HcalDepth1OverEcal_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr = Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HcalDepth1OverEcal_Vs_NVtx_SlopeErr/Ele_HcalDepth1OverEcal_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HcalDepth2OverEcal_Vs_NVtx_Graph->Fit(("Ele_HcalDepth2OverEcal_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea = Ele_HcalDepth2OverEcal_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr = Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HcalDepth2OverEcal_Vs_NVtx_SlopeErr/Ele_HcalDepth2OverEcal_Vs_NVtx_Slope,2));


  f1 = new TF1(("Ele_TrkIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_TrkIso03_Vs_NVtx_Graph->Fit(("Ele_TrkIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_TrkIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_TrkIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_TrkIso03_Vs_NVtx_EffArea = Ele_TrkIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_TrkIso03_Vs_NVtx_EffAreaErr = Ele_TrkIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_TrkIso03_Vs_NVtx_SlopeErr/Ele_TrkIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_EMIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_EMIso03_Vs_NVtx_Graph->Fit(("Ele_EMIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_EMIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_EMIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_EMIso03_Vs_NVtx_EffArea = Ele_EMIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_EMIso03_Vs_NVtx_EffAreaErr = Ele_EMIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_EMIso03_Vs_NVtx_SlopeErr/Ele_EMIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HadIso03_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HadIso03_Vs_NVtx_Graph->Fit(("Ele_HadIso03_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HadIso03_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HadIso03_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HadIso03_Vs_NVtx_EffArea = Ele_HadIso03_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HadIso03_Vs_NVtx_EffAreaErr = Ele_HadIso03_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HadIso03_Vs_NVtx_SlopeErr/Ele_HadIso03_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_TrkIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_TrkIso04_Vs_NVtx_Graph->Fit(("Ele_TrkIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_TrkIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_TrkIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_TrkIso04_Vs_NVtx_EffArea = Ele_TrkIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_TrkIso04_Vs_NVtx_EffAreaErr = Ele_TrkIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_TrkIso04_Vs_NVtx_SlopeErr/Ele_TrkIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_EMIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_EMIso04_Vs_NVtx_Graph->Fit(("Ele_EMIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_EMIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_EMIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_EMIso04_Vs_NVtx_EffArea = Ele_EMIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_EMIso04_Vs_NVtx_EffAreaErr = Ele_EMIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_EMIso04_Vs_NVtx_SlopeErr/Ele_EMIso04_Vs_NVtx_Slope,2));

  f1 = new TF1(("Ele_HadIso04_VS_NVtx"+label+"_Fit").c_str(), "pol1", 0.5, 20.5);
  Ele_HadIso04_Vs_NVtx_Graph->Fit(("Ele_HadIso04_VS_NVtx"+label+"_Fit").c_str(),"R");
  Double_t Ele_HadIso04_Vs_NVtx_Slope = f1->GetParameter(1);
  Double_t Ele_HadIso04_Vs_NVtx_SlopeErr = f1->GetParError(1);
  Double_t Ele_HadIso04_Vs_NVtx_EffArea = Ele_HadIso04_Vs_NVtx_Slope / Ele_Rho_Vs_NVtx_Slope;
  Double_t Ele_HadIso04_Vs_NVtx_EffAreaErr = Ele_HadIso04_Vs_NVtx_EffArea * TMath::Sqrt(pow( Ele_Rho_Vs_NVtx_SlopeErr/Ele_Rho_Vs_NVtx_Slope,2) + pow(Ele_HadIso04_Vs_NVtx_SlopeErr/Ele_HadIso04_Vs_NVtx_Slope,2));



  cout << endl << endl 
       << "***********************************************************************************"
       << endl << endl ;

  char buffer[200];
  char buffer2[200];

  if (isMC) {
    sprintf(buffer,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NPU_EffArea,Ele_ChargedIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NVtx_EffArea,Ele_ChargedIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NPU_EffArea,Ele_NeutralHadronIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NVtx_EffArea,Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NPU_EffArea,Ele_GammaIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NVtx_EffArea,Ele_GammaIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NPU_EffArea,Ele_ChargedIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NVtx_EffArea,Ele_ChargedIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NPU_EffArea,Ele_NeutralHadronIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NVtx_EffArea,Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NPU_EffArea,Ele_GammaIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NVtx_EffArea,Ele_GammaIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NPU_EffArea,Ele_NeutralHadronIso007_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NVtx_EffArea,Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso007_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HoverE_Vs_NPU_EffArea,Ele_HoverE_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HoverE_Vs_NVtx_EffArea,Ele_HoverE_Vs_NVtx_EffAreaErr);
    cout << "Ele_HoverE_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NPU_EffArea,Ele_HcalDepth1OverEcal_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth1OverEcal_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NPU_EffArea,Ele_HcalDepth2OverEcal_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth2OverEcal_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NPU_EffArea,Ele_TrkIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NVtx_EffArea,Ele_TrkIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_EMIso03_Vs_NPU_EffArea,Ele_EMIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso03_Vs_NVtx_EffArea,Ele_EMIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_HadIso03_Vs_NPU_EffArea,Ele_HadIso03_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso03_Vs_NVtx_EffArea,Ele_HadIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso03_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NPU_EffArea,Ele_TrkIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NVtx_EffArea,Ele_TrkIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_EMIso04_Vs_NPU_EffArea,Ele_EMIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso04_Vs_NVtx_EffArea,Ele_EMIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

    sprintf(buffer,"%.3f +/- %.3f",Ele_HadIso04_Vs_NPU_EffArea,Ele_HadIso04_Vs_NPU_EffAreaErr);
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso04_Vs_NVtx_EffArea,Ele_HadIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso04_Vs_NPU_EffArea : " << buffer << " " << buffer2 << endl;

  } else {

    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso03_Vs_NVtx_EffArea,Ele_ChargedIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso03_Vs_NVtx_EffArea,Ele_NeutralHadronIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso03_Vs_NVtx_EffArea,Ele_GammaIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_ChargedIso04_Vs_NVtx_EffArea,Ele_ChargedIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_ChargedIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso04_Vs_NVtx_EffArea,Ele_NeutralHadronIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIso04_Vs_NVtx_EffArea,Ele_GammaIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea,Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffAreaErr);
    cout << "Ele_GammaIsoVetoEtaStrip04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_NeutralHadronIso007_Vs_NVtx_EffArea,Ele_NeutralHadronIso007_Vs_NVtx_EffAreaErr);
    cout << "Ele_NeutralHadronIso007_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HoverE_Vs_NVtx_EffArea,Ele_HoverE_Vs_NVtx_EffAreaErr);
    cout << "Ele_HoverE_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth1OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth1OverEcal_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.5f +/- %.5f",Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea,Ele_HcalDepth2OverEcal_Vs_NVtx_EffAreaErr);
    cout << "Ele_HcalDepth2OverEcal_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso03_Vs_NVtx_EffArea,Ele_TrkIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso03_Vs_NVtx_EffArea,Ele_EMIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso03_Vs_NVtx_EffArea,Ele_HadIso03_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso03_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_TrkIso04_Vs_NVtx_EffArea,Ele_TrkIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_TrkIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_EMIso04_Vs_NVtx_EffArea,Ele_EMIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_EMIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
     
    sprintf(buffer2,"%.3f +/- %.3f",Ele_HadIso04_Vs_NVtx_EffArea,Ele_HadIso04_Vs_NVtx_EffAreaErr);
    cout << "Ele_HadIso04_Vs_NVtx_EffArea : " << buffer2 << endl;
    

    //*************************************************
    //Make Some Plots
    //*************************************************
    TCanvas *cv = 0;
    TPaveLabel *label = 0;
    
    cv = new TCanvas("cv","cv",800,600);
    Ele_Rho_Vs_NVtx_Graph->GetYaxis()->SetTitle("Energy Density (#rho) [GeV]");
    Ele_Rho_Vs_NVtx_Graph->GetXaxis()->SetTitle("Number of Reconstructed Vertices");
    Ele_Rho_Vs_NVtx_Graph->Draw("AP");
    Ele_Rho_Vs_NVtx_Graph->GetXaxis()->SetRangeUser(0,22);
    label = new TPaveLabel(5.0,20,20.0,22.5, Form("Slope = %.3f +/- %.3f", Ele_Rho_Vs_NVtx_Slope, Ele_Rho_Vs_NVtx_SlopeErr));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(0.6);
    label->Draw();
    cv->SaveAs("Ele_Rho_Vs_NVtx_Graph.eps");    

    cv = new TCanvas("cv","cv",800,600);
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetYaxis()->SetTitleSize(0.046);
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetYaxis()->SetTitleOffset(1.3);
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetYaxis()->SetTitle("Neutral Hadron Isolation (0.3 Cone) [GeV]");
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetXaxis()->SetTitle("Number of Reconstructed Vertices");
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->Draw("AP");
    Ele_NeutralHadronIso03_Vs_NVtx_Graph->GetXaxis()->SetRangeUser(0,22);
    label = new TPaveLabel(5.0,1.6,20.0,1.8, Form("Slope = %.3f +/- %.3f", Ele_NeutralHadronIso03_Vs_NVtx_Slope, Ele_NeutralHadronIso03_Vs_NVtx_SlopeErr));
    label->SetFillStyle(0);
    label->SetBorderSize(0);
    label->SetTextSize(0.6);
    label->Draw();
    cv->SaveAs("Ele_NeutralHadronIso03_Vs_NVtx_Graph.eps");
    
 
  }
   
}

void ComputeElectronIsolationEffectiveArea() {
  
//   ComputeElectronEffectiveArea("Data2011",kFALSE, 0);
   ComputeElectronEffectiveArea("Data2011",kFALSE, 1);
//   ComputeElectronEffectiveArea("Data2011",kFALSE, 2);
//    ComputeElectronEffectiveArea("Data2011",kFALSE, 3);
//   ComputeElectronEffectiveArea("Data2011",kFALSE, 4);

}
