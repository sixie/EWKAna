//root -l EWKAna/Hww/Acceptance/PlotJetVetoEfficiencySystematics.C+\(\"JetVetoEfficiencySystematics.root\"\)
//================================================================================================
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TStyle.h>                 // class to handle ROOT plotting style
#include <TCanvas.h>                // class for drawing
#include <TBenchmark.h>             // class to track macro running statistics
#include <iostream>                 // standard I/O
#include <iomanip>
#include <fstream>


// RooFit headers

#include "TFile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TMath.h"

#endif

#define LUMINOSITY 2.88 //(in pb^-1)
#define NBINSPASS 60
#define NBINSFAIL 24



void PrintJetBinFractions() {
  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  vector<TGraphAsymmErrors*> HwwZeroJetFraction;
  vector<TGraphAsymmErrors*> HwwOneJetFraction;
  vector<TGraphAsymmErrors*> HwwTwoJetFraction;

  HwwZeroJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("ZeroJetFractionVsJetPtCut_ggHww160_default"));
  HwwZeroJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("ZeroJetFractionVsJetPtCut_ggHww160_fscaleUp"));
  HwwZeroJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("ZeroJetFractionVsJetPtCut_ggHww160_fscaleDown"));
  HwwZeroJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("ZeroJetFractionVsJetPtCut_ggHww160_rscaleUp"));
  HwwZeroJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("ZeroJetFractionVsJetPtCut_ggHww160_rscaleDown"));
  HwwOneJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("OneJetFractionVsJetPtCut_ggHww160_default"));
  HwwOneJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("OneJetFractionVsJetPtCut_ggHww160_fscaleUp"));
  HwwOneJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("OneJetFractionVsJetPtCut_ggHww160_fscaleDown"));
  HwwOneJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("OneJetFractionVsJetPtCut_ggHww160_rscaleUp"));
  HwwOneJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("OneJetFractionVsJetPtCut_ggHww160_rscaleDown"));
  HwwTwoJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("TwoJetFractionVsJetPtCut_ggHww160_default"));
  HwwTwoJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("TwoJetFractionVsJetPtCut_ggHww160_fscaleUp"));
  HwwTwoJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("TwoJetFractionVsJetPtCut_ggHww160_fscaleDown"));
  HwwTwoJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("TwoJetFractionVsJetPtCut_ggHww160_rscaleUp"));
  HwwTwoJetFraction.push_back((TGraphAsymmErrors*)EffFile->Get("TwoJetFractionVsJetPtCut_ggHww160_rscaleDown"));

  for (int i=0; i<HwwZeroJetFraction.size(); ++i) {
    cout << i << endl;
    assert(HwwZeroJetFraction[i]);
    double x;
    double ZeroJetFraction;
    double OneJetFraction;
    double TwoJetFraction;
    HwwZeroJetFraction[i]->GetPoint(30, x,ZeroJetFraction);
    HwwOneJetFraction[i]->GetPoint(30, x,OneJetFraction);
    HwwTwoJetFraction[i]->GetPoint(30, x,TwoJetFraction);
    cout << HwwZeroJetFraction[i]->GetName() << " : " << ZeroJetFraction << " " << OneJetFraction << " " << TwoJetFraction << endl;
  }
  

}




void PlotJetVetoEfficiencySystematics(string inputFile = "" ) {
  PrintJetBinFractions();
  return;


  TFile *EffFile = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");

  vector<TH1D*> HwwJetVetoEfficiencyHists;
  vector<TH1D*> ZeeJetVetoEfficiencyHists;


//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_default_PowhegToMCAtNLOReweighted"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_default_PowhegToMCAtNLOReweighted"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_default_PowhegToMCAtNLOReweighted"));


//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_MCAtNLO_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_MCAtNLO_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_MCAtNLO_default"));


//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_default"));

//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleDown_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_reweighted_ggHww160_fscaleUp_rscaleUp"));

//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_MCAtNLO_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_MCAtNLO_fscaleUp_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_MCAtNLO_fscaleDown_rscaleDown"));


  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_default"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleUp"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_rscaleDown"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp_rscaleUp"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleUp_rscaleDown"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown_rscaleUp"));
  HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww160_fscaleDown_rscaleDown"));

//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleUp_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww200_fscaleDown_rscaleDown"));

//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleUp_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww250_fscaleDown_rscaleDown"));

//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_default"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleUp_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleUp_rscaleDown"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleDown_rscaleUp"));
//   HwwJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_ggHww400_fscaleDown_rscaleDown"));
   
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_default"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleUp"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleDown"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_rscaleUp"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_rscaleDown"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleUp_rscaleUp"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleUp_rscaleDown"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleDown_rscaleUp"));
  ZeeJetVetoEfficiencyHists.push_back((TH1D*)EffFile->Get("jetVetoEfficiency_Zee_fscaleDown_rscaleDown"));
  
  assert(HwwJetVetoEfficiencyHists.size() == ZeeJetVetoEfficiencyHists.size());
  for (int i=0; i<HwwJetVetoEfficiencyHists.size(); ++i) {
    cout << "Hww: " << i << endl;
    assert(HwwJetVetoEfficiencyHists[i]);
    cout << "zee: " << i << endl;
    assert(ZeeJetVetoEfficiencyHists[i]);
  }
  
  Double_t HwwXSDefault = 49.739376991;
  Double_t HwwXSScaleUp = 54.553365096;
  Double_t HwwXSScaleDown = 45.6830791089;


  //Construct Hww-Zee Jet Veto Efficiency Scale Factor
  TH1D* JetVetoEfficiencyScaleFactor = (TH1D*)HwwJetVetoEfficiencyHists[0]->Clone("HwwZllJetVetoEfficiencyScaleFactor");

  for (int i=1; i<JetVetoEfficiencyScaleFactor->GetXaxis()->GetNbins()+1; ++i) {

    Double_t scaleFactor = HwwJetVetoEfficiencyHists[0]->GetBinContent(i)/ZeeJetVetoEfficiencyHists[0]->GetBinContent(i);
    Double_t statError = scaleFactor * TMath::Sqrt(pow(HwwJetVetoEfficiencyHists[0]->GetBinError(i)/HwwJetVetoEfficiencyHists[0]->GetBinContent(i),2) + 
                                                   pow(ZeeJetVetoEfficiencyHists[0]->GetBinError(i)/ZeeJetVetoEfficiencyHists[0]->GetBinContent(i),2));
    Double_t sysErrorBand = 0;
    cout << "bin " << i << " " << scaleFactor << " : " ;
    for (int j=0; j < HwwJetVetoEfficiencyHists.size(); ++j) {
      for (int k=0; k < ZeeJetVetoEfficiencyHists.size(); ++k) {      
        Double_t tmp = HwwJetVetoEfficiencyHists[j]->GetBinContent(i)/ZeeJetVetoEfficiencyHists[k]->GetBinContent(i);
        if (sysErrorBand < fabs(tmp - scaleFactor)) {
          sysErrorBand = fabs(tmp - scaleFactor);
        }
        cout << j << " " << k << " " << HwwJetVetoEfficiencyHists[j]->GetBinContent(i) << "/" << ZeeJetVetoEfficiencyHists[k]->GetBinContent(i) << "=" << tmp << " \n";
      }
    }
    cout << endl;

    Double_t totalError = TMath::Sqrt(pow(statError,2) + pow(sysErrorBand,2));
    JetVetoEfficiencyScaleFactor->SetBinContent(i,scaleFactor);
    JetVetoEfficiencyScaleFactor->SetBinError(i,totalError);


    cout << "bin : " << i << " : " << scaleFactor << " " << statError << " " << sysErrorBand << " " << totalError << endl;

  }


  TH1D* HwwJetVetoCrossSection = (TH1D*)HwwJetVetoEfficiencyHists[0]->Clone("HwwJetVetoCrossSection");
  HwwJetVetoCrossSection->GetYaxis()->SetTitle("Cross Section (Jet Veto) [pb]");
  for (int i=1; i<JetVetoEfficiencyScaleFactor->GetXaxis()->GetNbins()+1; ++i) {

    Double_t XSDefault = HwwJetVetoEfficiencyHists[0]->GetBinContent(i) * HwwXSDefault;
    Double_t XSScaleUp = HwwJetVetoEfficiencyHists[5]->GetBinContent(i) * HwwXSScaleUp;
    Double_t XSScaleDown = HwwJetVetoEfficiencyHists[8]->GetBinContent(i) * HwwXSScaleDown;

    Double_t error = fabs(XSDefault - XSScaleDown);
    if (fabs(XSDefault - XSScaleUp) > error) error = fabs(XSDefault - XSScaleUp);

    HwwJetVetoCrossSection->SetBinContent(i,XSDefault);
    HwwJetVetoCrossSection->SetBinError(i,error);

    cout << i << " " << HwwJetVetoEfficiencyHists[0]->GetXaxis()->GetBinLowEdge(i) << " : " << HwwJetVetoEfficiencyHists[0]->GetBinContent(i) << " " << HwwJetVetoEfficiencyHists[5]->GetBinContent(i) << " " << HwwJetVetoEfficiencyHists[8]->GetBinContent(i) << endl;

  }





  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.73,0.25,0.93,0.40);   
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(ZeeJetVetoEfficiencyHists[0], "Default", "LP");   
  tmpLegend->AddEntry(ZeeJetVetoEfficiencyHists[5], "ScaleUp", "LP");   
  tmpLegend->AddEntry(ZeeJetVetoEfficiencyHists[7], "ScaleDown", "LP");   
  tmpLegend->Draw();

  ZeeJetVetoEfficiencyHists[0]->Draw("E");
  ZeeJetVetoEfficiencyHists[5]->SetLineColor(kRed);
  ZeeJetVetoEfficiencyHists[5]->Draw("Esame");
  ZeeJetVetoEfficiencyHists[7]->SetLineColor(kBlue);
  ZeeJetVetoEfficiencyHists[7]->Draw("Esame");

  cv->SaveAs("ZeeJetVetoEfficiencyScaleFactor.gif");


  tmpLegend->Clear();
  tmpLegend->AddEntry(HwwJetVetoEfficiencyHists[0], "Default", "LP");   
  tmpLegend->AddEntry(HwwJetVetoEfficiencyHists[5], "ScaleUp", "LP");   
  tmpLegend->AddEntry(HwwJetVetoEfficiencyHists[7], "ScaleDown", "LP");   
  tmpLegend->Draw();

  HwwJetVetoEfficiencyHists[0]->Draw("E");
  HwwJetVetoEfficiencyHists[5]->SetLineColor(kRed);
  HwwJetVetoEfficiencyHists[5]->Draw("Esame");
  HwwJetVetoEfficiencyHists[7]->SetLineColor(kBlue);
  HwwJetVetoEfficiencyHists[7]->Draw("Esame");

  cv->SaveAs("HwwJetVetoEfficiencyScaleFactor.gif");


  JetVetoEfficiencyScaleFactor->Draw("E3");
  JetVetoEfficiencyScaleFactor->SetMarkerSize(0);
  JetVetoEfficiencyScaleFactor->SetFillColor(kBlue);
  JetVetoEfficiencyScaleFactor->SetFillStyle(3001);
  JetVetoEfficiencyScaleFactor->GetYaxis()->SetTitleOffset(1.3);
  JetVetoEfficiencyScaleFactor->GetYaxis()->SetTitle("Jet Veto Efficiency Ratio");
  JetVetoEfficiencyScaleFactor->GetYaxis()->SetRangeUser(0.5, 0.8);
  JetVetoEfficiencyScaleFactor->GetXaxis()->SetRangeUser(20,50);  
  cv->SaveAs("jetVetoEfficiencySF_ScaleVariation_HWW160.gif");
  cv->SaveAs("jetVetoEfficiencySF_ScaleVariation_HWW160.pdf");


  HwwJetVetoCrossSection->Draw("E3");
  HwwJetVetoCrossSection->SetMaximum(60);
  HwwJetVetoCrossSection->SetMinimum(10);
  HwwJetVetoCrossSection->SetMarkerSize(0.0);
  HwwJetVetoCrossSection->SetFillColor(kBlue);
  HwwJetVetoCrossSection->SetFillStyle(3001);
  HwwJetVetoCrossSection->GetYaxis()->SetTitleOffset(1.3);
  HwwJetVetoCrossSection->GetYaxis()->SetTitle("Cross Section (after Jet Veto) [pb]");
  HwwJetVetoCrossSection->GetXaxis()->SetTitleOffset(1.05);
  HwwJetVetoCrossSection->GetXaxis()->SetTitle("GenJet Veto Threshold [GeV/c]");
  HwwJetVetoCrossSection->GetXaxis()->SetRangeUser(10,90);  
  TGraph *HwwJetVetoCrossSectionNNLO =  (TGraph*)EffFile->Get("HWWJetVetoCrossSection");
  HwwJetVetoCrossSectionNNLO->SetLineColor(kRed);
  HwwJetVetoCrossSectionNNLO->SetLineWidth(2);
  HwwJetVetoCrossSectionNNLO->Draw("L,same");

  tmpLegend = new TLegend(0.53,0.25,0.93,0.40);   
  tmpLegend->Clear();
  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(0);
  tmpLegend->AddEntry(HwwJetVetoCrossSection, "Powheg Band", "F");   
  tmpLegend->AddEntry(HwwJetVetoCrossSectionNNLO, "NNLO Band", "L");   
  tmpLegend->Draw();


  cv->SaveAs("HwwJetVetoCrossSection.gif");




//   EffFile->Close();

}


