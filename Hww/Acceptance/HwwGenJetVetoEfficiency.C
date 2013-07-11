//root -l EWKAna/Hww/Acceptance/HwwGenJetVetoEfficiency.C+\(\"HwwNtuple_f10-h160ww2l-gf-z2-v12_noskim_0000.root\",\"HWW\"\)
//root -l EWKAna/Hww/Acceptance/HwwGenJetVetoEfficiency.C+\(\"HwwNtuple_f10-ww2l-z2-v12_noskim_0000.root\",\"WW\"\)
//root -l EWKAna/Hww/Acceptance/HwwGenJetVetoEfficiency.C+\(\"/home/sixie/hist/HwwAcceptance/ggHWW160GenOnly_default.root\",\"ggHww160_default\"\)
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

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "TLegend.h"

#endif
 

//=== FUNCTION DECLARATIONS ======================================================================================
//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HwwNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}




//=== MAIN MACRO =================================================================================================

void HwwGenJetVetoEfficiency(const string inputFilename,  const string Label) 
{  
  gBenchmark->Start("WWTemplate");

  cout << "here1\n";
  string label = "";
  cout << "here11\n";
  if (Label != "") label = "_" + Label;
  cout << "here12\n";

  TFile *file = new TFile("JetVetoEfficiencySystematics.root", "UPDATE");
  cout << "here13\n";
  assert(file);
//   TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPtFineBinning_ggHww160_PowhegToNNLL");
  TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww160_PowhegToNNLL");
  cout << "here14\n";
//    TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww160_MCAtNLOToNNLL");
//  TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww180_PowhegToNNLL");//  
//  TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww200_PowhegToNNLL");
//   TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww220_PowhegToNNLL");
//  TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww250_PowhegToNNLL");
//  TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww300_PowhegToNNLL");
//   TH1D *kFactorHiggsPt = (TH1D*)file->Get("kFactorHiggsPt_ggHww400_PowhegToNNLL");
  assert(kFactorHiggsPt);
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    
  cout << "here2\n";

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *genMet = new TH1D("genMet", ";  Met [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonPtMin = new TH1D("leptonPtMin", ";  Lepton p_{T} [GeV/c]; Number of Events", 50, 0, 200);
  TH1D *leptonEta = new TH1D("leptonEta", "; Lepton #eta; Number of Events", 50, -5, 5);
  TH1D *leadingGenJetPt = new TH1D((string("leadingGenJetPt")+ label).c_str() , "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *leadingGenJetPt_reweighted = new TH1D((string("leadingGenJetPt_reweighted")+ label).c_str() , "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *secondGenJetPt_reweighted = new TH1D((string("leadingGenJetPt_reweighted")+ label).c_str() , "; Leading GenJet p_{T} [GeV/c]; Number of Events", 200, 0, 200);
  TH1D *nGenJets = new TH1D("nGenJets", "; Mass [GeV/c^{2}]; Number of Events", 10, -0.5, 9.5);
  TH1D *leptonDeltaPhi = new TH1D("leptonDeltaPhi", "; #Delta #phi (l,l); Number of Events", 50, 0, 180);


  vector<Double_t> NEvents_0Jets;
  vector<Double_t> NEvents_1Jets;
  vector<Double_t> NEvents_2Jets;
  for(UInt_t k=0 ; k < 200; ++k) {
    NEvents_0Jets.push_back(0); 
    NEvents_1Jets.push_back(0); 
    NEvents_2Jets.push_back(0);
  }


  Double_t NEventsGenerated = 0;
  Double_t NEventsAccepted = 0;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  mithep::TGenInfo *genInfo    = new mithep::TGenInfo();

  cout << "here3\n";

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Gen", &genInfo); TBranch *genInfoBr = eventTree->GetBranch("Gen");
   
  cout << "here4\n";

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    genInfoBr->GetEntry(ientry);
		
    NEventsGenerated++;

   mithep::FourVectorM lepton1;
    mithep::FourVectorM lepton2;
    lepton1.SetCoordinates(genInfo->pt_1, genInfo->eta_1, genInfo->phi_1, 0.51099892e-3 );
    lepton2.SetCoordinates(genInfo->pt_2, genInfo->eta_2, genInfo->phi_2, 0.51099892e-3 );
    mithep::FourVectorM dilepton = lepton1+lepton2;

    if (genInfo->pt_1 > 20.0 && genInfo->pt_2 > 20.0
        && fabs(genInfo->eta_1) < 2.5 && fabs(genInfo->eta_2) < 2.5
        && genInfo->met > 30.0
        && fabs(dilepton.M() - 91.2) > 15.0
      ) {
      Double_t scale = kFactorHiggsPt->GetBinContent(kFactorHiggsPt->GetXaxis()->FindFixBin(genInfo->ptBosonSystem));
      scale = 1;

      if (info->eventweight >= 0) {
        
        //don't use a scale factor for Higgs pt < 1.0, because HQT does negative.
        if (genInfo->ptBosonSystem < 1.25 ) scale = 1.0;

        //powheg
        if (genInfo->ptBosonSystem > 200) scale = 0.6;
        //MC@NLO
//          if (genInfo->ptBosonSystem > 200) scale = 1.5;        

         leadingGenJetPt->Fill(genInfo->jetpt_1,1.0);
         leadingGenJetPt_reweighted->Fill(genInfo->jetpt_1,scale);

         for(UInt_t k=0 ; k < 200; ++k) {
           if (genInfo->jetpt_1 > k && genInfo->jetpt_2 > k) {
             NEvents_2Jets[k] += scale;
           } else if (genInfo->jetpt_1 > k) {
             NEvents_1Jets[k] += scale; 
           } else {
             NEvents_0Jets[k] += scale;
           }                                
         }

      } else {       
        leadingGenJetPt_reweighted->Fill(genInfo->jetpt_1,-1*scale);
        leadingGenJetPt->Fill(genInfo->jetpt_1,-1*scale);
        
        for(UInt_t k=0 ; k < 200; ++k) {
          if (genInfo->jetpt_1 > k && genInfo->jetpt_2 > k) {
            NEvents_2Jets[k] += -1*scale;
          } else if (genInfo->jetpt_1 > k) {
            NEvents_1Jets[k] += -1*scale; 
          } else {
             NEvents_0Jets[k] += -1*scale;
          }                                
        }

      }
    }
  } //end loop over data     

  delete info;
  delete genInfo;

  //--------------------------------------------------------------------------------------------------------------
  // Compute Jet Veto Efficiency
  //==============================================================================================================
  TH1D *jetVetoEfficiency = (TH1D*)leadingGenJetPt->Clone((string("jetVetoEfficiency")+label).c_str()); 
  jetVetoEfficiency->GetYaxis()->SetTitle("Jet Veto Efficiency");
  jetVetoEfficiency->GetYaxis()->SetTitleOffset(1.5);
  for (int i=1; i<jetVetoEfficiency->GetXaxis()->GetNbins()+1; ++i) {
    Int_t n = leadingGenJetPt->Integral(1,i);
    Int_t d = leadingGenJetPt->Integral(1,leadingGenJetPt->GetXaxis()->GetNbins()+1);
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Int_t errorType = 2;
    mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);

    jetVetoEfficiency->SetBinContent(i,ratio);
    jetVetoEfficiency->SetBinError(i,(errLow+errHigh)/2); 
    cout << jetVetoEfficiency->GetXaxis()->GetBinLowEdge(i) << " : " << n << " / " << d << " = " 
         << ratio << " + " << errHigh << " - " << errLow << endl;
  }

  TH1D *jetVetoEfficiency_reweighted = (TH1D*)leadingGenJetPt_reweighted->Clone((string("jetVetoEfficiency_reweighted")+label).c_str()); 
  jetVetoEfficiency_reweighted->GetYaxis()->SetTitle("Jet Veto Efficiency");
  jetVetoEfficiency_reweighted->GetYaxis()->SetTitleOffset(1.5);
  for (int i=1; i<jetVetoEfficiency_reweighted->GetXaxis()->GetNbins()+1; ++i) {
    Int_t n = leadingGenJetPt_reweighted->Integral(1,i);
    Int_t d = leadingGenJetPt_reweighted->Integral(1,leadingGenJetPt_reweighted->GetXaxis()->GetNbins()+1);
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Int_t errorType = 2;
    mithep::MathUtils::CalcRatio(n , d, ratio, errLow, errHigh, errorType);

    jetVetoEfficiency_reweighted->SetBinContent(i,ratio);
    jetVetoEfficiency_reweighted->SetBinError(i,(errLow+errHigh)/2); 
    cout << jetVetoEfficiency_reweighted->GetXaxis()->GetBinLowEdge(i) << " : " << n << " / " << d << " = " 
         << ratio << " + " << errHigh << " - " << errLow << endl;
  }





  //*****************************************************************************************************
  //0,1,2 - Jet Fractions
  //*****************************************************************************************************
  const int nPoints = 100;
  double jetPtCut[nPoints];
  double jetPtCutError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    jetPtCut[i] = i;
    jetPtCutError[i] = 0;     
  }

  double ZeroJetFraction[nPoints];
  double ZeroJetFractionErrLow[nPoints];
  double ZeroJetFractionErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(NEvents_0Jets[i]);
    Double_t n2 = TMath::Nint(NEvents_0Jets[i]+NEvents_1Jets[i]+NEvents_2Jets[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ZeroJetFraction[i] = ratio;
    ZeroJetFractionErrLow[i] = errLow;
    ZeroJetFractionErrHigh[i] = errHigh;
  }
  TGraphAsymmErrors *ZeroJetFractionVsJetPtCut = new TGraphAsymmErrors (nPoints,  jetPtCut, ZeroJetFraction, jetPtCutError, jetPtCutError, ZeroJetFractionErrLow, ZeroJetFractionErrHigh);
  ZeroJetFractionVsJetPtCut->SetName(("ZeroJetFractionVsJetPtCut"+label).c_str());
  ZeroJetFractionVsJetPtCut->SetMarkerColor(kRed);
  ZeroJetFractionVsJetPtCut->GetXaxis()->SetTitleOffset(1.02);
  ZeroJetFractionVsJetPtCut->GetYaxis()->SetTitleOffset(1.05);



  double OneJetFraction[nPoints];
  double OneJetFractionErrLow[nPoints];
  double OneJetFractionErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(NEvents_1Jets[i]);
    Double_t n2 = TMath::Nint(NEvents_0Jets[i]+NEvents_1Jets[i]+NEvents_2Jets[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    OneJetFraction[i] = ratio;
    OneJetFractionErrLow[i] = errLow;
    OneJetFractionErrHigh[i] = errHigh;
  }
  TGraphAsymmErrors *OneJetFractionVsJetPtCut = new TGraphAsymmErrors (nPoints,  jetPtCut, OneJetFraction, jetPtCutError, jetPtCutError, OneJetFractionErrLow, OneJetFractionErrHigh);
  OneJetFractionVsJetPtCut->SetName(("OneJetFractionVsJetPtCut"+label).c_str());
  OneJetFractionVsJetPtCut->SetMarkerColor(kBlue);
  OneJetFractionVsJetPtCut->GetXaxis()->SetTitleOffset(1.02);
  OneJetFractionVsJetPtCut->GetYaxis()->SetTitleOffset(1.05);



  double TwoJetFraction[nPoints];
  double TwoJetFractionErrLow[nPoints];
  double TwoJetFractionErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(NEvents_2Jets[i]);
    Double_t n2 = TMath::Nint(NEvents_0Jets[i]+NEvents_1Jets[i]+NEvents_2Jets[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    TwoJetFraction[i] = ratio;
    TwoJetFractionErrLow[i] = errLow;
    TwoJetFractionErrHigh[i] = errHigh;
  }
  TGraphAsymmErrors *TwoJetFractionVsJetPtCut = new TGraphAsymmErrors (nPoints,  jetPtCut, TwoJetFraction, jetPtCutError, jetPtCutError, TwoJetFractionErrLow, TwoJetFractionErrHigh);
  TwoJetFractionVsJetPtCut->SetName(("TwoJetFractionVsJetPtCut"+label).c_str());
  TwoJetFractionVsJetPtCut->SetMarkerColor(kMagenta);
  TwoJetFractionVsJetPtCut->GetXaxis()->SetTitleOffset(1.02);
  TwoJetFractionVsJetPtCut->GetYaxis()->SetTitleOffset(1.05);




  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *tmpLegend = new TLegend(0.73,0.55,0.93,0.70); 
  

  tmpLegend->SetTextSize(0.03);
  tmpLegend->SetBorderSize(1);
  tmpLegend->AddEntry(leadingGenJetPt, "Powheg", "LP");   
  tmpLegend->AddEntry(leadingGenJetPt_reweighted, "MC@NLO(reweighted)", "LP");   
  cv->SetLogy();
  leadingGenJetPt->Draw();
  leadingGenJetPt->SetLineColor(kRed);
  leadingGenJetPt_reweighted->Draw("same");
  tmpLegend->Draw();
  cv->SaveAs("leadingGenJetPt.gif");

  jetVetoEfficiency->Draw();
  cv->SaveAs("jetVetoEfficiency.gif");
 
  cv->SetLogy(0);
  tmpLegend->Clear();
  tmpLegend->AddEntry(ZeroJetFractionVsJetPtCut, "0Jet Fraction", "LP");   
  tmpLegend->AddEntry(OneJetFractionVsJetPtCut, "1Jet Fraction", "LP");   
  tmpLegend->AddEntry(TwoJetFractionVsJetPtCut, ">2 Jet Fraction", "LP");   
  ZeroJetFractionVsJetPtCut->GetYaxis()->SetRangeUser(0.0,1.0);
  ZeroJetFractionVsJetPtCut->Draw("AP");
  OneJetFractionVsJetPtCut->Draw("Psame");
  TwoJetFractionVsJetPtCut->Draw("Psame");
  cv->SaveAs("NJetFractions.gif");


  file->WriteTObject(leadingGenJetPt, leadingGenJetPt->GetName(), "WriteDelete");
  file->WriteTObject(leadingGenJetPt_reweighted, leadingGenJetPt_reweighted->GetName(), "WriteDelete");
  file->WriteTObject(jetVetoEfficiency, jetVetoEfficiency->GetName(), "WriteDelete");
  file->WriteTObject(jetVetoEfficiency_reweighted, jetVetoEfficiency_reweighted->GetName(), "WriteDelete");
  file->WriteTObject(ZeroJetFractionVsJetPtCut, ZeroJetFractionVsJetPtCut->GetName(), "WriteDelete");
  file->WriteTObject(OneJetFractionVsJetPtCut, OneJetFractionVsJetPtCut->GetName(), "WriteDelete");
  file->WriteTObject(TwoJetFractionVsJetPtCut, TwoJetFractionVsJetPtCut->GetName(), "WriteDelete");
//   file->Close();

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
 


        
  gBenchmark->Show("plotZ");       
} 

