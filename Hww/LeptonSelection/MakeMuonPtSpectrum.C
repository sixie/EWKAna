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

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"

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
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void FillLeptonPtSpectrum()
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histRealMuPtSource = new TH1F("RealMuPtSource", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);
  TH1F *histRealMuPtTarget = new TH1F("RealMuPtTarget", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);
  TH1F *histFakeMuPtSource = new TH1F("FakeMuPtSource", "; p_{T} [GeV/c] ; Number of Events ",  50, 0 , 100);


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 


  
  //*****************************************************************************************
  //RealMuSourceTree
  //*****************************************************************************************
  TFile *RealMuSourceFile = new TFile("MuonSelectionTraining.Real.root", "READ");
  TTree *RealMuSourceTree = (TTree*)RealMuSourceFile->Get("Muons");
  RealMuSourceTree->SetBranchAddress( "weight", &fWeight);
  RealMuSourceTree->SetBranchAddress( "run", &fRunNumber);
  RealMuSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealMuSourceTree->SetBranchAddress( "event", &fEventNumber);
  RealMuSourceTree->SetBranchAddress( "pt", &fMuPt); 
  RealMuSourceTree->SetBranchAddress( "eta", &fMuEta); 
  RealMuSourceTree->SetBranchAddress( "phi", &fMuPhi); 
  RealMuSourceTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  RealMuSourceTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  RealMuSourceTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  RealMuSourceTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  RealMuSourceTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  RealMuSourceTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  RealMuSourceTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  RealMuSourceTree->SetBranchAddress( "D0", &fMuD0); 
  RealMuSourceTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  RealMuSourceTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  RealMuSourceTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  RealMuSourceTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  RealMuSourceTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  RealMuSourceTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  RealMuSourceTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  RealMuSourceTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  RealMuSourceTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  RealMuSourceTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  RealMuSourceTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  RealMuSourceTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  RealMuSourceTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  RealMuSourceTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  RealMuSourceTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  RealMuSourceTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  RealMuSourceTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  RealMuSourceTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  RealMuSourceTree->SetBranchAddress( "Rho", &fRho); 
  RealMuSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealMuSourceTree->GetEntries(); ientry++) {       	
    RealMuSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histRealMuPtSource->Fill(fMuPt);
  } 
  


  //*****************************************************************************************
  //RealEleTargetTree
  //*****************************************************************************************
  TFile *RealMuTargetFile = new TFile("MuonSelectionTraining.HWW115.root", "READ");
  TTree *RealMuTargetTree = (TTree*)RealMuTargetFile->Get("Muons");
  RealMuTargetTree->SetBranchAddress( "weight", &fWeight);
  RealMuTargetTree->SetBranchAddress( "run", &fRunNumber);
  RealMuTargetTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealMuTargetTree->SetBranchAddress( "event", &fEventNumber);
  RealMuTargetTree->SetBranchAddress( "pt", &fMuPt); 
  RealMuTargetTree->SetBranchAddress( "eta", &fMuEta); 
  RealMuTargetTree->SetBranchAddress( "phi", &fMuPhi); 
  RealMuTargetTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  RealMuTargetTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  RealMuTargetTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  RealMuTargetTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  RealMuTargetTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  RealMuTargetTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  RealMuTargetTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  RealMuTargetTree->SetBranchAddress( "D0", &fMuD0); 
  RealMuTargetTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  RealMuTargetTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  RealMuTargetTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  RealMuTargetTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  RealMuTargetTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  RealMuTargetTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  RealMuTargetTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  RealMuTargetTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  RealMuTargetTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  RealMuTargetTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  RealMuTargetTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  RealMuTargetTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  RealMuTargetTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  RealMuTargetTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  RealMuTargetTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  RealMuTargetTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  RealMuTargetTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  RealMuTargetTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  RealMuTargetTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  RealMuTargetTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  RealMuTargetTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  RealMuTargetTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  RealMuTargetTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  RealMuTargetTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  RealMuTargetTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  RealMuTargetTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  RealMuTargetTree->SetBranchAddress( "Rho", &fRho); 
  RealMuTargetTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealMuTargetTree->GetEntries(); ientry++) {       	
    RealMuTargetTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histRealMuPtTarget->Fill(fMuPt);
  } 


  //*****************************************************************************************
  //FakeMuSourceTree
  //*****************************************************************************************
  TFile *FakeMuSourceFile = new TFile("MuonSelectionTraining.Fake.root", "READ");
  TTree *FakeMuSourceTree = (TTree*)FakeMuSourceFile->Get("Muons");
  FakeMuSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeMuSourceTree->SetBranchAddress( "run", &fRunNumber);
  FakeMuSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeMuSourceTree->SetBranchAddress( "event", &fEventNumber);
  FakeMuSourceTree->SetBranchAddress( "pt", &fMuPt); 
  FakeMuSourceTree->SetBranchAddress( "eta", &fMuEta); 
  FakeMuSourceTree->SetBranchAddress( "phi", &fMuPhi); 
  FakeMuSourceTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  FakeMuSourceTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  FakeMuSourceTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  FakeMuSourceTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  FakeMuSourceTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  FakeMuSourceTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  FakeMuSourceTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  FakeMuSourceTree->SetBranchAddress( "D0", &fMuD0); 
  FakeMuSourceTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  FakeMuSourceTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  FakeMuSourceTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  FakeMuSourceTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  FakeMuSourceTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  FakeMuSourceTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  FakeMuSourceTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  FakeMuSourceTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  FakeMuSourceTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  FakeMuSourceTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  FakeMuSourceTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  FakeMuSourceTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  FakeMuSourceTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  FakeMuSourceTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  FakeMuSourceTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  FakeMuSourceTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  FakeMuSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeMuSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeMuSourceTree->GetEntries(); ientry++) {       	
    FakeMuSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    histFakeMuPtSource->Fill(fMuPt);
  } 



  TFile *file = new TFile("MuonPtSpectrum.root", "UPDATE");
  file->cd();
  file->WriteTObject(histRealMuPtSource, histRealMuPtSource->GetName(), "WriteDelete");  
  file->WriteTObject(histRealMuPtTarget, histRealMuPtTarget->GetName(), "WriteDelete");  
  file->WriteTObject(histFakeMuPtSource, histFakeMuPtSource->GetName(), "WriteDelete");  
  file->Close();
  delete file;



  gBenchmark->Show("WWTemplate");       
} 

//*************************************************************************************************
//Compute Reweight Factors
//*************************************************************************************************
void  MakeReweightFactors() {
  
  TFile* file = new TFile("MuonPtSpectrum.root","READ");
  TH1F *histRealMuPtSource = (TH1F*)file->Get("RealMuPtSource");
  TH1F *histRealMuPtTarget = (TH1F*)file->Get("RealMuPtTarget");
  TH1F *histFakeMuPtSource = (TH1F*)file->Get("FakeMuPtSource");
  TH1F *histFakeMuPtTarget = (TH1F*)file->Get("FakeMuPtTarget");

  NormalizeHist(histRealMuPtSource);
  NormalizeHist(histRealMuPtTarget);
  NormalizeHist(histFakeMuPtSource);
  NormalizeHist(histFakeMuPtTarget);


  TH1F *RealMuonPtReweightFactor = (TH1F*)histRealMuPtSource->Clone("RealMuonPtReweightFactor");
  RealMuonPtReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < RealMuonPtReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histRealMuPtSource->GetBinContent(a)>0) {
      RealMuonPtReweightFactor->SetBinContent(a,histRealMuPtTarget->GetBinContent(a) / histRealMuPtSource->GetBinContent(a));
      cout << a << " " << histRealMuPtTarget->GetBinContent(a) << " / " << histRealMuPtSource->GetBinContent(a) << " = " << histRealMuPtTarget->GetBinContent(a) / histRealMuPtSource->GetBinContent(a) << endl;
    } else {
      RealMuonPtReweightFactor->SetBinContent(a,1.0);
    }
  }

  TH1F *FakeMuonPtReweightFactor = (TH1F*)histFakeMuPtSource->Clone("FakeMuonPtReweightFactor");
  FakeMuonPtReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < FakeMuonPtReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histFakeMuPtSource->GetBinContent(a)>0) {
      FakeMuonPtReweightFactor->SetBinContent(a,histFakeMuPtTarget->GetBinContent(a) / histFakeMuPtSource->GetBinContent(a));
    } else {
      FakeMuonPtReweightFactor->SetBinContent(a,1.0);
    }
  }

  TFile *outfile = new TFile("MuonPtSpectrum.root", "UPDATE");
  outfile->cd();
  outfile->WriteTObject(RealMuonPtReweightFactor, RealMuonPtReweightFactor->GetName(), "WriteDelete");  
  outfile->WriteTObject(FakeMuonPtReweightFactor, FakeMuonPtReweightFactor->GetName(), "WriteDelete");  
  outfile->Close();
  delete outfile;

  TCanvas *cv = new TCanvas("cv","cv", 800, 600);
  histRealMuPtSource->GetXaxis()->SetTitleOffset(1.05);
  histRealMuPtSource->GetYaxis()->SetTitleOffset(1.4);
  histRealMuPtSource->Draw("hist");
  cv->SaveAs("SignalPt_Source.png");

  histRealMuPtTarget->GetXaxis()->SetTitleOffset(1.05);
  histRealMuPtTarget->GetYaxis()->SetTitleOffset(1.4);
  histRealMuPtTarget->Draw("hist");
  cv->SaveAs("SignalPt_Target.png");

  histFakeMuPtSource->GetXaxis()->SetTitleOffset(1.05);
  histFakeMuPtSource->GetYaxis()->SetTitleOffset(1.4);
  histFakeMuPtSource->Draw("hist");
  cv->SaveAs("BkgPt_Source.png");

  histFakeMuPtTarget->GetXaxis()->SetTitleOffset(1.05);
  histFakeMuPtTarget->GetYaxis()->SetTitleOffset(1.4);
  histFakeMuPtTarget->Draw("hist");
  cv->SaveAs("BkgPt_Target.png");


}


//*************************************************************************************************
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void UpdateNtupleWeights()
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TFile* inputfile = new TFile("MuonPtSpectrum.root","READ");
  TH1F *RealMuonPtReweightFactor = (TH1F*)inputfile->Get("RealMuonPtReweightFactor");
  RealMuonPtReweightFactor->SetDirectory(0);
  TH1F *FakeMuonPtReweightFactor = (TH1F*)inputfile->Get("FakeMuonPtReweightFactor");
  FakeMuonPtReweightFactor->SetDirectory(0);
  inputfile->Close();
  delete inputfile;


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  //*****************************************************************************************
  //RealMuSourceTree Output
  //*****************************************************************************************
  TFile *RealMuSourceOutputFile = new TFile("MuonSelectionTraining.Real.weighted.root", "RECREATE");
  TTree *RealMuSourceOutputTree = new TTree("Muons","Muons");
  RealMuSourceOutputTree->SetAutoFlush(0);

  RealMuSourceOutputTree->Branch("weight",&fWeight,"weight/F");
  RealMuSourceOutputTree->Branch("run",&fRunNumber,"run/i");
  RealMuSourceOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  RealMuSourceOutputTree->Branch("event",&fEventNumber,"event/i");
  RealMuSourceOutputTree->Branch("pt",&fMuPt,"pt/F"); 
  RealMuSourceOutputTree->Branch("eta",&fMuEta,"eta/F"); 
  RealMuSourceOutputTree->Branch("phi",&fMuPhi,"phi/F"); 
  RealMuSourceOutputTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  RealMuSourceOutputTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  RealMuSourceOutputTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  RealMuSourceOutputTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  RealMuSourceOutputTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  RealMuSourceOutputTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  RealMuSourceOutputTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  RealMuSourceOutputTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  RealMuSourceOutputTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  RealMuSourceOutputTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  RealMuSourceOutputTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  RealMuSourceOutputTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  RealMuSourceOutputTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  RealMuSourceOutputTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  RealMuSourceOutputTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  RealMuSourceOutputTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  RealMuSourceOutputTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  RealMuSourceOutputTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  RealMuSourceOutputTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  RealMuSourceOutputTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  RealMuSourceOutputTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  RealMuSourceOutputTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  RealMuSourceOutputTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  RealMuSourceOutputTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  RealMuSourceOutputTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  RealMuSourceOutputTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  RealMuSourceOutputTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  RealMuSourceOutputTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  RealMuSourceOutputTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  RealMuSourceOutputTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  RealMuSourceOutputTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  RealMuSourceOutputTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  RealMuSourceOutputTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  RealMuSourceOutputTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
  RealMuSourceOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  RealMuSourceOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 


  
  //*****************************************************************************************
  //RealMuSourceTree
  //*****************************************************************************************
  TFile *RealMuSourceFile = new TFile("MuonSelectionTraining.Real.root", "READ");
  TTree *RealMuSourceTree = (TTree*)RealMuSourceFile->Get("Muons");
  RealMuSourceTree->SetBranchAddress( "weight", &fWeight);
  RealMuSourceTree->SetBranchAddress( "run", &fRunNumber);
  RealMuSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealMuSourceTree->SetBranchAddress( "event", &fEventNumber);
  RealMuSourceTree->SetBranchAddress( "pt", &fMuPt); 
  RealMuSourceTree->SetBranchAddress( "eta", &fMuEta); 
  RealMuSourceTree->SetBranchAddress( "phi", &fMuPhi); 
  RealMuSourceTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  RealMuSourceTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  RealMuSourceTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  RealMuSourceTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  RealMuSourceTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  RealMuSourceTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  RealMuSourceTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  RealMuSourceTree->SetBranchAddress( "D0", &fMuD0); 
  RealMuSourceTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  RealMuSourceTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  RealMuSourceTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  RealMuSourceTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  RealMuSourceTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  RealMuSourceTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  RealMuSourceTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  RealMuSourceTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  RealMuSourceTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  RealMuSourceTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  RealMuSourceTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  RealMuSourceTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  RealMuSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  RealMuSourceTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  RealMuSourceTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  RealMuSourceTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  RealMuSourceTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  RealMuSourceTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  RealMuSourceTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  RealMuSourceTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  RealMuSourceTree->SetBranchAddress( "Rho", &fRho); 
  RealMuSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealMuSourceTree->GetEntries(); ientry++) {       	
    RealMuSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Muon " << ientry << endl;
    fWeight = RealMuonPtReweightFactor->GetBinContent(RealMuonPtReweightFactor->GetXaxis()->FindFixBin(fMuPt));
    RealMuSourceOutputTree->Fill();
  } 
  RealMuSourceOutputFile->Write();
  RealMuSourceOutputFile->Close();
  delete RealMuSourceFile;









  //*****************************************************************************************
  //FakeMuSourceTree Output
  //*****************************************************************************************
  TFile *FakeMuSourceOutputFile = new TFile("MuonSelectionTraining.Fake.weighted.root", "RECREATE");
  TTree *FakeMuSourceOutputTree = new TTree("Muons","Muons");
  FakeMuSourceOutputTree->SetAutoFlush(0);

  FakeMuSourceOutputTree->Branch("weight",&fWeight,"weight/F");
  FakeMuSourceOutputTree->Branch("run",&fRunNumber,"run/i");
  FakeMuSourceOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  FakeMuSourceOutputTree->Branch("event",&fEventNumber,"event/i");
  FakeMuSourceOutputTree->Branch("pt",&fMuPt,"pt/F"); 
  FakeMuSourceOutputTree->Branch("eta",&fMuEta,"eta/F"); 
  FakeMuSourceOutputTree->Branch("phi",&fMuPhi,"phi/F"); 
  FakeMuSourceOutputTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  FakeMuSourceOutputTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  FakeMuSourceOutputTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  FakeMuSourceOutputTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  FakeMuSourceOutputTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  FakeMuSourceOutputTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  FakeMuSourceOutputTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  FakeMuSourceOutputTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  FakeMuSourceOutputTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  FakeMuSourceOutputTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  FakeMuSourceOutputTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  FakeMuSourceOutputTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  FakeMuSourceOutputTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  FakeMuSourceOutputTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  FakeMuSourceOutputTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  FakeMuSourceOutputTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  FakeMuSourceOutputTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  FakeMuSourceOutputTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  FakeMuSourceOutputTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  FakeMuSourceOutputTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  FakeMuSourceOutputTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  FakeMuSourceOutputTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  FakeMuSourceOutputTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  FakeMuSourceOutputTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  FakeMuSourceOutputTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  FakeMuSourceOutputTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  FakeMuSourceOutputTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  FakeMuSourceOutputTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  FakeMuSourceOutputTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  FakeMuSourceOutputTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  FakeMuSourceOutputTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  FakeMuSourceOutputTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  FakeMuSourceOutputTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  FakeMuSourceOutputTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
  FakeMuSourceOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  FakeMuSourceOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 


  
  //*****************************************************************************************
  //FakeMuSourceTree
  //*****************************************************************************************
  TFile *FakeMuSourceFile = new TFile("MuonSelectionTraining.Fake.root", "READ");
  TTree *FakeMuSourceTree = (TTree*)FakeMuSourceFile->Get("Muons");
  FakeMuSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeMuSourceTree->SetBranchAddress( "run", &fRunNumber);
  FakeMuSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeMuSourceTree->SetBranchAddress( "event", &fEventNumber);
  FakeMuSourceTree->SetBranchAddress( "pt", &fMuPt); 
  FakeMuSourceTree->SetBranchAddress( "eta", &fMuEta); 
  FakeMuSourceTree->SetBranchAddress( "phi", &fMuPhi); 
  FakeMuSourceTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  FakeMuSourceTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  FakeMuSourceTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  FakeMuSourceTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  FakeMuSourceTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  FakeMuSourceTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  FakeMuSourceTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  FakeMuSourceTree->SetBranchAddress( "D0", &fMuD0); 
  FakeMuSourceTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  FakeMuSourceTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  FakeMuSourceTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  FakeMuSourceTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  FakeMuSourceTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  FakeMuSourceTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  FakeMuSourceTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  FakeMuSourceTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  FakeMuSourceTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  FakeMuSourceTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  FakeMuSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  FakeMuSourceTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  FakeMuSourceTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  FakeMuSourceTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  FakeMuSourceTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  FakeMuSourceTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  FakeMuSourceTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  FakeMuSourceTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  FakeMuSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeMuSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeMuSourceTree->GetEntries(); ientry++) {       	
    FakeMuSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Muon " << ientry << endl;
    fWeight = FakeMuonPtReweightFactor->GetBinContent(FakeMuonPtReweightFactor->GetXaxis()->FindFixBin(fMuPt));
    FakeMuSourceOutputTree->Fill();
  } 
  FakeMuSourceOutputFile->Write();
  FakeMuSourceOutputFile->Close();
  delete FakeMuSourceFile;



  gBenchmark->Show("WWTemplate");       
} 



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void MakeMuonPtSpectrum() {
  FillLeptonPtSpectrum();
  MakeReweightFactors();
  UpdateNtupleWeights();
  
}
