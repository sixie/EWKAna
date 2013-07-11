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
void DoFakeMuonPUReweighting()
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


  TH1F *histRhoRealMuonTarget = new TH1F("RhoRealMuonTarget", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);
  TH1F *histRhoFakeMuonSource = new TH1F("RhoFakeMuonSource", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);




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
  //FakeMuTree Output
  //*****************************************************************************************
  TFile *FakeMuOutputFile = new TFile("MuonSelectionTraining.Fake.PtAndPUWeighted.root", "RECREATE");
  TTree *FakeMuOutputTree = new TTree("Muons","Muons");
  FakeMuOutputTree->SetAutoFlush(0);

  FakeMuOutputTree->Branch("weight",&fWeight,"weight/F");
  FakeMuOutputTree->Branch("run",&fRunNumber,"run/i");
  FakeMuOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  FakeMuOutputTree->Branch("event",&fEventNumber,"event/i");
  FakeMuOutputTree->Branch("pt",&fMuPt,"pt/F"); 
  FakeMuOutputTree->Branch("eta",&fMuEta,"eta/F"); 
  FakeMuOutputTree->Branch("phi",&fMuPhi,"phi/F"); 
  FakeMuOutputTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  FakeMuOutputTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  FakeMuOutputTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  FakeMuOutputTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  FakeMuOutputTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  FakeMuOutputTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  FakeMuOutputTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  FakeMuOutputTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  FakeMuOutputTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  FakeMuOutputTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  FakeMuOutputTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  FakeMuOutputTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  FakeMuOutputTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  FakeMuOutputTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  FakeMuOutputTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  FakeMuOutputTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  FakeMuOutputTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  FakeMuOutputTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  FakeMuOutputTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  FakeMuOutputTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  FakeMuOutputTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  FakeMuOutputTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  FakeMuOutputTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  FakeMuOutputTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  FakeMuOutputTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  FakeMuOutputTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  FakeMuOutputTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  FakeMuOutputTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  FakeMuOutputTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  FakeMuOutputTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  FakeMuOutputTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  FakeMuOutputTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  FakeMuOutputTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  FakeMuOutputTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
  FakeMuOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  FakeMuOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 




  //*****************************************************************************************
  //Generate Target PU distribution from Real Muon Sample
  //*****************************************************************************************
  TFile *RealMuTargetFile = new TFile("MuonSelectionTraining.Real.weighted.root", "READ");
  TTree *RealMuTargetTree = (TTree*)RealMuTargetFile->Get("Muons");
  RealMuTargetTree->SetBranchAddress( "weight", &fWeight);
  RealMuTargetTree->SetBranchAddress( "Rho", &fRho); 
  RealMuTargetTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealMuTargetTree->GetEntries(); ientry++) {       	
    RealMuTargetTree->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Real Muon Sample " << ientry << endl;
    histRhoRealMuonTarget->Fill(fRho, fWeight);
  } 


  //*****************************************************************************************
  //Generate Source PU distribution from Fake Muon Sample
  //*****************************************************************************************
  TFile *FakeMuSourceFile = new TFile("MuonSelectionTraining.Fake.weighted.root", "READ");
  TTree *FakeMuSourceTree = (TTree*)FakeMuSourceFile->Get("Muons");
  FakeMuSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeMuSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeMuSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeMuSourceTree->GetEntries(); ientry++) {       	
    FakeMuSourceTree->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Fake Muon Sample " << ientry << endl;
    histRhoFakeMuonSource->Fill(fRho, fWeight);
  } 

  //*****************************************************************************************
  //Create PU reweight histogram
  //*****************************************************************************************
  NormalizeHist(histRhoRealMuonTarget);
  NormalizeHist(histRhoFakeMuonSource);

  TH1F *FakeMuonPUReweightFactor = (TH1F*)histRhoFakeMuonSource->Clone("FakeMuonPUReweightFactor");
  FakeMuonPUReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < FakeMuonPUReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histRhoFakeMuonSource->GetBinContent(a)>0) {
      FakeMuonPUReweightFactor->SetBinContent(a,histRhoRealMuonTarget->GetBinContent(a) / histRhoFakeMuonSource->GetBinContent(a));
    } else {
      FakeMuonPUReweightFactor->SetBinContent(a,1.0);
    }
  }

  
  //*****************************************************************************************
  //FakeMuSourceTree
  //*****************************************************************************************
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
    fWeight = fWeight * FakeMuonPUReweightFactor->GetBinContent(FakeMuonPUReweightFactor->GetXaxis()->FindFixBin(fRho));
    FakeMuOutputTree->Fill();
  } 
  FakeMuOutputFile->Write();
  FakeMuOutputFile->Close();
  delete FakeMuSourceFile;


  gBenchmark->Show("WWTemplate");       
} 



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void ReweightMuonPU() {

   DoFakeMuonPUReweighting();

}
