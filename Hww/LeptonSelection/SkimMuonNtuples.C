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
//Fill Lepton Pt Spectrum
//*************************************************************************************************
void SkimMuonNtuples(string label, Int_t option)
{  

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
  //RealMuTree Output
  //*****************************************************************************************
  TFile *RealMuOutputFile = new TFile(("MuonSelectionTraining.Real.weighted." + label + ".root").c_str(), "RECREATE");
  TTree *RealMuOutputTree = new TTree("Muons","Muons");
  RealMuOutputTree->SetAutoFlush(0);

  RealMuOutputTree->Branch("weight",&fWeight,"weight/F");
  RealMuOutputTree->Branch("run",&fRunNumber,"run/i");
  RealMuOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  RealMuOutputTree->Branch("event",&fEventNumber,"event/i");
  RealMuOutputTree->Branch("pt",&fMuPt,"pt/F"); 
  RealMuOutputTree->Branch("eta",&fMuEta,"eta/F"); 
  RealMuOutputTree->Branch("phi",&fMuPhi,"phi/F"); 
  RealMuOutputTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  RealMuOutputTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  RealMuOutputTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  RealMuOutputTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  RealMuOutputTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  RealMuOutputTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  RealMuOutputTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  RealMuOutputTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  RealMuOutputTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  RealMuOutputTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  RealMuOutputTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  RealMuOutputTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  RealMuOutputTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  RealMuOutputTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  RealMuOutputTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  RealMuOutputTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  RealMuOutputTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  RealMuOutputTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  RealMuOutputTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  RealMuOutputTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  RealMuOutputTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  RealMuOutputTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  RealMuOutputTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  RealMuOutputTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  RealMuOutputTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  RealMuOutputTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  RealMuOutputTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  RealMuOutputTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  RealMuOutputTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  RealMuOutputTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  RealMuOutputTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  RealMuOutputTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  RealMuOutputTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  RealMuOutputTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
  RealMuOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  RealMuOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 


  
  //*****************************************************************************************
  //RealMuTree
  //*****************************************************************************************
  TFile *RealMuFile = new TFile("MuonSelectionTraining.Real.weighted.root", "READ");
  TTree *RealMuTree = (TTree*)RealMuFile->Get("Muons");
  RealMuTree->SetBranchAddress( "weight", &fWeight);
  RealMuTree->SetBranchAddress( "run", &fRunNumber);
  RealMuTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealMuTree->SetBranchAddress( "event", &fEventNumber);
  RealMuTree->SetBranchAddress( "pt", &fMuPt); 
  RealMuTree->SetBranchAddress( "eta", &fMuEta); 
  RealMuTree->SetBranchAddress( "phi", &fMuPhi); 
  RealMuTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  RealMuTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  RealMuTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  RealMuTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  RealMuTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  RealMuTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  RealMuTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  RealMuTree->SetBranchAddress( "D0", &fMuD0); 
  RealMuTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  RealMuTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  RealMuTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  RealMuTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  RealMuTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  RealMuTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  RealMuTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  RealMuTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  RealMuTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  RealMuTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  RealMuTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  RealMuTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  RealMuTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  RealMuTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  RealMuTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  RealMuTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  RealMuTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  RealMuTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  RealMuTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  RealMuTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  RealMuTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  RealMuTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  RealMuTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  RealMuTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  RealMuTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  RealMuTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  RealMuTree->SetBranchAddress( "Rho", &fRho); 
  RealMuTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealMuTree->GetEntries(); ientry++) {       	
    RealMuTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Muon " << ientry << endl;

    Bool_t passCuts = kFALSE;
    
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fMuEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fMuPt > 14.5) ptBin = 1;
    if (fMuPt > 20.0) ptBin = 2;

    if (option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (option == 5) passCuts = (subdet == 1 && ptBin == 2);

    if (passCuts) {
      RealMuOutputTree->Fill();
    }
  } 
  RealMuOutputFile->Write();
  RealMuOutputFile->Close();
  delete RealMuFile;

  //*****************************************************************************************
  //FakeEleTree Output
  //*****************************************************************************************
  TFile *FakeMuOutputFile = new TFile(("MuonSelectionTraining.Fake.PtAndPUWeighted." + label + ".root").c_str(), "RECREATE");
  TTree *FakeMuOutputTree = new TTree("Muons","Muons");
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
  //FakeMuTargetTree
  //*****************************************************************************************
  TFile *FakeMuFile = new TFile("MuonSelectionTraining.Fake.PtAndPUWeighted.root", "READ");
  TTree *FakeMuTree = (TTree*)FakeMuFile->Get("Muons");
  FakeMuTree->SetBranchAddress( "weight", &fWeight);
  FakeMuTree->SetBranchAddress( "run", &fRunNumber);
  FakeMuTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeMuTree->SetBranchAddress( "event", &fEventNumber);
  FakeMuTree->SetBranchAddress( "pt", &fMuPt); 
  FakeMuTree->SetBranchAddress( "eta", &fMuEta); 
  FakeMuTree->SetBranchAddress( "phi", &fMuPhi); 
  FakeMuTree->SetBranchAddress( "pfiso", &fMuPFIso); 
  FakeMuTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  FakeMuTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  FakeMuTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  FakeMuTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  FakeMuTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  FakeMuTree->SetBranchAddress( "NMatches", &fMuNMatches); 
  FakeMuTree->SetBranchAddress( "D0", &fMuD0); 
  FakeMuTree->SetBranchAddress( "IP3d", &fMuIP3d); 
  FakeMuTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  FakeMuTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  FakeMuTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  FakeMuTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  FakeMuTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  FakeMuTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  FakeMuTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  FakeMuTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  FakeMuTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  FakeMuTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  FakeMuTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  FakeMuTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  FakeMuTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  FakeMuTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  FakeMuTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  FakeMuTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  FakeMuTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  FakeMuTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  FakeMuTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  FakeMuTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  FakeMuTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  FakeMuTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  FakeMuTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  FakeMuTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  FakeMuTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  FakeMuTree->SetBranchAddress( "Rho", &fRho); 
  FakeMuTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeMuTree->GetEntries(); ientry++) {       	
    FakeMuTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Muon " << ientry << endl;

    Bool_t passCuts = kFALSE;
    
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fMuEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fMuPt > 15.0) ptBin = 1;
    if (fMuPt > 20.0) ptBin = 2;

    if (option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (option == 2) passCuts = (subdet == 0 && ptBin == 1);
    if (option == 3) passCuts = (subdet == 1 && ptBin == 1);
    if (option == 4) passCuts = (subdet == 0 && ptBin == 2);
    if (option == 5) passCuts = (subdet == 1 && ptBin == 2);

    if (passCuts) {
      FakeMuOutputTree->Fill();
    }
  } 
  FakeMuOutputFile->Write();
  FakeMuOutputFile->Close();
  delete FakeMuFile;

  gBenchmark->Show("WWTemplate");       
} 



