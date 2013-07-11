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
void DoFakeElectronPUReweighting()
{  


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TFile* inputfile = new TFile("ElectronPtSpectrum.root","READ");
  TH1F *RealElectronPtReweightFactor = (TH1F*)inputfile->Get("RealElectronPtReweightFactor");
  RealElectronPtReweightFactor->SetDirectory(0);
  TH1F *FakeElectronPtReweightFactor = (TH1F*)inputfile->Get("FakeElectronPtReweightFactor");
  FakeElectronPtReweightFactor->SetDirectory(0);
  inputfile->Close();
  delete inputfile;


  TH1F *histRhoRealElectronTarget = new TH1F("RhoRealElectronTarget", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);
  TH1F *histRhoFakeElectronSource = new TH1F("RhoFakeElectronSource", "; #rho [GeV] ; Number of Events ",  50, 0 , 50);




  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fElePFIso; 
  
  //CutBased Variables
  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  //Additional Vars used in Likelihood
  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleCovIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;

  //Isolation Variables
  Float_t                 fEleStandardLikelihood; 
  Float_t                 fElePFMVA; 
  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 



  //*****************************************************************************************
  //FakeEleTree Output
  //*****************************************************************************************
  TFile *FakeEleOutputFile = new TFile("ElectronSelectionTraining.Fake.PtAndPUWeighted.root", "RECREATE");
  TTree *FakeEleOutputTree = new TTree("Electrons","Electrons");
  FakeEleOutputTree->SetAutoFlush(0);

  FakeEleOutputTree->Branch("weight",&fWeight,"weight/F");
  FakeEleOutputTree->Branch("run",&fRunNumber,"run/i");
  FakeEleOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  FakeEleOutputTree->Branch("event",&fEventNumber,"event/i");
  FakeEleOutputTree->Branch("pt",&fElePt,"pt/F"); 
  FakeEleOutputTree->Branch("eta",&fEleEta,"eta/F"); 
  FakeEleOutputTree->Branch("phi",&fElePhi,"phi/F"); 
  FakeEleOutputTree->Branch("scet",&fEleSCEt,"scet/F"); 
  FakeEleOutputTree->Branch("sceta",&fEleSCEta,"sceta/F"); 
  FakeEleOutputTree->Branch("scphi",&fEleSCPhi,"scphi/F"); 
  FakeEleOutputTree->Branch("pfiso",&fElePFIso,"pfiso/F"); 
  
  //CutBased Variables
  FakeEleOutputTree->Branch("SigmaIEtaIEta",&fEleSigmaIEtaIEta,"SigmaIEtaIEta/F"); 
  FakeEleOutputTree->Branch("DEtaIn",&fEleDEtaIn,"DEtaIn/F"); 
  FakeEleOutputTree->Branch("DPhiIn",&fEleDPhiIn,"DPhiIn/F"); 
  FakeEleOutputTree->Branch("HoverE",&fEleHoverE,"HoverE/F"); 
  FakeEleOutputTree->Branch("D0",&fEleD0,"D0/F");
  FakeEleOutputTree->Branch("DZ",&fEleDZ,"DZ/F"); 
  FakeEleOutputTree->Branch("FBrem",&fEleFBrem,"FBrem/F"); 
  FakeEleOutputTree->Branch("EOverP",&fEleEOverP,"EOverP/F"); 

  //Additional Vars used in Likelihood
  FakeEleOutputTree->Branch("ESeedClusterOverPout",&fEleESeedClusterOverPout,"ESeedClusterOverPout/F"); 
  FakeEleOutputTree->Branch("SigmaIPhiIPhi",&fEleSigmaIPhiIPhi,"SigmaIPhiIPhi/F"); 
  FakeEleOutputTree->Branch("NBrem",&fEleNBrem,"NBrem/F"); 
  FakeEleOutputTree->Branch("OneOverEMinusOneOverP",&fEleOneOverEMinusOneOverP,"OneOverEMinusOneOverP/F"); 
  FakeEleOutputTree->Branch("ESeedClusterOverPIn",&fEleESeedClusterOverPIn,"ESeedClusterOverPIn/F"); 
  FakeEleOutputTree->Branch("IP3d",&fEleIP3d,"IP3d/F"); 
  FakeEleOutputTree->Branch("IP3dSig",&fEleIP3dSig,"IP3dSig/F"); 

  FakeEleOutputTree->Branch("HcalDepth1OverEcal",&fEleHcalDepth1OverEcal,"HcalDepth1OverEcal/F"); 
  FakeEleOutputTree->Branch("HcalDepth2OverEcal",&fEleHcalDepth2OverEcal,"HcalDepth2OverEcal/F"); 
  FakeEleOutputTree->Branch("dEtaCalo",&fEledEtaCalo,"dEtaCalo/F"); 
  FakeEleOutputTree->Branch("dPhiCalo",&fEledPhiCalo,"dPhiCalo/F"); 
  FakeEleOutputTree->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
  FakeEleOutputTree->Branch("CovIEtaIPhi",&fEleCovIEtaIPhi,"CovIEtaIPhi/F"); 
  FakeEleOutputTree->Branch("SCEtaWidth",&fEleSCEtaWidth,"SCEtaWidth/F"); 
  FakeEleOutputTree->Branch("SCPhiWidth",&fEleSCPhiWidth,"SCPhiWidth/F"); 
  FakeEleOutputTree->Branch("GsfTrackChi2OverNdof",&fEleGsfTrackChi2OverNdof,"GsfTrackChi2OverNdof/F"); 
  FakeEleOutputTree->Branch("R9",&fEleR9,"R9/F"); 

  FakeEleOutputTree->Branch("SeedEMaxOverE",&fEleSeedEMaxOverE,"SeedEMaxOverE/F"); 
  FakeEleOutputTree->Branch("SeedETopOverE",&fEleSeedETopOverE,"SeedETopOverE/F"); 
  FakeEleOutputTree->Branch("SeedEBottomOverE",&fEleSeedEBottomOverE,"SeedEBottomOverE/F"); 
  FakeEleOutputTree->Branch("SeedELeftOverE",&fEleSeedELeftOverE,"SeedELeftOverE/F"); 
  FakeEleOutputTree->Branch("SeedERightOverE",&fEleSeedERightOverE,"SeedERightOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2ndOverE",&fEleSeedE2ndOverE,"SeedE2ndOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE,"SeedE2x5RightOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE,"SeedE2x5LeftOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x5BottomOverE",&fEleSeedE2x5BottomOverE,"SeedE2x5BottomOverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  FakeEleOutputTree->Branch("SeedE1x3OverE",&fEleSeedE1x3OverE,"SeedE1x3OverE/F"); 
  FakeEleOutputTree->Branch("SeedE3x1OverE",&fEleSeedE3x1OverE,"SeedE3x1OverE/F"); 
  FakeEleOutputTree->Branch("SeedE1x5OverE",&fEleSeedE1x5OverE,"SeedE1x5OverE/F"); 
  FakeEleOutputTree->Branch("SeedE2x2OverE",&fEleSeedE2x2OverE,"SeedE2x2OverE/F"); 
  FakeEleOutputTree->Branch("SeedE3x2OverE",&fEleSeedE3x2OverE,"SeedE3x2OverE/F"); 
  FakeEleOutputTree->Branch("SeedE3x3OverE",&fEleSeedE3x3OverE,"SeedE3x3OverE/F"); 
  FakeEleOutputTree->Branch("SeedE4x4OverE",&fEleSeedE4x4OverE,"SeedE4x4OverE/F"); 
  FakeEleOutputTree->Branch("SeedE5x5OverE",&fEleSeedE5x5OverE,"SeedE5x5OverE/F"); 

  //Isolation Variables
  FakeEleOutputTree->Branch("StandardLikelihood",&fEleStandardLikelihood,"StandardLikelihood/F"); 
  FakeEleOutputTree->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
  FakeEleOutputTree->Branch("ChargedIso03",&fEleChargedIso03,"ChargedIso03/F"); 
  FakeEleOutputTree->Branch("NeutralHadronIso03",&fEleNeutralHadronIso03,"NeutralHadronIso03/F"); 
  FakeEleOutputTree->Branch("GammaIso03",&fEleGammaIso03,"GammaIso03/F"); 
  FakeEleOutputTree->Branch("ChargedIso04",&fEleChargedIso04,"ChargedIso04/F"); 
  FakeEleOutputTree->Branch("NeutralHadronIso04",&fEleNeutralHadronIso04,"NeutralHadronIso04/F"); 
  FakeEleOutputTree->Branch("GammaIso04",&fEleGammaIso04,"GammaIso04/F"); 
  FakeEleOutputTree->Branch("ChargedIso04FromOtherVertices",&fEleChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  FakeEleOutputTree->Branch("NeutralHadronIso04_10Threshold",&fEleNeutralHadronIso04_10Threshold,"NeutralHadronIso04_10Threshold/F"); 
  FakeEleOutputTree->Branch("GammaIso04_10Threshold",&fEleGammaIso04_10Threshold,"GammaIso04_10Threshold/F"); 
  FakeEleOutputTree->Branch("TrkIso03",&fEleTrkIso03,"TrkIso03/F"); 
  FakeEleOutputTree->Branch("EMIso03",&fEleEMIso03,"EMIso03/F"); 
  FakeEleOutputTree->Branch("HadIso03",&fEleHadIso03,"HadIso03/F"); 
  FakeEleOutputTree->Branch("TrkIso04",&fEleTrkIso04,"TrkIso04/F"); 
  FakeEleOutputTree->Branch("EMIso04",&fEleEMIso04,"EMIso04/F"); 
  FakeEleOutputTree->Branch("HadIso04",&fEleHadIso04,"HadIso04/F"); 
  FakeEleOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  FakeEleOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 



  //*****************************************************************************************
  //Generate Target PU distribution from Real Electron Sample
  //*****************************************************************************************
  TFile *RealEleTargetFile = new TFile("ElectronSelectionTraining.Real.weighted.root", "READ");
  TTree *RealEleTargetTree = (TTree*)RealEleTargetFile->Get("Electrons");
  RealEleTargetTree->SetBranchAddress( "weight", &fWeight);
  RealEleTargetTree->SetBranchAddress( "Rho", &fRho); 
  RealEleTargetTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealEleTargetTree->GetEntries(); ientry++) {       	
    RealEleTargetTree->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Real Electron Sample " << ientry << endl;
    histRhoRealElectronTarget->Fill(fRho, fWeight);
  } 


  //*****************************************************************************************
  //Generate Source PU distribution from Fake Electron Sample
  //*****************************************************************************************
  TFile *FakeEleSourceFile = new TFile("ElectronSelectionTraining.Fake.weighted.root", "READ");
  TTree *FakeEleSourceTree = (TTree*)FakeEleSourceFile->Get("Electrons");
  FakeEleSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeEleSourceTree->GetEntries(); ientry++) {       	
    FakeEleSourceTree->GetEntry(ientry);    
    if (ientry % 100000 == 0) cout << "Fake Electron Sample " << ientry << endl;
    histRhoFakeElectronSource->Fill(fRho, fWeight);
  } 

  //*****************************************************************************************
  //Create PU reweight histogram
  //*****************************************************************************************
  NormalizeHist(histRhoRealElectronTarget);
  NormalizeHist(histRhoFakeElectronSource);

  TH1F *FakeElectronPUReweightFactor = (TH1F*)histRhoFakeElectronSource->Clone("FakeElectronPUReweightFactor");
  FakeElectronPUReweightFactor->SetBinContent(0,1.0);
  for(UInt_t a=1; a < FakeElectronPUReweightFactor->GetXaxis()->GetNbins()+2; ++a) {
    if (histRhoFakeElectronSource->GetBinContent(a)>0) {
      FakeElectronPUReweightFactor->SetBinContent(a,histRhoRealElectronTarget->GetBinContent(a) / histRhoFakeElectronSource->GetBinContent(a));
    } else {
      FakeElectronPUReweightFactor->SetBinContent(a,1.0);
    }
  }

  //*****************************************************************************************
  //FakeEleSourceTree
  //*****************************************************************************************
  FakeEleSourceTree->SetBranchAddress( "weight", &fWeight);
  FakeEleSourceTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleSourceTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleSourceTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleSourceTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleSourceTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleSourceTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleSourceTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleSourceTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleSourceTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleSourceTree->SetBranchAddress( "pfiso", &fElePFIso); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  FakeEleSourceTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  FakeEleSourceTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  FakeEleSourceTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  FakeEleSourceTree->SetBranchAddress( "D0", &fEleD0); 
  FakeEleSourceTree->SetBranchAddress( "DZ", &fEleDZ); 
  FakeEleSourceTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  FakeEleSourceTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  FakeEleSourceTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  FakeEleSourceTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  FakeEleSourceTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  FakeEleSourceTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleSourceTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  FakeEleSourceTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig);
  FakeEleSourceTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  FakeEleSourceTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  FakeEleSourceTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  FakeEleSourceTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  FakeEleSourceTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  FakeEleSourceTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  FakeEleSourceTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  FakeEleSourceTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  FakeEleSourceTree->SetBranchAddress( "R9", &fEleR9); 
  FakeEleSourceTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  FakeEleSourceTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  FakeEleSourceTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  FakeEleSourceTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleSourceTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleSourceTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  FakeEleSourceTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  FakeEleSourceTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  FakeEleSourceTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  FakeEleSourceTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  FakeEleSourceTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  FakeEleSourceTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleSourceTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeEleSourceTree->GetEntries(); ientry++) {       	
    FakeEleSourceTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Electron " << ientry << endl;
    fWeight = fWeight * FakeElectronPUReweightFactor->GetBinContent(FakeElectronPUReweightFactor->GetXaxis()->FindFixBin(fRho));
    FakeEleOutputTree->Fill();
  } 
  FakeEleOutputFile->Write();
  FakeEleOutputFile->Close();
  delete FakeEleSourceFile;


  gBenchmark->Show("WWTemplate");       
}



//*************************************************************************************************
//Main Function
//*************************************************************************************************
void ReweightElectronPU() {

   DoFakeElectronPUReweighting();

}
