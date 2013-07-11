//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0LowPt",0)'
//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1LowPt",1)'
//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2LowPt",2)'
//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet0HighPt",3)'
//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet1HighPt",4)'
//root -l -b -q EWKAna/Hww/LeptonSelection/SkimElectronNtuples.C+'("Subdet2HighPt",5)'

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
void SkimElectronNtuples(string label, Int_t option)
{  

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
  //RealEleTree Output
  //*****************************************************************************************
  TFile *RealEleOutputFile = new TFile(("ElectronSelectionTraining.Real.weighted." + label + ".root").c_str(), "RECREATE");
  TTree *RealEleOutputTree = new TTree("Electrons","Electrons");
  RealEleOutputTree->Branch("weight",&fWeight,"weight/F");
  RealEleOutputTree->Branch("run",&fRunNumber,"run/i");
  RealEleOutputTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  RealEleOutputTree->Branch("event",&fEventNumber,"event/i");
  RealEleOutputTree->Branch("pt",&fElePt,"pt/F"); 
  RealEleOutputTree->Branch("eta",&fEleEta,"eta/F"); 
  RealEleOutputTree->Branch("phi",&fElePhi,"phi/F"); 
  RealEleOutputTree->Branch("scet",&fEleSCEt,"scet/F"); 
  RealEleOutputTree->Branch("sceta",&fEleSCEta,"sceta/F"); 
  RealEleOutputTree->Branch("scphi",&fEleSCPhi,"scphi/F"); 
  RealEleOutputTree->Branch("pfiso",&fElePFIso,"pfiso/F"); 
  
  //CutBased Variables
  RealEleOutputTree->Branch("SigmaIEtaIEta",&fEleSigmaIEtaIEta,"SigmaIEtaIEta/F"); 
  RealEleOutputTree->Branch("DEtaIn",&fEleDEtaIn,"DEtaIn/F"); 
  RealEleOutputTree->Branch("DPhiIn",&fEleDPhiIn,"DPhiIn/F"); 
  RealEleOutputTree->Branch("HoverE",&fEleHoverE,"HoverE/F"); 
  RealEleOutputTree->Branch("D0",&fEleD0,"D0/F"); 
  RealEleOutputTree->Branch("DZ",&fEleDZ,"DZ/F"); 
  RealEleOutputTree->Branch("FBrem",&fEleFBrem,"FBrem/F"); 
  RealEleOutputTree->Branch("EOverP",&fEleEOverP,"EOverP/F"); 

  //Additional Vars used in Likelihood
  RealEleOutputTree->Branch("ESeedClusterOverPout",&fEleESeedClusterOverPout,"ESeedClusterOverPout/F"); 
  RealEleOutputTree->Branch("SigmaIPhiIPhi",&fEleSigmaIPhiIPhi,"SigmaIPhiIPhi/F"); 
  RealEleOutputTree->Branch("NBrem",&fEleNBrem,"NBrem/F"); 
  RealEleOutputTree->Branch("OneOverEMinusOneOverP",&fEleOneOverEMinusOneOverP,"OneOverEMinusOneOverP/F"); 
  RealEleOutputTree->Branch("ESeedClusterOverPIn",&fEleESeedClusterOverPIn,"ESeedClusterOverPIn/F"); 
  RealEleOutputTree->Branch("IP3d",&fEleIP3d,"IP3d/F"); 
  RealEleOutputTree->Branch("IP3dSig",&fEleIP3dSig,"IP3dSig/F"); 

  RealEleOutputTree->Branch("HcalDepth1OverEcal",&fEleHcalDepth1OverEcal,"HcalDepth1OverEcal/F"); 
  RealEleOutputTree->Branch("HcalDepth2OverEcal",&fEleHcalDepth2OverEcal,"HcalDepth2OverEcal/F"); 
  RealEleOutputTree->Branch("dEtaCalo",&fEledEtaCalo,"dEtaCalo/F"); 
  RealEleOutputTree->Branch("dPhiCalo",&fEledPhiCalo,"dPhiCalo/F"); 
  RealEleOutputTree->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
  RealEleOutputTree->Branch("CovIEtaIPhi",&fEleCovIEtaIPhi,"CovIEtaIPhi/F"); 
  RealEleOutputTree->Branch("SCEtaWidth",&fEleSCEtaWidth,"SCEtaWidth/F"); 
  RealEleOutputTree->Branch("SCPhiWidth",&fEleSCPhiWidth,"SCPhiWidth/F"); 
  RealEleOutputTree->Branch("GsfTrackChi2OverNdof",&fEleGsfTrackChi2OverNdof,"GsfTrackChi2OverNdof/F"); 
  RealEleOutputTree->Branch("R9",&fEleR9,"R9/F"); 

  RealEleOutputTree->Branch("SeedEMaxOverE",&fEleSeedEMaxOverE,"SeedEMaxOverE/F"); 
  RealEleOutputTree->Branch("SeedETopOverE",&fEleSeedETopOverE,"SeedETopOverE/F"); 
  RealEleOutputTree->Branch("SeedEBottomOverE",&fEleSeedEBottomOverE,"SeedEBottomOverE/F"); 
  RealEleOutputTree->Branch("SeedELeftOverE",&fEleSeedELeftOverE,"SeedELeftOverE/F"); 
  RealEleOutputTree->Branch("SeedERightOverE",&fEleSeedERightOverE,"SeedERightOverE/F"); 
  RealEleOutputTree->Branch("SeedE2ndOverE",&fEleSeedE2ndOverE,"SeedE2ndOverE/F"); 
  RealEleOutputTree->Branch("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE,"SeedE2x5RightOverE/F"); 
  RealEleOutputTree->Branch("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE,"SeedE2x5LeftOverE/F"); 
  RealEleOutputTree->Branch("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  RealEleOutputTree->Branch("SeedE2x5BottomOverE",&fEleSeedE2x5BottomOverE,"SeedE2x5BottomOverE/F"); 
  RealEleOutputTree->Branch("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  RealEleOutputTree->Branch("SeedE1x3OverE",&fEleSeedE1x3OverE,"SeedE1x3OverE/F"); 
  RealEleOutputTree->Branch("SeedE3x1OverE",&fEleSeedE3x1OverE,"SeedE3x1OverE/F"); 
  RealEleOutputTree->Branch("SeedE1x5OverE",&fEleSeedE1x5OverE,"SeedE1x5OverE/F"); 
  RealEleOutputTree->Branch("SeedE2x2OverE",&fEleSeedE2x2OverE,"SeedE2x2OverE/F"); 
  RealEleOutputTree->Branch("SeedE3x2OverE",&fEleSeedE3x2OverE,"SeedE3x2OverE/F"); 
  RealEleOutputTree->Branch("SeedE3x3OverE",&fEleSeedE3x3OverE,"SeedE3x3OverE/F"); 
  RealEleOutputTree->Branch("SeedE4x4OverE",&fEleSeedE4x4OverE,"SeedE4x4OverE/F"); 
  RealEleOutputTree->Branch("SeedE5x5OverE",&fEleSeedE5x5OverE,"SeedE5x5OverE/F"); 

  //Isolation Variables
  RealEleOutputTree->Branch("StandardLikelihood",&fEleStandardLikelihood,"StandardLikelihood/F"); 
  RealEleOutputTree->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
  RealEleOutputTree->Branch("ChargedIso03",&fEleChargedIso03,"ChargedIso03/F"); 
  RealEleOutputTree->Branch("NeutralHadronIso03",&fEleNeutralHadronIso03,"NeutralHadronIso03/F"); 
  RealEleOutputTree->Branch("GammaIso03",&fEleGammaIso03,"GammaIso03/F"); 
  RealEleOutputTree->Branch("ChargedIso04",&fEleChargedIso04,"ChargedIso04/F"); 
  RealEleOutputTree->Branch("NeutralHadronIso04",&fEleNeutralHadronIso04,"NeutralHadronIso04/F"); 
  RealEleOutputTree->Branch("GammaIso04",&fEleGammaIso04,"GammaIso04/F"); 
  RealEleOutputTree->Branch("ChargedIso04FromOtherVertices",&fEleChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  RealEleOutputTree->Branch("NeutralHadronIso04_10Threshold",&fEleNeutralHadronIso04_10Threshold,"NeutralHadronIso04_10Threshold/F"); 
  RealEleOutputTree->Branch("GammaIso04_10Threshold",&fEleGammaIso04_10Threshold,"GammaIso04_10Threshold/F"); 
  RealEleOutputTree->Branch("TrkIso03",&fEleTrkIso03,"TrkIso03/F"); 
  RealEleOutputTree->Branch("EMIso03",&fEleEMIso03,"EMIso03/F"); 
  RealEleOutputTree->Branch("HadIso03",&fEleHadIso03,"HadIso03/F"); 
  RealEleOutputTree->Branch("TrkIso04",&fEleTrkIso04,"TrkIso04/F"); 
  RealEleOutputTree->Branch("EMIso04",&fEleEMIso04,"EMIso04/F"); 
  RealEleOutputTree->Branch("HadIso04",&fEleHadIso04,"HadIso04/F"); 
  RealEleOutputTree->Branch("Rho",&fRho,"Rho/F"); 
  RealEleOutputTree->Branch("NVertices",&fNVertices,"NVertices/F"); 

  
  //*****************************************************************************************
  //RealEleTree
  //*****************************************************************************************
  TFile *RealEleFile = new TFile("ElectronSelectionTraining.Real.weighted.root", "READ");
  TTree *RealEleTree = (TTree*)RealEleFile->Get("Electrons");
  RealEleTree->SetBranchAddress( "weight", &fWeight);
  RealEleTree->SetBranchAddress( "run", &fRunNumber);
  RealEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  RealEleTree->SetBranchAddress( "event", &fEventNumber);
  RealEleTree->SetBranchAddress( "pt", &fElePt); 
  RealEleTree->SetBranchAddress( "eta", &fEleEta); 
  RealEleTree->SetBranchAddress( "phi", &fElePhi); 
  RealEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  RealEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  RealEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  RealEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  RealEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  RealEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  RealEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  RealEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  RealEleTree->SetBranchAddress( "D0", &fEleD0); 
  RealEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  RealEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  RealEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  RealEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  RealEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  RealEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  RealEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  RealEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  RealEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  RealEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig);
  RealEleTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  RealEleTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  RealEleTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  RealEleTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  RealEleTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  RealEleTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  RealEleTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  RealEleTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  RealEleTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  RealEleTree->SetBranchAddress( "R9", &fEleR9); 
  RealEleTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  RealEleTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  RealEleTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  RealEleTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  RealEleTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  RealEleTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  RealEleTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  RealEleTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  RealEleTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  RealEleTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  RealEleTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  RealEleTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  RealEleTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  RealEleTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  RealEleTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  RealEleTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  RealEleTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  RealEleTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  RealEleTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  RealEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  RealEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  RealEleTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  RealEleTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  RealEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  RealEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  RealEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  RealEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  RealEleTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  RealEleTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  RealEleTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  RealEleTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  RealEleTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  RealEleTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  RealEleTree->SetBranchAddress( "Rho", &fRho); 
  RealEleTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < RealEleTree->GetEntries(); ientry++) {       	
    RealEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Electron " << ientry << endl;

    Bool_t passCuts = kFALSE;
    
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

    if (option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);

    if (option == 221) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0 && fElePt <= 15);
    if (option == 222) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0 && fElePt > 15);

    if (passCuts) {
      RealEleOutputTree->Fill();
    }
  } 
  RealEleOutputFile->Write();
  RealEleOutputFile->Close();
  delete RealEleFile;

  //*****************************************************************************************
  //FakeEleTree Output
  //*****************************************************************************************
  TFile *FakeEleOutputFile = new TFile(("ElectronSelectionTraining.Fake.PtAndPUWeighted." + label + ".root").c_str(), "RECREATE");
  TTree *FakeEleOutputTree = new TTree("Electrons","Electrons");
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
  //FakeEleTargetTree
  //*****************************************************************************************
  TFile *FakeEleFile = new TFile("ElectronSelectionTraining.Fake.PtAndPUWeighted.root", "READ");
  TTree *FakeEleTree = (TTree*)FakeEleFile->Get("Electrons");
  FakeEleTree->SetBranchAddress( "weight", &fWeight);
  FakeEleTree->SetBranchAddress( "run", &fRunNumber);
  FakeEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  FakeEleTree->SetBranchAddress( "event", &fEventNumber);
  FakeEleTree->SetBranchAddress( "pt", &fElePt); 
  FakeEleTree->SetBranchAddress( "eta", &fEleEta); 
  FakeEleTree->SetBranchAddress( "phi", &fElePhi); 
  FakeEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  FakeEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  FakeEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  FakeEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  FakeEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  FakeEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  FakeEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  FakeEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  FakeEleTree->SetBranchAddress( "D0", &fEleD0); 
  FakeEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  FakeEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  FakeEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  FakeEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  FakeEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  FakeEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  FakeEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  FakeEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  FakeEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  FakeEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  FakeEleTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  FakeEleTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  FakeEleTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  FakeEleTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  FakeEleTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  FakeEleTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  FakeEleTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  FakeEleTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  FakeEleTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  FakeEleTree->SetBranchAddress( "R9", &fEleR9); 
  FakeEleTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  FakeEleTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  FakeEleTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  FakeEleTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  FakeEleTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  FakeEleTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  FakeEleTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  FakeEleTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  FakeEleTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  FakeEleTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  FakeEleTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  FakeEleTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  FakeEleTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  FakeEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  FakeEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  FakeEleTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  FakeEleTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  FakeEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  FakeEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  FakeEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  FakeEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  FakeEleTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  FakeEleTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  FakeEleTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  FakeEleTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  FakeEleTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  FakeEleTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  FakeEleTree->SetBranchAddress( "Rho", &fRho); 
  FakeEleTree->SetBranchAddress( "NVertices", &fNVertices); 

  for(UInt_t ientry=0; ientry < FakeEleTree->GetEntries(); ientry++) {       	
    FakeEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Electron " << ientry << endl;

    Bool_t passCuts = kFALSE;
    
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = -1;
    if (fElePt > 10.0) ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

    if (option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);

    if (option == 221) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0 && fElePt <= 15);
    if (option == 222) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0 && fElePt > 15);



    if (passCuts) {
      FakeEleOutputTree->Fill();
    }
    
  } //loop over electrons
  FakeEleOutputFile->Write();
  FakeEleOutputFile->Close();
  delete FakeEleFile;




  gBenchmark->Show("WWTemplate");       
} 



