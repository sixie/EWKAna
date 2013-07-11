//root -l /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/CompareElectronDistributions.C+'("/home/sixie/CMSSW_analysis/src/ElectronMVA/output/ElectronNtuple.Real.Subdet0LowPt.root","/home/sixie/CMSSW_analysis/src/ElectronSelectionTraining.HWW115.root","Subdet0LowPt",0)'


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
void CompareElectronDistributions(string ElectronFile1, string ElectronFile2, 
                                  string LegendLabel1, string LegendLabel2,
                                  string Label, Int_t Option)
{  

  string label = "";
  if (Label != "") label = "_" + Label;


  //*****************************************************************************************
  //Plotting Setup
  //*****************************************************************************************
  vector<Int_t> markers;
  vector<Int_t> colors;
  colors.push_back(kRed);     markers.push_back(20);
  colors.push_back(kBlue);    markers.push_back(21);
  colors.push_back(kMagenta); markers.push_back(22);
  colors.push_back(kCyan);    markers.push_back(34);
  colors.push_back(kBlack);   markers.push_back(29);
  colors.push_back(kGreen);   markers.push_back(33);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  

  //Input Variables
  TH1F *EleFBrem_Data = new TH1F(("EleFBrem_Data"+label).c_str(), "; FBrem ; Number of Events ",  200, -2 , 2);
  TH1F *EleDEta_Data = new TH1F(("EleDEta_Data"+label).c_str(), "; DEta ; Number of Events ",  200, -0.03 , 0.03);
  TH1F *EleDPhi_Data = new TH1F(("EleDPhi_Data"+label).c_str(), "; DPhi ; Number of Events ",  200, -0.3 , 0.3);
  TH1F *EleSigmaIEtaIEta_Data = new TH1F(("EleSigmaIEtaIEta_Data"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  200, 0 , 0.04);
  TH1F *EleSigmaIPhiIPhi_Data = new TH1F(("EleSigmaIPhiIPhi_Data"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  200, 0 , 0.1);
//   TH1F *EleOneOverEMinusOneOverP_Data = new TH1F(("EleOneOverEMinusOneOverP_Data"+label).c_str(), "; OneOverEMinusOneOverP ; Number of Events ",  200, -0.5 , 0.2);
  TH1F *EleOneOverEMinusOneOverP_Data = new TH1F(("EleOneOverEMinusOneOverP_Data"+label).c_str(), "; OneOverEMinusOneOverP ; Number of Events ",  200, -0.1 , 0.05);

  TH1F *EleNBrem_Data = new TH1F(("EleNBrem_Data"+label).c_str(), "; NBrem ; Number of Events ",  16, -0.5 , 15.5);
  TH1F *EleHoverE_Data = new TH1F(("EleHoverE_Data"+label).c_str(), "; HoverE ; Number of Events ",  200, -1 , 9.0);
  TH1F *EleD0_Data = new TH1F(("EleD0_Data"+label).c_str(), "; D0 ; Number of Events ",  200, -0.1 , 0.1);
  TH1F *EleEOverP_Data = new TH1F(("EleEOverP_Data"+label).c_str(), "; EOverP ; Number of Events ",  200, -5 , 5);
  TH1F *EleESeedClusterOverPout_Data = new TH1F(("EleESeedClusterOverPout_Data"+label).c_str(), "; ESeedClusterOverPout ; Number of Events ",  200, -5 , 5);
  TH1F *EleESeedClusterOverPIn_Data = new TH1F(("EleESeedClusterOverPIn_Data"+label).c_str(), "; ESeedClusterOverPIn ; Number of Events ",  200, -5 , 5);
  TH1F *EleIP3d_Data = new TH1F(("EleIP3d_Data"+label).c_str(), "; IP3d ; Number of Events ",  200, -0.1 , 0.1);
  TH1F *EleIP3dSig_Data = new TH1F(("EleIP3dSig_Data"+label).c_str(), "; IP3dSig ; Number of Events ",  200, -10, 10);


  TH1F *EleFBrem_MC = new TH1F(("EleFBrem_MC"+label).c_str(), "; FBrem ; Number of Events ",  200, -2 , 2);
  TH1F *EleDEta_MC = new TH1F(("EleDEta_MC"+label).c_str(), "; DEta ; Number of Events ",  200, -0.03 , 0.03);
  TH1F *EleDPhi_MC = new TH1F(("EleDPhi_MC"+label).c_str(), "; DPhi ; Number of Events ",  200, -0.3 , 0.3);
  TH1F *EleSigmaIEtaIEta_MC = new TH1F(("EleSigmaIEtaIEta_MC"+label).c_str(), "; SigmaIEtaIEta ; Number of Events ",  200, 0 , 0.04);
  TH1F *EleSigmaIPhiIPhi_MC = new TH1F(("EleSigmaIPhiIPhi_MC"+label).c_str(), "; SigmaIPhiIPhi ; Number of Events ",  200, 0 , 0.1);
//   TH1F *EleOneOverEMinusOneOverP_MC = new TH1F(("EleOneOverEMinusOneOverP_MC"+label).c_str(), "; OneOverEMinusOneOverP ; Number of Events ",  200, -0.5 , 0.2);

  TFile *fff = new TFile("/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronLikelihoodPdfs_MC.root","READ");
  TH1F *EleOneOverEMinusOneOverP_MC = (TH1F*)fff ->Get("OneOverEMinusOneOverPClass_hadrons_subdet0_ptbin0_class1;");


  TH1F *EleNBrem_MC = new TH1F(("EleNBrem_MC"+label).c_str(), "; NBrem ; Number of Events ",  16, -0.5 , 15.5);
  TH1F *EleHoverE_MC = new TH1F(("EleHoverE_MC"+label).c_str(), "; HoverE ; Number of Events ",  200, -1 , 9.0);
  TH1F *EleD0_MC = new TH1F(("EleD0_MC"+label).c_str(), "; D0 ; Number of Events ",  200, -0.1 , 0.1);
  TH1F *EleEOverP_MC = new TH1F(("EleEOverP_MC"+label).c_str(), "; EOverP ; Number of Events ",  200, -5 , 5);
  TH1F *EleESeedClusterOverPout_MC = new TH1F(("EleESeedClusterOverPout_MC"+label).c_str(), "; ESeedClusterOverPout ; Number of Events ",  200, -5 , 5);
  TH1F *EleESeedClusterOverPIn_MC = new TH1F(("EleESeedClusterOverPIn_MC"+label).c_str(), "; ESeedClusterOverPIn ; Number of Events ",  200, -5 , 5);
  TH1F *EleIP3d_MC = new TH1F(("EleIP3d_MC"+label).c_str(), "; IP3d ; Number of Events ",  200, -0.1 , 0.1);
  TH1F *EleIP3dSig_MC = new TH1F(("EleIP3dSig_MC"+label).c_str(), "; IP3dSig ; Number of Events ",  200, -10, 10);


  //MVA Output Variables
  TH1F *EleStandardLikelihood_Data = new TH1F(("EleStandardLikelihood_Data"+label).c_str(), "; StandardLikelihood ; Number of Events ",  200, -20 , 20);
  TH1F *ElePFMVA_Data = new TH1F(("ElePFMVA_Data"+label).c_str(), "; PFMVA ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDT_Data = new TH1F(("EleBDT_Data"+label).c_str(), "; BDT ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTG_Data = new TH1F(("EleBDTG_Data"+label).c_str(), "; BDTG ; Number of Events ",  200, -2 , 2);
  TH1F *EleTMVALikelihood_Data = new TH1F(("EleTMVALikelihood_Data"+label).c_str(), "; TMVALikelihood ; Number of Events ",  200, 0 , 1);
  TH1F *EleTMVALikelihoodD_Data = new TH1F(("EleTMVALikelihoodD_Data"+label).c_str(), "; TMVALikelihoodD ; Number of Events ",  200, -2 , 2);
  TH1F *EleNN_Data = new TH1F(("EleNN_Data"+label).c_str(), "; NN ; Number of Events ",  200, -2 , 2);

  TH1F *EleStandardLikelihood_MC = new TH1F(("EleStandardLikelihood_MC"+label).c_str(), "; StandardLikelihood ; Number of Events ",  200, -20 , 20);
  TH1F *ElePFMVA_MC = new TH1F(("ElePFMVA_MC"+label).c_str(), "; PFMVA ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDT_MC = new TH1F(("EleBDT_MC"+label).c_str(), "; BDT ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTG_MC = new TH1F(("EleBDTG_MC"+label).c_str(), "; BDTG ; Number of Events ",  200, -2 , 2);
  TH1F *EleTMVALikelihood_MC = new TH1F(("EleTMVALikelihood_MC"+label).c_str(), "; TMVALikelihood ; Number of Events ",  200, 0 , 1);
  TH1F *EleTMVALikelihoodD_MC = new TH1F(("EleTMVALikelihoodD_MC"+label).c_str(), "; TMVALikelihoodD ; Number of Events ",  200, -2 , 2);
  TH1F *EleNN_MC = new TH1F(("EleNN_MC"+label).c_str(), "; NN ; Number of Events ",  200, -2 , 2);


  TH1F *EleBDTGV1_Data = new TH1F(("EleBDTGV1_Data"+label).c_str(), "; BDTG V1 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV2_Data = new TH1F(("EleBDTGV2_Data"+label).c_str(), "; BDTG V2 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV3_Data = new TH1F(("EleBDTGV3_Data"+label).c_str(), "; BDTG V3 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV4_Data = new TH1F(("EleBDTGV4_Data"+label).c_str(), "; BDTG V4 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV5_Data = new TH1F(("EleBDTGV5_Data"+label).c_str(), "; BDTG V5 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV6_Data = new TH1F(("EleBDTGV6_Data"+label).c_str(), "; BDTG V6 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV7_Data = new TH1F(("EleBDTGV7_Data"+label).c_str(), "; BDTG V7 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV8_Data = new TH1F(("EleBDTGV8_Data"+label).c_str(), "; BDTG V8 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV9_Data = new TH1F(("EleBDTGV9_Data"+label).c_str(), "; BDTG V9 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV1_MC = new TH1F(("EleBDTGV1_MC"+label).c_str(), "; BDTG V1 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV2_MC = new TH1F(("EleBDTGV2_MC"+label).c_str(), "; BDTG V2 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV3_MC = new TH1F(("EleBDTGV3_MC"+label).c_str(), "; BDTG V3 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV4_MC = new TH1F(("EleBDTGV4_MC"+label).c_str(), "; BDTG V4 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV5_MC = new TH1F(("EleBDTGV5_MC"+label).c_str(), "; BDTG V5 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV6_MC = new TH1F(("EleBDTGV6_MC"+label).c_str(), "; BDTG V6 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV7_MC = new TH1F(("EleBDTGV7_MC"+label).c_str(), "; BDTG V7 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV8_MC = new TH1F(("EleBDTGV8_MC"+label).c_str(), "; BDTG V8 ; Number of Events ",  200, -2 , 2);
  TH1F *EleBDTGV9_MC = new TH1F(("EleBDTGV9_MC"+label).c_str(), "; BDTG V9 ; Number of Events ",  200, -2 , 2);

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

  //Isolation Variables
  Float_t                 fEleStandardLikelihood; 
  Float_t                 fElePFMVA; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  Float_t                 fEleBDT;
  Float_t                 fEleBDTG;
  Float_t                 fEleNN;
  Float_t                 fEleTMVALikelihood;
  Float_t                 fEleTMVALikelihoodD;

  Float_t                 fEleBDTGV1;
  Float_t                 fEleBDTGV2;
  Float_t                 fEleBDTGV3;
  Float_t                 fEleBDTGV4;
  Float_t                 fEleBDTGV5;
  Float_t                 fEleBDTGV6;
  Float_t                 fEleBDTGV7;
  Float_t                 fEleBDTGV8;
  Float_t                 fEleBDTGV9;

 
  //*****************************************************************************************
  //DataEleTree
  //*****************************************************************************************
  TFile *DataEleFile = new TFile(ElectronFile1.c_str(), "READ");
  TTree *DataEleTree = (TTree*)DataEleFile->Get("Electrons");
  DataEleTree->SetBranchAddress( "weight", &fWeight);
  DataEleTree->SetBranchAddress( "run", &fRunNumber);
  DataEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  DataEleTree->SetBranchAddress( "event", &fEventNumber);
  DataEleTree->SetBranchAddress( "pt", &fElePt); 
  DataEleTree->SetBranchAddress( "eta", &fEleEta); 
  DataEleTree->SetBranchAddress( "phi", &fElePhi); 
  DataEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  DataEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  DataEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  DataEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  DataEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  DataEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  DataEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  DataEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  DataEleTree->SetBranchAddress( "D0", &fEleD0); 
  DataEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  DataEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  DataEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  DataEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  DataEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  DataEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  DataEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  DataEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  DataEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  DataEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  DataEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  DataEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  DataEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  DataEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  DataEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  DataEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  DataEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  DataEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  DataEleTree->SetBranchAddress( "Rho", &fRho); 
  DataEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  DataEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  DataEleTree->SetBranchAddress( "BDTGV0", &fEleBDTG); 
  DataEleTree->SetBranchAddress( "NN", &fEleNN); 
  DataEleTree->SetBranchAddress( "LikelihoodV0", &fEleTMVALikelihood); 
  DataEleTree->SetBranchAddress( "LikelihoodD", &fEleTMVALikelihoodD); 
  DataEleTree->SetBranchAddress( "BDTGV1a", &fEleBDTGV1); 
  DataEleTree->SetBranchAddress( "BDTGV2a", &fEleBDTGV2); 
  DataEleTree->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  DataEleTree->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  DataEleTree->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  DataEleTree->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  DataEleTree->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  DataEleTree->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  DataEleTree->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 

  for(UInt_t ientry=0; ientry < DataEleTree->GetEntries(); ientry++) {       	
    DataEleTree->GetEntry(ientry);
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    if (!(fNVertices >= 3 && fNVertices <= 6)) continue;

    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 15.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply full isolation cut
    Bool_t passIsoCuts = kFALSE;
    if (ptBin == 0) passIsoCuts = ( fElePFIso  < 0.09 ); 
    if (ptBin == 1) passIsoCuts = ( fElePFIso  < 0.13 ); 
    if (!passIsoCuts) continue;

    //Fill Histograms
    EleFBrem_Data->Fill(fEleFBrem,fWeight);
    EleDEta_Data->Fill(fEleDEtaIn,fWeight);
    EleDPhi_Data->Fill(fEleDPhiIn,fWeight);
    EleSigmaIEtaIEta_Data->Fill(fEleSigmaIEtaIEta,fWeight);
    EleSigmaIPhiIPhi_Data->Fill(fEleSigmaIPhiIPhi,fWeight);
    EleOneOverEMinusOneOverP_Data->Fill(fEleOneOverEMinusOneOverP,fWeight);
    EleNBrem_Data->Fill(fEleNBrem,fWeight);

    EleStandardLikelihood_Data->Fill(fEleStandardLikelihood,fWeight);
    ElePFMVA_Data->Fill(fElePFMVA,fWeight);
    EleBDT_Data->Fill(fEleBDT,fWeight);
    EleBDTG_Data->Fill(fEleBDTG,fWeight);
    EleTMVALikelihood_Data->Fill(fEleTMVALikelihood,fWeight);
    EleTMVALikelihoodD_Data->Fill(fEleTMVALikelihoodD,fWeight);
    EleNN_Data->Fill(fEleNN,fWeight);

    EleBDTGV1_Data->Fill(fEleBDTGV1,fWeight);
    EleBDTGV2_Data->Fill(fEleBDTGV2,fWeight);
    EleBDTGV3_Data->Fill(fEleBDTGV3,fWeight);
    EleBDTGV4_Data->Fill(fEleBDTGV4,fWeight);
    EleBDTGV5_Data->Fill(fEleBDTGV5,fWeight);
    EleBDTGV6_Data->Fill(fEleBDTGV6,fWeight);
    EleBDTGV7_Data->Fill(fEleBDTGV7,fWeight);
    EleBDTGV8_Data->Fill(fEleBDTGV8,fWeight);
    EleBDTGV9_Data->Fill(fEleBDTGV9,fWeight);

  } 
  





  //*****************************************************************************************
  //MCEleTree
  //*****************************************************************************************
  TFile *MCEleFile = new TFile(ElectronFile2.c_str(), "READ");
  TTree *MCEleTree = (TTree*)MCEleFile->Get("Electrons");
  MCEleTree->SetBranchAddress( "weight", &fWeight);
  MCEleTree->SetBranchAddress( "run", &fRunNumber);
  MCEleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
  MCEleTree->SetBranchAddress( "event", &fEventNumber);
  MCEleTree->SetBranchAddress( "pt", &fElePt); 
  MCEleTree->SetBranchAddress( "eta", &fEleEta); 
  MCEleTree->SetBranchAddress( "phi", &fElePhi); 
  MCEleTree->SetBranchAddress( "scet", &fEleSCEt); 
  MCEleTree->SetBranchAddress( "sceta", &fEleSCEta); 
  MCEleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
  MCEleTree->SetBranchAddress( "pfiso", &fElePFIso); 
  MCEleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  MCEleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  MCEleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  MCEleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
  MCEleTree->SetBranchAddress( "D0", &fEleD0); 
  MCEleTree->SetBranchAddress( "DZ", &fEleDZ); 
  MCEleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
  MCEleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
  MCEleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  MCEleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  MCEleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
  MCEleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  MCEleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  MCEleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
  MCEleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  MCEleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  MCEleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
  MCEleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  MCEleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  MCEleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  MCEleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  MCEleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  MCEleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  MCEleTree->SetBranchAddress( "Rho", &fRho); 
  MCEleTree->SetBranchAddress( "NVertices", &fNVertices); 
  MCEleTree->SetBranchAddress( "BDT", &fEleBDT); 
  MCEleTree->SetBranchAddress( "BDTGV0", &fEleBDTG); 
  MCEleTree->SetBranchAddress( "NN", &fEleNN); 
  MCEleTree->SetBranchAddress( "LikelihoodV0", &fEleTMVALikelihood); 
  MCEleTree->SetBranchAddress( "LikelihoodD", &fEleTMVALikelihoodD); 
  MCEleTree->SetBranchAddress( "BDTGV1a", &fEleBDTGV1); 
  MCEleTree->SetBranchAddress( "BDTGV2a", &fEleBDTGV2); 
  MCEleTree->SetBranchAddress( "BDTGV3", &fEleBDTGV3); 
  MCEleTree->SetBranchAddress( "BDTGV4", &fEleBDTGV4); 
  MCEleTree->SetBranchAddress( "BDTGV5", &fEleBDTGV5); 
  MCEleTree->SetBranchAddress( "BDTGV6", &fEleBDTGV6); 
  MCEleTree->SetBranchAddress( "BDTGV7", &fEleBDTGV7); 
  MCEleTree->SetBranchAddress( "BDTGV8", &fEleBDTGV8); 
  MCEleTree->SetBranchAddress( "BDTGV9", &fEleBDTGV9); 


  for(UInt_t ientry=0; ientry < MCEleTree->GetEntries(); ientry++) {       	
    MCEleTree->GetEntry(ientry);
    if (!(fNVertices >= 3 && fNVertices <= 6)) continue;
    
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
    
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 15.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);

    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    //apply full isolation cut
    Bool_t passIsoCuts = kFALSE;
    if (ptBin == 0) passIsoCuts = ( fElePFIso  < 0.09 ); 
    if (ptBin == 1) passIsoCuts = ( fElePFIso  < 0.13 ); 
    if (!passIsoCuts) continue;

    //Fill Histograms
    EleFBrem_MC->Fill(fEleFBrem,fWeight);
    EleDEta_MC->Fill(fEleDEtaIn,fWeight);
    EleDPhi_MC->Fill(fEleDPhiIn,fWeight);
    EleSigmaIEtaIEta_MC->Fill(fEleSigmaIEtaIEta,fWeight);
    EleSigmaIPhiIPhi_MC->Fill(fEleSigmaIPhiIPhi,fWeight);
    EleOneOverEMinusOneOverP_MC->Fill(fEleOneOverEMinusOneOverP,fWeight);
    EleNBrem_MC->Fill(fEleNBrem,fWeight);

    EleStandardLikelihood_MC->Fill(fEleStandardLikelihood,fWeight);
    ElePFMVA_MC->Fill(fElePFMVA,fWeight);
    EleBDT_MC->Fill(fEleBDT,fWeight);
    EleBDTG_MC->Fill(fEleBDTG,fWeight);
    EleTMVALikelihood_MC->Fill(fEleTMVALikelihood,fWeight);
    EleTMVALikelihoodD_MC->Fill(fEleTMVALikelihoodD,fWeight);
    EleNN_MC->Fill(fEleNN,fWeight);
    EleBDTGV1_MC->Fill(fEleBDTGV1,fWeight);
    EleBDTGV2_MC->Fill(fEleBDTGV2,fWeight);
    EleBDTGV3_MC->Fill(fEleBDTGV3,fWeight);
    EleBDTGV4_MC->Fill(fEleBDTGV4,fWeight);
    EleBDTGV5_MC->Fill(fEleBDTGV5,fWeight);
    EleBDTGV6_MC->Fill(fEleBDTGV6,fWeight);
    EleBDTGV7_MC->Fill(fEleBDTGV7,fWeight);
    EleBDTGV8_MC->Fill(fEleBDTGV8,fWeight);
    EleBDTGV9_MC->Fill(fEleBDTGV9,fWeight);

  } //loop over electrons
  



  //*****************************************************************************************
  //Normalize
  //*****************************************************************************************
  NormalizeHist(EleFBrem_Data);
  NormalizeHist(EleDEta_Data);
  NormalizeHist(EleDPhi_Data);
  NormalizeHist(EleSigmaIEtaIEta_Data);
  NormalizeHist(EleSigmaIPhiIPhi_Data);
  NormalizeHist(EleOneOverEMinusOneOverP_Data);
  NormalizeHist(EleNBrem_Data);

  NormalizeHist(EleFBrem_MC);
  NormalizeHist(EleDEta_MC);
  NormalizeHist(EleDPhi_MC);
  NormalizeHist(EleSigmaIEtaIEta_MC);
  NormalizeHist(EleSigmaIPhiIPhi_MC);
  NormalizeHist(EleOneOverEMinusOneOverP_MC);
  NormalizeHist(EleNBrem_MC);

  NormalizeHist(EleStandardLikelihood_Data);
  NormalizeHist(ElePFMVA_Data);
  NormalizeHist(EleBDT_Data);
  NormalizeHist(EleBDTG_Data);
  NormalizeHist(EleTMVALikelihood_Data);
  NormalizeHist(EleTMVALikelihoodD_Data);
  NormalizeHist(EleNN_Data);

  NormalizeHist(EleStandardLikelihood_MC);
  NormalizeHist(ElePFMVA_MC);
  NormalizeHist(EleBDT_MC);
  NormalizeHist(EleBDTG_MC);
  NormalizeHist(EleTMVALikelihood_MC);
  NormalizeHist(EleTMVALikelihoodD_MC);
  NormalizeHist(EleNN_MC);

  NormalizeHist(EleBDTGV1_Data);
  NormalizeHist(EleBDTGV2_Data);
  NormalizeHist(EleBDTGV3_Data);
  NormalizeHist(EleBDTGV4_Data);
  NormalizeHist(EleBDTGV5_Data);
  NormalizeHist(EleBDTGV6_Data);
  NormalizeHist(EleBDTGV7_Data);
  NormalizeHist(EleBDTGV8_Data);
  NormalizeHist(EleBDTGV9_Data);
  NormalizeHist(EleBDTGV1_MC);
  NormalizeHist(EleBDTGV2_MC);
  NormalizeHist(EleBDTGV3_MC);
  NormalizeHist(EleBDTGV4_MC);
  NormalizeHist(EleBDTGV5_MC);
  NormalizeHist(EleBDTGV6_MC);
  NormalizeHist(EleBDTGV7_MC);
  NormalizeHist(EleBDTGV8_MC);
  NormalizeHist(EleBDTGV9_MC);

  //*****************************************************************************************
  //Plot Distributions
  //*****************************************************************************************
  TCanvas *cv;
  TLegend *legend;

  vector<TH1F*> hists;
  vector<string> histLabels;

  //*****************************************************************************************
  //FBrem
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleFBrem_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleFBrem_MC, LegendLabel2.c_str(), "L");

  EleFBrem_Data->SetLineColor(kBlue);
  EleFBrem_MC->SetLineColor(kRed);   

  EleFBrem_Data->SetMaximum(TMath::Max(EleFBrem_Data->GetMaximum(),EleFBrem_MC->GetMaximum()) * 1.2);
  EleFBrem_Data->GetYaxis()->SetTitleOffset(1.2);
  EleFBrem_Data->GetXaxis()->SetTitleOffset(1.05);
  EleFBrem_Data->Draw("hist");
  EleFBrem_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleFBrem" + label + ".gif").c_str());



  //*****************************************************************************************
  //DEta
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleDEta_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleDEta_MC, LegendLabel2.c_str(), "L");

  EleDEta_Data->SetLineColor(kBlue);
  EleDEta_MC->SetLineColor(kRed);   

  EleDEta_Data->SetMaximum(TMath::Max(EleDEta_Data->GetMaximum(),EleDEta_MC->GetMaximum()) * 1.2);
  EleDEta_Data->GetYaxis()->SetTitleOffset(1.2);
  EleDEta_Data->GetXaxis()->SetTitleOffset(1.05);
  EleDEta_Data->Draw("hist");
  EleDEta_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleDEta" + label + ".gif").c_str());



  //*****************************************************************************************
  //DPhi
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleDPhi_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleDPhi_MC, LegendLabel2.c_str(), "L");

  EleDPhi_Data->SetLineColor(kBlue);
  EleDPhi_MC->SetLineColor(kRed);   

  EleDPhi_Data->SetMaximum(TMath::Max(EleDPhi_Data->GetMaximum(),EleDPhi_MC->GetMaximum()) * 1.2);
  EleDPhi_Data->GetYaxis()->SetTitleOffset(1.2);
  EleDPhi_Data->GetXaxis()->SetTitleOffset(1.05);
  EleDPhi_Data->Draw("hist");
  EleDPhi_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleDPhi" + label + ".gif").c_str());


  //*****************************************************************************************
  //SigmaIEtaIEta
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  if(Option == 2 || Option == 5) legend = new TLegend(0.20,0.70,0.40,0.90);
  else legend = new TLegend(0.60,0.70,0.80,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleSigmaIEtaIEta_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleSigmaIEtaIEta_MC, LegendLabel2.c_str(), "L");

  EleSigmaIEtaIEta_Data->SetLineColor(kBlue);
  EleSigmaIEtaIEta_MC->SetLineColor(kRed);   

  EleSigmaIEtaIEta_Data->SetMaximum(TMath::Max(EleSigmaIEtaIEta_Data->GetMaximum(),EleSigmaIEtaIEta_MC->GetMaximum()) * 1.2);
  EleSigmaIEtaIEta_Data->GetYaxis()->SetTitleOffset(1.2);
  EleSigmaIEtaIEta_Data->GetXaxis()->SetTitleOffset(1.05);
  EleSigmaIEtaIEta_Data->Draw("hist");
  EleSigmaIEtaIEta_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleSigmaIEtaIEta" + label + ".gif").c_str());


  //*****************************************************************************************
  //SigmaIPhiIPhi
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.70,0.75,0.90,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleSigmaIPhiIPhi_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleSigmaIPhiIPhi_MC, LegendLabel2.c_str(), "L");

  EleSigmaIPhiIPhi_Data->SetLineColor(kBlue);
  EleSigmaIPhiIPhi_MC->SetLineColor(kRed);   

  EleSigmaIPhiIPhi_Data->SetMaximum(TMath::Max(EleSigmaIPhiIPhi_Data->GetMaximum(),EleSigmaIPhiIPhi_MC->GetMaximum()) * 1.2);
  EleSigmaIPhiIPhi_Data->GetYaxis()->SetTitleOffset(1.2);
  EleSigmaIPhiIPhi_Data->GetXaxis()->SetTitleOffset(1.05);
  EleSigmaIPhiIPhi_Data->Draw("hist");
  EleSigmaIPhiIPhi_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleSigmaIPhiIPhi" + label + ".gif").c_str());


  //*****************************************************************************************
  //OneOverEMinusOneOverP
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
//   cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleOneOverEMinusOneOverP_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleOneOverEMinusOneOverP_MC, LegendLabel2.c_str(), "L");

  EleOneOverEMinusOneOverP_Data->SetLineColor(kBlue);
  EleOneOverEMinusOneOverP_MC->SetLineColor(kRed);   

  EleOneOverEMinusOneOverP_Data->SetMaximum(TMath::Max(EleOneOverEMinusOneOverP_Data->GetMaximum(),EleOneOverEMinusOneOverP_MC->GetMaximum()) * 1.2);
  EleOneOverEMinusOneOverP_Data->GetYaxis()->SetTitleOffset(1.2);
  EleOneOverEMinusOneOverP_Data->GetXaxis()->SetTitleOffset(1.05);
  EleOneOverEMinusOneOverP_Data->Draw("hist");
  EleOneOverEMinusOneOverP_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleOneOverEMinusOneOverP" + label + ".gif").c_str());


  //*****************************************************************************************
  //NBrem
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleNBrem_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleNBrem_MC, LegendLabel2.c_str(), "L");

  EleNBrem_Data->SetLineColor(kBlue);
  EleNBrem_MC->SetLineColor(kRed);   

  EleNBrem_Data->SetMaximum(TMath::Max(EleNBrem_Data->GetMaximum(),EleNBrem_MC->GetMaximum()) * 1.2);
  EleNBrem_Data->GetYaxis()->SetTitleOffset(1.2);
  EleNBrem_Data->GetXaxis()->SetTitleOffset(1.05);
  EleNBrem_Data->Draw("hist");
  EleNBrem_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleNBrem" + label + ".gif").c_str());





  


  //*****************************************************************************************
  //StandardLikelihood
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleStandardLikelihood_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleStandardLikelihood_MC, LegendLabel2.c_str(), "L");

  EleStandardLikelihood_Data->SetLineColor(kBlue);
  EleStandardLikelihood_MC->SetLineColor(kRed);   

  EleStandardLikelihood_Data->SetMaximum(TMath::Max(EleStandardLikelihood_Data->GetMaximum(),EleStandardLikelihood_MC->GetMaximum()) * 1.2);
  EleStandardLikelihood_Data->GetYaxis()->SetTitleOffset(1.2);
  EleStandardLikelihood_Data->GetXaxis()->SetTitleOffset(1.05);
  EleStandardLikelihood_Data->Draw("hist");
  EleStandardLikelihood_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleStandardLikelihood" + label + ".gif").c_str());

  //*****************************************************************************************
  //PFMVA
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(ElePFMVA_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(ElePFMVA_MC, LegendLabel2.c_str(), "L");

  ElePFMVA_Data->SetLineColor(kBlue);
  ElePFMVA_MC->SetLineColor(kRed);   

  ElePFMVA_Data->SetMaximum(TMath::Max(ElePFMVA_Data->GetMaximum(),ElePFMVA_MC->GetMaximum()) * 1.2);
  ElePFMVA_Data->GetYaxis()->SetTitleOffset(1.2);
  ElePFMVA_Data->GetXaxis()->SetTitleOffset(1.05);
  ElePFMVA_Data->Draw("hist");
  ElePFMVA_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("ElePFMVA" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDT
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDT_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDT_MC, LegendLabel2.c_str(), "L");

  EleBDT_Data->SetLineColor(kBlue);
  EleBDT_MC->SetLineColor(kRed);   

  EleBDT_Data->SetMaximum(TMath::Max(EleBDT_Data->GetMaximum(),EleBDT_MC->GetMaximum()) * 1.2);
  EleBDT_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDT_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDT_Data->Draw("hist");
  EleBDT_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDT" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  if (Option >= 3) cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTG_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTG_MC, LegendLabel2.c_str(), "L");

  EleBDTG_Data->SetLineColor(kBlue);
  EleBDTG_MC->SetLineColor(kRed);   

  EleBDTG_Data->SetMaximum(TMath::Max(EleBDTG_Data->GetMaximum(),EleBDTG_MC->GetMaximum()) * 1.2);
  EleBDTG_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTG_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTG_Data->Draw("hist");
  EleBDTG_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTG" + label + ".gif").c_str());


  //*****************************************************************************************
  //TMVALikelihood
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleTMVALikelihood_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleTMVALikelihood_MC, LegendLabel2.c_str(), "L");

  EleTMVALikelihood_Data->SetLineColor(kBlue);
  EleTMVALikelihood_MC->SetLineColor(kRed);   

  EleTMVALikelihood_Data->SetMaximum(TMath::Max(EleTMVALikelihood_Data->GetMaximum(),EleTMVALikelihood_MC->GetMaximum()) * 1.2);
  EleTMVALikelihood_Data->GetYaxis()->SetTitleOffset(1.2);
  EleTMVALikelihood_Data->GetXaxis()->SetTitleOffset(1.05);
  EleTMVALikelihood_Data->Draw("hist");
  EleTMVALikelihood_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleTMVALikelihood" + label + ".gif").c_str());


  //*****************************************************************************************
  //TMVALikelihoodD
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleTMVALikelihoodD_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleTMVALikelihoodD_MC, LegendLabel2.c_str(), "L");

  EleTMVALikelihoodD_Data->SetLineColor(kBlue);
  EleTMVALikelihoodD_MC->SetLineColor(kRed);   

  EleTMVALikelihoodD_Data->SetMaximum(TMath::Max(EleTMVALikelihoodD_Data->GetMaximum(),EleTMVALikelihoodD_MC->GetMaximum()) * 1.2);
  EleTMVALikelihoodD_Data->GetYaxis()->SetTitleOffset(1.2);
  EleTMVALikelihoodD_Data->GetXaxis()->SetTitleOffset(1.05);
  EleTMVALikelihoodD_Data->Draw("hist");
  EleTMVALikelihoodD_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleTMVALikelihoodD" + label + ".gif").c_str());


  //*****************************************************************************************
  //NN
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleNN_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleNN_MC, LegendLabel2.c_str(), "L");

  EleNN_Data->SetLineColor(kBlue);
  EleNN_MC->SetLineColor(kRed);   

  EleNN_Data->SetMaximum(TMath::Max(EleNN_Data->GetMaximum(),EleNN_MC->GetMaximum()) * 1.2);
  EleNN_Data->GetYaxis()->SetTitleOffset(1.2);
  EleNN_Data->GetXaxis()->SetTitleOffset(1.05);
  EleNN_Data->Draw("hist");
  EleNN_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleNN" + label + ".gif").c_str());





  //*****************************************************************************************
  //BDTG V1
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  if (Option >= 3) cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV1_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV1_MC, LegendLabel2.c_str(), "L");

  EleBDTGV1_Data->SetLineColor(kBlue);
  EleBDTGV1_MC->SetLineColor(kRed);   

  EleBDTGV1_Data->GetXaxis()->SetTitle("BDTG Without IP Info");
  EleBDTGV1_Data->SetMaximum(TMath::Max(EleBDTGV1_Data->GetMaximum(),EleBDTGV1_MC->GetMaximum()) * 1.2);
  EleBDTGV1_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV1_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV1_Data->Draw("hist");
  EleBDTGV1_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV1" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V2
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  if (Option >= 3) cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV2_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV2_MC, LegendLabel2.c_str(), "L");

  EleBDTGV2_Data->SetLineColor(kBlue);
  EleBDTGV2_MC->SetLineColor(kRed);   

  EleBDTGV2_Data->SetMaximum(TMath::Max(EleBDTGV2_Data->GetMaximum(),EleBDTGV2_MC->GetMaximum()) * 1.2);
  EleBDTGV2_Data->GetXaxis()->SetTitle("BDTG With IP Info");
  EleBDTGV2_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV2_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV2_Data->Draw("hist");
  EleBDTGV2_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV2" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V3
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  if (Option >= 3) cv->SetLogy();
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV3_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV3_MC, LegendLabel2.c_str(), "L");

  EleBDTGV3_Data->SetLineColor(kBlue);
  EleBDTGV3_MC->SetLineColor(kRed);   

  EleBDTGV3_Data->SetMaximum(TMath::Max(EleBDTGV3_Data->GetMaximum(),EleBDTGV3_MC->GetMaximum()) * 1.2);
  EleBDTGV3_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV3_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV3_Data->Draw("hist");
  EleBDTGV3_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV3" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V4
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV4_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV4_MC, LegendLabel2.c_str(), "L");

  EleBDTGV4_Data->SetLineColor(kBlue);
  EleBDTGV4_MC->SetLineColor(kRed);   

  EleBDTGV4_Data->SetMaximum(TMath::Max(EleBDTGV4_Data->GetMaximum(),EleBDTGV4_MC->GetMaximum()) * 1.2);
  EleBDTGV4_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV4_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV4_Data->Draw("hist");
  EleBDTGV4_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV4" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V5
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV5_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV5_MC, LegendLabel2.c_str(), "L");

  EleBDTGV5_Data->SetLineColor(kBlue);
  EleBDTGV5_MC->SetLineColor(kRed);   

  EleBDTGV5_Data->SetMaximum(TMath::Max(EleBDTGV5_Data->GetMaximum(),EleBDTGV5_MC->GetMaximum()) * 1.2);
  EleBDTGV5_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV5_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV5_Data->Draw("hist");
  EleBDTGV5_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV5" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V6
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV6_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV6_MC, LegendLabel2.c_str(), "L");

  EleBDTGV6_Data->SetLineColor(kBlue);
  EleBDTGV6_MC->SetLineColor(kRed);   

  EleBDTGV6_Data->SetMaximum(TMath::Max(EleBDTGV6_Data->GetMaximum(),EleBDTGV6_MC->GetMaximum()) * 1.2);
  EleBDTGV6_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV6_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV6_Data->Draw("hist");
  EleBDTGV6_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV6" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V7
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV7_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV7_MC, LegendLabel2.c_str(), "L");

  EleBDTGV7_Data->SetLineColor(kBlue);
  EleBDTGV7_MC->SetLineColor(kRed);   

  EleBDTGV7_Data->SetMaximum(TMath::Max(EleBDTGV7_Data->GetMaximum(),EleBDTGV7_MC->GetMaximum()) * 1.2);
  EleBDTGV7_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV7_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV7_Data->Draw("hist");
  EleBDTGV7_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV7" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V8
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV8_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV8_MC, LegendLabel2.c_str(), "L");

  EleBDTGV8_Data->SetLineColor(kBlue);
  EleBDTGV8_MC->SetLineColor(kRed);   

  EleBDTGV8_Data->SetMaximum(TMath::Max(EleBDTGV8_Data->GetMaximum(),EleBDTGV8_MC->GetMaximum()) * 1.2);
  EleBDTGV8_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV8_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV8_Data->Draw("hist");
  EleBDTGV8_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV8" + label + ".gif").c_str());


  //*****************************************************************************************
  //BDTG V9
  //*****************************************************************************************
  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.20,0.70,0.40,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  legend->AddEntry(EleBDTGV9_Data, LegendLabel1.c_str(), "L");
  legend->AddEntry(EleBDTGV9_MC, LegendLabel2.c_str(), "L");

  EleBDTGV9_Data->SetLineColor(kBlue);
  EleBDTGV9_MC->SetLineColor(kRed);   

  EleBDTGV9_Data->SetMaximum(TMath::Max(EleBDTGV9_Data->GetMaximum(),EleBDTGV9_MC->GetMaximum()) * 1.2);
  EleBDTGV9_Data->GetYaxis()->SetTitleOffset(1.2);
  EleBDTGV9_Data->GetXaxis()->SetTitleOffset(1.05);
  EleBDTGV9_Data->Draw("hist");
  EleBDTGV9_MC->Draw("hist,same"); 
  legend->Draw();

  cv->SaveAs(("EleBDTGV9" + label + ".gif").c_str());





//   //*****************************************************************************************
//   //Save Histograms in file
//   //*****************************************************************************************
//   TFile *file = new TFile("ElectronIDMVAResults.root", "UPDATE");
//   file->cd();
//   file->WriteTObject(EleIDStandardLikelihood_Data, EleIDStandardLikelihood_Data->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDT_Data, EleIDBDT_Data->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDTG_Data, EleIDBDTG_Data->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDNN_Data, EleIDNN_Data->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihood_Data, EleIDTMVALikelihood_Data->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihoodD_Data, EleIDTMVALikelihoodD_Data->GetName(), "WriteDelete");  

//   file->WriteTObject(EleIDStandardLikelihood_MC, EleIDStandardLikelihood_MC->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDT_MC, EleIDBDT_MC->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDBDTG_MC, EleIDBDTG_MC->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDNN_MC, EleIDNN_MC->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihood_MC, EleIDTMVALikelihood_MC->GetName(), "WriteDelete");  
//   file->WriteTObject(EleIDTMVALikelihoodD_MC, EleIDTMVALikelihoodD_MC->GetName(), "WriteDelete"); 

//   file->WriteTObject(ROC_StandardLikelihood, ROC_StandardLikelihood->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_BDT, ROC_BDT->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_BDTG, ROC_BDTG->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_NN, ROC_NN->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihood, ROC_TMVALikelihood->GetName(), "WriteDelete");  
//   file->WriteTObject(ROC_TMVALikelihoodD, ROC_TMVALikelihoodD->GetName(), "WriteDelete");  
 
//   file->Close();
//   delete file;




  gBenchmark->Show("WWTemplate");       
} 

