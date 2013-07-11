//root -l EWKAna/Hww/LeptonSelection/MakeTagAndProbeIsolationDistributions.C+\(\"\"\)



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
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <TGraphAsymmErrors.h>      // 3D vector class
#include <TLegend.h>                // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum, Int_t triggerSelection);
Bool_t passElectronTagCuts(const mithep::TElectron *ele);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronProbeCuts(const mithep::TElectron *ele);
Bool_t passMuonTagCuts(const mithep::TMuon *mu);
Bool_t passMuonCuts(const mithep::TMuon *mu);
Bool_t passMuonProbeCuts(const mithep::TMuon *mu);
void WriteMassToFile( ofstream *file, Double_t mass);
void WriteMassToFile( ofstream *file, Double_t mass, Int_t evtNum);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2);
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

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}

//=== MAIN MACRO =================================================================================================

void MakeTagAndProbeIsolationDistributions(const string Label) {   

  gBenchmark->Start("ElectronTagAndProbe");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  string label = Label;
  if (Label != "") label = "_" + label;

  Double_t lumi;              // luminosity (pb^-1)
    

  //********************************************************
  // Define Bins
  //********************************************************
  vector<Double_t> Numerator_NVertex;
  vector<Double_t> Denominator_NVertex;

  for (UInt_t i=0 ; i < 20 ; ++i) {
    Numerator_NVertex.push_back(0);
    Denominator_NVertex.push_back(0);
  }

  vector<vector<string> > inputFiles;
  vector<string> processNames;
//         inputFiles.push_back(vector<string>());
//         inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full_TightPlusRecoTriggerSkim.root");
//         processNames.push_back("DataTagAndProbe");
   inputFiles.push_back(vector<string>());
 //   inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zeem20-v1g1-pu_noskim_normalized.root");
    inputFiles.back().push_back("/home/sixie/hist/HwwAnalysis/normalized/HwwAnalysis_p11-zmmm20-v1g1-pu_noskim_normalized.root");
   processNames.push_back("ZMCTagAndProbe");

  assert(processNames.size() == inputFiles.size());


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector<Double_t> ElectronHOverEEfficiencyBarrelNumerator;
  vector<Double_t> ElectronL1CorrectedHOverEEfficiencyBarrelNumerator;
  vector<Double_t> ElectronHOverEEfficiencyBarrelDenominator;
  vector<Double_t> ElectronHOverEEfficiencyEndcapNumerator;
  vector<Double_t> ElectronL1CorrectedHOverEEfficiencyEndcapNumerator;
  vector<Double_t> ElectronHOverEEfficiencyEndcapDenominator;

  vector<Double_t> ElectronIsolationEfficiencyBarrelNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyBarrelNumerator;
  vector<Double_t> ElectronIsolationEfficiencyBarrelDenominator;
  vector<Double_t> ElectronIsolationEfficiencyEndcapNumerator;
  vector<Double_t> ElectronL1CorrectedIsolationEfficiencyEndcapNumerator;
  vector<Double_t> ElectronIsolationEfficiencyEndcapDenominator;
  vector<Double_t> MuonIsolationEfficiencyNumerator;
  vector<Double_t> MuonL1CorrectedIsolationEfficiencyNumerator;
  vector<Double_t> MuonIsolationEfficiencyDenominator;
  for(UInt_t i=0 ; i<20; ++i) {
    ElectronHOverEEfficiencyBarrelNumerator.push_back(0.0);
    ElectronL1CorrectedHOverEEfficiencyBarrelNumerator.push_back(0.0);
    ElectronHOverEEfficiencyBarrelDenominator.push_back(0.0);
    ElectronHOverEEfficiencyEndcapNumerator.push_back(0.0);
    ElectronL1CorrectedHOverEEfficiencyEndcapNumerator.push_back(0.0);
    ElectronHOverEEfficiencyEndcapDenominator.push_back(0.0);
    
    ElectronIsolationEfficiencyBarrelNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyBarrelNumerator.push_back(0.0);
    ElectronIsolationEfficiencyBarrelDenominator.push_back(0.0);
    ElectronIsolationEfficiencyEndcapNumerator.push_back(0.0);
    ElectronL1CorrectedIsolationEfficiencyEndcapNumerator.push_back(0.0);
    ElectronIsolationEfficiencyEndcapDenominator.push_back(0.0);
    MuonIsolationEfficiencyNumerator.push_back(0.0);
    MuonL1CorrectedIsolationEfficiencyNumerator.push_back(0.0);
    MuonIsolationEfficiencyDenominator.push_back(0.0);
  }


  
  vector <TH1F*>  NVertices;
  vector<vector<TH1F*> >  Rho;
  vector<vector<TH1F*> >  Electron_caloIso_Barrel;
  vector<vector<TH1F*> >  Electron_caloIso_Endcap;
  vector<vector<TH1F*> >  Electron_ecalIso_Barrel;
  vector<vector<TH1F*> >  Electron_ecalIso_Endcap;
  vector<vector<TH1F*> >  Electron_hcalIso_Barrel;
  vector<vector<TH1F*> >  Electron_hcalIso_Endcap;
  vector<vector<TH1F*> >  Electron_trkIso_Barrel;
  vector<vector<TH1F*> >  Electron_trkIso_Endcap;
  vector<vector<TH1F*> >  Electron_ChargedIso_Barrel;
  vector<vector<TH1F*> >  Electron_ChargedIsoNoPU_Barrel;
  vector<vector<TH1F*> >  Electron_NeutralHadronIso_Barrel;
  vector<vector<TH1F*> >  Electron_GammaIso_Barrel;
  vector<vector<TH1F*> >  Electron_TotalPFIso_Barrel;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFIso_Barrel;
  vector<vector<TH1F*> >  Electron_ChargedIso_Endcap;
  vector<vector<TH1F*> >  Electron_ChargedIsoNoPU_Endcap;
  vector<vector<TH1F*> >  Electron_NeutralHadronIso_Endcap;
  vector<vector<TH1F*> >  Electron_GammaIso_Endcap;
  vector<vector<TH1F*> >  Electron_TotalPFIso_Endcap;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFIso_Endcap;
  vector<vector<TH1F*> >  Muon_caloIso_Barrel;
  vector<vector<TH1F*> >  Muon_caloIso_Endcap;
  vector<vector<TH1F*> >  Muon_ecalIso_Barrel;
  vector<vector<TH1F*> >  Muon_ecalIso_Endcap;
  vector<vector<TH1F*> >  Muon_hcalIso_Barrel;
  vector<vector<TH1F*> >  Muon_hcalIso_Endcap;
  vector<vector<TH1F*> >  Muon_trkIso_Barrel;
  vector<vector<TH1F*> >  Muon_trkIso_Endcap;
  vector<vector<TH1F*> >  Muon_ChargedIso_Barrel;
  vector<vector<TH1F*> >  Muon_ChargedIsoNoPU_Barrel;
  vector<vector<TH1F*> >  Muon_NeutralIso_Barrel;
  vector<vector<TH1F*> >  Muon_TotalPFIso_Barrel;
  vector<vector<TH1F*> >  Muon_ChargedIso_Endcap;
  vector<vector<TH1F*> >  Muon_ChargedIsoNoPU_Endcap;
  vector<vector<TH1F*> >  Muon_NeutralIso_Endcap;
  vector<vector<TH1F*> >  Muon_TotalPFIso_Endcap;

  vector<vector<TH1F*> >  Electron_caloIso04_Barrel;
  vector<vector<TH1F*> >  Electron_caloIso04_Endcap;
  vector<vector<TH1F*> >  Electron_ecalIso04_Barrel;
  vector<vector<TH1F*> >  Electron_ecalIso04_Endcap;
  vector<vector<TH1F*> >  Electron_hcalIso04_Barrel;
  vector<vector<TH1F*> >  Electron_hcalIso04_Endcap;
  vector<vector<TH1F*> >  Electron_trkIso04_Barrel;
  vector<vector<TH1F*> >  Electron_trkIso04_Endcap;
  vector<vector<TH1F*> >  Muon_caloIso05_Barrel;
  vector<vector<TH1F*> >  Muon_caloIso05_Endcap;
  vector<vector<TH1F*> >  Muon_ecalIso05_Barrel;
  vector<vector<TH1F*> >  Muon_ecalIso05_Endcap;
  vector<vector<TH1F*> >  Muon_hcalIso05_Barrel;
  vector<vector<TH1F*> >  Muon_hcalIso05_Endcap;
  vector<vector<TH1F*> >  Muon_trkIso05_Barrel;
  vector<vector<TH1F*> >  Muon_trkIso05_Endcap;

  vector<vector<TH1F*> >  Electron_RelIso_Barrel;
  vector<vector<TH1F*> >  Electron_TotalPFRelIso_Barrel;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFRelIso_Barrel;
  vector<vector<TH1F*> >  Electron_VertexSelectedPFRelIso_Barrel;
  vector<vector<TH1F*> >  Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
  vector<vector<TH1F*> >  Electron_RelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Electron_TotalPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Electron_RelIso_Endcap;
  vector<vector<TH1F*> >  Electron_TotalPFRelIso_Endcap;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFRelIso_Endcap;
  vector<vector<TH1F*> >  Electron_VertexSelectedPFRelIso_Endcap;
  vector<vector<TH1F*> >  Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
  vector<vector<TH1F*> >  Electron_RelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Electron_TotalPFRelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Muon_RelIso_Barrel;
  vector<vector<TH1F*> >  Muon_TotalPFRelIso_Barrel;
  vector<vector<TH1F*> >  Muon_VertexSelectedPFRelIso_Barrel;
  vector<vector<TH1F*> >  Muon_RelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Muon_TotalPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
  vector<vector<TH1F*> >  Muon_RelIso_Endcap;
  vector<vector<TH1F*> >  Muon_TotalPFRelIso_Endcap;
  vector<vector<TH1F*> >  Muon_VertexSelectedPFRelIso_Endcap;
  vector<vector<TH1F*> >  Muon_RelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Muon_TotalPFRelIsoRhoCorrected_Endcap;
  vector<vector<TH1F*> >  Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap;



 
//   vector <TH1F*>  Electron_relIso_Barrel;
//   vector <TH1F*>  Electron_relIso_Endcap;
//   vector <TH1F*>  Muon_relIso;
//   vector <TH1F*>  Electron_relIsoL1Corrected_Barrel;
//   vector <TH1F*>  Electron_relIsoL1Corrected_Endcap;
//   vector <TH1F*>  Muon_relIsoL1Corrected;
//   vector <TH1F*>  Electron_relIso04_Barrel;
//   vector <TH1F*>  Electron_relIso04_Endcap;
//   vector <TH1F*>  Muon_relIso05;
//   vector <TH1F*>  Electron_relIso04L1Corrected_Barrel;
//   vector <TH1F*>  Electron_relIso04L1Corrected_Endcap;
//   vector <TH1F*>  Muon_relIso05L1Corrected;

  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpNVertices = new TH1F((string("NVertices_")+processNames[q]+label).c_str(), "; Number of Vertices; Number of Events ", 20, -0.5, 19.5);
    
    vector<TH1F*>   tempVectorRho;
    vector<TH1F*>   tempVectorElectron_caloIso_Barrel;
    vector<TH1F*>   tempVectorElectron_caloIso_Endcap;
    vector<TH1F*>   tempVectorElectron_ecalIso_Barrel;
    vector<TH1F*>   tempVectorElectron_ecalIso_Endcap;
    vector<TH1F*>   tempVectorElectron_hcalIso_Barrel;
    vector<TH1F*>   tempVectorElectron_hcalIso_Endcap;
    vector<TH1F*>   tempVectorElectron_trkIso_Barrel;
    vector<TH1F*>   tempVectorElectron_trkIso_Endcap;
    vector<TH1F*>   tempVectorElectron_ChargedIso_Barrel;
    vector<TH1F*>   tempVectorElectron_ChargedIsoNoPU_Barrel;
    vector<TH1F*>   tempVectorElectron_NeutralHadronIso_Barrel;
    vector<TH1F*>   tempVectorElectron_GammaIso_Barrel;
    vector<TH1F*>   tempVectorElectron_TotalPFIso_Barrel;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFIso_Barrel;
    vector<TH1F*>   tempVectorElectron_ChargedIso_Endcap;
    vector<TH1F*>   tempVectorElectron_ChargedIsoNoPU_Endcap;
    vector<TH1F*>   tempVectorElectron_NeutralHadronIso_Endcap;
    vector<TH1F*>   tempVectorElectron_GammaIso_Endcap;
    vector<TH1F*>   tempVectorElectron_TotalPFIso_Endcap;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFIso_Endcap;
    vector<TH1F*>   tempVectorMuon_caloIso_Barrel;
    vector<TH1F*>   tempVectorMuon_caloIso_Endcap;
    vector<TH1F*>   tempVectorMuon_ecalIso_Barrel;
    vector<TH1F*>   tempVectorMuon_ecalIso_Endcap;
    vector<TH1F*>   tempVectorMuon_hcalIso_Barrel;
    vector<TH1F*>   tempVectorMuon_hcalIso_Endcap;
    vector<TH1F*>   tempVectorMuon_trkIso_Barrel;
    vector<TH1F*>   tempVectorMuon_trkIso_Endcap;
    vector<TH1F*>   tempVectorMuon_ChargedIso_Barrel;
    vector<TH1F*>   tempVectorMuon_ChargedIsoNoPU_Barrel;
    vector<TH1F*>   tempVectorMuon_NeutralIso_Barrel;
    vector<TH1F*>   tempVectorMuon_TotalPFIso_Barrel;
    vector<TH1F*>   tempVectorMuon_ChargedIso_Endcap;
    vector<TH1F*>   tempVectorMuon_ChargedIsoNoPU_Endcap;
    vector<TH1F*>   tempVectorMuon_NeutralIso_Endcap;
    vector<TH1F*>   tempVectorMuon_TotalPFIso_Endcap;


    vector<TH1F*>   tempVectorElectron_caloIso04_Barrel;
    vector<TH1F*>   tempVectorElectron_caloIso04_Endcap;
    vector<TH1F*>   tempVectorElectron_ecalIso04_Barrel;
    vector<TH1F*>   tempVectorElectron_ecalIso04_Endcap;
    vector<TH1F*>   tempVectorElectron_hcalIso04_Barrel;
    vector<TH1F*>   tempVectorElectron_hcalIso04_Endcap;
    vector<TH1F*>   tempVectorElectron_trkIso04_Barrel;
    vector<TH1F*>   tempVectorElectron_trkIso04_Endcap;
    vector<TH1F*>   tempVectorMuon_caloIso05_Barrel;
    vector<TH1F*>   tempVectorMuon_caloIso05_Endcap;
    vector<TH1F*>   tempVectorMuon_ecalIso05_Barrel;
    vector<TH1F*>   tempVectorMuon_ecalIso05_Endcap;
    vector<TH1F*>   tempVectorMuon_hcalIso05_Barrel;
    vector<TH1F*>   tempVectorMuon_hcalIso05_Endcap;
    vector<TH1F*>   tempVectorMuon_trkIso05_Barrel;
    vector<TH1F*>   tempVectorMuon_trkIso05_Endcap;

    vector<TH1F*>   tempVectorElectron_RelIso_Barrel;
    vector<TH1F*>   tempVectorElectron_TotalPFRelIso_Barrel;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFRelIso_Barrel;
    vector<TH1F*>   tempVectorElectron_VertexSelectedPFRelIso_Barrel;
    vector<TH1F*>   tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel;
    vector<TH1F*>   tempVectorElectron_RelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorElectron_TotalPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorElectron_RelIso_Endcap;
    vector<TH1F*>   tempVectorElectron_TotalPFRelIso_Endcap;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFRelIso_Endcap;
    vector<TH1F*>   tempVectorElectron_VertexSelectedPFRelIso_Endcap;
    vector<TH1F*>   tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap;
    vector<TH1F*>   tempVectorElectron_RelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorElectron_TotalPFRelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorMuon_RelIso_Barrel;
    vector<TH1F*>   tempVectorMuon_TotalPFRelIso_Barrel;
    vector<TH1F*>   tempVectorMuon_VertexSelectedPFRelIso_Barrel;
    vector<TH1F*>   tempVectorMuon_RelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorMuon_TotalPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel;
    vector<TH1F*>   tempVectorMuon_RelIso_Endcap;
    vector<TH1F*>   tempVectorMuon_TotalPFRelIso_Endcap;
    vector<TH1F*>   tempVectorMuon_VertexSelectedPFRelIso_Endcap;
    vector<TH1F*>   tempVectorMuon_RelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorMuon_TotalPFRelIsoRhoCorrected_Endcap;
    vector<TH1F*>   tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap;
 




    for (int n=0; n < 20; ++n) {     
      TH1F *tmpRho = new TH1F((string("Rho_") + IntToString(n) + "_" + processNames[q] + label).c_str(), "; Rho; Number of Events ", 2000, -20, 20);
      TH1F *tmpElectron_caloIso_Barrel = new TH1F((string("Electron_caloIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_caloIso_Endcap = new TH1F((string("Electron_caloIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ecalIso_Barrel = new TH1F((string("Electron_ecalIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ecalIso_Endcap = new TH1F((string("Electron_ecalIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_hcalIso_Barrel = new TH1F((string("Electron_hcalIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_hcalIso_Endcap = new TH1F((string("Electron_hcalIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_trkIso_Barrel = new TH1F((string("Electron_trkIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_trkIso_Endcap = new TH1F((string("Electron_trkIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ChargedIso_Barrel = new TH1F((string("Electron_ChargedIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ChargedIsoNoPU_Barrel = new TH1F((string("Electron_ChargedIsoNoPU_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ChargedIsoNoPU; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_NeutralHadronIso_Barrel = new TH1F((string("Electron_NeutralHadronIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_GammaIso_Barrel = new TH1F((string("Electron_GammaIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_TotalPFIso_Barrel = new TH1F((string("Electron_TotalPFIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_FootprintRemovedPFIso_Barrel = new TH1F((string("Electron_FootprintRemovedPFIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; FootprintRemovedPFIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ChargedIso_Endcap = new TH1F((string("Electron_ChargedIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ChargedIsoNoPU_Endcap = new TH1F((string("Electron_ChargedIsoNoPU_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ChargedIsoNoPU; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_NeutralHadronIso_Endcap = new TH1F((string("Electron_NeutralHadronIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; NeutralHadronIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_GammaIso_Endcap = new TH1F((string("Electron_GammaIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; GammaIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_TotalPFIso_Endcap = new TH1F((string("Electron_TotalPFIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_FootprintRemovedPFIso_Endcap = new TH1F((string("Electron_FootprintRemovedPFIso_Endcap_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; FootprintRemovedPFIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_caloIso_Barrel = new TH1F((string("Muon_caloIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_caloIso_Endcap = new TH1F((string("Muon_caloIso_Endcap_")+ IntToString(n) + "_" +processNames[q]+label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ecalIso_Barrel = new TH1F((string("Muon_ecalIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ecalIso_Endcap = new TH1F((string("Muon_ecalIso_Endcap_")+ IntToString(n) + "_" +processNames[q]+label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_hcalIso_Barrel = new TH1F((string("Muon_hcalIso_Barrel_")+ IntToString(n) + "_" +processNames[q] +label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_hcalIso_Endcap = new TH1F((string("Muon_hcalIso_Endcap_")+ IntToString(n) + "_" +processNames[q]+label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_trkIso_Barrel = new TH1F((string("Muon_trkIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_trkIso_Endcap = new TH1F((string("Muon_trkIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ChargedIso_Barrel = new TH1F((string("Muon_ChargedIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ChargedIsoNoPU_Barrel = new TH1F((string("Muon_ChargedIsoNoPU_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ChargedIsoNoPU; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_NeutralIso_Barrel = new TH1F((string("Muon_NeutralIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; NeutralIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_TotalPFIso_Barrel = new TH1F((string("Muon_TotalPFIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ChargedIso_Endcap = new TH1F((string("Muon_ChargedIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ChargedIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ChargedIsoNoPU_Endcap = new TH1F((string("Muon_ChargedIsoNoPU_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ChargedIsoNoPU; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_NeutralIso_Endcap = new TH1F((string("Muon_NeutralIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; NeutralIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_TotalPFIso_Endcap = new TH1F((string("Muon_TotalPFIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFIso; Fraction of Events ", 2000, -200.0, 200.0);

      TH1F *tmpElectron_caloIso04_Barrel = new TH1F((string("Electron_caloIso04_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_caloIso04_Endcap = new TH1F((string("Electron_caloIso04_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; caloIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ecalIso04_Barrel = new TH1F((string("Electron_ecalIso04_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_ecalIso04_Endcap = new TH1F((string("Electron_ecalIso04_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ecalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_hcalIso04_Barrel = new TH1F((string("Electron_hcalIso04_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_hcalIso04_Endcap = new TH1F((string("Electron_hcalIso04_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; hcalIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_trkIso04_Barrel = new TH1F((string("Electron_trkIso04_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpElectron_trkIso04_Endcap = new TH1F((string("Electron_trkIso04_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_caloIso05_Barrel = new TH1F((string("Muon_caloIso05_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; caloIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_caloIso05_Endcap = new TH1F((string("Muon_caloIso05_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; caloIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ecalIso05_Barrel = new TH1F((string("Muon_ecalIso05_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ecalIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_ecalIso05_Endcap = new TH1F((string("Muon_ecalIso05_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; ecalIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_hcalIso05_Barrel = new TH1F((string("Muon_hcalIso05_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; hcalIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_hcalIso05_Endcap = new TH1F((string("Muon_hcalIso05_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; hcalIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_trkIso05_Barrel = new TH1F((string("Muon_trkIso05_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso05; Fraction of Events ", 2000, -200.0, 200.0);
      TH1F *tmpMuon_trkIso05_Endcap = new TH1F((string("Muon_trkIso05_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; trkIso05; Fraction of Events ", 2000, -200.0, 200.0);


      TH1F *tmpElectron_RelIso_Barrel = new TH1F((string("Electron_RelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_TotalPFRelIso_Barrel = new TH1F((string("Electron_TotalPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_FootprintRemovedPFRelIso_Barrel = new TH1F((string("Electron_FootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; FootprintRemovedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedPFRelIso_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel = new TH1F((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedFootprintRemovedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_RelIsoRhoCorrected_Barrel = new TH1F((string("Electron_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_TotalPFRelIsoRhoCorrected_Barrel = new TH1F((string("Electron_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel = new TH1F((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; FootprintRemovedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel = new TH1F((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedFootprintRemovedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_RelIso_Endcap = new TH1F((string("Electron_RelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_TotalPFRelIso_Endcap = new TH1F((string("Electron_TotalPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_FootprintRemovedPFRelIso_Endcap = new TH1F((string("Electron_FootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; FootprintRemovedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedPFRelIso_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap = new TH1F((string("Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedFootprintRemovedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_RelIsoRhoCorrected_Endcap = new TH1F((string("Electron_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_TotalPFRelIsoRhoCorrected_Endcap = new TH1F((string("Electron_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap = new TH1F((string("Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; FootprintRemovedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap = new TH1F((string("Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedFootprintRemovedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_RelIso_Barrel = new TH1F((string("Muon_RelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_TotalPFRelIso_Barrel = new TH1F((string("Muon_TotalPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_VertexSelectedPFRelIso_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIso_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_RelIsoRhoCorrected_Barrel = new TH1F((string("Muon_RelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_TotalPFRelIsoRhoCorrected_Barrel = new TH1F((string("Muon_TotalPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel = new TH1F((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_RelIso_Endcap = new TH1F((string("Muon_RelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_TotalPFRelIso_Endcap = new TH1F((string("Muon_TotalPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_VertexSelectedPFRelIso_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIso_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIso; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_RelIsoRhoCorrected_Endcap = new TH1F((string("Muon_RelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; RelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_TotalPFRelIsoRhoCorrected_Endcap = new TH1F((string("Muon_TotalPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; TotalPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);
      TH1F *tmpMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap = new TH1F((string("Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap_") + IntToString(n) + "_" +processNames[q]+label).c_str(), "; VertexSelectedPFRelIsoRhoCorrected; Fraction of Events ", 2000, -4.0, 4.0);


      tempVectorRho.push_back(tmpRho);
      tempVectorElectron_caloIso_Barrel.push_back(tmpElectron_caloIso_Barrel);
      tempVectorElectron_caloIso_Endcap.push_back(tmpElectron_caloIso_Endcap);
      tempVectorElectron_ecalIso_Barrel.push_back(tmpElectron_ecalIso_Barrel);
      tempVectorElectron_ecalIso_Endcap.push_back(tmpElectron_ecalIso_Endcap);
      tempVectorElectron_hcalIso_Barrel.push_back(tmpElectron_hcalIso_Barrel);
      tempVectorElectron_hcalIso_Endcap.push_back(tmpElectron_hcalIso_Endcap);
      tempVectorElectron_trkIso_Barrel.push_back(tmpElectron_trkIso_Barrel);
      tempVectorElectron_trkIso_Endcap.push_back(tmpElectron_trkIso_Endcap);
      tempVectorElectron_ChargedIso_Barrel.push_back(tmpElectron_ChargedIso_Barrel);
      tempVectorElectron_ChargedIsoNoPU_Barrel.push_back(tmpElectron_ChargedIsoNoPU_Barrel);
      tempVectorElectron_NeutralHadronIso_Barrel.push_back(tmpElectron_NeutralHadronIso_Barrel);
      tempVectorElectron_GammaIso_Barrel.push_back(tmpElectron_GammaIso_Barrel);
      tempVectorElectron_TotalPFIso_Barrel.push_back(tmpElectron_TotalPFIso_Barrel);
      tempVectorElectron_FootprintRemovedPFIso_Barrel.push_back(tmpElectron_FootprintRemovedPFIso_Barrel);
      tempVectorElectron_ChargedIso_Endcap.push_back(tmpElectron_ChargedIso_Endcap);
      tempVectorElectron_ChargedIsoNoPU_Endcap.push_back(tmpElectron_ChargedIsoNoPU_Endcap);
      tempVectorElectron_NeutralHadronIso_Endcap.push_back(tmpElectron_NeutralHadronIso_Endcap);
      tempVectorElectron_GammaIso_Endcap.push_back(tmpElectron_GammaIso_Endcap);
      tempVectorElectron_TotalPFIso_Endcap.push_back(tmpElectron_TotalPFIso_Endcap);
      tempVectorElectron_FootprintRemovedPFIso_Endcap.push_back(tmpElectron_FootprintRemovedPFIso_Endcap);
      tempVectorMuon_caloIso_Barrel.push_back(tmpMuon_caloIso_Barrel);
      tempVectorMuon_caloIso_Endcap.push_back(tmpMuon_caloIso_Endcap);
      tempVectorMuon_ecalIso_Barrel.push_back(tmpMuon_ecalIso_Barrel);
      tempVectorMuon_ecalIso_Endcap.push_back(tmpMuon_ecalIso_Endcap);
      tempVectorMuon_hcalIso_Barrel.push_back(tmpMuon_hcalIso_Barrel);
      tempVectorMuon_hcalIso_Endcap.push_back(tmpMuon_hcalIso_Endcap);
      tempVectorMuon_trkIso_Barrel.push_back(tmpMuon_trkIso_Barrel);
      tempVectorMuon_trkIso_Endcap.push_back(tmpMuon_trkIso_Endcap);
      tempVectorMuon_ChargedIso_Barrel.push_back(tmpMuon_ChargedIso_Barrel);
      tempVectorMuon_ChargedIsoNoPU_Barrel.push_back(tmpMuon_ChargedIsoNoPU_Barrel);
      tempVectorMuon_NeutralIso_Barrel.push_back(tmpMuon_NeutralIso_Barrel);
      tempVectorMuon_TotalPFIso_Barrel.push_back(tmpMuon_TotalPFIso_Barrel);
      tempVectorMuon_ChargedIso_Endcap.push_back(tmpMuon_ChargedIso_Endcap);
      tempVectorMuon_ChargedIsoNoPU_Endcap.push_back(tmpMuon_ChargedIsoNoPU_Endcap);
      tempVectorMuon_NeutralIso_Endcap.push_back(tmpMuon_NeutralIso_Endcap);
      tempVectorMuon_TotalPFIso_Endcap.push_back(tmpMuon_TotalPFIso_Endcap);
      tempVectorElectron_caloIso04_Barrel.push_back(tmpElectron_caloIso04_Barrel);
      tempVectorElectron_caloIso04_Endcap.push_back(tmpElectron_caloIso04_Endcap);
      tempVectorElectron_ecalIso04_Barrel.push_back(tmpElectron_ecalIso04_Barrel);
      tempVectorElectron_ecalIso04_Endcap.push_back(tmpElectron_ecalIso04_Endcap);
      tempVectorElectron_hcalIso04_Barrel.push_back(tmpElectron_hcalIso04_Barrel);
      tempVectorElectron_hcalIso04_Endcap.push_back(tmpElectron_hcalIso04_Endcap);
      tempVectorElectron_trkIso04_Barrel.push_back(tmpElectron_trkIso04_Barrel);
      tempVectorElectron_trkIso04_Endcap.push_back(tmpElectron_trkIso04_Endcap);
      tempVectorMuon_caloIso05_Barrel.push_back(tmpMuon_caloIso05_Barrel);
      tempVectorMuon_caloIso05_Endcap.push_back(tmpMuon_caloIso05_Endcap);
      tempVectorMuon_ecalIso05_Barrel.push_back(tmpMuon_ecalIso05_Barrel);
      tempVectorMuon_ecalIso05_Endcap.push_back(tmpMuon_ecalIso05_Endcap);
      tempVectorMuon_hcalIso05_Barrel.push_back(tmpMuon_hcalIso05_Barrel);
      tempVectorMuon_hcalIso05_Endcap.push_back(tmpMuon_hcalIso05_Endcap);
      tempVectorMuon_trkIso05_Barrel.push_back(tmpMuon_trkIso05_Barrel);
      tempVectorMuon_trkIso05_Endcap.push_back(tmpMuon_trkIso05_Endcap);


    tempVectorElectron_RelIso_Barrel.push_back(tmpElectron_RelIso_Barrel);
    tempVectorElectron_TotalPFRelIso_Barrel.push_back(tmpElectron_TotalPFRelIso_Barrel);
    tempVectorElectron_FootprintRemovedPFRelIso_Barrel.push_back(tmpElectron_FootprintRemovedPFRelIso_Barrel);
    tempVectorElectron_VertexSelectedPFRelIso_Barrel.push_back(tmpElectron_VertexSelectedPFRelIso_Barrel);
    tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel.push_back(tmpElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel);
    tempVectorElectron_RelIsoRhoCorrected_Barrel.push_back(tmpElectron_RelIsoRhoCorrected_Barrel);
    tempVectorElectron_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpElectron_TotalPFRelIsoRhoCorrected_Barrel);
    tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel);
    tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tmpElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel);
    tempVectorElectron_RelIso_Endcap.push_back(tmpElectron_RelIso_Endcap);
    tempVectorElectron_TotalPFRelIso_Endcap.push_back(tmpElectron_TotalPFRelIso_Endcap);
    tempVectorElectron_FootprintRemovedPFRelIso_Endcap.push_back(tmpElectron_FootprintRemovedPFRelIso_Endcap);
    tempVectorElectron_VertexSelectedPFRelIso_Endcap.push_back(tmpElectron_VertexSelectedPFRelIso_Endcap);
    tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap.push_back(tmpElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap);
    tempVectorElectron_RelIsoRhoCorrected_Endcap.push_back(tmpElectron_RelIsoRhoCorrected_Endcap);
    tempVectorElectron_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpElectron_TotalPFRelIsoRhoCorrected_Endcap);
    tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap);
    tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tmpElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap);
    tempVectorMuon_RelIso_Barrel.push_back(tmpMuon_RelIso_Barrel);
    tempVectorMuon_TotalPFRelIso_Barrel.push_back(tmpMuon_TotalPFRelIso_Barrel);
    tempVectorMuon_VertexSelectedPFRelIso_Barrel.push_back(tmpMuon_VertexSelectedPFRelIso_Barrel);
    tempVectorMuon_RelIsoRhoCorrected_Barrel.push_back(tmpMuon_RelIsoRhoCorrected_Barrel);
    tempVectorMuon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tmpMuon_TotalPFRelIsoRhoCorrected_Barrel);
    tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tmpMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    tempVectorMuon_RelIso_Endcap.push_back(tmpMuon_RelIso_Endcap);
    tempVectorMuon_TotalPFRelIso_Endcap.push_back(tmpMuon_TotalPFRelIso_Endcap);
    tempVectorMuon_VertexSelectedPFRelIso_Endcap.push_back(tmpMuon_VertexSelectedPFRelIso_Endcap);
    tempVectorMuon_RelIsoRhoCorrected_Endcap.push_back(tmpMuon_RelIsoRhoCorrected_Endcap);
    tempVectorMuon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tmpMuon_TotalPFRelIsoRhoCorrected_Endcap);
    tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tmpMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);

    }

    NVertices.push_back(tmpNVertices);
    Rho.push_back(tempVectorRho);
    Electron_caloIso_Barrel.push_back(tempVectorElectron_caloIso_Barrel);
    Electron_caloIso_Endcap.push_back(tempVectorElectron_caloIso_Endcap);
    Electron_ecalIso_Barrel.push_back(tempVectorElectron_ecalIso_Barrel);
    Electron_ecalIso_Endcap.push_back(tempVectorElectron_ecalIso_Endcap);
    Electron_hcalIso_Barrel.push_back(tempVectorElectron_hcalIso_Barrel);
    Electron_hcalIso_Endcap.push_back(tempVectorElectron_hcalIso_Endcap);
    Electron_trkIso_Barrel.push_back(tempVectorElectron_trkIso_Barrel);
    Electron_trkIso_Endcap.push_back(tempVectorElectron_trkIso_Endcap);
    Electron_ChargedIso_Barrel.push_back(tempVectorElectron_ChargedIso_Barrel);
    Electron_ChargedIsoNoPU_Barrel.push_back(tempVectorElectron_ChargedIsoNoPU_Barrel);
    Electron_NeutralHadronIso_Barrel.push_back(tempVectorElectron_NeutralHadronIso_Barrel);
    Electron_GammaIso_Barrel.push_back(tempVectorElectron_GammaIso_Barrel);
    Electron_TotalPFIso_Barrel.push_back(tempVectorElectron_TotalPFIso_Barrel);
    Electron_FootprintRemovedPFIso_Barrel.push_back(tempVectorElectron_FootprintRemovedPFIso_Barrel);
    Electron_ChargedIso_Endcap.push_back(tempVectorElectron_ChargedIso_Endcap);
    Electron_ChargedIsoNoPU_Endcap.push_back(tempVectorElectron_ChargedIsoNoPU_Endcap);
    Electron_NeutralHadronIso_Endcap.push_back(tempVectorElectron_NeutralHadronIso_Endcap);
    Electron_GammaIso_Endcap.push_back(tempVectorElectron_GammaIso_Endcap);
    Electron_TotalPFIso_Endcap.push_back(tempVectorElectron_TotalPFIso_Endcap);
    Electron_FootprintRemovedPFIso_Endcap.push_back(tempVectorElectron_FootprintRemovedPFIso_Endcap);
    Muon_caloIso_Barrel.push_back(tempVectorMuon_caloIso_Barrel);
    Muon_caloIso_Endcap.push_back(tempVectorMuon_caloIso_Endcap);
    Muon_ecalIso_Barrel.push_back(tempVectorMuon_ecalIso_Barrel);
    Muon_ecalIso_Endcap.push_back(tempVectorMuon_ecalIso_Endcap);
    Muon_hcalIso_Barrel.push_back(tempVectorMuon_hcalIso_Barrel);
    Muon_hcalIso_Endcap.push_back(tempVectorMuon_hcalIso_Endcap);
    Muon_trkIso_Barrel.push_back(tempVectorMuon_trkIso_Barrel);
    Muon_trkIso_Endcap.push_back(tempVectorMuon_trkIso_Endcap);
    Muon_ChargedIso_Barrel.push_back(tempVectorMuon_ChargedIso_Barrel);
    Muon_ChargedIsoNoPU_Barrel.push_back(tempVectorMuon_ChargedIsoNoPU_Barrel);
    Muon_NeutralIso_Barrel.push_back(tempVectorMuon_NeutralIso_Barrel);
    Muon_TotalPFIso_Barrel.push_back(tempVectorMuon_TotalPFIso_Barrel);
    Muon_ChargedIso_Endcap.push_back(tempVectorMuon_ChargedIso_Endcap);
    Muon_ChargedIsoNoPU_Endcap.push_back(tempVectorMuon_ChargedIsoNoPU_Endcap);
    Muon_NeutralIso_Endcap.push_back(tempVectorMuon_NeutralIso_Endcap);
    Muon_TotalPFIso_Endcap.push_back(tempVectorMuon_TotalPFIso_Endcap);
    Electron_caloIso04_Barrel.push_back(tempVectorElectron_caloIso04_Barrel);
    Electron_caloIso04_Endcap.push_back(tempVectorElectron_caloIso04_Endcap);
    Electron_ecalIso04_Barrel.push_back(tempVectorElectron_ecalIso04_Barrel);
    Electron_ecalIso04_Endcap.push_back(tempVectorElectron_ecalIso04_Endcap);
    Electron_hcalIso04_Barrel.push_back(tempVectorElectron_hcalIso04_Barrel);
    Electron_hcalIso04_Endcap.push_back(tempVectorElectron_hcalIso04_Endcap);
    Electron_trkIso04_Barrel.push_back(tempVectorElectron_trkIso04_Barrel);
    Electron_trkIso04_Endcap.push_back(tempVectorElectron_trkIso04_Endcap);
    Muon_caloIso05_Barrel.push_back(tempVectorMuon_caloIso05_Barrel);
    Muon_caloIso05_Endcap.push_back(tempVectorMuon_caloIso05_Endcap);
    Muon_ecalIso05_Barrel.push_back(tempVectorMuon_ecalIso05_Barrel);
    Muon_ecalIso05_Endcap.push_back(tempVectorMuon_ecalIso05_Endcap);
    Muon_hcalIso05_Barrel.push_back(tempVectorMuon_hcalIso05_Barrel);
    Muon_hcalIso05_Endcap.push_back(tempVectorMuon_hcalIso05_Endcap);
    Muon_trkIso05_Barrel.push_back(tempVectorMuon_trkIso05_Barrel);
    Muon_trkIso05_Endcap.push_back(tempVectorMuon_trkIso05_Endcap);


    Electron_RelIso_Barrel.push_back(tempVectorElectron_RelIso_Barrel);
    Electron_TotalPFRelIso_Barrel.push_back(tempVectorElectron_TotalPFRelIso_Barrel);
    Electron_FootprintRemovedPFRelIso_Barrel.push_back(tempVectorElectron_FootprintRemovedPFRelIso_Barrel);
    Electron_VertexSelectedPFRelIso_Barrel.push_back(tempVectorElectron_VertexSelectedPFRelIso_Barrel);
    Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel.push_back(tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Barrel);
    Electron_RelIsoRhoCorrected_Barrel.push_back(tempVectorElectron_RelIsoRhoCorrected_Barrel);
    Electron_TotalPFRelIsoRhoCorrected_Barrel.push_back(tempVectorElectron_TotalPFRelIsoRhoCorrected_Barrel);
    Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Barrel);
    Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel.push_back(tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel);
    Electron_RelIso_Endcap.push_back(tempVectorElectron_RelIso_Endcap);
    Electron_TotalPFRelIso_Endcap.push_back(tempVectorElectron_TotalPFRelIso_Endcap);
    Electron_FootprintRemovedPFRelIso_Endcap.push_back(tempVectorElectron_FootprintRemovedPFRelIso_Endcap);
    Electron_VertexSelectedPFRelIso_Endcap.push_back(tempVectorElectron_VertexSelectedPFRelIso_Endcap);
    Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap.push_back(tempVectorElectron_VertexSelectedFootprintRemovedPFRelIso_Endcap);
    Electron_RelIsoRhoCorrected_Endcap.push_back(tempVectorElectron_RelIsoRhoCorrected_Endcap);
    Electron_TotalPFRelIsoRhoCorrected_Endcap.push_back(tempVectorElectron_TotalPFRelIsoRhoCorrected_Endcap);
    Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tempVectorElectron_FootprintRemovedPFRelIsoRhoCorrected_Endcap);
    Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tempVectorElectron_VertexSelectedPFRelIsoRhoCorrected_Endcap);
    Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap.push_back(tempVectorElectron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap);
    Muon_RelIso_Barrel.push_back(tempVectorMuon_RelIso_Barrel);
    Muon_TotalPFRelIso_Barrel.push_back(tempVectorMuon_TotalPFRelIso_Barrel);
    Muon_VertexSelectedPFRelIso_Barrel.push_back(tempVectorMuon_VertexSelectedPFRelIso_Barrel);
    Muon_RelIsoRhoCorrected_Barrel.push_back(tempVectorMuon_RelIsoRhoCorrected_Barrel);
    Muon_TotalPFRelIsoRhoCorrected_Barrel.push_back(tempVectorMuon_TotalPFRelIsoRhoCorrected_Barrel);
    Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel.push_back(tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Barrel);
    Muon_RelIso_Endcap.push_back(tempVectorMuon_RelIso_Endcap);
    Muon_TotalPFRelIso_Endcap.push_back(tempVectorMuon_TotalPFRelIso_Endcap);
    Muon_VertexSelectedPFRelIso_Endcap.push_back(tempVectorMuon_VertexSelectedPFRelIso_Endcap);
    Muon_RelIsoRhoCorrected_Endcap.push_back(tempVectorMuon_RelIsoRhoCorrected_Endcap);
    Muon_TotalPFRelIsoRhoCorrected_Endcap.push_back(tempVectorMuon_TotalPFRelIsoRhoCorrected_Endcap);
    Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap.push_back(tempVectorMuon_VertexSelectedPFRelIsoRhoCorrected_Endcap);


 //    TH1F *tmpElectron_relIsoL1Corrected_Barrel = new TH1F((string("Electron_relIsoL1Corrected_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpElectron_relIsoL1Corrected_Endcap = new TH1F((string("Electron_relIsoL1Corrected_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpMuon_relIsoL1Corrected = new TH1F((string("Muon_relIsoL1Corrected_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpElectron_relIso04_Barrel = new TH1F((string("Electron_relIso04_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpElectron_relIso04_Endcap = new TH1F((string("Electron_relIso04_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpMuon_relIso05 = new TH1F((string("Muon_relIso05_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpElectron_relIso04L1Corrected_Barrel = new TH1F((string("Electron_relIso04L1Corrected_Barrel_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpElectron_relIso04L1Corrected_Endcap = new TH1F((string("Electron_relIso04L1Corrected_Endcap_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     TH1F *tmpMuon_relIso05L1Corrected = new TH1F((string("Muon_relIso05L1Corrected_")+processNames[q]+label).c_str(), "; relIso; Fraction of Events ", 200, 0, 2.0);
//     Rho_NVertex2.push_back(tmpRho_NVertex2);
//     Rho_NVertex8.push_back(tmpRho_NVertex8);

//     Electron_relIso_Barrel.push_back(tmpElectron_relIso_Barrel);
//     Electron_relIso_Endcap.push_back(tmpElectron_relIso_Endcap);
//     Muon_relIso.push_back(tmpMuon_relIso);
//     Electron_relIsoL1Corrected_Barrel.push_back(tmpElectron_relIsoL1Corrected_Barrel);
//     Electron_relIsoL1Corrected_Endcap.push_back(tmpElectron_relIsoL1Corrected_Endcap);
//     Muon_relIsoL1Corrected.push_back(tmpMuon_relIsoL1Corrected);
//     Electron_relIso04_Barrel.push_back(tmpElectron_relIso04_Barrel);
//     Electron_relIso04_Endcap.push_back(tmpElectron_relIso04_Endcap);
//     Muon_relIso05.push_back(tmpMuon_relIso05);
//     Electron_relIso04L1Corrected_Barrel.push_back(tmpElectron_relIso04L1Corrected_Barrel);
//     Electron_relIso04L1Corrected_Endcap.push_back(tmpElectron_relIso04L1Corrected_Endcap);
//     Muon_relIso05L1Corrected.push_back(tmpMuon_relIso05L1Corrected);

  }


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
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  
    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("Cert_161079-161352_7TeV_PromptReco_Collisions11_JSON_noESpbl_v2.txt"); 
  hasJSON = kFALSE;

  for (int q = 0; q<inputFiles.size() ; ++q) { 
    for (int f = 0; f < inputFiles[q].size() ; ++f) {


      //********************************************************
      // Get Tree
      //********************************************************
      eventTree = getTreeFromFile(inputFiles[q][f].c_str(),"Events"); 

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");


      //*****************************************************************************************
      //Loop over Data Tree
      //*****************************************************************************************
      Double_t nsel=0, nselvar=0;
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);

        if (ientry % 100000 == 0 ) cout << "Event " << ientry << endl;

 //        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//         if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...    

        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear();    
        electronBr->GetEntry(ientry);
        muonArr->Clear();    
        muonBr->GetEntry(ientry);

        //********************************************************
        // TcMet
        //********************************************************
        TVector3 met;        
        if(info->tcMEx!=0 || info->tcMEy!=0) {       
          met.SetXYZ(info->tcMEx, info->tcMEy, 0);
        }

        //******************************************************************************
        //dilepton preselection
        //******************************************************************************
        if (!(electronArr->GetEntries() < 2 || muonArr->GetEntries() < 2)) continue;

        Double_t PUIsolationEnergy = info->PileupEnergyDensity * 3.14159 * pow(0.3,2) * 1.0 ;
        Double_t PUIsolationEnergy04Cone = info->PileupEnergyDensity * 3.14159 * pow(0.4,2) * 1.0 ;
        Double_t PUIsolationEnergy05Cone = info->PileupEnergyDensity * 3.14159 * pow(0.5,2) * 1.0 ;
        Double_t PUHOverE = info->PileupEnergyDensity * 3.14159 * pow(0.15,2) * 1.0 ;
        Int_t NVertex = info->nPV0; if (NVertex > 19) NVertex = 19;

       //******************************************************************************
        //loop over electron pairs
        //******************************************************************************
        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
          if ( !(
                 fabs(tag->eta) < 2.5
                 && tag->pt > 20.0
                 && passElectronTagCuts(tag)
                 )
            ) continue;
      

          for(Int_t j=0; j<electronArr->GetEntries(); j++) {
            if (i==j) continue;

            const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
            if ( !(
                   fabs(probe->eta) < 2.5
                   && probe->pt > 10.0
                   && passElectronProbeCuts(probe)
                   )
              ) continue;
        
 
            mithep::FourVectorM tagVector;
            mithep::FourVectorM probeVector;
            tagVector.SetCoordinates(tag->pt, tag->eta, tag->phi, 0.51099892e-3 );
            probeVector.SetCoordinates(probe->pt, probe->eta, probe->phi, 0.51099892e-3 );
            mithep::FourVectorM dilepton = tagVector+probeVector;
 
        
            if (dilepton.M() > 75 && dilepton.M() < 105
              ) {
              NVertices[q]->Fill(NVertex);
              Rho[q][NVertex]->Fill(info->PileupEnergyDensity);              


              Double_t relIso = (fabs(probe->eta) < 1.5)? (probe->trkIso03 + TMath::Max(probe->emIso03 - 1.0, 0.0) + probe->hadIso03)/probe->pt : (probe->trkIso03 + probe->emIso03  + probe->hadIso03)/probe->pt;
              Double_t relIsoL1Corrected = (fabs(probe->eta) < 1.5)? (probe->trkIso03 + TMath::Max(TMath::Max(probe->emIso03 - 1.0, 0.0) + probe->hadIso03 - PUIsolationEnergy,0.0)) / probe->pt:(probe->trkIso03 + TMath::Max(probe->emIso03 + probe->hadIso03 - PUIsolationEnergy,0.0)) / probe->pt;
              Double_t relIso04 = (fabs(probe->eta) < 1.5)? (probe->trkIso04 + TMath::Max(probe->emIso04 - 1.0, 0.0) + probe->hadIso04)/probe->pt : (probe->trkIso04 + probe->emIso04  + probe->hadIso04)/probe->pt;
              Double_t relIso04L1Corrected = (fabs(probe->eta) < 1.5)? (probe->trkIso04 + TMath::Max(TMath::Max(probe->emIso04 - 1.0, 0.0) + probe->hadIso04 - PUIsolationEnergy,0.0)) / probe->pt:(probe->trkIso04 + TMath::Max(probe->emIso04 + probe->hadIso04 - PUIsolationEnergy04Cone,0.0)) / probe->pt;

              if (fabs(probe->eta) < 1.5) {
                ElectronIsolationEfficiencyBarrelDenominator[NVertex]++;
                Electron_caloIso_Barrel[q][NVertex]->Fill(TMath::Min(double( TMath::Max(probe->emIso03 - 1.0, 0.0) + probe->hadIso03), double(199.99)));
                Electron_ecalIso_Barrel[q][NVertex]->Fill(TMath::Min(double( TMath::Max(probe->emIso03 - 1.0, 0.0)), double(199.99)));
                Electron_hcalIso_Barrel[q][NVertex]->Fill(TMath::Min(double( probe->hadIso03 ), double(199.99)));
                Electron_caloIso04_Barrel[q][NVertex]->Fill(TMath::Min( double(TMath::Max(probe->emIso04 - 1.0, 0.0) + probe->hadIso04), double(199.99)));
                Electron_ecalIso04_Barrel[q][NVertex]->Fill(TMath::Min( double(TMath::Max(probe->emIso04 - 1.0, 0.0) ), double(199.99)));
                Electron_hcalIso04_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->hadIso04), double(199.99)));
                Electron_trkIso_Barrel[q][NVertex]->Fill(TMath::Min(double( probe->trkIso03 ), double(199.99)));
                Electron_trkIso04_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->trkIso04 ) , double(199.99)));
                Electron_ChargedIso_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices ) , double(199.99)));
                Electron_ChargedIsoNoPU_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03FromOtherVertices ) , double(199.99)));
                Electron_NeutralHadronIso_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->NeutralHadronIso03) , double(199.99)));
                Electron_GammaIso_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->GammaIso03) , double(199.99)));
                Electron_TotalPFIso_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03) , double(199.99)));
                Electron_FootprintRemovedPFIso_Barrel[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip) , double(199.99)));
                if (relIso < 0.1) {
                  ElectronIsolationEfficiencyBarrelNumerator[NVertex]++;              
                }
                if (relIsoL1Corrected < 0.1) {
                  ElectronL1CorrectedIsolationEfficiencyBarrelNumerator[NVertex]++;              
                }  
                Electron_RelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((TMath::Max(probe->emIso03 - 1.0, 0.0) + probe->hadIso03 + probe->trkIso03)/probe->pt),double(3.99)));
                Electron_TotalPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03)/probe->pt),double(3.99)));
                Electron_FootprintRemovedPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip)/probe->pt),double(3.99)));
                Electron_VertexSelectedPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03+probe->GammaIso03)/probe->pt),double(3.99)));
                Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip)/probe->pt),double(3.99)));
                Electron_RelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((TMath::Max(probe->emIso03 - 1.0, 0.0) + probe->hadIso03 + probe->trkIso03 - info->PileupEnergyDensity*0.0798  )/probe->pt),double(3.99)));
                Electron_TotalPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03 - info->PileupEnergyDensity*0.2348)/probe->pt),double(3.99)));
                Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip - info->PileupEnergyDensity*0.1783)/probe->pt),double(3.99)));
                Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03+probe->GammaIso03 - info->PileupEnergyDensity*0.0686)/probe->pt),double(3.99)));
                Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip - info->PileupEnergyDensity*0.0686)/probe->pt),double(3.99)));


 
              } else {
                ElectronIsolationEfficiencyEndcapDenominator[NVertex]++;
                Electron_caloIso_Endcap[q][NVertex]->Fill(TMath::Min(double( probe->emIso03  + probe->hadIso03), double(199.99)));
                Electron_ecalIso_Endcap[q][NVertex]->Fill(TMath::Min(double( probe->emIso03 ), double(199.99)));
                Electron_hcalIso_Endcap[q][NVertex]->Fill(TMath::Min(double( probe->hadIso03 ), double(199.99)));
                Electron_caloIso04_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->emIso04  + probe->hadIso04), double(199.99)));
                Electron_ecalIso04_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->emIso04 ), double(199.99)));
                Electron_hcalIso04_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->hadIso04), double(199.99)));
                Electron_trkIso_Endcap[q][NVertex]->Fill(TMath::Min(double( probe->trkIso03 ), double(199.99)));
                Electron_trkIso04_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->trkIso04 ) , double(199.99)));
                Electron_ChargedIso_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices ) , double(199.99)));
                Electron_ChargedIsoNoPU_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03FromOtherVertices ) , double(199.99)));
                Electron_NeutralHadronIso_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->NeutralHadronIso03) , double(199.99)));
                Electron_GammaIso_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->GammaIso03) , double(199.99)));
                Electron_TotalPFIso_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03) , double(199.99)));
                Electron_FootprintRemovedPFIso_Endcap[q][NVertex]->Fill(TMath::Min( double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip) , double(199.99)));
                if (relIso < 0.1) {
                  ElectronIsolationEfficiencyEndcapNumerator[NVertex]++;              
                }
                if (relIsoL1Corrected < 0.1) {
                  ElectronL1CorrectedIsolationEfficiencyEndcapNumerator[NVertex]++;              
                }   
              }
              Electron_RelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03)/probe->pt),double(3.99)));
              Electron_TotalPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03)/probe->pt),double(3.99)));
              Electron_FootprintRemovedPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip)/probe->pt),double(3.99)));
              Electron_VertexSelectedPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03+probe->GammaIso03)/probe->pt),double(3.99)));
              Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip)/probe->pt),double(3.99)));
              Electron_RelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03 - info->PileupEnergyDensity*0.0919  )/probe->pt),double(3.99)));
              Electron_TotalPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03+probe->GammaIso03 - info->PileupEnergyDensity*0.3911)/probe->pt),double(3.99)));
              Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip - info->PileupEnergyDensity*0.2097)/probe->pt),double(3.99)));
              Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03+probe->GammaIso03 - info->PileupEnergyDensity*0.2215)/probe->pt),double(3.99)));
              Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralHadronIso03-probe->NeutralHadronIso007+probe->GammaIso03-probe->GammaIsoVetoEtaStrip - info->PileupEnergyDensity*0.2215)/probe->pt),double(3.99)));


//               Double_t HOverECut = (fabs(probe->eta) < 1.5)?0.04:0.025;
// //               Double_t HOverECut = (fabs(probe->eta) < 1.5)?0.025:0.025;
//               if (fabs(probe->eta) < 1.5) {

//                 ElectronHOverEEfficiencyBarrelDenominator[NVertex]++;
//                 if (relIso < 0.1) {
//                   if (probe->HoverE < HOverECut) {
//                     ElectronHOverEEfficiencyBarrelNumerator[NVertex]++;                     
//                   }
//                   if ((probe->HoverE*probe->e - PUHOverE)/probe->e  < HOverECut) {
//                     ElectronL1CorrectedHOverEEfficiencyBarrelNumerator[NVertex]++;                     
//                   }
//                 } 
//               } else {
//                 ElectronHOverEEfficiencyEndcapDenominator[NVertex]++;
//                 if (relIso < 0.1) {
//                   if (probe->HoverE < HOverECut) {
//                     ElectronHOverEEfficiencyEndcapNumerator[NVertex]++;                     
//                   }
//                   if ((probe->HoverE*probe->e - PUHOverE)/probe->e  < HOverECut) {
//                     ElectronL1CorrectedHOverEEfficiencyEndcapNumerator[NVertex]++;                     
//                   }
//                 } 
//               }
              
            } //passes T&P selection

              
          } //loop over probes
        } //loop over tags



        //******************************************************************************
        //loop over muon pairs
        //******************************************************************************
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *tag = (mithep::TMuon*)((*muonArr)[i]);
          if ( !(
                 fabs(tag->eta) < 2.5
                 && tag->pt > 20.0
                 && passMuonTagCuts(tag)
                 )
            ) continue;
          
          for(Int_t j=0; j<muonArr->GetEntries(); j++) {
            if (i==j) continue;
            
            const mithep::TMuon *probe = (mithep::TMuon*)((*muonArr)[j]);
            if ( !(
                   fabs(probe->eta) < 2.5
                   && probe->pt > 10.0
                   && passMuonProbeCuts(probe)
                   )
              ) continue;
            
            
            mithep::FourVectorM tagVector;
            mithep::FourVectorM probeVector;
            tagVector.SetCoordinates(tag->pt, tag->eta, tag->phi, 0.51099892e-3 );
            probeVector.SetCoordinates(probe->pt, probe->eta, probe->phi, 0.51099892e-3 );
            mithep::FourVectorM dilepton = tagVector+probeVector;
            
            if (dilepton.M() > 75 && dilepton.M() < 105
                
              ) {

              Rho[q][NVertex]->Fill(info->PileupEnergyDensity);

              MuonIsolationEfficiencyDenominator[NVertex]++;  

              if (fabs(probe->eta) < 1.5) {
                Muon_caloIso_Barrel[q][NVertex]->Fill( TMath::Min(double( probe->emIso03 + probe->hadIso03), double(199.99)));
                Muon_ecalIso_Barrel[q][NVertex]->Fill( TMath::Min(double( probe->emIso03 ), double(199.99)) );
                Muon_hcalIso_Barrel[q][NVertex]->Fill( TMath::Min(double( probe->hadIso03), double(199.99)));
                Muon_caloIso05_Barrel[q][NVertex]->Fill( TMath::Min( double(probe->emIso05 + probe->hadIso05), double(199.99)));
                Muon_ecalIso05_Barrel[q][NVertex]->Fill( TMath::Min( double(probe->emIso05 ), double(199.99)) );
                Muon_hcalIso05_Barrel[q][NVertex]->Fill( TMath::Min( double(probe->hadIso05), double(199.99)) );
                Muon_trkIso_Barrel[q][NVertex]->Fill( TMath::Min(double( probe->trkIso03 ), double(199.99)) );
                Muon_trkIso05_Barrel[q][NVertex]->Fill( probe->trkIso05 );
                Muon_ChargedIso_Barrel[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices),double(199.99)));
                Muon_ChargedIsoNoPU_Barrel[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03),double(199.99)));
                Muon_NeutralIso_Barrel[q][NVertex]->Fill(TMath::Min(double(probe->NeutralIso03) , double(199.99)));
                Muon_TotalPFIso_Barrel[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03) , double(199.99)));

                Muon_RelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03)/probe->pt) , double(199.99)));
                Muon_TotalPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03)/probe->pt) , double(199.99)));
                Muon_VertexSelectedPFRelIso_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralIso03)/probe->pt) , double(199.99)));
                Muon_RelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03 - info->PileupEnergyDensity*0.0964)/probe->pt) , double(199.99)));
                Muon_TotalPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03 - info->PileupEnergyDensity*0.2258)/probe->pt) , double(199.99)));
                Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralIso03 - info->PileupEnergyDensity*0.0602)/probe->pt) , double(199.99)));




              } else {
                Muon_caloIso_Endcap[q][NVertex]->Fill( TMath::Min(double( probe->emIso03 + probe->hadIso03), double(199.99)));
                Muon_ecalIso_Endcap[q][NVertex]->Fill( TMath::Min(double( probe->emIso03 ), double(199.99)) );
                Muon_hcalIso_Endcap[q][NVertex]->Fill( TMath::Min(double( probe->hadIso03), double(199.99)));
                Muon_caloIso05_Endcap[q][NVertex]->Fill( TMath::Min( double(probe->emIso05 + probe->hadIso05), double(199.99)));
                Muon_ecalIso05_Endcap[q][NVertex]->Fill( TMath::Min( double(probe->emIso05 ), double(199.99)) );
                Muon_hcalIso05_Endcap[q][NVertex]->Fill( TMath::Min( double(probe->hadIso05), double(199.99)) );
                Muon_trkIso_Endcap[q][NVertex]->Fill( TMath::Min(double( probe->trkIso03 ), double(199.99)) );
                Muon_trkIso05_Endcap[q][NVertex]->Fill( probe->trkIso05 );
                Muon_ChargedIso_Endcap[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices),double(199.99)));
                Muon_ChargedIsoNoPU_Endcap[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03),double(199.99)));
                Muon_NeutralIso_Endcap[q][NVertex]->Fill(TMath::Min(double(probe->NeutralIso03) , double(199.99)));
                Muon_TotalPFIso_Endcap[q][NVertex]->Fill(TMath::Min(double(probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03) , double(199.99)));
            
                Muon_RelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03)/probe->pt) , double(199.99)));
                Muon_TotalPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03)/probe->pt) , double(199.99)));
                Muon_VertexSelectedPFRelIso_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralIso03)/probe->pt) , double(199.99)));
                Muon_RelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->emIso03 + probe->hadIso03 + probe->trkIso03 - info->PileupEnergyDensity*0.0809)/probe->pt) , double(199.99)));
                Muon_TotalPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->ChargedIso03FromOtherVertices+probe->NeutralIso03 - info->PileupEnergyDensity*0.1907)/probe->pt) , double(199.99)));
                Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][NVertex]->Fill(TMath::Min(double((probe->ChargedIso03+probe->NeutralIso03 - info->PileupEnergyDensity*0.0464)/probe->pt) , double(199.99)));

              }


              if ((probe->trkIso03 + probe->emIso03 + probe->hadIso03) / probe->pt < 0.15) {
//               if ((probe->trkIso05 + probe->emIso05 + probe->hadIso05) / probe->pt < 0.15) {
                MuonIsolationEfficiencyNumerator[NVertex]++;
              }              
              if ((probe->trkIso03 + TMath::Max(probe->emIso03 + probe->hadIso03 - PUIsolationEnergy,0.0)) / probe->pt < 0.15) {
//               if ((probe->trkIso05 + TMath::Max(probe->emIso05 + probe->hadIso05 - PUIsolationEnergy05Cone,0.0)) / probe->pt < 0.15) {
                MuonL1CorrectedIsolationEfficiencyNumerator[NVertex]++;
              }
              
              
            } //passes T&P selection
            
            
          } //loop over probes
        } //loop over tags
                        
      } //end loop over data     

      delete info;
      delete electronArr;

    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make Efficiency plots
  //==============================================================================================================


  const int nPoints = ElectronIsolationEfficiencyBarrelNumerator.size();
  double NPileup[nPoints];
  double NPileupError[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    NPileup[i] = i;
    NPileupError[i] = 0.5;     
  }


  double ElectronIsolationBarrelEff[nPoints];
  double ElectronIsolationBarrelEffErrLow[nPoints];
  double ElectronIsolationBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationBarrelEff[i] = ratio;
    ElectronIsolationBarrelEffErrLow[i] = errLow;
    ElectronIsolationBarrelEffErrHigh[i] = errHigh;
    cout << "ElectronIsolationBarrelEff " << i << " : " << ElectronIsolationBarrelEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationBarrelEff, NPileupError, NPileupError, ElectronIsolationBarrelEffErrLow, ElectronIsolationBarrelEffErrHigh);
  ElectronIsolationBarrelEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedIsolationBarrelEff[nPoints];
  double ElectronL1CorrectedIsolationBarrelEffErrLow[nPoints];
  double ElectronL1CorrectedIsolationBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedIsolationBarrelEff[i] = ratio;
    ElectronL1CorrectedIsolationBarrelEffErrLow[i] = errLow;
    ElectronL1CorrectedIsolationBarrelEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedIsolationBarrelEff " << i << " : " << ElectronL1CorrectedIsolationBarrelEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedIsolationBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationBarrelEff,NPileupError, NPileupError, ElectronL1CorrectedIsolationBarrelEffErrLow, ElectronL1CorrectedIsolationBarrelEffErrHigh);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);



  double ElectronIsolationEndcapEff[nPoints];
  double ElectronIsolationEndcapEffErrLow[nPoints];
  double ElectronIsolationEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronIsolationEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronIsolationEndcapEff[i] = ratio;
    ElectronIsolationEndcapEffErrLow[i] = errLow;
    ElectronIsolationEndcapEffErrHigh[i] = errHigh;
    cout << "ElectronIsolationEndcapEff " << i << " : " << ElectronIsolationEndcapEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronIsolationEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronIsolationEndcapEff, NPileupError, NPileupError, ElectronIsolationEndcapEffErrLow, ElectronIsolationEndcapEffErrHigh);
  ElectronIsolationEndcapEffVsNPileup->SetMarkerColor(kRed);
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedIsolationEndcapEff[nPoints];
  double ElectronL1CorrectedIsolationEndcapEffErrLow[nPoints];
  double ElectronL1CorrectedIsolationEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedIsolationEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronIsolationEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedIsolationEndcapEff[i] = ratio;
    ElectronL1CorrectedIsolationEndcapEffErrLow[i] = errLow;
    ElectronL1CorrectedIsolationEndcapEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedIsolationEndcapEff " << i << " : " << ElectronL1CorrectedIsolationEndcapEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedIsolationEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedIsolationEndcapEff,NPileupError, NPileupError, ElectronL1CorrectedIsolationEndcapEffErrLow, ElectronL1CorrectedIsolationEndcapEffErrHigh);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);




  double ElectronHOverEBarrelEff[nPoints];
  double ElectronHOverEBarrelEffErrLow[nPoints];
  double ElectronHOverEBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronHOverEEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronHOverEBarrelEff[i] = ratio;
    ElectronHOverEBarrelEffErrLow[i] = errLow;
    ElectronHOverEBarrelEffErrHigh[i] = errHigh;
    cout << "ElectronHOverEBarrelEff " << i << " : " << ElectronHOverEBarrelEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronHOverEBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronHOverEBarrelEff, NPileupError, NPileupError, ElectronHOverEBarrelEffErrLow, ElectronHOverEBarrelEffErrHigh);
  ElectronHOverEBarrelEffVsNPileup->SetMarkerColor(kRed);
  ElectronHOverEBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronHOverEBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedHOverEBarrelEff[nPoints];
  double ElectronL1CorrectedHOverEBarrelEffErrLow[nPoints];
  double ElectronL1CorrectedHOverEBarrelEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedHOverEEfficiencyBarrelNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyBarrelDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedHOverEBarrelEff[i] = ratio;
    ElectronL1CorrectedHOverEBarrelEffErrLow[i] = errLow;
    ElectronL1CorrectedHOverEBarrelEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedHOverEBarrelEff " << i << " : " << ElectronL1CorrectedHOverEBarrelEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedHOverEBarrelEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedHOverEBarrelEff,NPileupError, NPileupError, ElectronL1CorrectedHOverEBarrelEffErrLow, ElectronL1CorrectedHOverEBarrelEffErrHigh);
  ElectronL1CorrectedHOverEBarrelEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedHOverEBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedHOverEBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);



  double ElectronHOverEEndcapEff[nPoints];
  double ElectronHOverEEndcapEffErrLow[nPoints];
  double ElectronHOverEEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronHOverEEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronHOverEEndcapEff[i] = ratio;
    ElectronHOverEEndcapEffErrLow[i] = errLow;
    ElectronHOverEEndcapEffErrHigh[i] = errHigh;
    cout << "ElectronHOverEEndcapEff " << i << " : " << ElectronHOverEEndcapEff[i] << endl;
  }
  TGraphAsymmErrors *ElectronHOverEEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronHOverEEndcapEff, NPileupError, NPileupError, ElectronHOverEEndcapEffErrLow, ElectronHOverEEndcapEffErrHigh);
  ElectronHOverEEndcapEffVsNPileup->SetMarkerColor(kRed);
  ElectronHOverEEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronHOverEEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double ElectronL1CorrectedHOverEEndcapEff[nPoints];
  double ElectronL1CorrectedHOverEEndcapEffErrLow[nPoints];
  double ElectronL1CorrectedHOverEEndcapEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(ElectronL1CorrectedHOverEEfficiencyEndcapNumerator[i]);
    Double_t n2 = TMath::Nint(ElectronHOverEEfficiencyEndcapDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    ElectronL1CorrectedHOverEEndcapEff[i] = ratio;
    ElectronL1CorrectedHOverEEndcapEffErrLow[i] = errLow;
    ElectronL1CorrectedHOverEEndcapEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "ElectronL1CorrectedHOverEEndcapEff " << i << " : " << ElectronL1CorrectedHOverEEndcapEff[i] << endl;

  }
  TGraphAsymmErrors *ElectronL1CorrectedHOverEEndcapEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, ElectronL1CorrectedHOverEEndcapEff,NPileupError, NPileupError, ElectronL1CorrectedHOverEEndcapEffErrLow, ElectronL1CorrectedHOverEEndcapEffErrHigh);
  ElectronL1CorrectedHOverEEndcapEffVsNPileup->SetMarkerColor(kBlue);
  ElectronL1CorrectedHOverEEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  ElectronL1CorrectedHOverEEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);



  double MuonIsolationEff[nPoints];
  double MuonIsolationEffErrLow[nPoints];
  double MuonIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonIsolationEff[i] = ratio;
    MuonIsolationEffErrLow[i] = errLow;
    MuonIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonIsolationEff " << i << " : " << MuonIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonIsolationEff, NPileupError, NPileupError, MuonIsolationEffErrLow, MuonIsolationEffErrHigh);
  MuonIsolationEffVsNPileup->SetMarkerColor(kRed);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);

  double MuonL1CorrectedIsolationEff[nPoints];
  double MuonL1CorrectedIsolationEffErrLow[nPoints];
  double MuonL1CorrectedIsolationEffErrHigh[nPoints];
  for (UInt_t i=0; i<nPoints; ++i) {
    
    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    
    Double_t n1 = TMath::Nint(MuonL1CorrectedIsolationEfficiencyNumerator[i]);
    Double_t n2 = TMath::Nint(MuonIsolationEfficiencyDenominator[i]);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    MuonL1CorrectedIsolationEff[i] = ratio;
    MuonL1CorrectedIsolationEffErrLow[i] = errLow;
    MuonL1CorrectedIsolationEffErrHigh[i] = errHigh;
    NPileup[i] = i;
    cout << "MuonL1CorrectedIsolationEff " << i << " : " << MuonL1CorrectedIsolationEff[i] << endl;
  }
  TGraphAsymmErrors *MuonL1CorrectedIsolationEffVsNPileup = new TGraphAsymmErrors (nPoints,  NPileup, MuonL1CorrectedIsolationEff, NPileupError, NPileupError,MuonL1CorrectedIsolationEffErrLow,MuonL1CorrectedIsolationEffErrHigh  );
  MuonL1CorrectedIsolationEffVsNPileup->SetMarkerColor(kBlue);
  MuonL1CorrectedIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.02);
  MuonL1CorrectedIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);




  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  Double_t ymin = 0.0;
  Double_t ymax = 1.1;


  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TLegend *legend = 0;

  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.2;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronIsolationBarrelEffVsNPileup, "NoCorrection", "LP");
//     legend->AddEntry(ElectronL1CorrectedIsolationBarrelEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronIsolationBarrelEffVsNPileup->SetTitle("");
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationBarrelEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronIsolationBarrelEffVsNPileup->Draw("AP");
//   ElectronL1CorrectedIsolationBarrelEffVsNPileup->Draw("Psame");
  ElectronIsolationBarrelEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs(("ElectronIsolationEfficiency_Barrel_vs_NVertices"+label+".gif").c_str());

  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.2;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronIsolationEndcapEffVsNPileup, "NoCorrection", "LP");
//     legend->AddEntry(ElectronL1CorrectedIsolationEndcapEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronIsolationEndcapEffVsNPileup->SetTitle("");
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronIsolationEndcapEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronIsolationEndcapEffVsNPileup->Draw("AP");
//   ElectronL1CorrectedIsolationEndcapEffVsNPileup->Draw("Psame");
  ElectronIsolationEndcapEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs(("ElectronIsolationEfficiency_Endcap_vs_NVertices"+label+".gif").c_str());


  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.2;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronHOverEBarrelEffVsNPileup, "NoCorrection", "LP");
//     legend->AddEntry(ElectronL1CorrectedHOverEBarrelEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronHOverEBarrelEffVsNPileup->SetTitle("");
  ElectronHOverEBarrelEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronHOverEBarrelEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronHOverEBarrelEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronHOverEBarrelEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronHOverEBarrelEffVsNPileup->Draw("AP");
//   ElectronL1CorrectedHOverEBarrelEffVsNPileup->Draw("Psame");
  ElectronHOverEBarrelEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs(("ElectronHOverEBarrelEfficiency_PUStudy_vs_NVertices"+label+".gif").c_str());

  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.0;
  ymax = 1.2;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(ElectronHOverEEndcapEffVsNPileup, "NoCorrection", "LP");
//     legend->AddEntry(ElectronL1CorrectedHOverEEndcapEffVsNPileup, "FastJetCorrected", "LP");
  }
  ElectronHOverEEndcapEffVsNPileup->SetTitle("");
  ElectronHOverEEndcapEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  ElectronHOverEEndcapEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  ElectronHOverEEndcapEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  ElectronHOverEEndcapEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  ElectronHOverEEndcapEffVsNPileup->Draw("AP");
//   ElectronL1CorrectedHOverEEndcapEffVsNPileup->Draw("Psame");
  ElectronHOverEEndcapEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);
  legend->Draw();
  cv->SaveAs(("ElectronHOverEEndcapEfficiency_PUStudy_vs_NVertices"+label+".gif").c_str());


  legend = new TLegend(0.20,0.25,0.43,0.40);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
  ymin = 0.8;
  ymax = 1.2;
  for (int q = 0; q<processNames.size() ; ++q) { 
    legend->AddEntry(MuonIsolationEffVsNPileup, "NoCorrection", "LP");
//     legend->AddEntry(MuonL1CorrectedIsolationEffVsNPileup, "FastJetCorrected", "LP");
  }
  MuonIsolationEffVsNPileup->SetTitle("");
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetYaxis()->SetTitle("Efficiency");
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitleOffset(1.05);
  MuonIsolationEffVsNPileup->GetXaxis()->SetTitle("Number of Reco Vertices (DA)");
  MuonIsolationEffVsNPileup->Draw("AP");
//   MuonL1CorrectedIsolationEffVsNPileup->Draw("Psame");
  MuonIsolationEffVsNPileup->GetYaxis()->SetRangeUser(ymin,ymax);

  legend->Draw();
  cv->SaveAs(("MuonIsolationEfficiency_PUStudy_vs_NVertices"+label+".gif").c_str());


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HwwSelectionPlots_LeptonEfficiency.root", "UPDATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
    file->WriteTObject(NVertices[q],NVertices[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(ElectronIsolationBarrelEffVsNPileup, ("ElectronIsolationBarrelEffVsNVertices_" + processNames[q]+label).c_str() , "WriteDelete") ;
//     file->WriteTObject(ElectronIsolationEndcapEffVsNPileup, ("ElectronIsolationEndcapEffVsNVertices_" + processNames[q]+label).c_str() , "WriteDelete") ;
     file->WriteTObject(MuonIsolationEffVsNPileup, ("MuonIsolationEffVsNVertices_" + processNames[q]+label).c_str() , "WriteDelete") ;


    for (int n=0; n < 20 ; ++n) {
      file->WriteTObject(Rho[q][n],Rho[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Barrel[q][n],Electron_caloIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso_Barrel[q][n],Electron_ecalIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso_Barrel[q][n],Electron_hcalIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso_Endcap[q][n],Electron_caloIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso_Endcap[q][n],Electron_ecalIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso_Endcap[q][n],Electron_hcalIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Barrel[q][n],Electron_trkIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso_Endcap[q][n],Electron_trkIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ChargedIso_Barrel[q][n],Electron_ChargedIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ChargedIsoNoPU_Barrel[q][n],Electron_ChargedIsoNoPU_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralHadronIso_Barrel[q][n],Electron_NeutralHadronIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_GammaIso_Barrel[q][n],Electron_GammaIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFIso_Barrel[q][n],Electron_TotalPFIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFIso_Barrel[q][n],Electron_FootprintRemovedPFIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ChargedIso_Endcap[q][n],Electron_ChargedIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ChargedIsoNoPU_Endcap[q][n],Electron_ChargedIsoNoPU_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_NeutralHadronIso_Endcap[q][n],Electron_NeutralHadronIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_GammaIso_Endcap[q][n],Electron_GammaIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFIso_Endcap[q][n],Electron_TotalPFIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFIso_Endcap[q][n],Electron_FootprintRemovedPFIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso_Barrel[q][n],Muon_caloIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso_Barrel[q][n],Muon_ecalIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso_Barrel[q][n],Muon_hcalIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso_Endcap[q][n],Muon_caloIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso_Endcap[q][n],Muon_ecalIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso_Endcap[q][n],Muon_hcalIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso_Barrel[q][n],Muon_trkIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso_Endcap[q][n],Muon_trkIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Barrel[q][n],Muon_ChargedIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIsoNoPU_Barrel[q][n],Muon_ChargedIsoNoPU_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Barrel[q][n],Muon_NeutralIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Barrel[q][n],Muon_TotalPFIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIso_Endcap[q][n],Muon_ChargedIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ChargedIsoNoPU_Endcap[q][n],Muon_ChargedIsoNoPU_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_NeutralIso_Endcap[q][n],Muon_NeutralIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFIso_Endcap[q][n],Muon_TotalPFIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Barrel[q][n],Electron_caloIso04_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso04_Barrel[q][n],Electron_ecalIso04_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso04_Barrel[q][n],Electron_hcalIso04_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_caloIso04_Endcap[q][n],Electron_caloIso04_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_ecalIso04_Endcap[q][n],Electron_ecalIso04_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_hcalIso04_Endcap[q][n],Electron_hcalIso04_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Barrel[q][n],Electron_trkIso04_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_trkIso04_Endcap[q][n],Electron_trkIso04_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso05_Barrel[q][n],Muon_caloIso05_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso05_Barrel[q][n],Muon_ecalIso05_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso05_Barrel[q][n],Muon_hcalIso05_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_caloIso05_Endcap[q][n],Muon_caloIso05_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_ecalIso05_Endcap[q][n],Muon_ecalIso05_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_hcalIso05_Endcap[q][n],Muon_hcalIso05_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso05_Barrel[q][n],Muon_trkIso05_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_trkIso05_Endcap[q][n],Muon_trkIso05_Endcap[q][n]->GetName(), "WriteDelete") ;

      file->WriteTObject(Electron_RelIso_Barrel[q][n],Electron_RelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Barrel[q][n],Electron_TotalPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFRelIso_Barrel[q][n],Electron_FootprintRemovedPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Barrel[q][n],Electron_VertexSelectedPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel[q][n],Electron_VertexSelectedFootprintRemovedPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIsoRhoCorrected_Barrel[q][n],Electron_RelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIsoRhoCorrected_Barrel[q][n],Electron_TotalPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel[q][n],Electron_FootprintRemovedPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][n],Electron_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel[q][n],Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIso_Endcap[q][n],Electron_RelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIso_Endcap[q][n],Electron_TotalPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFRelIso_Endcap[q][n],Electron_FootprintRemovedPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIso_Endcap[q][n],Electron_VertexSelectedPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap[q][n],Electron_VertexSelectedFootprintRemovedPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_RelIsoRhoCorrected_Endcap[q][n],Electron_RelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_TotalPFRelIsoRhoCorrected_Endcap[q][n],Electron_TotalPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap[q][n],Electron_FootprintRemovedPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][n],Electron_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap[q][n],Electron_VertexSelectedFootprintRemovedPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso_Barrel[q][n],Muon_RelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Barrel[q][n],Muon_TotalPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Barrel[q][n],Muon_VertexSelectedPFRelIso_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIsoRhoCorrected_Barrel[q][n],Muon_RelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIsoRhoCorrected_Barrel[q][n],Muon_TotalPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][n],Muon_VertexSelectedPFRelIsoRhoCorrected_Barrel[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIso_Endcap[q][n],Muon_RelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIso_Endcap[q][n],Muon_TotalPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIso_Endcap[q][n],Muon_VertexSelectedPFRelIso_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_RelIsoRhoCorrected_Endcap[q][n],Muon_RelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_TotalPFRelIsoRhoCorrected_Endcap[q][n],Muon_TotalPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;
      file->WriteTObject(Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][n],Muon_VertexSelectedPFRelIsoRhoCorrected_Endcap[q][n]->GetName(), "WriteDelete") ;


    }

//     file->WriteTObject(Electron_relIso_Barrel[q],Electron_relIso_Barrel[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIso_Endcap[q],Electron_relIso_Endcap[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Muon_relIso[q],Muon_relIso[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIsoL1Corrected_Barrel[q],Electron_relIsoL1Corrected_Barrel[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIsoL1Corrected_Endcap[q],Electron_relIsoL1Corrected_Endcap[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Muon_relIsoL1Corrected[q],Muon_relIsoL1Corrected[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIso04_Barrel[q],Electron_relIso04_Barrel[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIso04_Endcap[q],Electron_relIso04_Endcap[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Muon_relIso05[q],Muon_relIso05[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIso04L1Corrected_Barrel[q],Electron_relIso04L1Corrected_Barrel[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Electron_relIso04L1Corrected_Endcap[q],Electron_relIso04L1Corrected_Endcap[q]->GetName(), "WriteDelete") ;
//     file->WriteTObject(Muon_relIso05L1Corrected[q],Muon_relIso05L1Corrected[q]->GetName(), "WriteDelete") ;

  }

  file->Close();
  delete file;

  gBenchmark->Show("ElectronTagAndProbe");       


} 


void WriteMassToFile( ofstream *file, Double_t mass, Int_t eventNum) {

  (*file) << mass << " " << eventNum << endl;
}

void WriteMassToFile( ofstream *file, Double_t mass) {

  (*file) << left << mass <<  endl;
}



Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum, Int_t triggerSelection) {


  Bool_t pass = kFALSE;

  //it's electron data
  if (triggerSelection == 0) {
    if ( (ele->hltMatchBits & kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL) ) pass = kTRUE;
  } 
      
  return pass;

}


Bool_t passElectronTagCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (!passElectronCuts(ele)) pass = kFALSE;

  return pass;
}


Bool_t passElectronProbeCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
//             && ele->HoverE < 0.04
//             && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
//              && ele->HoverE < 0.025
//              && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
             && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }
  } else {
    pass = kFALSE;
    return pass;
  }

  return pass;
}

Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
            && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.025
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
             && fabs(ele->d0) < 0.02
          )
      ) {
      pass = kFALSE;
    }
  } else {
    pass = kFALSE;
    return pass;
  }

  return pass;
}

Bool_t passMuonTagCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;
  if (!passMuonCuts(mu)) pass = kFALSE;

//   //Match to HLT
//   if (!(mu->hltMatchBits & kHLT_Mu24 )) pass = kFALSE;


  return pass;
}


Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15

        && (mu->nSeg > 1 || mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  return pass;
}

Bool_t passMuonProbeCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02

        && (mu->nSeg > 1 || mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  return pass;
}


//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2)
{
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << mass;
  
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << endl;
  
}
