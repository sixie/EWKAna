//root -l -b -q EWKAna/Hww/FakeRate/ComputeMuonFakeRate_Data.C+\(\)
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
#include "EWKAna/Utils/LeptonIDCuts.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHiggs/Utils/interface/EfficiencyUtils.h"
#include "MitHiggs/Utils/interface/PlotUtils.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"

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


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu, Double_t fRho);
Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Int_t DenominatorType, string SampleType, Double_t fRho);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType, Double_t fRho);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);


//--------------------------------------------------------------------------------------------------
Bool_t passZVeto(TClonesArray *muonArr, Int_t DenominatorType, Double_t fRho)
{
  for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
    const mithep::TMuon* mu1 = (mithep::TMuon*)((*muonArr)[i]);
    if(mu1->pt        < 20)  continue;
    if(fabs(mu1->eta) > 2.4) continue;
    if(!passMuonDenominatorCuts(mu1, DenominatorType, fRho)) continue;
    for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) {
      const mithep::TMuon* mu2 = (mithep::TMuon*)((*muonArr)[j]);
      if(mu1->q == mu2->q)     continue;
      if(mu2->pt	< 20)  continue;
      if(fabs(mu2->eta) > 2.4) continue;
      if(!passMuonDenominatorCuts(mu2, DenominatorType, fRho)) continue;

      return kFALSE;  // Z candidate => fail Z veto
    }
  }
  
  return kTRUE;  // No Z candidate => pass Z veto
}


//--------------------------------------------------------------------------------------------------
void DoComputeMuonFakeRate(const string inputFilename,
                           const string label, 
                           const string outputFilename, 
                           const string smurfOutputFilename, Int_t Option);

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


//--------------------------------------------------------------------------------------------------
Double_t calcMt(const Double_t met, const Double_t metphi, const mithep::TMuon *muon)
{
  const Double_t m = 0.105659369;
  TLorentzVector vMuon; vMuon.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, m);
  TLorentzVector vMet;  vMet.SetPtEtaPhiM(met, 0, metphi, 0);
  Double_t et = (vMuon.E())*(vMuon.Pt())/(vMuon.P());
  
  return sqrt( (et+vMet.Perp())*(et+vMet.Perp()) - (vMuon.Px()+vMet.Px())*(vMuon.Px()+vMet.Px()) - (vMuon.Py()+vMet.Py())*(vMuon.Py()+vMet.Py()) );
}


//=== MAIN MACRO =================================================================================================
void ComputeMuonFakeRate_Data() {


//************************************
//Tests With Different WPs
//************************************



//************************************
//Real Stuff
//************************************
//   DoComputeMuonFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr_MuonFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV5.skim.root");

//   DoComputeMuonFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-pr_MuonFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.skim.prv2.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/2011Data/WWAnalysisSkimmed_r11a-smu-pr-v4_MuonFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.skim.prv4.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-smu-m10-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.skim.m10.root");


//   DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.SmurfV7.root", "FakeRates_Muon_SmurfV7.root", 0);
//   DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.CutBasedWithEACorrPFIso.root", "FakeRates_Muon_CutBasedWithEACorrPFIso.root", 0);


//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-m10-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.m10.skim.root", "FakeRates_SmurfV6.m10.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-pr-v4_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.prv4.skim.root", "FakeRates_SmurfV6.prv4.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-a05-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.a05.skim.root", "FakeRates_SmurfV6.a05.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-pr-v6_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.prv6.skim.root", "FakeRates_SmurfV6.prv6.root");
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-smu-pr-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.SmurfV6.r11b.skim.root", "FakeRates_SmurfV6.r11b.root");


  DoComputeMuonFakeRate("LIST","MuonFakeRate","MuonFakeRate.BDTGIDIsoCombined.root", "FakeRates_Muon_BDTGIDIsoCombined.root", 10);

//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-m10-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.BDTGV3.r11am10.root", "FakeRates_Muon_BDTGV3.r11am10.root", 10);
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-pr-v4_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.BDTGV3.r11aprv4.root", "FakeRates_Muon_BDTGV3.r11aprv4.root", 10);
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-a05-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.BDTGV3.r11aa05.root", "FakeRates_Muon_BDTGV3.r11aa05.root", 10);
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-o03-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.BDTGV3.r11ao03.root", "FakeRates_Muon_BDTGV3.r11ao03.root", 10);
//   DoComputeMuonFakeRate("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-smu-pr-v1_MuFakeRateTriggerAndDenominatorSkim.root","MuonFakeRate","MuonFakeRate.BDTGV3.r11b.root", "FakeRates_Muon_BDTGV3.r11b.root", 10);


}



void DoComputeMuonFakeRate(const string inputFilename,
                           const string label, 
                           const string outputFilename, 
                           const string smurfOutputFilename,
                           Int_t Option)
{  
  gBenchmark->Start("WWTemplate");


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  mithep::MuonIDMVA *muonIDBDTGIDIsoCombined = 0;
  if (Option >= 10 && Option <= 19) {

    muonIDBDTGIDIsoCombined = new mithep::MuonIDMVA();
    muonIDBDTGIDIsoCombined->Initialize("BDTG method",
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml", 
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml", 
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml", 
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml", 
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml", 
                                    "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml",
                                    mithep::MuonIDMVA::kV8);
    
  }
  

  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
//   ptbins.push_back(10);  
//   ptbins.push_back(80);  


  ptbins.push_back(10);  
  ptbins.push_back(11);  
  ptbins.push_back(12);  
  ptbins.push_back(13);  
  ptbins.push_back(14);  
  ptbins.push_back(15);  
  ptbins.push_back(16);  
  ptbins.push_back(17);  
  ptbins.push_back(18);  
  ptbins.push_back(19);  
  ptbins.push_back(20);  
  ptbins.push_back(22);  
  ptbins.push_back(24);  
  ptbins.push_back(26);  
  ptbins.push_back(28);
  ptbins.push_back(30);  
  ptbins.push_back(32.5);  
  ptbins.push_back(35);  
  ptbins.push_back(37.5);  
  ptbins.push_back(40);  
  ptbins.push_back(45);  
  ptbins.push_back(50);  
  ptbins.push_back(60);  
  ptbins.push_back(70);  
  ptbins.push_back(80);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(0.0);
  etabins.push_back(0.2);
  etabins.push_back(0.4);
  etabins.push_back(0.6);
  etabins.push_back(0.8);
  etabins.push_back(1.0);
  etabins.push_back(1.2);
  etabins.push_back(1.4);
  etabins.push_back(1.6);
  etabins.push_back(1.8); 
  etabins.push_back(2.0);
  etabins.push_back(2.25);
  etabins.push_back(2.5);

  vector<double> phibins;
  phibins.push_back(-3.25);
  phibins.push_back(-2.75);
  phibins.push_back(-2.25);
  phibins.push_back(-1.75);
  phibins.push_back(-1.25);
  phibins.push_back(-0.75);
  phibins.push_back(-0.25);
  phibins.push_back(0.25);
  phibins.push_back(0.75);
  phibins.push_back(1.25);
  phibins.push_back(1.75);
  phibins.push_back(2.25);
  phibins.push_back(2.75);
  phibins.push_back(3.25);


  vector<double> nvtxbins;
  nvtxbins.push_back(0);
  nvtxbins.push_back(2);
  nvtxbins.push_back(4);
  nvtxbins.push_back(6);
  nvtxbins.push_back(8);
  nvtxbins.push_back(10);
  nvtxbins.push_back(12);
  nvtxbins.push_back(14);
  nvtxbins.push_back(16);
  nvtxbins.push_back(20);
  nvtxbins.push_back(22);
  nvtxbins.push_back(24);


  vector<double> rhobins;
  rhobins.push_back(0);
  rhobins.push_back(2);
  rhobins.push_back(4);
  rhobins.push_back(6);
  rhobins.push_back(8);
  rhobins.push_back(10);
  rhobins.push_back(12);
  rhobins.push_back(14);
  rhobins.push_back(16);
  rhobins.push_back(18);
  rhobins.push_back(20);
  rhobins.push_back(22);
  rhobins.push_back(24);
  rhobins.push_back(26);
  rhobins.push_back(28);
  rhobins.push_back(30);



  vector<double> ptbins2D;
  ptbins2D.push_back(10);  
  ptbins2D.push_back(15);  
  ptbins2D.push_back(20);  
  ptbins2D.push_back(25);  
  ptbins2D.push_back(30);  
  ptbins2D.push_back(35);  
  vector<double> etabins2D;
  etabins2D.push_back(0.0);
  etabins2D.push_back(1.0);
  etabins2D.push_back(1.479);
  etabins2D.push_back(2.0);
  etabins2D.push_back(2.5);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *leadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH1F*> > > DenominatorVector_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Rho;
  vector<vector<vector<TH2F*> > > DenominatorVector_PtEta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH1F*> > > NumeratorVector_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Rho;
  vector<vector<vector<TH2F*> > > NumeratorVector_PtEta;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_Rho;


  vector<Double_t> denominatorType;
  denominatorType.push_back(1);
  denominatorType.push_back(2);
  denominatorType.push_back(3);
  vector<string> sampleLabel;
  sampleLabel.push_back("Mu8Sample");
  sampleLabel.push_back("Mu15Sample");
  sampleLabel.push_back("Mu8PtCombinedSample");
  sampleLabel.push_back("Mu8Jet40Sample");
  sampleLabel.push_back("PhotonJetsSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(0);
  ptThreshold.push_back(5);
  ptThreshold.push_back(10);
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH1F*> > tmpDenominatorVector_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Rho;
    vector<vector<TH2F*> > tmpDenominatorVector_PtEta;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH1F*> > tmpNumeratorVector_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Rho;
    vector<vector<TH2F*> > tmpNumeratorVector_PtEta;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_Rho;


    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH1F*>  tmptmpDenominatorVector_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Rho;
      vector<TH2F*>  tmptmpDenominatorVector_PtEta;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH1F*>  tmptmpNumeratorVector_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Rho;
      vector<TH2F*>  tmptmpNumeratorVector_PtEta;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho;

      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);

        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0.0, 3.0);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0.0, 3.0);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);


   
        TH1F *histDenominator_Pt10To20_Barrel_NVtx = new TH1F(("histDenominator_Pt10To20_Barrel_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Barrel_Rho = new TH1F(("histDenominator_Pt10To20_Barrel_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt10To20_Endcap_NVtx = new TH1F(("histDenominator_Pt10To20_Endcap_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Endcap_Rho = new TH1F(("histDenominator_Pt10To20_Endcap_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Barrel_NVtx = new TH1F(("histDenominator_Pt20ToInf_Barrel_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Barrel_Rho = new TH1F(("histDenominator_Pt20ToInf_Barrel_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Endcap_NVtx = new TH1F(("histDenominator_Pt20ToInf_Endcap_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Endcap_Rho = new TH1F(("histDenominator_Pt20ToInf_Endcap_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Barrel_NVtx = new TH1F(("histNumerator_Pt10To20_Barrel_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Barrel_Rho = new TH1F(("histNumerator_Pt10To20_Barrel_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Endcap_NVtx = new TH1F(("histNumerator_Pt10To20_Endcap_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Endcap_Rho = new TH1F(("histNumerator_Pt10To20_Endcap_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Barrel_NVtx = new TH1F(("histNumerator_Pt20ToInf_Barrel_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Barrel_Rho = new TH1F(("histNumerator_Pt20ToInf_Barrel_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Endcap_NVtx = new TH1F(("histNumerator_Pt20ToInf_Endcap_NVtx_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Endcap_Rho = new TH1F(("histNumerator_Pt20ToInf_Endcap_Rho_DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);



        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_NVtx.push_back(histDenominator_NVtx);
        tmptmpDenominatorVector_Rho.push_back(histDenominator_Rho);
        tmptmpDenominatorVector_PtEta.push_back(histDenominator_PtEta);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_NVtx.push_back(histNumerator_NVtx);
        tmptmpNumeratorVector_Rho.push_back(histNumerator_Rho);
        tmptmpNumeratorVector_PtEta.push_back(histNumerator_PtEta);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);
        tmptmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(histDenominator_Pt10To20_Barrel_NVtx);
        tmptmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(histDenominator_Pt10To20_Barrel_Rho);
        tmptmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(histDenominator_Pt10To20_Endcap_NVtx);
        tmptmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(histDenominator_Pt10To20_Endcap_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(histDenominator_Pt20ToInf_Barrel_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(histDenominator_Pt20ToInf_Barrel_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(histDenominator_Pt20ToInf_Endcap_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(histDenominator_Pt20ToInf_Endcap_Rho);
        tmptmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(histNumerator_Pt10To20_Barrel_NVtx);
        tmptmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(histNumerator_Pt10To20_Barrel_Rho);
        tmptmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(histNumerator_Pt10To20_Endcap_NVtx);
        tmptmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(histNumerator_Pt10To20_Endcap_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(histNumerator_Pt20ToInf_Barrel_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(histNumerator_Pt20ToInf_Barrel_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(histNumerator_Pt20ToInf_Endcap_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(histNumerator_Pt20ToInf_Endcap_Rho);

      }
        tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
        tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
        tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
        tmpDenominatorVector_NVtx.push_back(tmptmpDenominatorVector_NVtx);
        tmpDenominatorVector_Rho.push_back(tmptmpDenominatorVector_Rho);
        tmpDenominatorVector_PtEta.push_back(tmptmpDenominatorVector_PtEta);
        tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
        tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
        tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
        tmpNumeratorVector_NVtx.push_back(tmptmpNumeratorVector_NVtx);
        tmpNumeratorVector_Rho.push_back(tmptmpNumeratorVector_Rho);
        tmpNumeratorVector_PtEta.push_back(tmptmpNumeratorVector_PtEta);
        tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
        tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
        tmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_NVtx);
        tmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_Rho);
        tmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_NVtx);
        tmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_Rho);
        tmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
        tmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho);
        tmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
        tmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho);
        tmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_NVtx);
        tmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_Rho);
        tmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_NVtx);
        tmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_Rho);
        tmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
        tmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho);
        tmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
        tmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_NVtx.push_back(tmpDenominatorVector_NVtx);
    DenominatorVector_Rho.push_back(tmpDenominatorVector_Rho);
    DenominatorVector_PtEta.push_back(tmpDenominatorVector_PtEta);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_NVtx.push_back(tmpNumeratorVector_NVtx);
    NumeratorVector_Rho.push_back(tmpNumeratorVector_Rho);
    NumeratorVector_PtEta.push_back(tmpNumeratorVector_PtEta);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
    DenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmpDenominatorVector_Pt10To20_Barrel_NVtx);
    DenominatorVector_Pt10To20_Barrel_Rho.push_back(tmpDenominatorVector_Pt10To20_Barrel_Rho);
    DenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmpDenominatorVector_Pt10To20_Endcap_NVtx);
    DenominatorVector_Pt10To20_Endcap_Rho.push_back(tmpDenominatorVector_Pt10To20_Endcap_Rho);
    DenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
    DenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_Rho);
    DenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
    DenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_Rho);
    NumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmpNumeratorVector_Pt10To20_Barrel_NVtx);
    NumeratorVector_Pt10To20_Barrel_Rho.push_back(tmpNumeratorVector_Pt10To20_Barrel_Rho);
    NumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmpNumeratorVector_Pt10To20_Endcap_NVtx);
    NumeratorVector_Pt10To20_Endcap_Rho.push_back(tmpNumeratorVector_Pt10To20_Endcap_Rho);
    NumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
    NumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_Rho);
    NumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
    NumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_Rho);
  }

  ofstream eventListFile("eventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TMuon");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-m10-v1_MuFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-pr-v4_MuFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-a05-v1_MuFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-smu-o03-v1_MuFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-smu-pr-v1_MuFakeRateTriggerAndDenominatorSkim.root");    
  } else {
    inputfiles.push_back(inputFilename);
  }
  
  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *photonBr;


    //*****************************************************************************************
    //Loop over Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon", &electronArr); electronBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      if (Option >= 10 && Option < 20) {
        if (info->evtNum % 2 == 0) continue;
      }

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
      Double_t rho = 0;
      if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      jetArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);


      //********************************************************
      // TcMet
      //********************************************************
      TVector3 pfMet;        
      if(info->pfMEx!=0 || info->pfMEy!=0) {       
        pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
      }
      Double_t met = pfMet.Pt();

      Int_t NMuons = muonArr->GetEntries();
      
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
      
        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = mu->pt + mu->ChargedIso04 + mu->NeutralIso04_10Threshold;

//         cout << "Event: " << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
//         for(Int_t e=0; e<electronArr->GetEntries(); e++) {
//           const mithep::TMuon *muon = (mithep::TMuon*)((*electronArr)[e]);
//           cout << "Mu " << e << " : " << muon->pt << " " << muon->eta << " " << muon->phi << " " 
//                << mu->ChargedIso03 + mu->NeutralIso03_10Threshold << " "
//                << endl;
//         }
//         for(Int_t j=0; j<jetArr->GetEntries(); j++) {
//           const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
//           cout << "Jet " << j << " : " << jet->pt << " " << jet->eta << " " << jet->phi << endl;
//         }
        }


        //********************************************************
        // Photons
        //********************************************************
        Int_t NPhotons = 0;      
        Double_t photonEt = -1;
        for(Int_t p=0; p<photonArr->GetEntries(); p++) {
          const mithep::TPhoton *photon = (mithep::TPhoton*)((*photonArr)[p]);
        
          Bool_t isMu = kFALSE;
          for(Int_t e=0; e<electronArr->GetEntries(); e++) {
            const mithep::TMuon *tmpMu = (mithep::TMuon*)((*electronArr)[e]);   
            if (mithep::MathUtils::DeltaR(photon->phi, photon->eta, tmpMu->phi, tmpMu->eta) < 0.3) isMu = kTRUE;
          }
 


          //photon ID
          if ( photon->et > 20
               && !isMu
               && !photon->hasPixelSeed
               && photon->emIso03 < 2.0 + 0.006*photon->et
               && photon->hadIso03 < 2.0+0.0025*photon->et
               && photon->trkIso03Hollow < 1.5 + 0.001*photon->et
               && ( (fabs(photon->eta) < 1.5 && photon->sigiEtaiEta < 0.01) || (fabs(photon->eta) >= 1.5 && photon->sigiEtaiEta < 0.028) )
            ) {
            continue;
          }


          mithep::FourVectorM phFourVector;
          mithep::FourVectorM muFourVector;
          muFourVector.SetCoordinates(mu->pt, mu->eta, mu->phi, 0.51099892e-3 );
          phFourVector.SetCoordinates(photon->et, photon->eta, photon->phi, 0.0 );
          mithep::FourVectorM dilepton = phFourVector+muFourVector;
        
//         cout << "photon " << p << " : " << photon->et << " " << photon->eta << " " << photon->phi 
//              << "Mu : " << mu->pt << " " << mu->eta << " " << mu->phi << " : " 
//              << dilepton.M() << endl;
        
          if ( fabs(dilepton.M() - 91) > 20 ) {
            photonEt = photon->et;
            NPhotons++;
          }
        }
      
      
        for ( UInt_t denominatorTypeIndex = 0 ; denominatorTypeIndex < denominatorType.size() ; ++denominatorTypeIndex ) {
          for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
            for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
          
              //********************************************************
              // Event Selection Cuts
              //********************************************************

              if (NMuons > 1) continue;

              if (sampleLabel[sampleTypeIndex] == "Mu8Sample" || sampleLabel[sampleTypeIndex] == "Mu15Sample" || sampleLabel[sampleTypeIndex] == "Mu8PtCombinedSample" || sampleLabel[sampleTypeIndex] == "Mu8Jet40Sample") {
                if (met > 20) continue;

                if (calcMt(pfMet.Pt(), pfMet.Phi(), mu) > 20) continue;

                if (!passZVeto(muonArr, denominatorType[denominatorTypeIndex], rho)) continue;

                Bool_t passJetSelection = kFALSE;
                if (ptThreshold[ptThresholdIndex] == 0) passJetSelection = kTRUE;
                if (leadingJetPt > ptThreshold[ptThresholdIndex]) {
                  passJetSelection = kTRUE;               
                }
//               cout << leadingJetPt << " " << ptThreshold[ptThresholdIndex] << " " << passJetSelection << endl;
                if (!passJetSelection) continue;
              }
      
              if (sampleLabel[sampleTypeIndex] == "PhotonJetsSample") {
                if (NPhotons != 1) continue;
                if (!(photonEt > ptThreshold[ptThresholdIndex])) continue;  
              }

              if (!(mu->pt > 10 && mu->pt < 35)) continue;

              if (passMuonDenominatorCuts(info->triggerBits, mu, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex], rho)) {
              

//                 eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << mu->pt << " " << mu->eta << " " << mu->phi << " : " << passMuonNumeratorCuts(mu) << endl;
              

                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt,eventweight);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta),eventweight);
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi,eventweight);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, fabs(mu->eta), eventweight);
                
                if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                } 
              
                if (mu->pt < 30) {
                  Double_t iso = ( mu->ChargedIso03 + mu->NeutralIso03_10Threshold ) / mu->pt;
                  LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(leptonJetPt, 0.01),99.9),eventweight);              
                  DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(iso, 0.000001),0.9999999),eventweight); 
                }

              
//                 if (Option >= 10 && Option < 20) {
//                   cout << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << mu->pt << " " << mu->eta << " " << mu->phi << " : " 
// //                        << mu->tkNchi2 << " "
// //                        << mu->muNchi2 << " "
// //                        << mu->nValidHits << " "
// //                        << mu->nTkHits << " "
// //                        << mu->nPixHits << " "
// //                        << mu->nMatch << " "
// //                        << mu->d0 << " "
// //                        << mu->ip3d << " "
// //                        << mu->ip3dSig << " "
// //                        << mu->TrkKink << " "
// //                        << mu->SegmentCompatibility << " "
// //                        << mu->CaloCompatilibity << " "
// //                        << mu->HadEnergy/mu->pt << " "
// //                        << mu->HoEnergy/mu->pt << " "
// //                        << mu->EmEnergy/mu->pt << " "
// //                        << mu->HadS9Energy/mu->pt << " "
// //                        << mu->HoS9Energy/mu->pt << " "
// //                        << mu->EmS9Energy/mu->pt << " "                    
//                        << muonIDBDTGV3->MVAValue(mu->pt, mu->eta,
//                                                  mu->tkNchi2,
//                                                  mu->muNchi2,
//                                                  mu->nValidHits,
//                                                  mu->nTkHits,
//                                                  mu->nPixHits,
//                                                  mu->nMatch,
//                                                  mu->d0,
//                                                  mu->ip3d,
//                                                  mu->ip3dSig,
//                                                  mu->TrkKink,
//                                                  mu->SegmentCompatibility,
//                                                  mu->CaloCompatilibity,
//                                                  mu->HadEnergy/mu->pt,
//                                                  mu->HoEnergy/mu->pt,
//                                                  mu->EmEnergy/mu->pt,
//                                                  mu->HadS9Energy/mu->pt,
//                                                  mu->HoS9Energy/mu->pt,
//                                                  mu->EmS9Energy/mu->pt
//                          )
//                        << endl;
//                 }

                if (
                  (Option == 0 && passMuonNumeratorCuts(mu, rho))
//                   (Option == 1 && passMuonIDWithEACorrPFIso(mu, rho ))
                  ||
                  (Option >= 10 && Option < 20 && passMuonMVAIDIsoCombined(mu, muonIDBDTGIDIsoCombined->MVAValue(
                                                                               mu->pt, mu->eta,
                                                                               mu->tkNchi2,
                                                                               mu->muNchi2,
                                                                               mu->nValidHits,
                                                                               mu->nTkHits,
                                                                               mu->nPixHits,
                                                                               mu->nMatch,
                                                                               mu->d0,
                                                                               mu->ip3d,
                                                                               mu->ip3dSig,
                                                                               mu->TrkKink,
                                                                               mu->SegmentCompatibility,
                                                                               mu->CaloCompatilibity,
                                                                               (mu->HadEnergy - rho*MuonEffectiveArea(kMuHadEnergy,mu->eta))/mu->pt,
                                                                               (mu->HoEnergy - rho*MuonEffectiveArea(kMuHoEnergy,mu->eta))/mu->pt,
                                                                               (mu->EmEnergy - rho*MuonEffectiveArea(kMuEmEnergy,mu->eta))/mu->pt,
                                                                               (mu->HadS9Energy - rho*MuonEffectiveArea(kMuHadS9Energy,mu->eta))/mu->pt,
                                                                               (mu->HoS9Energy - rho*MuonEffectiveArea(kMuHoS9Energy,mu->eta))/mu->pt,
                                                                               (mu->EmS9Energy - rho*MuonEffectiveArea(kMuEmS9Energy,mu->eta))/mu->pt,
                                                                               (mu->ChargedIso03 - rho*MuonEffectiveArea(kMuChargedIso03,mu->eta))/mu->pt,
                                                                               (mu->NeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03,mu->eta))/mu->pt,
                                                                               (mu->ChargedIso04 - rho*MuonEffectiveArea(kMuChargedIso04,mu->eta))/mu->pt,
                                                                               (mu->NeutralIso04_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso04,mu->eta))/mu->pt,
                                                                               kFALSE
                                                                               ), rho)
                    )
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, eventweight);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta), eventweight);
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi, eventweight);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, fabs(mu->eta) , eventweight);   
     
                  if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                  if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->nPV0,eventweight);
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(rho,eventweight);
                  }
                  
                }
              }

            } //loop over denominator types
          } //loop over sample types
        } //loop over ptThresholds

      } //loop over muons

    } //end loop over data     
  } //end loop over files
  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_PtEta = 
          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D, ErrorType, kFALSE);
        

        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_nvtx = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_rho = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorM"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);


        TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_nvtx, efficiency_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        
        file->WriteTObject(efficiency_Pt10To20_Barrel_nvtx, efficiency_Pt10To20_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Barrel_rho, efficiency_Pt10To20_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_nvtx, efficiency_Pt10To20_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_rho, efficiency_Pt10To20_Endcap_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_nvtx, efficiency_Pt20ToInf_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_rho, efficiency_Pt20ToInf_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_nvtx, efficiency_Pt20ToInf_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_rho, efficiency_Pt20ToInf_Endcap_rho->GetName(), "WriteDelete");

        file->Close();
        
        if ((denominatorType[denominatorTypeIndex] == 1 || denominatorType[denominatorTypeIndex] == 2)  && sampleLabel[sampleTypeIndex] == "Mu8PtCombinedSample" ) {
          
          TH2F *eff = 0;
          TH2F *effErrorLow = 0;
          TH2F *effErrorHigh = 0;
          
          file = new TFile(smurfOutputFilename.c_str(), "UPDATE");

          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          "MuonFakeRate_M"+IntToString(denominatorType[denominatorTypeIndex])+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D, ErrorType, file);
        
          file->Close();
          delete file;
        }

      }
    }
  }
    
  gBenchmark->Show("WWTemplate");       
} 


// //SMURFV5
// Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu) {
  
//   Bool_t pass = kTRUE;

//   if (mu->pt < 10) pass = kFALSE;
//   if (fabs(mu->eta) > 2.4) pass = kFALSE;

//   Double_t iso = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
//   Double_t isoCutValue = 0;

//   if (fabs(mu->eta) < 1.479) {
//     if (mu->pt > 20) {
//       isoCutValue = 0.13;
//     } else {
//       isoCutValue = 0.06;
//     }
//   } else {
//     if (mu->pt > 20) {
//       isoCutValue = 0.09;
//     } else {
//       isoCutValue = 0.05;
//     }
//   } 

//   if (! 
//       ( mu->typeBits & kGlobal
//         && mu->typeBits & kTracker
//         && mu->nTkHits > 10
//         && mu->muNchi2 < 10.0
//         && (mu->qualityBits & kGlobalMuonPromptTight)
//         && fabs(mu->d0) < 0.02
//         && fabs(mu->dz) < 0.2
//         && iso / mu->pt < isoCutValue
//         && (mu->nMatch > 1 )
//         && (mu->nPixHits > 0)
//         && (mu->pterr / mu->pt < 0.1)
//         )
//     ) pass = kFALSE;

//   if (mu->pt < 20) {
//     if (!
//         ( fabs(mu->d0) < 0.01
//           )
//       ) {
//       pass = kFALSE;
//     }    
//   }

//   return pass;
// }

Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) > 2.4) pass = kFALSE;

  Double_t iso = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
  Double_t isoCutValue = 0;

  if (fabs(mu->eta) < 1.479) {
    if (mu->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (mu->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  } 

  if (! 
      ( (
          (Bool_t(mu->typeBits & kGlobal) 
           && mu->muNchi2 < 10.0
           && (mu->nValidHits > 0)
           && (mu->nMatch > 1 )
            )
          || 
          ( mu->typeBits & kTracker            
            && Bool_t(mu->qualityBits & kTMLastStationTight) 
            )
        )
                
        && mu->nTkHits > 10
        && (mu->nPixHits > 0)
        && fabs(mu->d0) < 0.02
        && fabs(mu->dz) < 0.1
        && iso / mu->pt < isoCutValue
        && (mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,mu->eta)) / mu->pt < 1.0          
        && (mu->pterr / mu->pt < 0.1)
        && (mu->TrkKink < 20)

        )
    ) pass = kFALSE;


  if (mu->pt < 20) {
    if (!
        ( fabs(mu->d0) < 0.01
          )
      ) {
      pass = kFALSE;
    }    
  }

  return pass;
}


Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Int_t DenominatorType, string SampleType, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  //match to HLT
  if (SampleType == "Mu8Sample") {
    if (!( (triggerBits & kHLT_Mu8)
            && (mu->hltMatchBits & kHLTObject_Mu8)
             
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu15Sample") {
    if (!( 
           (triggerBits & kHLT_Mu15)
            && (mu->hltMatchBits & kHLTObject_Mu15)
                       
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu8PtCombinedSample") {
    if (!( 
          ((triggerBits & kHLT_Mu15)
            && (mu->hltMatchBits & kHLTObject_Mu15)
            )
          ||
          ((triggerBits & kHLT_Mu8)
            && (mu->hltMatchBits & kHLTObject_Mu8)
            )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu8Jet40Sample") {
    if (!( (triggerBits & kHLT_Mu8_Jet40)
           && (mu->hltMatchBits & kHLTObject_Mu8)
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "PhotonJetsSample") {
    if (!(
          (triggerBits & kHLT_Mu8_Photon20_CaloIdVT_IsoT )
          && (mu->hltMatchBits & kHLTObject_Mu8)
          )
      ) pass = kFALSE;
  }

  if (!passMuonDenominatorCuts(mu,DenominatorType,fRho)) pass = kFALSE;


  return pass;
}



Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  if (DenominatorType == 1) {

    if (! 
        ( (
            (Bool_t(mu->typeBits & kGlobal) 
             && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              )
            || 
            ( mu->typeBits & kTracker            
              && Bool_t(mu->qualityBits & kTMLastStationTight) 
              )
          ) 
          && mu->typeBits & kTracker
          && mu->nTkHits > 10                   
          && (mu->nPixHits > 0)
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && (mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,mu->eta)) / mu->pt < 1.0          
          && (mu->pterr / mu->pt < 0.1)
          && (mu->TrkKink < 20)
          )
      ) pass = kFALSE;    
  }

  if (DenominatorType == 2) {
    if (! 
        ( (
            (Bool_t(mu->typeBits & kGlobal) 
             && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              )
            || 
            ( mu->typeBits & kTracker            
              && Bool_t(mu->qualityBits & kTMLastStationTight) 
              )
          )
          && mu->typeBits & kTracker
          && mu->nTkHits > 10
          && ( mu->nPixHits > 0)          
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && (mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,mu->eta)) / mu->pt < 0.4          
          && ( mu->pterr / mu->pt < 0.1)
          && (mu->TrkKink < 20)
          )
      ) pass = kFALSE;    
  }

  if (DenominatorType == 3) {
    if (! 
        ( (
            (Bool_t(mu->typeBits & kGlobal) 
             && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              )
            || 
            ( mu->typeBits & kTracker            
              && Bool_t(mu->qualityBits & kTMLastStationTight) 
              )
          )
          && mu->typeBits & kTracker
          && mu->nTkHits > 10
          && ( mu->nPixHits > 0)
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && mu->trkIso03 / mu->pt < 0.3
          && mu->emIso03 / mu->pt < 0.3
          && mu->hadIso03 / mu->pt < 0.3
          && ( mu->pterr / mu->pt < 0.1)
          && (mu->TrkKink < 20)
         )
      ) pass = kFALSE;    
  }


  if (DenominatorType == 100) {
    pass = kTRUE;
  }


  return pass;
}



//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs <<   runNum << " " ;
  ofs <<  lumiSec << " ";
  ofs << evtNum<< " ";
  ofs << mass<< " ";

//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
