//root -l EWKAna/Hww/FakeRate/ComputeMuonFakeRate_V4.C+\(\)
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


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType);
Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DoComputeMuonFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               Int_t SelectionType);

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
void ComputeMuonFakeRate_V4() {


   DoComputeMuonFakeRate("/home/sixie/hist/FakeRate/FakeRate_r10a-mu-d22_fakeskim.root","MuonFakeRate_v4_Mu9Jet30","MuonFakeRate_V4.root",2);

}



void DoComputeMuonFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               Int_t SelectionType) 
{  
  gBenchmark->Start("WWTemplate");


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  Double_t jetPtThreshold = 0;
  if (SelectionType == 1) {
    jetPtThreshold = 15.0;
  } else if (SelectionType == 2) {
    jetPtThreshold = 30.0;
  }else if (SelectionType == 3) {
    jetPtThreshold = 50.0;
  }else if (SelectionType == 4) {
    jetPtThreshold = 70.0;
  }else if (SelectionType == 5) {
    jetPtThreshold = 100.0;
  }


  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
  ptbins.push_back(10);  
  ptbins.push_back(12.5);  
  ptbins.push_back(15);  
  ptbins.push_back(20);  
  ptbins.push_back(25);  
  ptbins.push_back(30);  
  ptbins.push_back(35);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
  ptbins.push_back(80);  


  vector<double> etabins;
  etabins.push_back(-2.75);
  etabins.push_back(-2.25);
  etabins.push_back(-1.75);
  etabins.push_back(-1.25);
  etabins.push_back(-0.75);
  etabins.push_back(-0.25);
  etabins.push_back(0.25);
  etabins.push_back(0.75);
  etabins.push_back(1.25);
  etabins.push_back(1.75);
  etabins.push_back(2.25);
  etabins.push_back(2.75);

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


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *leadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  TH1F *denominator_Pt = new TH1F("denominator_Pt" , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
  TH1F *denominator_Eta = new TH1F("denominator_Eta" , "; #eta ; Number of Events ",  100, -3.0 , 3.0);
  TH1F *denominator_Phi = new TH1F("denominator_Phi" , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
  TH2F *denominator_PtEta = new TH2F("denominator_PtEta" , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, -3.5, 3.5);

  TH1F *numerator_Pt = new TH1F("numerator_Pt" , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
  TH1F *numerator_Eta = new TH1F("numerator_Eta" , "; #eta ; Number of Events ",  100, -3.0 , 3.0);
  TH1F *numerator_Phi = new TH1F("numerator_Phi" , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
  TH2F *numerator_PtEta = new TH2F("numerator_PtEta" , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, -3.5, 3.5);


  //N-1 histograms
  string tmplabel = label; if (tmplabel != "") tmplabel = "_"+label;

  TH1F *sigmaIEtaIEta_Barrel = new TH1F((string("sigmaIEtaIEta_Barrel")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 25, 0, 0.05);
  TH1F *DeltaEta_Barrel = new TH1F((string("DeltaEta_Barrel")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 25, -0.02, 0.02);
  TH1F *DeltaPhi_Barrel = new TH1F((string("DeltaPhi_Barrel")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 25, -0.5, 0.5);
  TH1F *HOverE_Barrel = new TH1F((string("HOverE_Barrel")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 25, 0, 2.0);
  TH1F *relIso_Barrel = new TH1F((string("relIso_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 25, 0, 1.0);
  TH1F *NExpectedHits_Barrel = new TH1F((string("NExpectedHits_Barrel")+tmplabel).c_str(), "; NExpectedHits; Fraction of Events ", 25, 0, 1.0);
  TH1F *d0_Barrel = new TH1F((string("d0_Barrel")+tmplabel).c_str(), "; d0; Fraction of Events ", 25, 0, 0.05);
  TH1F *sigmaIEtaIEta_Endcap = new TH1F((string("sigmaIEtaIEta_Endcap")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 25, 0, 0.1);
  TH1F *DeltaEta_Endcap = new TH1F((string("DeltaEta_Endcap")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 25, -0.02, 0.02);
  TH1F *DeltaPhi_Endcap = new TH1F((string("DeltaPhi_Endcap")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 25, -0.5, 0.5);
  TH1F *HOverE_Endcap = new TH1F((string("HOverE_Endcap")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 25, 0, 2.0);
  TH1F *relIso_Endcap = new TH1F((string("relIso_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 25, 0, 1.0);
  TH1F *NExpectedHits_Endcap = new TH1F((string("NExpectedHits_Endcap")+tmplabel).c_str(), "; NExpectedHits; Fraction of Events ", 25, 0, 1.0);
  TH1F *d0_Endcap = new TH1F((string("d0_Endcap")+tmplabel).c_str(), "; d0; Fraction of Events ", 25, 0, 0.05);

  TH1F *sigmaIEtaIEta_Barrel_NMinusOne = new TH1F((string("sigmaIEtaIEta_Barrel_NMinusOne")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 25, 0, 0.05);
  TH1F *DeltaEta_Barrel_NMinusOne = new TH1F((string("DeltaEta_Barrel_NMinusOne")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 25, -0.02, 0.02);
  TH1F *DeltaPhi_Barrel_NMinusOne = new TH1F((string("DeltaPhi_Barrel_NMinusOne")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 25, -0.5, 0.5);
  TH1F *HOverE_Barrel_NMinusOne = new TH1F((string("HOverE_Barrel_NMinusOne")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 25, 0, 2.0);
  TH1F *relIso_Barrel_NMinusOne = new TH1F((string("relIso_Barrel_NMinusOne")+tmplabel).c_str(), "; relIso; Fraction of Events ", 25, 0, 1.0);
  TH1F *NExpectedHits_Barrel_NMinusOne = new TH1F((string("NExpectedHits_Barrel_NMinusOne")+tmplabel).c_str(), "; NExpectedHits; Fraction of Events ", 25, 0, 1.0);
  TH1F *d0_Barrel_NMinusOne = new TH1F((string("d0_Barrel_NMinusOne")+tmplabel).c_str(), "; d0; Fraction of Events ", 25, 0, 0.05);
  TH1F *sigmaIEtaIEta_Endcap_NMinusOne = new TH1F((string("sigmaIEtaIEta_Endcap_NMinusOne")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 25, 0, 0.1);
  TH1F *DeltaEta_Endcap_NMinusOne = new TH1F((string("DeltaEta_Endcap_NMinusOne")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 25, -0.02, 0.02);
  TH1F *DeltaPhi_Endcap_NMinusOne = new TH1F((string("DeltaPhi_Endcap_NMinusOne")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 25, -0.5, 0.5);
  TH1F *HOverE_Endcap_NMinusOne = new TH1F((string("HOverE_Endcap_NMinusOne")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 25, 0, 2.0);
  TH1F *relIso_Endcap_NMinusOne = new TH1F((string("relIso_Endcap_NMinusOne")+tmplabel).c_str(), "; relIso; Fraction of Events ", 25, 0, 1.0);
  TH1F *NExpectedHits_Endcap_NMinusOne = new TH1F((string("NExpectedHits_Endcap_NMinusOne")+tmplabel).c_str(), "; NExpectedHits; Fraction of Events ", 25, 0, 1.0);
  TH1F *d0_Endcap_NMinusOne = new TH1F((string("d0_Endcap_NMinusOne")+tmplabel).c_str(), "; d0; Fraction of Events ", 25, 0, 0.05);



  ofstream eventListFile("eventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
   rlrm.AddJSONFile("merged_JsonReRecoSep17_JsonStreamExpressV2.txt"); 
   hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

  cout << "Total Events: " << eventTree->GetEntries() << endl;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
		
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear(); 
    jetArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    jetBr->GetEntry(ientry);


    //********************************************************
    // TcMet
    //********************************************************
    TVector3 met;        
    if(info->tcMEx!=0 || info->tcMEy!=0) {       
      met.SetXYZ(info->tcMEx, info->tcMEy, 0);
    }
	
    //********************************************************
    // Event Selection Cuts
    //********************************************************
    if (met.Pt() > 20) continue;
    if (muonArr->GetEntries() > 1) continue;

    Double_t tempLeadingJetPt = 0;
    for(Int_t i=0; i<jetArr->GetEntries(); i++) {
      const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
      if( jet->pt > tempLeadingJetPt) tempLeadingJetPt = jet->pt;
    }
    leadingJetPt->Fill(tempLeadingJetPt);
    
 

    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);

       //pass HLT selection
      if (!passHLT(info->triggerBits, info->runNum, SelectionType)) continue;
      
      //pass event selection
      Bool_t passJetSelection = kFALSE;
      for(Int_t i=0; i<jetArr->GetEntries(); i++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
        
        if (jet->pt > jetPtThreshold &&
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 0.3) {
          passJetSelection = kTRUE;
          break;
        }
      }
      if (!passJetSelection) continue;

      if (!passMuonDenominatorCuts(mu)) continue;

      denominator_Pt->Fill(mu->pt);
      denominator_Eta->Fill(mu->eta);
      denominator_Phi->Fill(mu->phi);
      denominator_PtEta->Fill(mu->pt, mu->eta);
      if (passMuonNumeratorCuts(mu)) {
        numerator_Pt->Fill(mu->pt);
        numerator_Eta->Fill(mu->eta);
        numerator_Phi->Fill(mu->phi);
        numerator_PtEta->Fill(mu->pt, mu->eta );   
      }           
    }

  } //end loop over data     


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;



  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  Int_t ErrorType = 2; //Clopper Pearson errors
  TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(numerator_Pt, denominator_Pt, label+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(numerator_Eta, denominator_Eta, label+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
  TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(numerator_Phi, denominator_Phi, label+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
  mithep::TH2DAsymErr *efficiency_PtEta = 
    mithep::EfficiencyUtils::createEfficiencyHist2D(numerator_PtEta, 
                                                    denominator_PtEta, label+"_PtEta", 
                                                    ptbins, etabins, ErrorType);


  //*****************************************************************************************
  //Draw Plots
  //*****************************************************************************************
  TCanvas *cv = MakeCanvas("cv", "cv", 800, 600);
  denominator_Pt->Draw();
  cv->SaveAs((label+"_DenominatorPt.gif").c_str());
  numerator_Pt->Draw();
  cv->SaveAs((label+"_NumeratorPt.gif").c_str());
  efficiency_pt->Draw("AP");
  cv->SaveAs((label+"_FakeRatePt.gif").c_str());
 

  //*****************************************************************************************
  //Save Efficiency Plots
  //*****************************************************************************************
  TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
  file->cd();

  file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
  file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");
 

  file->Close();
  delete file;


    
  gBenchmark->Show("WWTemplate");       
} 



Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Int_t SelectionType) {


  Bool_t pass = kFALSE;

  //it's electron data
  if (SelectionType < 10) {
    if ((runNum >= 136033) && (runNum <= 147116)) {
      if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
    } 
    if ((runNum >= 147196) && (runNum <= 999999)) {
      if ( (triggerBits & kHLT_Mu15) ) pass = kTRUE;
    }
  } else if (SelectionType >= 100 && SelectionType < 200) {
    if (SelectionType == 101) {
      if ( (triggerBits & kHLT_Jet15U) ) pass = kTRUE;
    } else if (SelectionType == 102) {
      if ( (triggerBits & kHLT_Jet30U) ) pass = kTRUE;
    }else if (SelectionType == 103) {
      if ( (triggerBits & kHLT_Jet50U) ) pass = kTRUE;
    }else if (SelectionType == 104) {
      if ( (triggerBits & kHLT_Jet70U) ) pass = kTRUE;
    }else if (SelectionType == 105) {
      if ( (triggerBits & kHLT_Jet100U) ) pass = kTRUE;
    }

  } else if (SelectionType >= 200 && SelectionType < 300) {
    if (SelectionType == 201) {
      if ( (triggerBits & kHLT_Photon20_L1R) ) pass = kTRUE;
    } else if (SelectionType == 202) {
      if ( (triggerBits & kHLT_Photon30_L1R) ) pass = kTRUE;
    }else if (SelectionType == 202) {
      if ( (triggerBits & kHLT_Photon40_L1R) ) pass = kTRUE;
    }else if (SelectionType == 202) {
      if ( (triggerBits & kHLT_Photon50_L1R) ) pass = kTRUE;
    }

  } else {
    cout << "Selection Type " << SelectionType << " not recognized\n";
  }

  return pass;

}



Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu) {
  
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




Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.2
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 1.0
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

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
