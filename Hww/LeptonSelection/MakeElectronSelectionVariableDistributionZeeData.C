//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\"/home/sixie/hist/WWAnalysis/HwwNtuple_f10-h160ww2l-gf-z2-v12_noskim_0000.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_p10-w1jets-0-100-v26_noskim.old.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\"HwwNtuple_w10-h130ww2l-gf-z2-v8-pu11_noskim_0000.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\"HwwNtuple_w10-wjetsl-z2-v8-pu11_noskim_0000.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysisSkimmed_mu_noskim.root\",\"\",kFALSE\)
//root -l EWKAna/Hww/LeptonSelection/ElectronSelection.C+\(\)
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
#include <TLegend.h>               // 3D vector class
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
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData);
Bool_t passTightElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Double_t pt2, Double_t eta2, Double_t phi2);
void runMakeElectronPlots(const string dataInputFilename,       
                          const string Label, Bool_t isMuonData);
void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname);

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

void MakeElectronSelectionVariableDistributionZeeData() 
{  
  
//    runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_ele-d22_noskim.root","Data_lowPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/WWAnalysisSkimmed_ele-d22_noskim.root","Data_highPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/HwwAnalysis_w10-zee-powheg-c10-v8-pu11_noskim.root","MC_lowPt",kFALSE);
//   runMakeElectronPlots("/home/sixie/hist/HwwAnalysis/HwwAnalysis_w10-zee-powheg-c10-v8-pu11_noskim.root","MCTagAndProbePU11_highPt",kFALSE);

//      runMakeElectronPlots("/home/sixie/hist/WWAnalysis/WWAnalysis_f10-zee-powheg-c10-v12-pu_noskim.root","MCTagAndProbe_lowPt",kFALSE);
  runMakeElectronPlots("/home/sixie/hist/WWAnalysis/WWAnalysis_f10-zee-powheg-c10-v12-pu_noskim.root","MCTagAndProbe_highPt",kFALSE);



}



void runMakeElectronPlots(const string dataInputFilename,       
                const string Label, Bool_t isMuonData) 
{  
  gBenchmark->Start("WWTemplate");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  Double_t lumi;              // luminosity (pb^-1)
    

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1D *dileptonMass = new TH1D("dileptonMass", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_ee = new TH1D("dileptonMass_ee", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_emu = new TH1D("dileptonMass_emu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  TH1D *dileptonMass_mumu = new TH1D("dileptonMass_mumu", "; Mass [GeV/c^{2}]; Number of Events", 50, 0, 200);
  Int_t Count_ee_BB = 0;
  Int_t Count_ee_BE = 0;
  Int_t Count_ee_EE = 0;
  Int_t Count_mm_BB = 0;
  Int_t Count_mm_BE = 0;
  Int_t Count_mm_EE = 0;
  Int_t Count_em_BB = 0;
  Int_t Count_em_BE = 0;
  Int_t Count_em_EE = 0;
  ofstream eventListFile("eventList.txt");

  TH1D *ptMin = new TH1D("ptMin", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  TH1D *ptMax = new TH1D("ptMax", "; p_{T} [GeV]; Number of Events", 100, 0, 100);
  string tmplabel = Label;
  if (tmplabel != "") tmplabel = "_Zee_" + Label;
  //electron variables
  TH1F *sigmaIEtaIEta_Barrel = new TH1F((string("sigmaIEtaIEta_Barrel")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.05);
  TH1F *DeltaEta_Barrel = new TH1F((string("DeltaEta_Barrel")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Barrel = new TH1F((string("DeltaPhi_Barrel")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.03, 0.03);
  TH1F *HOverE_Barrel = new TH1F((string("HOverE_Barrel")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *relIso_Barrel = new TH1F((string("relIso_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Barrel = new TH1F((string("fBrem_Barrel")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);

  TH1F *sigmaIEtaIEta_Endcap = new TH1F((string("sigmaIEtaIEta_Endcap")+tmplabel).c_str(), "; sigma ieta ieta; Fraction of Events ", 100, 0, 0.1);
  TH1F *DeltaEta_Endcap = new TH1F((string("DeltaEta_Endcap")+tmplabel).c_str(), "; deltaEta; Fraction of Events ", 100, -0.02, 0.02);
  TH1F *DeltaPhi_Endcap = new TH1F((string("DeltaPhi_Endcap")+tmplabel).c_str(), "; deltaPhi; Fraction of Events ", 100, -0.1, 0.1);
  TH1F *HOverE_Endcap = new TH1F((string("HOverE_Endcap")+tmplabel).c_str(), "; HOverE; Fraction of Events ", 100, 0, 0.2);
  TH1F *relIso_Endcap = new TH1F((string("relIso_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);
  TH1F *fBrem_Endcap = new TH1F((string("fBrem_Endcap")+tmplabel).c_str(), "; relIso; Fraction of Events ", 100, 0, 1.0);

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
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  

  infile = new TFile(dataInputFilename.c_str());
  assert(infile);

    
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kFALSE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("merged_JsonReRecoSep17_JsonStreamExpressV2.txt"); 
  

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(), "Events"); assert(eventTree);
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr); TBranch *muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("PFJet", &jetArr); TBranch *jetBr = eventTree->GetBranch("PFJet");
  

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
	
//     mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
//     if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
//     if (!passHLT(info->triggerBits, info->runNum, isMuonData)) continue;


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
    // TcMet
    //********************************************************

    Int_t NLeptons = 0;
    vector<Int_t> leptonType;
    vector<Int_t> leptonIndex;
    vector<Double_t> leptonPt;
    vector<Double_t> leptonEta;
    vector<Double_t> leptonPhi;

    Int_t NJets = 0;
    const mithep::TJet *leadingJet = 0;
    
     //******************************************************************************
    //loop over electron pairs
    //******************************************************************************
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele1 = (mithep::TElectron*)((*electronArr)[i]);
      if ( !(
             fabs(ele1->eta) < 2.5
             &&
             ele1->pt > 20.0
             &&
             passTightElectronCuts(ele1)
//              &&
//              ele1->isUsed == kTRUE
             )
        ) continue;
      

      for(Int_t j=0; j<electronArr->GetEntries(); j++) {
        if (i==j) continue;

        const mithep::TElectron *ele2 = (mithep::TElectron*)((*electronArr)[j]);
        if ( !(
               fabs(ele2->eta) < 2.5
               && ele2->pt > 10.0
               )
          ) continue;
        
        mithep::FourVectorM lepton1;
        mithep::FourVectorM lepton2;
        lepton1.SetCoordinates(ele1->pt, ele1->eta, ele1->phi, 0.51099892e-3 );
        lepton2.SetCoordinates(ele2->pt, ele2->eta, ele2->phi, 0.51099892e-3 );
        mithep::FourVectorM dilepton = lepton1+lepton2;
 
        ptMin->Fill(ele2->pt);
        if (ele2->pt >= 20 && ele2->pt < 30 ){ 
          if (fabs(ele2->eta) < 1.5) {
            if ((ele2->trkIso03 + TMath::Max(ele2->emIso03 - 1.0, 0.0) + ele2->hadIso03) / ele2->pt < 0.1) {
              //dileptonMass->Fill(dilepton.M());
            }
            if (ele2->sigiEtaiEta < 0.012) {
              dileptonMass->Fill(dilepton.M()); 
            }
          } else {
            if ((ele2->trkIso03 + ele2->emIso03  + ele2->hadIso03) / ele2->pt < 0.1) {
              //dileptonMass->Fill(dilepton.M());
            }
            if (ele2->sigiEtaiEta < 0.035) {
              dileptonMass->Fill(dilepton.M());
            }
          }          
        }

        if (ele2->pt >= 20 && ele2->pt < 30
            &&             
            dilepton.M() > 61 && dilepton.M() < 121
          ) {
          if (fabs(ele2->eta) < 1.5) {
            if ((ele2->trkIso03 + TMath::Max(ele2->emIso03 - 1.0, 0.0) + ele2->hadIso03) / ele2->pt < 0.1) {
              sigmaIEtaIEta_Barrel->Fill( ele2->sigiEtaiEta );
              DeltaEta_Barrel->Fill( ele2->deltaEtaIn );
              DeltaPhi_Barrel->Fill( ele2->deltaPhiIn );
              HOverE_Barrel->Fill( ele2->HoverE );
            }
            if (ele2->sigiEtaiEta < 0.012) {
              relIso_Barrel->Fill( (ele2->trkIso03 + TMath::Max(ele2->emIso03 - 1.0, 0.0) + ele2->hadIso03) / ele2->pt );
            }
          } else {
            if ((ele2->trkIso03 + ele2->emIso03  + ele2->hadIso03) / ele2->pt < 0.1) {
              sigmaIEtaIEta_Endcap->Fill( ele2->sigiEtaiEta );
              DeltaEta_Endcap->Fill( ele2->deltaEtaIn );
              DeltaPhi_Endcap->Fill( ele2->deltaPhiIn );
              HOverE_Endcap->Fill( ele2->HoverE );
            }
            if (ele2->sigiEtaiEta < 0.035) {
              relIso_Endcap->Fill( (ele2->trkIso03 + ele2->emIso03  + ele2->hadIso03) / ele2->pt );
            }
          }
        }
      }
    }


  } //end loop over data     

  infile->Close(); delete infile;
  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("WWSelectionPlotsFakePrediction.root", "UPDATE");


  file->WriteTObject(sigmaIEtaIEta_Barrel,sigmaIEtaIEta_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(DeltaEta_Barrel,DeltaEta_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(DeltaPhi_Barrel,DeltaPhi_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(HOverE_Barrel,HOverE_Barrel->GetName(), "WriteDelete");
  file->WriteTObject(relIso_Barrel,relIso_Barrel->GetName(), "WriteDelete");

  file->WriteTObject(sigmaIEtaIEta_Endcap,sigmaIEtaIEta_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(DeltaEta_Endcap,DeltaEta_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(DeltaPhi_Endcap,DeltaPhi_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(HOverE_Endcap,HOverE_Endcap->GetName(), "WriteDelete");
  file->WriteTObject(relIso_Endcap,relIso_Endcap->GetName(), "WriteDelete");
  file->Close();

//   //--------------------------------------------------------------------------------------------------------------
//   // Make plots
//   //==============================================================================================================
   TCanvas *cv = new TCanvas("cv","cv", 800,600);
   ptMin->Draw();
   cv->SaveAs("ptMin.gif");

   dileptonMass->Draw();
   cv->SaveAs("DileptonMass.gif");

//   ptMin->Draw();
//   cv->SaveAs("ptMin.gif");
//   ptMax->Draw();
//   cv->SaveAs("ptMax.gif");

//   sigmaIEtaIEta_Barrel_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Barrel_signal.gif");
//   DeltaEta_Barrel_signal->Draw();
//   cv->SaveAs("DeltaEta_Barrel_signal.gif");
//   DeltaPhi_Barrel_signal->Draw();
//   cv->SaveAs("DeltaPhi_Barrel_signal.gif");
//   HOverE_Barrel_signal->Draw();
//   cv->SaveAs("HOverE_Barrel_signal.gif");
//   relIso_Barrel_signal->Draw();
//   cv->SaveAs("relIso_Barrel_signal.gif");

//   sigmaIEtaIEta_Endcap_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Endcap_signal.gif");
//   DeltaEta_Endcap_signal->Draw();
//   cv->SaveAs("DeltaEta_Endcap_signal.gif");
//   DeltaPhi_Endcap_signal->Draw();
//   cv->SaveAs("DeltaPhi_Endcap_signal.gif");
//   HOverE_Endcap_signal->Draw();
//   cv->SaveAs("HOverE_Endcap_signal.gif");
//   relIso_Endcap_signal->Draw();
//   cv->SaveAs("relIso_Endcap_signal.gif");


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 


        
  gBenchmark->Show("plotZ");       
} 

void makePlot(string signalHistName, string signalLabel, string bkgHistName, string bkgLabel, string plotname) {  
  TFile *file = new TFile("WWSelectionPlotsFakePrediction.root", "READ");

  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TH1F *signal = 0;
  TH1F *bkg = 0;
  TLegend *legend = new TLegend(0.73,0.75,0.93,0.90);
  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);
    
 
  signal = (TH1F*)file->Get(signalHistName.c_str());
  bkg = (TH1F*)file->Get(bkgHistName.c_str());
  if (!signal) {cout << "Cannot load histogram " << signalHistName << endl;  assert(signal);}
  if (!bkg) {cout << "Cannot load histogram " << bkgHistName << endl;  assert(bkg);}

  Double_t norm = 0;
  norm = signal->Integral();
  for (int b=0; b<signal->GetXaxis()->GetNbins()+2; ++b) {
    signal->SetBinContent(b,signal->GetBinContent(b) / norm);
    signal->SetBinError(b,signal->GetBinError(b) / norm);
  }
  norm = bkg->Integral();
  for (int b=0; b<bkg->GetXaxis()->GetNbins()+2; ++b) {
    bkg->SetBinContent(b,bkg->GetBinContent(b) / norm);
    bkg->SetBinError(b,bkg->GetBinError(b) / norm);
  }

  Double_t maxY = signal->GetMaximum();
  if (bkg->GetMaximum() > maxY) maxY = bkg->GetMaximum();
  maxY = maxY * 1.2;
  signal->SetMaximum(maxY);
  bkg->SetMaximum(maxY);

  signal->SetLineColor(kRed);
  signal->SetMarkerColor(kRed);
  bkg->SetLineColor(kBlue);
  bkg->SetMarkerColor(kBlue);
  legend->Clear();
  legend->AddEntry(signal, signalLabel.c_str(),  "LP");
  legend->AddEntry(bkg,    bkgLabel.c_str(), "LP");

  signal->Draw();
  bkg->Draw("same");
  legend->Draw();
  cv->SaveAs((plotname+".gif").c_str());

//   sigmaIEtaIEta_Barrel_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Barrel_signal.gif");
//   DeltaEta_Barrel_signal->Draw();
//   cv->SaveAs("DeltaEta_Barrel_signal.gif");
//   DeltaPhi_Barrel_signal->Draw();
//   cv->SaveAs("DeltaPhi_Barrel_signal.gif");
//   HOverE_Barrel_signal->Draw();
//   cv->SaveAs("HOverE_Barrel_signal.gif");
//   relIso_Barrel_signal->Draw();
//   cv->SaveAs("relIso_Barrel_signal.gif");

//   sigmaIEtaIEta_Endcap_signal->Draw();
//   cv->SaveAs("sigmaIEtaIEta_Endcap_signal.gif");
//   DeltaEta_Endcap_signal->Draw();
//   cv->SaveAs("DeltaEta_Endcap_signal.gif");
//   DeltaPhi_Endcap_signal->Draw();
//   cv->SaveAs("DeltaPhi_Endcap_signal.gif");
//   HOverE_Endcap_signal->Draw();
//   cv->SaveAs("HOverE_Endcap_signal.gif");
//   relIso_Endcap_signal->Draw();
//   cv->SaveAs("relIso_Endcap_signal.gif");

  file->Close();
}


Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData) {


  Bool_t pass = kFALSE;

//   if (isMuonData) {
//     if ((runNum >= 132440) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
//     } 
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140401)) {
// //       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( (triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   } else {
//     //it's electron data
//     if ((runNum >= 132440) && (runNum <= 135058)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Photon10_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 135059) && (runNum <= 140041)) {
// //       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_LW_L1R) ) pass = kTRUE;
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 140042) && (runNum <= 141900)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//     } 
//     if ((runNum >= 141901) && (runNum <= 146427)) {
//       if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 146428) && (runNum <= 147119)) {
//       if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//     }
//     if ((runNum >= 147120) && (runNum <= 999999)) {
//       if ( !(triggerBits & kHLT_Mu15_v1) && (triggerBits & kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1) ) pass = kTRUE;
//     }
//   }

  return pass;

}



Bool_t passTightElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.003
            && fabs(ele->deltaPhiIn) < 0.01
            && ele->HoverE < 0.025
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.01
            && ele->nExpHitsInner <= 0
            && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
              && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.005
             && fabs(ele->deltaPhiIn) < 0.02
             && ele->HoverE < 0.01
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.01
             && ele->nExpHitsInner <= 0
             && !(fabs(ele->partnerDist) < 0.02 && fabs(ele->partnerDeltaCot) < 0.02)
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


Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //ECAL driven only
  if (!ele->isEcalDriven) {
    pass = kFALSE;
    return pass;
  }
  
  //Barrel 
  if (fabs(ele->eta) < 1.5) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.014          
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.5
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else if (fabs(ele->eta) > 1.5) {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.034
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.5
          )
      ) {
      pass = kFALSE;
    }
  }

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

        && (mu->nSeg >= 2)
        && (mu->nValidHits >= 1)
        && (mu->nPixHits >= 1)
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
