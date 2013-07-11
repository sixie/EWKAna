//root -l EWKAna/Hww/Efficiency/ElectronTagAndProbeMC.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_f10-zee-powheg-c10-v12_noskim.root\",\"\"\)



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
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TGraphAsymmErrors.h"

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
Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum, Bool_t isMC = kTRUE);
Bool_t passElectronTagCuts(const mithep::TElectron *ele);
Bool_t passTighterElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronCuts(const mithep::TElectron *ele);
void WriteMassToFile( ofstream *file, Double_t mass);
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
      cout << "Cannot get Directory ZeeAnalysisMod from file " << infname << endl;
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

void ElectronTagAndProbeMC(const string dataInputFilename, const string Label, 
                           Int_t ChargeSelection = 0) {   

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
  vector<string> ptbinLabel;
  vector<Double_t> ptbinLowEdge;
  vector<Double_t> ptbinUpEdge;
   ptbinLabel.push_back("Pt20ToInf"); ptbinLowEdge.push_back(20.0); ptbinUpEdge.push_back(14000.0);
//   ptbinLabel.push_back("Pt10To15"); ptbinLowEdge.push_back(10.0); ptbinUpEdge.push_back(15.0);
//   ptbinLabel.push_back("Pt15To20"); ptbinLowEdge.push_back(15.0); ptbinUpEdge.push_back(20.0);
//   ptbinLabel.push_back("Pt20To30"); ptbinLowEdge.push_back(20.0); ptbinUpEdge.push_back(30.0);
//   ptbinLabel.push_back("Pt30To40"); ptbinLowEdge.push_back(30.0); ptbinUpEdge.push_back(40.0);
//   ptbinLabel.push_back("Pt40To50"); ptbinLowEdge.push_back(40.0); ptbinUpEdge.push_back(50.0);
//   ptbinLabel.push_back("Pt50ToInf"); ptbinLowEdge.push_back(50.0); ptbinUpEdge.push_back(14000.0);
  vector<string> etabinLabel;
  vector<Double_t> etabinLowEdge;
  vector<Double_t> etabinUpEdge;
  etabinLabel.push_back("EE"); etabinLowEdge.push_back(1.5); etabinUpEdge.push_back(2.5);
  etabinLabel.push_back("EB"); etabinLowEdge.push_back(0); etabinUpEdge.push_back(1.5);
  vector<string> chargebinLabel;
  vector<Double_t> chargeSelection;
  chargebinLabel.push_back("All");  chargeSelection.push_back(0);
  chargebinLabel.push_back("Plus"); chargeSelection.push_back(1);
  chargebinLabel.push_back("Minus"); chargeSelection.push_back(-1);
  vector<vector<vector<Double_t> > > EventCount_TagPlusRecoFailWWTightIdIso_binned;
  vector<vector<vector<Double_t> > > EventCount_TagPlusRecoPassWWTightIdIso_binned;
  vector<vector<vector<Double_t> > > EventCount_TagPlusWWTightIdIsoPassHLT_binned;
  vector<vector<vector<Double_t> > > EventCount_TagPlusWWTightIdIsoFailHLT_binned;
  for (int k=0; k < chargebinLabel.size(); ++k) {
    vector<vector<Double_t> > tmp_tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned;
    vector<vector<Double_t> > tmp_tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned;
    vector<vector<Double_t> > tmp_tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned;
    vector<vector<Double_t> > tmp_tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned;    
    
    for (int i=0; i < etabinLabel.size(); ++i) {
      vector<Double_t> tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned;
      vector<Double_t> tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned;
      vector<Double_t> tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned;
      vector<Double_t> tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned;    
      
      for (int j=0; j < ptbinLabel.size(); ++j) {
        Double_t tmp_EventCount_TagPlusRecoFailWWTightIdIso = 0;
        tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned.push_back(tmp_EventCount_TagPlusRecoFailWWTightIdIso);
        
        Double_t tmp_EventCount_TagPlusRecoPassWWTightIdIso = 0;
        tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned.push_back(tmp_EventCount_TagPlusRecoPassWWTightIdIso);
        
        Double_t tmp_EventCount_TagPlusWWTightIdIsoFailHLT = 0;
        tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned.push_back(tmp_EventCount_TagPlusWWTightIdIsoFailHLT);
        
        Double_t tmp_EventCount_TagPlusWWTightIdIsoPassHLT = 0;
        tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned.push_back(tmp_EventCount_TagPlusWWTightIdIsoPassHLT);
        
        cout << "make bins : " << k << " " << i << " " << j << " " << endl;
        
      }
      
      tmp_tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned.push_back(tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned);
      tmp_tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned.push_back(tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned);
      tmp_tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned.push_back(tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned);
      tmp_tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned.push_back(tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned);
      
    }

    EventCount_TagPlusRecoFailWWTightIdIso_binned.push_back(tmp_tmp_EventCount_TagPlusRecoFailWWTightIdIso_binned);
    EventCount_TagPlusRecoPassWWTightIdIso_binned.push_back(tmp_tmp_EventCount_TagPlusRecoPassWWTightIdIso_binned);
    EventCount_TagPlusWWTightIdIsoPassHLT_binned.push_back(tmp_tmp_EventCount_TagPlusWWTightIdIsoPassHLT_binned);
    EventCount_TagPlusWWTightIdIsoFailHLT_binned.push_back(tmp_tmp_EventCount_TagPlusWWTightIdIsoFailHLT_binned);
    
  }

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  Double_t EventCount_TagPlusRecoFailWWTightIdIso = 0;
  Double_t EventCount_TagPlusRecoPassWWTightIdIso = 0;
  Double_t EventCount_TagPlusWWTightIdIsoFailHLT = 0;
  Double_t EventCount_TagPlusWWTightIdIsoPassHLT = 0;

  Int_t nbins = 20;
  TH1F *histogram = new TH1F ("histogram", "; xaxis; yaxis; ", 100, 0 , 200);
  TH1F *ProbePt = new TH1F ("ProbePt", "; xaxis; yaxis; ", 100, 0 , 200);
  TH1F *TagPt = new TH1F ("TagPt", "; xaxis; yaxis; ", 100, 0 , 200);
  TH1F *TagPlusRecoFailWWTightIdIso = new TH1F ("TagPlusRecoFailWWTightIdIso", "; xaxis; yaxis; ", nbins, 60 , 120);
  TH1F *TagPlusRecoPassWWTightIdIso = new TH1F ("TagPlusRecoPassWWTightIdIso", "; xaxis; yaxis; ", nbins, 60 , 120);
  TH1F *TagPlusWWTightIdIsoFailHLT = new TH1F ("TagPlusWWTightIdIsoFailHLT", "; xaxis; yaxis; ", nbins, 60 , 120);
  TH1F *TagPlusWWTightIdIsoPassHLT = new TH1F ("TagPlusWWTightIdIsoPassHLT", "; xaxis; yaxis; ", nbins, 60 , 120);


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
  

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(dataInputFilename.c_str(),"Events");
  
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);          TBranch *infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  

  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  Double_t nsel=0, nselvar=0;
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
    infoBr->GetEntry(ientry);
	
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear();    
    electronBr->GetEntry(ientry);
 

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
    if (electronArr->GetEntries() < 2) continue;


    //******************************************************************************
    //loop over electron pairs
    //******************************************************************************
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
      if ( !(
             fabs(tag->eta) < 2.5
             && 
             tag->pt > 20.0
             &&
             passElectronTagCuts(tag)
             )
        ) continue;
      

      for(Int_t j=0; j<electronArr->GetEntries(); j++) {
        if (i==j) continue;

        const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
        if ( !(
               fabs(probe->eta) < 2.5
               && probe->pt > 10.0
               && probe->isEcalDriven
               )
          ) continue;
        
        mithep::FourVectorM tagVector;
        mithep::FourVectorM probeVector;
        tagVector.SetCoordinates(tag->pt, tag->eta, tag->phi, 0.51099892e-3 );
        probeVector.SetCoordinates(probe->pt, probe->eta, probe->phi, 0.51099892e-3 );
        mithep::FourVectorM dilepton = tagVector+probeVector;
 
        if (dilepton.M() > 60 && dilepton.M() < 120
            
          ) {
          
          //For binned efficiencies
          for (int q=0; q < chargebinLabel.size(); ++q) {
            for (int e=0; e < etabinLabel.size(); ++e) {
              for (int p=0; p < ptbinLabel.size(); ++p) {
                
                //Require the probe is in the right pt and eta bin
                if ( 		  
                (probe->pt >= ptbinLowEdge[p] && probe->pt < ptbinUpEdge[p])
                &&
                (fabs(probe->eta) >= etabinLowEdge[e] && fabs(probe->eta) < etabinUpEdge[e]) 
                && 
                (probe->q == chargeSelection[q] || chargeSelection[q] == 0)
//                 &&
//                 met.Pt() < 20.0
                  ) {
                  
                  //*****************************************************************************************************
                  //Reco -> WWTight
                  //*****************************************************************************************************	    
                
                  if (passElectronCuts(probe)) {	 
                    EventCount_TagPlusRecoPassWWTightIdIso_binned[q][e][p]++;                  
                  } else {
                  EventCount_TagPlusRecoFailWWTightIdIso_binned[q][e][p]++;
                  }
                  
                  //*****************************************************************************************************
                  //WWTight -> HLT
                  //*****************************************************************************************************	    
                  if (passElectronCuts(probe)) {
                    if (passHLT(probe, info->runNum)) {	 		     
                      EventCount_TagPlusWWTightIdIsoPassHLT_binned[q][e][p]++;
                    } else {
                      EventCount_TagPlusWWTightIdIsoFailHLT_binned[q][e][p]++;
                    }
                  }	    	    	
                }
              }
            }
          }            

          //*****************************************************************************************************
          //Reco -> WWTight
          //*****************************************************************************************************	    
          if (passElectronCuts(probe)) {	 
            EventCount_TagPlusRecoPassWWTightIdIso++;
            TagPlusRecoPassWWTightIdIso->Fill(dilepton.M());	      		  	      
          } else {
            EventCount_TagPlusRecoFailWWTightIdIso++;
            TagPlusRecoFailWWTightIdIso->Fill(dilepton.M());	   
          }
	    
          //*****************************************************************************************************
          //WWTight -> HLT
          //*****************************************************************************************************	    
          if (passElectronCuts(probe)) {
            if (passHLT(probe, info->runNum)) {	 
              EventCount_TagPlusWWTightIdIsoPassHLT++;
              TagPlusWWTightIdIsoPassHLT->Fill(dilepton.M());	
              ProbePt->Fill(probe->pt);
            } else {
              EventCount_TagPlusWWTightIdIsoFailHLT++;
              TagPlusWWTightIdIsoPassHLT->Fill(dilepton.M());			
            }
          }

        } //passes T&P selection

              
      } //loop over probes
    } //loop over tags


  } //end loop over data     

  delete info;
  delete electronArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);
  TagPlusRecoPassWWTightIdIso->Draw();
  cv->SaveAs("TagPlusRecoPassWWTightIdIso.gif");
  TagPlusRecoFailWWTightIdIso->Draw();
  cv->SaveAs("TagPlusRecoFailWWTightIdIso.gif");
  TagPlusWWTightIdIsoPassHLT->Draw();
  cv->SaveAs("TagPlusWWTightIdIsoPassHLT.gif");
  TagPlusWWTightIdIsoFailHLT->Draw();
  cv->SaveAs("TagPlusWWTightIdIsoFailHLT.gif");

  ProbePt->Draw();
  cv->SaveAs("ProbePt.gif");
  TagPt->Draw();
  cv->SaveAs("TagPt.gif");

  //*****************************************************************************************
  //Summarize Efficiencies
  //*****************************************************************************************
  cout << "**********************************************************************\n";
  cout << "Summarize MC Efficiencies\n";
  TFile *file = new TFile("Efficiencies.root", "UPDATE");

  for (int q=0; q < chargebinLabel.size(); ++q) {    
    for (int e=0; e < etabinLabel.size(); ++e) {
      vector<Double_t> efficiency_RecoToWWTight;
      vector<Double_t> efficiency_RecoToWWTight_lowErr;
      vector<Double_t> efficiency_RecoToWWTight_highErr;
      vector<Double_t> efficiency_WWTightToHLT;
      vector<Double_t> efficiency_WWTightToHLT_lowErr;
      vector<Double_t> efficiency_WWTightToHLT_highErr;
      
      for (int p=0; p < ptbinLabel.size(); ++p) {
    
        cout << etabinLabel[e] << "  " << ptbinLabel[p] << "  Charge = " << chargeSelection[q] << endl;
        cout << endl;

        Int_t errorType = 2; //Clopper Pearson intervals
        Double_t ratio;
        Double_t errLow;
        Double_t errHigh;     
        Double_t n1;
        Double_t n2;

        n1 = TMath::Nint(EventCount_TagPlusRecoPassWWTightIdIso_binned[q][e][p]);
        n2 = TMath::Nint(EventCount_TagPlusRecoPassWWTightIdIso_binned[q][e][p] + EventCount_TagPlusRecoFailWWTightIdIso_binned[q][e][p]);
        if (n1 > n2) n1 = n2;
        mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);      
        cout << "Reco -> WWTight  Efficiency (MC) : " << n1 << " / " << n2 << " = " << ratio << " + " << errHigh << " - " << errLow << endl;
        cout << "\n";
        efficiency_RecoToWWTight.push_back(ratio);
        efficiency_RecoToWWTight_lowErr.push_back(errLow);
        efficiency_RecoToWWTight_highErr.push_back(errHigh);


        n1 = TMath::Nint(EventCount_TagPlusWWTightIdIsoPassHLT_binned[q][e][p]);
        n2 = TMath::Nint(EventCount_TagPlusWWTightIdIsoPassHLT_binned[q][e][p] + EventCount_TagPlusWWTightIdIsoFailHLT_binned[q][e][p]);
        if (n1 > n2) n1 = n2;
        mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);      
        cout << "WWTight  -> HLT : " << n1 << " / " << n2 << " = " << ratio << " + " << errHigh << " - " << errLow << endl;
        cout << "\n";

        efficiency_WWTightToHLT.push_back(ratio);
        efficiency_WWTightToHLT_lowErr.push_back(errLow);
        efficiency_WWTightToHLT_highErr.push_back(errHigh);

  
        cout << "\n\n\n";
      }

    //Make Efficiency Graphs
    const Int_t nbins = efficiency_RecoToWWTight.size();
    Double_t x[nbins];
    Double_t y[nbins];
    Double_t xErr[nbins];
    Double_t yErrLow[nbins];
    Double_t yErrHigh[nbins];

    //***********************************************************
    //Reco -> WW Tight    
    //***********************************************************
    for (int b=0; b< efficiency_RecoToWWTight.size(); ++b) {
      x[b] = (ptbinUpEdge[b] + ptbinLowEdge[b]) / 2;
      xErr[b] = fabs(ptbinUpEdge[b] - ptbinLowEdge[b]) / 2;
      if (b==efficiency_RecoToWWTight.size() - 1) {
        x[b] = (100 + ptbinLowEdge[b]) / 2;
        xErr[b] = fabs(100 - ptbinLowEdge[b]) / 2;
      }
      y[b] = efficiency_RecoToWWTight[b];
      yErrLow[b] = fabs(efficiency_RecoToWWTight_lowErr[b]);
      yErrHigh[b] = efficiency_RecoToWWTight_highErr[b];
    }
   
    TGraphAsymmErrors *efficiencyGraph_RecoToWWTight = new TGraphAsymmErrors(nbins, x, y, xErr, xErr, yErrLow,yErrHigh );
    efficiencyGraph_RecoToWWTight->SetName(("MCEfficiency_RecoToWWTight_" + chargebinLabel[q] + "_" + etabinLabel[e]).c_str());
    efficiencyGraph_RecoToWWTight->SetTitle("");
    efficiencyGraph_RecoToWWTight->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    efficiencyGraph_RecoToWWTight->GetYaxis()->SetTitle("Efficiency");
    efficiencyGraph_RecoToWWTight->SetMaximum(1.0);
    efficiencyGraph_RecoToWWTight->SetMinimum(0.0);
    efficiencyGraph_RecoToWWTight->SetMarkerSize(1);
    efficiencyGraph_RecoToWWTight->SetLineWidth(2);
    efficiencyGraph_RecoToWWTight->GetXaxis()->SetRangeUser(0,100);

    file->WriteTObject(efficiencyGraph_RecoToWWTight, efficiencyGraph_RecoToWWTight->GetName(), "WriteDelete");


    //***********************************************************
    //WW Tight -> HLT
    //***********************************************************
    for (int b=0; b< efficiency_WWTightToHLT.size(); ++b) {
      x[b] = (ptbinUpEdge[b] + ptbinLowEdge[b]) / 2;
      xErr[b] = fabs(ptbinUpEdge[b] - ptbinLowEdge[b]) / 2;
      if (b==efficiency_RecoToWWTight.size() - 1) {
        x[b] = (100 + ptbinLowEdge[b]) / 2;
        xErr[b] = fabs(100 - ptbinLowEdge[b]) / 2;
      }
      y[b] = efficiency_WWTightToHLT[b];
      yErrLow[b] = fabs(efficiency_WWTightToHLT_lowErr[b]);
      yErrHigh[b] = efficiency_WWTightToHLT_highErr[b];
    }      
    TGraphAsymmErrors *efficiencyGraph_WWTightToHLT = new TGraphAsymmErrors(nbins, x, y, xErr, xErr, yErrLow,yErrHigh );
    efficiencyGraph_WWTightToHLT->SetName(("MCEfficiency_WWTightToHLT_" + chargebinLabel[q] + "_" + etabinLabel[e]).c_str());
    efficiencyGraph_WWTightToHLT->SetTitle("");
    efficiencyGraph_WWTightToHLT->GetXaxis()->SetTitle("p_{T} [GeV/c]");
    efficiencyGraph_WWTightToHLT->GetYaxis()->SetTitle("Efficiency");
    efficiencyGraph_WWTightToHLT->SetMaximum(1.0);
    efficiencyGraph_WWTightToHLT->SetMinimum(0.0);
    efficiencyGraph_WWTightToHLT->SetMarkerSize(1);
    efficiencyGraph_WWTightToHLT->SetLineWidth(2);
    efficiencyGraph_WWTightToHLT->GetXaxis()->SetRangeUser(0,100);

    file->WriteTObject(efficiencyGraph_WWTightToHLT, efficiencyGraph_WWTightToHLT->GetName(), "WriteDelete");

    } //end for loop over eta bins
  } //end for loop over charge bins

  file->Close();

  Int_t errorType = 2; //Clopper Pearson intervals
  Double_t ratio;
  Double_t errLow;
  Double_t errHigh;     
  Double_t n1;
  Double_t n2;


  n1 = TMath::Nint(EventCount_TagPlusRecoPassWWTightIdIso);
  n2 = TMath::Nint(EventCount_TagPlusRecoPassWWTightIdIso + EventCount_TagPlusRecoFailWWTightIdIso);
  if (n1 > n2) n1 = n2;
  mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);      
  cout << "Reco -> WWTight  Efficiency (MC) : " << n1 << " / " << n2 << " = " << ratio << " + " << errHigh << " - " << errLow << endl;
  cout << "\n";



  n1 = TMath::Nint(EventCount_TagPlusWWTightIdIsoPassHLT);
  n2 = TMath::Nint(EventCount_TagPlusWWTightIdIsoPassHLT + EventCount_TagPlusWWTightIdIsoFailHLT);
  if (n1 > n2) n1 = n2;
  mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);      
  cout << "WWTight  -> HLT : " << n1 << " / " << n2 << " = " << ratio << " + " << errHigh << " - " << errLow << endl;
  cout << "\n";

  

  cout << "**********************************************************************\n";

      
  gBenchmark->Show("ElectronTagAndProbe");       


} 


void WriteMassToFile( ofstream *file, Double_t mass) {

  (*file) << left << mass << endl;
}



Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum, Bool_t isMC) {


  Bool_t pass = kFALSE;

  //it's electron data
  if (isMC) {
    if ( (ele->hltMatchBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
  } else {
    if ((runNum >= 136033) && (runNum <= 139980)) {
      if ( !(ele->hltMatchBits & kHLT_Mu9) && (ele->hltMatchBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 140058) && (runNum <= 141882)) {
      if ( !(ele->hltMatchBits & kHLT_Mu9) && (ele->hltMatchBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 141956) && (runNum <= 144114)) {
      if ( !(ele->hltMatchBits & kHLT_Mu9) && (ele->hltMatchBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
    } 
    if ((runNum >= 146428) && (runNum <= 147116)) {
      if ( !(ele->hltMatchBits & kHLT_Mu9) && (ele->hltMatchBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 147196) && (runNum <= 148058)) {
      if ( !(ele->hltMatchBits & kHLT_Mu15) && (ele->hltMatchBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
    }
    if ((runNum >= 148819) && (runNum <= 149442)) {
      if ( !(ele->hltMatchBits & kHLT_Mu15) && (ele->hltMatchBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
    } 
      
  }

  return pass;

}


Bool_t passElectronTagCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;
  if (!passElectronCuts(ele)) pass = kFALSE;

  return pass;
}

Bool_t passTighterElectronCuts(const mithep::TElectron *ele) {
  
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
            && fabs(ele->deltaPhiIn) < 0.025
            && ele->HoverE < 0.025
            && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.03
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
             && fabs(ele->deltaEtaIn) < 0.005
             && fabs(ele->deltaPhiIn) < 0.02
             && ele->HoverE < 0.025
             && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.02
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
