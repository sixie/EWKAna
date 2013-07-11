//root -l EWKAna/Hww/Efficiency/ElectronTagAndProbe.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10a-eg-s17_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Efficiency/ElectronTagAndProbe.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_r10b-el-pr-v2_noskim.root\",\"\"\)
//root -l EWKAna/Hww/Efficiency/ElectronTagAndProbe.C+\(\"/home/sixie/hist/WWAnalysis/WWAnalysis_ele_noskim.root\",\"\"\)



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
Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum);
Bool_t passElectronTagCuts(const mithep::TElectron *ele);
Bool_t passTighterElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronCuts(const mithep::TElectron *ele);
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


//=== MAIN MACRO =================================================================================================

void ElectronTagAndProbe(const string dataInputFilename, const string Label) {   

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
//    ptbinLabel.push_back("Pt10To15"); ptbinLowEdge.push_back(10.0); ptbinUpEdge.push_back(15.0);
//    ptbinLabel.push_back("Pt15To20"); ptbinLowEdge.push_back(15.0); ptbinUpEdge.push_back(20.0);
//   ptbinLabel.push_back("Pt20To30"); ptbinLowEdge.push_back(20.0); ptbinUpEdge.push_back(30.0);
//   ptbinLabel.push_back("Pt30To40"); ptbinLowEdge.push_back(30.0); ptbinUpEdge.push_back(40.0);
//   ptbinLabel.push_back("Pt40To50"); ptbinLowEdge.push_back(40.0); ptbinUpEdge.push_back(50.0);
//   ptbinLabel.push_back("Pt50ToInf"); ptbinLowEdge.push_back(50.0); ptbinUpEdge.push_back(14000.0);
  vector<string> etabinLabel;
  vector<Double_t> etabinLowEdge;
  vector<Double_t> etabinUpEdge;
  etabinLabel.push_back("EB"); etabinLowEdge.push_back(0); etabinUpEdge.push_back(1.5);
  etabinLabel.push_back("EE"); etabinLowEdge.push_back(1.5); etabinUpEdge.push_back(2.5);
  vector<vector<ofstream *> > outputFile_TagPlusRecoFailWWTightIdIso_binned;
  vector<vector<ofstream *> > outputFile_TagPlusRecoPassWWTightIdIso_binned;
  vector<vector<ofstream *> > outputFile_TagPlusWWTightIdIsoPassHLT_binned;
  vector<vector<ofstream *> > outputFile_TagPlusWWTightIdIsoFailHLT_binned;
  for (int i=0; i < int(etabinLabel.size()); ++i) {
    vector<ofstream *> tmp_outputFile_TagPlusRecoFailWWTightIdIso_binned;
    vector<ofstream *> tmp_outputFile_TagPlusRecoPassWWTightIdIso_binned;
    vector<ofstream *> tmp_outputFile_TagPlusWWTightIdIsoPassHLT_binned;
    vector<ofstream *> tmp_outputFile_TagPlusWWTightIdIsoFailHLT_binned;    
    
    for (int j=0; j < int(ptbinLabel.size()); ++j) {
      ofstream *tmp_outputFile_TagPlusRecoFailWWTightIdIso = new ofstream;
      tmp_outputFile_TagPlusRecoFailWWTightIdIso->open(("Mass_TagPlusRecoFailWWTightIdIso_" + etabinLabel[i] + "_" + ptbinLabel[j]+ label).c_str());
      tmp_outputFile_TagPlusRecoFailWWTightIdIso_binned.push_back(tmp_outputFile_TagPlusRecoFailWWTightIdIso);

      ofstream *tmp_outputFile_TagPlusRecoPassWWTightIdIso = new ofstream;
      tmp_outputFile_TagPlusRecoPassWWTightIdIso->open(("Mass_TagPlusRecoPassWWTightIdIso_" + etabinLabel[i] + "_" + ptbinLabel[j]+ label).c_str());
      tmp_outputFile_TagPlusRecoPassWWTightIdIso_binned.push_back(tmp_outputFile_TagPlusRecoPassWWTightIdIso);

      ofstream *tmp_outputFile_TagPlusWWTightIdIsoFailHLT = new ofstream;
      tmp_outputFile_TagPlusWWTightIdIsoFailHLT->open(("Mass_TagPlusWWTightIdIsoFailHLT_" + etabinLabel[i] + "_" + ptbinLabel[j]+ label).c_str());
      tmp_outputFile_TagPlusWWTightIdIsoFailHLT_binned.push_back(tmp_outputFile_TagPlusWWTightIdIsoFailHLT);

      ofstream *tmp_outputFile_TagPlusWWTightIdIsoPassHLT = new ofstream;
      tmp_outputFile_TagPlusWWTightIdIsoPassHLT->open(("Mass_TagPlusWWTightIdIsoPassHLT_" + etabinLabel[i] + "_" + ptbinLabel[j]+ label).c_str());
      tmp_outputFile_TagPlusWWTightIdIsoPassHLT_binned.push_back(tmp_outputFile_TagPlusWWTightIdIsoPassHLT);

      cout << "make bins : " << i << " " << j << endl;

    }

    outputFile_TagPlusRecoFailWWTightIdIso_binned.push_back(tmp_outputFile_TagPlusRecoFailWWTightIdIso_binned);
    outputFile_TagPlusRecoPassWWTightIdIso_binned.push_back(tmp_outputFile_TagPlusRecoPassWWTightIdIso_binned);
    outputFile_TagPlusWWTightIdIsoPassHLT_binned.push_back(tmp_outputFile_TagPlusWWTightIdIsoPassHLT_binned);
    outputFile_TagPlusWWTightIdIsoFailHLT_binned.push_back(tmp_outputFile_TagPlusWWTightIdIsoFailHLT_binned);
       
  }

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  ofstream *outputFile_TagPlusRecoFailWWTightIdIso = new ofstream;
  outputFile_TagPlusRecoFailWWTightIdIso->open(("Mass_TagPlusRecoFailWWTightIdIso" + label).c_str());
  ofstream *outputFile_TagPlusRecoPassWWTightIdIso = new ofstream;
  outputFile_TagPlusRecoPassWWTightIdIso->open(("Mass_TagPlusRecoPassWWTightIdIso" + label).c_str());
  ofstream *outputFile_TagPlusWWTightIdIsoFailHLT = new ofstream;
  outputFile_TagPlusWWTightIdIsoFailHLT->open(("Mass_TagPlusWWTightIdIsoFailHLT" + label).c_str());
  ofstream *outputFile_TagPlusWWTightIdIsoPassHLT = new ofstream;
  outputFile_TagPlusWWTightIdIsoPassHLT->open(("Mass_TagPlusWWTightIdIsoPassHLT" + label).c_str());

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
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_132440-147454_7TeV_StreamExpress_Collisions10_JSON.txt"); 
  rlrm.AddJSONFile("merged_JsonReRecoSep17_JsonStreamExpressV2.txt"); 
  hasJSON = kFALSE;

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
		
    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...    

//     if (info->runNum != 149294) continue;

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
          for (int e=0; e < etabinLabel.size(); ++e) {
            for (int p=0; p < ptbinLabel.size(); ++p) {
	      
              //Require the probe is in the right pt and eta bin
              if ( 		  
                (probe->pt >= ptbinLowEdge[p] && probe->pt < ptbinUpEdge[p])
                &&
                (fabs(probe->eta) >= etabinLowEdge[e] && fabs(probe->eta) < etabinUpEdge[e])                
                &&
                met.Pt() < 20.0
                ) {
                
                //*****************************************************************************************************
                //Reco -> WWTight
                //*****************************************************************************************************	    
                
                if (passElectronCuts(probe)) {	 
                  WriteMassToFile(outputFile_TagPlusRecoPassWWTightIdIso_binned[e][p], dilepton.M());
                } else {
                  WriteMassToFile(outputFile_TagPlusRecoFailWWTightIdIso_binned[e][p], dilepton.M());
                }
	
                //*****************************************************************************************************
                //WWTight -> HLT
                //*****************************************************************************************************	    
                if (passElectronCuts(probe)) {
                  if (passHLT(probe, info->runNum)) {	 		     
                    WriteMassToFile(outputFile_TagPlusWWTightIdIsoPassHLT_binned[e][p], dilepton.M());
                  } else {
                    WriteMassToFile(outputFile_TagPlusWWTightIdIsoPassHLT_binned[e][p], dilepton.M());
                  }
                }	    	    	
              }
            }
          }

          //*****************************************************************************************************
          //Reco -> WWTight
          //*****************************************************************************************************	    
	    
          if (passElectronCuts(probe)) {	 
            WriteMassToFile(outputFile_TagPlusRecoPassWWTightIdIso, dilepton.M());
            TagPlusRecoPassWWTightIdIso->Fill(dilepton.M());	      		  	      
          } else {
            WriteMassToFile(outputFile_TagPlusRecoFailWWTightIdIso, dilepton.M());
            TagPlusRecoFailWWTightIdIso->Fill(dilepton.M());	   
          }
	    
          //*****************************************************************************************************
          //WWTight -> HLT
          //*****************************************************************************************************	    
          if (passElectronCuts(probe)) {
            if (passHLT(probe, info->runNum)) {	 
              WriteMassToFile(outputFile_TagPlusWWTightIdIsoPassHLT, dilepton.M());
              TagPlusWWTightIdIsoPassHLT->Fill(dilepton.M());	
              ProbePt->Fill(probe->pt);
            } else {
              WriteMassToFile(outputFile_TagPlusWWTightIdIsoPassHLT, dilepton.M());
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

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 

      
  gBenchmark->Show("ElectronTagAndProbe");       


} 


void WriteMassToFile( ofstream *file, Double_t mass, Int_t eventNum) {

  (*file) << mass << " " << eventNum << endl;
}

void WriteMassToFile( ofstream *file, Double_t mass) {

  (*file) << left << mass <<  endl;
}



Bool_t passHLT(const mithep::TElectron *ele, Int_t runNum) {


  Bool_t pass = kFALSE;

  //it's electron data
  if ((runNum >= 136033) && (runNum <= 139980)) {
    if ( (ele->hltMatchBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
  } 
  if ((runNum >= 140058) && (runNum <= 141882)) {
    if ( (ele->hltMatchBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
  } 
  if ((runNum >= 141956) && (runNum <= 144114)) {
    if ( (ele->hltMatchBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
  } 
  if ((runNum >= 146428) && (runNum <= 147116)) {
    if ( (ele->hltMatchBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
  }
  if ((runNum >= 147196) && (runNum <= 148058)) {
    if ( (ele->hltMatchBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
  }
  if ((runNum >= 148819) && (runNum <= 149442)) {
    if ( (ele->hltMatchBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
  } 
      
  return pass;

}


Bool_t passElectronTagCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;
  if (!passTighterElectronCuts(ele)) pass = kFALSE;

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
