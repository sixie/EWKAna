//root -l EWKAna/HZZ4l/Selection/HZZ4lSelection.C+\(\"\"\)

//================================================================================================
//
// HZZ4l selection macro
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
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "EWKAna/HZZ4l/Selection/LeptonSelection.hh"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

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

void HZZ4lSelection(const string Label = "") 
{  
  gBenchmark->Start("HZZTemplate");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  Bool_t printDebug = kFALSE;
  Double_t lumi = 1092;              // luminosity (pb^-1)


  vector<vector<string> > inputFiles;
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/data/HZZNtuple_r11a-data.1092ipb.FourRecoLeptonSkim.root");

  inputFiles.push_back(vector<string>());
  inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/mc/HZZ4lAnalysis_s11-zz-v11-pu_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/mc/HZZ4lAnalysis_s11-wz-v11-pu_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/mc/HZZ4lAnalysis_s11-ttj-v11-pu_noskim_normalized.root");
//   inputFiles.push_back(vector<string>());
//   inputFiles.back().push_back("/data/blue/sixie/ntuples/HZZ4l/mc/HZZ4lAnalysis_s11-wwj-v11-pu_noskim_normalized.root");


  vector<string> processNames;
//   processNames.push_back("Data");
  processNames.push_back("ZZ");
//   processNames.push_back("WZ");
//   processNames.push_back("ttbar");
//   processNames.push_back("WW"); 

  assert(processNames.size() == inputFiles.size());

  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  vector <TH1F*>  fHZZ4lSelection; 
  vector <TH1F*>  fHZZTo4ESelection; 
  vector <TH1F*>  fHZZTo4MuSelection; 
  vector <TH1F*>  fHZZTo2E2MuSelection; 
  vector <TH1F*>  fZ1Mass; 
  vector <TH1F*>  fZ2Mass; 
  vector <TH1F*>  fZZMass; 
  for (int q=0; q<processNames.size() ; ++q) {
    TH1F *tmpHZZ4lSelection = new TH1F(("hHZZ4lSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4ESelection = new TH1F(("hHZZTo4ESelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo4MuSelection = new TH1F(("hHZZTo4MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpHZZTo2E2MuSelection = new TH1F(("hHZZTo2E2MuSelection"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 15, -1.5, 13.5);
    TH1F *tmpZ1Mass = new TH1F(("hZ1Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZ2Mass = new TH1F(("hZ2Mass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 200);
    TH1F *tmpZZMass = new TH1F(("hZZMass"+string("_")+processNames[q]).c_str(), ";Cut Number;Number of Events", 40, 0, 600);
    tmpHZZ4lSelection->Sumw2();
    tmpHZZTo4ESelection->Sumw2();
    tmpHZZTo4MuSelection->Sumw2();
    tmpHZZTo2E2MuSelection->Sumw2();   
    tmpZ1Mass->Sumw2();   
    tmpZ2Mass->Sumw2();   
    tmpZZMass->Sumw2();   
    
    fHZZ4lSelection.push_back(tmpHZZ4lSelection);
    fHZZTo4ESelection.push_back(tmpHZZTo4ESelection);
    fHZZTo4MuSelection.push_back(tmpHZZTo4MuSelection);
    fHZZTo2E2MuSelection.push_back(tmpHZZTo2E2MuSelection);
    fZ1Mass.push_back(tmpZ1Mass);
    fZ2Mass.push_back(tmpZ2Mass);
    fZZMass.push_back(tmpZZMass);
  }


  vector<string> CutLabel;
  CutLabel.push_back("Dilepton20/10");

  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TFile *inputFile=0;
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  
  for (int q = 0; q<inputFiles.size() ; ++q) { 
    for (int f = 0; f < inputFiles[q].size() ; ++f) {
      //********************************************************
      // Get Tree
      //********************************************************
      cout << "Reading File " << inputFiles[q][f] << endl;
      eventTree = getTreeFromFile(inputFiles[q][f].c_str(),"Events"); 
      TBranch *infoBr;
      TBranch *electronBr;
      TBranch *muonBr;
      TBranch *jetBr;


      //*****************************************************************************************
      //Loop over muon Data Tree
      //*****************************************************************************************
      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");
  
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
        infoBr->GetEntry(ientry);
        if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
		
        //********************************************************
        // Printdebug
        //********************************************************
        printDebug = kFALSE;
        if ((0 == 1) 
            || (info->evtNum == 9316)           
          ) printDebug = kTRUE;
        



        //********************************************************
        // Load the branches
        //********************************************************
        electronArr->Clear(); 
        muonArr->Clear(); 
        jetArr->Clear(); 
        electronBr->GetEntry(ientry);
        muonBr->GetEntry(ientry);
        jetBr->GetEntry(ientry);

        Double_t eventweight = info->eventweight * lumi;

        //********************************************************
        // Met
        //********************************************************
        TVector3 pfMet;        
        if(info->pfMEx!=0 || info->pfMEy!=0) {       
          pfMet.SetXYZ(info->pfMEx, info->pfMEy, 0);
        }

        //********************************************************
        // Lepton Selection
        //********************************************************
        Int_t NLeptons = 0;
        vector<Int_t> leptonType;
        vector<Int_t> leptonIndex;
        vector<Double_t> leptonPt;
        vector<Double_t> leptonEta;
        vector<Double_t> leptonPhi;
        vector<Int_t> leptonCharge;

        Int_t NJets = 0;
        const mithep::TJet *leadingJet = 0;
    
        for(Int_t i=0; i<muonArr->GetEntries(); i++) {
          const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
          if ( (0==0)
               &&
               passMuonCuts(mu)
               &&
               fabs(mu->eta) < 2.4
               && 
               mu->pt > 5.0
            ) {
            leptonPt.push_back(mu->pt);
            leptonEta.push_back(mu->eta);
            leptonPhi.push_back(mu->phi);
            leptonType.push_back(13);
            leptonIndex.push_back(i);  
            leptonCharge.push_back(mu->q);
          }
        }

        for(Int_t i=0; i<electronArr->GetEntries(); i++) {
          const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
          Bool_t isMuonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if ( leptonType[k] == 13 
                 && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
              ) {
              isMuonOverlap = kTRUE; 
              break;
            }
          }

          if ( (0==0)
               &&
               passElectronCuts(ele)
               &&
               fabs(ele->eta) < 2.5
               && 
               ele->pt > 7.0
               &&
               !isMuonOverlap
            ) {
            leptonPt.push_back(ele->pt);
            leptonEta.push_back(ele->eta);
            leptonPhi.push_back(ele->phi);
            leptonType.push_back(11);
            leptonIndex.push_back(i);
            leptonCharge.push_back(ele->q);
          }
        }

        //sort leptons
        Int_t tempType;
        Int_t tempIndex;
        Double_t tempPt;
        Double_t tempEta;
        Double_t tempPhi;
        Int_t tempCharge;
        for (int l=0; l<leptonIndex.size(); l++) {
          for (int k=0; k < leptonIndex.size() - 1; k++) {
            if (leptonPt[k+1] > leptonPt[k]) {
              tempType = leptonType[k];
              tempIndex = leptonIndex[k];
              tempPt = leptonPt[k];
              tempEta = leptonEta[k];
              tempPhi = leptonPhi[k];
              tempCharge = leptonCharge[k];
          
              leptonType[k] = leptonType[k+1];
              leptonIndex[k] = leptonIndex[k+1];
              leptonPt[k] = leptonPt[k+1];
              leptonEta[k] = leptonEta[k+1];
              leptonPhi[k] = leptonPhi[k+1];
              leptonCharge[k] = leptonCharge[k+1];

              leptonType[k+1] = tempType;
              leptonIndex[k+1] = tempIndex;
              leptonPt[k+1] = tempPt;
              leptonEta[k+1] = tempEta;
              leptonPhi[k+1] = tempPhi;
              leptonCharge[k+1] = tempCharge;
          
            }
          }
        }

        //******************************************************************************
        //Count Jets
        //******************************************************************************
        for(Int_t i=0; i<jetArr->GetEntries(); i++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

          if (jet->pt < 30.0) continue;

          Bool_t leptonOverlap = kFALSE;
          for (int k=0; k<leptonPt.size(); ++k) {
            if (leptonType[k] == 11) {
              const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[leptonIndex[k]]);
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta,tmpEle->phi, tmpEle->eta) < 0.5) {
                leptonOverlap = kTRUE;
              }
            } else {
              if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, leptonPhi[k],leptonEta[k]) < 0.5) {
                leptonOverlap = kTRUE;
              }
            }
          }

          if (!leptonOverlap) {
            if (jet->pt > 30.0 && fabs(jet->eta) < 5.0 ) {
              if (!leadingJet || jet->pt > leadingJet->pt) {
                leadingJet = jet;
              }
              NJets++;
            } 
          }
        }



        //******************************************************************************
        //Z1 Selection
        //******************************************************************************
        if (leptonPt.size() < 2) continue;

        Int_t Z1LeptonPlusIndex = -1;
        Int_t Z1LeptonMinusIndex = -1;
        Double_t BestZ1Mass = -1;
        for(int i = 0; i < leptonPt.size(); ++i) {
          for(int j = i+1; j < leptonPt.size(); ++j) {
            if (!(leptonPt[i] > 20.0 || leptonPt[j] > 20.0)) continue; //one lepton must be 20 GeV
            if (!(leptonPt[i] > 10.0 && leptonPt[j] > 10.0)) continue; //both leptons must be 10 GeV
            if (leptonCharge[i] == leptonCharge[j]) continue;          //require opp sign
            if (fabs(leptonType[i]) != fabs(leptonType[j])) continue;  //require same flavor

            //Make Z1 hypothesis
            mithep::FourVectorM leptonPlus;
            mithep::FourVectorM leptonMinus;

            Double_t Lepton1Mass = 0.51099892e-3;
            if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
            Double_t Lepton2Mass = 0.51099892e-3;
            if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

            if (leptonCharge[i] > 0 ) {
              leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
            } else {
              leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
            }
            mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
            
            if (BestZ1Mass < 0 || fabs(dilepton.M() - 91.1876) < fabs(BestZ1Mass - 91.1876)) {
              if (dilepton.M() > 60) {
                BestZ1Mass = dilepton.M();
                if (leptonCharge[i] > 0) {
                  Z1LeptonPlusIndex = i;
                  Z1LeptonMinusIndex = j;
                } else {
                  Z1LeptonPlusIndex = j;
                  Z1LeptonMinusIndex = i;
                }
              }
            }
          }
        }
        // stop if no Z1 candidate is found
        if (Z1LeptonPlusIndex == -1) continue;

        //******************************************************************************
        //Make Z1 Candidate FourVector
        //******************************************************************************
        mithep::FourVectorM Z1LeptonPlus;
        mithep::FourVectorM Z1LeptonMinus;
        if (leptonType[Z1LeptonPlusIndex] == 11) {
          Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 0.51099892e-3 );
        } else {
          Z1LeptonPlus.SetCoordinates(leptonPt[Z1LeptonPlusIndex], leptonEta[Z1LeptonPlusIndex], leptonPhi[Z1LeptonPlusIndex], 105.658369e-3 );
        }
        if (leptonType[Z1LeptonMinusIndex] == 11) {
          Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 0.51099892e-3 );
        } else {
          Z1LeptonMinus.SetCoordinates(leptonPt[Z1LeptonMinusIndex], leptonEta[Z1LeptonMinusIndex], leptonPhi[Z1LeptonMinusIndex], 105.658369e-3 );
        }
        mithep::FourVectorM Z1Candidate = Z1LeptonPlus+Z1LeptonMinus;
  

        //******************************************************************************
        //Z2 Selection
        //******************************************************************************
        Int_t Z2LeptonPlusIndex = -1;
        Int_t Z2LeptonMinusIndex = -1;
        Double_t BestZ2Mass = -1;
        for(int i = 0; i < leptonPt.size(); ++i) {
          for(int j = i+1; j < leptonPt.size(); ++j) {
            if (i == Z1LeptonPlusIndex || i == Z1LeptonMinusIndex) continue; //skip Z1 leptons
            if (j == Z1LeptonPlusIndex || j == Z1LeptonMinusIndex) continue; //skip Z1 leptons
            if (leptonCharge[i] == leptonCharge[j]) continue;         //require opp sign
            if (fabs(leptonType[i]) != fabs(leptonType[j])) continue; //require same flavor
            
            //Make Z2 hypothesis
            mithep::FourVectorM leptonPlus;
            mithep::FourVectorM leptonMinus;

            Double_t Lepton1Mass = 0.51099892e-3;
            if (leptonType[i] != 11) Lepton1Mass = 105.658369e-3;
            Double_t Lepton2Mass = 0.51099892e-3;
            if (leptonType[j] != 11) Lepton2Mass = 105.658369e-3;

            if (leptonCharge[i] > 0 ) {
              leptonPlus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
              leptonMinus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
            } else {
              leptonPlus.SetCoordinates(leptonPt[j], leptonEta[j], leptonPhi[j], Lepton2Mass );
              leptonMinus.SetCoordinates(leptonPt[i], leptonEta[i], leptonPhi[i], Lepton1Mass );
            }

            mithep::FourVectorM dilepton = leptonPlus+leptonMinus;
            mithep::FourVectorM fourLepton = Z1Candidate + dilepton;

            if (!(dilepton.M() > 12.0)) continue;
            if (!(fourLepton.M() > 100.0)) continue;
            
            //for 4e and 4mu, require at least 1 of the other opp sign lepton pairs have mass > 12
            if (fabs(leptonType[i]) == fabs(leptonType[Z1LeptonPlusIndex])) {
              mithep::FourVectorM pair1 = Z1LeptonPlus+leptonMinus;
              mithep::FourVectorM pair2 = Z1LeptonMinus+leptonPlus;
              if (!(pair1.M() > 12 || pair2.M() > 12)) continue;
            }
            
            //Disambiguiation is done by choosing the pair with the largest ptMax, and largest ptMin
            if (Z2LeptonPlusIndex < 0) {
              if (leptonCharge[i] > 0) {
                Z2LeptonPlusIndex = i;
                Z2LeptonMinusIndex = j;
              } else {
                Z2LeptonPlusIndex = j;
                Z2LeptonMinusIndex = i;
              }
            } else {
              Double_t BestPairPtMax = leptonPt[Z2LeptonPlusIndex];               
              Double_t BestPairPtMin = leptonPt[Z2LeptonMinusIndex]; 
              if (leptonPt[Z2LeptonMinusIndex] > BestPairPtMax) {
                BestPairPtMax = leptonPt[Z2LeptonMinusIndex];
                BestPairPtMin = leptonPt[Z2LeptonPlusIndex];
              }

              Double_t CurrentPairPtMax = leptonPt[i];               
              Double_t CurrentPairPtMin = leptonPt[j]; 
              if (leptonPt[j] > CurrentPairPtMax) {
                CurrentPairPtMax = leptonPt[j];
                CurrentPairPtMin = leptonPt[i];
              }

              if (CurrentPairPtMax > BestPairPtMax) {
                if (leptonCharge[i] > 0) {
                  Z2LeptonPlusIndex = i;
                  Z2LeptonMinusIndex = j;
                } else {
                  Z2LeptonPlusIndex = j;
                  Z2LeptonMinusIndex = i;
                }
              } else if (CurrentPairPtMax  == BestPairPtMax) {
                if (CurrentPairPtMin > BestPairPtMin) {
                  if (leptonCharge[i] > 0) {
                    Z2LeptonPlusIndex = i;
                    Z2LeptonMinusIndex = j;
                  } else {
                    Z2LeptonPlusIndex = j;
                    Z2LeptonMinusIndex = i;
                  }                  
                }
              }
            }            
          }
        }
 
        // stop if no Z1 candidate is found
        if (Z2LeptonPlusIndex == -1) continue;

        //******************************************************************************
        //Make Z2 Candidate FourVector
        //******************************************************************************
        mithep::FourVectorM Z2LeptonPlus;
        mithep::FourVectorM Z2LeptonMinus;
        if (leptonType[Z2LeptonPlusIndex] == 11) {
          Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 0.51099892e-3 );
        } else {
          Z2LeptonPlus.SetCoordinates(leptonPt[Z2LeptonPlusIndex], leptonEta[Z2LeptonPlusIndex], leptonPhi[Z2LeptonPlusIndex], 105.658369e-3 );
        }
        if (leptonType[Z2LeptonMinusIndex] == 11) {
          Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 0.51099892e-3 );
        } else {
          Z2LeptonMinus.SetCoordinates(leptonPt[Z2LeptonMinusIndex], leptonEta[Z2LeptonMinusIndex], leptonPhi[Z2LeptonMinusIndex], 105.658369e-3 );
        }
        mithep::FourVectorM Z2Candidate = Z2LeptonPlus+Z2LeptonMinus;
        mithep::FourVectorM ZZSystem = Z1Candidate + Z2Candidate;



        fZ1Mass[q]->Fill(Z1Candidate.M(),eventweight);
        fZ2Mass[q]->Fill(Z2Candidate.M(),eventweight);
        fZZMass[q]->Fill(ZZSystem.M(),eventweight);

      } //end loop over data     
    }
  }


  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  //--------------------------------------------------------------------------------------------------------------
  // Save Histograms;
  //============================================================================================================== 
  TFile *file = new TFile("HZZ4lSelectionPlots.root", "RECREATE");
  
  for (int q = 0; q<processNames.size() ; ++q) { 
    file->WriteTObject(fHZZ4lSelection[q], fHZZ4lSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHZZTo4ESelection[q], fHZZTo4ESelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHZZTo4MuSelection[q], fHZZTo4MuSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fHZZTo2E2MuSelection[q], fHZZTo2E2MuSelection[q]->GetName(), "WriteDelete");
    file->WriteTObject(fZ1Mass[q], fZ1Mass[q]->GetName(), "WriteDelete");
    file->WriteTObject(fZ2Mass[q], fZ2Mass[q]->GetName(), "WriteDelete");
    file->WriteTObject(fZZMass[q], fZZMass[q]->GetName(), "WriteDelete");
  }

  file->Close();
  delete file;
        
  gBenchmark->Show("WWTemplate");       
} 



