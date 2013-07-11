//root -l -b -q EWKAna/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodSwitches.h"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"

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
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType);
void MakeNtuple(const string inputFilename,  const string outputFilename);

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
void MakeWJetsMuonTrainingNtuple() {

  MakeNtuple("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-h115ww2l-gf-v11-pu_noskim_0000.root","MuonSelectionTraining.HWW115.root");
  MakeNtuple("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-h120ww2l-gf-v11-pu_noskim_0000.root","MuonSelectionTraining.HWW120.root");
  MakeNtuple("/home/sixie/hist/AllNtuple/HWWNtuple/mc/AllNtuple_HWWNtuple_s11-h130ww2l-gf-v11-pu_noskim_0000.root","MuonSelectionTraining.HWW130.root");


}


void MakeNtuple(const string inputFilename, const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 26.5;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  TTree *muTree = new TTree("Muons","Muons");
  muTree->SetAutoFlush(0);

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergyOverPt; 
  Float_t                 fMuHoEnergyOverPt; 
  Float_t                 fMuEmEnergyOverPt; 
  Float_t                 fMuHadS9EnergyOverPt; 
  Float_t                 fMuHoS9EnergyOverPt; 
  Float_t                 fMuEmS9EnergyOverPt; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 


  muTree->Branch("weight",&fWeight,"weight/F");
  muTree->Branch("run",&fRunNumber,"run/i");
  muTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  muTree->Branch("event",&fEventNumber,"event/i");
  muTree->Branch("pt",&fMuPt,"pt/F"); 
  muTree->Branch("eta",&fMuEta,"eta/F"); 
  muTree->Branch("phi",&fMuPhi,"phi/F"); 
  muTree->Branch("pfiso",&fMuPFIso,"pfiso/F"); 
  
  //CutBased Variables
  muTree->Branch("TkNchi2",&fMuTkNchi2,"TkNchi2/F"); 
  muTree->Branch("GlobalNchi2",&fMuGlobalNchi2,"GlobalNchi2/F"); 
  muTree->Branch("NValidHits",&fMuNValidHits,"NValidHits/F"); 
  muTree->Branch("NTrackerHits",&fMuNTrackerHits,"NTrackerHits/F"); 
  muTree->Branch("NPixelHits",&fMuNPixelHits,"NPixelHits/F"); 
  muTree->Branch("NMatches",&fMuNMatches,"NMatches/F"); 
  muTree->Branch("D0",&fMuD0,"D0/F"); 

  //Additional Vars used in Likelihood
  muTree->Branch("IP3d",&fMuIP3d,"IP3d/F"); 
  muTree->Branch("IP3dSig",&fMuIP3dSig,"IP3dSig/F"); 
  muTree->Branch("TrkKink",&fMuTrkKink,"TrkKink/F"); 
  muTree->Branch("GlobalKink",&fMuGlobalKink,"GlobalKink/F"); 
  muTree->Branch("SegmentCompatibility",&fMuSegmentCompatibility,"SegmentCompatibility/F"); 
  muTree->Branch("CaloCompatibility",&fMuCaloCompatibility,"CaloCompatibility/F"); 
  muTree->Branch("HadEnergyOverPt",&fMuHadEnergyOverPt,"HadEnergyOverPt/F"); 
  muTree->Branch("HoEnergyOverPt",&fMuHoEnergyOverPt,"HoEnergyOverPt/F"); 
  muTree->Branch("EmEnergyOverPt",&fMuEmEnergyOverPt,"EmEnergyOverPt/F"); 
  muTree->Branch("HadS9EnergyOverPt",&fMuHadS9EnergyOverPt,"HadS9EnergyOverPt/F"); 
  muTree->Branch("HoS9EnergyOverPt",&fMuHoS9EnergyOverPt,"HoS9EnergyOverPt/F"); 
  muTree->Branch("EmS9EnergyOverPt",&fMuEmS9EnergyOverPt,"EmS9EnergyOverPt/F"); 

  //Isolation Variables
  muTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  muTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  muTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  muTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  muTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  muTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  muTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  muTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  muTree->Branch("Rho",&fRho,"Rho/F"); 
  muTree->Branch("NVertices",&fNVertices,"NVertices/F"); 


  UInt_t NMuonsFilled = 0;
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
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  Int_t NEvents = 0;

  vector<string> inputfiles;
  if (inputFilename == "LIST") {
    inputfiles.push_back("/data/blue/sixie/ntuples/HWW/mc/HwwAnalysis_s11-h115ww2l-gf-v11-pu_noskim_normalized.root");
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
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      NEvents++;

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

      Int_t NElectrons = electronArr->GetEntries();

      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);

        //********************************************************
        // Select MC Truth Electrons
        //********************************************************
              
        //use only real electrons
        if (!((UInt_t(abs(max(0,mu->isMCReal))) & 2) == 2)) continue;
              
        if (!passMuonDenominatorCuts(mu, 2)) continue;

       //Fill These Electrons

          //Fill These Muons

          fWeight = 1.0;
          fRunNumber = info->runNum;
          fLumiSectionNumber = info->lumiSec;
          fEventNumber = info->evtNum;
          fMuPt = mu->pt; 
          fMuEta = mu->eta; 
          fMuPhi = mu->phi; 
          fMuPFIso = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
   
          //CutBased Variables
          fMuTkNchi2 = mu->tkNchi2 ; 
          fMuGlobalNchi2 = mu->muNchi2 ; 
          fMuNValidHits = mu->nValidHits; 
          fMuNTrackerHits = mu->nTkHits; 
          fMuNPixelHits = mu->nPixHits; 
          fMuNMatches = mu->nMatch ; 
          fMuD0 = mu->d0 ; 

          //Additional Vars 
          fMuIP3d = mu->ip3d ; 
          fMuIP3dSig = mu->ip3dSig ; 
          fMuTrkKink = mu->TrkKink ; 
          fMuGlobalKink = mu->GlobalKink ; 
          fMuSegmentCompatibility = mu->SegmentCompatibility ; 
          fMuCaloCompatibility = mu->CaloCompatilibity ; 
          fMuHadEnergyOverPt = mu->HadEnergy / mu->pt; 
          fMuHoEnergyOverPt = mu->HoEnergy / mu->pt; 
          fMuEmEnergyOverPt = mu->EmEnergy / mu->pt; 
          fMuHadS9EnergyOverPt = mu->HadS9Energy / mu->pt; 
          fMuHoS9EnergyOverPt = mu->HoS9Energy / mu->pt; 
          fMuEmS9EnergyOverPt = mu->EmS9Energy / mu->pt; 



          //Isolation Variables
          fMuChargedIso03 = mu->ChargedIso03 ; 
          fMuChargedIso03FromOtherVertices = mu->ChargedIso03FromOtherVertices ; 
          fMuNeutralIso03_05Threshold = mu->NeutralIso03_05Threshold ; 
          fMuNeutralIso03_10Threshold = mu->NeutralIso03_10Threshold ; 
          fMuChargedIso04 = mu->ChargedIso04 ; 
          fMuChargedIso04FromOtherVertices = mu->ChargedIso04FromOtherVertices ; 
          fMuNeutralIso04_05Threshold = mu->NeutralIso04_05Threshold ; 
          fMuNeutralIso04_10Threshold = mu->NeutralIso04_10Threshold ; 
          fRho = info->PileupEnergyDensity; 
          fNVertices = info->nPV0; 

          NMuonsFilled++;
          muTree->Fill();

      } //loop over electrons

    }

    cout << "Total Muons: " << NMuonsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Total Muons: " << NMuonsFilled << endl;
  outputFile->Write();
  outputFile->Close();

  gBenchmark->Show("WWTemplate");       
} 




//--------------------------------------------------------------------------------------------------
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType) {
  
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
          && (mu->ChargedIso03 + mu->NeutralIso03_10Threshold) / mu->pt < 1.0
          && (mu->pterr / mu->pt < 0.1)
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
          && (mu->ChargedIso03 + mu->NeutralIso03_10Threshold) / mu->pt < 0.4
          && ( mu->pterr / mu->pt < 0.1)
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
          )
      ) pass = kFALSE;    
  }


  if (DenominatorType == 100) {
    pass = kTRUE;
  }


  return pass;
}

