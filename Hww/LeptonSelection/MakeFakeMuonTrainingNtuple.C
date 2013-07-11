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
#include "EWKAna/Utils/LeptonIDCuts.hh"

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
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Double_t fRho);
Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Double_t fRho);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
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

//--------------------------------------------------------------------------------------------------
Bool_t passZVeto(TClonesArray *muonArr, Double_t fRho)
{
  for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
    const mithep::TMuon* mu1 = (mithep::TMuon*)((*muonArr)[i]);
    if(mu1->pt        < 20)  continue;
    if(fabs(mu1->eta) > 2.4) continue;
    if(!passMuonDenominatorCuts(mu1,  fRho)) continue;
    for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) {
      const mithep::TMuon* mu2 = (mithep::TMuon*)((*muonArr)[j]);
      if(mu1->q == mu2->q)     continue;
      if(mu2->pt	< 20)  continue;
      if(fabs(mu2->eta) > 2.4) continue;
      if(!passMuonDenominatorCuts(mu2,  fRho)) continue;

      return kFALSE;  // Z candidate => fail Z veto
    }
  }
  
  return kTRUE;  // No Z candidate => pass Z veto
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
void MakeFakeMuonTrainingNtuple() {

  MakeNtuple("LIST","MuonSelectionTraining.Fake.root");

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
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
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
  muTree->Branch("HadEnergy",&fMuHadEnergy,"HadEnergy/F"); 
  muTree->Branch("HoEnergy",&fMuHoEnergy,"HoEnergy/F"); 
  muTree->Branch("EmEnergy",&fMuEmEnergy,"EmEnergy/F"); 
  muTree->Branch("HadS9Energy",&fMuHadS9Energy,"HadS9Energy/F"); 
  muTree->Branch("HoS9Energy",&fMuHoS9Energy,"HoS9Energy/F"); 
  muTree->Branch("EmS9Energy",&fMuEmS9Energy,"EmS9Energy/F"); 

  //Isolation Variables
  muTree->Branch("ChargedIso03",&fMuChargedIso03,"ChargedIso03/F"); 
  muTree->Branch("ChargedIso03FromOtherVertices",&fMuChargedIso03FromOtherVertices,"ChargedIso03FromOtherVertices/F"); 
  muTree->Branch("NeutralIso03_05Threshold",&fMuNeutralIso03_05Threshold,"NeutralIso03_05Threshold/F"); 
  muTree->Branch("NeutralIso03_10Threshold",&fMuNeutralIso03_10Threshold,"NeutralIso03_10Threshold/F"); 
  muTree->Branch("ChargedIso04",&fMuChargedIso04,"ChargedIso04/F"); 
  muTree->Branch("ChargedIso04FromOtherVertices",&fMuChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  muTree->Branch("NeutralIso04_05Threshold",&fMuNeutralIso04_05Threshold,"NeutralIso04_05Threshold/F"); 
  muTree->Branch("NeutralIso04_10Threshold",&fMuNeutralIso04_10Threshold,"NeutralIso04_10Threshold/F"); 
  muTree->Branch("TrkIso03",&fMuTrkIso03,"TrkIso03/F"); 
  muTree->Branch("EMIso03",&fMuEMIso03,"EMIso03/F"); 
  muTree->Branch("HadIso03",&fMuHadIso03,"HadIso03/F"); 
  muTree->Branch("TrkIso05",&fMuTrkIso05,"TrkIso05/F"); 
  muTree->Branch("EMIso05",&fMuEMIso05,"EMIso05/F"); 
  muTree->Branch("HadIso05",&fMuHadIso05,"HadIso05/F"); 
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
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/data/smurf/data/Winter11_4700ipb/auxiliar/hww.Full2011.json"); 

  Int_t NEvents = 0;

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
    eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

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

      Int_t NMuons = muonArr->GetEntries();
      
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
 
        //find leading jet in the event
        Double_t leadingJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
        }

        //********************************************************
        // Event Selection Cuts
        //********************************************************
              
        //veto events with more than 1 reco muon
        if (NMuons > 1) continue;

        //met cut removed W events
        if (met > 20) continue;
        if (calcMt(pfMet.Pt(), pfMet.Phi(), mu) > 20) continue;

        if (!passZVeto(muonArr, info->PileupEnergyDensity)) continue;
        if (mu->pt > 35) continue;

//         //Jet Threshold Selection
//         Bool_t passJetSelection = kFALSE;
//         if (leadingJetPt > 0) {
//           passJetSelection = kTRUE;               
//         }
//         if (!passJetSelection) continue;

              
        //use V4 denominator, combined trigger sample
        if (passMuonDenominatorCuts(info->triggerBits, mu, info->PileupEnergyDensity)) {
              
          //Fill These Muons

          fWeight = 1.0;
          fRunNumber = info->runNum;
          fLumiSectionNumber = info->lumiSec;
          fEventNumber = info->evtNum;
          fMuPt = mu->pt; 
          fMuEta = mu->eta; 
          fMuPhi = mu->phi; 
          fMuPFIso = mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,mu->eta);
   
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
          fMuHadEnergy = mu->HadEnergy; 
          fMuHoEnergy = mu->HoEnergy; 
          fMuEmEnergy = mu->EmEnergy; 
          fMuHadS9Energy = mu->HadS9Energy; 
          fMuHoS9Energy = mu->HoS9Energy; 
          fMuEmS9Energy = mu->EmS9Energy; 

          //Isolation Variables
          fMuChargedIso03 = mu->ChargedIso03 ; 
          fMuChargedIso03FromOtherVertices = mu->ChargedIso03FromOtherVertices ; 
          fMuNeutralIso03_05Threshold = mu->NeutralIso03_05Threshold ; 
          fMuNeutralIso03_10Threshold = mu->NeutralIso03_10Threshold ; 
          fMuChargedIso04 = mu->ChargedIso04 ; 
          fMuChargedIso04FromOtherVertices = mu->ChargedIso04FromOtherVertices ; 
          fMuNeutralIso04_05Threshold = mu->NeutralIso04_05Threshold ; 
          fMuNeutralIso04_10Threshold = mu->NeutralIso04_10Threshold ; 
          fMuTrkIso03 = mu->trkIso03; 
          fMuEMIso03 = mu->emIso03; 
          fMuHadIso03 = mu->hadIso03; 
          fMuTrkIso05 = mu->trkIso05; 
          fMuEMIso05 = mu->emIso05; 
          fMuHadIso05 = mu->hadIso05; 
          fRho = info->PileupEnergyDensity; 
          fNVertices = info->nPV0; 

          NMuonsFilled++;
          muTree->Fill();

        }

      } //loop over muons

    } //end loop over data  

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


Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  if (!( 
        ((triggerBits & kHLT_Mu15)
         && (mu->hltMatchBits & kHLTObject_Mu15)
          )
        ||
        ((triggerBits & kHLT_Mu8)
         && (mu->hltMatchBits & kHLTObject_Mu8)
          )
        ||
        ((triggerBits & kHLT_Mu8_Jet40)
         && (mu->hltMatchBits & kHLTObject_Mu8)
          )
        )
    ) {
    pass = kFALSE;
  }
  
  pass = pass && passMuonDenominatorCuts(mu, fRho);

  return pass;
}

Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

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
        && mu->TrkKink  < 20
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
