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
Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele);
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


//=== MAIN MACRO =================================================================================================
void MakeFakeElectronTrainingNtuple() {

  MakeNtuple("LIST","ElectronSelectionTraining.Fakes.root");

}


void MakeNtuple(const string inputFilename, const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 26.5;

  //*****************************************************************************************
  //Setup Likelihood
  //*****************************************************************************************
  TFile *fileLH = TFile::Open("/home/sixie/CMSSW_analysis/src/MitPhysics/data/ElectronLikelihoodPdfs_MC.root");
  TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
  TDirectory *EElt15dir = fileLH->GetDirectory("/");
  TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
  TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
  TDirectory *EEgt15dir = fileLH->GetDirectory("/");

  LikelihoodSwitches defaultSwitches;
  defaultSwitches.m_useEoverP = false;
  defaultSwitches.m_useDeltaEta = true;
  defaultSwitches.m_useDeltaPhi = true;
  defaultSwitches.m_useHoverE = false;        
  defaultSwitches.m_useSigmaEtaEta = true;
  defaultSwitches.m_useSigmaPhiPhi = true;
  defaultSwitches.m_useFBrem = true;
  defaultSwitches.m_useOneOverEMinusOneOverP = true;
 
  ElectronLikelihood *LH = new ElectronLikelihood(&(*EB0lt15dir),&(*EB1lt15dir), &(*EElt15dir), 
                                                  &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir),
                                                  defaultSwitches,
                                                  std::string("class"),std::string("class"),true,true);

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************
  TFile *outputFile = new TFile(outputFilename.c_str(), "RECREATE");
  TTree *eleTree = new TTree("Electrons","Electrons");
  eleTree->SetAutoFlush(0);

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

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleCovIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;



  //Isolation Variables
  Float_t                 fEleStandardLikelihood;
  Float_t                 fElePFMVA;
  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 



  eleTree->Branch("weight",&fWeight,"weight/F");
  eleTree->Branch("run",&fRunNumber,"run/i");
  eleTree->Branch("lumi",&fLumiSectionNumber,"lumi/i");
  eleTree->Branch("event",&fEventNumber,"event/i");
  eleTree->Branch("pt",&fElePt,"pt/F"); 
  eleTree->Branch("eta",&fEleEta,"eta/F"); 
  eleTree->Branch("phi",&fElePhi,"phi/F"); 
  eleTree->Branch("scet",&fEleSCEt,"scet/F"); 
  eleTree->Branch("sceta",&fEleSCEta,"sceta/F"); 
  eleTree->Branch("scphi",&fEleSCPhi,"scphi/F"); 
  eleTree->Branch("pfiso",&fElePFIso,"pfiso/F"); 
  
  //CutBased Variables
  eleTree->Branch("SigmaIEtaIEta",&fEleSigmaIEtaIEta,"SigmaIEtaIEta/F"); 
  eleTree->Branch("DEtaIn",&fEleDEtaIn,"DEtaIn/F"); 
  eleTree->Branch("DPhiIn",&fEleDPhiIn,"DPhiIn/F"); 
  eleTree->Branch("HoverE",&fEleHoverE,"HoverE/F"); 
  eleTree->Branch("D0",&fEleD0,"D0/F"); 
  eleTree->Branch("DZ",&fEleDZ,"DZ/F"); 
  eleTree->Branch("FBrem",&fEleFBrem,"FBrem/F"); 
  eleTree->Branch("EOverP",&fEleEOverP,"EOverP/F"); 

  //Additional Vars used in Likelihood
  eleTree->Branch("ESeedClusterOverPout",&fEleESeedClusterOverPout,"ESeedClusterOverPout/F"); 
  eleTree->Branch("SigmaIPhiIPhi",&fEleSigmaIPhiIPhi,"SigmaIPhiIPhi/F"); 
  eleTree->Branch("NBrem",&fEleNBrem,"NBrem/F"); 
  eleTree->Branch("OneOverEMinusOneOverP",&fEleOneOverEMinusOneOverP,"OneOverEMinusOneOverP/F"); 
  eleTree->Branch("ESeedClusterOverPIn",&fEleESeedClusterOverPIn,"ESeedClusterOverPIn/F"); 
  eleTree->Branch("IP3d",&fEleIP3d,"IP3d/F"); 
  eleTree->Branch("IP3dSig",&fEleIP3dSig,"IP3dSig/F"); 

  eleTree->Branch("HcalDepth1OverEcal",&fEleHcalDepth1OverEcal,"HcalDepth1OverEcal/F"); 
  eleTree->Branch("HcalDepth2OverEcal",&fEleHcalDepth2OverEcal,"HcalDepth2OverEcal/F"); 
  eleTree->Branch("dEtaCalo",&fEledEtaCalo,"dEtaCalo/F"); 
  eleTree->Branch("dPhiCalo",&fEledPhiCalo,"dPhiCalo/F"); 
  eleTree->Branch("PreShowerOverRaw",&fElePreShowerOverRaw,"PreShowerOverRaw/F"); 
  eleTree->Branch("CovIEtaIPhi",&fEleCovIEtaIPhi,"CovIEtaIPhi/F"); 
  eleTree->Branch("SCEtaWidth",&fEleSCEtaWidth,"SCEtaWidth/F"); 
  eleTree->Branch("SCPhiWidth",&fEleSCPhiWidth,"SCPhiWidth/F"); 
  eleTree->Branch("GsfTrackChi2OverNdof",&fEleGsfTrackChi2OverNdof,"GsfTrackChi2OverNdof/F"); 
  eleTree->Branch("R9",&fEleR9,"R9/F"); 

  eleTree->Branch("SeedEMaxOverE",&fEleSeedEMaxOverE,"SeedEMaxOverE/F"); 
  eleTree->Branch("SeedETopOverE",&fEleSeedETopOverE,"SeedETopOverE/F"); 
  eleTree->Branch("SeedEBottomOverE",&fEleSeedEBottomOverE,"SeedEBottomOverE/F"); 
  eleTree->Branch("SeedELeftOverE",&fEleSeedELeftOverE,"SeedELeftOverE/F"); 
  eleTree->Branch("SeedERightOverE",&fEleSeedERightOverE,"SeedERightOverE/F"); 
  eleTree->Branch("SeedE2ndOverE",&fEleSeedE2ndOverE,"SeedE2ndOverE/F"); 
  eleTree->Branch("SeedE2x5RightOverE",&fEleSeedE2x5RightOverE,"SeedE2x5RightOverE/F"); 
  eleTree->Branch("SeedE2x5LeftOverE",&fEleSeedE2x5LeftOverE,"SeedE2x5LeftOverE/F"); 
  eleTree->Branch("SeedE2x5TopOverE",&fEleSeedE2x5TopOverE,"SeedE2x5TopOverE/F"); 
  eleTree->Branch("SeedE2x5BottomOverE",&fEleSeedE2x5BottomOverE,"SeedE2x5BottomOverE/F"); 
  eleTree->Branch("SeedE2x5MaxOverE",&fEleSeedE2x5MaxOverE,"SeedE2x5MaxOverE/F"); 
  eleTree->Branch("SeedE1x3OverE",&fEleSeedE1x3OverE,"SeedE1x3OverE/F"); 
  eleTree->Branch("SeedE3x1OverE",&fEleSeedE3x1OverE,"SeedE3x1OverE/F"); 
  eleTree->Branch("SeedE1x5OverE",&fEleSeedE1x5OverE,"SeedE1x5OverE/F"); 
  eleTree->Branch("SeedE2x2OverE",&fEleSeedE2x2OverE,"SeedE2x2OverE/F"); 
  eleTree->Branch("SeedE3x2OverE",&fEleSeedE3x2OverE,"SeedE3x2OverE/F"); 
  eleTree->Branch("SeedE3x3OverE",&fEleSeedE3x3OverE,"SeedE3x3OverE/F"); 
  eleTree->Branch("SeedE4x4OverE",&fEleSeedE4x4OverE,"SeedE4x4OverE/F"); 
  eleTree->Branch("SeedE5x5OverE",&fEleSeedE5x5OverE,"SeedE5x5OverE/F"); 

  //Isolation Variables
  eleTree->Branch("StandardLikelihood",&fEleStandardLikelihood,"StandardLikelihood/F"); 
  eleTree->Branch("PFMVA",&fElePFMVA,"PFMVA/F"); 
  eleTree->Branch("ChargedIso03",&fEleChargedIso03,"ChargedIso03/F"); 
  eleTree->Branch("NeutralHadronIso03",&fEleNeutralHadronIso03,"NeutralHadronIso03/F"); 
  eleTree->Branch("GammaIso03",&fEleGammaIso03,"GammaIso03/F"); 
  eleTree->Branch("ChargedIso04",&fEleChargedIso04,"ChargedIso04/F"); 
  eleTree->Branch("NeutralHadronIso04",&fEleNeutralHadronIso04,"NeutralHadronIso04/F"); 
  eleTree->Branch("GammaIso04",&fEleGammaIso04,"GammaIso04/F"); 
  eleTree->Branch("ChargedIso04FromOtherVertices",&fEleChargedIso04FromOtherVertices,"ChargedIso04FromOtherVertices/F"); 
  eleTree->Branch("NeutralHadronIso04_10Threshold",&fEleNeutralHadronIso04_10Threshold,"NeutralHadronIso04_10Threshold/F"); 
  eleTree->Branch("GammaIso04_10Threshold",&fEleGammaIso04_10Threshold,"GammaIso04_10Threshold/F"); 
  eleTree->Branch("TrkIso03",&fEleTrkIso03,"TrkIso03/F"); 
  eleTree->Branch("EMIso03",&fEleEMIso03,"EMIso03/F"); 
  eleTree->Branch("HadIso03",&fEleHadIso03,"HadIso03/F"); 
  eleTree->Branch("TrkIso04",&fEleTrkIso04,"TrkIso04/F"); 
  eleTree->Branch("EMIso04",&fEleEMIso04,"EMIso04/F"); 
  eleTree->Branch("HadIso04",&fEleHadIso04,"HadIso04/F"); 
  eleTree->Branch("Rho",&fRho,"Rho/F"); 
  eleTree->Branch("NVertices",&fNVertices,"NVertices/F"); 

  UInt_t NElectronsFilled = 0;
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
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-m10-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-pr-v4_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-a05-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11a-del-o03-v1_EleFakeRateTriggerAndDenominatorSkim.root");
    inputfiles.push_back("/home/sixie/hist/AllNtuple/HWWNtuple/data/AllNtuple_HWWNtuple_r11b-del-pr-v1_EleFakeRateTriggerAndDenominatorSkim.root");
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

      Int_t NElectrons = electronArr->GetEntries();
      
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
 
        //find leading jet in the event
        Double_t leadingJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
        }

        //********************************************************
        // Event Selection Cuts
        //********************************************************
              
        //veto events with more than 1 reco electron
        if (NElectrons > 1) continue;

        //met cut removed W events
        if (met > 20) continue;

//         //Jet Threshold Selection
//         Bool_t passJetSelection = kFALSE;
//         if (leadingJetPt > 0) {
//           passJetSelection = kTRUE;               
//         }
//         if (!passJetSelection) continue;
              
        //protect against pathologies
        if (TMath::IsNaN(ele->sigiPhiiPhi)) {
          cout << "Pathological SigmaIPhiIPhi : " 
               << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
          continue;
        }
       
        //use V4 denominator, combined trigger sample
        if (passElectronDenominatorCuts(info->triggerBits, ele)) {
              
          //Fill These Electrons

          fWeight = 1.0;
          fRunNumber = info->runNum;
          fLumiSectionNumber = info->lumiSec;
          fEventNumber = info->evtNum;
          fElePt = ele->pt; 
          fEleEta = ele->eta; 
          fElePhi = ele->phi; 
          fEleSCEt = ele->scEt; 
          fEleSCEta = ele->scEta; 
          fEleSCPhi = ele->scPhi; 
          fElePFIso = ( ele->ChargedIso04 + ele->NeutralHadronIso04_10Threshold
            + ele->GammaIso04_10Threshold - ele->GammaIsoVetoEtaStrip04_10Threshold
                        - ele->ChargedEMIsoVetoEtaStrip04 - ele->NeutralHadronIso007_10Threshold) / ele->pt; 
  
          fEleSigmaIEtaIEta = ele->sigiEtaiEta; 
          fEleDEtaIn = ele->deltaEtaIn; 
          fEleDPhiIn = ele->deltaPhiIn; 
          fEleHoverE = ele->HoverE; 
          fEleD0 = ele->d0; 
          fEleDZ = ele->dz; 
          fEleFBrem = ele->fBrem; 
          fEleEOverP = ele->EOverP; 

          fEleESeedClusterOverPout = ele->ESeedClusterOverPout; 
          fEleSigmaIPhiIPhi = TMath::Sqrt(ele->sigiPhiiPhi); 
          fEleNBrem = ele->nBrem; 
          fEleOneOverEMinusOneOverP = (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p);

          fEleESeedClusterOverPIn = ele->ESeedClusterOverPIn; 
          fEleIP3d = ele->ip3d; 
          fEleIP3dSig = ele->ip3dSig; 

          fEleHcalDepth1OverEcal = ele->HcalDepth1OverEcal;
          fEleHcalDepth2OverEcal = ele->HcalDepth2OverEcal;
          fEledEtaCalo = ele->dEtaCalo;
          fEledPhiCalo = ele->dPhiCalo;
          fElePreShowerOverRaw = ele->PreShowerOverRaw;
          fEleCovIEtaIPhi = ele->CovIEtaIPhi;
          fEleSCEtaWidth = ele->SCEtaWidth;
          fEleSCPhiWidth = ele->SCPhiWidth;
          fEleGsfTrackChi2OverNdof = ele->GsfTrackChi2OverNdof;
          fEleR9 = ele->R9;
          
          fEleSeedEMaxOverE = ele->SeedEMaxOverE;
          fEleSeedETopOverE = ele->SeedETopOverE;
          fEleSeedEBottomOverE = ele->SeedEBottomOverE;
          fEleSeedELeftOverE = ele->SeedELeftOverE;
          fEleSeedERightOverE = ele->SeedERightOverE;
          fEleSeedE2ndOverE = ele->SeedE2ndOverE;
          fEleSeedE2x5RightOverE = ele->SeedE2x5RightOverE;
          fEleSeedE2x5LeftOverE = ele->SeedE2x5LeftOverE;
          fEleSeedE2x5TopOverE = ele->SeedE2x5TopOverE;
          fEleSeedE2x5BottomOverE = ele->SeedE2x5BottomOverE;
          fEleSeedE2x5MaxOverE = ele->SeedE2x5MaxOverE;
          fEleSeedE1x3OverE = ele->SeedE1x3OverE;
          fEleSeedE3x1OverE = ele->SeedE3x1OverE;
          fEleSeedE1x5OverE = ele->SeedE1x5OverE;
          fEleSeedE2x2OverE = ele->SeedE2x2OverE;
          fEleSeedE3x2OverE = ele->SeedE3x2OverE;
          fEleSeedE3x3OverE = ele->SeedE3x3OverE;
          fEleSeedE4x4OverE = ele->SeedE4x4OverE;
          fEleSeedE5x5OverE = ele->SeedE5x5OverE;

          fElePFMVA = ele->mva; 
          fEleChargedIso03 = ele->ChargedIso03; 
          fEleNeutralHadronIso03 = ele->NeutralHadronIso03_05Threshold - ele->NeutralHadronIso007_05Threshold; 
          fEleGammaIso03 = ele->GammaIso03_05Threshold - ele->GammaIsoVetoEtaStrip03_05Threshold; 
          fEleChargedIso04 = ele->ChargedIso04; 
          fEleNeutralHadronIso04 = ele->NeutralHadronIso04_05Threshold - ele->NeutralHadronIso007_05Threshold; 
          fEleGammaIso04 = ele->GammaIso04_05Threshold - ele->GammaIsoVetoEtaStrip04_05Threshold; 
          fEleChargedIso04FromOtherVertices = ele->ChargedIso04FromOtherVertices; 
          fEleNeutralHadronIso04_10Threshold = ele->NeutralHadronIso04_10Threshold; 
          fEleGammaIso04_10Threshold = ele->GammaIso04_10Threshold; 
          fEleTrkIso03 = ele->trkIso03; 
          fEleEMIso03 = ele->emIso03; 
          fEleHadIso03 = ele->hadIso03; 
          fEleTrkIso04 = ele->trkIso04; 
          fEleEMIso04 = ele->emIso04; 
          fEleHadIso04 = ele->hadIso04; 
          fRho = info->PileupEnergyDensity; 
          fNVertices = info->nPV0; 

          //likelihood value
          mithep::FourVectorM tmpEleSC;
          tmpEleSC.SetCoordinates(ele->scEt, ele->scEta, ele->scPhi, 0.51099892e-3 );
            
          LikelihoodMeasurements measurements;
          measurements.pt = ele->pt;
          if (ele->isEB && (fabs(ele->eta)<1.0)) measurements.subdet = 0;
          else if (ele->isEB) measurements.subdet = 1;
          else measurements.subdet = 2;
          measurements.deltaPhi = ele->deltaPhiIn;
          measurements.deltaEta = ele->deltaEtaIn;
          measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout;
          measurements.eSuperClusterOverP = ele->EOverP;
          measurements.hadronicOverEm = ele->HoverE;
          measurements.sigmaIEtaIEta = ele->sigiEtaiEta;
          measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->sigiPhiiPhi);
          measurements.fBrem = ele->fBrem;
          measurements.nBremClusters = ele->nBrem;
          measurements.OneOverEMinusOneOverP = (1.0/(ele->EOverP*ele->p)) - 1.0 / ele->p;

          double likelihood = LH->result(measurements);
          fEleStandardLikelihood = 0; 
          if (likelihood <= 0) {
            fEleStandardLikelihood = -20.0;
          } else if (likelihood == 1) {
            fEleStandardLikelihood = 20.0;
          } else {
            fEleStandardLikelihood = TMath::Log(likelihood / (1.0-likelihood));
          }

          NElectronsFilled++;
          eleTree->Fill();

        }

      } //loop over electrons

    } //end loop over data  

     cout << "Total Electrons: " << NElectronsFilled << endl;

  } //end loop over files

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;


  cout << "Total Electrons: " << NElectronsFilled << endl;
  outputFile->Write();
  outputFile->Close();

  gBenchmark->Show("WWTemplate");       
} 

Bool_t passConversionVeto(Int_t isConv) {
 
  Int_t tmp0 = floor(double(isConv) / 2.0);
  Int_t tmp1 = floor(double(tmp0) / 2.0);
  Int_t tmp2 = floor(double(tmp1) / 2.0);
  Int_t tmp3 = floor(double(tmp2) / 2.0);
  Int_t tmp4 = floor(double(tmp3) / 2.0);
  Int_t tmp5 = floor(double(tmp4) / 2.0);
  Int_t tmp6 = floor(double(tmp5) / 2.0);
  Int_t tmp7 = floor(double(tmp6) / 2.0);
  Int_t tmp8 = floor(double(tmp7) / 2.0);
  Int_t tmp9 = floor(double(tmp8) / 2.0);
  Int_t tmp10 = floor(double(tmp9) / 2.0);
  Int_t tmp11 = floor(double(tmp10) / 2.0);
  Int_t tmp12 = floor(double(tmp11) / 2.0);

  Bool_t pass;
  pass =  tmp9 % 2;
  return pass; 
}

Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  if (!(
        ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
         && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL)
        ||
        ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
         && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL)
        ||
        ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
         && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL)
        ||
        ((triggerBits & kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL)
         && ele->hltMatchBits & kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL
          )
        )
    ) {
    pass = kFALSE;
  }

  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->dz) < 0.1
            && fabs(ele->d0) < 0.02
            && (ele->trkIso03) / ele->pt < 0.2
            && (ele->emIso03) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20
              
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
             && ele->HoverE < 0.10
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->dz) < 0.1
             && fabs(ele->d0) < 0.02
             && (ele->trkIso03) / ele->pt < 0.2
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20
          )
      ) {
      pass = kFALSE;
    }
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
