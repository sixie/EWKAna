//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwNtupler.C+\(\"0000\",\"noskim\",\"r11a-dmu-pr-v1_1\",\"cern/filefi/020\",\"/home/mitprod/catalog\",\"HwwNtuple\",10000,-1\) 
//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwNtupler.C+\(\"0000\",\"noskim\",\"p11-h130ww2l-gf-v1g1-pu\",\"cern/filefi/020\",\"/home/mitprod/catalog\",\"HwwNtuple\",-1,1100\) 
//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwNtupler.C+\(\"0000\",\"noskim\",\"f10-zmm-powheg-c10-v12\",\"cern/filefi/018\",\"/home/mitprod/catalog\",\"HwwNtuple\",10000,2\) 
//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwNtupler.C+\(\"0000\",\"noskim\",\"w10-wjets-z2-v8-pu11-2l\",\"cern/filler/020\",\"/home/ceballos/catalog\",\"HwwAnalysis\",-1,333\)   

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "EWKAna/Ntupler/interface/HwwNtuplerMod.hh"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"
#include "MitPhysics/Mods/interface/PartonFlavorHistoryMod.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/JetCol.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/PFMetCol.h"
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitPhysics/SelMods/interface/GenericSelMod.h"
#include "MitHiggs/EwkMods/interface/SignalFilterForJettinessMod.h"
#endif

using namespace mithep;

//==================================================================================================
/*
 * Triggers of interest
 *

  kHLT_Mu9                                   = 0x0000001,
  kHLT_Mu11                                  = 0x0000002,
  kHLT_Mu13_v1                               = 0x0000004,
  kHLT_Mu15_v1                               = 0x0000008,
  kHLT_DoubleMu5_v1                          = 0x0000010,
  kHLT_Jet15U                                = 0x0000020,
  kHLT_Jet30U                                = 0x0000040,
  kHLT_Jet50U                                = 0x0000080,
  kHLT_Photon10_L1R                          = 0x0000100,
  kHLT_Photon10_Cleaned_L1R                  = 0x0000200,
  kHLT_Photon15_L1R                          = 0x0000400,
  kHLT_Photon15_Cleaned_L1R                  = 0x0000800,
  kHLT_Photon20_L1R                          = 0x0001000,
  kHLT_Photon20_Cleaned_L1R                  = 0x0002000,
  kHLT_Photon30_Cleaned_L1R                  = 0x0004000,
  kHLT_Ele15_SW_L1R                          = 0x0008000,
  kHLT_Ele15_LW_L1R                          = 0x0010000,
  kHLT_Ele15_SW_CaloEleId_L1R                = 0x0020000,
  kHLT_Ele17_SW_L1R                          = 0x0040000,
  kHLT_Ele17_SW_CaloEleId_L1R                = 0x0080000,       
  kHLT_Ele17_SW_TightEleId_L1R               = 0x0100000,
  kHLT_Ele17_SW_TightCaloEleId_Ele8HE_L1R_v1 = 0x0200000,
  kHLT_Ele17_SW_TightCaloEleId_SC8HE_L1R_v1  = 0x0400000,
  kHLT_Ele17_SW_TighterEleIdIsol_L1R_v2      = 0x0800000,
  kHLT_DoubleEle10_SW_L1R                    = 0x1000000,
  kHLT_DoubleEle15_SW_L1R_v1                 = 0x2000000

  */
    
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -l -q -b runZeeNtupler.C+\(\"0000\",\"p10-zee-v26\",\"cern/filler/014a\",\"/home/ceballos/catalog\",1,0,1,-1,0,1\)
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runHwwNtupler(
  const char *fileset  = "",
  const char *skim         = "noskim",
  const char *dataset    = "s8-ttbar-id9",
  const char *book       = "mit/filler/006",
  const char *catalogDir = "/home/mitprod/catalog",
  const char *outputName = "HwwHiggsNtupleMaker",
  int   nEvents          = -1,
  int   sampleID         = -1
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
 
  Bool_t useGen = kFALSE;
  Bool_t applyWWFilter = kFALSE; 
  Int_t fsrmode = kFALSE;
  Bool_t useHlt = kTRUE;
  Bool_t skipHLTFail = kFALSE;
  Int_t runSkim = 0; if (sampleID == -1) runSkim = 1; if (sampleID == -2) runSkim = 2; if (sampleID == -3) runSkim = 3;
  Bool_t computePDFWeights = kFALSE; 
  string pdfSetName = "";
  if (sampleID == 10010) {
    runSkim = 0;
    skipHLTFail = kFALSE;
    computePDFWeights = kTRUE;
    pdfSetName = "cteq66.LHgrid";
  } else if (sampleID == 10011) {
    runSkim = 0;
    skipHLTFail = kFALSE;
    computePDFWeights = kTRUE;
    pdfSetName = "MSTW2008nlo68cl.LHgrid";
  } else if (sampleID == 10012) {
    runSkim = 0;
    skipHLTFail = kFALSE;
    computePDFWeights = kTRUE;
    pdfSetName = "NNPDF20_100.LHgrid";
  } else if (sampleID == 22) {
    applyWWFilter = kTRUE;
  }

  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // supercluster kinematics
  const Double_t SCEtMin  = 15;
  const Double_t SCEtMax  = 7000;
  const Double_t SCEtaMin = -3;
  const Double_t SCEtaMax =  3;
  
  // muon kinematics
  const Double_t MuonPtMin  = 3;
  const Double_t MuonPtMax  = 7000;
  const Double_t MuonEtaMin = -3;
  const Double_t MuonEtaMax =  3;

  // jet requirements
  const Double_t jetPtMin = 7;

  //for MC
  if ( sampleID >= 0 ) useGen = kTRUE; 

  //------------------------------------------------------------------------------------------------
  // generator information
  //------------------------------------------------------------------------------------------------
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(15.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(10.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(!useGen);
  generatorMod->SetFillHist(useGen);
  generatorMod->SetApplyISRFilter(kFALSE);
  generatorMod->SetApplyWWFilter(applyWWFilter);

  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  HLTMod *hltmod = new HLTMod;

  hltmod->AddTrigger("HLT_Mu24_v1",0,999999);
  hltmod->AddTrigger("HLT_DoubleMu7_v1",0,999999);
  hltmod->SetAbortIfNotAccepted(kFALSE);
  hltmod->SetTrigObjsName("myhltobjs");


  //------------------------------------------------------------------------------------------------
  // PV filter selection
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);
  if (sampleID == 333)   goodPVFilterMod->SetVertexesName("DAPrimaryVertexes");
  else goodPVFilterMod->SetVertexesName("PrimaryVertexes");

 //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------

  MuonIDMod *muTightId = new MuonIDMod;  
  muTightId->SetClassType ("GlobalTracker");
  muTightId->SetIDType    ("WWMuIdV2");
  muTightId->SetIsoType   ("PFIso");
  muTightId->SetApplyD0Cut(kTRUE);
  muTightId->SetApplyDZCut(kTRUE);
  muTightId->SetWhichVertex(0);
  muTightId->SetOutputName("SkimMuons");

  ElectronIDMod       *electronTightId       = new ElectronIDMod;
  electronTightId->SetIDType(TString("VBTFWorkingPointLowPtId"));
  electronTightId->SetIsoType(TString("PFIso"));
  electronTightId->SetApplyConversionFilterType1(kTRUE);
  electronTightId->SetApplyConversionFilterType2(kFALSE);
  electronTightId->SetChargeFilter(kFALSE);
  electronTightId->SetApplyD0Cut(kTRUE);
  electronTightId->SetApplyDZCut(kTRUE);
  electronTightId->SetWhichVertex(0);
  electronTightId->SetNExpectedHitsInnerCut(0);
  electronTightId->SetOutputName("SkimElectrons");




  MuonIDMod *muDenominator = new MuonIDMod;  
  muDenominator->SetClassType ("GlobalTracker");
  muDenominator->SetIDType    ("WWMuIdV2");
  muDenominator->SetIsoType   ("PFIso");
  muDenominator->SetApplyD0Cut(kTRUE);
  muDenominator->SetApplyDZCut(kTRUE);
  muDenominator->SetD0Cut(0.20);
  muDenominator->SetDZCut(0.10);
  muDenominator->SetCombIsoCut(1.00); 
  muDenominator->SetOutputName("SkimDenominatorMuons");
  muDenominator->SetWhichVertex(0);

  ElectronIDMod *electronDenominator = new ElectronIDMod;
  electronDenominator->SetIDType("VBTFWorkingPointFakeableId");
  electronDenominator->SetIsoType("TrackJura");
  electronDenominator->SetTrackIsoCut(0.2);
  electronDenominator->SetEcalJurIsoCut(0.2);
  electronDenominator->SetHcalIsoCut(0.2);
  electronDenominator->SetApplyConversionFilterType1(kTRUE);
  electronDenominator->SetApplyConversionFilterType2(kFALSE);
  electronDenominator->SetChargeFilter              (kFALSE);
  electronDenominator->SetApplyD0Cut                (kTRUE);
  electronDenominator->SetApplyDZCut                (kTRUE);
  electronDenominator->SetNExpectedHitsInnerCut     (0);
  electronDenominator->SetD0Cut(0.02);
  electronDenominator->SetDZCut(0.10);
  electronDenominator->SetOutputName("SkimDenominatorElectrons");


  MergeLeptonsMod *mergedTight = new MergeLeptonsMod;
  mergedTight->SetMuonsName    (muTightId->GetOutputName());
  mergedTight->SetElectronsName(electronTightId->GetOutputName());
  mergedTight->SetOutputName("mergedTightLeptons");

  MergeLeptonsMod *mergedLoose = new MergeLeptonsMod;
  mergedLoose->SetMuonsName    (muDenominator->GetOutputName());
  mergedLoose->SetElectronsName(electronDenominator->GetOutputName());
  mergedLoose->SetOutputName("mergedLooseLeptons");

  GenericSelMod<mithep::Particle> *selModTight = new GenericSelMod<mithep::Particle>;
  selModTight->SetPtMin(0.0);
  selModTight->SetMinCounts(1);
  selModTight->SetColName(mergedTight->GetOutputName());

  GenericSelMod<mithep::Particle> *selModLoose = new GenericSelMod<mithep::Particle>;
  selModLoose->SetPtMin(0.0);
  selModLoose->SetMinCounts(1);
  selModLoose->SetColName(mergedLoose->GetOutputName());

  GenericSelMod<mithep::Particle> *selModDoubleLoose = new GenericSelMod<mithep::Particle>;
  selModDoubleLoose->SetPtMin(0.0);
  selModDoubleLoose->SetMinCounts(2);
  selModDoubleLoose->SetColName(mergedLoose->GetOutputName());



  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("PubPFMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromFile("MitPhysics/data/START41_V0_AK5PF_L1FastJet.txt"); 
  jetCorr->AddCorrectionFromFile("MitPhysics/data/START41_V0_AK5PF_L2Relative.txt"); 
  jetCorr->AddCorrectionFromFile("MitPhysics/data/START41_V0_AK5PF_L3Absolute.txt");
  if(!useGen){ 
    jetCorr->AddCorrectionFromFile("MitPhysics/data/START41_V0_AK5PF_L2L3Residual.txt");
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod *muonID = new MuonIDMod;
  muonID->SetClassType("GlobalTracker");
  muonID->SetIDType("WWMuIdV2");
  muonID->SetIsoType("PFIso");
  muonID->SetApplyD0Cut(kTRUE);
  muonID->SetApplyDZCut(kTRUE);
  muonID->SetWhichVertex(0);

//   TFile *fileLH = TFile::Open("MitPhysics/data/ElectronLikelihoodPdfs_MC.root");
//   TDirectory *EB0lt15dir = fileLH->GetDirectory("/");
//   TDirectory *EB1lt15dir = fileLH->GetDirectory("/");
//   TDirectory *EElt15dir = fileLH->GetDirectory("/");
//   TDirectory *EB0gt15dir = fileLH->GetDirectory("/");
//   TDirectory *EB1gt15dir = fileLH->GetDirectory("/");
//   TDirectory *EEgt15dir = fileLH->GetDirectory("/");
//   //fileLH->Close();
//   //delete fileLH;
//   LikelihoodSwitches defaultSwitches;
//   defaultSwitches.m_useFBrem = true;
//   defaultSwitches.m_useEoverP = false;
//   defaultSwitches.m_useOneOverEMinusOneOverP = true;
//   defaultSwitches.m_useSigmaPhiPhi = true;
//   defaultSwitches.m_useHoverE = false;        
//   ElectronLikelihood *LH = new ElectronLikelihood(&(*EB0lt15dir),&(*EB1lt15dir), &(*EElt15dir), &(*EB0gt15dir), &(*EB1gt15dir), &(*EEgt15dir),
//                               defaultSwitches, std::string("class"),std::string("class"),true,true);

//   ElectronIDMod *electronID = new ElectronIDMod;
//   electronID->SetIDLikelihoodCut(1.0);
//   electronID->SetLH(LH);
//   electronID->SetIDType("Likelihood");
//   //electronID->SetIDType("VBTFWorkingPointLowPtId");
//   electronID->SetIsoType("PFIso");
//   electronID->SetApplyConversionFilterType1(kTRUE);
//   electronID->SetApplyConversionFilterType2(kFALSE);
//   electronID->SetChargeFilter(kFALSE);
//   electronID->SetApplyD0Cut(kTRUE);
//   electronID->SetApplyDZCut(kTRUE);
//   electronID->SetWhichVertex(0);
//   electronID->SetNExpectedHitsInnerCut(0);

  ElectronIDMod *electronID = new ElectronIDMod;
  electronID->SetIDType("VBTFWorkingPointLowPtId");
  electronID->SetIsoType("PFIso");
  electronID->SetApplyConversionFilterType1(kTRUE);
  electronID->SetApplyConversionFilterType2(kFALSE);
  electronID->SetChargeFilter(kFALSE);
  electronID->SetApplyD0Cut(kTRUE);
  electronID->SetApplyDZCut(kTRUE);
  electronID->SetWhichVertex(0);
  electronID->SetNExpectedHitsInnerCut(0);

 // Object ID and Cleaning Sequence
  PhotonIDMod         *photonID      = new PhotonIDMod;
  TauIDMod            *tauID         = new TauIDMod;
  JetIDMod            *jetID         = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(30.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetJetEEMFractionMinCut(0.0);
  jetID->SetOutputName("GoodJets");
  jetID->SetApplyBetaCut(kFALSE);

  ElectronCleaningMod *electronCleaning = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning   = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning      = new TauCleaningMod;
  JetCleaningMod      *jetCleaning      = new JetCleaningMod;
  jetCleaning->SetGoodJetsName("GoodJets");
  jetCleaning->SetCleanJetsName("CleanJets");

  JetIDMod            *jetIDNoPtCut     = new JetIDMod;
  jetIDNoPtCut->SetInputName(jetCorr->GetOutputName());
  jetIDNoPtCut->SetPtCut(0.0);
  jetIDNoPtCut->SetEtaMaxCut(5.0);
  jetIDNoPtCut->SetJetEEMFractionMinCut(0.0);
  jetIDNoPtCut->SetOutputName("GoodJetsNoPtCut");
  jetIDNoPtCut->SetApplyBetaCut(kFALSE);


  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *leptonMerger = new MergeLeptonsMod;
  leptonMerger->SetMuonsName(muonID->GetOutputName());
  leptonMerger->SetElectronsName(electronCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // Produce candidates for Jettiness calculation
  //------------------------------------------------------------------------------------------------
  SignalFilterForJettinessMod * signalFilter = new SignalFilterForJettinessMod;
  signalFilter->SetMuonsCollectionName(muonID->GetOutputName());
  signalFilter->SetElectronsCollectionName(electronCleaning->GetOutputName());
  signalFilter->SetLeptonCollectionName(leptonMerger->GetOutputName());
  signalFilter->SetPFJetCollectionName(jetCleaning->GetOutputName());
  signalFilter->SetMetCollectionName(pubPFMet->GetOutputName());
  signalFilter->SetOutputJetsPtCut(25);
  signalFilter->SetOutputTrackPtCut(0.5);
  signalFilter->SetOutputPFCandidatePtCut(0.5);



  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString("/data/blue/sixie/temp/") + TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  printf("\nRoot output: %s\n\n",rootFile.Data());  

  //------------------------------------------------------------------------------------------------
  //
  // setup ntupler module
  //
  //------------------------------------------------------------------------------------------------
  HwwNtuplerMod *mymod = new HwwNtuplerMod;
  mymod->SetOutputName((rootFile+TString("_copy")).Data());          // output ntuple file name
  mymod->SetUseGen(useGen);              // look at generator information (must set to kFALSE if MCParticle collection do not exist)
  mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  mymod->SetFSRMode(fsrmode);
  mymod->SetMuonPtMin(MuonPtMin);
  mymod->SetMuonPtMax(MuonPtMax);
  mymod->SetMuonEtaMin(MuonEtaMin);
  mymod->SetMuonEtaMax(MuonEtaMax);
  mymod->SetCleanJetsName(jetCleaning->GetOutputName());
  mymod->SetCleanJetsNoPtCutName(jetIDNoPtCut->GetOutputName());
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetComputePDFWeights(computePDFWeights);
  mymod->SetPDFName(pdfSetName.c_str());
  if (sampleID > 20000) {
    mymod->SetConversionName("MvfConversions");
    mymod->SetReadPileupInfo(kFALSE);
  }

 //Single Lepton Triggers
  mymod->AddTrigger("HLT_Mu24_v1",                                   kHLT_Mu24, kHLTObject_Mu24 , "hltSingleMu24L3Filtered24" );
  mymod->AddTrigger("HLT_Mu24_v2",                                   kHLT_Mu24, kHLTObject_Mu24 , "hltSingleMu24L3Filtered24" );
  mymod->AddTrigger("HLT_Mu24_v3",                                   kHLT_Mu24, kHLTObject_Mu24 , "hltSingleMu24L3Filtered24" );
  mymod->AddTrigger("HLT_Mu30_v1",                                   kHLT_Mu30, kHLTObject_Mu30 , "hltSingleMu30L3Filtered30" );
  mymod->AddTrigger("HLT_Mu30_v2",                                   kHLT_Mu30, kHLTObject_Mu30 , "hltSingleMu30L3Filtered30" );
  mymod->AddTrigger("HLT_Mu30_v3",                                   kHLT_Mu30, kHLTObject_Mu30 , "hltSingleMu30L3Filtered30" );
  mymod->AddTrigger("HLT_IsoMu12_v1",                                kHLT_IsoMu12, kHLTObject_IsoMu12, "hltSingleMuIsoL3IsoFiltered12" );
  mymod->AddTrigger("HLT_IsoMu12_v2",                                kHLT_IsoMu12, kHLTObject_IsoMu12, "hltSingleMuIsoL3IsoFiltered12" );
  mymod->AddTrigger("HLT_IsoMu12_v3",                                kHLT_IsoMu12, kHLTObject_IsoMu12, "hltSingleMuIsoL3IsoFiltered12" );
  mymod->AddTrigger("HLT_IsoMu12_v4",                                kHLT_IsoMu12, kHLTObject_IsoMu12, "hltSingleMuIsoL3IsoFiltered12" );
  mymod->AddTrigger("HLT_IsoMu12_v5",                                kHLT_IsoMu12, kHLTObject_IsoMu12, "hltSingleMuIsoL3IsoFiltered12" );
  mymod->AddTrigger("HLT_IsoMu17_v5",                                kHLT_IsoMu17, kHLTObject_IsoMu17, "hltSingleMuIsoL3IsoFiltered17" );
  mymod->AddTrigger("HLT_IsoMu17_v6",                                kHLT_IsoMu17, kHLTObject_IsoMu17, "hltSingleMuIsoL3IsoFiltered17" );
  mymod->AddTrigger("HLT_IsoMu17_v7",                                kHLT_IsoMu17, kHLTObject_IsoMu17, "hltSingleMuIsoL3IsoFiltered17" );
  mymod->AddTrigger("HLT_IsoMu17_v8",                                kHLT_IsoMu17, kHLTObject_IsoMu17, "hltSingleMuIsoL3IsoFiltered17" );
  mymod->AddTrigger("HLT_IsoMu17_v9",                                kHLT_IsoMu17, kHLTObject_IsoMu17, "hltSingleMuIsoL3IsoFiltered17" );
  mymod->AddTrigger("HLT_IsoMu24_v1",                                kHLT_IsoMu24, kHLTObject_IsoMu24, "hltSingleMuIsoL3IsoFiltered24" );
  mymod->AddTrigger("HLT_IsoMu24_v2",                                kHLT_IsoMu24, kHLTObject_IsoMu24, "hltSingleMuIsoL3IsoFiltered24" );
  mymod->AddTrigger("HLT_IsoMu24_v3",                                kHLT_IsoMu24, kHLTObject_IsoMu24, "hltSingleMuIsoL3IsoFiltered24" );
  mymod->AddTrigger("HLT_IsoMu24_v4",                                kHLT_IsoMu24, kHLTObject_IsoMu24, "hltSingleMuIsoL3IsoFiltered24" );
  mymod->AddTrigger("HLT_IsoMu24_v5",                                kHLT_IsoMu24, kHLTObject_IsoMu24, "hltSingleMuIsoL3IsoFiltered24" );
  mymod->AddTrigger("HLT_IsoMu30_v1",                                kHLT_IsoMu30, kHLTObject_IsoMu30, "hltSingleMuIsoL3IsoFiltered30" );
  mymod->AddTrigger("HLT_IsoMu30_v2",                                kHLT_IsoMu30, kHLTObject_IsoMu30, "hltSingleMuIsoL3IsoFiltered30" );
  mymod->AddTrigger("HLT_IsoMu30_v3",                                kHLT_IsoMu30, kHLTObject_IsoMu30, "hltSingleMuIsoL3IsoFiltered30" );
  mymod->AddTrigger("HLT_IsoMu30_v4",                                kHLT_IsoMu30, kHLTObject_IsoMu30, "hltSingleMuIsoL3IsoFiltered30" );
  mymod->AddTrigger("HLT_IsoMu30_v5",                                kHLT_IsoMu30, kHLTObject_IsoMu30, "hltSingleMuIsoL3IsoFiltered30" );



  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4", kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, kHLTObject_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter");


  //Main Dielectron Triggers
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1", kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter" );
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2", kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"); 
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3", kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"); 
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4", kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"); 
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5", kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter"); 
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v1", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter" );
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5", kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL , "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter");



  //Main Dimuon Triggers
  mymod->AddTrigger("HLT_DoubleMu7_v1",                          kHLT_DoubleMu7, kHLTObject_Mu7, "hltDiMuonL3PreFiltered7" );
  mymod->AddTrigger("HLT_DoubleMu7_v2",                          kHLT_DoubleMu7, kHLTObject_Mu7, "hltDiMuonL3PreFiltered7" );
  mymod->AddTrigger("HLT_DoubleMu7_v3",                          kHLT_DoubleMu7, kHLTObject_Mu7, "hltDiMuonL3PreFiltered7" );
  mymod->AddTrigger("HLT_Mu13_Mu8_v1",                           kHLT_Mu17_Mu8,  kHLTObject_Mu13, "hltSingleMu13L3Filtered13", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  mymod->AddTrigger("HLT_Mu13_Mu8_v2",                           kHLT_Mu17_Mu8, kHLTObject_Mu13, "hltSingleMu13L3Filtered13", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  mymod->AddTrigger("HLT_Mu17_Mu8_v1",                           kHLT_Mu17_Mu8, kHLTObject_Mu17, "hltSingleMu13L3Filtered17", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  mymod->AddTrigger("HLT_Mu17_Mu8_v2",                           kHLT_Mu17_Mu8, kHLTObject_Mu17, "hltSingleMu13L3Filtered17", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  //Note: This guy was put in with new menu in run 166346, but that menu was immediately switched off due to problems.
  mymod->AddTrigger("HLT_Mu13_Mu8_v3",                           kHLT_Mu17_Mu8, kHLTObject_Mu13, "hltSingleMu13L3Filtered13", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  mymod->AddTrigger("HLT_Mu17_Mu8_v3",                           kHLT_Mu17_Mu8, kHLTObject_Mu17, "hltSingleMu13L3Filtered17", kHLTObject_Mu8, "hltDiMuonL3PreFiltered8");
  


  //Main EMu Triggers
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1Mu3EG5L3Filtered17", kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter"  );
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1Mu3EG5L3Filtered17", kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1MuOpenEG5L3Filtered17", kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1MuOpenEG5L3Filtered17", kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v5",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1MuOpenEG5L3Filtered17", kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1Mu3EG5L3Filtered8",kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter" );
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1Mu3EG5L3Filtered8",kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1MuOpenEG5L3Filtered8",kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1MuOpenEG5L3Filtered8",kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v5",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1MuOpenEG5L3Filtered8",kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");



  //T&P triggers
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v1",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdL_CaloIsoVL, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLTObject_SC17, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter"   );
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdL_CaloIsoVL, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLTObject_SC17, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter" );
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdL_CaloIsoVL, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLTObject_SC17, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter" );
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1",     kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PixelMatchFilter", kHLTObject_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter" );
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v2",     kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, kHLTObject_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PixelMatchFilter", kHLTObject_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter" );
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLTObject_SC8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter" );
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLTObject_SC8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLTObject_SC8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLTObject_SC8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v5", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLTObject_SC8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v2", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLTObject_Ele8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v3", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, kHLTObject_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT , "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLTObject_Ele8, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter");



  //Other triggers
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_v1",    kHLT_Photon26_CaloIdL_IsoVL_Photon18, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLTrackIsoFilter");
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_v2",    kHLT_Photon26_CaloIdL_IsoVL_Photon18, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLTrackIsoFilter");
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_v3",    kHLT_Photon26_CaloIdL_IsoVL_Photon18, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLTrackIsoFilter");
//   mymod->AddTrigger("HLT_Photon26_IsoVL_Photon18_v1",            kHLT_Photon26_IsoVL_Photon18);
//   mymod->AddTrigger("HLT_Photon26_IsoVL_Photon18_IsoVL_v1",                 kHLT_Photon26_IsoVL_Photon18_IsoVL );
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v1",  kHLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLHcalIsoLastFilter");
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v2",  kHLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLHcalIsoLastFilter");
//   mymod->AddTrigger("HLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL_v3",  kHLT_Photon26_CaloIdL_IsoVL_Photon18_CaloIdL_IsoVL, kHLTObject_Photon26_CaloIdL_IsoVL, "hltEG26CaloIdLIsoVLHcalIsoLastFilter");
//   mymod->AddTrigger("HLT_TripleEle10_CaloIdL_TrkIdVL_v1",                    kHLT_TripleEle10_CaloIdL_TrkIdVL);
//   mymod->AddTrigger("HLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10_v1",              kHLT_DoubleEle10_CaloIdL_TrkIdVL_Ele10);
//   mymod->AddTrigger("HLT_DoubleMu5_Ele8_v2",                                 kHLT_DoubleMu5_Ele8);
//   mymod->AddTrigger("HLT_Mu5_DoubleEle8_v2",                                 kHLT_Mu5_DoubleEle8);
//   mymod->AddTrigger("HLT_TripleMu5_v2",                                      kHLT_TripleMu5);

  //Fake Rate triggers
  mymod->AddTrigger("HLT_Mu8_v1",                                kHLT_Mu8, kHLTObject_Mu8, "hltSingleMu8L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_v2",                                kHLT_Mu8, kHLTObject_Mu8, "hltSingleMu8L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_v3",                                kHLT_Mu8, kHLTObject_Mu8, "hltSingleMu8L3Filtered8");
  mymod->AddTrigger("HLT_Mu15_v1",                               kHLT_Mu15, kHLTObject_Mu15 ,  "hltL3Muon15");
  mymod->AddTrigger("HLT_Mu15_v2",                               kHLT_Mu15, kHLTObject_Mu15 , "hltL3Muon15");
  mymod->AddTrigger("HLT_Mu15_v3",                               kHLT_Mu15, kHLTObject_Mu15 , "hltSingleMu15L3Filtered15");
  mymod->AddTrigger("HLT_Mu15_v4",                               kHLT_Mu15, kHLTObject_Mu15 , "hltSingleMu15L3Filtered15");
  mymod->AddTrigger("HLT_Mu8_Jet40_v1",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8" );
  mymod->AddTrigger("HLT_Mu8_Jet40_v2",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Jet40_v3",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Jet40_v4",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Jet40_v5",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Jet40_v6",                          kHLT_Mu8_Jet40,kHLTObject_Mu8,"hltL3Mu8Jet20L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,kHLTObject_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Mu8,"hltSingleMu8EG5L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,kHLTObject_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8,"hltSingleMu8EG5L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,kHLTObject_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8,"hltSingleMu8EG5L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v4",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,kHLTObject_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8,"hltSingleMu8EG5L3Filtered8");
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v5",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,kHLTObject_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8,"hltSingleMu8EG5L3Filtered8");
  mymod->AddTrigger("HLT_Ele8_v1",                               kHLT_Ele8,kHLTObject_Ele8 , "hltEle8PixelMatchFilter" );
  mymod->AddTrigger("HLT_Ele8_v2",                               kHLT_Ele8,kHLTObject_Ele8 , "hltEle8PixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_v3",                               kHLT_Ele8,kHLTObject_Ele8 , "hltEle8PixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_v4",                               kHLT_Ele8,kHLTObject_Ele8 , "hltEle8PixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_v5",                               kHLT_Ele8,kHLTObject_Ele8 , "hltEle8PixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",             kHLT_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",             kHLT_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",             kHLT_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v4",             kHLT_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v5",             kHLT_Ele8_CaloIdL_CaloIsoVL, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, kHLTObject_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",            kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",            kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",            kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v4",            kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v5",            kHLT_Ele17_CaloIdL_CaloIsoVL, kHLTObject_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, kHLTObject_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter");
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , kHLTObject_Photon20_CaloIdVT_IsoT , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLTObject_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" );
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , kHLTObject_Photon20_CaloIdVT_IsoT, "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLTObject_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" );
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , kHLTObject_Photon20_CaloIdVT_IsoT, "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLTObject_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" );
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , kHLTObject_Photon20_CaloIdVT_IsoT, "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLTObject_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" );
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , kHLTObject_Photon20_CaloIdVT_IsoT, "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLTObject_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" );



//   mymod->AddTrigger("HLT_Jet30_v1", kHLT_Jet30 , kHLTObject_Jet, "hltSingleJet30");
//   mymod->AddTrigger("HLT_Jet30_v2", kHLT_Jet30 , kHLTObject_Jet, "hltSingleJet30");
//   mymod->AddTrigger("HLT_Jet60_v1", kHLT_Jet60 , kHLTObject_Jet, "hltSingleJet60Regional");
//   mymod->AddTrigger("HLT_Jet60_v2", kHLT_Jet60 , kHLTObject_Jet, "hltSingleJet60Regional");
//   mymod->AddTrigger("HLT_Jet110_v1", kHLT_Jet110 , kHLTObject_Jet, "hltSingleJet110Regional");
//   mymod->AddTrigger("HLT_Jet110_v2", kHLT_Jet110 , kHLTObject_Jet, "hltSingleJet110Regional");


  //Add L1SeedModules
  mymod->AddL1SeedModule("hltL1sL1SingleEG20",    kL1_SingleEG20 );
  mymod->AddL1SeedModule("hltL1sL1SingleEG12",    kL1_SingleEG12 );
  mymod->AddL1SeedModule("hltL1sL1SingleEG5",    kL1_SingleEG5 );
  mymod->AddL1SeedModule("hltL1sL1SingleMu3",    kL1_SingleMu3 );



  //Add L1Triggers
  mymod->AddL1Trigger("L1_SingleEG5",    kL1_SingleEG5 );
  mymod->AddL1Trigger("L1_SingleEG12", kL1_SingleEG12  );                  
  mymod->AddL1Trigger("L1_SingleEG20", kL1_SingleEG20  );                  
  mymod->AddL1Trigger("L1_SingleEG30", kL1_SingleEG30  );                  
  mymod->AddL1Trigger("L1_SingleMu3", kL1_SingleMu3 );                     
  mymod->AddL1Trigger("L1_SingleMu7 ", kL1_SingleMu7 );                    
  mymod->AddL1Trigger("L1_SingleMu10",kL1_SingleMu10 );                    
  mymod->AddL1Trigger("L1_SingleMu12", kL1_SingleMu12 );                   
  mymod->AddL1Trigger("L1_SingleMu20", kL1_SingleMu20 );
  mymod->AddL1Trigger("L1_SingleMu25", kL1_SingleMu25 );
  mymod->AddL1Trigger("L1_MuOpen_EG12", kL1_MuOpen_EG12 );                 
  mymod->AddL1Trigger("L1_MuOpen_EG5", kL1_MuOpen_EG5 );                   
  mymod->AddL1Trigger("L1_Mu3_EG5",   kL1_Mu3_EG5  );                    
  mymod->AddL1Trigger("L1_Mu5_EG12",   kL1_Mu5_EG12  );                    
  mymod->AddL1Trigger("L1_Mu7_EG5",      kL1_Mu7_EG5 );                    
  mymod->AddL1Trigger("L1_Mu12_EG5",  kL1_Mu12_EG5 );                      
  mymod->AddL1Trigger("L1_DoubleMu0", kL1_DoubleMu0 );                     
  mymod->AddL1Trigger("L1_DoubleMu5",  kL1_DoubleMu5 );                    
  mymod->AddL1Trigger("L1_DoubleEG3", kL1_DoubleEG3  );                    
  mymod->AddL1Trigger("L1_DoubleEG5", kL1_DoubleEG5  );                    
  mymod->AddL1Trigger("L1_DoubleEG_12_5", kL1_DoubleEG_12_5  );            
  mymod->AddL1Trigger("L1_Mu3_Jet20_Central", kL1_Mu3_Jet20_Central  );    
  mymod->AddL1Trigger("L1_EG5_Jet36_deltaPhi1",  kL1_EG5_Jet36_deltaPhi1 );


  mymod->SetPrintHLT(kTRUE); // print HLT table at start of analysis?
  



  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------

   hltmod->Add(goodPVFilterMod);

  //------------------------------------------------------------------------------------------------
  // Run Lepton + Denominator Skim
  //------------------------------------------------------------------------------------------------
  //   goodPVFilterMod->Add(muonID);

  //loose+loose
  if (runSkim == 1) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModDoubleLoose);
    selModDoubleLoose->Add(muonID); 
  } 
  
  //loose
  else if (runSkim == 2) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModLoose);
    selModLoose->Add(muonID); 
  }
  //tight+loose
  else if (runSkim == 3) {
    goodPVFilterMod->Add(muTightId);
    muTightId->Add(electronTightId);
    electronTightId->Add(muDenominator);
    muDenominator->Add(electronDenominator);
    electronDenominator->Add(mergedTight);
    mergedTight->Add(mergedLoose);
    mergedLoose->Add(selModTight);
    selModTight->Add(selModDoubleLoose); 
    selModDoubleLoose->Add(muonID);
  } else {
   goodPVFilterMod->Add(muonID);
  }

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  muonID->Add(electronID);
  electronID->Add(photonID);
  photonID->Add(tauID);
  tauID->Add(pubJet);
  pubJet->Add(pubPFMet); 
  pubPFMet->Add(jetCorr); 
  jetCorr->Add(jetID);
  jetID->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning->Add(tauCleaning);
  tauCleaning->Add(jetCleaning);
  jetCleaning->Add(jetIDNoPtCut);
  jetIDNoPtCut->Add(leptonMerger);
  leptonMerger->Add(signalFilter);
  signalFilter->Add(mymod);



  //------------------------------------------------------------------------------------------------
  //
  // setup analysis object
  //
  //------------------------------------------------------------------------------------------------
  Analysis *ana = new Analysis;
  ana->SetUseHLT(useHlt);
  if(nEvents >= 0) 
    ana->SetProcessNEvents(nEvents);
  if (useGen) {
    ana->AddSuperModule(generatorMod);
    generatorMod->Add(goodPVFilterMod);
    cout << "use gen\n";
  } else {
    ana->AddSuperModule(hltmod);
//     ana->AddSuperModule(goodPVFilterMod);    
    cout << "not gen\n";
  }


  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
//   printf("\nRely on Catalog: %s\n",catalogDir);
//   printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
//   Catalog *c = new Catalog(catalogDir);
//   TString skimdataset = TString(dataset)+TString("/") +TString(skim);
//   Dataset *d = NULL;
//   if (TString(skim).CompareTo("noskim") == 0)
//     d = c->FindDataset(book,dataset,fileset);
//   else 
//     d = c->FindDataset(book,skimdataset.Data(),fileset);
//   ana->AddDataset(d);

//     ana->AddFile("/data/smurf/sixie/BAMBU/020/r11-express-*.root");
//     ana->AddFile("/server/03a/mitprod/skimExpress/skimExpress_doublelepton_*.root");
//    ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/p11-h160ww2l-gf-v1g1-pu/3EABE6D3-F150-E011-A1E9-00A0D1EE89E0.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/r11a-mueg-pr-v1/D2EFD897-095A-E011-94A8-0030487CBD0A.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/r11a-smu-pr-v1/1A4E4D44-F556-E011-B246-000423D98B6C.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/r11a-sel-pr-v2/AE30DEFF-F772-E011-8B67-003048F0258C.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/r11a-sel-pr-v2/2277B418-7975-E011-98D4-0019DB29C5FC.root");
  ana->AddFile("/castor/cern.ch/user/p/paus/filefi/020/p11-vvj-v1g1-pu_1//3E906C34-6D5A-E011-9A08-0024E87663EE.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName(rootFile.Data());
  ana->SetCacheSize(0);


  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

