//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwFakeRateNtupler.C+\(\"0000\",\"noskim\",\"r11a-dmu-pr-v1\",\"cern/filefi/020\",\"/home/mitprod/catalog\",\"HwwNtuple\",10000,-1\) 
//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwFakeRateNtupler.C+\(\"0000\",\"noskim\",\"f10-h160ww2l-gf-z2-v12\",\"cern/filler/015\",\"/home/ceballos/catalog\",\"HwwNtuple\",-1,2\) 
//root -l -q -b $CMSSW_BASE/src/EWKAna/Ntupler/macros/runHwwFakeRateNtupler.C+\(\"0000\",\"noskim\",\"f10-zmm-powheg-c10-v12\",\"cern/filefi/015\",\"/home/mitprod/catalog\",\"HwwNtuple\",10000,2\) 

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
void runHwwFakeRateNtupler(
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
  Bool_t fakeRateSkim = kTRUE;
  Bool_t computePDFWeights = kFALSE; 
  string pdfSetName = "";

  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // supercluster kinematics
  const Double_t SCEtMin  = 10;
  const Double_t SCEtMax  = 7000;
  const Double_t SCEtaMin = -3;
  const Double_t SCEtaMax =  3;
  
  // muon kinematics
  const Double_t MuonPtMin  = 3;
  const Double_t MuonPtMax  = 7000;
  const Double_t MuonEtaMin = -3;
  const Double_t MuonEtaMax =  3;

  // jet requirements
  const Double_t jetEtMin = 15;

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
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetVertexesName("DAPrimaryVertexes");
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);

 //------------------------------------------------------------------------------------------------
  // organize selection
  //------------------------------------------------------------------------------------------------

  MuonIDMod *muId = new MuonIDMod;  
  muId->SetPtMin     (17.0);
  muId->SetApplyD0Cut(kTRUE);
  muId->SetClassType ("Global");
  muId->SetIDType    ("WWMuId");
  muId->SetIsoType   ("TrackCaloSliding");
  muId->SetOutputName("SkimMuons");

  MuonIDMod *muIdDenominator = new MuonIDMod;  
  muIdDenominator ->SetPtMin     (10.0);
  muIdDenominator ->SetApplyD0Cut(kFALSE);
  muIdDenominator ->SetClassType ("Global");
  muIdDenominator ->SetIDType    ("NoId");
  muIdDenominator ->SetIsoType   ("NoIso");
  muIdDenominator ->SetOutputName("DenominatorMuons");

  ElectronIDMod *electronsIsolated = new ElectronIDMod;
  electronsIsolated->SetInputName                 (Names::gkElectronBrn);
  electronsIsolated->SetPtMin                     (10.0);
  electronsIsolated->SetApplyConversionFilterType1(kFALSE);
  electronsIsolated->SetApplyConversionFilterType2(kFALSE);
  electronsIsolated->SetChargeFilter              (kFALSE);
  electronsIsolated->SetApplySpikeRemoval         (kTRUE);
  electronsIsolated->SetApplyD0Cut                (kFALSE);
  electronsIsolated->SetNExpectedHitsInnerCut     (999);
  electronsIsolated->SetIDType                    ("NoId");
  electronsIsolated->SetIsoType                   ("TrackJuraSliding");
  electronsIsolated->SetOutputName("ElectronsIsolated");

  ElectronIDMod       *electronVBTF80NonIsolated    = new ElectronIDMod;
  electronVBTF80NonIsolated->SetPtMin(10.0);
  electronVBTF80NonIsolated->SetEtaMax(999.);
  electronVBTF80NonIsolated->SetIDType(TString("VBTFWorkingPoint80Id"));
  electronVBTF80NonIsolated->SetIsoType(TString("NoIso"));
  electronVBTF80NonIsolated->SetNExpectedHitsInnerCut(0);
  electronVBTF80NonIsolated->SetApplyConversionFilterType1(kFALSE);
  electronVBTF80NonIsolated->SetApplyConversionFilterType2(kTRUE);
  electronVBTF80NonIsolated->SetApplyD0Cut                (kFALSE);
  electronVBTF80NonIsolated->SetChargeFilter              (kFALSE);
  electronVBTF80NonIsolated->SetApplySpikeRemoval         (kTRUE);
  electronVBTF80NonIsolated->SetApplyTriggerMatching      (kFALSE);
  electronVBTF80NonIsolated->SetOutputName("ElectronsVBTF80NonIsolated");

  ElectronIDMod       *electronVBTF90PartiallyIsolated    = new ElectronIDMod;
  electronVBTF90PartiallyIsolated->SetPtMin(10.0);
  electronVBTF90PartiallyIsolated->SetEtaMax(999.);
  electronVBTF90PartiallyIsolated->SetIDType(TString("VBTFWorkingPoint90Id"));
  electronVBTF90PartiallyIsolated->SetIsoType(TString("Custom"));
  electronVBTF90PartiallyIsolated->SetNExpectedHitsInnerCut(0);
  electronVBTF90PartiallyIsolated->SetApplyConversionFilterType1(kFALSE);
  electronVBTF90PartiallyIsolated->SetApplyConversionFilterType2(kTRUE);
  electronVBTF90PartiallyIsolated->SetApplyD0Cut                (kFALSE);
  electronVBTF90PartiallyIsolated->SetChargeFilter              (kFALSE);
  electronVBTF90PartiallyIsolated->SetApplySpikeRemoval         (kTRUE);
  electronVBTF90PartiallyIsolated->SetApplyTriggerMatching      (kFALSE);
  electronVBTF90PartiallyIsolated->SetOutputName("ElectronsVBTF90PartiallyIsolated");


  //------------------------------------------------------------------------------------------------
  // publisher Mod
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubJet->SetInputName("AKt5PFJets");
  pubJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("pubPFMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->AddCorrectionFromFile("MitPhysics/data/START38_V13_AK5PF_L2Relative.txt"); 
  jetCorr->AddCorrectionFromFile("MitPhysics/data/START38_V13_AK5PF_L3Absolute.txt");
  if(!useGen){ 
    jetCorr->AddCorrectionFromFile("MitPhysics/data/START38_V13_AK5PF_L2L3Residual.txt");
  }
  jetCorr->SetInputName(pubJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonID        = new MuonIDMod;  
  muonID->SetClassType("Global");
  muonID->SetIDType("WWMuId");
  muonID->SetIsoType("TrackCaloSliding");
  muonID->SetApplyD0Cut(kTRUE);

  ElectronIDMod       *electronID    = new ElectronIDMod;
  electronID->SetIDType("VBTFWorkingPoint80Id");
  electronID->SetIsoType("TrackJuraSliding");
  electronID->SetApplyConversionFilterType1(kFALSE);
  electronID->SetApplyConversionFilterType2(kTRUE);
  electronID->SetChargeFilter(kFALSE);
  electronID->SetApplyD0Cut(kTRUE);
  electronID->SetNExpectedHitsInnerCut(0);

  PhotonIDMod         *photonID      = new PhotonIDMod;
  TauIDMod            *tauID         = new TauIDMod;
  JetIDMod            *jetID         = new JetIDMod;
  jetID->SetInputName(jetCorr->GetOutputName());
  jetID->SetPtCut(25.0);
  jetID->SetEtaMaxCut(5.0);
  jetID->SetJetEEMFractionMinCut(0.0);
  jetID->SetOutputName("GoodJets");

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

  JetCleaningMod      *jetCleaningNoPtCut = new JetCleaningMod;
  jetCleaningNoPtCut->SetGoodJetsName("GoodJetsNoPtCut");
  jetCleaningNoPtCut->SetCleanJetsName("CleanJetsNoPtCut");

  // get the leptons
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
  TString rootFile = TString(outputName);
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
  mymod->SetFakeRateSkim(fakeRateSkim);            // keep all electrons even if not used in a pre-selected dielectron
  mymod->SetFSRMode(fsrmode);
  mymod->SetSCEtMin(SCEtMin);
  mymod->SetSCEtMax(SCEtMax);
  mymod->SetSCEtaMin(SCEtaMin);
  mymod->SetSCEtaMax(SCEtaMax);
  mymod->SetMuonPtMin(MuonPtMin);
  mymod->SetMuonPtMax(MuonPtMax);
  mymod->SetMuonEtaMin(MuonEtaMin);
  mymod->SetMuonEtaMax(MuonEtaMax);
  mymod->SetCleanJetsName(jetCleaning->GetOutputName());
  mymod->SetCleanJetsNoPtCutName(jetCleaningNoPtCut->GetOutputName());
  mymod->SetJetEtMin(jetEtMin);
  mymod->SetComputePDFWeights(computePDFWeights);
  mymod->SetPDFName(pdfSetName.c_str());

  mymod->AddTrigger("HLT_Mu8_v1",                                kHLT_Mu8);
  mymod->AddTrigger("HLT_Mu15_v1",                               kHLT_Mu15);
  mymod->AddTrigger("HLT_Mu15_v2",                               kHLT_Mu15);
  mymod->AddTrigger("HLT_Mu8_Jet40_v1",                          kHLT_Mu8_Jet40);
  mymod->AddTrigger("HLT_Mu8_Jet40_v2",                          kHLT_Mu8_Jet40);
  mymod->AddTrigger("HLT_Mu8_Jet40_v3",                          kHLT_Mu8_Jet40);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",         kHLT_Mu8_Photon20_CaloIdVT_IsoT);
  mymod->AddTrigger("HLT_Ele8_v1",                               kHLT_Ele8);
  mymod->AddTrigger("HLT_Ele8_v2",                               kHLT_Ele8);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",             kHLT_Ele8_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",             kHLT_Ele8_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",            kHLT_Ele17_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",            kHLT_Ele17_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
  mymod->AddTrigger("HLT_Jet30_v1", kHLT_Jet30);
  mymod->AddTrigger("HLT_Jet60_v1", kHLT_Jet60);
  mymod->AddTrigger("HLT_Jet110_v1", kHLT_Jet110);
  
  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  

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
    
  } else {
    ana->AddSuperModule(goodPVFilterMod);    
  }


  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n\n",book,dataset,skim,fileset);
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);

//    ana->AddFile("/castor/cern.ch/user/p/paus/filefi/015/r10b-mu-pr-v2/C86AB5EF-6FE0-DF11-ABDF-000423D98B6C.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/015/r10b-mu-pr-v2/6225A696-84DA-DF11-9EF7-0030487CAF0E.root");
//   ana->AddFile("/castor/cern.ch/user/p/paus/filefi/018/r10a-eg-d22/A01B80F5-500F-E011-AA56-003048C693D6.root");
//    ana->AddFile("/server/03a/mitprod/skimExpress/skimExpress_singlelepton_*.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName(rootFile.Data());

  

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  goodPVFilterMod->Add(electronsIsolated);
  electronsIsolated->Add(electronVBTF80NonIsolated);
  electronVBTF80NonIsolated->Add(electronVBTF90PartiallyIsolated);
  electronVBTF90PartiallyIsolated->Add(muIdDenominator);
  muIdDenominator->Add(muonID);  
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
  jetIDNoPtCut->Add(jetCleaningNoPtCut);
  jetCleaningNoPtCut->Add(leptonMerger);
  leptonMerger->Add(signalFilter);
  signalFilter->Add(mymod);

  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

