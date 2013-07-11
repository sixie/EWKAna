//foreach f ( /home/sixie/hist/AllNtuple/HWWNtuple/data/unprocessed/AllNtuple_HWWNtuple_r11a-del-pr-v4_noskim_????.root )
//root -l -b -q EWKAna/Hww/Skim/SkimTightPlusRecoPerFile.C+\(\"$f\",\"$f.TightPlusRecoSkimmed.root\"\)
//end



#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TH1F.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#endif

// define structures to read in ntuple
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#include "EWKAna/Utils/LeptonIDCuts.hh"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"

//=== FUNCTION DECLARATIONS ======================================================================================
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);

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

Double_t getNEvents(string filename) {

  //Get Number of Events in the Sample
  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny("AnaFwkMod");
  if (!dir) {
    cout << "Could not find directory AnaFwkMod"
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *hist = (TH1F*)dir->Get("hDAllEvents");
  if (!hist) {
    cout << "Could not find histogram hDEvents in directory AnaFwkMod"
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }  
  return hist->Integral();
}

// Main macro function
//--------------------------------------------------------------------------------------------------
void SkimTightPlusRecoPerFile(string inputFilename, string outputFilename) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  

  mithep::MuonIDMVA *muonIDMVA = new mithep::MuonIDMVA();
  muonIDMVA->Initialize("BDTG method",
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin0_IDIsoCombined_BDTG.weights.xml", 
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin0_IDIsoCombined_BDTG.weights.xml", 
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin1_IDIsoCombined_BDTG.weights.xml", 
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin1_IDIsoCombined_BDTG.weights.xml", 
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/BarrelPtBin2_IDIsoCombined_BDTG.weights.xml", 
                        "/home/sixie/CMSSW_analysis/src/MitPhysics/data/MuonMVAWeights/EndcapPtBin2_IDIsoCombined_BDTG.weights.xml",
                        mithep::MuonIDMVA::kIDIsoCombinedDetIso);
  
  mithep::ElectronIDMVA *electronIDMVA = new mithep::ElectronIDMVA();
  electronIDMVA->Initialize("BDTG method",
                            "MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_IDIsoCombined_BDTG.weights.xml", 
                            "MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_IDIsoCombined_BDTG.weights.xml", 
                            "MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_IDIsoCombined_BDTG.weights.xml", 
                            "MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_IDIsoCombined_BDTG.weights.xml", 
                            "MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_IDIsoCombined_BDTG.weights.xml", 
                            "MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_IDIsoCombined_BDTG.weights.xml",
                            mithep::ElectronIDMVA::kIDIsoCombined);



  // Don't write TObject part of the objects
  mithep::TEventInfo::Class()->IgnoreTObjectStreamer();
  mithep::TMuon::Class()->IgnoreTObjectStreamer();
  mithep::TElectron::Class()->IgnoreTObjectStreamer();
  mithep::TJet::Class()->IgnoreTObjectStreamer();
  mithep::TPhoton::Class()->IgnoreTObjectStreamer();

  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr     = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr    = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr     = new TClonesArray("mithep::TPhoton");
    
  UInt_t nInputEvts = 0;
  UInt_t nPassEvts  = 0;
  UInt_t nEventsTotal = 0;

  TFile* outfile = new TFile(outputFilename.c_str(), "RECREATE");
  
  //
  // Initialize data trees and structs
  // 
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Info",     &info);
  outEventTree->Branch("Muon",     &muonArr);
  outEventTree->Branch("Electron", &electronArr);
  outEventTree->Branch("PFJet",    &jetArr);
  outEventTree->Branch("Photon",   &photonArr);

  cout << "Skimming " << inputFilename << "..." << endl;
  TTree *eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  nEventsTotal += getNEvents(inputFilename.c_str()); 
  assert(eventTree);
    
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");
  eventTree->SetBranchAddress("Photon",    &photonArr);  TBranch *photonBr   = eventTree->GetBranch("Photon");
     
  for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) { 
    infoBr->GetEntry(ientry);
    muonArr->Clear();     muonBr->GetEntry(ientry);
    electronArr->Clear(); electronBr->GetEntry(ientry);      
    jetArr->Clear(); jetBr->GetEntry(ientry);
    photonArr->Clear(); photonBr->GetEntry(ientry);
     
    if (ientry % 100000 == 0) cout << "Events: " << ientry << endl;

    nInputEvts++;
      
    Bool_t keep = kTRUE;
    
    Double_t rho = 0;
    if (!(TMath::IsNaN(info->PileupEnergyDensity) || isinf(info->PileupEnergyDensity))) rho = info->PileupEnergyDensity;

    Int_t NTightMuons = 0;
    Int_t NRecoMuons = 0;
    for(Int_t i=0; i<muonArr->GetEntries(); i++) {
      const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
      if ( fabs(mu->eta) < 2.4
           && 
           mu->pt > 10.0
        ) {
        NRecoMuons++;

//         if ( passMuonCuts(mu)) NTightMuons++;
        if ( passMuonMVAIDIsoCombined(mu, muonIDMVA->MVAValue(
                                          mu->pt, mu->eta,
                                          mu->tkNchi2,
                                          mu->muNchi2,
                                          mu->nValidHits,
                                          mu->nTkHits,
                                          mu->nPixHits,
                                          mu->nMatch,
                                          mu->d0,
                                          mu->ip3d,
                                          mu->ip3dSig,
                                          mu->TrkKink,
                                          mu->SegmentCompatibility,
                                          mu->CaloCompatilibity,
                                          (mu->HadEnergy - rho*MuonEffectiveArea(kMuHadEnergy,mu->eta))/mu->pt,
                                          (mu->HoEnergy - rho*MuonEffectiveArea(kMuHoEnergy,mu->eta))/mu->pt,
                                          (mu->EmEnergy - rho*MuonEffectiveArea(kMuEmEnergy,mu->eta))/mu->pt,
                                          (mu->HadS9Energy - rho*MuonEffectiveArea(kMuHadS9Energy,mu->eta))/mu->pt,
                                          (mu->HoS9Energy - rho*MuonEffectiveArea(kMuHoS9Energy,mu->eta))/mu->pt,
                                          (mu->EmS9Energy - rho*MuonEffectiveArea(kMuEmS9Energy,mu->eta))/mu->pt,
                                          (mu->trkIso03 - rho*MuonEffectiveArea(kMuTrkIso03,mu->eta))/mu->pt,
                                          (mu->emIso03 - rho*MuonEffectiveArea(kMuEMIso03,mu->eta))/mu->pt,
                                          (mu->hadIso03 - rho*MuonEffectiveArea(kMuHadIso03,mu->eta))/mu->pt,
                                          (mu->trkIso05 - rho*MuonEffectiveArea(kMuTrkIso05,mu->eta))/mu->pt,
                                          (mu->emIso05 - rho*MuonEffectiveArea(kMuEMIso05,mu->eta))/mu->pt,
                                          (mu->hadIso05 - rho*MuonEffectiveArea(kMuHadIso05,mu->eta))/mu->pt
                                        ), rho)) NTightMuons++;
      }
    }
    Int_t NTightElectrons = 0;
    Int_t NRecoElectrons = 0;
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      if ( fabs(ele->eta) < 2.5
           && 
           ele->pt > 10.0
        ) {
        NRecoElectrons++;

//         if (passElectronCuts(ele)) NTightElectrons++;
        if (passElectronMVAIDIsoCombined(ele, 
                                         electronIDMVA->MVAValue(
                                           ele->pt,ele->scEta,rho,
                                           ele->sigiEtaiEta, 
                                           ele->deltaEtaIn,
                                           ele->deltaPhiIn, 
                                           ele->HoverE,
                                           ele->d0,
                                           ele->dz, 
                                           ele->fBrem,
                                           ele->EOverP,
                                           ele->ESeedClusterOverPout,
                                           TMath::Sqrt(ele->sigiPhiiPhi),
                                           ele->nBrem,
                                           (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p), 
                                           ele->ESeedClusterOverPIn,
                                           ele->ip3d,
                                           ele->ip3dSig,
                                           ele->GsfTrackChi2OverNdof,
                                           ele->dEtaCalo,
                                           ele->dPhiCalo,
                                           ele->R9,
                                           ele->SCEtaWidth,
                                           ele->SCPhiWidth,
                                           ele->CovIEtaIPhi,
                                           ele->PreShowerOverRaw,
                                           ele->ChargedIso03,
                                           (ele->NeutralHadronIso03_05Threshold - ele->NeutralHadronIso007_05Threshold),
                                           (ele->GammaIso03_05Threshold - ele->GammaIsoVetoEtaStrip03_05Threshold),
                                           ele->ChargedIso04 ,
                                           (ele->NeutralHadronIso04_05Threshold - ele->NeutralHadronIso007_05Threshold),
                                           (ele->GammaIso04_05Threshold - ele->GammaIsoVetoEtaStrip04_05Threshold) 
                                           ), 
                                         0)) NTightElectrons++;



      }
    }

    //One Tight AND Two Loose
    if (!( NTightElectrons + NTightMuons >= 1  
           && NRecoElectrons + NRecoMuons >= 2
          ) ) {
      keep = kFALSE;
    }
      
    if(keep) {
      outEventTree->Fill();
      nPassEvts++;
    }
  }
  outfile->Write();
  outfile->Close();
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
    
  std::cout << outputFilename << " created!" << std::endl;
  std::cout << " >>> Total Number of Events: " << nEventsTotal << std::endl;
  std::cout << " >>> Events processed: " << nInputEvts << std::endl;
  std::cout << " >>>   Events passing: " << nPassEvts << std::endl;

  gBenchmark->Show("SkimNtuples");
}  

Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (fabs(ele->eta) > 2.5) pass = kFALSE;

  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->NeutralHadronIso007_10Threshold;
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
    +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
    -ele->NeutralHadronIso007_10Threshold;
  Double_t iso = iso04;
  Double_t isoCutValue = 0;

  if (fabs(ele->scEta) < 1.479) {
    if (ele->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.13;
    }
  } else {
    if (ele->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.09;
    }
  }


  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.10
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }
  } 


  if (ele->pt < 20) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.03
              && ele->HoverE < 0.025
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else  {
      if (! (  (0==0)
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
            )
        ) {
        pass = kFALSE;
      }
    } 

    if (ele->fBrem <= 0.15) {
      if (fabs(ele->scEta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EOverP > 0.95 )) pass = kFALSE;
      }
    }
  }



  return pass;
}

Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  Double_t iso = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
  Double_t isoCutValue = 0;
  if (fabs(mu->eta) < 1.479) {
    if (mu->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (mu->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  }

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
                
        && mu->nTkHits > 10
        && (mu->nPixHits > 0)
        && fabs(mu->d0) < 0.02
        && fabs(mu->dz) < 0.1
        && iso / mu->pt < isoCutValue
        && (mu->pterr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  if (mu->pt < 20) {
    if (!
        ( 
          fabs(mu->d0) < 0.01
          )
      ) pass = kFALSE;    
  }
  
  return pass;
}
