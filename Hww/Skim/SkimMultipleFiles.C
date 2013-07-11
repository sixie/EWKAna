//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-dmu-pr-v4.txt\",\"WWAnalysisSkimmed_r11a-dmu-pr-v4_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-del-pr-v4.txt\",\"WWAnalysisSkimmed_r11a-del-pr-v4_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-mueg-pr-v4.txt\",\"WWAnalysisSkimmed_r11a-mueg-pr-v4_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-smu-pr-v4.txt\",\"WWAnalysisSkimmed_r11a-smu-pr-v4_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-sel-pr-v4.txt\",\"WWAnalysisSkimmed_r11a-sel-pr-v4_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)

//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-dmu-m10-v1.txt\",\"WWAnalysisSkimmed_r11a-dmu-m10-v1_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-del-m10-v1.txt\",\"WWAnalysisSkimmed_r11a-del-m10-v1_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-mueg-m10-v1.txt\",\"WWAnalysisSkimmed_r11a-mueg-m10-v1_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-smu-m10-v1.txt\",\"WWAnalysisSkimmed_r11a-smu-m10-v1_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)
//root -l -b -q EWKAna/Hww/Skim/SkimMultipleFiles.C+\(\"EWKAna/Hww/Skim/SkimInput_r11a-sel-m10-v1.txt\",\"WWAnalysisSkimmed_r11a-sel-m10-v1_TightPlusDenominatorNoTriggerSkim.root\",kTRUE,kFALSE\)

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

//=== FUNCTION DECLARATIONS ======================================================================================
Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData, Bool_t isMC);
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);

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
void SkimMultipleFiles(string inputFilename, string outputFilename, Bool_t isMuonData, Bool_t isMC) 
{
  gBenchmark->Start("SkimNtuples");
    
  TTree::SetMaxTreeSize(kMaxLong64);
  
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


  // list input ntuple files to be skimmed
  vector<string> infilenames;  
  ifstream ifs;
  ifs.open(inputFilename.c_str()); 
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) { infilenames.push_back(line); }
  ifs.close();

  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Skimming " << infilenames[ifile] << "..." << endl;
    TTree *eventTree = getTreeFromFile(infilenames[ifile].c_str(),"Events"); 
    nEventsTotal += getNEvents(infilenames[ifile].c_str()); 
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
//       if (!passHLT(info->triggerBits, info->runNum, isMuonData, isMC)) keep = kFALSE;
//       if (!(info->triggerBits & kHLT_Ele10_SW_L1R)) keep = kFALSE;
    
      Int_t NRecoLeptons = 0;
      Int_t NMuons = 0;
      Int_t NDenominatorMuons = 0;
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
        if ( passMuonCuts(mu)
             &&
             fabs(mu->eta) < 2.4
             && 
             mu->pt > 10.0
          ) {
          NMuons++;
        }
        if ( fabs(mu->eta) < 2.4 && mu->pt > 10.0) {
          NRecoLeptons++;
        }
        if ( passMuonDenominatorCuts(mu)
             &&
             fabs(mu->eta) < 2.4
             && 
             mu->pt > 10.0
          ) {
          NDenominatorMuons++;
        }
      }
      Int_t NElectrons = 0;
      Int_t NDenominatorElectrons = 0;
     for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
        if ( passElectronCuts(ele)
             &&
             fabs(ele->eta) < 2.5
             && 
             ele->pt > 10.0
          ) {
          NElectrons++;
        }
        if ( fabs(ele->eta) < 2.5 && ele->pt > 10.0 ) {
          NRecoLeptons++;
        }
        if ( passElectronDenominatorCuts(ele)
             &&
             fabs(ele->eta) < 2.5
             && 
             ele->pt > 10.0
          ) {
          NDenominatorElectrons++;
        }
      }

//      //One Tight OR Two Loose
//      if (!(NDenominatorMuons+NDenominatorElectrons >= 2 || NElectrons+NMuons >= 1 )) {
//        keep = kFALSE;
//      }
     //One Tight OR Two Loose
     if (!(NRecoLeptons >= 2 || NElectrons+NMuons >= 1 )) {
       keep = kFALSE;
     }




      if(keep) {
        outEventTree->Fill();
        nPassEvts++;
      }
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

Bool_t passHLT(UInt_t triggerBits, Int_t runNum, Bool_t isMuonData, Bool_t isMC) {


  Bool_t pass = kFALSE;
//   if (isMC) {
//     if (triggerBits & kHLT_Mu9) pass = kTRUE;
//     if (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) pass = kTRUE;
//     if (triggerBits & kHLT_Ele15_LW_L1R) pass = kTRUE;
//   } else {
//     if (isMuonData) {
//       if ((runNum >= 136033) && (runNum <= 147116)) {
//         if ( (triggerBits & kHLT_Mu9) ) pass = kTRUE;
//       } 
//       if ((runNum >= 136033) && (runNum <= 139980)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 140058) && (runNum <= 141882)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 141956) && (runNum <= 144114)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 146428) && (runNum <= 147116)) {
//         if ( (triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 149442)) {
//         if ( (triggerBits & kHLT_Mu15) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 148058)) {
//         if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 148819) && (runNum <= 149442)) {
//         if ( (triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
//       }
//     } else {
//       pass = kTRUE;
//       //it's electron data
//       if ((runNum >= 136033) && (runNum <= 139980)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele10_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 140058) && (runNum <= 141882)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 141956) && (runNum <= 144114)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele15_SW_CaloEleId_L1R) ) pass = kTRUE;
//       } 
//       if ((runNum >= 146428) && (runNum <= 147116)) {
//         if ( !(triggerBits & kHLT_Mu9) && (triggerBits & kHLT_Ele17_SW_CaloEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 147196) && (runNum <= 148058)) {
//         if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TightEleId_L1R) ) pass = kTRUE;
//       }
//       if ((runNum >= 148819) && (runNum <= 149442)) {
//         if ( !(triggerBits & kHLT_Mu15) && (triggerBits & kHLT_Ele17_SW_TighterEleIdIsol_L1R) ) pass = kTRUE;
//       }
//     }
//   }
  return pass;

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

Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (!(fabs(ele->eta) <= 2.5)) pass = kFALSE;

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




Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) <= 2.5)) pass = kFALSE;

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
            && (ele->trkIso03) / ele->pt < 0.20
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

Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) > 2.4) pass = kFALSE;

  Double_t iso04 = mu->ChargedIso04 + mu->NeutralIso04_10Threshold;
  Double_t iso03 = mu->ChargedIso03 + mu->NeutralIso03_10Threshold;
  Double_t iso = iso03;
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
        ( fabs(mu->d0) < 0.01
          )
      ) {
      pass = kFALSE;
    }    
  }

  return pass;
}


Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

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

  return pass;
}
