#include "EWKAna/Ntupler/interface/HwwGenNtuplerMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/CaloJet.h"
#include "MitAna/DataTree/interface/TrackJet.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/L1TriggerMask.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitAna/DataTree/interface/CaloMet.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/DecayParticle.h"
#include "MitAna/DataTree/interface/StableData.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
// #include "LHAPDF/LHAPDF.h"
#include <vector>

//#define __DATA_36X__

using namespace mithep;

ClassImp(mithep::HwwGenNtuplerMod)

HwwGenNtuplerMod::HwwGenNtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
  fOutputName    ("ntuple.root"),
  fPDFName       ("cteq66.LHgrid"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fParticles     (0),
  fMCEvtInfo     (0),
  fGenJets       (0),
  fUseGen        (kTRUE),
  fComputePDFWeights(kFALSE),
  fEventTree     (0)
{
  // Constructor

  // Don't write TObject part of the objects
  TEventInfo::Class()->IgnoreTObjectStreamer();
  TGenInfo::Class()->IgnoreTObjectStreamer();
  TElectron::Class()->IgnoreTObjectStreamer();
  TMuon::Class()->IgnoreTObjectStreamer();
  TJet::Class()->IgnoreTObjectStreamer();
  TPhoton::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
HwwGenNtuplerMod::~HwwGenNtuplerMod()
{
  // Destructor
}
	//--------------------------------------------------------------------------------------------------      
void HwwGenNtuplerMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::BeginRun()
{

}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,      fParticles); 
  ReqBranch(fMCEvtInfoName, fMCEvtInfo);
  ReqBranch(Names::gkGenJetBrn , fGenJets); 
  ReqBranch("PFMet",     fPFMet);


  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");

  //
  // Initialize data trees and structs
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  fEventTree->Branch("Gen",&fGenInfo);
  

//   if (fComputePDFWeights) {
//  //    gSystem->Setenv("LHAPATH","/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/lhapdf/5.6.0-cms4/share/lhapdf/PDFsets/");
//     LHAPDF::setVerbosity(LHAPDF::SILENT);
//     LHAPDF::initPDFSet(fPDFName.Data());
//     LHAPDF::getDescription();
    
//     //Branches for PDFset weights
//     fEventTree->Branch("NPDFMembers",&fNPDFMembers,"NPDFMembers/I");  
//     fEventTree->Branch("PDFWeights",fPDFWeights,"PDFWeights[NPDFMembers]/F");
//   }


  //Tau Decay Branches
  
  fEventTree->Branch("NGenTaus",&fNGenTaus,"NGenTaus/I");  

  fEventTree->Branch("GenTauMass",fGenTauMass,"GenTauMass[NGenTaus]/F");
  fEventTree->Branch("GenTauPt",fGenTauPt,"GenTauPt[NGenTaus]/F");
  fEventTree->Branch("GenTauEta",fGenTauEta,"GenTauEta[NGenTaus]/F");
  fEventTree->Branch("GenTauPhi",fGenTauPhi,"GenTauPhi[NGenTaus]/F");

  fEventTree->Branch("GenTauDecayLeptonType",fGenTauDecayLeptonType,"GenTauDecayLeptonType[NGenTaus]/I");
  fEventTree->Branch("GenTauDecayLeptonPt",fGenTauDecayLeptonPt,"GenTauDecayLeptonPt[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayLeptonEta",fGenTauDecayLeptonEta,"GenTauDecayLeptonEta[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayLeptonPhi",fGenTauDecayLeptonPhi,"GenTauDecayLeptonPhi[NGenTaus]/F");

  fEventTree->Branch("GenTauDecayTauNeutrinoType",fGenTauDecayTauNeutrinoType,"GenTauDecayTauNeutrinoType[NGenTaus]/I");
  fEventTree->Branch("GenTauDecayTauNeutrinoPt",fGenTauDecayTauNeutrinoPt,"GenTauDecayTauNeutrinoPt[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayTauNeutrinoEta",fGenTauDecayTauNeutrinoEta,"GenTauDecayTauNeutrinoEta[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayTauNeutrinoPhi",fGenTauDecayTauNeutrinoPhi,"GenTauDecayTauNeutrinoPhi[NGenTaus]/F");

  fEventTree->Branch("GenTauDecayLeptonNeutrinoPt",fGenTauDecayLeptonNeutrinoPt,"GenTauDecayLeptonNeutrinoPt[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayLeptonNeutrinoEta",fGenTauDecayLeptonNeutrinoEta,"GenTauDecayLeptonNeutrinoEta[NGenTaus]/F");
  fEventTree->Branch("GenTauDecayLeptonNeutrinoPhi",fGenTauDecayLeptonNeutrinoPhi,"GenTauDecayLeptonNeutrinoPhi[NGenTaus]/F");

  fEventTree->Branch("NLeptons", &fNLeptons,"NLeptons/I");  
  fEventTree->Branch("LeptonType",fLeptonType,"LeptonType[NLeptons]/I");
  fEventTree->Branch("LeptonPt"  ,fLeptonPt  ,"LeptonPt[NLeptons]/F");
  fEventTree->Branch("LeptonEta" ,fLeptonEta ,"LeptonEta[NLeptons]/F");
  fEventTree->Branch("LeptonPhi" ,fLeptonPhi ,"LeptonPhi[NLeptons]/F");

  fEventTree->Branch("PFMet"   ,&fNtuplePFMet ,"PFMet/F");
  fEventTree->Branch("PFMetPhi",&fNtuplePFMetPhi,"PFMetPhi/F");


//   AddOutput(fEventTree);



}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  


}  

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::Process()
{
  //
  // Load branches
  //
  LoadBranch(fPartName);
  LoadBranch(fMCEvtInfoName);
  LoadBranch(Names::gkGenJetBrn);


 // GeneratorTools::PrintHepMCTable(fParticles, true, -1);     
  

 //***********************************************************************************************
  //Some MC Specific Requirements
  //***********************************************************************************************
  const MCParticleCol *GenLeptons = 0;
  GenLeptons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCLeptonsName,0);
  const MCParticleCol *GenNeutrinos = GetObjThisEvt<MCParticleCol>(ModNames::gkMCNeutrinosName,0);
  const MCParticleCol *GenTaus = GetObjThisEvt<MCParticleCol>(ModNames::gkMCTausName,0);
  const MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCBosonsName);    

  ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >(FindObjThisEvt(ModNames::gkCleanMuonsName));
  ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >(FindObjThisEvt(ModNames::gkCleanElectronsName));

  //Fill Gen Info
  FillGenInfo(GenLeptons, GenNeutrinos, GenBosons);  // fill the data structure


//   //************************************************************************************************
//   //Fill PDF information
//   //************************************************************************************************
//   if (fComputePDFWeights ) {
//     //cout << "Event " << fMCEvtInfo->Weight() << " " << fMCEvtInfo->Id1() << " " << fMCEvtInfo->Id2() << " " << fMCEvtInfo->X1() << " " << fMCEvtInfo->X2() << endl;    
//     Double_t Q    = fMCEvtInfo->Scale();
//     Int_t    id1  = fMCEvtInfo->Id1();
//     Double_t x1   = fMCEvtInfo->X1();
//     Double_t pdf1 = fMCEvtInfo->Pdf1();
//     Int_t    id2  = fMCEvtInfo->Id2();
//     Double_t x2   = fMCEvtInfo->X2();
//     Double_t pdf2 = fMCEvtInfo->Pdf2();
    
//     UInt_t nmembers = LHAPDF::numberPDF() + 1;
//     fNPDFMembers = nmembers;
//     //    cout << "PDFset has " << nmembers << " members\n";
      
//     //don't use default pdf numbers. use values of pdf member 0 as the default
//     fPDFWeights[0] = 1;
//     LHAPDF::usePDFMember(0);
//     pdf1 =  LHAPDF::xfx(x1, Q, id1)/x1;
//     pdf2 = LHAPDF::xfx(x2, Q, id2)/x2;   
      
//     for (UInt_t i=1; i<nmembers; ++i) {
//       LHAPDF::usePDFMember(i);
//       Double_t newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
//       Double_t newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
//       Double_t TheWeight = newpdf1/pdf1*newpdf2/pdf2;
        
// //             cout << i << " --> " << newpdf1 << " " << newpdf2 << " | " 
// //                  << pdf1 << " "   << pdf2 << " | "
// //                  << x1   << " "   << x2   << " | "
// //                  << id1  << " "   << id2  << " | "
// //                  << Q    << " : " <<  TheWeight << endl;
        
//       fPDFWeights[i] = TheWeight;
//     } 
//   } // end if comp



  //Fill GenTau Info
  fNGenTaus = 0;

  for (UInt_t i = 0; i < fParticles->GetEntries() ; ++i) {
    if (fParticles->At(i)->Is(MCParticle::kTau) && fParticles->At(i)->Status() == 2) {
      const MCParticle *tau = fParticles->At(i);
      const MCParticle *tauNeutrino = fParticles->At(i)->FindDaughter(MCParticle::kTauNu);
      const MCParticle *electron = fParticles->At(i)->FindDaughter(MCParticle::kEl);
      const MCParticle *electronNeutrino = fParticles->At(i)->FindDaughter(MCParticle::kElNu);
      const MCParticle *muon = fParticles->At(i)->FindDaughter(MCParticle::kMu);
      const MCParticle *muonNeutrino = fParticles->At(i)->FindDaughter(MCParticle::kMuNu);

      if (tauNeutrino && ( (electron && electronNeutrino) || (muon && muonNeutrino))) {

 //        cout << i << " : " << fParticles->At(i)->PdgId() << " " << fParticles->At(i)->Pt() << " " << fParticles->At(i)->Eta() << " " << fParticles->At(i)->Phi() << " " << fParticles->At(i)->Mass() << endl;
//         cout << "tauNu: " << tauNeutrino->PdgId() << " " << tauNeutrino->Pt() << " " << tauNeutrino->Eta() << " " << tauNeutrino->Phi() << " " << tauNeutrino->Mass() << endl;
//         if (electron && electronNeutrino) {
//           cout << "electron: " << electron->PdgId() << " " << electron->Pt() << " " << electron->Eta() << " " << electron->Phi() << " " << electron->Mass() << endl;
//           cout << "electronNu: " << electronNeutrino->PdgId() << " " << electronNeutrino->Pt() << " " << electronNeutrino->Eta() << " " << electronNeutrino->Phi() << " " << electronNeutrino->Mass() << endl;
//         }
        
//         if (muon && muonNeutrino) {
//           cout << "muon: " << muon->PdgId() << " " << muon->Pt() << " " << muon->Eta() << " " << muon->Phi() << " " << muon->Mass() << endl;
//           cout << "muonNu: " << muonNeutrino->PdgId() << " " << muonNeutrino->Pt() << " " << muonNeutrino->Eta() << " " << muonNeutrino->Phi() << " " << muonNeutrino->Mass() << endl;
//         }

//         cout << "GENTABLE\n";
//         GeneratorTools::PrintHepMCTable(fParticles, true, -1); 
        
        //we don't find both: sanity check
        if (!((electron && electronNeutrino) && (muon && muonNeutrino))) {

        //Found the leptonic Tau
        fGenTauMass[fNGenTaus] = tau->Mass();
        fGenTauPt[fNGenTaus]   = tau->Pt();
        fGenTauEta[fNGenTaus]  = tau->Eta();
        fGenTauPhi[fNGenTaus]  = tau->Phi();

        fGenTauDecayTauNeutrinoType[fNGenTaus] = tauNeutrino->PdgId();
        fGenTauDecayTauNeutrinoPt[fNGenTaus]   = tauNeutrino->Pt();
        fGenTauDecayTauNeutrinoEta[fNGenTaus]  = tauNeutrino->Eta();
        fGenTauDecayTauNeutrinoPhi[fNGenTaus]  = tauNeutrino->Phi();


        if (muon && muonNeutrino) {
          fGenTauDecayLeptonType[fNGenTaus] = muon->PdgId();
          fGenTauDecayLeptonPt[fNGenTaus]   = muon->Pt();
          fGenTauDecayLeptonEta[fNGenTaus]  = muon->Eta();
          fGenTauDecayLeptonPhi[fNGenTaus]  = muon->Phi();
          fGenTauDecayLeptonNeutrinoType[fNGenTaus] = muonNeutrino->PdgId();
          fGenTauDecayLeptonNeutrinoPt[fNGenTaus]   = muonNeutrino->Pt();
          fGenTauDecayLeptonNeutrinoEta[fNGenTaus]  = muonNeutrino->Eta();
          fGenTauDecayLeptonNeutrinoPhi[fNGenTaus]  = muonNeutrino->Phi();
        }

        if (electron && electronNeutrino) {
          fGenTauDecayLeptonType[fNGenTaus] = electron->PdgId();
          fGenTauDecayLeptonPt[fNGenTaus]   = electron->Pt();
          fGenTauDecayLeptonEta[fNGenTaus]  = electron->Eta();
          fGenTauDecayLeptonPhi[fNGenTaus]  = electron->Phi();
          fGenTauDecayLeptonNeutrinoType[fNGenTaus] = electronNeutrino->PdgId();
          fGenTauDecayLeptonNeutrinoPt[fNGenTaus]   = electronNeutrino->Pt();
          fGenTauDecayLeptonNeutrinoEta[fNGenTaus]  = electronNeutrino->Eta();
          fGenTauDecayLeptonNeutrinoPhi[fNGenTaus]  = electronNeutrino->Phi();
        }

        fNGenTaus++;

        } else {
          cout << "ERROR: Found both electron and muon from tau decay!\n";
        }

      } //found leptonic tau decay
    } //if it's a tau
  } //loop over gen particles



  //Fill Electrons and Muons and MET for Tau study
  fNLeptons = 0;
  for (UInt_t i = 0; i < CleanMuons->GetEntries() ; ++i) {
    fLeptonType[fNLeptons] = 13 * Int_t(CleanMuons->At(i)->Charge());
    //cout << CleanMuons->At(i)->Charge() << " " << fLeptonType[fNLeptons] << endl;
    fLeptonPt[fNLeptons] = CleanMuons->At(i)->Pt();
    fLeptonEta[fNLeptons] = CleanMuons->At(i)->Eta();
    fLeptonPhi[fNLeptons] = CleanMuons->At(i)->Phi();
    fNLeptons++;
  }
  for (UInt_t i = 0; i < CleanElectrons->GetEntries() ; ++i) {
    fLeptonType[fNLeptons] = 11 * Int_t(CleanElectrons->At(i)->Charge());
    //cout << CleanElectrons->At(i)->Charge() << " " << fLeptonType[fNLeptons] << endl;
    fLeptonPt[fNLeptons] = CleanElectrons->At(i)->Pt();
    fLeptonEta[fNLeptons] = CleanElectrons->At(i)->Eta();
    fLeptonPhi[fNLeptons] = CleanElectrons->At(i)->Phi();
    fNLeptons++;
  }
  fNtuplePFMet    = fPFMet->At(0)->Pt();
  fNtuplePFMetPhi = fPFMet->At(0)->Phi();

  //If fill gen only, then skip the rest of the event.
  fEventInfo.eventweight  = fMCEvtInfo->Weight();
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();

  if (fNGenTaus >= 2) {
    fEventTree->Fill();
  }

}

//--------------------------------------------------------------------------------------------------
void HwwGenNtuplerMod::FillGenInfo(const MCParticleCol *GenLeptons, const MCParticleCol *GenNeutrinos, 
                                const MCParticleCol *GenBosons)
{
  // N.B. TGenInfo::nGenPart and TGenInfo::nGenCh is not filled here
  
//   cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
//   GeneratorTools::PrintHepMCTable(fParticles, kTRUE, -1);

  fGenInfo.id_1     = fMCEvtInfo->Id1();
  fGenInfo.id_2     = fMCEvtInfo->Id2();
  fGenInfo.x_1      = fMCEvtInfo->X1();
  fGenInfo.x_2      = fMCEvtInfo->X2();
  fGenInfo.pdf_1    = fMCEvtInfo->Pdf1();
  fGenInfo.pdf_2    = fMCEvtInfo->Pdf2();
  fGenInfo.scalePdf = fMCEvtInfo->ScalePdf();
  fGenInfo.scale    = fMCEvtInfo->Scale();
  fGenInfo.weight   = fMCEvtInfo->Weight();
  
  fGenInfo.mass = 0;
  fGenInfo.pt = 0;
  fGenInfo.y = 0;
  fGenInfo.phi = 0;
  fGenInfo.pt_1 = 0;
  fGenInfo.eta_1 = 0;
  fGenInfo.phi_1 = 0;
  fGenInfo.pt_2 = 0;
  fGenInfo.eta_2 = 0;
  fGenInfo.phi_2 = 0;
  fGenInfo.met = 0;
  fGenInfo.vmass_1 = 0;
  fGenInfo.vpt_1 = 0;
  fGenInfo.vy_1 = 0;
  fGenInfo.veta_1 = 0;
  fGenInfo.vphi_1 = 0;
  fGenInfo.vmass_2 = 0;
  fGenInfo.vpt_2 = 0;
  fGenInfo.vy_2 = 0;
  fGenInfo.veta_2 = 0;
  fGenInfo.vphi_2 = 0;
  fGenInfo.ptBosonSystem = 0;

  if (GenLeptons->GetEntries() > 0) {
    fGenInfo.pt_1 = GenLeptons->At(0)->Pt();
    fGenInfo.eta_1 = GenLeptons->At(0)->Eta();
    fGenInfo.phi_1 = GenLeptons->At(0)->Phi();
   }
  if (GenLeptons->GetEntries() > 1) {
    fGenInfo.pt_2 = GenLeptons->At(1)->Pt();
    fGenInfo.eta_2 = GenLeptons->At(1)->Eta();
    fGenInfo.phi_2 = GenLeptons->At(1)->Phi();
  }

  if (GenNeutrinos->GetEntries() >= 2) {
    CompositeParticle *neutrinoSystem = new CompositeParticle;
    neutrinoSystem->AddDaughter(GenNeutrinos->At(0));
    neutrinoSystem->AddDaughter(GenNeutrinos->At(1));
    fGenInfo.met = neutrinoSystem->Pt();
    fGenInfo.metPhi = neutrinoSystem->Phi();
//     cout << "Neutrino 1 : " << GenNeutrinos->At(0)->Pt() << " " << GenNeutrinos->At(0)->Eta() << " " << GenNeutrinos->At(0)->Phi() << endl;
//     cout << "Neutrino 2 : " << GenNeutrinos->At(1)->Pt() << " " << GenNeutrinos->At(1)->Eta() << " " << GenNeutrinos->At(1)->Phi() << endl;
//     cout << "Met: " << neutrinoSystem->Pt() << " " << neutrinoSystem->Phi() << endl;
    delete neutrinoSystem;
  }

  const MCParticle *GenWBoson1 = 0;
  const MCParticle *GenWBoson2 = 0;
  const MCParticle *GenHiggsBoson = 0;
  for (UInt_t i=0; i< GenBosons->GetEntries(); ++i) {
    if (GenBosons->At(i)->PdgId() == 24) {
      if (!GenWBoson1) {
        GenWBoson1 = GenBosons->At(i);       
      } 
    }
    if (GenBosons->At(i)->PdgId() == -24) {
      if (!GenWBoson2) {
        GenWBoson2 = GenBosons->At(i);       
      }
    }
    if (GenBosons->At(i)->PdgId() == 25) {
      GenHiggsBoson = GenBosons->At(i);
    }
  }

//   if (GenWBoson1) {
//     cout << "W1 : " << GenWBoson1->Pt() << " " << GenWBoson1->Eta() << " " << GenWBoson1->Phi() << " : " << GenWBoson1->PdgId() << endl;
//   } 
//   if (GenWBoson2) {
//     cout << "W2 : " << GenWBoson2->Pt() << " " << GenWBoson2->Eta() << " " << GenWBoson2->Phi() << " : " << GenWBoson2->PdgId() << endl;
//   } 
//   if (GenHiggsBoson) {
//     cout << "Higgs : " << GenHiggsBoson->Pt() << " " << GenHiggsBoson->Eta() << " " << GenHiggsBoson->Phi() << " : " << GenHiggsBoson->PdgId() << endl;
//   } 


  if (GenWBoson1) {
    fGenInfo.vmass_1 = GenWBoson1->Mass();
    fGenInfo.vpt_1 = GenWBoson1->Pt();
    fGenInfo.vy_1 = GenWBoson1->Rapidity();
    fGenInfo.veta_1 = GenWBoson1->Eta();
    fGenInfo.vphi_1 = GenWBoson1->Phi();
  }
  if (GenWBoson2) {
    fGenInfo.vmass_2 = GenWBoson2->Mass();
    fGenInfo.vpt_2 = GenWBoson2->Pt();
    fGenInfo.vy_2 = GenWBoson2->Rapidity();
    fGenInfo.veta_2 = GenWBoson2->Eta();
    fGenInfo.vphi_2 = GenWBoson2->Phi();

    CompositeParticle *BosonSystem = new CompositeParticle;
    BosonSystem->AddDaughter(GenWBoson1);
    BosonSystem->AddDaughter(GenWBoson2);
//     cout << "WW System: " << BosonSystem->Pt() << " " << BosonSystem->Eta() << " " << BosonSystem->Phi() << endl;
    fGenInfo.ptBosonSystem = BosonSystem->Pt();
    fGenInfo.mass = BosonSystem->Mass();
    fGenInfo.pt = BosonSystem->Pt();
    fGenInfo.y = BosonSystem->Rapidity();
    fGenInfo.phi = BosonSystem->Phi();
    fGenInfo.ptBosonSystem = BosonSystem->Pt();
    delete BosonSystem;
  }

  //take higgs pt from the higgs particle directly
  if (GenHiggsBoson) {
    fGenInfo.mass = GenHiggsBoson->Mass();
    fGenInfo.pt = GenHiggsBoson->Pt();
    fGenInfo.y = GenHiggsBoson->Rapidity();
    fGenInfo.phi = GenHiggsBoson->Phi();
  }


//   if (GenLeptons->GetEntries() > 0) 
//     cout << "Lepton1: " << fGenInfo.pt_1 << " " << fGenInfo.eta_1 << " " << fGenInfo.phi_1 << endl;
//   if (GenLeptons->GetEntries() > 1)
//     cout << "Lepton2: " << fGenInfo.pt_2 << " " << fGenInfo.eta_2 << " " << fGenInfo.phi_2 << endl;
//   cout << "BosonSystemPt: " << fGenInfo.ptBosonSystem  << endl;


  if (fGenJets) {
    Int_t nGenJets = 0;
    Int_t genJetsIndex = 0;
    fGenInfo.jetpt_1 = 0;
    fGenInfo.jetpt_2 = 0;
    fGenInfo.jetpt_3 = 0;
    fGenInfo.jetpt_4 = 0;
    for (UInt_t i=0; i<fGenJets->GetEntries(); ++i) {
      
      if (fabs(fGenJets->At(i)->Eta()) < 5.0 && fGenJets->At(i)->Pt() > 10.0  
          && (GenLeptons->GetEntries() < 1 || MathUtils::DeltaR(*fGenJets->At(i), *GenLeptons->At(0)) > 0.3)
          && (GenLeptons->GetEntries() < 2 || MathUtils::DeltaR(*fGenJets->At(i), *GenLeptons->At(1)) > 0.3)
        ) {

//         cout << "GenJet " << i << " : " << fGenJets->At(i)->Pt() << " " << fGenJets->At(i)->Eta() << " " << fGenJets->At(i)->Phi() << endl;
//         for (UInt_t k=0; k < fParticles->GetEntries() ; ++k) {
//           if (MathUtils::DeltaR(*fGenJets->At(i), *fParticles->At(k)) < 0.5) {
//             cout << "GenJet Particle: " << k << " : " << fParticles->At(k)->PdgId() << " " << fParticles->At(k)->Status() << " " << fParticles->At(k)->Pt() << " " << fParticles->At(k)->Eta() << " " << fParticles->At(k)->Phi() << endl;
//           }
//         }

        fGenInfo.jetpt_1 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 1) fGenInfo.jetpt_2 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 2) fGenInfo.jetpt_3 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 3) fGenInfo.jetpt_4 = fGenJets->At(i)->Pt();
        if (fGenJets->At(i)->Pt() > 30.0) nGenJets++;
        genJetsIndex++;        
      }
      fGenInfo.nGenJets = nGenJets;   
    }
  
  }
}
