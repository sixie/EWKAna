#include "EWKAna/Ntupler/interface/HwwNtuplerMod.hh"
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
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Utils/interface/GeneratorTools.h"
#include "MitPhysics/Utils/interface/JetTools.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"
#include "MitPhysics/Utils/interface/MetTools.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include "LHAPDF/LHAPDF.h"
#include <vector>

//#define __DATA_36X__

using namespace mithep;

ClassImp(mithep::HwwNtuplerMod)

HwwNtuplerMod::HwwNtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
  fOutputName    ("ntuple.root"),
  fPDFName       ("cteq66.LHgrid"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fElectronName  (Names::gkElectronBrn),
  fMuonName     (Names::gkMuonBrn),
  fTrackName     (Names::gkTrackBrn),
  fCaloTowerName (Names::gkCaloTowerBrn),
  fPrimVtxName   (Names::gkPVBeamSpotBrn),
  fBeamSpotName  ("BeamSpot"),
  fCaloJetName   ("AKt5Jets"),
  fCleanJetsName ("NoDefaultNameSet"),
  fCleanJetsNoPtCutName("NoDefaultNameSet"),
  fTrigMaskName  (Names::gkHltBitBrn),
  fCaloMetName   (Names::gkCaloMetBrn),
  fTCMetName     ("TCMet"),
  fPFMetName     ("PFMet"),
  fConversionName(Names::gkMvfConversionBrn),
  fParticles     (0),
  fMCEvtInfo     (0),
  fGenJets       (0),
  fElectrons     (0),
  fTracks        (0),
  fPrimVerts     (0),
  fBeamSpot      (0),
  fCaloJets      (0),
  fTrackJets     (0),
  fPFJets        (0),
  fPFCandidates  (0),
  fPhotons       (0),
  fTrigMask      (0),
  fCaloMet       (0),
  fTCMet         (0),
  fPFMet         (0),
  fConversions   (0),  
  fUseGen        (kTRUE),
  fPrintTable    (kFALSE),
  fSkipIfHLTFail (kFALSE),
  fFakeRateSkim  (kFALSE),
  fFillGenOnly   (kFALSE),
  fComputePDFWeights(kFALSE),
  fReadPileupInfo(kTRUE),
  fSCEtMin       (5),
  fSCEtMax       (14000),
  fSCEtaMin      (-3),
  fSCEtaMax      (3),
  fMuonPtMin     (5),
  fMuonPtMax     (14000),
  fMuonEtaMin    (-3),
  fMuonEtaMax    (3),
  fMassMin       (0),
  fMassMax       (1000),
  fJetPtMin      (10),
  fPhotonEtMin   (10), 
  fEventTree     (0),
  fFSRMode       (0)
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
HwwNtuplerMod::~HwwNtuplerMod()
{
  // Destructor
}
	//--------------------------------------------------------------------------------------------------      
void HwwNtuplerMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::BeginRun()
{
//   if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); GetL1AlgoTable()->Print(); }
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,      fParticles); 
  ReqBranch(fMCEvtInfoName, fMCEvtInfo);
  ReqBranch(Names::gkGenJetBrn , fGenJets); 
  ReqBranch(fElectronName,  fElectrons);
  ReqBranch(fMuonName,      fMuons);
  ReqBranch(fTrackName,     fTracks);
  ReqBranch(fCaloTowerName, fCaloTowers);
  ReqBranch(Names::gkPFCandidatesBrn,     fPFCandidates);
  ReqBranch(fBeamSpotName,  fBeamSpot);
  ReqBranch(fTrigMaskName,  fTrigMask);
  ReqBranch("L1AlgoBitsBeforeMask",  fL1TrigMask);
  ReqBranch(fTCMetName,     fTCMet);
  ReqBranch(fPFMetName,     fPFMet);
  ReqBranch(fConversionName,fConversions);
  ReqBranch(fCaloJetName,   fCaloJets);
  ReqBranch(Names::gkPhotonBrn,    fPhotons);
  ReqBranch(Names::gkPFCandidatesBrn, fPFCandidates);
  cout << "ReadPileup: " << fReadPileupInfo << endl;
  cout << "UseGen : " << fUseGen << endl;
  if (fReadPileupInfo) {
    cout << "useGen: " << fUseGen << endl;
    if (Bool_t(fUseGen)) ReqBranch(Names::gkPileupInfoBrn,   fPileupInfo);
    ReqBranch(Names::gkPileupEnergyDensityBrn, fPileupEnergyDensity);
  }

  //
  // Set up arrays
  //
  fElectronArr   = new TClonesArray("mithep::TElectron");   assert(fElectronArr);
  fMuonArr       = new TClonesArray("mithep::TMuon");       assert(fMuonArr);
  fPFJetArr      = new TClonesArray("mithep::TJet");        assert(fPFJetArr);
  fPhotonArr     = new TClonesArray("mithep::TPhoton");     assert(fPhotonArr);

  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");

  //
  // Initialize data trees and structs
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(fUseGen)
    fEventTree->Branch("Gen",&fGenInfo);
  
  fEventTree->Branch("Electron",   &fElectronArr);
  fEventTree->Branch("Muon",       &fMuonArr);
  fEventTree->Branch("PFJet",      &fPFJetArr);
  fEventTree->Branch("Photon",     &fPhotonArr);

  if (fComputePDFWeights) {
 //    gSystem->Setenv("LHAPATH","/afs/cern.ch/cms/sw/slc4_ia32_gcc345/external/lhapdf/5.6.0-cms4/share/lhapdf/PDFsets/");
    LHAPDF::setVerbosity(LHAPDF::SILENT);
    LHAPDF::initPDFSet(fPDFName.Data());
    LHAPDF::getDescription();
    
    //Branches for PDFset weights
    fEventTree->Branch("NPDFMembers",&fNPDFMembers,"NPDFMembers/I");  
    fEventTree->Branch("PDFWeights",fPDFWeights,"PDFWeights[NPDFMembers]/F");
  }

//   AddOutput(fEventTree);

  fMuonTools = new MuonTools();
  

}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();

  fOutputFile->cd();
  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fElectronArr;
  delete fMuonArr;
  delete fPFJetArr;
  delete fPhotonArr;

}  

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::Process()
{
  //
  // Load branches
  //
  if(fUseGen) LoadBranch(fPartName);
  if(fUseGen) LoadBranch(fMCEvtInfoName);
  if(fUseGen) LoadBranch(Names::gkGenJetBrn);
  if (!fFillGenOnly) {
    LoadBranch(fElectronName);
    LoadBranch(fMuonName);
    LoadBranch(fTrackName);
    LoadBranch(fCaloTowerName);
    LoadBranch(Names::gkPFCandidatesBrn);
    LoadBranch(fBeamSpotName);
    LoadBranch(fTrigMaskName);
    LoadBranch("L1AlgoBitsBeforeMask");
    LoadBranch(fTCMetName);
    LoadBranch(fPFMetName); 
    LoadBranch(fConversionName);
    LoadBranch(fCaloJetName);
    LoadBranch(Names::gkPhotonBrn);
    LoadBranch(Names::gkPFCandidatesBrn);
    if (fReadPileupInfo) {
      if (Bool_t(fUseGen)) LoadBranch(Names::gkPileupInfoBrn);
      LoadBranch(Names::gkPileupEnergyDensityBrn);
    }
  }


 //***********************************************************************************************
  //Some MC Specific Requirements
  //***********************************************************************************************
  const MCParticleCol *GenLeptons = 0;
  GenLeptons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCLeptonsName,0);
  const MCParticleCol *GenNeutrinos = GetObjThisEvt<MCParticleCol>(ModNames::gkMCNeutrinosName,0);
  const MCParticleCol *GenTaus = GetObjThisEvt<MCParticleCol>(ModNames::gkMCTausName,0);
  const MCParticleCol *GenPhotons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCPhotonsName,0);  
  const MCParticleCol *GenBosons = GetObjThisEvt<MCParticleCol>(ModNames::gkMCBosonsName,0);

  //exclude Wgamma , Zgamma
  Bool_t FoundZ = false;
  Bool_t FoundZGamma = false;
  Bool_t FoundWGammaEvent = false;
  Bool_t FoundW = false;


  if (fUseGen) {

    for (UInt_t i=0; i<GenBosons->GetEntries(); i++) {  
      if (GenBosons->At(i)->AbsPdgId() == 23) {
        FoundZ = true;
      }
    }
    
    for (UInt_t i=0; i<GenPhotons->GetEntries(); i++) {  
      if (GenPhotons->At(i)->Pt() > 10.0) {
        
        //ISR Photon
        if ( (GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->IsParton())
            || (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
                && GenPhotons->At(i)->Mother()->Status() ==3  
                && GenPhotons->At(i)->Mother()->Mother() 
                && GenPhotons->At(i)->Mother()->Mother()->IsParton())
          ) {
          FoundZGamma = true;    
        }
        
        //Pythia FSR
        if (GenPhotons->At(i)->Mother() && (GenPhotons->At(i)->Mother()->Status() == 3 || GenPhotons->At(i)->Mother()->Status() == 2)
            && (GenPhotons->At(i)->Mother()->AbsPdgId() == 11 
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 13
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 15)
          ) {
          CompositeParticle *object = new CompositeParticle();
          object->AddDaughter(GenPhotons->At(i));
          object->AddDaughter(GenPhotons->At(i)->Mother());
          if(object->Mass() > 1.0) FoundZGamma = true;
          delete object;          
        }      
      }
    }
    
    //For WJets sample: Remove FSR W+gamma Events.
   
    for (UInt_t i=0; i<GenBosons->GetEntries(); i++) {  
      if (GenBosons->At(i)->AbsPdgId() == 24) {
        FoundW = true;
      }
    }

    for (UInt_t i=0; i<GenPhotons->GetEntries(); i++) {  
      if (GenPhotons->At(i)->Pt() > 10.0) {
        
        //ISR Photon
        if ((GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->IsParton())
            || (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
                && GenPhotons->At(i)->Mother()->Status() ==3  
                && GenPhotons->At(i)->Mother()->Mother() 
                && GenPhotons->At(i)->Mother()->Mother()->IsParton())
          ) {
          FoundWGammaEvent = true;    
        }
        
        //WWgamma vertex
        if ((GenPhotons->At(i)->Mother() && GenPhotons->At(i)->Mother()->AbsPdgId() == 24) 
            || 
            (GenPhotons->At(i)->Mother()->AbsPdgId() == 22 
             && GenPhotons->At(i)->Mother()->Status() == 3
             && GenPhotons->At(i)->Mother()->Mother()
             && GenPhotons->At(i)->Mother()->Mother()->AbsPdgId() == 24
              )
          ) {
          FoundWGammaEvent = true;
        }
        
        //Pythia FSR
        if (GenPhotons->At(i)->Mother() && (GenPhotons->At(i)->Mother()->Status() == 3 || GenPhotons->At(i)->Mother()->Status() == 2)
            && (GenPhotons->At(i)->Mother()->AbsPdgId() == 11 
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 13
                || GenPhotons->At(i)->Mother()->AbsPdgId() == 15)
          ) {
          CompositeParticle *object = new CompositeParticle();
          object->AddDaughter(GenPhotons->At(i));
          object->AddDaughter(GenPhotons->At(i)->Mother());
          if(object->Mass() > 1.0) FoundWGammaEvent = true;
          delete object;          
        }
      }
    } //for all gen photons   
   
    //Fill Gen Info
    FillGenInfo(GenLeptons, GenNeutrinos, GenBosons);  // fill the data structure



  }

  fPrintDebug = kFALSE;
  if ( (GetEventHeader()->RunNum() == 162926 && GetEventHeader()->LumiSec() == 647 && GetEventHeader()->EvtNum() == 389987568 )
       || 
       (GetEventHeader()->RunNum() == 162926 && GetEventHeader()->LumiSec() == 647 && GetEventHeader()->EvtNum() == 390151968 )


       || GetEventHeader()->EvtNum() == 7755124
       || GetEventHeader()->EvtNum() == 19668349
       || GetEventHeader()->EvtNum() == 215345879
       || GetEventHeader()->EvtNum() == 288792650
       
       || GetEventHeader()->EvtNum() == 6086587
       || GetEventHeader()->EvtNum() == 71683939
       || GetEventHeader()->EvtNum() == 251947246
       || GetEventHeader()->EvtNum() == 295317759



    ) {
    fPrintDebug = kTRUE;    
  }
  if (fPrintDebug) {
    cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
  }

  //************************************************************************************************
  //Fill PDF information
  //************************************************************************************************
  if (fComputePDFWeights && fUseGen) {
    //cout << "Event " << fMCEvtInfo->Weight() << " " << fMCEvtInfo->Id1() << " " << fMCEvtInfo->Id2() << " " << fMCEvtInfo->X1() << " " << fMCEvtInfo->X2() << endl;    
    Double_t Q    = fMCEvtInfo->Scale();
    Int_t    id1  = fMCEvtInfo->Id1();
    Double_t x1   = fMCEvtInfo->X1();
    Double_t pdf1 = fMCEvtInfo->Pdf1();
    Int_t    id2  = fMCEvtInfo->Id2();
    Double_t x2   = fMCEvtInfo->X2();
    Double_t pdf2 = fMCEvtInfo->Pdf2();
    
    UInt_t nmembers = LHAPDF::numberPDF() + 1;
    fNPDFMembers = nmembers;
    //    cout << "PDFset has " << nmembers << " members\n";
      
    //don't use default pdf numbers. use values of pdf member 0 as the default
    fPDFWeights[0] = 1;
    LHAPDF::usePDFMember(0);
    pdf1 =  LHAPDF::xfx(x1, Q, id1)/x1;
    pdf2 = LHAPDF::xfx(x2, Q, id2)/x2;   
      
    for (UInt_t i=1; i<nmembers; ++i) {
      LHAPDF::usePDFMember(i);
      Double_t newpdf1 = LHAPDF::xfx(x1, Q, id1)/x1;
      Double_t newpdf2 = LHAPDF::xfx(x2, Q, id2)/x2;
      Double_t TheWeight = newpdf1/pdf1*newpdf2/pdf2;
        
//             cout << i << " --> " << newpdf1 << " " << newpdf2 << " | " 
//                  << pdf1 << " "   << pdf2 << " | "
//                  << x1   << " "   << x2   << " | "
//                  << id1  << " "   << id2  << " | "
//                  << Q    << " : " <<  TheWeight << endl;
        
      fPDFWeights[i] = TheWeight;
    } 
  } // end if comp

  //If fill gen only, then skip the rest of the event.
  if (fFillGenOnly){ 
    fEventInfo.eventweight  = fMCEvtInfo->Weight();
    fEventInfo.runNum       = GetEventHeader()->RunNum();
    fEventInfo.evtNum       = GetEventHeader()->EvtNum();
    fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
    fEventTree->Fill();

    return;
  }


  //************************************************************************************************
  //Obtain all the good objects from the event cleaning module
  //************************************************************************************************
//   ObjArray<Muon> *CleanMuons = dynamic_cast<ObjArray<Muon>* >(FindObjThisEvt(ModNames::gkCleanMuonsName));
//   ObjArray<Electron> *CleanElectrons = dynamic_cast<ObjArray<Electron>* >(FindObjThisEvt(ModNames::gkCleanElectronsName));
  ObjArray<Photon> *CleanPhotons = dynamic_cast<ObjArray<Photon>* >(FindObjThisEvt(ModNames::gkCleanPhotonsName));

//   JetOArr *CleanJets    = GetObjThisEvt<JetOArr>(fCleanJetsName);
  ObjArray<Jet> *CleanJetsNoPtCut = dynamic_cast<ObjArray<Jet>* >
    (FindObjThisEvt(fCleanJetsNoPtCutName.Data()));



  //************************************************************************************************
  //Fake Rate Skim
  //Keep only events with one denominator electron (v2,v3,or v4) and a jet with pt > 15
  //************************************************************************************************
  if (fFakeRateSkim) {
    Bool_t passSkim = kFALSE;
    ElectronOArr *electronsIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsIsolated");
    ElectronOArr *electronsVBTF80NonIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsVBTF80NonIsolated");
    ElectronOArr *electronsVBTF90PartiallyIsolated    = GetObjThisEvt<ElectronOArr>("ElectronsVBTF90PartiallyIsolated");
    MuonOArr *muonsDenominator    = GetObjThisEvt<MuonOArr>("DenominatorMuons");

    Bool_t foundJet = kFALSE;
    if (CleanJetsNoPtCut) {
      for (UInt_t k=0; k<CleanJetsNoPtCut->GetEntries() ; ++k) {

        if (CleanJetsNoPtCut->At(k)->Pt() < 15.0) continue;

        Bool_t hasOverlap = kFALSE;
        if (electronsIsolated) {
          for (UInt_t l = 0; l < electronsIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (electronsVBTF80NonIsolated) {
          for (UInt_t l = 0; l < electronsVBTF80NonIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsVBTF80NonIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (electronsVBTF90PartiallyIsolated) {
          for (UInt_t l = 0; l < electronsVBTF90PartiallyIsolated->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*electronsVBTF90PartiallyIsolated->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }
        if (!hasOverlap) {
          foundJet = kTRUE;
          break;
        }
      }
    }

    Int_t NRecoElectrons = 0;
    for(UInt_t k=0; k<fElectrons->GetEntries(); ++k) {
      if (fElectrons->At(k)->Pt() > 10.0) NRecoElectrons++;
    } 

    //if we didn't find an electron and a jet then we don't save the event
    if ( NRecoElectrons + electronsIsolated->GetEntries() + 
         electronsVBTF80NonIsolated->GetEntries() + 
         electronsVBTF90PartiallyIsolated->GetEntries() >= 1 && foundJet         
      ) {
      passSkim = kTRUE;
    }


    if (CleanJetsNoPtCut) {
      for (UInt_t k=0; k<CleanJetsNoPtCut->GetEntries() ; ++k) {
        
        if (CleanJetsNoPtCut->At(k)->Pt() < 15.0) continue;
        
        Bool_t hasOverlap = kFALSE;
        if (muonsDenominator) {
          for (UInt_t l = 0; l < muonsDenominator->GetEntries() ; ++l) {
            if (MathUtils::DeltaR(*muonsDenominator->At(l),*CleanJetsNoPtCut->At(k)) < 0.3)
              hasOverlap = kTRUE;
          }
        }        
        if (!hasOverlap) {
          foundJet = kTRUE;
          break;
        }
      }
    }
    
    //if we didn't find a muon and a jet then we don't save the event
    if (muonsDenominator->GetEntries() >= 1 && foundJet) {
      passSkim = kTRUE;
    }

    if (!passSkim) {
      return;
    }

  }




  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  //
  ULong_t trigbits=0;
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();

    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
       
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }

      if (fPrintDebug)
      {
        cout << "Trig : " << trigname->GetName() << " " << fTrigMask->At(trigname->Id())  << " "
	     << fTriggerIdsv[itrig] << " ::: " << trigbits << " " 
	     << endl;
      }
  
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0))
    return;
  
  ULong_t l1trigbits = 0;
  if (fL1TrigMask) {
    const TriggerTable *l1Table = GetL1AlgoTable();
    assert(l1Table);
    for(UInt_t itrig=0; itrig<fL1TriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = l1Table->Get(fL1TriggerNamesv[itrig].Data());
      if(!trigname) continue;


      if (fPrintDebug)
      {
        cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
        cout << "L1Trig : " << trigname->GetName() << " " << trigname->Id() << " ";
        cout << fL1TrigMask->At(trigname->Id()) << endl;
      }
      
      if(fL1TrigMask->At(trigname->Id())) { l1trigbits |= fL1TriggerIdsv[itrig]; }
    }  
  }

  



  IncNEventsProcessed();
  
  fElectronArr->Clear();
  fMuonArr->Clear();
  fPFJetArr->Clear();
  fPhotonArr->Clear();
   

  //
  // Get beam spot. If no beam spot information is available, default the coordinates to 99999
  //
  Double_t bsx=99999, bsy=99999, bsz=99999;
  if(fBeamSpot) {
    if(fBeamSpot->GetEntries() > 1) 
      std::cout << "********** More than 1 beam spot! **********" << std::endl;
    const BeamSpot *bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }


  //
  // Get primary vertex (for corrected d0)
  // Take the first primary vertex listed (should be ordered by sum-pT as in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  fPrimVerts = GetObjThisEvt<VertexOArr>(ModNames::gkGoodVertexesName);
  const Vertex* primaryVertex = 0;
  Bool_t hasGoodPV = kFALSE;   
  if (fPrimVerts->At(0)) {
    primaryVertex = fPrimVerts->At(0);
    hasGoodPV = kTRUE;
  } 


  for(UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    if (fPrintDebug) {
      cout << "Vertex " << i << " : " << fPrimVerts->At(i)->Z() << " " << fPrimVerts->At(i)->X() << " " << fPrimVerts->At(i)->Y() << " " 
           <<  fPrimVerts->At(i)->NTracks() << " " 
           << fPrimVerts->At(i)->NTracksFit() << " "
           << fPrimVerts->At(i)->Ndof() << " "
           << fPrimVerts->At(i)->Position().Rho() << " "
           << fPrimVerts->At(i)->Position().Z() << " "
           << endl;
    }
  }


  
  //-----------------------------------------
  //
  // Loop through electrons.
  //
  //-----------------------------------------
  vector<const Electron*> elev;  // array of pointers to preselected electrons ... 
  ElectronTools eleTools;        // helper class for electron ID decisions

  assert(fElectrons);                
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
    
    if (fPrintDebug) {
      cout << "ele " << i << " " << ele->Pt() << " " << ele->Eta() << " " << ele->Phi() << endl;
    }

    if (ele->SCluster()->Et()  < fSCEtMin  || ele->SCluster()->Et()  > fSCEtMax ) continue;
    if (ele->SCluster()->Eta() < fSCEtaMin || ele->SCluster()->Eta() > fSCEtaMax) continue;
    
    if (fPrintDebug) {
      cout << "pass \n";
    }


    elev.push_back(ele);       
  }
  
  for(UInt_t i=0; i<elev.size(); i++) {  

    //********************************************************************************
    //Match MC Truth
    //********************************************************************************
    Int_t isMCMatched = 0;
    Bool_t GenEleMatch = kFALSE;
    Bool_t GenMuMatch = kFALSE;
    Bool_t GenTauMatch = kFALSE;
    Bool_t GenPhotonMatch = kFALSE;
    if (GenLeptons) {
      for (UInt_t l=0; l<GenLeptons->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenLeptons->At(l), *(elev[i])) < 0.3  ) {
          if (GenLeptons->At(l)->AbsPdgId() == 11) {
            GenEleMatch = kTRUE;
          } 
          if (GenLeptons->At(l)->AbsPdgId() == 13) {
            GenMuMatch = kTRUE;
          }           
        }
      }
    }
    if (GenTaus) {
      for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenTaus->At(l), *(elev[i])) < 0.3 ) {
          GenTauMatch = kTRUE;
          break;
        }
      }
    }
    if (GenPhotons) {
      for (UInt_t l=0; l < GenPhotons->GetEntries(); l++) {  
        if (GenPhotons->At(l)->Pt() > 10.0 &&  MathUtils::DeltaR(*GenPhotons->At(l), *(elev[i])) < 0.3 ) {
        
          //ISR Photon
          if ( (GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->IsParton())
               || (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
                   && GenPhotons->At(l)->Mother()->Status() ==3  
                   && GenPhotons->At(l)->Mother()->Mother() 
                   && GenPhotons->At(l)->Mother()->Mother()->IsParton())
            ) {
            GenPhotonMatch = kTRUE;    
          }
        
          //WWgamma vertex
          if ((GenPhotons->At(l)->Mother() && GenPhotons->At(l)->Mother()->AbsPdgId() == 24) 
              || 
              (GenPhotons->At(l)->Mother()->AbsPdgId() == 22 
               && GenPhotons->At(l)->Mother()->Status() == 3
               && GenPhotons->At(l)->Mother()->Mother()
               && GenPhotons->At(l)->Mother()->Mother()->AbsPdgId() == 24
                )
            ) {
            GenPhotonMatch = true;
          }

          //Pythia FSR
          if (GenPhotons->At(l)->Mother() && (GenPhotons->At(l)->Mother()->Status() == 3 || GenPhotons->At(l)->Mother()->Status() == 2)
              && (GenPhotons->At(l)->Mother()->AbsPdgId() == 11 
                  || GenPhotons->At(l)->Mother()->AbsPdgId() == 13
                  || GenPhotons->At(l)->Mother()->AbsPdgId() == 15)
            ) {
            CompositeParticle *object = new CompositeParticle();
            object->AddDaughter(GenPhotons->At(l));
            object->AddDaughter(GenPhotons->At(l)->Mother());
            if(object->Mass() > 1.0) GenPhotonMatch = kTRUE;
            delete object;          
          }      
        }
      }
    }


    if (GenEleMatch) isMCMatched += 1;
    if (GenMuMatch) isMCMatched += 2;
    if (GenTauMatch) isMCMatched += 4;
    if (GenPhotonMatch) isMCMatched += 8;

    
  
    //********************************************************************************
    //Fill Electron
    //********************************************************************************
    FillElectron(elev[i], isMCMatched, primaryVertex);  // fill electron data object    

//      if (elev[i]->Pt() > 10 && !isMCMatched && (elev[i]->TrackIsolationDr03() + elev[i]->EcalRecHitIsoDr03() + elev[i]->HcalTowerSumEtDr03()) / elev[i]->Pt() < 0.1) {
//        cout << "ELE : " << elev[i]->Pt() << " " << elev[i]->Eta() << " " << elev[i]->Phi() << endl;
//        cout << elev[i]->TrackIsolationDr03() << " " << elev[i]->EcalRecHitIsoDr03() << " " << elev[i]->HcalTowerSumEtDr03() << " ";
//        cout << elev[i]->BestTrk()->D0Corrected(primaryVertex) << " " << elev[i]->BestTrk()->DzCorrected(primaryVertex);
//        cout << endl;
       
//        GeneratorTools::PrintNearbyParticles(fParticles, elev[i]->Eta(), elev[i]->Phi(), 0.5);
//        cout << "FULL GEN TABLE\n";
//        GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//      }

  }
  
  
  //-----------------------------------------
  //
  // Loop through muons (and general tracks if desired).
  // If desired, tracks not matched to muon system information but passing 
  // kinematic preselection will be included in the muon TTree. 
  // Tracks from tracker info only will have typeBits=0
  //
  //-----------------------------------------
  vector<const Muon*> muonv;    // array of pointers to preselected muons ... 
  assert(fMuons);
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    
    // Use tracker tracks for kinematics when available
    const Track *muTrk=0;
    if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk(); }
    else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
          
    if((muTrk->Pt()  < fMuonPtMin)  || (muTrk->Pt()  > fMuonPtMax))  continue;  // muon pT cut
    if((muTrk->Eta() < fMuonEtaMin) || (muTrk->Eta() > fMuonEtaMax)) continue;  // muon eta cut
    
    muonv.push_back(mu);        
  }

  for(UInt_t i=0; i<muonv.size(); i++) {    

    //********************************************************************************
    //Match MC Truth
    //********************************************************************************
    Int_t isMCMatched = 0;
    Bool_t GenEleMatch = kFALSE;
    Bool_t GenMuMatch = kFALSE;
    Bool_t GenTauMatch = kFALSE;
    if (GenLeptons) {
      for (UInt_t l=0; l<GenLeptons->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenLeptons->At(l), *(muonv[i])) < 0.3 ) {
          if (GenLeptons->At(l)->AbsPdgId() == 11) {
            GenEleMatch = kTRUE;
          } 
          if (GenLeptons->At(l)->AbsPdgId() == 13) {
            GenMuMatch = kTRUE;
          }           
        }
      }
    }
    if (GenTaus) {
      for (UInt_t l=0; l<GenTaus->GetEntries(); ++l) {
        if (MathUtils::DeltaR(*GenTaus->At(l), *(muonv[i])) < 0.3 && (GenTaus->At(l)->AbsPdgId() == 13 || GenTaus->At(l)->AbsPdgId() == 15)) {
          GenTauMatch = kTRUE;
          break;
        }
      }
    }
    if (GenEleMatch) isMCMatched += 1;
    if (GenMuMatch) isMCMatched += 2;
    if (GenTauMatch) isMCMatched += 4;




    // fill muon data object
    FillMuon(muonv[i], isMCMatched, primaryVertex);

    if (fPrintDebug) {
      cout << "MUON : " << muonv[i]->Pt() << " " << muonv[i]->Eta() << " " << muonv[i]->Phi() << endl;
      cout << muonv[i]->IsoR03SumPt() << " " << muonv[i]->IsoR03EmEt() << " " << muonv[i]->IsoR03HadEt() << " ";
      
      const Track *muTrk=0;
      if(muonv[i]->HasTrackerTrk())         { muTrk = muonv[i]->TrackerTrk(); }
      else if(muonv[i]->HasStandaloneTrk()) { muTrk = muonv[i]->StandaloneTrk(); } 
      cout << muTrk->D0Corrected(*primaryVertex) << " " << muTrk->DzCorrected(*primaryVertex);
      cout << endl;
    }

//       GeneratorTools::PrintNearbyParticles(fParticles, muonv[i]->Eta(), muonv[i]->Phi(), 0.5);
//       cout << "FULL GEN TABLE\n";
//       GeneratorTools::PrintHepMCTable(fParticles,kTRUE,-1);
//     }

  }

  //-----------------------------------------
  //
  //Fill tracks for muon tag and probe
  //
  //-----------------------------------------
  assert(fTracks);
  for(UInt_t i=0; i<fTracks->GetEntries(); ++i) {
    const Track *track = fTracks->At(i);

    if((track->Pt()  < 9.0)  || (track->Pt()  > 7000.0))  continue;  // pT cut
    if((track->Eta() < fMuonEtaMin) || (track->Eta() > fMuonEtaMax)) continue;  // eta cut
      
    // Check that the track is not associated with a muon.
    // If it is, skip to next track...
    Bool_t isMuon = kFALSE;
    for(UInt_t j=0; j<fMuons->GetEntries(); ++j) {
      if(track == (fMuons->At(j)->TrackerTrk())) isMuon = kTRUE;
    }
    if(isMuon) continue;
  
    FillMuon(track, primaryVertex);
  }



  //
  // Loop through jets
  //
  assert(CleanJetsNoPtCut);
  for(UInt_t i=0; i<CleanJetsNoPtCut->GetEntries(); ++i) {
    const PFJet *jet = (PFJet*)CleanJetsNoPtCut->At(i);

    if (fPrintDebug) {
      cout << "Jet " << i << " : " <<  CleanJetsNoPtCut->At(i)->Pt() << " : " << jet->Et() << " " << jet->Pt() << " " << jet->Eta() << " " << jet->Phi() << " : " 
           << jet->L1OffsetCorrectionScale() << " " << jet->L2RelativeCorrectionScale() << " " << jet->L3AbsoluteCorrectionScale() << " "  
           << jet->L4EMFCorrectionScale() << " " << jet->L5FlavorCorrectionScale() << " " << jet->L6LSBCorrectionScale() << " " << jet->L7PartonCorrectionScale() << " " 
           << jet->CombinedCorrectionScale() << " " 
           << " === " << jet->RawMom().Pt() << endl;
    }

    if(CleanJetsNoPtCut->At(i)->RawMom().Pt() > fJetPtMin ) { 
      FillJet(jet, primaryVertex); 
    }
  }


//   //************************************************************************************************
//   // Compute Jettiness Variables
//   //************************************************************************************************
//   // signal 
//   CompositeParticleOArr *Signals = dynamic_cast<CompositeParticleOArr* >(FindObjThisEvt("signalsForJettiness"));
//   // signal jets
//   JetOArr *signalJets = dynamic_cast<JetOArr* >(FindObjThisEvt("signalJetsForJettiness"));
//   // radiation jets
// //   JetOArr *radiationJets = dynamic_cast<JetOArr* >(FindObjThisEvt("radiationJetsForJettiness"));
//   // radiation tracks
//   TrackOArr *radiationTracks = dynamic_cast<TrackOArr* > (FindObjThisEvt("radiationTracksForJettiness"));
//   // pf candidates tracks
//   PFCandidateOArr *radiationPFCandidates = dynamic_cast<PFCandidateOArr* > (FindObjThisEvt("radiationPFCandidatesForJettiness"));

  double jettinessTracks = -999;
  double jettinessPFCands = -999;
  double beamThrustTracks = -999;
  double beamThrustPFCands = -999;
//   if ( Signals->Entries() > 0 ) {
//     CompositeParticle * signal = Signals->At(0);
//     double beta = sqrt(1-signal->Mass()*signal->Mass()/(signal->E()*signal->E()));
//     double costheta = signal->Pz()/signal->P();
//     double rapidity = 0.5*log( (1+beta*costheta) / (1-beta*costheta) );
//     double Q = sqrt(signal->Mass()*signal->Mass()+signal->Pt()*signal->Pt());
//     jettinessTracks = JetTools::NJettiness(radiationTracks,signalJets, Q,rapidity);
//     jettinessPFCands = JetTools::NJettiness(radiationPFCandidates,signalJets, Q,rapidity);
//     beamThrustTracks = JetTools::NJettiness(radiationTracks,signalJets, Q, 0);
//     beamThrustPFCands = JetTools::NJettiness(radiationPFCandidates,signalJets, Q, 0);
//   }
  
  //************************************************************************************************
  // Loop through photons
  //************************************************************************************************
//   assert(fPhotons);  
//   for(UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
  if (CleanPhotons) {
    for(UInt_t i=0; i<CleanPhotons->GetEntries(); ++i) {
      const Photon *pho= fPhotons->At(i);
      if(pho->HasPixelSeed()) continue;  // require NOT pixel seeded
      if(pho->Pt() >  10) { FillPhoton(pho); }
    }
  }

  //************************************************************************************************
  // Compute variations of MET
  //************************************************************************************************
  double MET_trk_X     = 0. ;
  double MET_trk_Y     = 0. ;
  double MET_trk_sumPt = 0. ;
  double MET_trkplusneu_X = 0. ;
  double MET_trkplusneu_Y = 0. ;
  double MET_trkplusneu_sumPt = 0. ;
  double MET_neu_noFwd_X = 0. ;
  double MET_neu_noFwd_Y = 0. ;
  double MET_neu_noFwd_sumPt = 0. ;

  for (UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    // charged

    if ((fPFCandidates->At(i)->HasTrackerTrk() && fabs(fPFCandidates->At(i)->TrackerTrk()->DzCorrected(*primaryVertex)) < 0.1) ||
        (fPFCandidates->At(i)->HasGsfTrk()     && fabs(fPFCandidates->At(i)->GsfTrk()->DzCorrected(*primaryVertex)    ) < 0.1)) {      
      MET_trk_X -= fPFCandidates->At(i)->Px();
      MET_trk_Y -= fPFCandidates->At(i)->Py();
      MET_trk_sumPt += fPFCandidates->At(i)->Pt();
      MET_trkplusneu_X -= fPFCandidates->At(i)->Px();
      MET_trkplusneu_Y -= fPFCandidates->At(i)->Py();
      MET_trkplusneu_sumPt += fPFCandidates->At(i)->Pt();
    }

    //neutral
    if (fPFCandidates->At(i)->PFType()== PFCandidate::eNeutralHadron || fPFCandidates->At(i)->PFType()== PFCandidate::eGamma) {
      //Not forward
      if (fabs(fPFCandidates->At(i)->Eta()) < 3.0 && fPFCandidates->At(i)->Pt() > 2.0) {
        MET_neu_noFwd_X -= fPFCandidates->At(i)->Px();
        MET_neu_noFwd_Y -= fPFCandidates->At(i)->Py();
        MET_neu_noFwd_sumPt += fPFCandidates->At(i)->Pt();
      }
      
      //trk+neutral
      if (fabs(fPFCandidates->At(i)->Eta()) < 5.0 && fPFCandidates->At(i)->Pt() > 2.0 ) {
        MET_trkplusneu_X -= fPFCandidates->At(i)->Px();
        MET_trkplusneu_Y -= fPFCandidates->At(i)->Py();
        MET_trkplusneu_sumPt += fPFCandidates->At(i)->Pt();
      }
    }

  }

//   cout << "OLD Met calculation: " << endl;
//   MetTools metTools(CleanMuons, CleanElectrons, fPFCandidates, primaryVertex, 0.1, 8.0, 5.0);
//   MetTools metToolsNoFwd(CleanMuons, CleanElectrons, fPFCandidates, primaryVertex, 0.1, 2.0, 3.0);


  //
  // Fill event info tree
  //
  fEventInfo.eventweight  = 1.0;
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nTracks0     = fTracks->GetEntries();
  fEventInfo.nLeptons0    = fElectrons->GetEntries();
  fEventInfo.nCaloTowers0 = fCaloTowers->GetEntries();
  fEventInfo.nPV0         = fPrimVerts->GetEntries();
  fEventInfo.triggerBits  = trigbits;
  fEventInfo.l1triggerBits  = l1trigbits;
  fEventInfo.pvx          = primaryVertex->X();
  fEventInfo.pvy          = primaryVertex->Y();
  fEventInfo.pvz          = primaryVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.tcMEx        = fTCMet->At(0)->Mex();
  fEventInfo.tcMEy        = fTCMet->At(0)->Mey();
  fEventInfo.tcSumET      = fTCMet->At(0)->SumEt();
  fEventInfo.pfMEx        = fPFMet->At(0)->Mex();
  fEventInfo.pfMEy        = fPFMet->At(0)->Mey();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.pfTrackMEx   = MET_trk_X;
  fEventInfo.pfTrackMEy   = MET_trk_Y;
  fEventInfo.pfTrackSumET = MET_trk_sumPt;
  fEventInfo.pfNeutralMEx   = MET_trkplusneu_X;
  fEventInfo.pfNeutralMEy   = MET_trkplusneu_Y;
  fEventInfo.pfNeutralSumET = MET_trkplusneu_sumPt;
  fEventInfo.pfNeutralNoFwdMEx   = MET_neu_noFwd_X;
  fEventInfo.pfNeutralNoFwdMEy   = MET_neu_noFwd_Y;
  fEventInfo.pfNeutralNoFwdSumET = MET_neu_noFwd_sumPt;
  if (fReadPileupInfo) {
    if (fUseGen) {
      Int_t NPU = 0;
      Int_t NPU_PlusOne = 0;
      Int_t NPU_MinusOne = 0;
      for (UInt_t k=0; k < fPileupInfo->GetEntries() ; ++k) {
        if (fPileupInfo->At(k)->GetBunchCrossing() == 0)  NPU          = fPileupInfo->At(k)->GetPU_NumInteractions();
        if (fPileupInfo->At(k)->GetBunchCrossing() == 1)  NPU_PlusOne  = fPileupInfo->At(k)->GetPU_NumInteractions();
        if (fPileupInfo->At(k)->GetBunchCrossing() == -1) NPU_MinusOne = fPileupInfo->At(k)->GetPU_NumInteractions();
      }
      fEventInfo.nPUEvents   = NPU;
      fEventInfo.nPUMinusOne = NPU_MinusOne;
      fEventInfo.nPUPlusOne  = NPU_PlusOne;
    }
    fEventInfo.PileupEnergyDensity = fPileupEnergyDensity->At(0)->Rho();
    fEventInfo.PileupEnergyDensityHighEta = fPileupEnergyDensity->At(0)->RhoHighEta();
  } else {
    fEventInfo.nPUEvents = 0;
    fEventInfo.PileupEnergyDensity = 0;
    fEventInfo.PileupEnergyDensityHighEta =0;
  }

  if ( FoundZGamma || FoundWGammaEvent) {
    fEventInfo.VGammaEvent = kTRUE;
  } else {
    fEventInfo.VGammaEvent = kFALSE;
  }
  fEventInfo.hasGoodPV    = hasGoodPV;
  fEventInfo.ZeroJettinessPFCandidates = jettinessPFCands;
  fEventInfo.ZeroJettinessTracks = jettinessTracks;
  fEventInfo.ZeroJettinessCalotowers = -999;
  fEventInfo.BeamThrustPFCandidates = beamThrustPFCands;
  fEventInfo.BeamThrustTracks = beamThrustTracks;
  fEventInfo.BeamThrustCalotowers = -999;


  fEventTree->Fill();


  //********************************************************************************************
  //Debug
  //********************************************************************************************

  if ( fPrintDebug   
//        ||
//       (GetEventHeader()->RunNum() == 1 && GetEventHeader()->LumiSec() == 127 && GetEventHeader()->EvtNum() == 90831)
//        ||
//     (GetEventHeader()->RunNum() == 148031 && GetEventHeader()->EvtNum() == 507254451)
//              ||
//     (GetEventHeader()->RunNum() == 148822 && GetEventHeader()->EvtNum() == 215107719)
//              ||
//     (GetEventHeader()->RunNum() == 148862 && GetEventHeader()->EvtNum() == 655863681)
//              ||
//     (GetEventHeader()->RunNum() == 149181 && GetEventHeader()->EvtNum() == 469493528)
//       (GetEventHeader()->RunNum() == 142933 && GetEventHeader()->LumiSec() == 926 && GetEventHeader()->EvtNum() == 565024903)
//       ||
//       (GetEventHeader()->RunNum() == 142933 && GetEventHeader()->LumiSec() == 939 && GetEventHeader()->EvtNum() == 572606305) 




    ) {

    cout << "DEBUG!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
    cout << "Event " << GetEventHeader()->RunNum() << " " << GetEventHeader()->LumiSec() << " " << GetEventHeader()->EvtNum() << "\n";
    
    cout << "Met : " << fPFMet->At(0)->Pt() << " " << fPFMet->At(0)->Mex() << " " << fPFMet->At(0)->Mey() << " : " << endl;
    cout << "TrackMet : " << TMath::Sqrt(pow(MET_trk_X,2) + pow(MET_trk_Y,2))  << " " << MET_trk_X << " " << MET_trk_Y << endl;
   

    for (UInt_t i=0; i<fElectrons->GetEntries(); ++i) {

//       const Electron *ele = fElectrons->At(i);
      cout << "All Electron" << i << " : " << fElectrons->At(i)->Pt() << " " << fElectrons->At(i)->Eta() << " " 
           << fElectrons->At(i)->Phi() << " "  
           << fElectrons->At(i)->IsEB() << " "  
           << fElectrons->At(i)->Charge() << " " 
           << " : " 
           << fElectrons->At(i)->CoviEtaiEta() << " " 
           << fElectrons->At(i)->HadronicOverEm() << " " 
           << fElectrons->At(i)->EcalRecHitIsoDr03() << " " 
           << fElectrons->At(i)->TrackIsolationDr03() << " " 
           << fElectrons->At(i)->HcalTowerSumEtDr03() << " " 
           << IsolationTools::PFElectronIsolation(fElectrons->At(i), fPFCandidates, primaryVertex, 0.1, 1.0, 0.4, 0.0) << " " 
           << fElectrons->At(i)->BestTrk()->NExpectedHitsInner() << " "
           << fElectrons->At(i)->SCluster()->Seed()->EMax() << " " 
           << fElectrons->At(i)->SCluster()->Seed()->E3x3() << " "
           << fElectrons->At(i)->FBrem() << " "
           << fElectrons->At(i)->ESuperClusterOverP() << " "
           << fElectrons->At(i)->SCluster()->Eta() << " : "
           << fElectrons->At(i)->BestTrk()->D0Corrected(*primaryVertex) << " "
           << fElectrons->At(i)->BestTrk()->DzCorrected(*primaryVertex) << " "
           << " : "
           << ElectronTools::PassCustomID(fElectrons->At(i), ElectronTools::kVBTFWorkingPointLowPtId) << " " 
           << endl;  
      
      
      for (UInt_t p=0; p<fPFCandidates->GetEntries(); ++p) {
        if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
          if ( (fElectrons->At(i)->TrackerTrk() == fPFCandidates->At(p)->TrackerTrk()) or
               (fElectrons->At(i)->HasGsfTrk() and fElectrons->At(i)->GsfTrk() == fPFCandidates->At(p)->GsfTrk()) ) {
            cout << "Matched PFCandidate: " << fPFCandidates->At(p)->Pt() << " " << fPFCandidates->At(p)->Eta() << " " 
                 << fPFCandidates->At(p)->Phi() << " : " << fPFCandidates->At(i)->HasTrackerTrk() << " " << fPFCandidates->At(i)->HasGsfTrk() 
                 << endl;
          }
        }
      }
      
    }

    for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {

//       const Muon *mu = fMuons->At(i);
      cout << "Muon" << i << " : " << fMuons->At(i)->Pt() << " " << fMuons->At(i)->Eta() << " " << fMuons->At(i)->Phi() << " "
	   << fMuons->At(i)->Charge() << " " 
	   << fMuons->At(i)->BestTrk()->NHits() << " " 
	   << fMuons->At(i)->NMatches() << " " 
	   << fMuons->At(i)->NSegments() << " " 
	   << fMuons->At(i)->BestTrk()->NPixelHits() << " " 
	   << fMuons->At(i)->Quality().Quality(MuonQuality::GlobalMuonPromptTight) << " " 
	   << fMuons->At(i)->BestTrk()->PtErr() << " " ;
      if (fMuons->At(i)->GlobalTrk()) {
        cout << fMuons->At(i)->GlobalTrk()->Chi2()/fMuons->At(i)->GlobalTrk()->Ndof() << " ";
      } 
      cout   << fMuons->At(i)->IsoR03SumPt() << " " 
             << fMuons->At(i)->IsoR03EmEt() << " " 
             << fMuons->At(i)->IsoR03HadEt() << " " 	  
             << fMuons->At(i)->BestTrk()->D0Corrected(*primaryVertex) << " "
             << fMuons->At(i)->BestTrk()->DzCorrected(*primaryVertex) << " "
             << endl;   

      for (UInt_t p=0; p<fPFCandidates->GetEntries(); ++p) {
        if (fPFCandidates->At(i)->HasTrackerTrk() || fPFCandidates->At(i)->HasGsfTrk()) {
          if ( (fMuons->At(i)->TrackerTrk() == fPFCandidates->At(p)->TrackerTrk()) ) {
            cout << "Matched PFCandidate: " << fPFCandidates->At(p)->Pt() << " " << fPFCandidates->At(p)->Eta() << " " 
                 << fPFCandidates->At(p)->Phi() << " : " << fPFCandidates->At(i)->HasTrackerTrk() << " " << fPFCandidates->At(i)->HasGsfTrk() 
                 << endl;
          }
        }
      }

     
    }
  }


}

//--------------------------------------------------------------------------------------------------
  void HwwNtuplerMod::FillGenInfo(const MCParticleCol *GenLeptons, const MCParticleCol *GenNeutrinos, 
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
//     fGenInfo.ptBosonSystem = BosonSystem->Pt();
//     fGenInfo.mass = BosonSystem->Mass();
//     fGenInfo.pt = BosonSystem->Pt();
//     fGenInfo.y = BosonSystem->Rapidity();
//     fGenInfo.phi = BosonSystem->Phi();
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
//              cout << "GenJet " << i << " : " << fGenJets->At(i)->Pt() << " " << fGenJets->At(i)->Eta() << " " << fGenJets->At(i)->Phi() << endl;

        if (genJetsIndex == 0) {
          fGenInfo.jetpt_1 = fGenJets->At(i)->Pt();
 //          for (UInt_t k=0; k < fParticles->GetEntries() ; ++k) {
//             if (MathUtils::DeltaR(*fGenJets->At(i), *fParticles->At(k)) < 0.5) {
//               cout << "GenJet Particle: " << k << " : " << fParticles->At(k)->PdgId() << " " << fParticles->At(k)->Status() << " " << fParticles->At(k)->Pt() << " " << fParticles->At(k)->Eta() << " " << fParticles->At(k)->Phi() << endl;
//             }
//           }
        }
        if (genJetsIndex == 1) fGenInfo.jetpt_2 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 2) fGenInfo.jetpt_3 = fGenJets->At(i)->Pt();
        if (genJetsIndex == 3) fGenInfo.jetpt_4 = fGenJets->At(i)->Pt();
        if (fGenJets->At(i)->Pt() > 25.0) nGenJets++;
        genJetsIndex++;
      }
    }
    fGenInfo.nGenJets = nGenJets;

   

  }

  
}


//--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillElectron(const Electron *ele, Int_t isMCMatched, const Vertex* primaryVertex
  )
{
  assert(ele);
 

  //********************************************************************************
  //PF Isolation
  //********************************************************************************
  Double_t Electron_ChargedIso03 = 0;
  Double_t Electron_ChargedIso03FromOtherVertices = 0;
  Double_t Electron_NeutralHadronIso03_01Threshold = 0;
  Double_t Electron_GammaIso03_01Threshold = 0;
  Double_t Electron_NeutralHadronIso03_05Threshold = 0;
  Double_t Electron_GammaIso03_05Threshold = 0;
  Double_t Electron_NeutralHadronIso03_10Threshold = 0;
  Double_t Electron_GammaIso03_10Threshold = 0;
  Double_t Electron_NeutralHadronIso03_15Threshold = 0;
  Double_t Electron_GammaIso03_15Threshold = 0;

  Double_t Electron_ChargedIso04 = 0;
  Double_t Electron_ChargedIso04FromOtherVertices = 0;
  Double_t Electron_NeutralHadronIso04_01Threshold = 0;
  Double_t Electron_GammaIso04_01Threshold = 0;
  Double_t Electron_NeutralHadronIso04_05Threshold = 0;
  Double_t Electron_GammaIso04_05Threshold = 0;
  Double_t Electron_NeutralHadronIso04_10Threshold = 0;
  Double_t Electron_GammaIso04_10Threshold = 0;
  Double_t Electron_NeutralHadronIso04_15Threshold = 0;
  Double_t Electron_GammaIso04_15Threshold = 0;

  Double_t Electron_ChargedEMIsoVetoEtaStrip03 = 0;
  Double_t Electron_ChargedEMIsoVetoEtaStrip04 = 0;
  Double_t Electron_NeutralHadronIso007_01Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip03_01Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip04_01Threshold = 0;
  Double_t Electron_NeutralHadronIso007_05Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip03_05Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip04_05Threshold = 0;
  Double_t Electron_NeutralHadronIso007_10Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip03_10Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip04_10Threshold = 0;
  Double_t Electron_NeutralHadronIso007_15Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip03_15Threshold = 0;
  Double_t Electron_GammaIsoVetoEtaStrip04_15Threshold = 0;

  const PFCandidate *pfEleMatch = 0;

  Double_t zElectron = 0.0;
  if(ele->BestTrk()) zElectron = ele->BestTrk()->DzCorrected(*primaryVertex);
    
  for (UInt_t p=0; p<fPFCandidates->GetEntries();p++) {   
    const PFCandidate *pf = fPFCandidates->At(p);
      
    //find matching PFElectron
    if (!pfEleMatch) {
      if(
        (pf->GsfTrk() && ele->GsfTrk() &&
         pf->GsfTrk() == ele->GsfTrk())
        ||
        (pf->TrackerTrk() && ele->TrackerTrk() &&
         pf->TrackerTrk() == ele->TrackerTrk())
        ) {
        pfEleMatch = pf;
      }
    }
  
    if (fPrintDebug) {
      if (pf->HasTrackerTrk() || pf->HasGsfTrk()) {
        if (pf->TrackerTrk() == ele->TrackerTrk()) {
          cout << "PFMatch TrkTrack : " << pf->Pt() << " " << pf->Eta() << " " << pf->Phi() << " : " << pf->HasTrackerTrk() << " " << ele->HasTrackerTrk() << endl;
        }
        if(
          (pf->GsfTrk() && ele->GsfTrk() &&
           pf->GsfTrk() == ele->GsfTrk())
          ) {
        cout << "PFMatch GsfTrack : " << pf->Pt() << " " << pf->Eta() << " " << pf->Phi() << endl;
        }              
      }      
    }

  
    //exclude the electron itself
    if(pf->GsfTrk() && ele->GsfTrk() &&
       pf->GsfTrk() == ele->GsfTrk()) continue;
    if(pf->TrackerTrk() && ele->TrackerTrk() &&
       pf->TrackerTrk() == ele->TrackerTrk()) continue;      

    

    Double_t dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
      

    if (fPrintDebug) {
      cout << "PFCandidates For Electron " << ele->Pt() << " " << ele->Eta() << " " << ele->Phi() << "  : " << pf->Pt() << " " << pf->Eta() << " " << pf->Phi() << " : " << pf->PFType() << " " << dr << " ";

      if (dr < 0.4) {
        cout << "INCONE " ;
      }

    }




    if(pf->BestTrk()) {

      //remove charged particles from other vertices
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*primaryVertex) - zElectron);

      if (fPrintDebug) cout << " : " << deltaZ << " "; 

      if (deltaZ <= 0.1) {
        if (dr < 0.4) {
          Electron_ChargedIso04 += pf->Pt();          
        }
        if (dr < 0.3) {
          Electron_ChargedIso03 += pf->Pt();
        }
        if (pf->PFType() == PFCandidate::eElectron) {
          if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.3) {
            Electron_ChargedEMIsoVetoEtaStrip03 += pf->Pt();
          }
          if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.4) {
            Electron_ChargedEMIsoVetoEtaStrip04 += pf->Pt();
          }
        } 
      } else {
        if (dr < 0.4) {
          Electron_ChargedIso04FromOtherVertices += pf->Pt();          
        }
        if (dr < 0.3) {
          Electron_ChargedIso03FromOtherVertices += pf->Pt();
        }
      }
    } else if (pf->PFType() == PFCandidate::eGamma) {
      if (dr < 0.4) {
        if (pf->Pt() > 0.1) Electron_GammaIso04_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_GammaIso04_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_GammaIso04_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_GammaIso04_15Threshold += pf->Pt();
      }
      if (dr < 0.3) {
        if (pf->Pt() > 0.1) Electron_GammaIso03_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_GammaIso03_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_GammaIso03_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_GammaIso03_15Threshold += pf->Pt();
      }
      if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.3) {
        if (pf->Pt() > 0.1) Electron_GammaIsoVetoEtaStrip03_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_GammaIsoVetoEtaStrip03_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_GammaIsoVetoEtaStrip03_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_GammaIsoVetoEtaStrip03_15Threshold += pf->Pt();
      }
      if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.4) {
        if (pf->Pt() > 0.1) Electron_GammaIsoVetoEtaStrip04_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_GammaIsoVetoEtaStrip04_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_GammaIsoVetoEtaStrip04_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_GammaIsoVetoEtaStrip04_15Threshold += pf->Pt();
      }
    } else {
      if (dr < 0.4) {
        if (pf->Pt() > 0.1) Electron_NeutralHadronIso04_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_NeutralHadronIso04_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_NeutralHadronIso04_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_NeutralHadronIso04_15Threshold += pf->Pt();
      }
      if (dr < 0.3) {
        if (pf->Pt() > 0.1) Electron_NeutralHadronIso03_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_NeutralHadronIso03_05Threshold += pf->Pt();        
        if (pf->Pt() > 1.0) Electron_NeutralHadronIso03_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_NeutralHadronIso03_15Threshold += pf->Pt();
      }
      if (dr < 0.07) {
        if (pf->Pt() > 0.1) Electron_NeutralHadronIso007_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Electron_NeutralHadronIso007_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Electron_NeutralHadronIso007_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Electron_NeutralHadronIso007_15Threshold += pf->Pt();
      }
    }
    if (fPrintDebug) cout << endl;
  }


  TClonesArray &rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const Int_t index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron *pElectron = (TElectron*)rElectronArr[index];


  pElectron->pt              = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  if (pfEleMatch) {
    pElectron->pfPt              = pfEleMatch->Pt();
    pElectron->pfEta             = pfEleMatch->Eta();
    pElectron->pfPhi             = pfEleMatch->Phi();
  } else {
    pElectron->pfPt              = -999;
    pElectron->pfEta             = 0;
    pElectron->pfPhi             = 0;
  }
  pElectron->isEB            = ele->IsEB();
  pElectron->p               = ele->BestTrk()->P();
  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->trkIso04        = ele->TrackIsolationDr04();
  pElectron->emIso04         = ele->EcalRecHitIsoDr04();
  pElectron->hadIso04        = ele->HcalTowerSumEtDr04();

  pElectron->ChargedIso03 = Electron_ChargedIso03;
  pElectron->ChargedIso03FromOtherVertices = Electron_ChargedIso03FromOtherVertices;
  pElectron->NeutralHadronIso03_01Threshold = Electron_NeutralHadronIso03_01Threshold;
  pElectron->GammaIso03_01Threshold = Electron_GammaIso03_01Threshold;
  pElectron->NeutralHadronIso03_05Threshold = Electron_NeutralHadronIso03_05Threshold;
  pElectron->GammaIso03_05Threshold = Electron_GammaIso03_05Threshold;
  pElectron->NeutralHadronIso03_10Threshold = Electron_NeutralHadronIso03_10Threshold;
  pElectron->GammaIso03_10Threshold = Electron_GammaIso03_10Threshold;
  pElectron->NeutralHadronIso03_15Threshold = Electron_NeutralHadronIso03_15Threshold;
  pElectron->GammaIso03_15Threshold = Electron_GammaIso03_15Threshold;

  pElectron->ChargedIso04 = Electron_ChargedIso04;
  pElectron->ChargedIso04FromOtherVertices = Electron_ChargedIso04FromOtherVertices;
  pElectron->NeutralHadronIso04_01Threshold = Electron_NeutralHadronIso04_01Threshold;
  pElectron->GammaIso04_01Threshold = Electron_GammaIso04_01Threshold;
  pElectron->NeutralHadronIso04_05Threshold = Electron_NeutralHadronIso04_05Threshold;
  pElectron->GammaIso04_05Threshold = Electron_GammaIso04_05Threshold;
  pElectron->NeutralHadronIso04_10Threshold = Electron_NeutralHadronIso04_10Threshold;
  pElectron->GammaIso04_10Threshold = Electron_GammaIso04_10Threshold;
  pElectron->NeutralHadronIso04_15Threshold = Electron_NeutralHadronIso04_15Threshold;
  pElectron->GammaIso04_15Threshold = Electron_GammaIso04_15Threshold;

  pElectron->ChargedEMIsoVetoEtaStrip03 = Electron_ChargedEMIsoVetoEtaStrip03;
  pElectron->ChargedEMIsoVetoEtaStrip04 = Electron_ChargedEMIsoVetoEtaStrip04;
  pElectron->NeutralHadronIso007_01Threshold = Electron_NeutralHadronIso007_01Threshold ;
  pElectron->GammaIsoVetoEtaStrip03_01Threshold = Electron_GammaIsoVetoEtaStrip03_01Threshold;
  pElectron->GammaIsoVetoEtaStrip04_01Threshold = Electron_GammaIsoVetoEtaStrip04_01Threshold ;
  pElectron->NeutralHadronIso007_05Threshold = Electron_NeutralHadronIso007_05Threshold ;
  pElectron->GammaIsoVetoEtaStrip03_05Threshold = Electron_GammaIsoVetoEtaStrip03_05Threshold;
  pElectron->GammaIsoVetoEtaStrip04_05Threshold = Electron_GammaIsoVetoEtaStrip04_05Threshold ;
  pElectron->NeutralHadronIso007_10Threshold = Electron_NeutralHadronIso007_10Threshold ;
  pElectron->GammaIsoVetoEtaStrip03_10Threshold = Electron_GammaIsoVetoEtaStrip03_10Threshold;
  pElectron->GammaIsoVetoEtaStrip04_10Threshold = Electron_GammaIsoVetoEtaStrip04_10Threshold ;
  pElectron->NeutralHadronIso007_15Threshold = Electron_NeutralHadronIso007_15Threshold ;
  pElectron->GammaIsoVetoEtaStrip03_15Threshold = Electron_GammaIsoVetoEtaStrip03_15Threshold;
  pElectron->GammaIsoVetoEtaStrip04_15Threshold = Electron_GammaIsoVetoEtaStrip04_15Threshold ;



  pElectron->d0              = ele->BestTrk()->D0Corrected(*primaryVertex);
  pElectron->d0Err           = ele->BestTrk()->D0Err();
  pElectron->dz              = ele->BestTrk()->DzCorrected(*primaryVertex);
  pElectron->scEt            = ele->SCluster()->Et();
  pElectron->scEta           = ele->SCluster()->Eta();
  pElectron->scPhi           = ele->SCluster()->Phi();
  pElectron->fBrem           = ele->FBrem();
  pElectron->nBrem           = ele->NumberOfClusters() - 1;
  pElectron->EOverP          = ele->ESuperClusterOverP();
  pElectron->ESeedClusterOverPIn  = ele->ESeedClusterOverPIn();
  pElectron->ESeedClusterOverPout = ele->ESeedClusterOverPout();
  pElectron->HoverE          = ele->HadronicOverEm();
  pElectron->deltaEtaIn      = ele->DeltaEtaSuperClusterTrackAtVtx();
  pElectron->deltaPhiIn      = ele->DeltaPhiSuperClusterTrackAtVtx();
  pElectron->sigiEtaiEta     = ele->CoviEtaiEta();
  pElectron->sigiPhiiPhi     = ele->SCluster()->Seed()->CoviPhiiPhi();
  pElectron->nExpHitsInner   = ele->CorrectedNExpectedHitsInner();
  pElectron->partnerDeltaCot = ele->ConvPartnerDCotTheta();
  pElectron->partnerDist     = ele->ConvPartnerDist();
  pElectron->partnerRadius   = ele->ConvPartnerRadius();  
  pElectron->passSuperTightId      = ElectronTools::PassTightId(ele, fPrimVerts, fConversions, 2);
  pElectron->passHyperTight1Id     = ElectronTools::PassTightId(ele, fPrimVerts, fConversions, 3);
  pElectron->passHyperTight2Id     = ElectronTools::PassTightId(ele, fPrimVerts, fConversions, 4);
  pElectron->passHyperTight3Id     = ElectronTools::PassTightId(ele, fPrimVerts, fConversions, 5);
  pElectron->passHyperTight4Id     = ElectronTools::PassTightId(ele, fPrimVerts, fConversions, 6);
  pElectron->passCustomTightId  = ElectronTools::PassCustomID(ele, ElectronTools::kCustomIdTight);
  pElectron->likelihood      = ele->IDLikelihood();
  pElectron->mva             = ele->Mva();
  pElectron->q               = ele->Charge();

  //Conversion Veto
  Int_t passConvVeto = 0;
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kFALSE, kTRUE)) {
    passConvVeto += 1;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kFALSE, kFALSE)) {
    passConvVeto += 2;
  }
  //attempts matching to all conversions without regard to double-counting)
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, 2.0, kTRUE, kFALSE)) {
  //also attempt matching through ckf track in addition to gsf track
    passConvVeto += 4;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, -1.0, kTRUE, kFALSE)) {
  //relax lxy cut
    passConvVeto += 8;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-6, -999.9, kTRUE, kFALSE)) {
  //remove lxy cut
    passConvVeto += 16;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 1, 1e-10, -999.9, kTRUE, kFALSE)) {
  //relax probability cut
    passConvVeto += 32;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 2, 1e-10, -999.9, kTRUE, kFALSE)) {
  //relax hits before vtx cut
    passConvVeto += 64;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 999, 1e-10, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 128;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kFALSE, kTRUE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 256;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kFALSE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 512;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, 2.0, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 1024;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, -1.0, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 2048;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-6, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 4096;
  }
  if (ElectronTools::PassConversionFilter(ele, fConversions, fBeamSpot->At(0), 0, 1e-10, -999.9, kTRUE, kFALSE)) {
  //remove hits before vtx cut (this will probably be too tight)
    passConvVeto += 8192;
  }


  pElectron->isConv          = passConvVeto;

  pElectron->isEcalDriven    = ele->IsEcalDriven();
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi(), GetEventHeader()->RunNum(), GetEventHeader()->EvtNum());  
  pElectron->l1TriggerMatchBits = MatchL1(ele->SCluster()->Eta(),ele->SCluster()->Phi());

  pElectron->isMCReal          = isMCMatched;
  pElectron->ip3d              = ele->Ip3dPV();
  pElectron->ip3dSig           = ele->Ip3dPVSignificance();
  pElectron->isTrackerDriven   = ele->IsTrackerDriven();

  pElectron->HcalDepth1OverEcal = ele->HcalDepth1OverEcal();
  pElectron->HcalDepth2OverEcal = ele->HcalDepth2OverEcal();
  pElectron->dEtaCalo = ele->DeltaEtaSeedClusterTrackAtCalo();
  pElectron->dPhiCalo = ele->DeltaPhiSeedClusterTrackAtCalo();
  pElectron->PreShowerOverRaw = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();
  pElectron->CovIEtaIPhi = ele->SCluster()->Seed()->CoviEtaiPhi();
  pElectron->SCEtaWidth = ele->SCluster()->EtaWidth();
  pElectron->SCPhiWidth = ele->SCluster()->PhiWidth();
  pElectron->GsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
  
  pElectron->R9 = ele->SCluster()->R9();
  pElectron->SeedEMaxOverE = ele->SCluster()->Seed()->EMax() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedETopOverE = ele->SCluster()->Seed()->ETop() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedEBottomOverE = ele->SCluster()->Seed()->EBottom() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedELeftOverE = ele->SCluster()->Seed()->ELeft() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedERightOverE = ele->SCluster()->Seed()->ERight() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2ndOverE = ele->SCluster()->Seed()->E2nd() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x5RightOverE = ele->SCluster()->Seed()->E2x5Right() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x5LeftOverE = ele->SCluster()->Seed()->E2x5Left() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x5TopOverE = ele->SCluster()->Seed()->E2x5Top() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x5BottomOverE = ele->SCluster()->Seed()->E2x5Bottom() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x5MaxOverE = ele->SCluster()->Seed()->E2x5Max() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE1x3OverE = ele->SCluster()->Seed()->E1x3() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE3x1OverE = ele->SCluster()->Seed()->E3x1() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE1x5OverE = ele->SCluster()->Seed()->E1x5() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE2x2OverE = ele->SCluster()->Seed()->E2x2() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE3x2OverE = ele->SCluster()->Seed()->E3x2() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE3x3OverE = ele->SCluster()->Seed()->E3x3() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE4x4OverE = ele->SCluster()->Seed()->E4x4() / ele->SCluster()->Seed()->Energy();;
  pElectron->SeedE5x5OverE = ele->SCluster()->Seed()->E5x5() / ele->SCluster()->Seed()->Energy();;

}
    
  //--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillMuon(const Muon *mu, Int_t isMCMatched, const Vertex* primaryVertex)
{
  // N.B. TMuon::isUsed is not to be filled here
  
  assert(mu);
 

  //********************************************************************************
  //PF Isolation
  //********************************************************************************
  Double_t Muon_ChargedIso03 = 0;
  Double_t Muon_ChargedIso03FromOtherVertices = 0;
  Double_t Muon_NeutralIso03_01Threshold = 0;
  Double_t Muon_NeutralIso03_05Threshold = 0;
  Double_t Muon_NeutralIso03_10Threshold = 0;
  Double_t Muon_NeutralIso03_15Threshold = 0;

  Double_t Muon_ChargedIso04 = 0;
  Double_t Muon_ChargedIso04FromOtherVertices = 0;
  Double_t Muon_NeutralIso04_01Threshold = 0;
  Double_t Muon_NeutralIso04_05Threshold = 0;
  Double_t Muon_NeutralIso04_10Threshold = 0;
  Double_t Muon_NeutralIso04_15Threshold = 0;

  Double_t zMuon = 0.0;
  if(mu->BestTrk()) zMuon = mu->BestTrk()->DzCorrected(*primaryVertex);

  const PFCandidate *pfMuonMatched = 0;

  for (UInt_t p=0; p<fPFCandidates->GetEntries();p++) {   
    const PFCandidate *pf = fPFCandidates->At(p);
      
    if (!pfMuonMatched) {
      if(pf->TrackerTrk() && mu->TrackerTrk() &&
         pf->TrackerTrk() == mu->TrackerTrk()) {
        pfMuonMatched = pf;
      }           
    }

    //exclude the muon itself
    if(pf->TrackerTrk() && mu->TrackerTrk() &&
       pf->TrackerTrk() == mu->TrackerTrk()) continue;      

    Double_t dr = MathUtils::DeltaR(mu->Mom(), pf->Mom());
      
    if(pf->BestTrk()) {

      //remove charged particles from other vertices
      Double_t deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*primaryVertex) - zMuon);
      if (deltaZ <= 0.1) {
        if (dr < 0.4) {
          Muon_ChargedIso04 += pf->Pt();          
        }
        if (dr < 0.3) {
          Muon_ChargedIso03 += pf->Pt();
        }
      } else {
        if (dr < 0.4) {
          Muon_ChargedIso04FromOtherVertices += pf->Pt();          
        }
        if (dr < 0.3) {
          Muon_ChargedIso03FromOtherVertices += pf->Pt();
        }
      }
    } else {
      if (dr < 0.4) {
        if (pf->Pt() > 0.1) Muon_NeutralIso04_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Muon_NeutralIso04_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Muon_NeutralIso04_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Muon_NeutralIso04_15Threshold += pf->Pt();
      }
      if (dr < 0.3) {
        if (pf->Pt() > 0.1) Muon_NeutralIso03_01Threshold += pf->Pt();
        if (pf->Pt() > 0.5) Muon_NeutralIso03_05Threshold += pf->Pt();
        if (pf->Pt() > 1.0) Muon_NeutralIso03_10Threshold += pf->Pt();
        if (pf->Pt() > 1.5) Muon_NeutralIso03_15Threshold += pf->Pt();
      }
    }
  }

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); } 
  assert(muTrk);                 
  
  pMuon->pt        = muTrk->Pt();
  pMuon->eta       = muTrk->Eta();
  pMuon->phi       = muTrk->Phi();
  if (pfMuonMatched) {
    pMuon->pfPt        = pfMuonMatched->Pt();
    pMuon->pfEta       = pfMuonMatched->Eta();
    pMuon->pfPhi       = pfMuonMatched->Phi();  
  } else {
    pMuon->pfPt        = -999;
    pMuon->pfEta       = 0;
    pMuon->pfPhi       = 0;
  }
  pMuon->pterr     = muTrk->PtErr();
  pMuon->trkIso03  = mu->IsoR03SumPt();
  pMuon->emIso03   = mu->IsoR03EmEt();
  pMuon->hadIso03  = mu->IsoR03HadEt();
  pMuon->hoIso03   = mu->IsoR03HoEt();
  pMuon->trkIso05  = mu->IsoR05SumPt();
  pMuon->emIso05   = mu->IsoR05EmEt();
  pMuon->hadIso05  = mu->IsoR05HadEt();
  pMuon->hoIso05   = mu->IsoR05HoEt();

  pMuon->ChargedIso03 = Muon_ChargedIso03;
  pMuon->ChargedIso03FromOtherVertices = Muon_ChargedIso03FromOtherVertices;
  pMuon->NeutralIso03_01Threshold = Muon_NeutralIso03_01Threshold ;
  pMuon->NeutralIso03_05Threshold = Muon_NeutralIso03_05Threshold ;
  pMuon->NeutralIso03_10Threshold = Muon_NeutralIso03_10Threshold ;
  pMuon->NeutralIso03_15Threshold = Muon_NeutralIso03_15Threshold ;

  pMuon->ChargedIso04 = Muon_ChargedIso04;
  pMuon->ChargedIso04FromOtherVertices = Muon_ChargedIso04FromOtherVertices;
  pMuon->NeutralIso04_01Threshold = Muon_NeutralIso04_01Threshold ;
  pMuon->NeutralIso04_05Threshold = Muon_NeutralIso04_05Threshold ;
  pMuon->NeutralIso04_10Threshold = Muon_NeutralIso04_10Threshold ;
  pMuon->NeutralIso04_15Threshold = Muon_NeutralIso04_15Threshold ;
   
  pMuon->d0        = muTrk->D0Corrected(*primaryVertex);
  pMuon->d0Err     = muTrk->D0Err();
  pMuon->dz        = muTrk->DzCorrected(*primaryVertex);
  pMuon->tkNchi2   = muTrk->RChi2();
  
  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  //
  // NOTE:
  // It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA if
  // chamber matches for the TK and hits for the SA share DetIDs
  //   (see hypernews thread: https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html)
  //	      
  pMuon->typeBits = 0;
  if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  
  pMuon->nTkHits      = muTrk->NHits();
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nSeg         = mu->NSegments();
  pMuon->nMatch       = mu->NMatches();
  pMuon->hltMatchBits = MatchHLTMuon(muTrk->Pt(),muTrk->Eta(),muTrk->Phi(),GetEventHeader()->RunNum(), GetEventHeader()->EvtNum() );
  
  pMuon->staPt  = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Phi() : -999;

  pMuon->isMCReal = isMCMatched;
  pMuon->ip3d              = mu->Ip3dPV();
  pMuon->ip3dSig           = mu->Ip3dPVSignificance();

  pMuon->SegmentCompatibility     = fMuonTools->GetSegmentCompatability(mu);
  pMuon->CaloCompatilibity        = fMuonTools->GetCaloCompatability(mu, kTRUE, kTRUE);

  pMuon->TrkKink           = mu->TrkKink();
  pMuon->GlobalKink        = mu->GlbKink();
  pMuon->HadEnergy         = mu->HadEnergy();
  pMuon->HadS9Energy       = mu->HadS9Energy();
  pMuon->HoEnergy          = mu->HoEnergy();
  pMuon->HoS9Energy        = mu->HoS9Energy();
  pMuon->EmEnergy          = mu->EmEnergy();
  pMuon->EmS9Energy        = mu->EmS9Energy();

  pMuon->Station0NSegments    = mu->GetNSegments(0);
  pMuon->Station0TrackDist    = mu->GetTrackDist(0);
  pMuon->Station0TrackDistErr = mu->GetTrackDistErr(0);
  pMuon->Station0dX           = mu->GetDX(0);
  pMuon->Station0dY           = mu->GetDY(0);
  pMuon->Station0PullX        = mu->GetPullX(0);
  pMuon->Station0PullY        = mu->GetPullY(0);

  pMuon->Station1NSegments    = mu->GetNSegments(1);
  pMuon->Station1TrackDist    = mu->GetTrackDist(1);
  pMuon->Station1TrackDistErr = mu->GetTrackDistErr(1);
  pMuon->Station1dX           = mu->GetDX(1);
  pMuon->Station1dY           = mu->GetDY(1);
  pMuon->Station1PullX        = mu->GetPullX(1);
  pMuon->Station1PullY        = mu->GetPullY(1);

  pMuon->Station2NSegments    = mu->GetNSegments(2);
  pMuon->Station2TrackDist    = mu->GetTrackDist(2);
  pMuon->Station2TrackDistErr = mu->GetTrackDistErr(2);
  pMuon->Station2dX           = mu->GetDX(2);
  pMuon->Station2dY           = mu->GetDY(2);
  pMuon->Station2PullX        = mu->GetPullX(2);
  pMuon->Station2PullY        = mu->GetPullY(2);

  pMuon->Station3NSegments    = mu->GetNSegments(3);
  pMuon->Station3TrackDist    = mu->GetTrackDist(3);
  pMuon->Station3TrackDistErr = mu->GetTrackDistErr(3);
  pMuon->Station3dX           = mu->GetDX(3);
  pMuon->Station3dY           = mu->GetDY(3);
  pMuon->Station3PullX        = mu->GetPullX(3);
  pMuon->Station3PullY        = mu->GetPullY(3);

  pMuon->Station4NSegments    = mu->GetNSegments(4);
  pMuon->Station4TrackDist    = mu->GetTrackDist(4);
  pMuon->Station4TrackDistErr = mu->GetTrackDistErr(4);
  pMuon->Station4dX           = mu->GetDX(4);
  pMuon->Station4dY           = mu->GetDY(4);
  pMuon->Station4PullX        = mu->GetPullX(4);
  pMuon->Station4PullY        = mu->GetPullY(4);

  pMuon->Station5NSegments    = mu->GetNSegments(5);
  pMuon->Station5TrackDist    = mu->GetTrackDist(5);
  pMuon->Station5TrackDistErr = mu->GetTrackDistErr(5);
  pMuon->Station5dX           = mu->GetDX(5);
  pMuon->Station5dY           = mu->GetDY(5);
  pMuon->Station5PullX        = mu->GetPullX(5);
  pMuon->Station5PullY        = mu->GetPullY(5);

  pMuon->Station6NSegments    = mu->GetNSegments(6);
  pMuon->Station6TrackDist    = mu->GetTrackDist(6);
  pMuon->Station6TrackDistErr = mu->GetTrackDistErr(6);
  pMuon->Station6dX           = mu->GetDX(6);
  pMuon->Station6dY           = mu->GetDY(6);
  pMuon->Station6PullX        = mu->GetPullX(6);
  pMuon->Station6PullY        = mu->GetPullY(6);

  pMuon->Station7NSegments    = mu->GetNSegments(7);
  pMuon->Station7TrackDist    = mu->GetTrackDist(7);
  pMuon->Station7TrackDistErr = mu->GetTrackDistErr(7);
  pMuon->Station7dX           = mu->GetDX(7);
  pMuon->Station7dY           = mu->GetDY(7);
  pMuon->Station7PullX        = mu->GetPullX(7);
  pMuon->Station7PullY        = mu->GetPullY(7);


}




void HwwNtuplerMod::FillMuon(const Track *mu, const Vertex* primaryVertex)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];                 
  
  pMuon->pt           = mu->Pt();
  pMuon->eta          = mu->Eta();
  pMuon->phi          = mu->Phi();
  pMuon->pfPt         = -999;
  pMuon->pfEta        = 0;
  pMuon->pfPhi        = 0;
  pMuon->pterr        = mu->PtErr();
  pMuon->trkIso03     = 0;
  pMuon->emIso03      = 0;
  pMuon->hadIso03     = 0;
  pMuon->hoIso03      = 0;
  pMuon->trkIso05     = 0;
  pMuon->emIso05      = 0;
  pMuon->hadIso05     = 0;
  pMuon->hoIso05      = 0;

  pMuon->ChargedIso03 = 0;
  pMuon->ChargedIso03FromOtherVertices = 0;
  pMuon->NeutralIso03_01Threshold = 0;
  pMuon->NeutralIso03_05Threshold = 0;
  pMuon->NeutralIso03_10Threshold = 0;
  pMuon->NeutralIso03_15Threshold = 0;

  pMuon->ChargedIso04 = 0;
  pMuon->ChargedIso04FromOtherVertices = 0;
  pMuon->NeutralIso04_01Threshold = 0;
  pMuon->NeutralIso04_05Threshold = 0;
  pMuon->NeutralIso04_10Threshold = 0;
  pMuon->NeutralIso04_15Threshold = 0;
 
  pMuon->d0           = mu->D0Corrected(*primaryVertex);
  pMuon->d0Err        = mu->D0Err();
  pMuon->dz           = mu->DzCorrected(*primaryVertex);
  pMuon->tkNchi2      = mu->RChi2();
  pMuon->muNchi2      = mu->RChi2();
  pMuon->q            = mu->Charge();
  pMuon->nValidHits   = 0;  
  pMuon->qualityBits  = 0;      
  pMuon->typeBits     = 0; 
  pMuon->nTkHits      = mu->NHits();
  pMuon->nPixHits     = mu->NPixelHits();
  pMuon->nSeg         = 0;
  pMuon->nMatch       = 0;
  pMuon->hltMatchBits = MatchHLT(mu->Pt(),mu->Eta(),mu->Phi());
  
  pMuon->staPt  = -1;  
  pMuon->staEta = -999;
  pMuon->staPhi = -999;
  

  pMuon->isMCReal = 0;
  pMuon->ip3d              = 0;
  pMuon->ip3dSig           = 0;

  pMuon->SegmentCompatibility     = 0;
  pMuon->CaloCompatilibity        = 0;

  pMuon->TrkKink           = 0;
  pMuon->GlobalKink        = 0;
  pMuon->HadEnergy         = 0;
  pMuon->HadS9Energy       = 0;
  pMuon->HoEnergy          = 0;
  pMuon->HoS9Energy        = 0;
  pMuon->EmEnergy          = 0;
  pMuon->EmS9Energy        = 0;

  pMuon->Station0NSegments    = 0;
  pMuon->Station0TrackDist    = 0;
  pMuon->Station0TrackDistErr = 0;
  pMuon->Station0dX           = 0;
  pMuon->Station0dY           = 0;
  pMuon->Station0PullX        = 0;
  pMuon->Station0PullY        = 0;

  pMuon->Station1NSegments    = 0;
  pMuon->Station1TrackDist    = 0;
  pMuon->Station1TrackDistErr = 0;
  pMuon->Station1dX           = 0;
  pMuon->Station1dY           = 0;
  pMuon->Station1PullX        = 0;
  pMuon->Station1PullY        = 0;

  pMuon->Station2NSegments    = 0;
  pMuon->Station2TrackDist    = 0;
  pMuon->Station2TrackDistErr = 0;
  pMuon->Station2dX           = 0;
  pMuon->Station2dY           = 0;
  pMuon->Station2PullX        = 0;
  pMuon->Station2PullY        = 0;

  pMuon->Station3NSegments    = 0;
  pMuon->Station3TrackDist    = 0;
  pMuon->Station3TrackDistErr = 0;
  pMuon->Station3dX           = 0;
  pMuon->Station3dY           = 0;
  pMuon->Station3PullX        = 0;
  pMuon->Station3PullY        = 0;

  pMuon->Station4NSegments    = 0;
  pMuon->Station4TrackDist    = 0;
  pMuon->Station4TrackDistErr = 0;
  pMuon->Station4dX           = 0;
  pMuon->Station4dY           = 0;
  pMuon->Station4PullX        = 0;
  pMuon->Station4PullY        = 0;

  pMuon->Station5NSegments    = 0;
  pMuon->Station5TrackDist    = 0;
  pMuon->Station5TrackDistErr = 0;
  pMuon->Station5dX           = 0;
  pMuon->Station5dY           = 0;
  pMuon->Station5PullX        = 0;
  pMuon->Station5PullY        = 0;

  pMuon->Station6NSegments    = 0;
  pMuon->Station6TrackDist    = 0;
  pMuon->Station6TrackDistErr = 0;
  pMuon->Station6dX           = 0;
  pMuon->Station6dY           = 0;
  pMuon->Station6PullX        = 0;
  pMuon->Station6PullY        = 0;

  pMuon->Station7NSegments    = 0;
  pMuon->Station7TrackDist    = 0;
  pMuon->Station7TrackDistErr = 0;
  pMuon->Station7dX           = 0;
  pMuon->Station7dY           = 0;
  pMuon->Station7PullX        = 0;
  pMuon->Station7PullY        = 0;


}




void HwwNtuplerMod::FillJet(const PFJet *jet, const Vertex* primaryVertex)
{
  TClonesArray &rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const Int_t index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) TJet();
  TJet *pPFJet = (TJet*)rPFJetArr[index]; 
   
  pPFJet->pt           = jet->Pt();
  pPFJet->eta          = jet->Eta();
  pPFJet->phi          = jet->Phi();
  pPFJet->mass         = jet->Mass();
  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());  
  pPFJet->TrackCountingHighEffBJetTagsDisc = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->TrackCountingHighPurBJetTagsDisc  = jet->TrackCountingHighPurBJetTagsDisc();
  pPFJet->GhostTrackBJetTagsDisc = jet->GhostTrackBJetTagsDisc();
  pPFJet->SoftElectronByPtBJetTagsDisc = jet->SoftElectronByPtBJetTagsDisc();
  pPFJet->SoftElectronByIP3dBJetTagsDisc = jet->SoftElectronByIP3dBJetTagsDisc();
  pPFJet->SoftMuonByPtBJetTagsDisc = jet->SoftMuonByPtBJetTagsDisc();
  pPFJet->SoftMuonByIP3dBJetTagsDisc = jet->SoftMuonByIP3dBJetTagsDisc();
  pPFJet->SoftMuonBJetTagsDisc = jet->SoftMuonBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexHighPurBJetTagsDisc = jet->SimpleSecondaryVertexHighPurBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexHighEffBJetTagsDisc = jet->SimpleSecondaryVertexHighEffBJetTagsDisc();
  pPFJet->SimpleSecondaryVertexBJetTagsDisc = jet->SimpleSecondaryVertexBJetTagsDisc();
  pPFJet->CombinedSecondaryVertexBJetTagsDisc = jet->CombinedSecondaryVertexBJetTagsDisc();
  pPFJet->CombinedSecondaryVertexMVABJetTagsDisc = jet->CombinedSecondaryVertexMVABJetTagsDisc();
  pPFJet->PassBetaVertexAssociationCut = JetTools::PassBetaVertexAssociationCut(jet, fPrimVerts->At(0), fPrimVerts,0.2);
  pPFJet->NConstituents = jet->NConstituents();
  pPFJet->NeutralHadronFraction = jet->NeutralHadronEnergy() / jet->E();
  pPFJet->NeutralEMFraction = jet->NeutralEmEnergy() / jet->E();
  pPFJet->ChargedHadronFraction = jet->ChargedHadronEnergy() / jet->E();
  pPFJet->ChargedEMFraction = jet->ChargedEmEnergy() / jet->E();
  pPFJet->ChargedMultiplicity = jet->ChargedMultiplicity();
  pPFJet->JetProbabilityBJetTagsDisc = jet->JetProbabilityBJetTagsDisc();
  pPFJet->JetBProbabilityBJetTagsDisc = jet->JetBProbabilityBJetTagsDisc();

  pPFJet->JetArea = jet->JetArea();

  Double_t ChargedPtSqrSum = 0;
  Double_t DzSum = 0;
  Double_t DzPtSqrWeightedSum = 0;
  Int_t NChargedPFCandidates = 0;
  for (UInt_t i = 0; i < jet->NPFCands(); ++i) {
    //Charged PF Cand
    if (jet->PFCand(i)->HasTrackerTrk() || jet->PFCand(i)->HasGsfTrk()) {
      Double_t dz = 0;
      ChargedPtSqrSum += pow(jet->PFCand(i)->Pt(),2);
      NChargedPFCandidates++;
      if (jet->PFCand(i)->HasTrackerTrk()) {
        dz = jet->PFCand(i)->TrackerTrk()->DzCorrected(*primaryVertex);
      } else {
        dz = jet->PFCand(i)->GsfTrk()->DzCorrected(*primaryVertex);
      }

      DzSum += dz;
      DzPtSqrWeightedSum += dz * pow(jet->PFCand(i)->Pt(),2);
    }
  }
  
  pPFJet->DzAvg = DzSum / NChargedPFCandidates;
  pPFJet->DzPtSqrWeightedAvg = DzPtSqrWeightedSum / ChargedPtSqrSum;

  pPFJet->rawPt = jet->RawMom().Pt();

}
   

  //--------------------------------------------------------------------------------------------------
void HwwNtuplerMod::FillPhoton(const Photon *pho)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pPhoton = (TPhoton*)rPhotonArr[index];
  
  pPhoton->et		  = pho->Et(); 
  pPhoton->eta  	  = pho->Eta();
  pPhoton->phi  	  = pho->Phi();
  pPhoton->scEt		  = pho->SCluster()->Et(); 
  pPhoton->scEta  	  = pho->SCluster()->Eta();
  pPhoton->scPhi  	  = pho->SCluster()->Phi();
  pPhoton->trkIso03Hollow = pho->HollowConeTrkIsoDr03(); 
  pPhoton->trkIso03Solid  = pho->SolidConeTrkIsoDr03();
  pPhoton->emIso03        = pho->EcalRecHitIsoDr03();
  pPhoton->hadIso03	  = pho->HcalTowerSumEtDr03(); 
  pPhoton->HoverE	  = pho->HadOverEm();
  pPhoton->R9		  = pho->R9();
  pPhoton->sigiEtaiEta    = pho->CoviEtaiEta();
  pPhoton->hltMatchBits   = MatchHLT(pho->SCluster()->Eta(),pho->SCluster()->Phi());
  pPhoton->scID           = pho->SCluster()->GetUniqueID();
  pPhoton->hasPixelSeed   = pho->HasPixelSeed();
}
 

//--------------------------------------------------------------------------------------------------
ULong_t HwwNtuplerMod::MatchL1(const Double_t eta, const Double_t phi)
{

  ULong_t bits = 0;
  
  const Double_t L1MatchDR = 0.5;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);

    for(ULong_t iseed=0; iseed<fL1SeedModuleNamesv.size(); iseed++) {

      for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {

        const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
        if(!trigname) continue;
   
        const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
        if(!list) continue;
        TIter iter(list->MakeIterator());
        const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
        while(to) {         
        
          
        
          if(to->IsL1()) {
            Bool_t match = true;
            
            //match L1 seed module name
            if (!(string(to->ModuleName()) == string(fL1SeedModuleNamesv[iseed].Data()) )) { match = false; }

//             if (match) cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " " << endl;
        


            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > L1MatchDR) { match = false; }

            // set appropriate bits
            if(match) { 
              bits |= fL1SeedModuleIdsv[iseed]; 
            }
          }
          to = dynamic_cast<const TriggerObject*>(iter.Next());
        }
      }
    }
  }
  
  return bits;
}


ULong_t HwwNtuplerMod::MatchHLT(const Double_t eta, const Double_t phi, Double_t runNum, Double_t evtNum)
{



//   if (fPrintDebug) {
//     cout << "Match To HLT Object : " << eta << " " << phi << endl;
//   }

  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.1;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      bool doDebug = kFALSE;

      if (doDebug) cout << "DEBUG: Match To HLT Object : " << eta << " " << phi << endl;

      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        
        if (doDebug || fPrintDebug) {         
          cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " : " 
               << MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) << " --> " << hltMatchR << " : " 
               << fFirstTriggerObjectIdsv[itrig] << " " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " "  
               << (fFirstTriggerObjectIdsv[itrig] != 0) << " " 
               << (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "") << " " 
               << endl;          
        }
        
        if(to->IsHLT()) {
          
          //Match First Trigger Object
          if (fFirstTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }

            // set appropriate bits
            if(match) { 
              bits |= fFirstTriggerObjectIdsv[itrig]; 
              if (doDebug ||fPrintDebug) {
                cout << "Matched First : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
              }
            }
          }
          
          //Match Second Trigger Object
          if (fSecondTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fSecondTriggerObjectModuleNamesv[itrig].Data()) || string(fSecondTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fSecondTriggerObjectIdsv[itrig];               
              if (doDebug ||fPrintDebug) {
                cout << "Matched Second : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
              }
            }
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }
    }
  }
  
  return bits;
}

ULong_t HwwNtuplerMod::MatchHLTMuon(const Double_t pt, const Double_t eta, const Double_t phi, Double_t runNum, Double_t evtNum)
{
  ULong_t bits = 0;
  
//   if (fPrintDebug) {
//     cout << "Match HLT Mu : " << pt << " " << eta << " " << phi << endl;
//   }


  const Double_t hltMatchR = 0.2;
//   const Double_t hltMatchPtFrac = 10;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(ULong_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {  

//         if (fPrintDebug) {
    
//           cout << to->Pt() << " " << to->Eta() << " " << to->Phi() << " : " << to->IsHLT() << " " << to->IsL1() << " " << to->TrigName() << " " << to->FilterName() << " " << to->ModuleName() << " " << to->TagName() << " " << endl;
          
//         }
        
        if(to->IsHLT()) {
          
          //Match First Trigger Object
          if (fFirstTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fFirstTriggerObjectModuleNamesv[itrig].Data()) || string(fFirstTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
//             // pT matching
//             if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fFirstTriggerObjectIdsv[itrig]; 
//               if (fPrintDebug) {
//                 cout << "Matched First : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
//               }              
            }
          }
          
          //Match Second Trigger Object
          if (fSecondTriggerObjectIdsv[itrig] != 0 && 
              (string(to->ModuleName()) == string(fSecondTriggerObjectModuleNamesv[itrig].Data()) || string(fSecondTriggerObjectModuleNamesv[itrig].Data()) == "")
            ) {
            Bool_t match = true;
            
            // eta-phi matching
            if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) { match = false; }
            
//             // pT matching
//             if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) { match = false; }
            
            // set appropriate bits
            if(match) { 
              bits |= fSecondTriggerObjectIdsv[itrig]; 
//               if (fPrintDebug) {
//                 cout << "Matched Second : " << fFirstTriggerObjectModuleNamesv[itrig].Data() << " : " << endl; 
//               }
            }
          }

        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}


// //--------------------------------------------------------------------------------------------------
// Bool_t HwwNtuplerMod::IsConversion(const Electron *ele) 
// {
//   const Track* trk = ele->GsfTrk();
//   Bool_t isConv = kFALSE;

//   for (UInt_t k=0; k<fConversions->GetEntries(); ++k) {
//     const DecayParticle *decaypart = fConversions->At(k);
//     if(decaypart->Prob() < 0.005) continue;
//     if(decaypart->Lxy()  < 0)     continue;
//     if(decaypart->Lz()   < 0)     continue;
//     for(UInt_t i=0;i<decaypart->NDaughters();i++){
//       const Track *t = dynamic_cast<const ChargedParticle*>(decaypart->Daughter(i))->Trk();
//       if(t==trk) isConv = kTRUE;
//     }
//     if(!isConv) continue;
//     for(UInt_t i=0;i<decaypart->NDaughters();i++){
//       const Track *t       = dynamic_cast<const ChargedParticle*>(decaypart->Daughter(i))->Trk();
//       const StableData *sd = dynamic_cast<const StableData*>(decaypart->DaughterDat(i));
//       if(sd->NWrongHits()!=0) isConv = kFALSE;
//       if (i==0) {
//         if(t->NHits()<8)    isConv = kFALSE;
//         if(t->Prob()<0.005) isConv = kFALSE;
//       }
//     }
//   }
//   return isConv;
// }

