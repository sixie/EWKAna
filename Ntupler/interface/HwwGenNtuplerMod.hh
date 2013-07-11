#ifndef EWKANA_NTUPLER_HWWGENNTUPLERMOD_H
#define EWKANA_NTUPLER_HWWGENNTUPLERMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"

#include "EWKAnaDefs.hh"
#include "TEventInfo.hh"
#include "TGenInfo.hh"
#include <TClonesArray.h>
#include "TDielectron.hh"
#include "TElectron.hh"
#include "TJet.hh"
#include "TMuon.hh"
#include "TPhoton.hh"

#include <vector>

class TTree;
class TFile;
class TString;

namespace mithep
{
  class HwwGenNtuplerMod : public BaseMod
  {    
    public:
      HwwGenNtuplerMod(const char *name="HwwGenNtuplerMod", const char *title="BAMBU to ntuple");
      ~HwwGenNtuplerMod();	

      void SetOutputName(const char *f)  { fOutputName = f; }
      
      void SetPDFName(const char *s)     { fPDFName            = s;    }
      void SetComputePDFWeights(Bool_t b) { fComputePDFWeights = b;    }

            
    protected:
      void Begin();
      void BeginRun();
      void EndRun();
      void SlaveBegin();
      void SlaveTerminate();
      void Terminate();
      void Process();

      // Fill generator info data object
      void FillGenInfo(const MCParticleCol *GenLeptons, const MCParticleCol *GenNeutrinos, 
                       const MCParticleCol *GenBosons);
      
      
      
      TFile                  *fOutputFile;      // output file handle
      TString                 fOutputName;      // output file name
      
      TString                 fPDFName;         //PDF name
      TString                 fPartName;        // MC particle collection name
      TString                 fMCEvtInfoName;   // MC event info name


      const MCParticleCol    *fParticles;       // MC particle collection handle
      const MCEventInfo      *fMCEvtInfo;       // MC event info handle
      const GenJetCol        *fGenJets;         // MC particle collection handle
      const PFMetCol         *fPFMet;           // particle flow MET handle

      Bool_t                  fPrintDebug;
      Bool_t                  fUseGen;          // flag whether to look at generator info
      Bool_t                  fComputePDFWeights; //compute weights for different pdfs?

      
      TTree*                  fEventTree;       // event tree
            
      TEventInfo              fEventInfo;       // general event information
      TGenInfo                fGenInfo;         // generator information
 
      Int_t                   fNPDFMembers;     // Count How many PDF Members we have
      Float_t                 fPDFWeights[100]; // Save weights for each of the PDFs

      Int_t                   fNGenTaus;        
      Float_t                 fGenTauMass[10]; 
      Float_t                 fGenTauPt[10]; 
      Float_t                 fGenTauEta[10];  
      Float_t                 fGenTauPhi[10];  
      Int_t                   fGenTauDecayLeptonType[10]; 
      Float_t                 fGenTauDecayLeptonPt[10]; 
      Float_t                 fGenTauDecayLeptonEta[10];  
      Float_t                 fGenTauDecayLeptonPhi[10];  
      Int_t                   fGenTauDecayTauNeutrinoType[10]; 
      Float_t                 fGenTauDecayTauNeutrinoPt[10]; 
      Float_t                 fGenTauDecayTauNeutrinoEta[10];  
      Float_t                 fGenTauDecayTauNeutrinoPhi[10];  
      Int_t                   fGenTauDecayLeptonNeutrinoType[10]; 
      Float_t                 fGenTauDecayLeptonNeutrinoPt[10]; 
      Float_t                 fGenTauDecayLeptonNeutrinoEta[10];  
      Float_t                 fGenTauDecayLeptonNeutrinoPhi[10];  
      Int_t                   fNLeptons;
      Int_t                   fLeptonType[10]; 
      Float_t                 fLeptonPt[10]; 
      Float_t                 fLeptonEta[10];  
      Float_t                 fLeptonPhi[10];  
      Float_t                 fNtuplePFMet;  
      Float_t                 fNtuplePFMetPhi;  


    ClassDef(HwwGenNtuplerMod,1)
  };
}

#endif    
