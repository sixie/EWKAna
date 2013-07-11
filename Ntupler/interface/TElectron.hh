#ifndef EWKANA_NTUPLER_TELECTRON_HH
#define EWKANA_NTUPLER_TELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TElectron : public TObject
  {
    public:
      TElectron(){}
      ~TElectron(){}
    
      Float_t pt, eta, phi, p;     // kinematics      
      Float_t pfPt, pfEta, pfPhi;  // kinematics      
      Bool_t  isEB;                // is in barrel
      Float_t trkIso03;            // track isolation
      Float_t emIso03;             // ECAL-based isolation
      Float_t hadIso03;            // HCAL-based isolation
      Float_t trkIso04;            // track isolation
      Float_t emIso04;             // ECAL-based isolation
      Float_t hadIso04;            // HCAL-based isolation

      Float_t ChargedIso03;
      Float_t ChargedIso03FromOtherVertices;

      Float_t NeutralHadronIso03_01Threshold;
      Float_t GammaIso03_01Threshold;
      Float_t NeutralHadronIso03_05Threshold;
      Float_t GammaIso03_05Threshold;
      Float_t NeutralHadronIso03_10Threshold;
      Float_t GammaIso03_10Threshold;
      Float_t NeutralHadronIso03_15Threshold;
      Float_t GammaIso03_15Threshold;

      Float_t ChargedIso04;
      Float_t ChargedIso04FromOtherVertices;
      Float_t NeutralHadronIso04_01Threshold;
      Float_t GammaIso04_01Threshold;
      Float_t NeutralHadronIso04_05Threshold;
      Float_t GammaIso04_05Threshold;
      Float_t NeutralHadronIso04_10Threshold;
      Float_t GammaIso04_10Threshold;
      Float_t NeutralHadronIso04_15Threshold;
      Float_t GammaIso04_15Threshold;

      Float_t ChargedEMIsoVetoEtaStrip03;      
      Float_t ChargedEMIsoVetoEtaStrip04;
      Float_t NeutralHadronIso007_01Threshold;
      Float_t GammaIsoVetoEtaStrip03_01Threshold;      
      Float_t GammaIsoVetoEtaStrip04_01Threshold;      
      Float_t NeutralHadronIso007_05Threshold;
      Float_t GammaIsoVetoEtaStrip03_05Threshold;      
      Float_t GammaIsoVetoEtaStrip04_05Threshold;      
      Float_t NeutralHadronIso007_10Threshold;
      Float_t GammaIsoVetoEtaStrip03_10Threshold;      
      Float_t GammaIsoVetoEtaStrip04_10Threshold;      
      Float_t NeutralHadronIso007_15Threshold;
      Float_t GammaIsoVetoEtaStrip03_15Threshold;      
      Float_t GammaIsoVetoEtaStrip04_15Threshold;      



      Float_t d0, d0Err, dz;       // impact parameter
      Float_t scEt, scEta, scPhi;  // supercluster
      Float_t fBrem, EOverP;       // fBrem, EOverP
      Int_t   nBrem;               // Number of Brems
      Float_t ESeedClusterOverPIn; // ESeedClusterOverPout     
      Float_t ESeedClusterOverPout;// ESeedClusterOverPout     
      Float_t HoverE;              // H / E
      Float_t deltaEtaIn;          // eta difference between track (at vertex) and SC
      Float_t deltaPhiIn;          // phi difference between track (at vertex) and SC
      Float_t sigiEtaiEta;         // eta-width of shower in number of crystals
      Float_t sigiPhiiPhi;         // phi-width of shower in number of crystals
      Float_t nExpHitsInner;       // number of hits expected before first hit
      Float_t partnerDeltaCot;     // cot(theta) difference with conversion partner track       
      Float_t partnerDist;         // distance in x-y plane to nearest conversion partner track
      Float_t partnerRadius;       // radius of helix intersection with conversion partner track
      Int_t   passSuperTightId;    // Bit for CiC ID
      Int_t   passHyperTight1Id;   // Bit for CiC ID
      Int_t   passHyperTight2Id;   // Bit for CiC ID
      Int_t   passHyperTight3Id;   // Bit for CiC ID
      Int_t   passHyperTight4Id;   // Bit for CiC ID
      Bool_t  passCustomTightId;   // Old Hww cuts
      Float_t likelihood;          // likelihood
      Float_t mva;                 // mva
      Int_t q;                     // charge
      ULong_t hltMatchBits;         // bits for matching with HLT primitives
      UInt_t l1TriggerMatchBits;   // bits for matching with L1 Seeds
      Int_t  isConv;               // is conversion? (vertexing method)
      Bool_t isEcalDriven;         // is ECAL seeded electron?
      Int_t  isMCReal;        
      Float_t ip3d, ip3dSig;       //IP3D PV
      Bool_t isTrackerDriven;      // is TrackerDriven electron


      Float_t HcalDepth1OverEcal;
      Float_t HcalDepth2OverEcal;
      Float_t dEtaCalo;
      Float_t dPhiCalo;
      Float_t PreShowerOverRaw;
      Float_t CovIEtaIPhi;
      Float_t SCEtaWidth;
      Float_t SCPhiWidth;
      Float_t  GsfTrackChi2OverNdof;
  
      Float_t R9;
      Float_t SeedEMaxOverE;
      Float_t SeedETopOverE;
      Float_t SeedEBottomOverE;
      Float_t SeedELeftOverE;
      Float_t SeedERightOverE;
      Float_t SeedE2ndOverE;
      Float_t SeedE2x5RightOverE;
      Float_t SeedE2x5LeftOverE;
      Float_t SeedE2x5TopOverE;
      Float_t SeedE2x5BottomOverE;
      Float_t SeedE2x5MaxOverE;
      Float_t SeedE1x3OverE;
      Float_t SeedE3x1OverE;
      Float_t SeedE1x5OverE;
      Float_t SeedE2x2OverE;
      Float_t SeedE3x2OverE;
      Float_t SeedE3x3OverE;
      Float_t SeedE4x4OverE;
      Float_t SeedE5x5OverE;

    ClassDef(TElectron,1)
  };
}
#endif
