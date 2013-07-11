#ifndef EWKANA_NTUPLER_TDIELECTRON_HH
#define EWKANA_NTUPLER_TDIELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TDielectron : public TObject
  {
    public:
      TDielectron(){}
      ~TDielectron(){} 
   
      Float_t mass, pt, y, phi;  // dielectron kinematics
  
      // leading electron
      Float_t pt_1, eta_1, phi_1, e_1;
      Float_t trkIso03_1;
      Float_t emIso03_1;
      Float_t hadIso03_1;
      Float_t d0_1, d0Err_1, dz_1;
      Float_t scEt_1, scEta_1, scPhi_1;
      Float_t HoverE_1;
      Float_t deltaEtaIn_1;
      Float_t deltaPhiIn_1;
      Float_t sigiEtaiEta_1;
      Float_t nExpHitsInner_1;
      Float_t partnerDeltaCot_1;     
      Float_t partnerDist_1;
      Float_t partnerRadius_1;      
      Int_t q_1;
      UInt_t hltMatchBits_1;
      Bool_t isConv_1;
      Bool_t isEcalDriven_1;

      // lagging electron
      Float_t pt_2, eta_2, phi_2, e_2;
      Float_t trkIso03_2;
      Float_t emIso03_2;
      Float_t hadIso03_2;
      Float_t d0_2, d0Err_2, dz_2;
      Float_t scEt_2, scEta_2, scPhi_2;
      Float_t HoverE_2;
      Float_t deltaEtaIn_2;
      Float_t deltaPhiIn_2;
      Float_t sigiEtaiEta_2;
      Float_t nExpHitsInner_2;
      Float_t partnerDeltaCot_2;     
      Float_t partnerDist_2;
      Float_t partnerRadius_2;      
      Int_t q_2;
      UInt_t hltMatchBits_2;
      Bool_t isConv_2;
      Bool_t isEcalDriven_2;
    
    ClassDef(TDielectron,1)
  };
}
#endif
