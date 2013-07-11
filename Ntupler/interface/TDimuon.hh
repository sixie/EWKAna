#ifndef EWKANA_NTUPLER_TDIMUON_HH
#define EWKANA_NTUPLER_TDIMUON_HH

#include <TObject.h>

namespace mithep 
{
  class TDimuon : public TObject
  {
    public:
      TDimuon(){}
      ~TDimuon(){}

      Float_t mass, pt, y, phi;  // dimuon kinematics
  
      // leading pT muon
      Float_t pt_1, eta_1, phi_1;
      Float_t staMass_1, staPt_1, staEta_1, staPhi_1;
      Float_t trkIso03_1;
      Float_t emIso03_1;
      Float_t hadIso03_1;
      Float_t hoIso03_1;
      Float_t d0_1, d0Err_1, dz_1;
      Float_t caloComp_1;
      Float_t segComp_1;
      Float_t tkNchi2_1;
      Float_t muNchi2_1;
      Int_t q_1;
      Int_t lastHit_1;
      Int_t nValidHits_1;
      UInt_t qualityBits_1;
      UInt_t typeBits_1;
      UInt_t nTkHits_1;
      UInt_t nPixHits_1;
      UInt_t nSeg_1;
      UInt_t hltMatchBits_1;
	 
      // lagging pT muon
      Float_t pt_2, eta_2, phi_2;
      Float_t staMass_2, staPt_2, staEta_2, staPhi_2;
      Float_t trkIso03_2;
      Float_t emIso03_2;
      Float_t hadIso03_2;
      Float_t hoIso03_2;
      Float_t d0_2, d0Err_2, dz_2;
      Float_t caloComp_2;
      Float_t segComp_2;
      Float_t tkNchi2_2;
      Float_t muNchi2_2;
      Int_t q_2;
      Int_t lastHit_2;
      Int_t nValidHits_2;
      UInt_t qualityBits_2;
      UInt_t typeBits_2;
      UInt_t nTkHits_2;
      UInt_t nPixHits_2;
      UInt_t nSeg_2;
      UInt_t hltMatchBits_2;

    ClassDef(TDimuon,1)
  };
}
#endif
