#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/EWKAnaDefs.hh"
#include "Common/MyTools.hh"
#include "MitPhysics/ElectronLikelihood/interface/LikelihoodMeasurements.h"
#include "MitPhysics/ElectronLikelihood/interface/ElectronLikelihood.h"
#include <cassert>
#include "TMath.h"

// void computeTrkMet(const mithep::TMuon *mu1, const mithep::TMuon *mu2, 
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met=0, Double_t *out_metphi=0, Double_t *out_sumet=0);

// void computeTrkMet(const mithep::TMuon *mu, const mithep::TElectron *ele, 
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met=0, Double_t *out_metphi=0, Double_t *out_sumet=0);
		   
// void computeTrkMet(const mithep::TElectron *ele1, const mithep::TElectron *ele2,
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met=0, Double_t *out_metphi=0, Double_t *out_sumet=0);

Double_t Likelihood(const mithep::TElectron *ele, ElectronLikelihood *LH);
Bool_t passConversionVeto(Int_t isConv);
Bool_t passMuonID(const mithep::TMuon *muon);
Bool_t passMuonIDWithEACorrPFIso(const mithep::TMuon *muon, Double_t fRho);
Bool_t passMuonMVAIDIsoCombined(const mithep::TMuon *muon, Double_t mvaValue, Double_t fRho);
Bool_t passMuonIsoOnly(const mithep::TMuon *muon);
Bool_t passMuonPFIso03Only( const mithep::TMuon *muon, Double_t fRho );
Bool_t passMuonPFIso04Only( const mithep::TMuon *muon, Double_t fRho );
Bool_t passEleID(const mithep::TElectron *electron);
Bool_t passElectronLH(const mithep::TElectron *ele, Double_t LHValue);
Bool_t passElectronMVA(const mithep::TElectron *ele, Double_t mvaValue);
Bool_t passElectronMVAWithEACorrPFIso(const mithep::TElectron *ele, Double_t mvaValue, Double_t Rho);
Bool_t passElectronMVAIDIsoCombined(const mithep::TElectron *ele, Double_t mvaValue, Double_t Rho);
Bool_t passElectronProbeForMVAID(const mithep::TElectron *electron);
Bool_t passElectronPFIsoOnly(const mithep::TElectron *electron);
Bool_t passElectronEACorrPFIsoOnly(const mithep::TElectron *electron, Double_t Rho);
Bool_t passElectronMVAIDOnly(const mithep::TElectron *electron, Double_t mvaValue);
Bool_t passMuonDenominatorM2(const mithep::TMuon *muon, Double_t fRho);


Bool_t isSoftMuon(const mithep::TMuon *muon);

Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver=1);
Bool_t isEleFO(const mithep::TElectron *electron);

Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);


//*******************************************
//=== Effective Area Pileup Corrections  ====
//*******************************************
enum { kMuChargedIso03, kMuNeutralIso03, kMuChargedIso04, kMuNeutralIso04, 
       kMuHadEnergy, kMuHoEnergy, kMuEmEnergy, kMuHadS9Energy, kMuHoS9Energy, kMuEmS9Energy,
       kMuTrkIso03, kMuEMIso03, kMuHadIso03, 
       kMuTrkIso05, kMuEMIso05, kMuHadIso05 
};
Double_t MuonEffectiveArea(UInt_t type, Double_t Eta) {

  Double_t EffectiveArea = 0;
  //inclusive
//   if (type == kMuChargedIso03) EffectiveArea = 0.000;
//   if (type == kMuNeutralIso03) EffectiveArea = 0.077;
//   if (type == kMuChargedIso04) EffectiveArea = 0.000;
//   if (type == kMuNeutralIso04) EffectiveArea = 0.159;
//   if (type == kMuHadEnergy)    EffectiveArea = 0.014;
//   if (type == kMuHoEnergy)     EffectiveArea = 0.000;
//   if (type == kMuEmEnergy)     EffectiveArea = 0.000;
//   if (type == kMuHadS9Energy)  EffectiveArea = 0.053;
//   if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
//   if (type == kMuEmS9Energy)   EffectiveArea = 0.000;

  if (fabs(Eta) < 1.0) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.080;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.163;
    if (type == kMuHadEnergy)    EffectiveArea = 0.000;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.016;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.080;
    if (type == kMuHadIso03)     EffectiveArea = 0.025;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.290;
    if (type == kMuHadIso05)     EffectiveArea = 0.091;
  } else if (fabs(Eta) >= 1.0 && fabs(Eta) < 1.479 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.083;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.168;
    if (type == kMuHadEnergy)    EffectiveArea = 0.005;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.041;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.043;
    if (type == kMuHadIso03)     EffectiveArea = 0.028;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.184;
    if (type == kMuHadIso05)     EffectiveArea = 0.106;
  } else if (fabs(Eta) >= 1.479 && fabs(Eta) < 2.0 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.060;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.131;
    if (type == kMuHadEnergy)    EffectiveArea = 0.020;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.072;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.036;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.124;
    if (type == kMuHadIso05)     EffectiveArea = 0.140;
  } else if (fabs(Eta) >= 2.0 && fabs(Eta) < 2.25 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.066;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.149;
    if (type == kMuHadEnergy)    EffectiveArea = 0.056;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.148;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.025;
    if (type == kMuHadIso03)     EffectiveArea = 0.050;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.120;
    if (type == kMuHadIso05)     EffectiveArea = 0.186;
  } else if (fabs(Eta) >= 2.25 && fabs(Eta) < 2.4 ) {
    if (type == kMuChargedIso03) EffectiveArea = 0.000;
    if (type == kMuNeutralIso03) EffectiveArea = 0.098;
    if (type == kMuChargedIso04) EffectiveArea = 0.000;
    if (type == kMuNeutralIso04) EffectiveArea = 0.200;
    if (type == kMuHadEnergy)    EffectiveArea = 0.093;
    if (type == kMuHoEnergy)     EffectiveArea = 0.000;
    if (type == kMuEmEnergy)     EffectiveArea = 0.000;
    if (type == kMuHadS9Energy)  EffectiveArea = 0.260;
    if (type == kMuHoS9Energy)   EffectiveArea = 0.000;
    if (type == kMuEmS9Energy)   EffectiveArea = 0.000;
    if (type == kMuTrkIso03)     EffectiveArea = 0.000;
    if (type == kMuEMIso03)      EffectiveArea = 0.027;
    if (type == kMuHadIso03)     EffectiveArea = 0.060;
    if (type == kMuTrkIso05)     EffectiveArea = 0.000;
    if (type == kMuEMIso05)      EffectiveArea = 0.139;
    if (type == kMuHadIso05)     EffectiveArea = 0.228;
  }
  return EffectiveArea;

}

enum { kEleChargedIso03, kEleNeutralHadronIso03, kEleGammaIso03, kEleGammaIsoVetoEtaStrip03, kEleChargedIso04, kEleNeutralHadronIso04, kEleGammaIso04, kEleGammaIsoVetoEtaStrip04, kEleNeutralHadronIso007, kEleHoverE, kEleHcalDepth1OverEcal, kEleHcalDepth2OverEcal };
Double_t ElectronEffectiveArea(UInt_t type, Double_t SCEta) {

  Double_t EffectiveArea = 0;

//Inclusive Effective Areas
//   if (type == kEleChargedIso03) EffectiveArea = 0.000;
//   if (type == kEleNeutralHadronIso03) EffectiveArea = 0.021;
//   if (type == kEleGammaIso03) EffectiveArea = 0.153;
//   if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.117;
//   if (type == kEleChargedIso04) EffectiveArea = 0.000;
//   if (type == kEleNeutralHadronIso04) EffectiveArea = 0.044;
//   if (type == kEleGammaIso04) EffectiveArea = 0.183;
//   if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.118;
//   if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
//   if (type == kEleHoverE) EffectiveArea = 0.00025;
//   if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00023;
//   if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;

  if (fabs(SCEta) < 1.0) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.017;
    if (type == kEleGammaIso03) EffectiveArea = 0.045;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.014;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.034;
    if (type == kEleGammaIso04) EffectiveArea = 0.079;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.014;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00016;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(SCEta) >= 1.0 && fabs(SCEta) < 1.479 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.025;
    if (type == kEleGammaIso03) EffectiveArea = 0.052;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.030;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.050;
    if (type == kEleGammaIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.030;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00022;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00000;    
  } else if (fabs(SCEta) >= 1.479 && fabs(SCEta) < 2.0 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.030;
    if (type == kEleGammaIso03) EffectiveArea = 0.170;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.134;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.060;
    if (type == kEleGammaIso04) EffectiveArea = 0.187;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.134;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00030;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00026;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00002;        
  } else if (fabs(SCEta) >= 2.0 && fabs(SCEta) < 2.25 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.022;
    if (type == kEleGammaIso03) EffectiveArea = 0.623;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 0.516;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.055;
    if (type == kEleGammaIso04) EffectiveArea = 0.659;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 0.517;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00054;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00045;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00003;
  } else if (fabs(SCEta) >= 2.25 && fabs(SCEta) < 2.5 ) {
    if (type == kEleChargedIso03) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso03) EffectiveArea = 0.018;
    if (type == kEleGammaIso03) EffectiveArea = 1.198;
    if (type == kEleGammaIsoVetoEtaStrip03) EffectiveArea = 1.049;
    if (type == kEleChargedIso04) EffectiveArea = 0.000;
    if (type == kEleNeutralHadronIso04) EffectiveArea = 0.073;
    if (type == kEleGammaIso04) EffectiveArea = 1.258;
    if (type == kEleGammaIsoVetoEtaStrip04) EffectiveArea = 1.051;
    if (type == kEleNeutralHadronIso007) EffectiveArea = 0.000;
    if (type == kEleHoverE) EffectiveArea = 0.00082;
    if (type == kEleHcalDepth1OverEcal) EffectiveArea = 0.00066;
    if (type == kEleHcalDepth2OverEcal) EffectiveArea = 0.00004;
  }
  
  
  return EffectiveArea;

}


//=== FUNCTION DEFINITIONS ======================================================================================

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


// //--------------------------------------------------------------------------------------------------
// void computeTrkMet(const mithep::TMuon *mu1, const mithep::TMuon *mu2,
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met, Double_t *out_metphi, Double_t *out_sumet)
// {
//   TLorentzVector met;
//   met.SetPxPyPzE(in_met*cos(in_metphi) + mu1->pfPx + mu2->pfPx - (mu1->pt)*cos(mu1->phi) - (mu2->pt)*cos(mu2->phi),
//                  in_met*sin(in_metphi) + mu1->pfPy + mu2->pfPy - (mu1->pt)*sin(mu1->phi) - (mu2->pt)*sin(mu2->phi),
// 		 0,0);
//   Double_t sumet = in_sumet - sqrt((mu1->pfPx)*(mu1->pfPx) + (mu1->pfPy)*(mu1->pfPy)) - sqrt((mu1->pfPx)*(mu1->pfPx) + (mu1->pfPy)*(mu1->pfPy)) 
//                    + mu1->pt + mu2->pt;
  
//   if(out_met)    *out_met    = met.Pt();
//   if(out_metphi) *out_metphi = met.Phi();
//   if(out_sumet)  *out_sumet  = sumet;
// }

// void computeTrkMet(const mithep::TMuon *mu, const mithep::TElectron *ele, 
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met, Double_t *out_metphi, Double_t *out_sumet)
// {
//   TLorentzVector met;
//   met.SetPxPyPzE(in_met*cos(in_metphi) + mu->pfPx + ele->pfPx - (mu->pt)*cos(mu->phi) - (ele->pt)*cos(ele->phi),
//                  in_met*sin(in_metphi) + mu->pfPy + ele->pfPy - (mu->pt)*sin(mu->phi) - (ele->pt)*sin(ele->phi),
// 		 0,0);
//   Double_t sumet = in_sumet - sqrt((mu->pfPx)*(mu->pfPx) + (mu->pfPy)*(mu->pfPy)) - sqrt((ele->pfPx)*(ele->pfPx) + (ele->pfPy)*(ele->pfPy)) 
//                    + mu->pt + ele->pt;
  
//   if(out_met)    *out_met    = met.Pt();
//   if(out_metphi) *out_metphi = met.Phi();
//   if(out_sumet)  *out_sumet  = sumet;
// }

// void computeTrkMet(const mithep::TElectron *ele1, const mithep::TElectron *ele2, 
//                    const Double_t in_met, const Double_t in_metphi, const Double_t in_sumet,
//                    Double_t *out_met, Double_t *out_metphi, Double_t *out_sumet)
// {
//   TLorentzVector met;
//   met.SetPxPyPzE(in_met*cos(in_metphi) + ele1->pfPx + ele2->pfPx - (ele1->pt)*cos(ele1->phi) - (ele2->pt)*cos(ele2->phi),
//                  in_met*sin(in_metphi) + ele1->pfPy + ele2->pfPy - (ele1->pt)*sin(ele1->phi) - (ele2->pt)*sin(ele2->phi),
// 		 0,0);
//   Double_t sumet = in_sumet - sqrt((ele1->pfPx)*(ele1->pfPx) + (ele1->pfPy)*(ele1->pfPy)) - sqrt((ele2->pfPx)*(ele2->pfPx) + (ele2->pfPy)*(ele2->pfPy)) 
//                    + ele1->pt + ele2->pt;
  
//   if(out_met)    *out_met    = met.Pt();
//   if(out_metphi) *out_metphi = met.Phi();
//   if(out_sumet)  *out_sumet  = sumet;
// }

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const mithep::TMuon *muon)
{

  Double_t pfIso03 = muon->ChargedIso03 + muon->NeutralIso03_10Threshold;

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(muon->TrkKink        >= 20)   return kFALSE;

/*
if(muon->muNchi2    > 10)        return kFALSE;
if(muon->nMatch     < 2)         return kFALSE;
if(muon->nValidHits < 1)         return kFALSE;
if(!(muon->typeBits & kTracker)) return kFALSE;
if(!(muon->typeBits & kGlobal))  return kFALSE;
*/
  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(muon->pt>20) {
    if(fabs(muon->d0)>0.02)   return kFALSE;
    if(fabs(muon->eta)<1.479) return (pfIso03<0.13*(muon->pt));
    else                      return (pfIso03<0.09*(muon->pt));
  } else {
    if(fabs(muon->d0)>0.01)   return kFALSE;
    if(fabs(muon->eta)<1.479) return (pfIso03<0.06*(muon->pt));
    else                      return (pfIso03<0.05*(muon->pt));
  }
}

//--------------------------------------------------------------------------------------------------
Bool_t passMuonIDWithEACorrPFIso(const mithep::TMuon *muon, Double_t fRho)
{

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(muon->TrkKink        >= 20)   return kFALSE;

/*
if(muon->muNchi2    > 10)        return kFALSE;
if(muon->nMatch     < 2)         return kFALSE;
if(muon->nValidHits < 1)         return kFALSE;
if(!(muon->typeBits & kTracker)) return kFALSE;
if(!(muon->typeBits & kGlobal))  return kFALSE;
*/
  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(muon->pt>20) {
    if(fabs(muon->d0)>0.02)   return kFALSE;
  } else {
    if(fabs(muon->d0)>0.01)   return kFALSE;
  }

  Bool_t passIso = passMuonPFIso03Only(muon, fRho);
  if (!passIso) return kFALSE;

  return kTRUE;

}

//--------------------------------------------------------------------------------------------------
Bool_t passMuonCutBasedIDOnly(const mithep::TMuon *muon)
{

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(muon->TrkKink        >= 20)   return kFALSE;

  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(muon->pt>20) {
    if(fabs(muon->d0)>0.02)   return kFALSE;
  } else {
    if(fabs(muon->d0)>0.01)   return kFALSE;
  }
  return kTRUE;
}



//--------------------------------------------------------------------------------------------------
Bool_t passMuonIsoOnly(const mithep::TMuon *muon)
{

  Double_t pfIso03 = muon->ChargedIso03 + muon->NeutralIso03_10Threshold;

  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(muon->pt>20) {
    if(fabs(muon->eta)<1.479) return (pfIso03<0.13*(muon->pt));
    else                      return (pfIso03<0.09*(muon->pt));
  } else {
    if(fabs(muon->eta)<1.479) return (pfIso03<0.06*(muon->pt));
    else                      return (pfIso03<0.05*(muon->pt));
  }
}


Bool_t passMuonPFIso04Only( const mithep::TMuon *muon , Double_t fRho) {

  Int_t subdet = 0;
  if (fabs(muon->eta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (muon->pt > 14.5) ptBin = 1;
  if (muon->pt > 20.0) ptBin = 2;

  Int_t PFIsoBin = -1;
  if (subdet == 0 && ptBin == 0) PFIsoBin = 0;
  if (subdet == 1 && ptBin == 0) PFIsoBin = 1;
  if (subdet == 0 && ptBin == 1) PFIsoBin = 2;
  if (subdet == 1 && ptBin == 1) PFIsoBin = 3;
  if (subdet == 0 && ptBin == 2) PFIsoBin = 4;
  if (subdet == 1 && ptBin == 2) PFIsoBin = 5;

  Double_t PFIsoCut = -999;
  if (PFIsoBin == 0) PFIsoCut = 0.089;
  if (PFIsoBin == 1) PFIsoCut = 0.083;
  if (PFIsoBin == 2) PFIsoCut = 0.081;
  if (PFIsoBin == 3) PFIsoCut = 0.069;
  if (PFIsoBin == 4) PFIsoCut = 0.183;
  if (PFIsoBin == 5) PFIsoCut = 0.120;

  Double_t pfiso = muon->ChargedIso04 + muon->NeutralIso04_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso04, muon->eta);

  if (pfiso/muon->pt < PFIsoCut) return kTRUE;
  return kFALSE;
}

Bool_t passMuonPFIso03Only( const mithep::TMuon *muon, Double_t fRho ) {

  Int_t subdet = 0;
  if (fabs(muon->eta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (muon->pt > 14.5) ptBin = 1;
  if (muon->pt > 20.0) ptBin = 2;

  Int_t PFIsoBin = -1;
  if (subdet == 0 && ptBin == 0) PFIsoBin = 0;
  if (subdet == 1 && ptBin == 0) PFIsoBin = 1;
  if (subdet == 0 && ptBin == 1) PFIsoBin = 2;
  if (subdet == 1 && ptBin == 1) PFIsoBin = 3;
  if (subdet == 0 && ptBin == 2) PFIsoBin = 4;
  if (subdet == 1 && ptBin == 2) PFIsoBin = 5;

  Double_t PFIsoCut = -999;
  if (PFIsoBin == 0) PFIsoCut = 0.047;
  if (PFIsoBin == 1) PFIsoCut = 0.042;
  if (PFIsoBin == 2) PFIsoCut = 0.049;
  if (PFIsoBin == 3) PFIsoCut = 0.038;
  if (PFIsoBin == 4) PFIsoCut = 0.128;
  if (PFIsoBin == 5) PFIsoCut = 0.082;

  Double_t pfiso = muon->ChargedIso03 + muon->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03, muon->eta);
  if (pfiso/muon->pt < PFIsoCut) return kTRUE;
  return kFALSE;

}




Bool_t passMuonMVASameCutBasedSig(const mithep::TMuon *mu, Double_t mvaValue, Double_t fRho) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(mu->eta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (mu->pt > 14.5) ptBin = 1;
  if (mu->pt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
  //same Eff as LP2011 using HWW MC (no H/E, proper E/P)
  if (MVABin == 0) MVACut = -0.3642;
  if (MVABin == 1) MVACut = 0.103;
  if (MVABin == 2) MVACut = -0.3054;
  if (MVABin == 3) MVACut = 0.0398;
  if (MVABin == 4) MVACut = 0.6838;
  if (MVABin == 5) MVACut = 0.7954;

  //Isolation
  Double_t pfIso03 = mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03, mu->eta);

  //Explicitly Apply M2 Denominator Cuts
  Bool_t pass = kTRUE;
  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) >= 2.4) pass = kFALSE;
  
  if (! ( (0==0)
          &&
          (
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
          && (pfIso03 / mu->pt) < 0.4
          && ( mu->pterr / mu->pt < 0.1)
          
          && mvaValue > MVACut
          && (pfIso03 / mu->pt) < 0.4
        )
    ) {
    pass = kFALSE;
  }

  return pass;
}

//Half Cut-Based Bkg
Bool_t passMuonMVAIDIsoCombined(const mithep::TMuon *mu, Double_t mvaValue, Double_t fRho) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(mu->eta) < 1.479) subdet = 0;
  else subdet = 1;
  Int_t ptBin = 0;
  if (mu->pt > 14.5) ptBin = 1;
  if (mu->pt > 20.0) ptBin = 2;

  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 0 && ptBin == 1) MVABin = 2;
  if (subdet == 1 && ptBin == 1) MVABin = 3;
  if (subdet == 0 && ptBin == 2) MVABin = 4;
  if (subdet == 1 && ptBin == 2) MVABin = 5;

  Double_t MVACut = -999;
//   //50% bkg eff as cut-based (using V8 - PFISO)
//   if (MVABin == 0) MVACut = -0.377;
//   if (MVABin == 1) MVACut = -0.0902;
//   if (MVABin == 2) MVACut = -0.221;
//   if (MVABin == 3) MVACut = -0.0154;
//   if (MVABin == 4) MVACut = 0.459;
//   if (MVABin == 5) MVACut = 0.9158;

//   //same signal eff as cut-based (using V8 - PFISO)
//   if (MVABin == 0) MVACut = -0.5514;
//   if (MVABin == 1) MVACut = -0.303;
//   if (MVABin == 2) MVACut = -0.4562;
//   if (MVABin == 3) MVACut = -0.269;
//   if (MVABin == 4) MVACut = 0.1726;
//   if (MVABin == 5) MVACut = 0.801;

  //same signal eff as cut-based (using V10 - Detector Based Iso)
  if (MVABin == 0) MVACut = -0.5618;
  if (MVABin == 1) MVACut = -0.3002;
  if (MVABin == 2) MVACut = -0.4642;
  if (MVABin == 3) MVACut = -0.2478;
  if (MVABin == 4) MVACut = 0.1706;
  if (MVABin == 5) MVACut = 0.8146;

  //Isolation
  Double_t pfIso03 = mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03, mu->eta);

  //Explicitly Apply M2 Denominator Cuts
  Bool_t pass = kTRUE;
  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) >= 2.4) pass = kFALSE;
  
  if (! ( (0==0)
          &&
          (
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
          && (pfIso03 / mu->pt) < 0.4
          && ( mu->pterr / mu->pt < 0.1)
          
          && mvaValue > MVACut
          && (pfIso03 / mu->pt) < 0.4
        )
    ) {
    pass = kFALSE;
  }

  return pass;
}


//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const mithep::TElectron *electron)
{
  Double_t pfIso04 = electron->ChargedIso04+electron->NeutralHadronIso04_10Threshold
    +electron->GammaIso04_10Threshold-electron->GammaIsoVetoEtaStrip04_10Threshold
    -electron->ChargedEMIsoVetoEtaStrip04-electron->NeutralHadronIso007_10Threshold;

  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0)           return kFALSE;
  if(!passConversionVeto(electron->isConv)) return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(pfIso04 > 0.13*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
      if(electron->HoverE	    > 0.04)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.025) return kFALSE;    
    }
  
  } else {
    // endcap
    if(pfIso04 > 0.09*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;      
    }
  }

/*** NOTE: change to SC eta! ***/  
  if(electron->pt < 20)
    return ((electron->fBrem>0.15) || (fabs(electron->scEta)<1 && electron->EOverP>0.95));

  return kTRUE;
}


Double_t Likelihood(const mithep::TElectron *ele, ElectronLikelihood *LH) {


  //Get likelihood value
  LikelihoodMeasurements measurements;
  measurements.pt = ele->pt;
  if (ele->isEB && (fabs(ele->eta)<1.0)) measurements.subdet = 0;
  else if (ele->isEB) measurements.subdet = 1;
  else measurements.subdet = 2;
  measurements.deltaPhi = ele->deltaPhiIn;
  measurements.deltaEta = ele->deltaEtaIn;
  measurements.eSeedClusterOverPout = ele->ESeedClusterOverPout;
  measurements.eSuperClusterOverP = ele->EOverP;
  measurements.hadronicOverEm = ele->HoverE;
  measurements.sigmaIEtaIEta = ele->sigiEtaiEta;
  measurements.sigmaIPhiIPhi = TMath::Sqrt(ele->sigiPhiiPhi);
  measurements.fBrem = ele->fBrem;
  measurements.nBremClusters = ele->nBrem;
  measurements.OneOverEMinusOneOverP = (1.0/(ele->EOverP*ele->p)) - 1.0 / ele->p;
  double likelihood = LH->result(measurements);
  Double_t likelihoodValue = 0; 
  if (likelihood <= 0) {
    likelihoodValue = -20.0;
  } else if (likelihood == 1) {
    likelihoodValue = 20.0;
  } else {
    likelihoodValue = TMath::Log(likelihood / (1.0-likelihood));
  }

  return likelihoodValue;
}


Bool_t passElectronLH(const mithep::TElectron *ele, ElectronLikelihood *LH) {

  Double_t LHValue = Likelihood(ele,LH);

  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
      +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
      -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_10Threshold;
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold;
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
  
  Double_t LHCutValue = 0;
  //NEW Tight WP
    if(ele->pt > 20){
      if(fabs(ele->scEta) < 1.479){
        if(ele->nBrem == 0)           LHCutValue = 3.5;
	else                          LHCutValue = 4.0;
      }
      else  {                                
        if(ele->nBrem == 0)           LHCutValue = 4.0;
	else                          LHCutValue = 4.0;
      }
    }
    else {
      if(fabs(ele->scEta) < 1.479){
        if(ele->nBrem == 0)           LHCutValue =  4.0;
	else                          LHCutValue =  4.5;
      }
      else  {                                
        if(ele->nBrem == 0)           LHCutValue =  4.0;
	else                          LHCutValue =  4.0;
      }
    }


  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && LHValue > LHCutValue

            && ( ele->trkIso03 ) / ele->pt < 0.2
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
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
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && LHValue > LHCutValue

             && (ele->trkIso03 ) / ele->pt < 0.2
             && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}






Bool_t passElectronMVA(const mithep::TElectron *ele, Double_t mvaValue, Int_t Option) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 1.0) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20.0) ptBin = 1;
  
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;


  //Isolation
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
      +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
      -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_10Threshold;
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold;
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

  Double_t MVACut = -999;

  //WP with same Eff as LP2011 Cut based
  if (Option == 1) {
    if (MVABin == 0) MVACut = 0.133;
    if (MVABin == 1) MVACut = 0.465;
    if (MVABin == 2) MVACut = 0.518; 
    if (MVABin == 3) MVACut = 0.942;
    if (MVABin == 4) MVACut = 0.947;
    if (MVABin == 5) MVACut = 0.878 ;
  } else if (Option == 2) {
    if (MVABin == 0) MVACut = 0.139;
    if (MVABin == 1) MVACut = 0.525;
    if (MVABin == 2) MVACut = 0.543; 
    if (MVABin == 3) MVACut = 0.947;
    if (MVABin == 4) MVACut = 0.950;
    if (MVABin == 5) MVACut = 0.884;
  }
  
  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && mvaValue > MVACut

            && ( ele->trkIso03 ) / ele->pt < 0.2
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
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
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && mvaValue > MVACut

             && (ele->trkIso03 ) / ele->pt < 0.2
             && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}


Bool_t passElectronMVAWithEACorrPFIso(const mithep::TElectron *ele, Double_t mvaValue, Int_t Option, Double_t Rho) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 1.0) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20.0) ptBin = 1;
  
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;


  //Isolation
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_05Threshold
    +ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold
    -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_05Threshold
    - Rho*ElectronEffectiveArea(kEleNeutralHadronIso03,ele->scEta ) 
    - Rho*ElectronEffectiveArea(kEleGammaIso03,ele->scEta)
    + Rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,ele->scEta)
    + Rho*ElectronEffectiveArea(kEleNeutralHadronIso007,ele->scEta);
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_05Threshold
    +ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_05Threshold
    - Rho*ElectronEffectiveArea(kEleNeutralHadronIso04,ele->scEta)
    - Rho*ElectronEffectiveArea(kEleGammaIso04,ele->scEta)
    + Rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,ele->scEta)
    + Rho*ElectronEffectiveArea(kEleNeutralHadronIso007,ele->scEta);
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

  Double_t MVACut = -999;

  //WP with same Eff as LP2011 Cut based
  if (Option == 1) {
    if (MVABin == 0) MVACut = 0.133;
    if (MVABin == 1) MVACut = 0.465;
    if (MVABin == 2) MVACut = 0.518; 
    if (MVABin == 3) MVACut = 0.942;
    if (MVABin == 4) MVACut = 0.947;
    if (MVABin == 5) MVACut = 0.878 ;
  } else if (Option == 2) {
    if (MVABin == 0) MVACut = 0.139;
    if (MVABin == 1) MVACut = 0.525;
    if (MVABin == 2) MVACut = 0.543; 
    if (MVABin == 3) MVACut = 0.947;
    if (MVABin == 4) MVACut = 0.950;
    if (MVABin == 5) MVACut = 0.884;
  }
  
  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && iso / ele->pt < isoCutValue
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && mvaValue > MVACut

            && ( ele->trkIso03 ) / ele->pt < 0.2
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
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
             && iso / ele->pt < isoCutValue
             && ele->nExpHitsInner <= 0
             && passConversionVeto(ele->isConv)
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && mvaValue > MVACut

             && (ele->trkIso03 ) / ele->pt < 0.2
             && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}



Bool_t passElectronMVAIDOnly(const mithep::TElectron *ele, Double_t mvaValue, Int_t Option) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 1.0) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20.0) ptBin = 1;
  
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;

  //WP with same Eff as LP2011 Cut based
  if (Option == 1) {
    if (MVABin == 0) MVACut = 0.133;
    if (MVABin == 1) MVACut = 0.465;
    if (MVABin == 2) MVACut = 0.518; 
    if (MVABin == 3) MVACut = 0.942;
    if (MVABin == 4) MVACut = 0.947;
    if (MVABin == 5) MVACut = 0.878 ;
  } else if (Option == 2) {
    if (MVABin == 0) MVACut = 0.139;
    if (MVABin == 1) MVACut = 0.525;
    if (MVABin == 2) MVACut = 0.543; 
    if (MVABin == 3) MVACut = 0.947;
    if (MVABin == 4) MVACut = 0.950;
    if (MVABin == 5) MVACut = 0.884;
  }
  
  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && mvaValue > MVACut

            && ( ele->trkIso03 ) / ele->pt < 0.2
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
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
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && mvaValue > MVACut

             && (ele->trkIso03 ) / ele->pt < 0.2
             && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}





Bool_t passElectronMVAIDIsoCombined(const mithep::TElectron *ele, Double_t mvaValue, Int_t Option) {

  //Find MVA Bin
  Int_t subdet = 0;
  if (fabs(ele->scEta) < 1.0) subdet = 0;
  else if (fabs(ele->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (ele->pt > 20.0) ptBin = 1;
  
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -999;

  //WP with same Eff as LP2011 Cut based
  if (MVABin == 0) MVACut = 0.4202;
  if (MVABin == 1) MVACut = 0.6206;
  if (MVABin == 2) MVACut = 0.619; 
  if (MVABin == 3) MVACut = 0.959;
  if (MVABin == 4) MVACut = 0.9586;
  if (MVABin == 5) MVACut = 0.9278;
  
  //Explicitly Apply V4 Denominator Cuts

  Bool_t pass = kTRUE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && passConversionVeto(ele->isConv)
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
            && mvaValue > MVACut

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
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
             && mvaValue > MVACut

             && (ele->trkIso03) / ele->pt < 0.20
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20

          )
      ) {
      pass = kFALSE;
    }
  } 


  return pass;
}





Bool_t passElectronPFIsoOnly(const mithep::TElectron *ele) {

  //Isolation
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_10Threshold
      +ele->GammaIso03_10Threshold-ele->GammaIsoVetoEtaStrip03_10Threshold
      -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_10Threshold;
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_10Threshold
    +ele->GammaIso04_10Threshold-ele->GammaIsoVetoEtaStrip04_10Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_10Threshold;
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


  Bool_t pass = kTRUE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  if (! ( (0==0)
          && iso / ele->pt < isoCutValue
          && ( ele->trkIso03 ) / ele->pt < 0.2
          && (ele->emIso03) / ele->pt < 0.20
          && (ele->hadIso03) / ele->pt < 0.20
        )
    ) {
    pass = kFALSE;
  }      

  return pass;
}

Bool_t passElectronEACorrPFIsoOnly(const mithep::TElectron *ele, Double_t Rho) {

  //Isolation
  Double_t iso03 = ele->ChargedIso03+ele->NeutralHadronIso03_05Threshold
    +ele->GammaIso03_05Threshold-ele->GammaIsoVetoEtaStrip03_05Threshold
    -ele->ChargedEMIsoVetoEtaStrip03-ele->NeutralHadronIso007_05Threshold 
    - Rho*ElectronEffectiveArea(kEleNeutralHadronIso03, ele->scEta) 
    - Rho*ElectronEffectiveArea(kEleGammaIso03, ele->scEta)
    + Rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03, ele->scEta)
    + Rho*ElectronEffectiveArea(kEleNeutralHadronIso007, ele->scEta);
  Double_t iso04 = ele->ChargedIso04+ele->NeutralHadronIso04_05Threshold
    +ele->GammaIso04_05Threshold-ele->GammaIsoVetoEtaStrip04_05Threshold
    -ele->ChargedEMIsoVetoEtaStrip04-ele->NeutralHadronIso007_05Threshold
    - Rho*ElectronEffectiveArea(kEleNeutralHadronIso04, ele->scEta)
    - Rho*ElectronEffectiveArea(kEleGammaIso04, ele->scEta)
    + Rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04, ele->scEta)
    + Rho*ElectronEffectiveArea(kEleNeutralHadronIso007, ele->scEta);
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


  Bool_t pass = kTRUE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;
  
  if (! ( (0==0)
          && iso / ele->pt < isoCutValue
          && ( ele->trkIso03 ) / ele->pt < 0.2
          && (ele->emIso03) / ele->pt < 0.20
          && (ele->hadIso03) / ele->pt < 0.20
        )
    ) {
    pass = kFALSE;
  }      

  return pass;
}


//--------------------------------------------------------------------------------------------------
//Passes conversion veto, d0 , isolation
Bool_t passElectronProbeForMVAID(const mithep::TElectron *electron)
{

  Double_t pfIso04 = electron->ChargedIso04+electron->NeutralHadronIso04_10Threshold
    +electron->GammaIso04_10Threshold-electron->GammaIsoVetoEtaStrip04_10Threshold
    -electron->ChargedEMIsoVetoEtaStrip04-electron->NeutralHadronIso007_10Threshold;

  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel     
    if(pfIso04 > 0.13*(electron->pt)) return kFALSE;
    
  } else {
    // endcap 
    if(pfIso04 > 0.09*(electron->pt)) return kFALSE;    
  }
  
  return kTRUE;
}




//--------------------------------------------------------------------------------------------------
Bool_t isSoftMuon(const mithep::TMuon *muon)
{
  if(muon->nTkHits  < 11)  return kFALSE;
  if(fabs(muon->d0) > 0.2) return kFALSE;
  if(fabs(muon->dz) > 0.1) return kFALSE;

  if(!(muon->typeBits & kTracker)) return kFALSE;  

  if(!(muon->qualityBits & kTMLastStationAngTight)) return kFALSE;
	  
  Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  if(muon->pt>20 && iso<0.1) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver)
{

  Double_t pfIso03 = muon->ChargedIso03 + muon->NeutralIso03_10Threshold;

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(muon->TrkKink        >= 20)   return kFALSE;

/*
if(muon->muNchi2    > 10)        return kFALSE;
if(muon->nMatch     < 2)         return kFALSE;
if(muon->nValidHits < 1)         return kFALSE;
if(!(muon->typeBits & kTracker)) return kFALSE;
if(!(muon->typeBits & kGlobal))  return kFALSE;
*/
  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(fabs(muon->d0) > 0.2) return kFALSE;
  
  if(ver==1) return (pfIso03/muon->pt<1.0);
  if(ver==2) return (pfIso03/muon->pt<0.4);
  if(ver==3) return (muon->trkIso03/muon->pt<0.2 && muon->emIso03/muon->pt<0.2 && muon->hadIso03/muon->pt<0.2);
  
  return kFALSE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isEleFO(const mithep::TElectron *electron)
{
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;

  //IP cut
  if(fabs(electron->d0) > 0.02) return kFALSE;
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {  
    // barrel
    if(electron->sigiEtaiEta      > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.15)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE	          > 0.12)  return kFALSE;

    if(electron->trkIso03                         > 0.2*(electron->pt)) return kFALSE;
    if(TMath::Max(electron->emIso03-1,Float_t(0)) > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03                         > 0.2*(electron->pt)) return kFALSE;
        
  } else {
    // endcap
    if(electron->sigiEtaiEta	  > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.10)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.009) return kFALSE;
    if(electron->HoverE	          > 0.10)  return kFALSE;

    if(electron->trkIso03 > 0.2*(electron->pt)) return kFALSE;
    if(electron->emIso03  > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03 > 0.2*(electron->pt)) return kFALSE;
  }
    
  return kTRUE;
}
#endif

//--------------------------------------------------------------------------------------------------
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi) 
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = toolbox::deltaPhi(lepPhi,metPhi);
  if(dphi > 0.5*pi)
    return met;
    
  return met*sin(dphi);
}

Bool_t passMuonDenominatorM2(const mithep::TMuon *mu, Double_t fRho) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

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
        && (mu->ChargedIso03 + mu->NeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03, mu->eta)) / mu->pt < 0.4          
        && ( mu->pterr / mu->pt < 0.1)
        && (mu->TrkKink < 20)
        )
    ) pass = kFALSE;    



  return pass;
}
