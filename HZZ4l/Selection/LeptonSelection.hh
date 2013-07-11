#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"

Int_t Classify(const mithep::TElectron *ele);
bool compute_cut(double x, double et, double cut_min, double cut_max, bool gtn = false);

Bool_t passConversionVeto(Int_t isConv);
Bool_t passMuonCuts(const mithep::TMuon *muon, Int_t Type);
Bool_t passElectronCuts(const mithep::TElectron *electron);
Bool_t isSoftMuon(const mithep::TMuon *muon);


//=== FUNCTION DEFINITIONS ======================================================================================

Bool_t passConversionVeto(Int_t isConv) {
 
  Bool_t pass = kFALSE;

//   Int_t tmp0 = floor(double(isConv) / 2.0);
//   Int_t tmp1 = floor(double(tmp0) / 2.0);
//   Int_t tmp2 = floor(double(tmp1) / 2.0);
//   Int_t tmp3 = floor(double(tmp2) / 2.0);
//   Int_t tmp4 = floor(double(tmp3) / 2.0);
//   Int_t tmp5 = floor(double(tmp4) / 2.0);
//   Int_t tmp6 = floor(double(tmp5) / 2.0);
//   Int_t tmp7 = floor(double(tmp6) / 2.0);
//   Int_t tmp8 = floor(double(tmp7) / 2.0);
//   Int_t tmp9 = floor(double(tmp8) / 2.0);
//   Int_t tmp10 = floor(double(tmp9) / 2.0);
//   Int_t tmp11 = floor(double(tmp10) / 2.0);
//   Int_t tmp12 = floor(double(tmp11) / 2.0);  
//   pass =  tmp9 % 2;

  pass = ( (UInt_t(isConv) & 512) == 512);

  return pass; 
}

Int_t Classify(const mithep::TElectron *ele) {
  
  double eta    = fabs(ele->eta);
  double eOverP = ele->EOverP;
  double fBrem  = ele->fBrem;

  int cat = -1;
  if (ele->isEB == kTRUE) {
    if ((fBrem >= 0.12) and (eOverP > 0.9) and (eOverP < 1.2))
      cat = 0;
    else if (((eta >  .445   and eta <  .45  ) or
  	      (eta >  .79    and eta <  .81  ) or
  	      (eta > 1.137   and eta < 1.157 ) or
  	      (eta > 1.47285 and eta < 1.4744)))
      cat = 6;
    else if (!ele->isEcalDriven)
      cat = 8;
    else if (fBrem < 0.12)
      cat = 1;
    else
      cat = 2;
  } else {
    if ((fBrem >= 0.2) and (eOverP > 0.82) and (eOverP < 1.22))
      cat = 3;
    else if (eta > 1.5 and eta <  1.58)
      cat = 7;
    else if (!ele->isEcalDriven)
      cat = 8;
    else if (fBrem < 0.2)
      cat = 4;
    else
      cat = 5;
  }

  return cat;
}

//--------------------------------------------------------------------------------------------------
bool compute_cut(double x, double et, double cut_min, double cut_max, bool gtn) {

  float et_min = 10;
  float et_max = 40;

  bool accept = false;
  float cut = cut_max; //  the cut at et=40 GeV

  if(et < et_max) {
    cut = cut_min + (1/et_min - 1/et)*(cut_max - cut_min)/(1/et_min - 1/et_max);
  } 
  
  if(et < et_min) {
    cut = cut_min;
  } 

  if(gtn) {   // useful for e/p cut which is gt
    accept = (x >= cut);
  } 
  else {
    accept = (x <= cut);
  }

  return accept;
}



//--------------------------------------------------------------------------------------------------
Bool_t passMuonCuts(const mithep::TMuon *muon)
{
  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->muNchi2	  > 10)    return kFALSE;
  if(muon->nMatch 	  < 2)     return kFALSE;
  if(muon->nValidHits	  < 1)     return kFALSE;
  if(muon->pterr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
 
  //isolation
  Double_t relIso03 = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  Double_t relIso05 = (muon->trkIso05 + muon->emIso05 + muon->hadIso05)/muon->pt; 
  if (relIso03 + relIso05 > 0.35) return kFALSE;
  
  //IP cut
  if (fabs(muon->d0) > 0.02) return kFALSE;



//*/
}

//--------------------------------------------------------------------------------------------------
Bool_t passElectronCuts(const mithep::TElectron *electron)
{
  //IP
  if(fabs(electron->d0) > 0.02) return kFALSE;
  if(fabs(electron->dz) > 0.1)  return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0)       return kFALSE;
  if(electron->partnerDeltaCot < 0.02 
     && electron->partnerDist < 0.02)   return kFALSE;
  //trackfit conversion rejection
  //if(! passConversionVeto(electron->isConv)) return kFALSE;
  

  //isolation
  Double_t relIso03 = (electron->trkIso03 + electron->emIso03 + electron->hadIso03)/electron->pt;
  Double_t relIso04 = (electron->trkIso04 + electron->emIso04 + electron->hadIso04)/electron->pt;
  if (relIso03 + relIso04 > 0.35) return kFALSE;


  //CiC tight (ID only == bit 1)
  if (!(UInt_t(electron->passSuperTightId) & 1 == 1)) return kFALSE;


  //VBTF80 selection
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
    if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
    if(electron->HoverE	    > 0.04)  return kFALSE;    
  } else {
    // endcap
    if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE	    > 0.10)  return kFALSE;
  }
  
  
  return kTRUE;
}


#endif

