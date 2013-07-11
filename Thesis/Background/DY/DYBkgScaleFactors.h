Double_t DYBkgScaleFactor(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselection[3] = { 2.58852, 2.68315, 3.27757  };
  Double_t DYBkgScaleFactorHiggsSelection[3][18] = { 
    { 3.95802,4.61216,3.80891,3.41355,3.49018,2.82142,2.62148,2.91899,2.86186,2.58023,3.31926,1.1371,1.37976,0.852867,2.7417,1.10374,0.532186,2.85512},
    { 4.09989,4.49503,4.35328,2.84783,2.90038,3.04521,2.849,2.71709,3.25583,3.44848,3.1455,3.89679,5.94272,6.34336,4.4965,4.48271,1.44758,2.81446},
    { 2.31245,2.86642,3.11711,3.35054,3.58409,3.94173,3.51267,3.95617,4.30867,4.53273,3.1952,2.8723,3.18828,3.81704,3.34463,3.80673,1.87365,2.83505} };
  if(mH == 0) return DYBkgScaleFactorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselection[jetBin];
  }
}

Double_t DYBkgScaleFactorKappa(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[18] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190,200,250,300};
  Double_t DYBkgScaleFactorWWPreselectionKappa[3] = { 1.25702, 1.16728, 1.11468  };
  Double_t DYBkgScaleFactorHiggsSelectionKappa[3][18] = { 
    { 1.5558,1.55414,1.55471,1.506,1.52613,1.71353,1.53319,1.5262,1.51919,1.52623,2.10923,2.06111,2.02321,2.73401,1.68552,2.86116,2.74988,4.85643},
    { 1.24061,1.23017,1.22405,1.25255,1.28519,1.31632,1.34992,1.35261,1.34394,1.34162,1.49355,1.41494,1.39579,1.39409,1.24076,1.19959,1.46577,1.47737},
    { 1.31823,1.31334,1.30907,1.32662,1.3531,1.37317,1.40058,1.43173,1.42997,1.42955,1.27953,1.34484,1.41726,1.2742,1.37997,1.27351,1.60898,1.26474} };
  if(mH == 0) return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 18 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return DYBkgScaleFactorHiggsSelectionKappa[jetBin][massIndex];
  } else {
    return DYBkgScaleFactorWWPreselectionKappa[jetBin];
  }
}

