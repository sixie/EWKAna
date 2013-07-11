Double_t WWBkgScaleFactorCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15882,1.15882,1.16394,1.17263,1.15137,1.1544,1.15083,1.14006,1.13984},
    { 1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23408,1.23407,1.24242,1.26125,1.26628,1.26833,1.25867,1.24801,1.24782} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorHiggsSelection[2][15] = { 
    { 1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.15883,1.14287,1.15406},
    { 1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.23405,1.21773,1.26472} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaCutBased(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07678,1.0762,1.08008,1.07992,1.08039,1.08126,1.08174,},
    { 1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13058,1.13059,1.12983,1.12779,1.13189,1.13176,1.13281,1.13396,1.13441} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

Double_t WWBkgScaleFactorKappaMVA(Int_t mH, Int_t jetBin) {
assert(jetBin >= 0 && jetBin <= 1);
  Int_t mHiggs[15] = {115,118,120,122,124,126,128,130,135,140,150,160,170,180,190};
  Double_t WWBkgScaleFactorKappaHiggsSelection[2][15] = { 
    { 1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.07721,1.08211,1.08674},
    { 1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13059,1.13642,1.14047} };
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 15 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return WWBkgScaleFactorKappaHiggsSelection[jetBin][massIndex];
  } else {
    return 1.0;
  }
}

