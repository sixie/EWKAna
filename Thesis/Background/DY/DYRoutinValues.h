Double_t RoutinValue(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinValuesWWPreselection[3] = { 0.195597, 0.164994, 0.161433  };
  Double_t RoutinValuesHiggsSelection[3][12] = { 
    { 0.151596,0.151596,0.507958,0.507958,0.325019,0.651159,0.67763,0.587311,0.388038,0.177564,0.0457175,0.259502},
    { 0.0794746,0.0794746,0.188197,0.188197,0.10524,0.231166,0.211369,0.194161,0.171163,0.124037,0.0568555,0.08387},
    { 0.0763378,0.0763378,0.14521,0.14521,0.0931249,0.169452,0.146787,0.143813,0.114942,0.0927258,0.0623241,0.109334} };
  if(mH == 0) return RoutinValuesWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinValuesHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinValuesWWPreselection[jetBin];
  }
}

Double_t RoutinStatError(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinStatErrorWWPreselection[3] = { 0.0192559, 0.00981837, 0.0149001  };
  Double_t RoutinStatErrorHiggsSelection[3][12] = { 
    { 0.024669,0.024669,0.086434,0.086434,0.0886281,0.266728,0.303059,0.274616,0.130994,0.0594005,0.0126074,0.0714244},
    { 0.00779779,0.00779779,0.0180301,0.0180301,0.0146813,0.0373394,0.0363474,0.0298794,0.0204553,0.014479,0.00702086,0.0121658},
    { 0.0117457,0.0117457,0.0220529,0.0220529,0.0193253,0.0408327,0.0379793,0.0335558,0.0220736,0.0169654,0.0106849,0.0204282} };
  if(mH == 0) return RoutinStatErrorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinStatErrorHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinStatErrorWWPreselection[jetBin];
  }
}

Double_t RoutinSystError(Int_t mH, Int_t jetBin) {
  Int_t mHiggs[12] = {115,120,130,140,150,160,170,180,190,200,250,300};
  Double_t RoutinSystErrorWWPreselection[3] = { 0.0145007, 0.022322, 0.0180462  };
  Double_t RoutinSystErrorHiggsSelection[3][12] = { 
    { 0.0762337,0.0762337,0.20281,0.20281,0.337348,0.245088,0.0931776,0.190155,0.117772,0.302086,0.0108007,0.144759},
    { 0.0132844,0.0132844,0.0540898,0.0540898,0.0474966,0.0655354,0.0591825,0.0561268,0.0256975,0.0121744,0.0245898,0.0365242},
    { 0.0221026,0.0221026,0.0564368,0.0564368,0.0230271,0.0995228,0.119595,0.0632419,0.0806375,0.051229,0.0575905,0.0620268} };
  if(mH == 0) return RoutinSystErrorWWPreselection[jetBin];
  Int_t massIndex = -1;
  for (UInt_t m=0; m < 12 ; ++m) {
    if (mH == mHiggs[m]) massIndex = m;
  }
  if (massIndex >= 0) {
    return RoutinSystErrorHiggsSelection[jetBin][massIndex];
  } else {
    return RoutinSystErrorWWPreselection[jetBin];
  }
}

