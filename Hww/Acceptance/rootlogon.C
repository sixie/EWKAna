{
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }      
    
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
    
    gSystem->Load("$CMSSW_BASE/lib/slc5_ia32_gcc434/libEWKAnaNtupler.so");
  
    gROOT->Macro("$CMSSW_BASE/src/Common/CPlot.cc+");
    gROOT->Macro("$CMSSW_BASE/src/Common/MitStyleRemix.cc+");  
  }  
  
  TString rfitpath("/afs/cern.ch/cms/sw/slc5_ia32_gcc434/lcg/roofit/5.26.00-cms3/include");
  TString path = gSystem->GetIncludePath();
  path += "-I. -I$ROOTSYS/src -I";
  path += rfitpath;
  gSystem->SetIncludePath(path.Data());
       
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
