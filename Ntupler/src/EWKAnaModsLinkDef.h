#ifndef EWKANA_NTUPLER_LINKDEF_H
#define EWKANA_NTUPLER_LINKDEF_H
#include "EWKAna/Ntupler/interface/HwwNtuplerMod.hh"
#include "EWKAna/Ntupler/interface/HwwGenNtuplerMod.hh"
#include "EWKAna/Ntupler/interface/BambuGenDumperMod.hh"
#include "EWKAna/Ntupler/interface/TEventInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TDimuon.hh"
#include "EWKAna/Ntupler/interface/TMuon.hh"
#include "EWKAna/Ntupler/interface/TDielectron.hh"
#include "EWKAna/Ntupler/interface/TElectron.hh"
#include "EWKAna/Ntupler/interface/TJet.hh"
#include "EWKAna/Ntupler/interface/TPhoton.hh"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::HwwNtuplerMod+;
#pragma link C++ class mithep::HwwGenNtuplerMod+;
#pragma link C++ class mithep::BambuGenDumperMod+;
#pragma link C++ class mithep::TEventInfo+;
#pragma link C++ class mithep::TGenInfo+;
#pragma link C++ class mithep::TDimuon+;
#pragma link C++ class mithep::TMuon+;
#pragma link C++ class mithep::TDielectron+;
#pragma link C++ class mithep::TElectron+;
#pragma link C++ class mithep::TJet+;
#pragma link C++ class mithep::TPhoton+;
#endif
