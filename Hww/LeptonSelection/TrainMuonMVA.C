// @(#)root/tmva $Id: TrainMuonMVA.C,v 1.3 2012/02/07 13:35:47 sixie Exp $
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <set>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TChainElement.h"

#include "TMVAGui.C"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "EWKAna/Utils/LeptonIDCuts.hh"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void TrainMuonMVA(
 TString sigInputFile    = "MuonSelectionTraining.Real.weighted.root",
 TString bgdInputFile    = "MuonSelectionTraining.Fake.weighted.root",
 TString label           = "default",
 string  versionLabel    = "",
 Bool_t  buildCommittee  = kFALSE,
 Int_t Option            = -1,
// TString myMethodList    = "Likelihood,BDT,BDTG"
 TString myMethodList    = "Likelihood"
 )
{
  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the
  // corresponding lines from .rootrc

  // methods to be processed can be given as an argument; use format:
  //
  // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
  //
  // if you like to use a method via the plugin mechanism, we recommend using
  //
  // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
  // (an example is given for using the BDT as plugin (see below),
  // but of course the real application is when you write your own
  // method based)

   
  //-------------------------------------------------------------------------------------
  // define event selection
  //
  // current selection is: WW selection, except anti b-tagging and soft muon veto
  //
  // The event selection is applied below, see: 
  // SIGNAL EVENT SELECTION and BACKGROUND EVENT SELECTION
  // if you want to apply any additional selection this needs to be implemented below
  //-------------------------------------------------------------------------------------

  //-----------------------------------------------------
  // choose which variables to include in MVA training
  //-----------------------------------------------------
  
  std::map<std::string,int> mvaVar;
  mvaVar[ "TkNchi2" ]              = 0;  
  mvaVar[ "GlobalNchi2" ]          = 0; 
  mvaVar[ "NValidHits" ]           = 0;  
  mvaVar[ "NTrackerHits" ]         = 0; 
  mvaVar[ "NPixelHits" ]           = 0;  
  mvaVar[ "NMatches" ]             = 0; 
  mvaVar[ "D0" ]                   = 0; 
  mvaVar[ "IP3d" ]                 = 0;
  mvaVar[ "IP3dSig" ]              = 0;
  mvaVar[ "TrkKink" ]              = 0;  
  mvaVar[ "GlobalKink" ]           = 0;  
  mvaVar[ "SegmentCompatibility" ] = 0;  
  mvaVar[ "CaloCompatibility" ]    = 0;  
  mvaVar[ "HadEnergyOverPt" ]      = 0;  
  mvaVar[ "HoEnergyOverPt" ]       = 0;  
  mvaVar[ "EmEnergyOverPt" ]       = 0;  
  mvaVar[ "HadS9EnergyOverPt" ]    = 0;  
  mvaVar[ "HoS9EnergyOverPt" ]     = 0;  
  mvaVar[ "EmS9EnergyOverPt" ]     = 0;  
  mvaVar[ "ChargedIso03OverPt" ]   = 0;  
  mvaVar[ "NeutralIso03OverPt" ]   = 0;  
  mvaVar[ "ChargedIso04OverPt" ]   = 0;  
  mvaVar[ "NeutralIso04OverPt" ]   = 0;  
  mvaVar[ "TrkIso03OverPt" ]       = 0;  
  mvaVar[ "EMIso03OverPt" ]        = 0;  
  mvaVar[ "HadIso03OverPt" ]       = 0;  
  mvaVar[ "TrkIso05OverPt" ]       = 0;  
  mvaVar[ "EMIso05OverPt" ]        = 0;  
  mvaVar[ "HadIso05OverPt" ]       = 0;  
  mvaVar[ "NVtx" ]                 = 0;  
 
  if (versionLabel == "V0") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 0;
    mvaVar[ "IP3dSig" ]              = 0;
    mvaVar[ "TrkKink" ]              = 0;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 0;  
    mvaVar[ "CaloCompatibility" ]    = 0;  
    mvaVar[ "HadEnergyOverPt" ]      = 0;  
    mvaVar[ "EmEnergyOverPt" ]       = 0;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 0;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 0;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V1") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 0;  
    mvaVar[ "CaloCompatibility" ]    = 0;  
    mvaVar[ "HadEnergyOverPt" ]      = 0;  
    mvaVar[ "EmEnergyOverPt" ]       = 0;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 0;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 0;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V2") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 0;  
    mvaVar[ "HadEnergyOverPt" ]      = 0;  
    mvaVar[ "EmEnergyOverPt" ]       = 0;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 0;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 0;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V3") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 0;  
    mvaVar[ "EmEnergyOverPt" ]       = 0;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 0;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 0;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V4") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V5") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 1;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V6") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 1;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 1;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V7") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 1;  
    mvaVar[ "NeutralIso03OverPt" ]   = 1;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V8") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 1;  
    mvaVar[ "NeutralIso03OverPt" ]   = 1;  
    mvaVar[ "ChargedIso04OverPt" ]   = 1;  
    mvaVar[ "NeutralIso04OverPt" ]   = 1;  
    mvaVar[ "TrkIso03OverPt" ]       = 0;  
    mvaVar[ "EMIso03OverPt" ]        = 0;  
    mvaVar[ "HadIso03OverPt" ]       = 0;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V9") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 1;  
    mvaVar[ "EMIso03OverPt" ]        = 1;  
    mvaVar[ "HadIso03OverPt" ]       = 1;  
    mvaVar[ "TrkIso05OverPt" ]       = 0;  
    mvaVar[ "EMIso05OverPt" ]        = 0;  
    mvaVar[ "HadIso05OverPt" ]       = 0;  
  }

  if (versionLabel == "V10") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 1;  
    mvaVar[ "EMIso03OverPt" ]        = 1;  
    mvaVar[ "HadIso03OverPt" ]       = 1;  
    mvaVar[ "TrkIso05OverPt" ]       = 1;  
    mvaVar[ "EMIso05OverPt" ]        = 1;  
    mvaVar[ "HadIso05OverPt" ]       = 1;  
  }

  if (versionLabel == "V11") {
    mvaVar[ "TkNchi2" ]              = 1;  
    mvaVar[ "GlobalNchi2" ]          = 1; 
    mvaVar[ "NValidHits" ]           = 1;  
    mvaVar[ "NTrackerHits" ]         = 1; 
    mvaVar[ "NPixelHits" ]           = 1;  
    mvaVar[ "NMatches" ]             = 1; 
    mvaVar[ "D0" ]                   = 1; 
    mvaVar[ "IP3d" ]                 = 1;
    mvaVar[ "IP3dSig" ]              = 1;
    mvaVar[ "TrkKink" ]              = 1;  
    mvaVar[ "GlobalKink" ]           = 0;  
    mvaVar[ "SegmentCompatibility" ] = 1;  
    mvaVar[ "CaloCompatibility" ]    = 1;  
    mvaVar[ "HadEnergyOverPt" ]      = 1;  
    mvaVar[ "EmEnergyOverPt" ]       = 1;  
    mvaVar[ "HadS9EnergyOverPt" ]    = 1;  
    mvaVar[ "EmS9EnergyOverPt" ]     = 1;  
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
    mvaVar[ "TrkIso03OverPt" ]       = 1;  
    mvaVar[ "EMIso03OverPt" ]        = 1;  
    mvaVar[ "HadIso03OverPt" ]       = 1;  
    mvaVar[ "TrkIso05OverPt" ]       = 1;  
    mvaVar[ "EMIso05OverPt" ]        = 1;  
    mvaVar[ "HadIso05OverPt" ]       = 1;  
    mvaVar[ "NVtx" ]                 = 1;  
  }



  //---------------------------------
  //choose bkg samples to include
  //---------------------------------
  
  TChain *chbackground = new TChain("Muons");
  chbackground->Add(bgdInputFile);

  //---------------------------------
  //choose signal sample to include
  //---------------------------------

  TChain *chsignal = new TChain("Muons");
  chsignal->Add(sigInputFile);

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // --- Cut optimisation
  Use["Cuts"]            = 1;
  Use["CutsD"]           = 1;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  // 
  // --- 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 1;
  Use["LikelihoodD"]     = 1; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // --- Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 1;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 1;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 1; // k-nearest neighbour method
  //
  // --- Linear Discriminant Analysis
  Use["LD"]              = 1; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // --- Function Discriminant analysis
  Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  //
  // --- Support Vector Machine 
  Use["SVM"]             = 1;
  // 
  // --- Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 1; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  // 
  // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 1;
  //
  // --- multi-output MVA's
  Use["multi_BDTG"]      = 1;
  Use["multi_MLP"]       = 1;
  Use["multi_FDA_GA"]    = 0;
  //
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  // --- Here the preparation phase begins

  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  TString Label;
  if (buildCommittee) Label = label + "_C";
  else Label = label;
  TString outfileName = Label + ".root";

  TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is 
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( Label.Data(), outputFile,
                                              "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );


  if (mvaVar["TkNchi2"])              factory->AddVariable( "TkNchi2",              "TkNchi2",              "", 'F' );
  if (mvaVar["GlobalNchi2"])          factory->AddVariable( "GlobalNchi2",          "GlobalNchi2",          "", 'F' );
  if (mvaVar["NValidHits"])           factory->AddVariable( "NValidHits",           "NValidHits",           "", 'F' );
  if (mvaVar["NTrackerHits"])         factory->AddVariable( "NTrackerHits",         "NTrackerHits",         "", 'F' );
  if (mvaVar["NPixelHits"])           factory->AddVariable( "NPixelHits",           "NPixelHits",           "", 'F' );
  if (mvaVar["NMatches"])             factory->AddVariable( "NMatches",             "NMatches",             "", 'F' );
  if (mvaVar["D0"])                   factory->AddVariable( "D0",                   "D0",                   "", 'F' );
  if (mvaVar["IP3d"])                 factory->AddVariable( "IP3d",                 "IP3d",                 "", 'F' );
  if (mvaVar["IP3dSig"])              factory->AddVariable( "IP3dSig",              "IP3dSig",              "", 'F' );
  if (mvaVar["TrkKink"])              factory->AddVariable( "TrkKink",              "TrkKink",              "", 'F' );
  if (mvaVar["GlobalKink"])           factory->AddVariable( "GlobalKink",           "GlobalKink",           "", 'I' );
  if (mvaVar["SegmentCompatibility"]) factory->AddVariable( "SegmentCompatibility", "SegmentCompatibility", "", 'F' );
  if (mvaVar["CaloCompatibility"])    factory->AddVariable( "CaloCompatibility",    "CaloCompatibility",    "", 'F' );
  if (mvaVar["HadEnergyOverPt"])      factory->AddVariable( "HadEnergyOverPt",      "HadEnergyOverPt",      "", 'F' );
  if (mvaVar["HoEnergyOverPt"])       factory->AddVariable( "HoEnergyOverPt",       "HoEnergyOverPt",       "", 'F' );
  if (mvaVar["EmEnergyOverPt"])       factory->AddVariable( "EmEnergyOverPt",       "EmEnergyOverPt",       "", 'F' );
  if (mvaVar["HadS9EnergyOverPt"])    factory->AddVariable( "HadS9EnergyOverPt",    "HadS9EnergyOverPt",    "", 'F' );
  if (mvaVar["HoS9EnergyOverPt"])     factory->AddVariable( "HoS9EnergyOverPt",     "HoS9EnergyOverPt",     "", 'F' );
  if (mvaVar["EmS9EnergyOverPt"])     factory->AddVariable( "EmS9EnergyOverPt",     "EmS9EnergyOverPt",     "", 'F' );
  if (mvaVar["ChargedIso03OverPt"])   factory->AddVariable( "ChargedIso03OverPt",   "ChargedIso03OverPt",   "", 'F' );
  if (mvaVar["NeutralIso03OverPt"])   factory->AddVariable( "NeutralIso03OverPt",   "NeutralIso03OverPt",   "", 'F' );
  if (mvaVar["ChargedIso04OverPt"])   factory->AddVariable( "ChargedIso04OverPt",   "ChargedIso04OverPt",   "", 'F' );
  if (mvaVar["NeutralIso04OverPt"])   factory->AddVariable( "NeutralIso04OverPt",   "NeutralIso04OverPt",   "", 'F' );
  if (mvaVar["TrkIso03OverPt"])       factory->AddVariable( "TrkIso03OverPt",       "TrkIso03OverPt",       "", 'F' );
  if (mvaVar["EMIso03OverPt"])        factory->AddVariable( "EMIso03OverPt",        "EMIso03OverPt",        "", 'F' );
  if (mvaVar["HadIso03OverPt"])       factory->AddVariable( "HadIso03OverPt",       "HadIso03OverPt",       "", 'F' );
  if (mvaVar["TrkIso05OverPt"])       factory->AddVariable( "TrkIso05OverPt",       "TrkIso05OverPt",       "", 'F' );
  if (mvaVar["EMIso05OverPt"])        factory->AddVariable( "EMIso05OverPt",        "EMIso05OverPt",        "", 'F' );
  if (mvaVar["HadIso05OverPt"])       factory->AddVariable( "HadIso05OverPt",       "HadIso05OverPt",       "", 'F' );
  if (mvaVar["NVtx"])                 factory->AddVariable( "NVtx",                 "NVtx",                 "", 'F' );

  int nVariablesTemp = 0;

  if (mvaVar["TkNchi2"])              { cout << "Adding variable to MVA training: TkNchi2"              << endl; nVariablesTemp++; }
  if (mvaVar["GlobalNchi2"])          { cout << "Adding variable to MVA training: GlobalNchi2"          << endl; nVariablesTemp++; }
  if (mvaVar["NValidHits"])           { cout << "Adding variable to MVA training: NValidHits"           << endl; nVariablesTemp++; }
  if (mvaVar["NTrackerHits"])         { cout << "Adding variable to MVA training: NTrackerHits"         << endl; nVariablesTemp++; }
  if (mvaVar["NPixelHits"])           { cout << "Adding variable to MVA training: NPixelHits"           << endl; nVariablesTemp++; }
  if (mvaVar["NMatches"])             { cout << "Adding variable to MVA training: NMatches"             << endl; nVariablesTemp++; }
  if (mvaVar["D0"])                   { cout << "Adding variable to MVA training: D0"                   << endl; nVariablesTemp++; }
  if (mvaVar["IP3d"])                 { cout << "Adding variable to MVA training: IP3d"                 << endl; nVariablesTemp++; }
  if (mvaVar["IP3dSig"])              { cout << "Adding variable to MVA training: IP3dSig"              << endl; nVariablesTemp++; }
  if (mvaVar["TrkKink"])              { cout << "Adding variable to MVA training: TrkKink"              << endl; nVariablesTemp++; }
  if (mvaVar["GlobalKink"])           { cout << "Adding variable to MVA training: GlobalKink"           << endl; nVariablesTemp++; }
  if (mvaVar["SegmentCompatibility"]) { cout << "Adding variable to MVA training: SegmentCompatibility" << endl; nVariablesTemp++; }
  if (mvaVar["CaloCompatibility"])    { cout << "Adding variable to MVA training: CaloCompatibility"    << endl; nVariablesTemp++; }
  if (mvaVar["HadEnergyOverPt"])      { cout << "Adding variable to MVA training: HadEnergyOverPt"      << endl; nVariablesTemp++; }
  if (mvaVar["HoEnergyOverPt"])       { cout << "Adding variable to MVA training: HoEnergyOverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["EmEnergyOverPt"])       { cout << "Adding variable to MVA training: EmEnergyOverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["HadS9EnergyOverPt"])    { cout << "Adding variable to MVA training: HadS9EnergyOverPt"    << endl; nVariablesTemp++; }
  if (mvaVar["HoS9EnergyOverPt"])     { cout << "Adding variable to MVA training: HoS9EnergyOverPt"     << endl; nVariablesTemp++; }
  if (mvaVar["EmS9EnergyOverPt"])     { cout << "Adding variable to MVA training: EmS9EnergyOverPt"     << endl; nVariablesTemp++; }
  if (mvaVar["ChargedIso03OverPt"])   { cout << "Adding variable to MVA training: ChargedIso03OverPt"   << endl; nVariablesTemp++; }
  if (mvaVar["NeutralIso03OverPt"])   { cout << "Adding variable to MVA training: NeutralIso03OverPt"   << endl; nVariablesTemp++; }
  if (mvaVar["ChargedIso04OverPt"])   { cout << "Adding variable to MVA training: ChargedIso04OverPt"   << endl; nVariablesTemp++; }
  if (mvaVar["NeutralIso04OverPt"])   { cout << "Adding variable to MVA training: NeutralIso04OverPt"   << endl; nVariablesTemp++; }
  if (mvaVar["TrkIso03OverPt"])       { cout << "Adding variable to MVA training: TrkIso03OverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["EMIso03OverPt"])        { cout << "Adding variable to MVA training: EMIso03OverPt"        << endl; nVariablesTemp++; }
  if (mvaVar["HadIso03OverPt"])       { cout << "Adding variable to MVA training: HadIso03OverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["TrkIso05OverPt"])       { cout << "Adding variable to MVA training: TrkIso05OverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["EMIso05OverPt"])        { cout << "Adding variable to MVA training: EMIso05OverPt"        << endl; nVariablesTemp++; }
  if (mvaVar["HadIso05OverPt"])       { cout << "Adding variable to MVA training: HadIso05OverPt"       << endl; nVariablesTemp++; }
  if (mvaVar["NVtx"])                 { cout << "Adding variable to MVA training: NVtx"                 << endl; nVariablesTemp++; }
                                    
  const unsigned int nVariables = nVariablesTemp;
  cout << "Using " << nVariables << " variables for MVA training" << endl;

  // You can add so-called "Spectator variables", which are not used in the MVA training,
  // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
  // input variables, the response values of all trained MVAs, and the spectator variables
  //factory->AddSpectator( "njets",  "Event Type", "units", 'F' );
  //factory->AddSpectator( "mydilmass := ",  "Spectator 2", "units", 'F' );

  TTree *signal     = (TTree*) chsignal;
  TTree *background = (TTree*) chbackground;
   
  std::cout << "--- TMVAClassification       : Using bkg input files: -------------------" <<  std::endl;

  TObjArray *listOfBkgFiles = chbackground->GetListOfFiles();
  TIter bkgFileIter(listOfBkgFiles);
  TChainElement* currentBkgFile = 0;

  while((currentBkgFile = (TChainElement*)bkgFileIter.Next())) {
    std::cout << currentBkgFile->GetTitle() << std::endl;
  }

  std::cout << "--- TMVAClassification       : Using sig input files: -------------------" <<  std::endl;
   
  TObjArray *listOfSigFiles = chsignal->GetListOfFiles();
  TIter sigFileIter(listOfSigFiles);
  TChainElement* currentSigFile = 0;

  while((currentSigFile = (TChainElement*)sigFileIter.Next())) {
    std::cout << currentSigFile->GetTitle() << std::endl;
  }
                                              
  // global event weights per tree (see below for setting event-wise weights)
  //Double_t signalWeight     = 1.0;
  //Double_t backgroundWeight = 1.0;
   
  // You can add an arbitrary number of signal or background trees
  //factory->AddSignalTree    ( signal,     signalWeight     );
  //factory->AddBackgroundTree( background, backgroundWeight );
      
  // To give different trees for training and testing, do as follows:
  //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
  //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

  // Use the following code instead of the above two or four lines to add signal and background
  // training and test events "by hand"
  // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
  //      variable definition, but simply compute the expression before adding the event
  //
  
  std::vector<Double_t> vars( nVariables );

  //Variables
  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fMuPt; 
  Float_t                 fMuEta; 
  Float_t                 fMuPhi; 
  Float_t                 fMuPFIso; 
  
  //CutBased Variables
  Float_t                 fMuTkNchi2; 
  Float_t                 fMuGlobalNchi2; 
  Float_t                 fMuNValidHits; 
  Float_t                 fMuNTrackerHits; 
  Float_t                 fMuNPixelHits; 
  Float_t                 fMuNMatches; 
  Float_t                 fMuD0; 

  //Additional Vars used in Likelihood
  Float_t                 fMuIP3d; 
  Float_t                 fMuIP3dSig; 
  Float_t                 fMuTrkKink; 
  Float_t                 fMuGlobalKink; 
  Float_t                 fMuSegmentCompatibility; 
  Float_t                 fMuCaloCompatibility; 
  Float_t                 fMuHadEnergy; 
  Float_t                 fMuHoEnergy; 
  Float_t                 fMuEmEnergy; 
  Float_t                 fMuHadS9Energy; 
  Float_t                 fMuHoS9Energy; 
  Float_t                 fMuEmS9Energy; 

  //Isolation Variables
  Float_t                 fMuChargedIso03; 
  Float_t                 fMuChargedIso03FromOtherVertices; 
  Float_t                 fMuNeutralIso03_05Threshold; 
  Float_t                 fMuNeutralIso03_10Threshold; 
  Float_t                 fMuChargedIso04; 
  Float_t                 fMuChargedIso04FromOtherVertices; 
  Float_t                 fMuNeutralIso04_05Threshold; 
  Float_t                 fMuNeutralIso04_10Threshold; 
  Float_t                 fMuTrkIso03; 
  Float_t                 fMuEMIso03; 
  Float_t                 fMuHadIso03; 
  Float_t                 fMuTrkIso05; 
  Float_t                 fMuEMIso05; 
  Float_t                 fMuHadIso05; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 


  signal->SetBranchAddress( "weight", &fWeight);
  signal->SetBranchAddress( "run", &fRunNumber);
  signal->SetBranchAddress( "lumi", &fLumiSectionNumber);
  signal->SetBranchAddress( "event", &fEventNumber);
  signal->SetBranchAddress( "pt", &fMuPt); 
  signal->SetBranchAddress( "eta", &fMuEta); 
  signal->SetBranchAddress( "phi", &fMuPhi); 
  signal->SetBranchAddress( "pfiso", &fMuPFIso); 
  signal->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  signal->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  signal->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  signal->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  signal->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  signal->SetBranchAddress( "NMatches", &fMuNMatches); 
  signal->SetBranchAddress( "D0", &fMuD0); 
  signal->SetBranchAddress( "IP3d", &fMuIP3d); 
  signal->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  signal->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  signal->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  signal->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  signal->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  signal->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  signal->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  signal->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  signal->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  signal->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  signal->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  signal->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  signal->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  signal->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  signal->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  signal->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  signal->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  signal->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  signal->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  signal->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  signal->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  signal->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  signal->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  signal->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  signal->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  signal->SetBranchAddress( "Rho", &fRho); 
  signal->SetBranchAddress( "NVertices", &fNVertices); 

  int nsigtrain = 0;
  int nsigtest  = 0;
  int nbkgtrain = 0;
  int nbkgtest  = 0;

  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);

    //--------------------------------------------------
    // SIGNAL EVENT SELECTION
    //--------------------------------------------------

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    if( fMuPt < 10.0                     ) continue; 
    if( fabs(fMuEta) > 2.4               ) continue; 
    if( fMuPt > 35.0                     ) continue; // too high pt muon are contaminated with signal
    if( (fMuChargedIso03 + fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03,fMuEta))/fMuPt > 0.4 ) continue;

    int varCounter = 0;

    if (mvaVar["TkNchi2"])              vars[varCounter++] = fMuTkNchi2;
    if (mvaVar["GlobalNchi2"])          vars[varCounter++] = fMuGlobalNchi2;
    if (mvaVar["NValidHits"])           vars[varCounter++] = fMuNValidHits;
    if (mvaVar["NTrackerHits"])         vars[varCounter++] = fMuNTrackerHits;
    if (mvaVar["NPixelHits"])           vars[varCounter++] = fMuNPixelHits;
    if (mvaVar["NMatches"])             vars[varCounter++] = fMuNMatches;
    if (mvaVar["D0"])                   vars[varCounter++] = fMuD0;
    if (mvaVar["IP3d"])                 vars[varCounter++] = fMuIP3d;
    if (mvaVar["IP3dSig"])              vars[varCounter++] = fMuIP3dSig;
    if (mvaVar["TrkKink"])              vars[varCounter++] = fMuTrkKink;
    if (mvaVar["GlobalKink"])           vars[varCounter++] = fMuGlobalKink;
    if (mvaVar["SegmentCompatibility"]) vars[varCounter++] = fMuSegmentCompatibility;
    if (mvaVar["CaloCompatibility"])    vars[varCounter++] = fMuCaloCompatibility;
    if (mvaVar["HadEnergyOverPt"])      vars[varCounter++] = (fMuHadEnergy - rho*MuonEffectiveArea(kMuHadEnergy,fMuEta))/fMuPt;
    if (mvaVar["HoEnergyOverPt"])       vars[varCounter++] = (fMuHoEnergy - rho*MuonEffectiveArea(kMuHoEnergy,fMuEta))/fMuPt;
    if (mvaVar["EmEnergyOverPt"])       vars[varCounter++] = (fMuEmEnergy - rho*MuonEffectiveArea(kMuEmEnergy,fMuEta))/fMuPt;
    if (mvaVar["HadS9EnergyOverPt"])    vars[varCounter++] = (fMuHadS9Energy - rho*MuonEffectiveArea(kMuHadS9Energy,fMuEta))/fMuPt;
    if (mvaVar["HoS9EnergyOverPt"])     vars[varCounter++] = (fMuHoS9Energy - rho*MuonEffectiveArea(kMuHoS9Energy,fMuEta))/fMuPt;
    if (mvaVar["EmS9EnergyOverPt"])     vars[varCounter++] = (fMuEmS9Energy - rho*MuonEffectiveArea(kMuEmS9Energy,fMuEta))/fMuPt;
    if (mvaVar["ChargedIso03OverPt"])   vars[varCounter++] = (fMuChargedIso03 - rho*MuonEffectiveArea(kMuChargedIso03,fMuEta))/fMuPt;
    if (mvaVar["NeutralIso03OverPt"])   vars[varCounter++] = (fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03,fMuEta))/fMuPt;
    if (mvaVar["ChargedIso04OverPt"])   vars[varCounter++] = (fMuChargedIso04 - rho*MuonEffectiveArea(kMuChargedIso04,fMuEta))/fMuPt;
    if (mvaVar["NeutralIso04OverPt"])   vars[varCounter++] = (fMuNeutralIso04_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso04,fMuEta))/fMuPt;
    if (versionLabel == "V11") {
      if (mvaVar["TrkIso03OverPt"])       vars[varCounter++] = (fMuTrkIso03)/fMuPt;
      if (mvaVar["EMIso03OverPt"])        vars[varCounter++] = (fMuEMIso03)/fMuPt;
      if (mvaVar["HadIso03OverPt"])       vars[varCounter++] = (fMuHadIso03)/fMuPt;
      if (mvaVar["TrkIso05OverPt"])       vars[varCounter++] = (fMuTrkIso05)/fMuPt;
      if (mvaVar["EMIso05OverPt"])        vars[varCounter++] = (fMuEMIso05)/fMuPt;
      if (mvaVar["HadIso05OverPt"])       vars[varCounter++] = (fMuHadIso05)/fMuPt;
    } else {
      if (mvaVar["TrkIso03OverPt"])       vars[varCounter++] = (fMuTrkIso03 - rho*MuonEffectiveArea(kMuTrkIso03,fMuEta))/fMuPt;
      if (mvaVar["EMIso03OverPt"])        vars[varCounter++] = (fMuEMIso03 - rho*MuonEffectiveArea(kMuEMIso03,fMuEta))/fMuPt;
      if (mvaVar["HadIso03OverPt"])       vars[varCounter++] = (fMuHadIso03 - rho*MuonEffectiveArea(kMuHadIso03,fMuEta))/fMuPt;
      if (mvaVar["TrkIso05OverPt"])       vars[varCounter++] = (fMuTrkIso05 - rho*MuonEffectiveArea(kMuTrkIso05,fMuEta))/fMuPt;
      if (mvaVar["EMIso05OverPt"])        vars[varCounter++] = (fMuEMIso05 - rho*MuonEffectiveArea(kMuEMIso05,fMuEta))/fMuPt;
      if (mvaVar["HadIso05OverPt"])       vars[varCounter++] = (fMuHadIso05 - rho*MuonEffectiveArea(kMuHadIso05,fMuEta))/fMuPt;
    }
    if (mvaVar["NVtx"])                 vars[varCounter++] = fNVertices;

    //This bin has too many events in the training sample, causes program to crash
    if (Option == 4) {
      if ( fEventNumber % 8 == 0 ) {
        factory->AddSignalTrainingEvent( vars, fWeight );
        nsigtrain++;
      }
      else if ( fEventNumber % 8 == 1 ) {
        factory->AddSignalTestEvent    ( vars, fWeight );
        nsigtest++;
      }
    } else {
      if ( fEventNumber % 2 == 0 ) {
        factory->AddSignalTrainingEvent( vars, fWeight );
        nsigtrain++;
      }
      else if ( fEventNumber % 2 == 1 ) {
        factory->AddSignalTestEvent    ( vars, fWeight );
        nsigtest++;
      }
    }
  }

  background->SetBranchAddress( "weight", &fWeight);
  background->SetBranchAddress( "run", &fRunNumber);
  background->SetBranchAddress( "lumi", &fLumiSectionNumber);
  background->SetBranchAddress( "event", &fEventNumber);
  background->SetBranchAddress( "pt", &fMuPt); 
  background->SetBranchAddress( "eta", &fMuEta); 
  background->SetBranchAddress( "phi", &fMuPhi); 
  background->SetBranchAddress( "pfiso", &fMuPFIso); 
  background->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
  background->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
  background->SetBranchAddress( "NValidHits", &fMuNValidHits); 
  background->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
  background->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
  background->SetBranchAddress( "NMatches", &fMuNMatches); 
  background->SetBranchAddress( "D0", &fMuD0); 
  background->SetBranchAddress( "IP3d", &fMuIP3d); 
  background->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
  background->SetBranchAddress( "TrkKink", &fMuTrkKink); 
  background->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
  background->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
  background->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
  background->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
  background->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
  background->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
  background->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
  background->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
  background->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
  background->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
  background->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
  background->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
  background->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
  background->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
  background->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
  background->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
  background->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
  background->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
  background->SetBranchAddress( "EMIso03", &fMuEMIso03); 
  background->SetBranchAddress( "HadIso03", &fMuHadIso03); 
  background->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
  background->SetBranchAddress( "EMIso05", &fMuEMIso05); 
  background->SetBranchAddress( "HadIso05", &fMuHadIso05); 
  background->SetBranchAddress( "Rho", &fRho); 
  background->SetBranchAddress( "NVertices", &fNVertices); 

  cout << "Add background events" << endl;
  cout << "Added " << nsigtrain << " training events" << endl;
  cout << "Added " << nsigtest  << " test events" << endl;
  
  for (UInt_t i=0; i<background->GetEntries(); i++) {
    
    background->GetEntry(i);

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

    //--------------------------------------------------
    // BACKGROUND EVENT SELECTION
    //--------------------------------------------------
    if( fMuPt < 10.0                     ) continue; 
    if( fabs(fMuEta) > 2.4               ) continue; 
    if( fMuPt > 35.0                     ) continue; // too high pt muon are contaminated with signal
    if( (fMuChargedIso03 + fMuNeutralIso03_05Threshold - fRho*MuonEffectiveArea(kMuNeutralIso03,fMuEta))/fMuPt > 0.4 ) continue;

    int varCounter = 0;

    if (mvaVar["TkNchi2"])              vars[varCounter++] = fMuTkNchi2;
    if (mvaVar["GlobalNchi2"])          vars[varCounter++] = fMuGlobalNchi2;
    if (mvaVar["NValidHits"])           vars[varCounter++] = fMuNValidHits;
    if (mvaVar["NTrackerHits"])         vars[varCounter++] = fMuNTrackerHits;
    if (mvaVar["NPixelHits"])           vars[varCounter++] = fMuNPixelHits;
    if (mvaVar["NMatches"])             vars[varCounter++] = fMuNMatches;
    if (mvaVar["D0"])                   vars[varCounter++] = fMuD0;
    if (mvaVar["IP3d"])                 vars[varCounter++] = fMuIP3d;
    if (mvaVar["IP3dSig"])              vars[varCounter++] = fMuIP3dSig;
    if (mvaVar["TrkKink"])              vars[varCounter++] = fMuTrkKink;
    if (mvaVar["GlobalKink"])           vars[varCounter++] = fMuGlobalKink;
    if (mvaVar["SegmentCompatibility"]) vars[varCounter++] = fMuSegmentCompatibility;
    if (mvaVar["CaloCompatibility"])    vars[varCounter++] = fMuCaloCompatibility;
    if (mvaVar["HadEnergyOverPt"])      vars[varCounter++] = (fMuHadEnergy - rho*MuonEffectiveArea(kMuHadEnergy,fMuEta))/fMuPt;
    if (mvaVar["HoEnergyOverPt"])       vars[varCounter++] = (fMuHoEnergy - rho*MuonEffectiveArea(kMuHoEnergy,fMuEta))/fMuPt;
    if (mvaVar["EmEnergyOverPt"])       vars[varCounter++] = (fMuEmEnergy - rho*MuonEffectiveArea(kMuEmEnergy,fMuEta))/fMuPt;
    if (mvaVar["HadS9EnergyOverPt"])    vars[varCounter++] = (fMuHadS9Energy - rho*MuonEffectiveArea(kMuHadS9Energy,fMuEta))/fMuPt;
    if (mvaVar["HoS9EnergyOverPt"])     vars[varCounter++] = (fMuHoS9Energy - rho*MuonEffectiveArea(kMuHoS9Energy,fMuEta))/fMuPt;
    if (mvaVar["EmS9EnergyOverPt"])     vars[varCounter++] = (fMuEmS9Energy - rho*MuonEffectiveArea(kMuEmS9Energy,fMuEta))/fMuPt;
    if (mvaVar["ChargedIso03OverPt"])   vars[varCounter++] = (fMuChargedIso03 - rho*MuonEffectiveArea(kMuChargedIso03,fMuEta))/fMuPt;
    if (mvaVar["NeutralIso03OverPt"])   vars[varCounter++] = (fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03,fMuEta))/fMuPt;
    if (mvaVar["ChargedIso04OverPt"])   vars[varCounter++] = (fMuChargedIso04 - rho*MuonEffectiveArea(kMuChargedIso04,fMuEta))/fMuPt;
    if (mvaVar["NeutralIso04OverPt"])   vars[varCounter++] = (fMuNeutralIso04_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso04,fMuEta))/fMuPt;
    if (versionLabel == "V11") {
      if (mvaVar["TrkIso03OverPt"])       vars[varCounter++] = (fMuTrkIso03)/fMuPt;
      if (mvaVar["EMIso03OverPt"])        vars[varCounter++] = (fMuEMIso03)/fMuPt;
      if (mvaVar["HadIso03OverPt"])       vars[varCounter++] = (fMuHadIso03)/fMuPt;
      if (mvaVar["TrkIso05OverPt"])       vars[varCounter++] = (fMuTrkIso05)/fMuPt;
      if (mvaVar["EMIso05OverPt"])        vars[varCounter++] = (fMuEMIso05)/fMuPt;
      if (mvaVar["HadIso05OverPt"])       vars[varCounter++] = (fMuHadIso05)/fMuPt;
    } else {
      if (mvaVar["TrkIso03OverPt"])       vars[varCounter++] = (fMuTrkIso03 - rho*MuonEffectiveArea(kMuTrkIso03,fMuEta))/fMuPt;
      if (mvaVar["EMIso03OverPt"])        vars[varCounter++] = (fMuEMIso03 - rho*MuonEffectiveArea(kMuEMIso03,fMuEta))/fMuPt;
      if (mvaVar["HadIso03OverPt"])       vars[varCounter++] = (fMuHadIso03 - rho*MuonEffectiveArea(kMuHadIso03,fMuEta))/fMuPt;
      if (mvaVar["TrkIso05OverPt"])       vars[varCounter++] = (fMuTrkIso05 - rho*MuonEffectiveArea(kMuTrkIso05,fMuEta))/fMuPt;
      if (mvaVar["EMIso05OverPt"])        vars[varCounter++] = (fMuEMIso05 - rho*MuonEffectiveArea(kMuEMIso05,fMuEta))/fMuPt;
      if (mvaVar["HadIso05OverPt"])       vars[varCounter++] = (fMuHadIso05 - rho*MuonEffectiveArea(kMuHadIso05,fMuEta))/fMuPt;
    }
    if (mvaVar["NVtx"])                 vars[varCounter++] = fNVertices;
 
    if ( fEventNumber % 2 == 0 ){
      factory->AddBackgroundTrainingEvent( vars, fWeight );
      nbkgtrain++;
    }
    else{
      factory->AddBackgroundTestEvent    ( vars, fWeight );
      nbkgtest++;
    }
  }
  
  cout << "Done adding background" << endl;
  cout << "Added " << nbkgtrain << " training events" << endl;
  cout << "Added " << nbkgtest  << " test events" << endl;
  

  // --- end ------------------------------------------------------------
  //
  // --- end of tree registration 
   
  // Set individual event weights (the variables must exist in the original TTree)
 //  factory->SetSignalWeightExpression    ("weight");
 //  factory->SetBackgroundWeightExpression("weight");
  cout << "Done setting weights" << endl;

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  // Tell the factory how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used 
  // for training, and the other half for testing:
  //    factory->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  // To also specify the number of testing events, use:
  //    factory->PrepareTrainingAndTestTree( mycut,
  //                                         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );
   
  //Use random splitting
  factory->PrepareTrainingAndTestTree( mycuts, mycutb,
                                       "nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

  //Use alternate splitting 
  //(this is preferable since its easier to track which events were used for training, but the job crashes! need to fix this...)
  //factory->PrepareTrainingAndTestTree( mycuts, mycutb,
  //                                     "nTrain_Signal=0:nTrain_Background=0:SplitMode=Alternate:NormMode=NumEvents:!V" );

  // ---- Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( TMVA::Types::kCuts, "Cuts",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

  if (Use["CutsD"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsD",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

  if (Use["CutsPCA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsPCA",
                         "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

  if (Use["CutsGA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsGA",
                         "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  if (Use["CutsSA"])
    factory->BookMethod( TMVA::Types::kCuts, "CutsSA",
                         "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
//     factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
//                          "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=100:NSmoothBkg[0]=100:NSmoothBkg[1]=10:NSmooth=0:Nbins=200" );
    factory->BookMethod( TMVA::Types::kLikelihood, "Likelihood",
                         "H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=30:NSmoothBkg[0]=30:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=40:Nbins=200" );

  // Decorrelated likelihood
  if (Use["LikelihoodD"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodD",
                         "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=40:VarTransform=Decorrelate" );

  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodPCA",
                         "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" ); 

  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodKDE",
                         "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" ); 

  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"])
    factory->BookMethod( TMVA::Types::kLikelihood, "LikelihoodMIX",
                         "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" ); 

  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERS",
                         "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

  if (Use["PDERSD"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERSD",
                         "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

  if (Use["PDERSPCA"])
    factory->BookMethod( TMVA::Types::kPDERS, "PDERSPCA",
                         "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"])
    factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoam",
                         "H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0333:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

  if (Use["PDEFoamBoost"])
    factory->BookMethod( TMVA::Types::kPDEFoam, "PDEFoamBoost",
                         "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod( TMVA::Types::kKNN, "KNN",
                         "H:CreateMVAPdfs:nkNN=31:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"])
    factory->BookMethod( TMVA::Types::kHMatrix, "HMatrix", "!H:!V" );

  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"])
    factory->BookMethod( TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=40:NsmoothMVAPdf=10" );

  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"])
    factory->BookMethod( TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"])
    factory->BookMethod( TMVA::Types::kFisher, "BoostedFisher", 
                         "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MC",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GA",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );

  if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_SA",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  if (Use["FDA_MT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_GAMT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  if (Use["FDA_MCMT"])
    factory->BookMethod( TMVA::Types::kFDA, "FDA_MCMT",
                         "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( TMVA::Types::kMLP, "MLP", "H:!V:CreateMVAPdfs:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  if (Use["MLPBFGS"])
    factory->BookMethod( TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

  if (Use["MLPBNN"])
    factory->BookMethod( TMVA::Types::kMLP, "MLPBNN", "H:!V:CreateMVAPdfs:NeuronType=tanh:VarTransform=N:NCycles=500:HiddenLayers=N+3:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

  // CF(Clermont-Ferrand)ANN
  if (Use["CFMlpANN"])
    factory->BookMethod( TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...  

  // Tmlp(Root)ANN
  if (Use["TMlpANN"])
    factory->BookMethod( TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod( TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
//     factory->BookMethod( TMVA::Types::kBDT, "BDTG",
//                          "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=2000:NNodesMax=5" );


//For NEW Version
//Difference seems to be UsePoissonNvars == True, which is not implemented in old version
//     if (versionLabel == "V2") {
//       //turn off UsePoissonNvars (V2)
//       factory->BookMethod( TMVA::Types::kBDT, "BDTG",
//                            "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:UseBaggedGrad:GradBaggingFraction=0.5:nCuts=2000:NNodesMax=5:!UsePoissonNvars" );
//     }
//   //JOSH Suggestions (V4)
//     if (versionLabel == "V4") {
//       factory->BookMethod( TMVA::Types::kBDT, "BDTG",
//                            "!H:!V:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:MaxDepth=6" );
//     }

//For TMVA 4.0.7 VErsion
    // Don't use bagging
//     if (versionLabel == "V0") {
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:CreateMVAPdfs:NTrees=1000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning" );
//     }
  

  if (Use["BDT"])  // Adaptive Boost
    factory->BookMethod(TMVA::Types::kBDT,"BDT",
			"!H:!V:CreateMVAPdfs:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=2000:PruneMethod=CostComplexity:PruneStrength=20.0");

  if (Use["BDTB"]) // Bagging
    factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                         "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=2000:PruneMethod=NoPruning" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod(TMVA::Types::kBDT,"BDTD",
                    "!H:!V:NTrees=400:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=2000:PruneMethod=CostComplexity:PruneStrength=25.0:VarTransform=Decorrelate");

  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( TMVA::Types::kRuleFit, "RuleFit",
                         "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

   
  // For an example of the category classifier usage, see: TMVAClassificationCategory

  // --------------------------------------------------------------------------------------------------

  // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

  // factory->OptimizeAllMethods("SigEffAt001","Scan");
  // factory->OptimizeAllMethods("ROCIntegral","GA");

  // --------------------------------------------------------------------------------------------------

  // ---- Now you can tell the factory to train, test, and evaluate the MVAs
  
  // Train MVAs using the set of training events
  factory->TrainAllMethods();
  
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  
  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  
  delete factory;

  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVAGui( outfileName );
}
