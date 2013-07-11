// @(#)root/tmva $Id: TrainElectronMVA.C,v 1.5 2012/02/27 09:26:38 sixie Exp $
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
#include "TMVA/Config.h"
#endif

typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

void TrainElectronMVA(
 TString sigInputFile    = "ElectronSelectionTraining.Real.weighted.root",
 TString bgdInputFile    = "ElectronSelectionTraining.Fake.weighted.root",
 TString label           = "default",
 string  versionLabel    = "",
 Bool_t  buildCommittee  = kFALSE,
 Int_t Option            = -1,
//  TString myMethodList    = "Likelihood,BDT,BDTG"
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
  mvaVar[ "SigmaIEtaIEta" ]         = 0;  
  mvaVar[ "DEtaIn" ]                = 0; 
  mvaVar[ "DPhiIn" ]                = 0;  
  mvaVar[ "D0" ]                    = 0;  
  mvaVar[ "DZ" ]                    = 0; 
  mvaVar[ "FBrem" ]                 = 0; 
  mvaVar[ "EOverP" ]                = 0;
  mvaVar[ "ESeedClusterOverPout" ]  = 0;
  mvaVar[ "SigmaIPhiIPhi" ]         = 0;  
  mvaVar[ "NBrem" ]                 = 0;  
  mvaVar[ "OneOverEMinusOneOverP" ] = 0;  
  mvaVar[ "ESeedClusterOverPIn" ]   = 0;  
  mvaVar[ "IP3d" ]                  = 0;  
  mvaVar[ "IP3dSig" ]               = 0;  
  mvaVar[ "StandardLikelihood" ]    = 0;  
  mvaVar[ "PFMVA" ]                 = 0;  
  mvaVar[ "TMVAKNN" ]               = 0;  
  mvaVar[ "TMVAMLP" ]               = 0;  
  mvaVar[ "TMVAMLPBNN" ]            = 0;  
  mvaVar[ "TMVABDTG" ]              = 0;  
  mvaVar[ "TMVABDT" ]               = 0;  
  mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
  mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
  mvaVar[ "dPhiCalo" ]              = 0;
  mvaVar[ "R9" ]                    = 0; //ShowerShape
  mvaVar[ "SCEtaWidth" ]            = 0;
  mvaVar[ "SCPhiWidth" ]            = 0;
  mvaVar[ "CovIEtaIPhi" ]           = 0;
  mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
  mvaVar[ "HoverE" ]                = 0; //HOverE
  mvaVar[ "HcalDepth1OverEcal" ]    = 0;
  mvaVar[ "HcalDepth2OverEcal" ]    = 0;
  mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
  mvaVar[ "SeedETopOverE" ]         = 0;
  mvaVar[ "SeedEBottomOverE" ]      = 0;
  mvaVar[ "SeedELeftOverE" ]        = 0;
  mvaVar[ "SeedERightOverE" ]       = 0;
  mvaVar[ "SeedE2ndOverE" ]         = 0;
  mvaVar[ "SeedE2x5RightOverE" ]    = 0;
  mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
  mvaVar[ "SeedE2x5TopOverE" ]      = 0;
  mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
  mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
  mvaVar[ "SeedE1x3OverE" ]         = 0;
  mvaVar[ "SeedE3x1OverE" ]         = 0;
  mvaVar[ "SeedE1x5OverE" ]         = 0;
  mvaVar[ "SeedE2x2OverE" ]         = 0;
  mvaVar[ "SeedE3x2OverE" ]         = 0;
  mvaVar[ "SeedE3x3OverE" ]         = 0;
  mvaVar[ "SeedE4x4OverE" ]         = 0;
  mvaVar[ "SeedE5x5OverE" ]         = 0;
  mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
  mvaVar[ "NeutralHadronIso03" ]    = 0; 
  mvaVar[ "GammaIso03" ]            = 0; 
  mvaVar[ "ChargedIso04" ]          = 0; 
  mvaVar[ "NeutralHadronIso04" ]    = 0; 
  mvaVar[ "GammaIso04" ]            = 0; 
  mvaVar[ "TrkIso03" ]              = 0; 
  mvaVar[ "EMIso03" ]               = 0; 
  mvaVar[ "HadIso03" ]              = 0; 
  mvaVar[ "TrkIso04" ]              = 0; 
  mvaVar[ "EMIso04" ]               = 0; 
  mvaVar[ "HadIso04" ]              = 0; 


  if (versionLabel == "V0") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }


  if (versionLabel == "V1") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }



  if (versionLabel == "V2") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V3") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V4") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 1; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 1;
    mvaVar[ "HcalDepth2OverEcal" ]    = 1;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V5") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 1; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 1;
    mvaVar[ "SeedEBottomOverE" ]      = 1;
    mvaVar[ "SeedELeftOverE" ]        = 1;
    mvaVar[ "SeedERightOverE" ]       = 1;
    mvaVar[ "SeedE2ndOverE" ]         = 1;
    mvaVar[ "SeedE2x5RightOverE" ]    = 1;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 1;
    mvaVar[ "SeedE2x5TopOverE" ]      = 1;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 1;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 1;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 1;
    mvaVar[ "SeedE3x2OverE" ]         = 1;
    mvaVar[ "SeedE3x3OverE" ]         = 1;
    mvaVar[ "SeedE4x4OverE" ]         = 1;
    mvaVar[ "SeedE5x5OverE" ]         = 1;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V6") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 1; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 1;
    mvaVar[ "HcalDepth2OverEcal" ]    = 1;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V7") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 1; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 1;
    mvaVar[ "HcalDepth2OverEcal" ]    = 1;
    mvaVar[ "SeedEMaxOverE" ]         = 1; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 1;
    mvaVar[ "SeedEBottomOverE" ]      = 1;
    mvaVar[ "SeedELeftOverE" ]        = 1;
    mvaVar[ "SeedERightOverE" ]       = 1;
    mvaVar[ "SeedE2ndOverE" ]         = 1;
    mvaVar[ "SeedE2x5RightOverE" ]    = 1;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 1;
    mvaVar[ "SeedE2x5TopOverE" ]      = 1;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 1;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 1;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 1;
    mvaVar[ "SeedE3x2OverE" ]         = 1;
    mvaVar[ "SeedE3x3OverE" ]         = 1;
    mvaVar[ "SeedE4x4OverE" ]         = 1;
    mvaVar[ "SeedE5x5OverE" ]         = 1;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V9") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V10") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V11") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 0; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 0; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 0;
    mvaVar[ "R9" ]                    = 0; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 0;
    mvaVar[ "SCPhiWidth" ]            = 0;
    mvaVar[ "CovIEtaIPhi" ]           = 0;
    mvaVar[ "PreShowerOverRaw" ]      = 0; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 1; 
    mvaVar[ "NeutralHadronIso04" ]    = 1; 
    mvaVar[ "GammaIso04" ]            = 1; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V12") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 1; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 1;
    mvaVar[ "HcalDepth2OverEcal" ]    = 1;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 1; 
    mvaVar[ "NeutralHadronIso04" ]    = 1; 
    mvaVar[ "GammaIso04" ]            = 1; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V13") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 1; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 1;
    mvaVar[ "HcalDepth2OverEcal" ]    = 1;
    mvaVar[ "SeedEMaxOverE" ]         = 1; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 1;
    mvaVar[ "SeedEBottomOverE" ]      = 1;
    mvaVar[ "SeedELeftOverE" ]        = 1;
    mvaVar[ "SeedERightOverE" ]       = 1;
    mvaVar[ "SeedE2ndOverE" ]         = 1;
    mvaVar[ "SeedE2x5RightOverE" ]    = 1;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 1;
    mvaVar[ "SeedE2x5TopOverE" ]      = 1;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 1;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 1;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 1;
    mvaVar[ "SeedE3x2OverE" ]         = 1;
    mvaVar[ "SeedE3x3OverE" ]         = 1;
    mvaVar[ "SeedE4x4OverE" ]         = 1;
    mvaVar[ "SeedE5x5OverE" ]         = 1;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 1; 
    mvaVar[ "NeutralHadronIso04" ]    = 1; 
    mvaVar[ "GammaIso04" ]            = 1; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  if (versionLabel == "V14") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 1; 
    mvaVar[ "NeutralHadronIso04" ]    = 1; 
    mvaVar[ "GammaIso04" ]            = 1; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }

  //We currently use V15
  if (versionLabel == "V15") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 1; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 1; 
    mvaVar[ "GammaIso03" ]            = 1; 
    mvaVar[ "ChargedIso04" ]          = 1; 
    mvaVar[ "NeutralHadronIso04" ]    = 1; 
    mvaVar[ "GammaIso04" ]            = 1; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }


  if (versionLabel == "V16") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 1; 
    mvaVar[ "EMIso03" ]               = 1; 
    mvaVar[ "HadIso03" ]              = 1; 
    mvaVar[ "TrkIso04" ]              = 0; 
    mvaVar[ "EMIso04" ]               = 0; 
    mvaVar[ "HadIso04" ]              = 0; 
  }


  if (versionLabel == "V17") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 0; 
    mvaVar[ "EMIso03" ]               = 0; 
    mvaVar[ "HadIso03" ]              = 0; 
    mvaVar[ "TrkIso04" ]              = 1; 
    mvaVar[ "EMIso04" ]               = 1; 
    mvaVar[ "HadIso04" ]              = 1; 
  }

  if (versionLabel == "V18") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 0;  
    mvaVar[ "OneOverEMinusOneOverP" ] = 1;  
    mvaVar[ "ESeedClusterOverPIn" ]   = 1;  
    mvaVar[ "IP3d" ]                  = 1;  
    mvaVar[ "IP3dSig" ]               = 1;  
    mvaVar[ "StandardLikelihood" ]    = 0;  
    mvaVar[ "PFMVA" ]                 = 0;  
    mvaVar[ "TMVAKNN" ]               = 0;  
    mvaVar[ "TMVAMLP" ]               = 0;  
    mvaVar[ "TMVAMLPBNN" ]            = 0;  
    mvaVar[ "TMVABDTG" ]              = 0;  
    mvaVar[ "TMVABDT" ]               = 0;  
    mvaVar[ "GsfTrackChi2OverNdof" ]  = 1; //TrackQuality
    mvaVar[ "dEtaCalo" ]              = 1; //Track Matching
    mvaVar[ "dPhiCalo" ]              = 1;
    mvaVar[ "R9" ]                    = 1; //ShowerShape
    mvaVar[ "SCEtaWidth" ]            = 1;
    mvaVar[ "SCPhiWidth" ]            = 1;
    mvaVar[ "CovIEtaIPhi" ]           = 1;
    mvaVar[ "PreShowerOverRaw" ]      = 1; //Preshower
    mvaVar[ "HoverE" ]                = 0; //HOverE
    mvaVar[ "HcalDepth1OverEcal" ]    = 0;
    mvaVar[ "HcalDepth2OverEcal" ]    = 0;
    mvaVar[ "SeedEMaxOverE" ]         = 0; //Single Crystal Information
    mvaVar[ "SeedETopOverE" ]         = 0;
    mvaVar[ "SeedEBottomOverE" ]      = 0;
    mvaVar[ "SeedELeftOverE" ]        = 0;
    mvaVar[ "SeedERightOverE" ]       = 0;
    mvaVar[ "SeedE2ndOverE" ]         = 0;
    mvaVar[ "SeedE2x5RightOverE" ]    = 0;
    mvaVar[ "SeedE2x5LeftOverE" ]     = 0;
    mvaVar[ "SeedE2x5TopOverE" ]      = 0;
    mvaVar[ "SeedE2x5BottomOverE" ]   = 0;
    mvaVar[ "SeedE2x5MaxOverE" ]      = 0;
    mvaVar[ "SeedE1x3OverE" ]         = 0;
    mvaVar[ "SeedE3x1OverE" ]         = 0;
    mvaVar[ "SeedE1x5OverE" ]         = 0;
    mvaVar[ "SeedE2x2OverE" ]         = 0;
    mvaVar[ "SeedE3x2OverE" ]         = 0;
    mvaVar[ "SeedE3x3OverE" ]         = 0;
    mvaVar[ "SeedE4x4OverE" ]         = 0;
    mvaVar[ "SeedE5x5OverE" ]         = 0;
    mvaVar[ "ChargedIso03" ]          = 0; //Isolation 
    mvaVar[ "NeutralHadronIso03" ]    = 0; 
    mvaVar[ "GammaIso03" ]            = 0; 
    mvaVar[ "ChargedIso04" ]          = 0; 
    mvaVar[ "NeutralHadronIso04" ]    = 0; 
    mvaVar[ "GammaIso04" ]            = 0; 
    mvaVar[ "TrkIso03" ]              = 1; 
    mvaVar[ "EMIso03" ]               = 1; 
    mvaVar[ "HadIso03" ]              = 1; 
    mvaVar[ "TrkIso04" ]              = 1; 
    mvaVar[ "EMIso04" ]               = 1; 
    mvaVar[ "HadIso04" ]              = 1; 
  }


  //---------------------------------
  //choose bkg samples to include
  //---------------------------------
  
  TChain *chbackground = new TChain("Electrons");
  chbackground->Add(bgdInputFile);

  //---------------------------------
  //choose signal sample to include
  //---------------------------------

  TChain *chsignal = new TChain("Electrons");
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

  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
//   (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";


  //Allow more variables
  (TMVA::gConfig().GetVariablePlotting()).fMaxNumOfAllowedVariablesForScatterPlots = 50;

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  //factory->AddVariable( "myvar1 := var1+var2", 'F' );
  //factory->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
  //factory->AddVariable( "var3",                "Variable 3", "units", 'F' );
  //factory->AddVariable( "var4",                "Variable 4", "units", 'F' );

  if (mvaVar["SigmaIEtaIEta"])         factory->AddVariable( "SigmaIEtaIEta",           "SigmaIEtaIEta",         "", 'F' );
  if (mvaVar["DEtaIn"])                factory->AddVariable( "DEtaIn",                  "DEtaIn",                "", 'F' );
  if (mvaVar["DPhiIn"])                factory->AddVariable( "DPhiIn",                  "DPhiIn",                "", 'F' );
  if (mvaVar["D0"])                    factory->AddVariable( "D0",                      "D0",                    "", 'F' );
  if (mvaVar["DZ"])                    factory->AddVariable( "DZ",                      "DZ",                    "", 'F' );
  if (mvaVar["FBrem"])                 factory->AddVariable( "FBrem",                   "FBrem",                 "", 'F' );
  if (mvaVar["EOverP"])                factory->AddVariable( "EOverP",                  "EOverP",                "", 'F' );
  if (mvaVar["ESeedClusterOverPout"])  factory->AddVariable( "ESeedClusterOverPout",    "ESeedClusterOverPout",  "", 'F' );
  if (mvaVar["SigmaIPhiIPhi"])         factory->AddVariable( "SigmaIPhiIPhi",           "SigmaIPhiIPhi",         "", 'F' );
  if (mvaVar["NBrem"])                 factory->AddVariable( "NBrem",                   "NBrem",                 "", 'I' );
  if (mvaVar["OneOverEMinusOneOverP"]) factory->AddVariable( "OneOverEMinusOneOverP",   "OneOverEMinusOneOverP", "", 'F' );
  if (mvaVar["ESeedClusterOverPIn"])   factory->AddVariable( "ESeedClusterOverPIn",     "ESeedClusterOverPIn",   "", 'F' );
  if (mvaVar["IP3d"])                  factory->AddVariable( "IP3d",                    "IP3d",                  "", 'F' );
  if (mvaVar["IP3dSig"])               factory->AddVariable( "IP3dSig",                 "IP3dSig",               "", 'F' );
  if (mvaVar["StandardLikelihood"])    factory->AddVariable( "StandardLikelihood",      "StandardLikelihood",    "", 'F' );
  if (mvaVar["PFMVA"])                 factory->AddVariable( "PFMVA",                   "PFMVA",                 "", 'F' );
  if (mvaVar["TMVAKNN"])               factory->AddVariable( "TMVAKNN",                 "TMVAKNN",               "", 'F' );
  if (mvaVar["TMVAMLP"])               factory->AddVariable( "TMVAMLP",                 "TMVAMLP",               "", 'F' );
  if (mvaVar["TMVAMLPBNN"])            factory->AddVariable( "TMVAMLPBNN",              "TMVAMLPBNN",            "", 'F' );
  if (mvaVar["TMVABDTG"])              factory->AddVariable( "TMVABDTG",                "TMVABDTG",              "", 'F' );
  if (mvaVar["TMVABDT"])               factory->AddVariable( "TMVABDT",                 "TMVABDT",               "", 'F' );
  if (mvaVar["GsfTrackChi2OverNdof"])  factory->AddVariable( "GsfTrackChi2OverNdof",    "GsfTrackChi2OverNdof",  "", 'F' );
  if (mvaVar["dEtaCalo"])              factory->AddVariable( "dEtaCalo",                "dEtaCalo",              "", 'F' );
  if (mvaVar["dPhiCalo"])              factory->AddVariable( "dPhiCalo",                "dPhiCalo",              "", 'F' );
  if (mvaVar["R9"])                    factory->AddVariable( "R9",                      "R9",                    "", 'F' );
  if (mvaVar["SCEtaWidth"])            factory->AddVariable( "SCEtaWidth",              "SCEtaWidth",            "", 'F' );
  if (mvaVar["SCPhiWidth"])            factory->AddVariable( "SCPhiWidth",              "SCPhiWidth",            "", 'F' );
  if (mvaVar["CovIEtaIPhi"])           factory->AddVariable( "CovIEtaIPhi",             "CovIEtaIPhi",           "", 'F' );
  if (Option == 2 || Option == 5) {
    if (mvaVar["PreShowerOverRaw"])      factory->AddVariable( "PreShowerOverRaw",        "PreShowerOverRaw",      "", 'F' );
  }
  if (mvaVar["HoverE"])                factory->AddVariable( "HoverE",                  "HoverE",                "", 'F' );
  if (mvaVar["HcalDepth1OverEcal"])    factory->AddVariable( "HcalDepth1OverEcal",      "HcalDepth1OverEcal",    "", 'F' );
  if (Option == 1 || Option == 2 || Option == 4 || Option == 5 ) {
    if (mvaVar["HcalDepth2OverEcal"])    factory->AddVariable( "HcalDepth2OverEcal",      "HcalDepth2OverEcal",    "", 'F' );
  }
  if (mvaVar["SeedEMaxOverE"])         factory->AddVariable( "SeedEMaxOverE",           "SeedEMaxOverE",         "", 'F' );
  if (mvaVar["SeedETopOverE"])         factory->AddVariable( "SeedETopOverE",           "SeedETopOverE",         "", 'F' );
  if (mvaVar["SeedEBottomOverE"])      factory->AddVariable( "SeedEBottomOverE",        "SeedEBottomOverE",      "", 'F' );
  if (mvaVar["SeedELeftOverE"])        factory->AddVariable( "SeedELeftOverE",          "SeedELeftOverE",        "", 'F' );
  if (mvaVar["SeedERightOverE"])       factory->AddVariable( "SeedERightOverE",         "SeedERightOverE",       "", 'F' );
  if (mvaVar["SeedE2ndOverE"])         factory->AddVariable( "SeedE2ndOverE",           "SeedE2ndOverE",         "", 'F' );
  if (mvaVar["SeedE2x5RightOverE"])    factory->AddVariable( "SeedE2x5RightOverE",      "SeedE2x5RightOverE",    "", 'F' );
  if (mvaVar["SeedE2x5LeftOverE"])     factory->AddVariable( "SeedE2x5LeftOverE",       "SeedE2x5LeftOverE",     "", 'F' );
  if (mvaVar["SeedE2x5TopOverE"])      factory->AddVariable( "SeedE2x5TopOverE",        "SeedE2x5TopOverE",      "", 'F' );
  if (mvaVar["SeedE2x5BottomOverE"])   factory->AddVariable( "SeedE2x5BottomOverE",     "SeedE2x5BottomOverE",   "", 'F' );
  if (mvaVar["SeedE2x5MaxOverE"])      factory->AddVariable( "SeedE2x5MaxOverE",        "SeedE2x5MaxOverE",      "", 'F' );
  if (mvaVar["SeedE1x3OverE"])         factory->AddVariable( "SeedE1x3OverE",           "SeedE1x3OverE",         "", 'F' );
  if (mvaVar["SeedE3x1OverE"])         factory->AddVariable( "SeedE3x1OverE",           "SeedE3x1OverE",         "", 'F' );
  if (mvaVar["SeedE1x5OverE"])         factory->AddVariable( "SeedE1x5OverE",           "SeedE1x5OverE",         "", 'F' );
  if (mvaVar["SeedE2x2OverE"])         factory->AddVariable( "SeedE2x2OverE",           "SeedE2x2OverE",         "", 'F' );
  if (mvaVar["SeedE3x2OverE"])         factory->AddVariable( "SeedE3x2OverE",           "SeedE3x2OverE",         "", 'F' );
  if (mvaVar["SeedE3x3OverE"])         factory->AddVariable( "SeedE3x3OverE",           "SeedE3x3OverE",         "", 'F' );
  if (mvaVar["SeedE4x4OverE"])         factory->AddVariable( "SeedE4x4OverE",           "SeedE4x4OverE",         "", 'F' );
  if (mvaVar["SeedE5x5OverE"])         factory->AddVariable( "SeedE5x5OverE",           "SeedE5x5OverE",         "", 'F' );
  if (mvaVar["ChargedIso03"])          factory->AddVariable( "ChargedIso03",            "ChargedIso03",          "", 'F' );
  if (mvaVar["NeutralHadronIso03"])    factory->AddVariable( "NeutralHadronIso03",      "NeutralHadronIso03",    "", 'F' );
  if (mvaVar["GammaIso03"])            factory->AddVariable( "GammaIso03",              "GammaIso03",            "", 'F' );
  if (mvaVar["ChargedIso04"])          factory->AddVariable( "ChargedIso04",            "ChargedIso04",          "", 'F' );
  if (mvaVar["NeutralHadronIso04"])    factory->AddVariable( "NeutralHadronIso04",      "NeutralHadronIso04",    "", 'F' );
  if (mvaVar["GammaIso04"])            factory->AddVariable( "GammaIso04",              "GammaIso04",            "", 'F' );
  if (mvaVar["TrkIso03"])              factory->AddVariable( "TrkIso03",                "TrkIso03",              "", 'F' );
  if (mvaVar["EMIso03"])               factory->AddVariable( "EMIso03",                 "EMIso03",               "", 'F' );
  if (mvaVar["HadIso03"])              factory->AddVariable( "HadIso03",                "HadIso03",              "", 'F' );
  if (mvaVar["TrkIso04"])              factory->AddVariable( "TrkIso04",                "TrkIso04",              "", 'F' );
  if (mvaVar["EMIso04"])               factory->AddVariable( "EMIso04",                 "EMIso04",               "", 'F' );
  if (mvaVar["HadIso04"])              factory->AddVariable( "HadIso04",                "HadIso04",              "", 'F' );

  int nVariablesTemp = 0;

  if (mvaVar["SigmaIEtaIEta"])          { cout << "Adding variable to MVA training: SigmaIEtaIEta"         << endl; nVariablesTemp++; }
  if (mvaVar["DEtaIn"])                 { cout << "Adding variable to MVA training: DEtaIn"                << endl; nVariablesTemp++; }
  if (mvaVar["DPhiIn"])                 { cout << "Adding variable to MVA training: DPhiIn"                << endl; nVariablesTemp++; }
  if (mvaVar["D0"])                     { cout << "Adding variable to MVA training: D0"                    << endl; nVariablesTemp++; }
  if (mvaVar["DZ"])                     { cout << "Adding variable to MVA training: DZ"                    << endl; nVariablesTemp++; }
  if (mvaVar["FBrem"])                  { cout << "Adding variable to MVA training: FBrem"                 << endl; nVariablesTemp++; }
  if (mvaVar["EOverP"])                 { cout << "Adding variable to MVA training: EOverP"                << endl; nVariablesTemp++; }
  if (mvaVar["ESeedClusterOverPout"])   { cout << "Adding variable to MVA training: ESeedClusterOverPout"  << endl; nVariablesTemp++; }
  if (mvaVar["SigmaIPhiIPhi"])          { cout << "Adding variable to MVA training: SigmaIPhiIPhi"         << endl; nVariablesTemp++; }
  if (mvaVar["NBrem"])                  { cout << "Adding variable to MVA training: NBrem"                 << endl; nVariablesTemp++; }
  if (mvaVar["OneOverEMinusOneOverP"])  { cout << "Adding variable to MVA training: OneOverEMinusOneOverP" << endl; nVariablesTemp++; }
  if (mvaVar["ESeedClusterOverPIn"])    { cout << "Adding variable to MVA training: ESeedClusterOverPIn"   << endl; nVariablesTemp++; }
  if (mvaVar["IP3d"])                   { cout << "Adding variable to MVA training: IP3d"                  << endl; nVariablesTemp++; }
  if (mvaVar["IP3dSig"])                { cout << "Adding variable to MVA training: IP3dSig"               << endl; nVariablesTemp++; }
  if (mvaVar["PtLessThan15"])           { cout << "Adding variable to MVA training: PtLessThan15"          << endl; nVariablesTemp++; }
  if (mvaVar["StandardLikelihood"])     { cout << "Adding variable to MVA training: StandardLikelihood"    << endl; nVariablesTemp++; }
  if (mvaVar["PFMVA"])                  { cout << "Adding variable to MVA training: PFMVA"                 << endl; nVariablesTemp++; }
  if (mvaVar["TMVAKNN"])                { cout << "Adding variable to MVA training: TMVAKNN"               << endl; nVariablesTemp++; }
  if (mvaVar["TMVAMLP"])                { cout << "Adding variable to MVA training: TMVAMLP"               << endl; nVariablesTemp++; }
  if (mvaVar["TMVAMLPBNN"])             { cout << "Adding variable to MVA training: TMVAMLPBNN"            << endl; nVariablesTemp++; }
  if (mvaVar["TMVABDTG"])               { cout << "Adding variable to MVA training: TMVABDTG"              << endl; nVariablesTemp++; }
  if (mvaVar["TMVABDT"])                { cout << "Adding variable to MVA training: TMVABDT"               << endl; nVariablesTemp++; }
  if (mvaVar["GsfTrackChi2OverNdof"])   { cout << "Adding variable to MVA training: GsfTrackChi2OverNdof"  << endl; nVariablesTemp++; }
  if (mvaVar["dEtaCalo"])               { cout << "Adding variable to MVA training: dEtaCalo"              << endl; nVariablesTemp++; }
  if (mvaVar["dPhiCalo"])               { cout << "Adding variable to MVA training: dPhiCalo"              << endl; nVariablesTemp++; }
  if (mvaVar["R9"])                     { cout << "Adding variable to MVA training: R9"                    << endl; nVariablesTemp++; }
  if (mvaVar["SCEtaWidth"])             { cout << "Adding variable to MVA training: SCEtaWidth"            << endl; nVariablesTemp++; }
  if (mvaVar["SCPhiWidth"])             { cout << "Adding variable to MVA training: SCPhiWidth"            << endl; nVariablesTemp++; }
  if (mvaVar["CovIEtaIPhi"])            { cout << "Adding variable to MVA training: CovIEtaIPhi"           << endl; nVariablesTemp++; }
  if (Option == 2 || Option == 5) {
    if (mvaVar["PreShowerOverRaw"])       { cout << "Adding variable to MVA training: PreShowerOverRaw"      << endl; nVariablesTemp++; }
  }
  if (mvaVar["HoverE"])                 { cout << "Adding variable to MVA training: HoverE"                << endl; nVariablesTemp++; }
  if (mvaVar["HcalDepth1OverEcal"])     { cout << "Adding variable to MVA training: HcalDepth1OverEcal"    << endl; nVariablesTemp++; }
  if (Option == 1 || Option == 2 || Option == 4 || Option == 5 ) {
    if (mvaVar["HcalDepth2OverEcal"])     { cout << "Adding variable to MVA training: HcalDepth2OverEcal"    << endl; nVariablesTemp++; }
  }
  if (mvaVar["SeedEMaxOverE"])          { cout << "Adding variable to MVA training: SeedEMaxOverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedETopOverE"])          { cout << "Adding variable to MVA training: SeedETopOverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedEBottomOverE"])       { cout << "Adding variable to MVA training: SeedEBottomOverE"      << endl; nVariablesTemp++; }
  if (mvaVar["SeedELeftOverE"])         { cout << "Adding variable to MVA training: SeedELeftOverE"        << endl; nVariablesTemp++; }
  if (mvaVar["SeedERightOverE"])        { cout << "Adding variable to MVA training: SeedERightOverE"       << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2ndOverE"])          { cout << "Adding variable to MVA training: SeedE2ndOverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x5RightOverE"])     { cout << "Adding variable to MVA training: SeedE2x5RightOverE"    << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x5LeftOverE"])      { cout << "Adding variable to MVA training: SeedE2x5LeftOverE"     << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x5TopOverE"])       { cout << "Adding variable to MVA training: SeedE2x5TopOverE"      << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x5BottomOverE"])    { cout << "Adding variable to MVA training: SeedE2x5BottomOverE"   << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x5MaxOverE"])       { cout << "Adding variable to MVA training: SeedE2x5MaxOverE"      << endl; nVariablesTemp++; }
  if (mvaVar["SeedE1x3OverE"])          { cout << "Adding variable to MVA training: SeedE1x3OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE3x1OverE"])          { cout << "Adding variable to MVA training: SeedE3x1OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE1x5OverE"])          { cout << "Adding variable to MVA training: SeedE1x5OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE2x2OverE"])          { cout << "Adding variable to MVA training: SeedE2x2OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE3x2OverE"])          { cout << "Adding variable to MVA training: SeedE3x2OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE3x3OverE"])          { cout << "Adding variable to MVA training: SeedE3x3OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE4x4OverE"])          { cout << "Adding variable to MVA training: SeedE4x4OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["SeedE5x5OverE"])          { cout << "Adding variable to MVA training: SeedE5x5OverE"         << endl; nVariablesTemp++; }
  if (mvaVar["ChargedIso03"])           { cout << "Adding variable to MVA training: ChargedIso03"          << endl; nVariablesTemp++; }
  if (mvaVar["NeutralHadronIso03"])     { cout << "Adding variable to MVA training: NeutralHadronIso03"    << endl; nVariablesTemp++; }
  if (mvaVar["GammaIso03"])             { cout << "Adding variable to MVA training: GammaIso03"            << endl; nVariablesTemp++; }
  if (mvaVar["ChargedIso04"])           { cout << "Adding variable to MVA training: ChargedIso04"          << endl; nVariablesTemp++; }
  if (mvaVar["NeutralHadronIso04"])     { cout << "Adding variable to MVA training: NeutralHadronIso04"    << endl; nVariablesTemp++; }
  if (mvaVar["GammaIso04"])             { cout << "Adding variable to MVA training: GammaIso04"            << endl; nVariablesTemp++; }
  if (mvaVar["TrkIso03"])               { cout << "Adding variable to MVA training: TrkIso03"              << endl; nVariablesTemp++; }
  if (mvaVar["EMIso03"])                { cout << "Adding variable to MVA training: EMIso03"               << endl; nVariablesTemp++; }
  if (mvaVar["HadIso03"])               { cout << "Adding variable to MVA training: HadIso03"              << endl; nVariablesTemp++; }
  if (mvaVar["TrkIso04"])               { cout << "Adding variable to MVA training: TrkIso04"              << endl; nVariablesTemp++; }
  if (mvaVar["EMIso04"])                { cout << "Adding variable to MVA training: EMIso04"               << endl; nVariablesTemp++; }
  if (mvaVar["HadIso04"])               { cout << "Adding variable to MVA training: HadIso04"              << endl; nVariablesTemp++; }

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

  Float_t                 fWeight;
  UInt_t                  fRunNumber;
  UInt_t                  fLumiSectionNumber;
  UInt_t                  fEventNumber;
  Float_t                 fElePt; 
  Float_t                 fEleEta; 
  Float_t                 fElePhi; 
  Float_t                 fEleSCEt; 
  Float_t                 fEleSCEta; 
  Float_t                 fEleSCPhi; 
  Float_t                 fElePFIso; 

  Float_t                 fEleSigmaIEtaIEta; 
  Float_t                 fEleDEtaIn; 
  Float_t                 fEleDPhiIn; 
  Float_t                 fEleHoverE; 
  Float_t                 fEleD0; 
  Float_t                 fEleDZ; 
  Float_t                 fEleFBrem; 
  Float_t                 fEleEOverP; 

  Float_t                 fEleESeedClusterOverPout; 
  Float_t                 fEleSigmaIPhiIPhi; 
  Float_t                 fEleNBrem; 
  Float_t                 fEleOneOverEMinusOneOverP; 
  Float_t                 fEleESeedClusterOverPIn; 
  Float_t                 fEleIP3d; 
  Float_t                 fEleIP3dSig; 

  Float_t                 fEleHcalDepth1OverEcal;
  Float_t                 fEleHcalDepth2OverEcal;
  Float_t                 fEledEtaCalo;
  Float_t                 fEledPhiCalo;
  Float_t                 fElePreShowerOverRaw;
  Float_t                 fEleCovIEtaIPhi;
  Float_t                 fEleSCEtaWidth;
  Float_t                 fEleSCPhiWidth;
  Float_t                 fEleGsfTrackChi2OverNdof;
  Float_t                 fEleR9;

  Float_t                 fEleSeedEMaxOverE;
  Float_t                 fEleSeedETopOverE;
  Float_t                 fEleSeedEBottomOverE;
  Float_t                 fEleSeedELeftOverE;
  Float_t                 fEleSeedERightOverE;
  Float_t                 fEleSeedE2ndOverE;
  Float_t                 fEleSeedE2x5RightOverE;
  Float_t                 fEleSeedE2x5LeftOverE;
  Float_t                 fEleSeedE2x5TopOverE;
  Float_t                 fEleSeedE2x5BottomOverE;
  Float_t                 fEleSeedE2x5MaxOverE;
  Float_t                 fEleSeedE1x3OverE;
  Float_t                 fEleSeedE3x1OverE;
  Float_t                 fEleSeedE1x5OverE;
  Float_t                 fEleSeedE2x2OverE;
  Float_t                 fEleSeedE3x2OverE;
  Float_t                 fEleSeedE3x3OverE;
  Float_t                 fEleSeedE4x4OverE;
  Float_t                 fEleSeedE5x5OverE;

  Float_t                 fEleStandardLikelihood; 
  Float_t                 fElePFMVA; 
  Float_t                 fEleKNN; 
  Float_t                 fEleMLP; 
  Float_t                 fEleMLPBNN; 
  Float_t                 fEleBDTG; 
  Float_t                 fEleBDT; 

  Float_t                 fEleChargedIso03; 
  Float_t                 fEleNeutralHadronIso03; 
  Float_t                 fEleGammaIso03; 
  Float_t                 fEleChargedIso04; 
  Float_t                 fEleNeutralHadronIso04; 
  Float_t                 fEleGammaIso04; 
  Float_t                 fEleChargedIso04FromOtherVertices; 
  Float_t                 fEleNeutralHadronIso04_10Threshold; 
  Float_t                 fEleGammaIso04_10Threshold; 
  Float_t                 fEleTrkIso03; 
  Float_t                 fEleEMIso03; 
  Float_t                 fEleHadIso03; 
  Float_t                 fEleTrkIso04; 
  Float_t                 fEleEMIso04; 
  Float_t                 fEleHadIso04; 
  Float_t                 fRho; 
  Float_t                 fNVertices; 

  signal->SetBranchAddress( "weight", &fWeight);
  signal->SetBranchAddress( "run", &fRunNumber);
  signal->SetBranchAddress( "lumi", &fLumiSectionNumber);
  signal->SetBranchAddress( "event", &fEventNumber);
  signal->SetBranchAddress( "pt", &fElePt); 
  signal->SetBranchAddress( "eta", &fEleEta); 
  signal->SetBranchAddress( "phi", &fElePhi); 
  signal->SetBranchAddress( "scet", &fEleSCEt); 
  signal->SetBranchAddress( "sceta", &fEleSCEta); 
  signal->SetBranchAddress( "scphi", &fEleSCPhi); 
  signal->SetBranchAddress( "pfiso", &fElePFIso); 
  signal->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  signal->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  signal->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  signal->SetBranchAddress( "HoverE", &fEleHoverE); 
  signal->SetBranchAddress( "D0", &fEleD0); 
  signal->SetBranchAddress( "DZ", &fEleDZ); 
  signal->SetBranchAddress( "FBrem", &fEleFBrem); 
  signal->SetBranchAddress( "EOverP", &fEleEOverP); 
  signal->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  signal->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  signal->SetBranchAddress( "NBrem", &fEleNBrem); 
  signal->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  signal->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  signal->SetBranchAddress( "IP3d", &fEleIP3d); 
  signal->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  signal->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  signal->SetBranchAddress( "PFMVA", &fElePFMVA); 
  signal->SetBranchAddress( ("KNN"+versionLabel).c_str(), &fEleKNN); 
  signal->SetBranchAddress( ("MLP"+versionLabel).c_str(), &fEleMLP); 
  signal->SetBranchAddress( ("MLPBNN"+versionLabel).c_str(), &fEleMLPBNN); 
  signal->SetBranchAddress( ("BDTG"+versionLabel).c_str(), &fEleBDTG); 
  signal->SetBranchAddress( ("BDT"+versionLabel).c_str(), &fEleBDT); 
  signal->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  signal->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  signal->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  signal->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  signal->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  signal->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  signal->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  signal->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  signal->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  signal->SetBranchAddress( "R9", &fEleR9); 
  signal->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  signal->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  signal->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  signal->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  signal->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  signal->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  signal->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  signal->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  signal->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  signal->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  signal->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  signal->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  signal->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  signal->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  signal->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  signal->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  signal->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  signal->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  signal->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  signal->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  signal->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  signal->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  signal->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  signal->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  signal->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  signal->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  signal->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  signal->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  signal->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  signal->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  signal->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  signal->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  signal->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  signal->SetBranchAddress( "HadIso04", &fEleHadIso04); 
  signal->SetBranchAddress( "Rho", &fRho); 
  signal->SetBranchAddress( "NVertices", &fNVertices); 

  int nsigtrain = 0;
  int nsigtest  = 0;
  int nbkgtrain = 0;
  int nbkgtest  = 0;


  for (UInt_t i=0; i<signal->GetEntries(); i++) {
    
    signal->GetEntry(i);

    Double_t rho = 0;
    if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;
    
   //--------------------------------------------------
    // SIGNAL EVENT SELECTION
    //--------------------------------------------------
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

    Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0 && fEleFBrem < 0.1);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0 && fEleFBrem < 0.1);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    


    if( fElePt < 10.0                    ) continue; // cut on dilepton mass
    if( fabs(fEleEta) > 2.5              ) continue; // cut on opposite-sign leptons
    if( fabs(fEleDZ) > 0.1               ) continue; // cut on low dilepton mass
    if( fabs(fEleD0) > 0.02    	         ) continue; // cut on leading lepton pt
    if( fElePt > 35.0                    ) continue; 


    int varCounter = 0;

    if (mvaVar["SigmaIEtaIEta"])          vars[varCounter++] = fEleSigmaIEtaIEta;
    if (mvaVar["DEtaIn"])                 vars[varCounter++] = fEleDEtaIn;
    if (mvaVar["DPhiIn"])                 vars[varCounter++] = fEleDPhiIn;
    if (mvaVar["D0"])                     vars[varCounter++] = fEleD0;
    if (mvaVar["DZ"])                     vars[varCounter++] = fEleDZ;
    if (mvaVar["FBrem"])                  vars[varCounter++] = fEleFBrem;
    if (mvaVar["EOverP"])                 vars[varCounter++] = fEleEOverP;
    if (mvaVar["ESeedClusterOverPout"])   vars[varCounter++] = fEleESeedClusterOverPout;
    if (mvaVar["SigmaIPhiIPhi"])          vars[varCounter++] = fEleSigmaIPhiIPhi;
    if (mvaVar["NBrem"])                  vars[varCounter++] = fEleNBrem;
    if (mvaVar["OneOverEMinusOneOverP"])  vars[varCounter++] = fEleOneOverEMinusOneOverP;
    if (mvaVar["ESeedClusterOverPIn"])    vars[varCounter++] = fEleESeedClusterOverPIn;
    if (mvaVar["IP3d"])                   vars[varCounter++] = fEleIP3d;
    if (mvaVar["IP3dSig"])                vars[varCounter++] = fEleIP3dSig;
    if (mvaVar["StandardLikelihood"])     vars[varCounter++] = fEleStandardLikelihood;
    if (mvaVar["PFMVA"])                  vars[varCounter++] = fElePFMVA;
    if (mvaVar["TMVAKNN"])                vars[varCounter++] = fEleKNN;
    if (mvaVar["TMVAMLP"])                vars[varCounter++] = fEleMLP;
    if (mvaVar["TMVAMLPBNN"])             vars[varCounter++] = fEleMLPBNN;
    if (mvaVar["TMVABDTG"])               vars[varCounter++] = fEleBDTG;
    if (mvaVar["TMVABDT"])                vars[varCounter++] = fEleBDT;
    if (mvaVar["GsfTrackChi2OverNdof"])   vars[varCounter++] = fEleGsfTrackChi2OverNdof;
    if (mvaVar["dEtaCalo"])               vars[varCounter++] = fEledEtaCalo;
    if (mvaVar["dPhiCalo"])               vars[varCounter++] = fEledPhiCalo;
    if (mvaVar["R9"])                     vars[varCounter++] = fEleR9;
    if (mvaVar["SCEtaWidth"])             vars[varCounter++] = fEleSCEtaWidth;
    if (mvaVar["SCPhiWidth"])             vars[varCounter++] = fEleSCPhiWidth;
    if (mvaVar["CovIEtaIPhi"])            vars[varCounter++] = fEleCovIEtaIPhi;
    if (Option == 2 || Option == 5) {
      if (mvaVar["PreShowerOverRaw"])       vars[varCounter++] = fElePreShowerOverRaw;
    }
    if (mvaVar["HoverE"])                 vars[varCounter++] = fEleHoverE - rho*ElectronEffectiveArea(kEleHoverE, fEleEta);
    if (mvaVar["HcalDepth1OverEcal"])     vars[varCounter++] = fEleHcalDepth1OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth1OverEcal, fEleEta);
    if (Option == 1 || Option == 2 || Option == 4 || Option == 5 ) {
      if (mvaVar["HcalDepth2OverEcal"])     vars[varCounter++] = fEleHcalDepth2OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth2OverEcal, fEleEta);
    }
    if (mvaVar["SeedEMaxOverE"])          vars[varCounter++] = fEleSeedEMaxOverE;
    if (mvaVar["SeedETopOverE"])          vars[varCounter++] = fEleSeedETopOverE;
    if (mvaVar["SeedEBottomOverE"])       vars[varCounter++] = fEleSeedEBottomOverE;
    if (mvaVar["SeedELeftOverE"])         vars[varCounter++] = fEleSeedELeftOverE;
    if (mvaVar["SeedERightOverE"])        vars[varCounter++] = fEleSeedERightOverE;
    if (mvaVar["SeedE2ndOverE"])          vars[varCounter++] = fEleSeedE2ndOverE;
    if (mvaVar["SeedE2x5RightOverE"])     vars[varCounter++] = fEleSeedE2x5RightOverE;
    if (mvaVar["SeedE2x5LeftOverE"])      vars[varCounter++] = fEleSeedE2x5LeftOverE;
    if (mvaVar["SeedE2x5TopOverE"])       vars[varCounter++] = fEleSeedE2x5TopOverE;
    if (mvaVar["SeedE2x5BottomOverE"])    vars[varCounter++] = fEleSeedE2x5BottomOverE;
    if (mvaVar["SeedE2x5MaxOverE"])       vars[varCounter++] = fEleSeedE2x5MaxOverE;
    if (mvaVar["SeedE1x3OverE"])          vars[varCounter++] = fEleSeedE1x3OverE;
    if (mvaVar["SeedE3x1OverE"])          vars[varCounter++] = fEleSeedE3x1OverE;
    if (mvaVar["SeedE1x5OverE"])          vars[varCounter++] = fEleSeedE1x5OverE;
    if (mvaVar["SeedE2x2OverE"])          vars[varCounter++] = fEleSeedE2x2OverE;
    if (mvaVar["SeedE3x2OverE"])          vars[varCounter++] = fEleSeedE3x2OverE;
    if (mvaVar["SeedE3x3OverE"])          vars[varCounter++] = fEleSeedE3x3OverE;
    if (mvaVar["SeedE4x4OverE"])          vars[varCounter++] = fEleSeedE4x4OverE;
    if (mvaVar["SeedE5x5OverE"])          vars[varCounter++] = fEleSeedE5x5OverE;

    if (versionLabel == "V15") {
      if (mvaVar["ChargedIso03"])           vars[varCounter++] = (fEleChargedIso03- rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta))/fElePt;
      if (mvaVar["NeutralHadronIso03"])     vars[varCounter++] = (fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
      if (mvaVar["GammaIso03"])             vars[varCounter++] = (fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta))/fElePt;
      if (mvaVar["ChargedIso04"])           vars[varCounter++] = (fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta))/fElePt;
      if (mvaVar["NeutralHadronIso04"])     vars[varCounter++] = (fEleNeutralHadronIso04- rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
      if (mvaVar["GammaIso04"])             vars[varCounter++] = (fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta))/fElePt;      
    } else {
      if (mvaVar["ChargedIso03"])           vars[varCounter++] = fEleChargedIso03- rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta);
      if (mvaVar["NeutralHadronIso03"])     vars[varCounter++] = fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
      if (mvaVar["GammaIso03"])             vars[varCounter++] = fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta);
      if (mvaVar["ChargedIso04"])           vars[varCounter++] = fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta);
      if (mvaVar["NeutralHadronIso04"])     vars[varCounter++] = fEleNeutralHadronIso04- rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
      if (mvaVar["GammaIso04"])             vars[varCounter++] = fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta);
    }

    if (mvaVar["TrkIso03"])             vars[varCounter++] = (fEleTrkIso03 - rho*ElectronEffectiveArea(kEleTrkIso03,fEleEta))/fElePt;
    if (mvaVar["EMIso03"])              vars[varCounter++] = (fEleEMIso03 - rho*ElectronEffectiveArea(kEleEMIso03,fEleEta))/fElePt;
    if (mvaVar["HadIso03"])             vars[varCounter++] = (fEleHadIso03 - rho*ElectronEffectiveArea(kEleHadIso03,fEleEta))/fElePt;
    if (mvaVar["TrkIso04"])             vars[varCounter++] = (fEleTrkIso04 - rho*ElectronEffectiveArea(kEleTrkIso04,fEleEta))/fElePt;
    if (mvaVar["EMIso04"])              vars[varCounter++] = (fEleEMIso04 - rho*ElectronEffectiveArea(kEleEMIso04,fEleEta))/fElePt;
    if (mvaVar["HadIso04"])             vars[varCounter++] = (fEleHadIso04 - rho*ElectronEffectiveArea(kEleHadIso04,fEleEta))/fElePt;



//     cout << "Signal " << i << " : varCounter: " << varCounter << " --> " << nVariables << endl;
//     for (UInt_t k = 0; k < varCounter ; ++k) {
//       cout << vars[k] << endl;
//     }

    if ( fEventNumber % 2 == 0 ){
//       cout << "bkg hovere:" << fEleHoverE << " " << rho << " " << ElectronEffectiveArea(kEleHoverE, fEleEta) << " : " << fEleHoverE - rho*ElectronEffectiveArea(kEleHoverE, fEleEta) << endl;

      factory->AddSignalTrainingEvent( vars, fWeight );
//       factory->AddSignalTrainingEvent( vars, 1.0 );
      nsigtrain++;
    }
    else{
      factory->AddSignalTestEvent    ( vars, fWeight );
//       factory->AddSignalTestEvent    ( vars, 1.0 );
      nsigtest++;
    }
  }

  background->SetBranchAddress( "weight", &fWeight);
  background->SetBranchAddress( "run", &fRunNumber);
  background->SetBranchAddress( "lumi", &fLumiSectionNumber);
  background->SetBranchAddress( "event", &fEventNumber);
  background->SetBranchAddress( "pt", &fElePt); 
  background->SetBranchAddress( "eta", &fEleEta); 
  background->SetBranchAddress( "phi", &fElePhi); 
  background->SetBranchAddress( "scet", &fEleSCEt); 
  background->SetBranchAddress( "sceta", &fEleSCEta); 
  background->SetBranchAddress( "scphi", &fEleSCPhi); 
  background->SetBranchAddress( "pfiso", &fElePFIso); 
  background->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
  background->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
  background->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
  background->SetBranchAddress( "HoverE", &fEleHoverE); 
  background->SetBranchAddress( "D0", &fEleD0); 
  background->SetBranchAddress( "DZ", &fEleDZ); 
  background->SetBranchAddress( "FBrem", &fEleFBrem); 
  background->SetBranchAddress( "EOverP", &fEleEOverP); 
  background->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
  background->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
  background->SetBranchAddress( "NBrem", &fEleNBrem); 
  background->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
  background->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
  background->SetBranchAddress( "IP3d", &fEleIP3d); 
  background->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
  background->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
  background->SetBranchAddress( "PFMVA", &fElePFMVA); 
  background->SetBranchAddress( ("KNN"+versionLabel).c_str(), &fEleKNN); 
  background->SetBranchAddress( ("MLP"+versionLabel).c_str(), &fEleMLP); 
  background->SetBranchAddress( ("MLPBNN"+versionLabel).c_str(), &fEleMLPBNN); 
  background->SetBranchAddress( ("BDTG"+versionLabel).c_str(), &fEleBDTG); 
  background->SetBranchAddress( ("BDT"+versionLabel).c_str(), &fEleBDT); 
  background->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
  background->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
  background->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
  background->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
  background->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
  background->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
  background->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
  background->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
  background->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
  background->SetBranchAddress( "R9", &fEleR9); 
  background->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
  background->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
  background->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
  background->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
  background->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
  background->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
  background->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
  background->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
  background->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
  background->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
  background->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
  background->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
  background->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
  background->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
  background->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
  background->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
  background->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
  background->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
  background->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
  background->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
  background->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
  background->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
  background->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
  background->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
  background->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
  background->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
  background->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
  background->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
  background->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
  background->SetBranchAddress( "EMIso03", &fEleEMIso03); 
  background->SetBranchAddress( "HadIso03", &fEleHadIso03); 
  background->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
  background->SetBranchAddress( "EMIso04", &fEleEMIso04); 
  background->SetBranchAddress( "HadIso04", &fEleHadIso04); 
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
    Int_t subdet = 0;
    if (fabs(fEleEta) < 1.0) subdet = 0;
    else if (fabs(fEleEta) < 1.479) subdet = 1;
    else subdet = 2;
    Int_t ptBin = 0;
    if (fElePt > 20.0) ptBin = 1;

     Bool_t passCuts = kFALSE;
    if (Option == 0) passCuts = (subdet == 0 && ptBin == 0);
    if (Option == 1) passCuts = (subdet == 1 && ptBin == 0);
    if (Option == 2) passCuts = (subdet == 2 && ptBin == 0);
    if (Option == 3) passCuts = (subdet == 0 && ptBin == 1);
    if (Option == 4) passCuts = (subdet == 1 && ptBin == 1);
    if (Option == 5) passCuts = (subdet == 2 && ptBin == 1);    
    if (Option == 10) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem == 0 && fEleFBrem < 0.1);
    if (Option == 11) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 12) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem == 0);
    if (Option == 13) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem == 0 && fEleFBrem < 0.1);
    if (Option == 14) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 15) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem == 0);
    if (Option == 20) passCuts = (subdet == 0 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 21) passCuts = (subdet == 1 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 22) passCuts = (subdet == 2 && ptBin == 0 && fEleNBrem > 0);
    if (Option == 23) passCuts = (subdet == 0 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 24) passCuts = (subdet == 1 && ptBin == 1 && fEleNBrem > 0);
    if (Option == 25) passCuts = (subdet == 2 && ptBin == 1 && fEleNBrem > 0);
    if (Option == -1) passCuts = kTRUE;
    if (!passCuts) continue;    

    if( fElePt < 10.0                    ) continue; // cut on dilepton mass
    if( fabs(fEleEta) > 2.5              ) continue; // cut on opposite-sign leptons
    if( fabs(fEleDZ) > 0.1               ) continue; // cut on low dilepton mass
    if( fabs(fEleD0) > 0.02    	         ) continue; // cut on leading lepton pt
    if( fElePt > 35.0                    ) continue; 

//     if (!(fEleSigmaIEtaIEta >= 0.0 && fEleSigmaIEtaIEta <= 0.035)) continue;
//     if (!(fEleDEtaIn >= -0.015 && fEleDEtaIn <= 0.015)) continue;
//     if (!(fEleDPhiIn >= -0.1 && fEleDPhiIn <= 0.1)) continue;
//     if (!(fEleFBrem >= -0.2 && fEleFBrem <= 1.0)) continue;
// //     if (!(fEleSigmaIPhiIPhi >= 0.0 && fEleSigmaIPhiIPhi <= 0.03)) continue;
//     if (!(fEleOneOverEMinusOneOverP >= -1.0 && fEleOneOverEMinusOneOverP <= 0.05)) continue;

    int varCounter = 0;
    
    if (mvaVar["SigmaIEtaIEta"])          vars[varCounter++] = fEleSigmaIEtaIEta;
    if (mvaVar["DEtaIn"])                 vars[varCounter++] = fEleDEtaIn;
    if (mvaVar["DPhiIn"])                 vars[varCounter++] = fEleDPhiIn;
    if (mvaVar["D0"])                     vars[varCounter++] = fEleD0;
    if (mvaVar["DZ"])                     vars[varCounter++] = fEleDZ;
    if (mvaVar["FBrem"])                  vars[varCounter++] = fEleFBrem;
    if (mvaVar["EOverP"])                 vars[varCounter++] = fEleEOverP;
    if (mvaVar["ESeedClusterOverPout"])   vars[varCounter++] = fEleESeedClusterOverPout;
    if (mvaVar["SigmaIPhiIPhi"])          vars[varCounter++] = fEleSigmaIPhiIPhi;
    if (mvaVar["NBrem"])                  vars[varCounter++] = fEleNBrem;
    if (mvaVar["OneOverEMinusOneOverP"])  vars[varCounter++] = fEleOneOverEMinusOneOverP;
    if (mvaVar["ESeedClusterOverPIn"])    vars[varCounter++] = fEleESeedClusterOverPIn;
    if (mvaVar["IP3d"])                   vars[varCounter++] = fEleIP3d;
    if (mvaVar["IP3dSig"])                vars[varCounter++] = fEleIP3dSig;
    if (mvaVar["StandardLikelihood"])     vars[varCounter++] = fEleStandardLikelihood;
    if (mvaVar["PFMVA"])                  vars[varCounter++] = fElePFMVA;
    if (mvaVar["TMVAKNN"])                vars[varCounter++] = fEleKNN;
    if (mvaVar["TMVAMLP"])                vars[varCounter++] = fEleMLP;
    if (mvaVar["TMVAMLPBNN"])             vars[varCounter++] = fEleMLPBNN;
    if (mvaVar["TMVABDTG"])               vars[varCounter++] = fEleBDTG;
    if (mvaVar["TMVABDT"])                vars[varCounter++] = fEleBDT;
    if (mvaVar["GsfTrackChi2OverNdof"])   vars[varCounter++] = fEleGsfTrackChi2OverNdof;
    if (mvaVar["dEtaCalo"])               vars[varCounter++] = fEledEtaCalo;
    if (mvaVar["dPhiCalo"])               vars[varCounter++] = fEledPhiCalo;
    if (mvaVar["R9"])                     vars[varCounter++] = fEleR9;
    if (mvaVar["SCEtaWidth"])             vars[varCounter++] = fEleSCEtaWidth;
    if (mvaVar["SCPhiWidth"])             vars[varCounter++] = fEleSCPhiWidth;
    if (mvaVar["CovIEtaIPhi"])            vars[varCounter++] = fEleCovIEtaIPhi;
    if (Option == 2 || Option == 5) {
      if (mvaVar["PreShowerOverRaw"])       vars[varCounter++] = fElePreShowerOverRaw;
    }
    if (mvaVar["HoverE"])                 vars[varCounter++] = fEleHoverE - rho*ElectronEffectiveArea(kEleHoverE, fEleEta);
    if (mvaVar["HcalDepth1OverEcal"])     vars[varCounter++] = fEleHcalDepth1OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth1OverEcal, fEleEta);
    if (Option == 1 || Option == 2 || Option == 4 || Option == 5 ) {
      if (mvaVar["HcalDepth2OverEcal"])     vars[varCounter++] = fEleHcalDepth2OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth2OverEcal, fEleEta);
    }
    if (mvaVar["SeedEMaxOverE"])          vars[varCounter++] = fEleSeedEMaxOverE;
    if (mvaVar["SeedETopOverE"])          vars[varCounter++] = fEleSeedETopOverE;
    if (mvaVar["SeedEBottomOverE"])       vars[varCounter++] = fEleSeedEBottomOverE;
    if (mvaVar["SeedELeftOverE"])         vars[varCounter++] = fEleSeedELeftOverE;
    if (mvaVar["SeedERightOverE"])        vars[varCounter++] = fEleSeedERightOverE;
    if (mvaVar["SeedE2ndOverE"])          vars[varCounter++] = fEleSeedE2ndOverE;
    if (mvaVar["SeedE2x5RightOverE"])     vars[varCounter++] = fEleSeedE2x5RightOverE;
    if (mvaVar["SeedE2x5LeftOverE"])      vars[varCounter++] = fEleSeedE2x5LeftOverE;
    if (mvaVar["SeedE2x5TopOverE"])       vars[varCounter++] = fEleSeedE2x5TopOverE;
    if (mvaVar["SeedE2x5BottomOverE"])    vars[varCounter++] = fEleSeedE2x5BottomOverE;
    if (mvaVar["SeedE2x5MaxOverE"])       vars[varCounter++] = fEleSeedE2x5MaxOverE;
    if (mvaVar["SeedE1x3OverE"])          vars[varCounter++] = fEleSeedE1x3OverE;
    if (mvaVar["SeedE3x1OverE"])          vars[varCounter++] = fEleSeedE3x1OverE;
    if (mvaVar["SeedE1x5OverE"])          vars[varCounter++] = fEleSeedE1x5OverE;
    if (mvaVar["SeedE2x2OverE"])          vars[varCounter++] = fEleSeedE2x2OverE;
    if (mvaVar["SeedE3x2OverE"])          vars[varCounter++] = fEleSeedE3x2OverE;
    if (mvaVar["SeedE3x3OverE"])          vars[varCounter++] = fEleSeedE3x3OverE;
    if (mvaVar["SeedE4x4OverE"])          vars[varCounter++] = fEleSeedE4x4OverE;
    if (mvaVar["SeedE5x5OverE"])          vars[varCounter++] = fEleSeedE5x5OverE;
    if (versionLabel == "V15") {
      if (mvaVar["ChargedIso03"])           vars[varCounter++] = (fEleChargedIso03- rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta))/fElePt;
      if (mvaVar["NeutralHadronIso03"])     vars[varCounter++] = (fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
      if (mvaVar["GammaIso03"])             vars[varCounter++] = (fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta))/fElePt;
      if (mvaVar["ChargedIso04"])           vars[varCounter++] = (fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta))/fElePt;
      if (mvaVar["NeutralHadronIso04"])     vars[varCounter++] = (fEleNeutralHadronIso04- rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
      if (mvaVar["GammaIso04"])             vars[varCounter++] = (fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta))/fElePt;      
    } else {
      if (mvaVar["ChargedIso03"])           vars[varCounter++] = fEleChargedIso03- rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta);
      if (mvaVar["NeutralHadronIso03"])     vars[varCounter++] = fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
      if (mvaVar["GammaIso03"])             vars[varCounter++] = fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta);
      if (mvaVar["ChargedIso04"])           vars[varCounter++] = fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta);
      if (mvaVar["NeutralHadronIso04"])     vars[varCounter++] = fEleNeutralHadronIso04- rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
      if (mvaVar["GammaIso04"])             vars[varCounter++] = fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta);
    }

    if (mvaVar["TrkIso03"])             vars[varCounter++] = (fEleTrkIso03 - rho*ElectronEffectiveArea(kEleTrkIso03,fEleEta))/fElePt;
    if (mvaVar["EMIso03"])              vars[varCounter++] = (fEleEMIso03 - rho*ElectronEffectiveArea(kEleEMIso03,fEleEta))/fElePt;
    if (mvaVar["HadIso03"])             vars[varCounter++] = (fEleHadIso03 - rho*ElectronEffectiveArea(kEleHadIso03,fEleEta))/fElePt;
    if (mvaVar["TrkIso04"])             vars[varCounter++] = (fEleTrkIso04 - rho*ElectronEffectiveArea(kEleTrkIso04,fEleEta))/fElePt;
    if (mvaVar["EMIso04"])              vars[varCounter++] = (fEleEMIso04 - rho*ElectronEffectiveArea(kEleEMIso04,fEleEta))/fElePt;
    if (mvaVar["HadIso04"])             vars[varCounter++] = (fEleHadIso04 - rho*ElectronEffectiveArea(kEleHadIso04,fEleEta))/fElePt;


//     cout << "Signal " << i << " : varCounter: " << varCounter << " --> " << nVariables << endl;
//     cout << i << " : varCounter: " << varCounter << " : " << fEleChargedIso03 - rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta) << endl;
//     for (UInt_t k = 0; k < varCounter ; ++k) {
//       cout << vars[k] << endl;
//     }

    if ( fEventNumber % 2 == 0 ){
//       cout << "bkg hovere:" << fEleHoverE << " " << rho << " " << ElectronEffectiveArea(kEleHoverE, fEleEta) << " : " << fEleHoverE - rho*ElectronEffectiveArea(kEleHoverE, fEleEta) << endl;
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
                         "H:!V:CreateMVAPdfs:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=30:NSmoothBkg[0]=30:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=40:Nbins=200" );

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
    factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                         "!H:!V:CreateMVAPdfs:CreateMVAPdfs:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:NNodesMax=5:UseNvars=4:PruneStrength=0.0:PruneMethod=nopruning:MaxDepth=6" );
  

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
