//root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root","output/ElectronNtuple.Fake.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG,MLPBNN")'
//root -l -q -b /home/sixie/CMSSW_analysis/src/EWKAna/Hww/LeptonSelection/EvaluateElectronMVA.C+'("ElectronSelectionTraining.Real.weighted.Subdet0LowPt.root","output/ElectronNtuple.Real.Subdet0LowPt.root","Subdet0LowPt","Likelihood,LikelihoodD,BDT,BDTG,MLPBNN")'





/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TChain.h"
#include "TChainElement.h"
#include "TBranch.h"
#include "TRandom.h"

#include "Math/LorentzVector.h"
#include "Math/VectorUtil.h"
#include "TMath.h"
#include "EWKAna/Utils/LeptonIDCuts.hh"

#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector; 

//--------------------------------------------------------------------

void fillUnderOverFlow(TH1F *h1, float value, float weight = 1)
{
  float min = h1->GetXaxis()->GetXmin();
  float max = h1->GetXaxis()->GetXmax();

  if (value > max) value = h1->GetBinCenter(h1->GetNbinsX());
  if (value < min) value = h1->GetBinCenter(1);

  h1->Fill(value, weight);
}

//--------------------------------------------------------------------

void EvaluateElectronMVA(
string inputFile      = "ElectronSelectionTraining.Fake.weighted.Subdet0LowPt.root", 
string outputFile     = "output/ElectronNtuple.Fake.Subdet0LowPt.root",
TString label         = "Subdet0LowPt",
string versionLabel   = "V0",
Bool_t  buildCommittee  = kFALSE,
Bool_t evaluateAllBins = kFALSE,
Int_t Option            = -1,
TString myMethodList  = "Likelihood,LikelihoodD,BDT,BDTG,MLPBNN"
) {   
#ifdef __CINT__
  gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

  //--------------------------------------------------------------------
  // path to weights dir (this is where MVA training info is stored)
  // output root file will be stored at [path]/output
  //--------------------------------------------------------------------
  
  vector<string> samples;
  samples.push_back(inputFile);

  //--------------------------------------------------------------------------------
  // IMPORTANT: set the following variables to the same set used for MVA training!!!
  //--------------------------------------------------------------------------------

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


  if (versionLabel == "V8") {
    mvaVar[ "SigmaIEtaIEta" ]         = 1;  
    mvaVar[ "DEtaIn" ]                = 1; 
    mvaVar[ "DPhiIn" ]                = 1;  
    mvaVar[ "D0" ]                    = 1;  
    mvaVar[ "DZ" ]                    = 0; 
    mvaVar[ "FBrem" ]                 = 1; 
    mvaVar[ "EOverP" ]                = 1;
    mvaVar[ "ESeedClusterOverPout" ]  = 1;
    mvaVar[ "SigmaIPhiIPhi" ]         = 1;  
    mvaVar[ "NBrem" ]                 = 1;  
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

  //---------------------------------------------------------------
  // specifies the selection applied to events in the training
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
  std::cout << "==> Start TMVAClassificationApplication" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
        std::cout << "Method \"" << regMethod 
                  << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
        for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
          std::cout << it->first << " ";
        }
        std::cout << std::endl;
        return;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------

  const unsigned int nsamples = samples.size();
  
  for( unsigned int i = 0 ; i < nsamples ; ++i ){

    // --- Create the Reader object

    TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

    TMVA::Reader  *fTMVAReader[6];
    for(UInt_t j=0; j<6; ++j) {
      fTMVAReader[j] = new TMVA::Reader( "!Color:!Silent:Error" );  
    }

    // Create a set of variables and declare them to the reader
    // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
    //    Float_t var1, var2;
    //    Float_t var3, var4;
    //    reader->AddVariable( "myvar1 := var1+var2", &var1 );
    //    reader->AddVariable( "myvar2 := var1-var2", &var2 );
    //    reader->AddVariable( "var3",                &var3 );
    //    reader->AddVariable( "var4",                &var4 );

    Float_t                 varEleSigmaIEtaIEta; 
    Float_t                 varEleDEtaIn; 
    Float_t                 varEleDPhiIn; 
    Float_t                 varEleHoverE; 
    Float_t                 varEleD0; 
    Float_t                 varEleDZ; 
    Float_t                 varEleFBrem; 
    Float_t                 varEleEOverP; 
    Float_t                 varEleESeedClusterOverPout; 
    Float_t                 varEleSigmaIPhiIPhi; 
    Float_t                 varEleNBrem; 
    Float_t                 varEleOneOverEMinusOneOverP; 
    Float_t                 varEleESeedClusterOverPIn; 
    Float_t                 varEleIP3d; 
    Float_t                 varEleIP3dSig; 
    Float_t                 varEleStandardLikelihood; 
    Float_t                 varElePFMVA; 
    Float_t                 varEleTMVAKNN; 
    Float_t                 varEleTMVAMLP; 
    Float_t                 varEleTMVAMLPBNN; 
    Float_t                 varEleTMVABDTG; 
    Float_t                 varEleTMVABDT; 
    Float_t                 varEleHcalDepth1OverEcal;
    Float_t                 varEleHcalDepth2OverEcal;
    Float_t                 varEledEtaCalo;
    Float_t                 varEledPhiCalo;
    Float_t                 varElePreShowerOverRaw;
    Float_t                 varEleCovIEtaIPhi;
    Float_t                 varEleSCEtaWidth;
    Float_t                 varEleSCPhiWidth;
    Float_t                 varEleGsfTrackChi2OverNdof;
    Float_t                 varEleR9;
    Float_t                 varEleSeedEMaxOverE;
    Float_t                 varEleSeedETopOverE;
    Float_t                 varEleSeedEBottomOverE;
    Float_t                 varEleSeedELeftOverE;
    Float_t                 varEleSeedERightOverE;
    Float_t                 varEleSeedE2ndOverE;
    Float_t                 varEleSeedE2x5RightOverE;
    Float_t                 varEleSeedE2x5LeftOverE;
    Float_t                 varEleSeedE2x5TopOverE;
    Float_t                 varEleSeedE2x5BottomOverE;
    Float_t                 varEleSeedE2x5MaxOverE;
    Float_t                 varEleSeedE1x3OverE;
    Float_t                 varEleSeedE3x1OverE;
    Float_t                 varEleSeedE1x5OverE;
    Float_t                 varEleSeedE2x2OverE;
    Float_t                 varEleSeedE3x2OverE;
    Float_t                 varEleSeedE3x3OverE;
    Float_t                 varEleSeedE4x4OverE;
    Float_t                 varEleSeedE5x5OverE;
    Float_t                 varEleChargedIso03; 
    Float_t                 varEleNeutralHadronIso03; 
    Float_t                 varEleGammaIso03; 
    Float_t                 varEleChargedIso04; 
    Float_t                 varEleNeutralHadronIso04; 
    Float_t                 varEleGammaIso04; 
    Float_t                 varEleTrkIso03; 
    Float_t                 varEleEMIso03; 
    Float_t                 varEleHadIso03; 
    Float_t                 varEleTrkIso04; 
    Float_t                 varEleEMIso04; 
    Float_t                 varEleHadIso04; 

    if (mvaVar["SigmaIEtaIEta"])           reader->AddVariable( "SigmaIEtaIEta",         &varEleSigmaIEtaIEta         );
    if (mvaVar["DEtaIn"])                  reader->AddVariable( "DEtaIn",                &varEleDEtaIn                );
    if (mvaVar["DPhiIn"])                  reader->AddVariable( "DPhiIn",                &varEleDPhiIn                );
    if (mvaVar["D0"])                      reader->AddVariable( "D0",                    &varEleD0                    );
    if (mvaVar["DZ"])                      reader->AddVariable( "DZ",                    &varEleDZ                    );
    if (mvaVar["FBrem"])                   reader->AddVariable( "FBrem",                 &varEleFBrem                 );
    if (mvaVar["EOverP"])                  reader->AddVariable( "EOverP",                &varEleEOverP                );
    if (mvaVar["ESeedClusterOverPout"])    reader->AddVariable( "ESeedClusterOverPout",  &varEleESeedClusterOverPout  );
    if (mvaVar["SigmaIPhiIPhi"])           reader->AddVariable( "SigmaIPhiIPhi",         &varEleSigmaIPhiIPhi         );
    if (mvaVar["NBrem"])                   reader->AddVariable( "NBrem",                 &varEleNBrem                 );
    if (mvaVar["OneOverEMinusOneOverP"])   reader->AddVariable( "OneOverEMinusOneOverP", &varEleOneOverEMinusOneOverP );
    if (mvaVar["ESeedClusterOverPIn"])     reader->AddVariable( "ESeedClusterOverPIn",   &varEleESeedClusterOverPIn   );
    if (mvaVar["IP3d"])                    reader->AddVariable( "IP3d",                  &varEleIP3d                  );
    if (mvaVar["IP3dSig"])                 reader->AddVariable( "IP3dSig",               &varEleIP3dSig               );
    if (mvaVar["StandardLikelihood"])      reader->AddVariable( "StandardLikelihood",    &varEleStandardLikelihood    );
    if (mvaVar["PFMVA"])                   reader->AddVariable( "PFMVA",                 &varElePFMVA                 );
    if (mvaVar["TMVAKNN"])                 reader->AddVariable( "TMVAKNN",               &varEleTMVAKNN               );
    if (mvaVar["TMVAMLP"])                 reader->AddVariable( "TMVAMLP",               &varEleTMVAMLP               );
    if (mvaVar["TMVAMLPBNN"])              reader->AddVariable( "TMVAMLPBNN",            &varEleTMVAMLPBNN            );
    if (mvaVar["TMVABDTG"])                reader->AddVariable( "TMVABDTG",              &varEleTMVABDTG              );
    if (mvaVar["TMVABDT"])                 reader->AddVariable( "TMVABDT",               &varEleTMVABDT               );
    if (mvaVar["GsfTrackChi2OverNdof"])    reader->AddVariable( "GsfTrackChi2OverNdof",  &varEleGsfTrackChi2OverNdof  );
    if (mvaVar["dEtaCalo"])                reader->AddVariable( "dEtaCalo",              &varEledEtaCalo              );
    if (mvaVar["dPhiCalo"])                reader->AddVariable( "dPhiCalo",              &varEledPhiCalo              );
    if (mvaVar["R9"])                      reader->AddVariable( "R9",                    &varEleR9                    );
    if (mvaVar["SCEtaWidth"])              reader->AddVariable( "SCEtaWidth",            &varEleSCEtaWidth            );
    if (mvaVar["SCPhiWidth"])              reader->AddVariable( "SCPhiWidth",            &varEleSCPhiWidth            );
    if (mvaVar["CovIEtaIPhi"])             reader->AddVariable( "CovIEtaIPhi",           &varEleCovIEtaIPhi           );
    if (Option == 2 || Option == 5) {
      if (mvaVar["PreShowerOverRaw"])      reader->AddVariable( "PreShowerOverRaw",      &varElePreShowerOverRaw      );
    }
    if (mvaVar["HoverE"])                  reader->AddVariable( "HoverE",                &varEleHoverE                );
    if (mvaVar["HcalDepth1OverEcal"])      reader->AddVariable( "HcalDepth1OverEcal",    &varEleHcalDepth1OverEcal    );
    if (Option == 1 || Option == 2 || Option == 4 || Option == 5 ) {
      if (mvaVar["HcalDepth2OverEcal"])    reader->AddVariable( "HcalDepth2OverEcal",    &varEleHcalDepth2OverEcal    );
    }
    if (mvaVar["SeedEMaxOverE"])           reader->AddVariable( "SeedEMaxOverE",         &varEleSeedEMaxOverE         );
    if (mvaVar["SeedETopOverE"])           reader->AddVariable( "SeedETopOverE",         &varEleSeedETopOverE         );
    if (mvaVar["SeedEBottomOverE"])        reader->AddVariable( "SeedEBottomOverE",      &varEleSeedEBottomOverE      );
    if (mvaVar["SeedELeftOverE"])          reader->AddVariable( "SeedELeftOverE",        &varEleSeedELeftOverE        );
    if (mvaVar["SeedERightOverE"])         reader->AddVariable( "SeedERightOverE",       &varEleSeedERightOverE       );
    if (mvaVar["SeedE2ndOverE"])           reader->AddVariable( "SeedE2ndOverE",         &varEleSeedE2ndOverE         );
    if (mvaVar["SeedE2x5RightOverE"])      reader->AddVariable( "SeedE2x5RightOverE",    &varEleSeedE2x5RightOverE    );
    if (mvaVar["SeedE2x5LeftOverE"])       reader->AddVariable( "SeedE2x5LeftOverE",     &varEleSeedE2x5LeftOverE     );
    if (mvaVar["SeedE2x5TopOverE"])        reader->AddVariable( "SeedE2x5TopOverE",      &varEleSeedE2x5TopOverE      );
    if (mvaVar["SeedE2x5BottomOverE"])     reader->AddVariable( "SeedE2x5BottomOverE",   &varEleSeedE2x5BottomOverE   );
    if (mvaVar["SeedE2x5MaxOverE"])        reader->AddVariable( "SeedE2x5MaxOverE",      &varEleSeedE2x5MaxOverE      );
    if (mvaVar["SeedE1x3OverE"])           reader->AddVariable( "SeedE1x3OverE",         &varEleSeedE1x3OverE         );
    if (mvaVar["SeedE3x1OverE"])           reader->AddVariable( "SeedE3x1OverE",         &varEleSeedE3x1OverE         );
    if (mvaVar["SeedE1x5OverE"])           reader->AddVariable( "SeedE1x5OverE",         &varEleSeedE1x5OverE         );
    if (mvaVar["SeedE2x2OverE"])           reader->AddVariable( "SeedE2x2OverE",         &varEleSeedE2x2OverE         );
    if (mvaVar["SeedE3x2OverE"])           reader->AddVariable( "SeedE3x2OverE",         &varEleSeedE3x2OverE         );
    if (mvaVar["SeedE3x3OverE"])           reader->AddVariable( "SeedE3x3OverE",         &varEleSeedE3x3OverE         );
    if (mvaVar["SeedE4x4OverE"])           reader->AddVariable( "SeedE4x4OverE",         &varEleSeedE4x4OverE         );
    if (mvaVar["SeedE5x5OverE"])           reader->AddVariable( "SeedE5x5OverE",         &varEleSeedE5x5OverE         );
    if (mvaVar["ChargedIso03"])            reader->AddVariable( "ChargedIso03",          &varEleChargedIso03          );
    if (mvaVar["NeutralHadronIso03"])      reader->AddVariable( "NeutralHadronIso03",    &varEleNeutralHadronIso03    );
    if (mvaVar["GammaIso03"])              reader->AddVariable( "GammaIso03",            &varEleGammaIso03            );
    if (mvaVar["ChargedIso04"])            reader->AddVariable( "ChargedIso04",          &varEleChargedIso04          );
    if (mvaVar["NeutralHadronIso04"])      reader->AddVariable( "NeutralHadronIso04",    &varEleNeutralHadronIso04    );
    if (mvaVar["GammaIso04"])              reader->AddVariable( "GammaIso04",            &varEleGammaIso04            );
    if (mvaVar["TrkIso03"])                reader->AddVariable( "TrkIso03",              &varEleTrkIso03              );
    if (mvaVar["EMIso03"])                 reader->AddVariable( "EMIso03",               &varEleEMIso03               );
    if (mvaVar["HadIso03"])                reader->AddVariable( "HadIso03",              &varEleHadIso03              );
    if (mvaVar["TrkIso04"])                reader->AddVariable( "TrkIso04",              &varEleTrkIso04              );
    if (mvaVar["EMIso04"])                 reader->AddVariable( "EMIso04",               &varEleEMIso04               );
    if (mvaVar["HadIso04"])                reader->AddVariable( "HadIso04",              &varEleHadIso04              );

    for(UInt_t j=0; j<6; ++j) {
      if (mvaVar["SigmaIEtaIEta"])           fTMVAReader[j]->AddVariable( "SigmaIEtaIEta",         &varEleSigmaIEtaIEta         );
      if (mvaVar["DEtaIn"])                  fTMVAReader[j]->AddVariable( "DEtaIn",                &varEleDEtaIn                );
      if (mvaVar["DPhiIn"])                  fTMVAReader[j]->AddVariable( "DPhiIn",                &varEleDPhiIn                );
      if (mvaVar["D0"])                      fTMVAReader[j]->AddVariable( "D0",                    &varEleD0                    );
      if (mvaVar["DZ"])                      fTMVAReader[j]->AddVariable( "DZ",                    &varEleDZ                    );
      if (mvaVar["FBrem"])                   fTMVAReader[j]->AddVariable( "FBrem",                 &varEleFBrem                 );
      if (mvaVar["EOverP"])                  fTMVAReader[j]->AddVariable( "EOverP",                &varEleEOverP                );
      if (mvaVar["ESeedClusterOverPout"])    fTMVAReader[j]->AddVariable( "ESeedClusterOverPout",  &varEleESeedClusterOverPout  );
      if (mvaVar["SigmaIPhiIPhi"])           fTMVAReader[j]->AddVariable( "SigmaIPhiIPhi",         &varEleSigmaIPhiIPhi         );
      if (mvaVar["NBrem"])                   fTMVAReader[j]->AddVariable( "NBrem",                 &varEleNBrem                 );
      if (mvaVar["OneOverEMinusOneOverP"])   fTMVAReader[j]->AddVariable( "OneOverEMinusOneOverP", &varEleOneOverEMinusOneOverP );
      if (mvaVar["ESeedClusterOverPIn"])     fTMVAReader[j]->AddVariable( "ESeedClusterOverPIn",   &varEleESeedClusterOverPIn   );
      if (mvaVar["IP3d"])                    fTMVAReader[j]->AddVariable( "IP3d",                  &varEleIP3d                  );
      if (mvaVar["IP3dSig"])                 fTMVAReader[j]->AddVariable( "IP3dSig",               &varEleIP3dSig               );
      if (mvaVar["StandardLikelihood"])      fTMVAReader[j]->AddVariable( "StandardLikelihood",    &varEleStandardLikelihood    );
      if (mvaVar["PFMVA"])                   fTMVAReader[j]->AddVariable( "PFMVA",                 &varElePFMVA                 );
      if (mvaVar["TMVAKNN"])                 fTMVAReader[j]->AddVariable( "TMVAKNN",               &varEleTMVAKNN               );
      if (mvaVar["TMVAMLP"])                 fTMVAReader[j]->AddVariable( "TMVAMLP",               &varEleTMVAMLP               );
      if (mvaVar["TMVAMLPBNN"])              fTMVAReader[j]->AddVariable( "TMVAMLPBNN",            &varEleTMVAMLPBNN            );
      if (mvaVar["TMVABDTG"])                fTMVAReader[j]->AddVariable( "TMVABDTG",              &varEleTMVABDTG              );
      if (mvaVar["TMVABDT"])                 fTMVAReader[j]->AddVariable( "TMVABDT",               &varEleTMVABDT               );
      if (mvaVar["GsfTrackChi2OverNdof"])    fTMVAReader[j]->AddVariable( "GsfTrackChi2OverNdof",  &varEleGsfTrackChi2OverNdof  );
      if (mvaVar["dEtaCalo"])                fTMVAReader[j]->AddVariable( "dEtaCalo",              &varEledEtaCalo              );
      if (mvaVar["dPhiCalo"])                fTMVAReader[j]->AddVariable( "dPhiCalo",              &varEledPhiCalo              );
      if (mvaVar["R9"])                      fTMVAReader[j]->AddVariable( "R9",                    &varEleR9                    );
      if (mvaVar["SCEtaWidth"])              fTMVAReader[j]->AddVariable( "SCEtaWidth",            &varEleSCEtaWidth            );
      if (mvaVar["SCPhiWidth"])              fTMVAReader[j]->AddVariable( "SCPhiWidth",            &varEleSCPhiWidth            );
      if (mvaVar["CovIEtaIPhi"])             fTMVAReader[j]->AddVariable( "CovIEtaIPhi",           &varEleCovIEtaIPhi           );
      if (j == 2 || j == 5) {
        if (mvaVar["PreShowerOverRaw"])      fTMVAReader[j]->AddVariable( "PreShowerOverRaw",      &varElePreShowerOverRaw      );
      }
      if (mvaVar["HoverE"])                  fTMVAReader[j]->AddVariable( "HoverE",                &varEleHoverE                );
      if (mvaVar["HcalDepth1OverEcal"])      fTMVAReader[j]->AddVariable( "HcalDepth1OverEcal",    &varEleHcalDepth1OverEcal    );
      if (j == 1 || j == 2 || j == 4 || j == 5 ) {
        if (mvaVar["HcalDepth2OverEcal"])    fTMVAReader[j]->AddVariable( "HcalDepth2OverEcal",    &varEleHcalDepth2OverEcal    );
      }
      if (mvaVar["SeedEMaxOverE"])           fTMVAReader[j]->AddVariable( "SeedEMaxOverE",         &varEleSeedEMaxOverE         );
      if (mvaVar["SeedETopOverE"])           fTMVAReader[j]->AddVariable( "SeedETopOverE",         &varEleSeedETopOverE         );
      if (mvaVar["SeedEBottomOverE"])        fTMVAReader[j]->AddVariable( "SeedEBottomOverE",      &varEleSeedEBottomOverE      );
      if (mvaVar["SeedELeftOverE"])          fTMVAReader[j]->AddVariable( "SeedELeftOverE",        &varEleSeedELeftOverE        );
      if (mvaVar["SeedERightOverE"])         fTMVAReader[j]->AddVariable( "SeedERightOverE",       &varEleSeedERightOverE       );
      if (mvaVar["SeedE2ndOverE"])           fTMVAReader[j]->AddVariable( "SeedE2ndOverE",         &varEleSeedE2ndOverE         );
      if (mvaVar["SeedE2x5RightOverE"])      fTMVAReader[j]->AddVariable( "SeedE2x5RightOverE",    &varEleSeedE2x5RightOverE    );
      if (mvaVar["SeedE2x5LeftOverE"])       fTMVAReader[j]->AddVariable( "SeedE2x5LeftOverE",     &varEleSeedE2x5LeftOverE     );
      if (mvaVar["SeedE2x5TopOverE"])        fTMVAReader[j]->AddVariable( "SeedE2x5TopOverE",      &varEleSeedE2x5TopOverE      );
      if (mvaVar["SeedE2x5BottomOverE"])     fTMVAReader[j]->AddVariable( "SeedE2x5BottomOverE",   &varEleSeedE2x5BottomOverE   );
      if (mvaVar["SeedE2x5MaxOverE"])        fTMVAReader[j]->AddVariable( "SeedE2x5MaxOverE",      &varEleSeedE2x5MaxOverE      );
      if (mvaVar["SeedE1x3OverE"])           fTMVAReader[j]->AddVariable( "SeedE1x3OverE",         &varEleSeedE1x3OverE         );
      if (mvaVar["SeedE3x1OverE"])           fTMVAReader[j]->AddVariable( "SeedE3x1OverE",         &varEleSeedE3x1OverE         );
      if (mvaVar["SeedE1x5OverE"])           fTMVAReader[j]->AddVariable( "SeedE1x5OverE",         &varEleSeedE1x5OverE         );
      if (mvaVar["SeedE2x2OverE"])           fTMVAReader[j]->AddVariable( "SeedE2x2OverE",         &varEleSeedE2x2OverE         );
      if (mvaVar["SeedE3x2OverE"])           fTMVAReader[j]->AddVariable( "SeedE3x2OverE",         &varEleSeedE3x2OverE         );
      if (mvaVar["SeedE3x3OverE"])           fTMVAReader[j]->AddVariable( "SeedE3x3OverE",         &varEleSeedE3x3OverE         );
      if (mvaVar["SeedE4x4OverE"])           fTMVAReader[j]->AddVariable( "SeedE4x4OverE",         &varEleSeedE4x4OverE         );
      if (mvaVar["SeedE5x5OverE"])           fTMVAReader[j]->AddVariable( "SeedE5x5OverE",         &varEleSeedE5x5OverE         );
      if (mvaVar["ChargedIso03"])            fTMVAReader[j]->AddVariable( "ChargedIso03",          &varEleChargedIso03          );
      if (mvaVar["NeutralHadronIso03"])      fTMVAReader[j]->AddVariable( "NeutralHadronIso03",    &varEleNeutralHadronIso03    );
      if (mvaVar["GammaIso03"])              fTMVAReader[j]->AddVariable( "GammaIso03",            &varEleGammaIso03            );
      if (mvaVar["ChargedIso04"])            fTMVAReader[j]->AddVariable( "ChargedIso04",          &varEleChargedIso04          );
      if (mvaVar["NeutralHadronIso04"])      fTMVAReader[j]->AddVariable( "NeutralHadronIso04",    &varEleNeutralHadronIso04    );
      if (mvaVar["GammaIso04"])              fTMVAReader[j]->AddVariable( "GammaIso04",            &varEleGammaIso04            );
      if (mvaVar["TrkIso03"])                fTMVAReader[j]->AddVariable( "TrkIso03",              &varEleTrkIso03              );
      if (mvaVar["EMIso03"])                 fTMVAReader[j]->AddVariable( "EMIso03",               &varEleEMIso03               );
      if (mvaVar["HadIso03"])                fTMVAReader[j]->AddVariable( "HadIso03",              &varEleHadIso03              );
      if (mvaVar["TrkIso04"])                fTMVAReader[j]->AddVariable( "TrkIso04",              &varEleTrkIso04              );
      if (mvaVar["EMIso04"])                 fTMVAReader[j]->AddVariable( "EMIso04",               &varEleEMIso04               );
      if (mvaVar["HadIso04"])                fTMVAReader[j]->AddVariable( "HadIso04",              &varEleHadIso04              );
    }

    // Spectator variables declared in the training have to be added to the reader, too
    //    Float_t spec1,spec2;
    //    reader->AddSpectator( "spec1 := var1*2",   &spec1 );
    //    reader->AddSpectator( "spec2 := var1*3",   &spec2 );

    //Float_t Category_cat1, Category_cat2, Category_cat3;
    if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      //       reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      //       reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      //       reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
    }

    // --- Book the MVA methods

    //--------------------------------------------------------------------------------------
    // tell evaluateMVA_smurf_hww where to find the weights dir, which contains the trained MVA's. 
    // In this example, the weights dir is located at [path]/[dir]
    // and the output root file is written to [path]/[output]
    //--------------------------------------------------------------------------------------

    TString dir    = "weights/";
    TString outdir = "output/";
    TString prefix = label;

    // Book method(s)
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
        TString methodName = TString(it->first) + TString(" method");
        if (!evaluateAllBins) reader->BookMVA( methodName, dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml")); 
        
        if (evaluateAllBins) {
          for(UInt_t j=0; j<6; ++j) {
            TString weightfile;
            if (j==0) weightfile = dir + TString("Subdet0LowPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;
            if (j==1) weightfile = dir + TString("Subdet1LowPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;
            if (j==2) weightfile = dir + TString("Subdet2LowPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;
            if (j==3) weightfile = dir + TString("Subdet0HighPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;
            if (j==4) weightfile = dir + TString("Subdet1HighPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;
            if (j==5) weightfile = dir + TString("Subdet2HighPt_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml") ;

            cout << "j = " << j << " : " << weightfile << endl;            
            fTMVAReader[j]->BookMVA(methodName , weightfile );
          }
        }
      }
    }


    //*****************************************************************************************
    //Prepare Tree Variables
    //*****************************************************************************************
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
    Float_t                 fEleTMVAKNN; 
    Float_t                 fEleTMVAMLP; 
    Float_t                 fEleTMVAMLPBNN; 
    Float_t                 fEleTMVABDTG; 
    Float_t                 fEleTMVABDT; 
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


    //*****************************************************************************************
    //Prepare Output Tree
    //*****************************************************************************************
    TFile *tmpfile = new TFile(inputFile.c_str(), "READ");
    TTree *tmptree = (TTree*)tmpfile->Get("Electrons");


    TFile *EleOutputFile = new TFile(outputFile.c_str(), "RECREATE");
    TTree *EleOutputTree = tmptree->CloneTree(-1, "fast");

    tmpfile->Close();
    delete tmpfile;

    //Add new branches for MVA output
    Float_t fEleBDT         = -9999;
    Float_t fEleBDTG        = -9999;
    Float_t fEleMLPBNN      = -9999;
    Float_t fEleMLP         = -9999;
    Float_t fEleKNN         = -9999;
    Float_t fEleLikelihood  = -9999;
    Float_t fEleCombinedMVA = -9999;

    TBranch* branchEleBDT          = 0;
    TBranch* branchEleBDTG         = 0;
    TBranch* branchEleMLPBNN       = 0;
    TBranch* branchEleMLP          = 0;
    TBranch* branchEleKNN          = 0;
    TBranch* branchEleLikelihood   = 0;
    TBranch* branchEleCombinedMVA  = 0;

    string branchLabel = versionLabel;
    if (buildCommittee) branchLabel = versionLabel + "_C";

    if(Use["BDT"])         branchEleBDT         = EleOutputTree->Branch( ("BDT"+branchLabel).c_str()         , &fEleBDT         , ("BDT"        +branchLabel+"/F").c_str() );
    if(Use["BDTG"])        branchEleBDTG        = EleOutputTree->Branch( ("BDTG"+branchLabel).c_str()        , &fEleBDTG        , ("BDTG"       +branchLabel+"/F").c_str() );
    if(Use["MLPBNN"])      branchEleMLPBNN      = EleOutputTree->Branch( ("MLPBNN"+branchLabel).c_str()      , &fEleMLPBNN      , ("MLPBNN"     +branchLabel+"/F").c_str() );
    if(Use["MLP"])         branchEleMLP         = EleOutputTree->Branch( ("MLP"+branchLabel).c_str()         , &fEleMLP         , ("MLP"        +branchLabel+"/F").c_str() );
    if(Use["KNN"])         branchEleKNN         = EleOutputTree->Branch( ("KNN"+branchLabel).c_str()         , &fEleKNN         , ("KNN"        +branchLabel+"/F").c_str() );
    if(Use["Likelihood"])  branchEleLikelihood  = EleOutputTree->Branch( ("Likelihood"+branchLabel).c_str()  , &fEleLikelihood  , ("Likelihood" +branchLabel+"/F").c_str() );

    if(Use["BDT"])         branchEleBDT         -> SetTitle(("BDT"         + branchLabel + " Output").c_str());
    if(Use["BDTG"])        branchEleBDTG        -> SetTitle(("BDTG"        + branchLabel + " Output").c_str());
    if(Use["MLPBNN"])      branchEleMLPBNN      -> SetTitle(("MLPBNN"      + branchLabel + " Output").c_str());
    if(Use["MLP"])         branchEleMLP         -> SetTitle(("MLP"         + branchLabel + " Output").c_str());
    if(Use["KNN"])         branchEleKNN         -> SetTitle(("KNN"         + branchLabel + " Output").c_str());
    if(Use["Likelihood"])  branchEleLikelihood  -> SetTitle(("Likelihood"  + branchLabel + " Output").c_str());

//     branchEleCombinedMVA  = EleOutputTree->Branch( ("CombinedMVA"+branchLabel).c_str()  , &fEleCombinedMVA  , ("CombinedMVA" +branchLabel+"/F").c_str() );
//     branchEleCombinedMVA  -> SetTitle((string("CombinedMVA")  + branchLabel + " Output").c_str());
 
    //*****************************************************************************************
    //Prepare Input Tree
    //*****************************************************************************************
    TFile *EleFile = new TFile(inputFile.c_str(), "READ");
    TTree *EleTree = (TTree*)EleFile->Get("Electrons");
    EleTree->SetBranchAddress( "weight", &fWeight);
    EleTree->SetBranchAddress( "run", &fRunNumber);
    EleTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
    EleTree->SetBranchAddress( "event", &fEventNumber);
    EleTree->SetBranchAddress( "pt", &fElePt); 
    EleTree->SetBranchAddress( "eta", &fEleEta); 
    EleTree->SetBranchAddress( "phi", &fElePhi); 
    EleTree->SetBranchAddress( "scet", &fEleSCEt); 
    EleTree->SetBranchAddress( "sceta", &fEleSCEta); 
    EleTree->SetBranchAddress( "scphi", &fEleSCPhi); 
    EleTree->SetBranchAddress( "pfiso", &fElePFIso); 
    EleTree->SetBranchAddress( "SigmaIEtaIEta", &fEleSigmaIEtaIEta); 
    EleTree->SetBranchAddress( "DEtaIn", &fEleDEtaIn); 
    EleTree->SetBranchAddress( "DPhiIn", &fEleDPhiIn); 
    EleTree->SetBranchAddress( "HoverE", &fEleHoverE); 
    EleTree->SetBranchAddress( "D0", &fEleD0); 
    EleTree->SetBranchAddress( "DZ", &fEleDZ); 
    EleTree->SetBranchAddress( "FBrem", &fEleFBrem); 
    EleTree->SetBranchAddress( "EOverP", &fEleEOverP); 
    EleTree->SetBranchAddress( "ESeedClusterOverPout", &fEleESeedClusterOverPout); 
    EleTree->SetBranchAddress( "SigmaIPhiIPhi", &fEleSigmaIPhiIPhi); 
    EleTree->SetBranchAddress( "NBrem", &fEleNBrem); 
    EleTree->SetBranchAddress( "OneOverEMinusOneOverP", &fEleOneOverEMinusOneOverP); 
    EleTree->SetBranchAddress( "ESeedClusterOverPIn", &fEleESeedClusterOverPIn); 
    EleTree->SetBranchAddress( "IP3d", &fEleIP3d); 
    EleTree->SetBranchAddress( "IP3dSig", &fEleIP3dSig); 
    EleTree->SetBranchAddress( "StandardLikelihood", &fEleStandardLikelihood); 
    EleTree->SetBranchAddress( "PFMVA", &fElePFMVA); 
    EleTree->SetBranchAddress( ("KNN"+versionLabel).c_str(), &fEleTMVAKNN); 
    EleTree->SetBranchAddress( ("MLP"+versionLabel).c_str(), &fEleTMVAMLP); 
    EleTree->SetBranchAddress( ("MLPBNN"+versionLabel).c_str(), &fEleTMVAMLPBNN); 
    EleTree->SetBranchAddress( ("BDTG"+versionLabel).c_str(), &fEleTMVABDTG); 
    EleTree->SetBranchAddress( ("BDT"+versionLabel).c_str(), &fEleTMVABDT); 
    EleTree->SetBranchAddress( "HcalDepth1OverEcal", &fEleHcalDepth1OverEcal); 
    EleTree->SetBranchAddress( "HcalDepth2OverEcal", &fEleHcalDepth2OverEcal); 
    EleTree->SetBranchAddress( "dEtaCalo", &fEledEtaCalo); 
    EleTree->SetBranchAddress( "dPhiCalo", &fEledPhiCalo); 
    EleTree->SetBranchAddress( "PreShowerOverRaw", &fElePreShowerOverRaw); 
    EleTree->SetBranchAddress( "CovIEtaIPhi", &fEleCovIEtaIPhi); 
    EleTree->SetBranchAddress( "SCEtaWidth", &fEleSCEtaWidth); 
    EleTree->SetBranchAddress( "SCPhiWidth", &fEleSCPhiWidth); 
    EleTree->SetBranchAddress( "GsfTrackChi2OverNdof", &fEleGsfTrackChi2OverNdof); 
    EleTree->SetBranchAddress( "R9", &fEleR9); 
    EleTree->SetBranchAddress( "SeedEMaxOverE", &fEleSeedEMaxOverE); 
    EleTree->SetBranchAddress( "SeedETopOverE", &fEleSeedETopOverE); 
    EleTree->SetBranchAddress( "SeedEBottomOverE", &fEleSeedEBottomOverE); 
    EleTree->SetBranchAddress( "SeedELeftOverE", &fEleSeedELeftOverE); 
    EleTree->SetBranchAddress( "SeedERightOverE", &fEleSeedERightOverE); 
    EleTree->SetBranchAddress( "SeedE2ndOverE", &fEleSeedE2ndOverE); 
    EleTree->SetBranchAddress( "SeedE2x5RightOverE", &fEleSeedE2x5RightOverE); 
    EleTree->SetBranchAddress( "SeedE2x5LeftOverE", &fEleSeedE2x5LeftOverE); 
    EleTree->SetBranchAddress( "SeedE2x5TopOverE", &fEleSeedE2x5TopOverE); 
    EleTree->SetBranchAddress( "SeedE2x5BottomOverE", &fEleSeedE2x5BottomOverE); 
    EleTree->SetBranchAddress( "SeedE2x5MaxOverE", &fEleSeedE2x5MaxOverE); 
    EleTree->SetBranchAddress( "SeedE1x3OverE", &fEleSeedE1x3OverE); 
    EleTree->SetBranchAddress( "SeedE3x1OverE", &fEleSeedE3x1OverE); 
    EleTree->SetBranchAddress( "SeedE1x5OverE", &fEleSeedE1x5OverE); 
    EleTree->SetBranchAddress( "SeedE2x2OverE", &fEleSeedE2x2OverE); 
    EleTree->SetBranchAddress( "SeedE3x2OverE", &fEleSeedE3x2OverE); 
    EleTree->SetBranchAddress( "SeedE3x3OverE", &fEleSeedE3x3OverE); 
    EleTree->SetBranchAddress( "SeedE4x4OverE", &fEleSeedE4x4OverE); 
    EleTree->SetBranchAddress( "SeedE5x5OverE", &fEleSeedE5x5OverE); 
    EleTree->SetBranchAddress( "ChargedIso03", &fEleChargedIso03); 
    EleTree->SetBranchAddress( "NeutralHadronIso03", &fEleNeutralHadronIso03); 
    EleTree->SetBranchAddress( "GammaIso03", &fEleGammaIso03); 
    EleTree->SetBranchAddress( "ChargedIso04", &fEleChargedIso04); 
    EleTree->SetBranchAddress( "NeutralHadronIso04", &fEleNeutralHadronIso04); 
    EleTree->SetBranchAddress( "GammaIso04", &fEleGammaIso04); 
    EleTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fEleChargedIso04FromOtherVertices); 
    EleTree->SetBranchAddress( "NeutralHadronIso04_10Threshold", &fEleNeutralHadronIso04_10Threshold); 
    EleTree->SetBranchAddress( "GammaIso04_10Threshold", &fEleGammaIso04_10Threshold); 
    EleTree->SetBranchAddress( "TrkIso03", &fEleTrkIso03); 
    EleTree->SetBranchAddress( "EMIso03", &fEleEMIso03); 
    EleTree->SetBranchAddress( "HadIso03", &fEleHadIso03); 
    EleTree->SetBranchAddress( "TrkIso04", &fEleTrkIso04); 
    EleTree->SetBranchAddress( "EMIso04", &fEleEMIso04); 
    EleTree->SetBranchAddress( "HadIso04", &fEleHadIso04); 
    EleTree->SetBranchAddress( "Rho", &fRho); 
    EleTree->SetBranchAddress( "NVertices", &fNVertices); 


    // --- Event loop

    // Prepare the event tree
    // - here the variable names have to corresponds to your tree
    // - you can use the same variables as above which is slightly faster,
    //   but of course you can use different ones and copy the values inside the event loop
    //
  
    // Efficiency calculator for cut method
    Int_t    nSelCutsGA = 0;
    Double_t effS       = 0.7;

    std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

    std::cout << "--- Processing: " << EleTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    int npass   = 0;
    float yield = 0.;

   for (Long64_t ievt=0; ievt<EleTree->GetEntries();ievt++) {

      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      EleTree->GetEntry(ievt);

      Double_t rho = 0;
      if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

     //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------
      varEleSigmaIEtaIEta         = fEleSigmaIEtaIEta; 
      varEleDEtaIn                = fEleDEtaIn; 
      varEleDPhiIn                = fEleDPhiIn; 
      varEleD0                    = fEleD0; 
      varEleDZ                    = fEleDZ; 
      varEleFBrem                 = fEleFBrem; 
      varEleEOverP                = fEleEOverP; 
      varEleESeedClusterOverPout  = fEleESeedClusterOverPout; 
      varEleSigmaIPhiIPhi         = fEleSigmaIPhiIPhi; 
      varEleNBrem                 = fEleNBrem; 
      varEleOneOverEMinusOneOverP = fEleOneOverEMinusOneOverP; 
      varEleESeedClusterOverPIn   = fEleESeedClusterOverPIn; 
      varEleIP3d                  = fEleIP3d; 
      varEleIP3dSig               = fEleIP3dSig; 
      varEleStandardLikelihood    = fEleStandardLikelihood; 
      varElePFMVA                 = fElePFMVA; 
      varEleTMVAKNN               = fEleTMVAKNN; 
      varEleTMVAMLP               = fEleTMVAMLP; 
      varEleTMVAMLPBNN            = fEleTMVAMLPBNN; 
      varEleTMVABDTG              = fEleTMVABDTG; 
      varEleTMVABDT               = fEleTMVABDT; 
      varEleGsfTrackChi2OverNdof  = fEleGsfTrackChi2OverNdof;
      varEledEtaCalo              = fEledEtaCalo;
      varEledPhiCalo              = fEledPhiCalo;
      varEleR9                    = fEleR9;
      varEleSCEtaWidth            = fEleSCEtaWidth;
      varEleSCPhiWidth            = fEleSCPhiWidth;
      varEleCovIEtaIPhi           = fEleCovIEtaIPhi;
      varElePreShowerOverRaw      = fElePreShowerOverRaw;
      varEleHoverE                = fEleHoverE - rho*ElectronEffectiveArea(kEleHoverE,fEleEta);
      varEleHcalDepth1OverEcal    = fEleHcalDepth1OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth1OverEcal,fEleEta);
      varEleHcalDepth2OverEcal    = fEleHcalDepth2OverEcal - rho*ElectronEffectiveArea(kEleHcalDepth2OverEcal,fEleEta);
      varEleSeedEMaxOverE         = fEleSeedEMaxOverE;
      varEleSeedETopOverE         = fEleSeedETopOverE;
      varEleSeedEBottomOverE      = fEleSeedEBottomOverE;
      varEleSeedELeftOverE        = fEleSeedELeftOverE;
      varEleSeedERightOverE       = fEleSeedERightOverE;
      varEleSeedE2ndOverE         = fEleSeedE2ndOverE;
      varEleSeedE2x5RightOverE    = fEleSeedE2x5RightOverE;
      varEleSeedE2x5LeftOverE     = fEleSeedE2x5LeftOverE;
      varEleSeedE2x5TopOverE      = fEleSeedE2x5TopOverE;
      varEleSeedE2x5BottomOverE   = fEleSeedE2x5BottomOverE;
      varEleSeedE2x5MaxOverE      = fEleSeedE2x5MaxOverE;
      varEleSeedE1x3OverE         = fEleSeedE1x3OverE;
      varEleSeedE3x1OverE         = fEleSeedE3x1OverE;
      varEleSeedE1x5OverE         = fEleSeedE1x5OverE;
      varEleSeedE2x2OverE         = fEleSeedE2x2OverE;
      varEleSeedE3x2OverE         = fEleSeedE3x2OverE;
      varEleSeedE3x3OverE         = fEleSeedE3x3OverE;
      varEleSeedE4x4OverE         = fEleSeedE4x4OverE;
      varEleSeedE5x5OverE         = fEleSeedE5x5OverE;
      if (versionLabel == "V15") {
        varEleChargedIso03          = (fEleChargedIso03 - rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta))/fElePt;
        varEleNeutralHadronIso03    = (fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
        varEleGammaIso03            = (fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta))/fElePt;
        varEleChargedIso04          = (fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta))/fElePt;
        varEleNeutralHadronIso04    = (fEleNeutralHadronIso04 - rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta))/fElePt;
        varEleGammaIso04            = (fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta))/fElePt;
      } else {
        varEleChargedIso03          = fEleChargedIso03 - rho*ElectronEffectiveArea(kEleChargedIso03,fEleEta);
        varEleNeutralHadronIso03    = fEleNeutralHadronIso03 - rho*ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
        varEleGammaIso03            = fEleGammaIso03 - rho*ElectronEffectiveArea(kEleGammaIso03,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta);
        varEleChargedIso04          = fEleChargedIso04 - rho*ElectronEffectiveArea(kEleChargedIso04,fEleEta);
        varEleNeutralHadronIso04    = fEleNeutralHadronIso04 - rho*ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta) + rho*ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta);
        varEleGammaIso04            = fEleGammaIso04 - rho*ElectronEffectiveArea(kEleGammaIso04,fEleEta) + rho*ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta);
      }
      varEleTrkIso03              = (fEleTrkIso03 - rho*ElectronEffectiveArea(kEleTrkIso03,fEleEta))/fElePt;
      varEleEMIso03               = (fEleEMIso03 - rho*ElectronEffectiveArea(kEleEMIso03,fEleEta))/fElePt;
      varEleHadIso03              = (fEleHadIso03 - rho*ElectronEffectiveArea(kEleHadIso03,fEleEta))/fElePt;
      varEleTrkIso04              = (fEleTrkIso04 - rho*ElectronEffectiveArea(kEleTrkIso04,fEleEta))/fElePt;
      varEleEMIso04               = (fEleEMIso04 - rho*ElectronEffectiveArea(kEleEMIso04,fEleEta))/fElePt;
      varEleHadIso04              = (fEleHadIso04 - rho*ElectronEffectiveArea(kEleHadIso04,fEleEta))/fElePt;

      // --- Return the MVA outputs and weights
      Int_t subdet = 0;
      if (fabs(fEleSCEta) < 1.0) subdet = 0;
      else if (fabs(fEleSCEta) < 1.479) subdet = 1;
      else subdet = 2;
      Int_t ptBin = 0;
      if (fElePt > 20.0) ptBin = 1;
      
      Int_t MVABin = -1;
      if (subdet == 0 && ptBin == 0) MVABin = 0;
      if (subdet == 1 && ptBin == 0) MVABin = 1;
      if (subdet == 2 && ptBin == 0) MVABin = 2;
      if (subdet == 0 && ptBin == 1) MVABin = 3;
      if (subdet == 1 && ptBin == 1) MVABin = 4;
      if (subdet == 2 && ptBin == 1) MVABin = 5;
      assert(MVABin >= 0 && MVABin <= 5);

      TMVA::Reader  *tmpReader = reader;
      if (evaluateAllBins) tmpReader = fTMVAReader[MVABin];      
      
      if (Use["BDT"]){
        fEleBDT = tmpReader->EvaluateMVA( "BDT method" );
        branchEleBDT->Fill();
      }
      if (Use["BDTG"]){
        fEleBDTG  = tmpReader->EvaluateMVA( "BDTG method" );
        branchEleBDTG->Fill();

//         cout << fEventNumber << " " << MVABin << " : " << fElePt << " " << fEleBDTG << endl;

      }
      if (Use["MLPBNN"]){
        fEleMLPBNN  = tmpReader->EvaluateMVA( "MLPBNN method" );
        branchEleMLPBNN->Fill();
      }
      if (Use["MLP"]){
        fEleMLP  = tmpReader->EvaluateMVA( "MLP method" );
        branchEleMLP->Fill();
      }
      if (Use["KNN"]){
        fEleKNN  = tmpReader->EvaluateMVA( "KNN method" );
        branchEleKNN->Fill();
      }
       if (Use["Likelihood"]){
        Double_t LH = tmpReader->EvaluateMVA( "Likelihood method" );
        fEleLikelihood = LH;
//         double newLik = 0.0;
//         if     (LH<=0) newLik = -20.0;
//         else if(LH>=1) newLik =  20.0;
//         else                   newLik = log(LH/(1.0-LH));        
//         fEleLikelihood = newLik;
        branchEleLikelihood->Fill();
      }

//        //Combined MVA
//       fEleCombinedMVA  = tmpReader->GetProba( "BDTG method" ) 
//         * tmpReader->GetProba( "BDT method" ) 
//         * tmpReader->GetProba( "MLP method" ) 
//         * tmpReader->GetProba( "MLPBNN method" )  
//         * tmpReader->GetProba( "KNN method" );
//       branchEleCombinedMVA->Fill();
 


//        cout << "Test Event : " << fRunNumber << " " << fEventNumber << " : " << fElePt << " " << fEleEta << " " << fElePhi << " : ";
//        cout << varEleSigmaIEtaIEta         << " "
//             << varEleDEtaIn                << " "
//             << varEleDPhiIn                << " "
//             << varEleD0                    << " "
//             << varEleDZ                    << " "
//             << varEleFBrem                 << " "
//             << varEleEOverP                << " "
//             << varEleESeedClusterOverPout  << " "
//             << varEleSigmaIPhiIPhi         << " "
//             << varEleNBrem                 << " "
//             << varEleOneOverEMinusOneOverP << " "
//             << varEleESeedClusterOverPIn   << " "
//             << varEleIP3d                  << " "
//             << varEleIP3dSig               << " "
//             << varEleGsfTrackChi2OverNdof  << " "
//             << varEledEtaCalo              << " "
//             << varEledPhiCalo              << " "
//             << varEleR9                    << " "
//             << varEleSCEtaWidth            << " "
//             << varEleSCPhiWidth            << " "
//             << varEleCovIEtaIPhi           << " "
//             << varElePreShowerOverRaw      << " "

//             << "   :   " 
//             << fEleChargedIso03 << " "
//             << fEleNeutralHadronIso03 << " "
//             << fEleGammaIso03  << " "
//             << fEleChargedIso04 << " "
//             << fEleNeutralHadronIso04  << " "
//             << fEleGammaIso04  << " "
//             << "   -   "         
//             << ElectronEffectiveArea(kEleChargedIso03,fEleEta) << " "
//             <<  ElectronEffectiveArea(kEleNeutralHadronIso03,fEleEta) << " "
//             << ElectronEffectiveArea(kEleGammaIso03,fEleEta)  << " "
//             <<  ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip03,fEleEta) << " "
//             <<  ElectronEffectiveArea(kEleChargedIso04,fEleEta) << " "
//             << ElectronEffectiveArea(kEleNeutralHadronIso04,fEleEta)  << " "
//             <<    ElectronEffectiveArea(kEleNeutralHadronIso007,fEleEta) << " "
//             <<  ElectronEffectiveArea(kEleGammaIso04,fEleEta) << " "
//             <<    ElectronEffectiveArea(kEleGammaIsoVetoEtaStrip04,fEleEta) << " "
//             << "   :   " 
//             << varEleChargedIso03          << " "
//             << varEleNeutralHadronIso03    << " "
//             << varEleGammaIso03            << " "
//             << varEleChargedIso04          << " "
//             << varEleNeutralHadronIso04    << " "
//             << varEleGammaIso04            << " "
//             << "   :   Bin"
//             << MVABin << " : " 
//             << fEleBDTG
//             << endl;


    } // End main loop


    std::cout << npass << " events passing selection, yield " << yield << std::endl;
 
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    //Write Output file
    EleOutputFile->Write();
    EleOutputFile->Close();

    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done with sample " << samples.at(i) << endl << std::endl;
  } 
}
