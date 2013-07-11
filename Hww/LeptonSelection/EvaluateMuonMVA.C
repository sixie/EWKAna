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

void EvaluateMuonMVA(
string inputFile      = "MuonSelectionTraining.Fake.weighted.Subdet0LowPt.root", 
string outputFile     = "output/MuonNtuple.Fake.Subdet0LowPt.root",
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
  //Default Variables to Use
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
    mvaVar[ "ChargedIso03OverPt" ]   = 0;  
    mvaVar[ "NeutralIso03OverPt" ]   = 0;  
    mvaVar[ "ChargedIso04OverPt" ]   = 0;  
    mvaVar[ "NeutralIso04OverPt" ]   = 0;  
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

    Float_t                 varMuTkNchi2; 
    Float_t                 varMuGlobalNchi2; 
    Float_t                 varMuNValidHits; 
    Float_t                 varMuNTrackerHits; 
    Float_t                 varMuNPixelHits; 
    Float_t                 varMuNMatches; 
    Float_t                 varMuD0; 
    Float_t                 varMuIP3d; 
    Float_t                 varMuIP3dSig; 
    Float_t                 varMuTrkKink; 
    Float_t                 varMuGlobalKink; 
    Float_t                 varMuSegmentCompatibility; 
    Float_t                 varMuCaloCompatibility; 
    Float_t                 varMuHadEnergyOverPt; 
    Float_t                 varMuHoEnergyOverPt; 
    Float_t                 varMuEmEnergyOverPt; 
    Float_t                 varMuHadS9EnergyOverPt; 
    Float_t                 varMuHoS9EnergyOverPt; 
    Float_t                 varMuEmS9EnergyOverPt; 
    Float_t                 varMuChargedIso03OverPt; 
    Float_t                 varMuNeutralIso03OverPt; 
    Float_t                 varMuChargedIso04OverPt; 
    Float_t                 varMuNeutralIso04OverPt; 
    Float_t                 varMuTrkIso03OverPt; 
    Float_t                 varMuEMIso03OverPt; 
    Float_t                 varMuHadIso03OverPt; 
    Float_t                 varMuTrkIso05OverPt; 
    Float_t                 varMuEMIso05OverPt; 
    Float_t                 varMuHadIso05OverPt; 
    Float_t                 varMuNVtx; 
   
    if (mvaVar["TkNchi2"])               reader->AddVariable( "TkNchi2",              &varMuTkNchi2              );
    if (mvaVar["GlobalNchi2"])           reader->AddVariable( "GlobalNchi2",          &varMuGlobalNchi2          );
    if (mvaVar["NValidHits"])            reader->AddVariable( "NValidHits",           &varMuNValidHits           );
    if (mvaVar["NTrackerHits"])          reader->AddVariable( "NTrackerHits",         &varMuNTrackerHits         );
    if (mvaVar["NPixelHits"])            reader->AddVariable( "NPixelHits",           &varMuNPixelHits           );
    if (mvaVar["NMatches"])              reader->AddVariable( "NMatches",             &varMuNMatches             );
    if (mvaVar["D0"])                    reader->AddVariable( "D0",                   &varMuD0                   );
    if (mvaVar["IP3d"])                  reader->AddVariable( "IP3d",                 &varMuIP3d                 );
    if (mvaVar["IP3dSig"])               reader->AddVariable( "IP3dSig",              &varMuIP3dSig              );
    if (mvaVar["TrkKink"])               reader->AddVariable( "TrkKink",              &varMuTrkKink              );
    if (mvaVar["GlobalKink"])            reader->AddVariable( "GlobalKink",           &varMuGlobalKink           );
    if (mvaVar["SegmentCompatibility"])  reader->AddVariable( "SegmentCompatibility", &varMuSegmentCompatibility );
    if (mvaVar["CaloCompatibility"])     reader->AddVariable( "CaloCompatibility",    &varMuCaloCompatibility    );
    if (mvaVar["HadEnergyOverPt"])       reader->AddVariable( "HadEnergyOverPt",      &varMuHadEnergyOverPt      );
    if (mvaVar["HoEnergyOverPt"])        reader->AddVariable( "HoEnergyOverPt",       &varMuHoEnergyOverPt       );
    if (mvaVar["EmEnergyOverPt"])        reader->AddVariable( "EmEnergyOverPt",       &varMuEmEnergyOverPt       );
    if (mvaVar["HadS9EnergyOverPt"])     reader->AddVariable( "HadS9EnergyOverPt",    &varMuHadS9EnergyOverPt    );
    if (mvaVar["HoS9EnergyOverPt"])      reader->AddVariable( "HoS9EnergyOverPt",     &varMuHoS9EnergyOverPt     );
    if (mvaVar["EmS9EnergyOverPt"])      reader->AddVariable( "EmS9EnergyOverPt",     &varMuEmS9EnergyOverPt     );
    if (mvaVar["ChargedIso03OverPt"])    reader->AddVariable( "ChargedIso03OverPt",   &varMuChargedIso03OverPt   );
    if (mvaVar["NeutralIso03OverPt"])    reader->AddVariable( "NeutralIso03OverPt",   &varMuNeutralIso03OverPt   );
    if (mvaVar["ChargedIso04OverPt"])    reader->AddVariable( "ChargedIso04OverPt",   &varMuChargedIso04OverPt   );
    if (mvaVar["NeutralIso04OverPt"])    reader->AddVariable( "NeutralIso04OverPt",   &varMuNeutralIso04OverPt   );
    if (mvaVar["TrkIso03OverPt"])        reader->AddVariable( "TrkIso03OverPt",       &varMuTrkIso03OverPt       );
    if (mvaVar["EMIso03OverPt"])         reader->AddVariable( "EMIso03OverPt",        &varMuEMIso03OverPt        );
    if (mvaVar["HadIso03OverPt"])        reader->AddVariable( "HadIso03OverPt",       &varMuHadIso03OverPt       );
    if (mvaVar["TrkIso05OverPt"])        reader->AddVariable( "TrkIso05OverPt",       &varMuTrkIso05OverPt       );
    if (mvaVar["EMIso05OverPt"])         reader->AddVariable( "EMIso05OverPt",        &varMuEMIso05OverPt        );
    if (mvaVar["HadIso05OverPt"])        reader->AddVariable( "HadIso05OverPt",       &varMuHadIso05OverPt       );
    if (mvaVar["NVtx"])                  reader->AddVariable( "NVtx",                 &varMuNVtx                 );

    for(UInt_t j=0; j<6; ++j) {
      if (mvaVar["TkNchi2"])               fTMVAReader[j]->AddVariable( "TkNchi2",              &varMuTkNchi2              );
      if (mvaVar["GlobalNchi2"])           fTMVAReader[j]->AddVariable( "GlobalNchi2",          &varMuGlobalNchi2          );
      if (mvaVar["NValidHits"])            fTMVAReader[j]->AddVariable( "NValidHits",           &varMuNValidHits           );
      if (mvaVar["NTrackerHits"])          fTMVAReader[j]->AddVariable( "NTrackerHits",         &varMuNTrackerHits         );
      if (mvaVar["NPixelHits"])            fTMVAReader[j]->AddVariable( "NPixelHits",           &varMuNPixelHits           );
      if (mvaVar["NMatches"])              fTMVAReader[j]->AddVariable( "NMatches",             &varMuNMatches             );
      if (mvaVar["D0"])                    fTMVAReader[j]->AddVariable( "D0",                   &varMuD0                   );
      if (mvaVar["IP3d"])                  fTMVAReader[j]->AddVariable( "IP3d",                 &varMuIP3d                 );
      if (mvaVar["IP3dSig"])               fTMVAReader[j]->AddVariable( "IP3dSig",              &varMuIP3dSig              );
      if (mvaVar["TrkKink"])               fTMVAReader[j]->AddVariable( "TrkKink",              &varMuTrkKink              );
      if (mvaVar["GlobalKink"])            fTMVAReader[j]->AddVariable( "GlobalKink",           &varMuGlobalKink           );
      if (mvaVar["SegmentCompatibility"])  fTMVAReader[j]->AddVariable( "SegmentCompatibility", &varMuSegmentCompatibility );
      if (mvaVar["CaloCompatibility"])     fTMVAReader[j]->AddVariable( "CaloCompatibility",    &varMuCaloCompatibility    );
      if (mvaVar["HadEnergyOverPt"])       fTMVAReader[j]->AddVariable( "HadEnergyOverPt",      &varMuHadEnergyOverPt      );
      if (mvaVar["HoEnergyOverPt"])        fTMVAReader[j]->AddVariable( "HoEnergyOverPt",       &varMuHoEnergyOverPt       );
      if (mvaVar["EmEnergyOverPt"])        fTMVAReader[j]->AddVariable( "EmEnergyOverPt",       &varMuEmEnergyOverPt       );
      if (mvaVar["HadS9EnergyOverPt"])     fTMVAReader[j]->AddVariable( "HadS9EnergyOverPt",    &varMuHadS9EnergyOverPt    );
      if (mvaVar["HoEnergyOverPt"])        fTMVAReader[j]->AddVariable( "HoS9EnergyOverPt",     &varMuHoS9EnergyOverPt     );
      if (mvaVar["EmS9EnergyOverPt"])      fTMVAReader[j]->AddVariable( "EmS9EnergyOverPt",     &varMuEmS9EnergyOverPt     );
      if (mvaVar["ChargedIso03OverPt"])    fTMVAReader[j]->AddVariable( "ChargedIso03OverPt",   &varMuChargedIso03OverPt   );
      if (mvaVar["NeutralIso03OverPt"])    fTMVAReader[j]->AddVariable( "NeutralIso03OverPt",   &varMuNeutralIso03OverPt   );
      if (mvaVar["ChargedIso04OverPt"])    fTMVAReader[j]->AddVariable( "ChargedIso04OverPt",   &varMuChargedIso04OverPt   );
      if (mvaVar["NeutralIso04OverPt"])    fTMVAReader[j]->AddVariable( "NeutralIso04OverPt",   &varMuNeutralIso04OverPt   );
      if (mvaVar["TrkIso03OverPt"])        fTMVAReader[j]->AddVariable( "TrkIso03OverPt",       &varMuTrkIso03OverPt       );
      if (mvaVar["EMIso03OverPt"])         fTMVAReader[j]->AddVariable( "EMIso03OverPt",        &varMuEMIso03OverPt        );
      if (mvaVar["HadIso03OverPt"])        fTMVAReader[j]->AddVariable( "HadIso03OverPt",       &varMuHadIso03OverPt       );
      if (mvaVar["TrkIso05OverPt"])        fTMVAReader[j]->AddVariable( "TrkIso05OverPt",       &varMuTrkIso05OverPt       );
      if (mvaVar["EMIso05OverPt"])         fTMVAReader[j]->AddVariable( "EMIso05OverPt",        &varMuEMIso05OverPt        );
      if (mvaVar["HadIso05OverPt"])        fTMVAReader[j]->AddVariable( "HadIso05OverPt",       &varMuHadIso05OverPt       );
      if (mvaVar["NVtx"])                  fTMVAReader[j]->AddVariable( "NVtx",                 &varMuNVtx                 );
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
            if (i==0) weightfile = dir + TString("BarrelPtBin0") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            if (j==1) weightfile = dir + TString("EndcapPtBin0") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            if (j==2) weightfile = dir + TString("BarrelPtBin1") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            if (j==3) weightfile = dir + TString("EndcapPtBin1") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            if (j==4) weightfile = dir + TString("BarrelPtBin2") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
            if (j==5) weightfile = dir + TString("EndcapPtBin2") + TString("_") + prefix + TString("_") + TString(it->first) + TString(".weights.xml");

            cout << "j = " << j << " : " << weightfile << endl;            
            fTMVAReader[j]->BookMVA(methodName , weightfile );
 
          }
        }
      }
    }


    //*****************************************************************************************
    //Prepare Tree Variables
    //*****************************************************************************************
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


    //*****************************************************************************************
    //Prepare Output Tree
    //*****************************************************************************************
    TFile *tmpfile = new TFile(inputFile.c_str(), "READ");
    TTree *tmptree = (TTree*)tmpfile->Get("Muons");


    TFile *MuOutputFile = new TFile(outputFile.c_str(), "RECREATE");
    TTree *MuOutputTree = tmptree->CloneTree(-1, "fast");

    tmpfile->Close();
    delete tmpfile;

    //Add new branches for MVA output
    Float_t fMuBDT         = -9999;
    Float_t fMuBDTG        = -9999;
    Float_t fMuMLPBNN      = -9999;
    Float_t fMuMLP         = -9999;
    Float_t fMuKNN         = -9999;
    Float_t fMuLikelihood  = -9999;
    Float_t fMuCombinedMVA = -9999;

    TBranch* branchMuBDT          = 0;
    TBranch* branchMuBDTG         = 0;
    TBranch* branchMuMLPBNN       = 0;
    TBranch* branchMuMLP          = 0;
    TBranch* branchMuKNN          = 0;
    TBranch* branchMuLikelihood   = 0;
    TBranch* branchMuCombinedMVA  = 0;

    string branchLabel = versionLabel;
    if (buildCommittee) branchLabel = versionLabel + "_C";

    if(Use["BDT"])         branchMuBDT         = MuOutputTree->Branch( ("BDT"+branchLabel).c_str()         , &fMuBDT         , ("BDT"        +branchLabel+"/F").c_str() );
    if(Use["BDTG"])        branchMuBDTG        = MuOutputTree->Branch( ("BDTG"+branchLabel).c_str()        , &fMuBDTG        , ("BDTG"       +branchLabel+"/F").c_str() );
    if(Use["MLPBNN"])      branchMuMLPBNN      = MuOutputTree->Branch( ("MLPBNN"+branchLabel).c_str()      , &fMuMLPBNN      , ("MLPBNN"     +branchLabel+"/F").c_str() );
    if(Use["MLP"])         branchMuMLP         = MuOutputTree->Branch( ("MLP"+branchLabel).c_str()         , &fMuMLP         , ("MLP"        +branchLabel+"/F").c_str() );
    if(Use["KNN"])         branchMuKNN         = MuOutputTree->Branch( ("KNN"+branchLabel).c_str()         , &fMuKNN         , ("KNN"        +branchLabel+"/F").c_str() );
    if(Use["Likelihood"])  branchMuLikelihood  = MuOutputTree->Branch( ("Likelihood"+branchLabel).c_str()  , &fMuLikelihood  , ("Likelihood" +branchLabel+"/F").c_str() );

    if(Use["BDT"])         branchMuBDT         -> SetTitle(("BDT"         + branchLabel + " Output").c_str());
    if(Use["BDTG"])        branchMuBDTG        -> SetTitle(("BDTG"        + branchLabel + " Output").c_str());
    if(Use["MLPBNN"])      branchMuMLPBNN      -> SetTitle(("MLPBNN"      + branchLabel + " Output").c_str());
    if(Use["MLP"])         branchMuMLP         -> SetTitle(("MLP"         + branchLabel + " Output").c_str());
    if(Use["KNN"])         branchMuKNN         -> SetTitle(("KNN"         + branchLabel + " Output").c_str());
    if(Use["Likelihood"])  branchMuLikelihood  -> SetTitle(("Likelihood"  + branchLabel + " Output").c_str());

 
    //*****************************************************************************************
    //Prepare Input Tree
    //*****************************************************************************************
    TFile *MuFile = new TFile(inputFile.c_str(), "READ");
    TTree *MuTree = (TTree*)MuFile->Get("Muons");
    MuTree->SetBranchAddress( "weight", &fWeight);
    MuTree->SetBranchAddress( "run", &fRunNumber);
    MuTree->SetBranchAddress( "lumi", &fLumiSectionNumber);
    MuTree->SetBranchAddress( "event", &fEventNumber);
    MuTree->SetBranchAddress( "pt", &fMuPt); 
    MuTree->SetBranchAddress( "eta", &fMuEta); 
    MuTree->SetBranchAddress( "phi", &fMuPhi); 
    MuTree->SetBranchAddress( "pfiso", &fMuPFIso); 
    MuTree->SetBranchAddress( "TkNchi2", &fMuTkNchi2); 
    MuTree->SetBranchAddress( "GlobalNchi2", &fMuGlobalNchi2); 
    MuTree->SetBranchAddress( "NValidHits", &fMuNValidHits); 
    MuTree->SetBranchAddress( "NTrackerHits", &fMuNTrackerHits); 
    MuTree->SetBranchAddress( "NPixelHits", &fMuNPixelHits); 
    MuTree->SetBranchAddress( "NMatches", &fMuNMatches); 
    MuTree->SetBranchAddress( "D0", &fMuD0); 
    MuTree->SetBranchAddress( "IP3d", &fMuIP3d); 
    MuTree->SetBranchAddress( "IP3dSig", &fMuIP3dSig); 
    MuTree->SetBranchAddress( "TrkKink", &fMuTrkKink); 
    MuTree->SetBranchAddress( "GlobalKink", &fMuGlobalKink); 
    MuTree->SetBranchAddress( "SegmentCompatibility", &fMuSegmentCompatibility); 
    MuTree->SetBranchAddress( "CaloCompatibility", &fMuCaloCompatibility); 
    MuTree->SetBranchAddress( "HadEnergy", &fMuHadEnergy); 
    MuTree->SetBranchAddress( "HoEnergy", &fMuHoEnergy); 
    MuTree->SetBranchAddress( "EmEnergy", &fMuEmEnergy); 
    MuTree->SetBranchAddress( "HadS9Energy", &fMuHadS9Energy); 
    MuTree->SetBranchAddress( "HoS9Energy", &fMuHoS9Energy); 
    MuTree->SetBranchAddress( "EmS9Energy", &fMuEmS9Energy); 
    MuTree->SetBranchAddress( "ChargedIso03", &fMuChargedIso03); 
    MuTree->SetBranchAddress( "ChargedIso03FromOtherVertices", &fMuChargedIso03FromOtherVertices); 
    MuTree->SetBranchAddress( "NeutralIso03_05Threshold", &fMuNeutralIso03_05Threshold); 
    MuTree->SetBranchAddress( "NeutralIso03_10Threshold", &fMuNeutralIso03_10Threshold); 
    MuTree->SetBranchAddress( "ChargedIso04", &fMuChargedIso04); 
    MuTree->SetBranchAddress( "ChargedIso04FromOtherVertices", &fMuChargedIso04FromOtherVertices); 
    MuTree->SetBranchAddress( "NeutralIso04_05Threshold", &fMuNeutralIso04_05Threshold); 
    MuTree->SetBranchAddress( "NeutralIso04_10Threshold", &fMuNeutralIso04_10Threshold); 
    MuTree->SetBranchAddress( "TrkIso03", &fMuTrkIso03); 
    MuTree->SetBranchAddress( "EMIso03", &fMuEMIso03); 
    MuTree->SetBranchAddress( "HadIso03", &fMuHadIso03); 
    MuTree->SetBranchAddress( "TrkIso05", &fMuTrkIso05); 
    MuTree->SetBranchAddress( "EMIso05", &fMuEMIso05); 
    MuTree->SetBranchAddress( "HadIso05", &fMuHadIso05); 
    MuTree->SetBranchAddress( "Rho", &fRho); 
    MuTree->SetBranchAddress( "NVertices", &fNVertices); 
 

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

    std::cout << "--- Processing: " << MuTree->GetEntries() << " events" << std::endl;
    TStopwatch sw;
    sw.Start();

    int npass   = 0;
    float yield = 0.;

   for (Long64_t ievt=0; ievt<MuTree->GetEntries();ievt++) {

      if (ievt%10000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      MuTree->GetEntry(ievt);

      Double_t rho = 0;
      if (!(TMath::IsNaN(fRho) || isinf(fRho))) rho = fRho;

      //--------------------------------------------------------
      // important: here we associate branches to MVA variables
      //--------------------------------------------------------
      varMuTkNchi2 		  = fMuTkNchi2; 
      varMuGlobalNchi2 	          = fMuGlobalNchi2; 
      varMuNValidHits  	          = fMuNValidHits; 
      varMuNTrackerHits 	  = fMuNTrackerHits; 
      varMuNPixelHits  	          = fMuNPixelHits; 
      varMuNMatches 		  = fMuNMatches; 
      varMuD0  		          = fMuD0; 
      varMuIP3d 		  = fMuIP3d; 
      varMuIP3dSig 		  = fMuIP3dSig; 
      varMuTrkKink 		  = fMuTrkKink; 
      varMuGlobalKink  	          = fMuGlobalKink; 
      varMuSegmentCompatibility   = fMuSegmentCompatibility; 
      varMuCaloCompatibility 	  = fMuCaloCompatibility; 
      varMuHadEnergyOverPt 	  = (fMuHadEnergy - rho*MuonEffectiveArea(kMuHadEnergy,fMuEta))/fMuPt;
      varMuHoEnergyOverPt 	  = (fMuHoEnergy - rho*MuonEffectiveArea(kMuHoEnergy,fMuEta))/fMuPt;
      varMuEmEnergyOverPt 	  = (fMuEmEnergy - rho*MuonEffectiveArea(kMuEmEnergy,fMuEta))/fMuPt;
      varMuHadS9EnergyOverPt 	  = (fMuHadS9Energy - rho*MuonEffectiveArea(kMuHadS9Energy,fMuEta))/fMuPt;
      varMuHoS9EnergyOverPt 	  = (fMuHoS9Energy - rho*MuonEffectiveArea(kMuHoS9Energy,fMuEta))/fMuPt;
      varMuEmS9EnergyOverPt 	  = (fMuEmS9Energy - rho*MuonEffectiveArea(kMuEmS9Energy,fMuEta))/fMuPt;
      varMuChargedIso03OverPt     = (fMuChargedIso03 - rho*MuonEffectiveArea(kMuChargedIso03,fMuEta))/fMuPt;
      varMuNeutralIso03OverPt     = (fMuNeutralIso03_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso03,fMuEta))/fMuPt;
      varMuChargedIso04OverPt     = (fMuChargedIso04 - rho*MuonEffectiveArea(kMuChargedIso04,fMuEta))/fMuPt;
      varMuNeutralIso04OverPt     = (fMuNeutralIso04_05Threshold - rho*MuonEffectiveArea(kMuNeutralIso04,fMuEta))/fMuPt;

      if (versionLabel == "V11") {
        varMuTrkIso03OverPt         = (fMuTrkIso03)/fMuPt;
        varMuEMIso03OverPt          = (fMuEMIso03)/fMuPt;
        varMuHadIso03OverPt         = (fMuHadIso03)/fMuPt;
        varMuTrkIso05OverPt         = (fMuTrkIso05)/fMuPt;
        varMuEMIso05OverPt          = (fMuEMIso05)/fMuPt;
        varMuHadIso05OverPt         = (fMuHadIso05)/fMuPt;
      } else {
        varMuTrkIso03OverPt         = (fMuTrkIso03 - rho*MuonEffectiveArea(kMuTrkIso03,fMuEta))/fMuPt;
        varMuEMIso03OverPt          = (fMuEMIso03 - rho*MuonEffectiveArea(kMuEMIso03,fMuEta))/fMuPt;
        varMuHadIso03OverPt         = (fMuHadIso03 - rho*MuonEffectiveArea(kMuHadIso03,fMuEta))/fMuPt;
        varMuTrkIso05OverPt         = (fMuTrkIso05 - rho*MuonEffectiveArea(kMuTrkIso05,fMuEta))/fMuPt;
        varMuEMIso05OverPt          = (fMuEMIso05 - rho*MuonEffectiveArea(kMuEMIso05,fMuEta))/fMuPt;
        varMuHadIso05OverPt         = (fMuHadIso05 - rho*MuonEffectiveArea(kMuHadIso05,fMuEta))/fMuPt;
      }
      varMuNVtx                   = fNVertices;

    // --- Return the MVA outputs and weights
    //classify by eta and pt bins
    Int_t subdet = 0;
    if (fabs(fMuEta) < 1.479) subdet = 0;
    else subdet = 1;
    Int_t ptBin = 0;
    if (fMuPt > 14.5) ptBin = 1;
    if (fMuPt > 20.0) ptBin = 2;
      
      Int_t MVABin = -1;
      if (subdet == 0 && ptBin == 0) MVABin = 0;
      if (subdet == 1 && ptBin == 0) MVABin = 1;
      if (subdet == 0 && ptBin == 1) MVABin = 2;
      if (subdet == 1 && ptBin == 1) MVABin = 3;
      if (subdet == 0 && ptBin == 2) MVABin = 4;
      if (subdet == 1 && ptBin == 2) MVABin = 5;
      assert(MVABin >= 0 && MVABin <= 5);

      TMVA::Reader  *tmpReader = reader;
      if (evaluateAllBins) tmpReader = fTMVAReader[MVABin];      
      
       if (Use["BDT"]){
        fMuBDT = tmpReader->EvaluateMVA( "BDT method" );
        branchMuBDT->Fill();
      }
      if (Use["BDTG"]){
        fMuBDTG  = tmpReader->EvaluateMVA( "BDTG method" );
        branchMuBDTG->Fill();
      }
      if (Use["MLPBNN"]){
        fMuMLPBNN  = tmpReader->EvaluateMVA( "MLPBNN method" );
        branchMuMLPBNN->Fill();
      }
      if (Use["MLP"]){
        fMuMLP  = tmpReader->EvaluateMVA( "MLP method" );
        branchMuMLP->Fill();
      }
      if (Use["KNN"]){
        fMuKNN  = tmpReader->EvaluateMVA( "KNN method" );
        branchMuKNN->Fill();
      }
       if (Use["Likelihood"]){
        Double_t LH = tmpReader->EvaluateMVA( "Likelihood method" );
        fMuLikelihood = LH;
//         double newLik = 0.0;
//         if     (LH<=0) newLik = -20.0;
//         else if(LH>=1) newLik =  20.0;
//         else                   newLik = log(LH/(1.0-LH));        
//         fMuLikelihood = newLik;
        branchMuLikelihood->Fill();
      }

//        //Combined MVA
//       fMuCombinedMVA  = tmpReader->GetProba( "BDTG method" ) 
//         * tmpReader->GetProba( "BDT method" ) 
//         * tmpReader->GetProba( "MLP method" ) 
//         * tmpReader->GetProba( "MLPBNN method" )  
//         * tmpReader->GetProba( "KNN method" );
//       branchMuCombinedMVA->Fill();
 
    } // End main loop


    std::cout << npass << " events passing selection, yield " << yield << std::endl;
 
    // Get elapsed time
    sw.Stop();
    std::cout << "--- End of event loop: "; sw.Print();

    //Write Output file
    MuOutputFile->Write();
    MuOutputFile->Close();

    delete reader;
    
    std::cout << "==> TMVAClassificationApplication is done with sample " << samples.at(i) << endl << std::endl;
  } 
}
