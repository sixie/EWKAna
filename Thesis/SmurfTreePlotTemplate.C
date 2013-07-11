#include "/home/ceballos/releases/CMSSW_4_2_2/src/Smurf/Core/SmurfTree.h"
#include "/home/ceballos/HiggsMVA/NeuralNetworkMaker-3x/factors.h"
#include "/home/ceballos/HiggsMVA/NeuralNetworkMaker-3x/trilepton.h"
#include "Smurf/HWW/LeptonScaleLookup.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <fstream>
#include "TLegend.h"
#include "TPaveText.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TSystem.h"

const int verboseLevel =   1;

//------------------------------------------------------------------------------
// optimalCuts_42x
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void SmurfTreePlotTemplate
(
  int     mH  ,
  string  Label,
  TString bgdInputFile    = "ntuples_41x/background_mconly.root",
  TString signalInputFile = "ntuples_41x/hww160.root",
  TString dataInputFile   = "ntuples_41x/data2l.root"
  )
{

  string label = Label;
  if (Label != "") label = "_"+Label;

  SmurfTree sigEvent;
  sigEvent.LoadTree(signalInputFile,-1);
  sigEvent.InitTree(0);

  SmurfTree bgdEvent;
  bgdEvent.LoadTree(bgdInputFile,-1);
  bgdEvent.InitTree(0);

  SmurfTree dataEvent;
  dataEvent.LoadTree(dataInputFile,-1);
  dataEvent.InitTree(0);

  int channel = mH;
  int binc = -1;
  if     (mH == 115) binc = 0;
  else if(mH == 120) binc = 1;
  else if(mH == 130) binc = 2;
  else if(mH == 140) binc = 3;
  else if(mH == 150) binc = 4;
  else if(mH == 160) binc = 5;
  else if(mH == 170) binc = 6;
  else if(mH == 180) binc = 7;
  else if(mH == 190) binc = 8;
  else if(mH == 200) binc = 9;
  else if(mH == 210) binc = 10;
  else if(mH == 220) binc = 11;
  else if(mH == 230) binc = 12;
  else if(mH == 250) binc = 13;
  else if(mH == 300) binc = 14;
  else if(mH == 350) binc = 15;
  else if(mH == 400) binc = 16;
  else if(mH == 450) binc = 17;
  else if(mH == 500) binc = 18;
  else if(mH == 550) binc = 19;
  else if(mH == 600) binc = 20;
  else               mH   = 115;

  TFile *fHiggsPtKFactorFile = TFile::Open("/data/smurf/sixie/KFactors/ggHWW_KFactors_PowhegToHQT.root");
  TH1D *HiggsPtKFactor;
  char kfactorHistName[100];
  sprintf(kfactorHistName, "KFactor_PowhegToHQT_mH%d", mH);
  HiggsPtKFactor = (TH1D*)(fHiggsPtKFactorFile->Get(kfactorHistName));
  if (HiggsPtKFactor) {
    HiggsPtKFactor->SetDirectory(0);
  }
  assert(HiggsPtKFactor);
  fHiggsPtKFactorFile->Close();
  delete fHiggsPtKFactorFile;


  double cutMassLow[21]       = { 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12};
  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};
  double cutMetLow[21]        = { 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35};
  double cutMetLowEM[21]      = { 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20};



  //***********************************************************************************************
  //Define Histograms & Yields
  //***********************************************************************************************
  TH1D *Met_HWW160 = new TH1D("Met_HWW160", ";Missing Transverse Energy [GeV/c];Fraction of Events", 60, 0, 120);
  TH1D *ProjectedMet_HWW160 = new TH1D("ProjectedMet_HWW160", ";Projected Missing Transverse Energy [GeV/c];Fraction of Events", 60, 0, 120);
  TH1D *Met_Z = new TH1D("Met_Z", ";Missing Transverse Energy [GeV/c];Fraction of Events", 60, 0, 120);
  TH1D *ProjectedMet_Zdl = new TH1D("ProjectedMet_Z", ";Projected Missing Transverse Energy [GeV/c];Fraction of Events", 60, 0, 120);

  Met->Sumw2();
  ProjectedMet->Sumw2();


  //***********************************************************************************************
  //Signal
  //***********************************************************************************************
  for (int i=0; i<nSig; ++i) {

    if (i%10000 == 0 && verboseLevel > 0)
      printf("--- reading event %5d of %5d\n",i,sigEvent.tree_->GetEntries());
    sigEvent.tree_->GetEntry(i);

    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( Njet3 != nJetsType          					 ) continue; // select n-jet type events
    if( dilep->mass() > dilmass_cut 					 ) continue; // cut on dilepton mass
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if( TMath::Min(pmet,pTrackMet) <= 20                                 ) continue; // cut on pmet for all lepton-pair flavors
    if( TMath::Min(pmet,pTrackMet) <= 40 && 
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on pmet for ee/mm lepton-pair
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair
    if( (cuts & SmurfTree::ExtraLeptonVeto) != SmurfTree::ExtraLeptonVeto) continue; // cut on dileptons
    if( (cuts & patternTopTag) == patternTopTag                          ) continue; // cut on btagging

    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1

    double add = 1.0;
    add = add*nPUScaleFactor(fhDPUS4,npu);

    add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
        					      fabs(lep2->eta()), lep2->pt(), 
        					      TMath::Abs( lid1), TMath::Abs(lid2));

    if (processId == 10010) {
      float theMax = 0.00;
      add = add * HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(TMath::Max(higgsPt, theMax)));
      // add = add * enhancementFactor(mH,false);
    }
    double myWeight = scaleFactorLum * scale1fb * add;

    int nSigBin = -1;
    // GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
    if     (processId==121 ||
    	    processId==122)   nSigBin = 1;
    else if(processId==24)    nSigBin = 2;
    else if(processId==26)    nSigBin = 3;
    else if(processId==10001) nSigBin = 4;
    else if(processId==10010) nSigBin = 5;
    else  {return;}

    bool passAllCuts = dilep->mass()         < cutMassHigh[channel] &&
        	       mt		     > cutMTLow[channel] &&
        	       mt		     < cutMTHigh[channel] &&
        	       lep1->pt()	     > cutPtMaxLow[channel] &&
        	       lep2->pt()	     > cutPtMinLow[channel] &&
        	       dPhi*180.0/TMath::Pi()< cutDeltaphilHigh[channel];
    if(nJetsType == 2){
      int centrality = 0;
      if(((jet1->eta()-lep1->eta() > 0 && jet2->eta()-lep1->eta() < 0) ||
          (jet2->eta()-lep1->eta() > 0 && jet1->eta()-lep1->eta() < 0)) &&
         ((jet1->eta()-lep2->eta() > 0 && jet2->eta()-lep2->eta() < 0) ||
          (jet2->eta()-lep2->eta() > 0 && jet1->eta()-lep2->eta() < 0))) centrality = 1; 
      passAllCuts = (*jet1+*jet2).M() > 450. &&
                    TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 &&
		    (mH > 200 || dilep->mass() < 100.) &&
		    centrality == 1;
    }
    if(passAllCuts == true) {

      
    }



  }







  TFile* file = new TFile("DYBkgPlots.root", "UPDATE");
  file->WriteTObject(hPtMax, hPtMax->GetName(), "WriteDelete");
  file->WriteTObject(hPtMin, hPtMin->GetName(), "WriteDelete");
  file->WriteTObject(hDPhi, hDPhi->GetName(), "WriteDelete");
  file->WriteTObject(hMet, hMet->GetName(), "WriteDelete");
  file->WriteTObject(hMt, hMt->GetName(), "WriteDelete");
  file->WriteTObject(hDileptonMass, hDileptonMass->GetName(), "WriteDelete");
  file->WriteTObject(hDileptonPt, hDileptonPt->GetName(), "WriteDelete");
  file->Close();


  return;

}
