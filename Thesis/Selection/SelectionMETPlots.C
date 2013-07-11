#include "Smurf/Core/SmurfTree.h"
#include "factors.h"
#include "Smurf/Core/LeptonScaleLookup.h"
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
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
#include "TArrow.h"
#include "TLine.h"
#include "TPaveLabel.h"


void NormalizeHist(TH1D *hist) {
  Double_t norm = 0;

  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


//------------------------------------------------------------------------------
// optimalCuts_42x
//------------------------------------------------------------------------------
// GF  == 10010, WBF == 10001, WH == 26, ZH == 24, ttH=121/122
void PlotProjectedMetSelection
(
  int     mH  ,
  string  Label,
  TString signalInputFile = "ntuples_41x/hww160.root",
  TString bgdInputFile    = "ntuples_41x/background_mconly.root"
  )
{

  string label = Label;
  if (Label != "") label = "_"+Label;

  TChain *chsignal = new TChain("tree");
  chsignal->Add("/data/smurf/data/LP2011/mitf_noweights/hww160.root");
  TTree *signal     = (TTree*) chsignal;

  TChain *chbackground = new TChain("tree");
  chbackground->Add("/data/smurf/data/LP2011/mitf_noweights/dymm_nometcut.root");
  chbackground->Add("/data/smurf/data/LP2011/mitf_noweights/dyee_nometcut.root");
  chbackground->Add("/data/smurf/data/LP2011/mitf_noweights/dytt.root");
  TTree *background = (TTree*) chbackground;


  Double_t scaleFactorLum = 1.0;

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

  TFile *fLeptonEffFile = TFile::Open("/data/smurf/sixie/data/Thesis/auxiliar/Winter11_4700ipb/efficiency_results_v7_42x_Full2011_4700ipb.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;

  LeptonScaleLookup trigLookup("/data/smurf/sixie/data/Thesis/auxiliar/Winter11_4700ipb/efficiency_results_v7_42x_Full2011_4700ipb.root");

  TFile *fPUS4File = TFile::Open("/data/smurf/data/LP2011/auxiliar/puWeights_PU4_68mb.root");
  TH1D *fhDPUS4 = (TH1D*)(fPUS4File->Get("puWeights"));
  assert(fhDPUS4);
  fhDPUS4->SetDirectory(0);
  delete fPUS4File;

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
  TH1D *MinMet_HWW160 = new TH1D("MinMet_HWW160", ";Min(E_{T}^{miss}, track-E_{T}^{miss}) [GeV/c];Fraction of Events", 40, 0, 100);
  TH1D *MinProjectedMet_HWW160 = new TH1D("MinProjectedMet_HWW160", ";Min(Projected E_{T}^{miss}, Projected trackMET) [GeV/c];Fraction of Events", 40, 0, 100);
  TH1D *DPhiDileptonJet_HWW160 = new TH1D("DPhiDileptonJet_HWW160", ";#Delta#phi(dilepton system, jet);Fraction of Events", 40, 0, 180);
  TH1D *MinMet_Z = new TH1D("MinMet_Z", ";Min(E_{T}^{miss}, track-E_{T}^{miss}) [GeV/c] ;Fraction of Events", 40, 0, 100);
  TH1D *MinProjectedMet_Z = new TH1D("MinProjectedMet_Z", ";Min(proj. E_{T}^{miss}, proj. track-E_{T}^{miss}) [GeV/c];Fraction of Events", 40, 0, 100);
  TH1D *DPhiDileptonJet_Z = new TH1D("DPhiDileptonJet_Z", ";#Delta#phi(dilepton system, jet);Fraction of Events", 40, 0, 180);
  MinMet_HWW160->Sumw2();
  MinProjectedMet_HWW160->Sumw2();
  DPhiDileptonJet_HWW160->Sumw2();
  MinMet_Z->Sumw2();
  MinProjectedMet_Z->Sumw2();
  DPhiDileptonJet_Z->Sumw2();


  TH1D *NSoftMuons_HWW160 = new TH1D("NSoftMuons_HWW160", ";Number of Soft Muons; Fraction of Events", 10, -0.5, 9.5);
  TH1D *MaxBTag_HWW160 = new TH1D("MaxBTag_HWW160", ";Maximum b-tag discriminator; Fraction of Events", 40, -20, 20);
  TH1D *NJets_HWW160 = new TH1D("NJets_HWW160", ";Number of Jets ( p_{T} > 30 GeV ); Fraction of Events", 10, -0.5, 9.5);
  TH1D *NSoftMuons_Top = new TH1D("NSoftMuons_Top", ";Number of Soft Muons; Fraction of Events", 10, -0.5, 9.5);
  TH1D *MaxBTag_Top = new TH1D("MaxBTag_Top", ";Maximum b-tag discriminator; Fraction of Events", 40, -20, 20);
  TH1D *NJets_Top = new TH1D("NJets_Top", ";Number of Jets ( p_{T} > 30 GeV ); Fraction of Events", 10, -0.5, 9.5);
  NSoftMuons_HWW160->Sumw2();
  MaxBTag_HWW160->Sumw2();
  NJets_HWW160->Sumw2();
  NSoftMuons_Top->Sumw2();
  MaxBTag_Top->Sumw2();
  NJets_Top->Sumw2();


  //----------------------------------------------------------------------------
  UInt_t          cuts;
  UInt_t          dstype;
  UInt_t          nvtx;
  UInt_t          npu;
  UInt_t          njets;
  UInt_t          run;
  UInt_t          event;
  Float_t         scale1fb;
  LorentzVector*  lep1  = 0;
  LorentzVector*  lep2  = 0;
  LorentzVector*  jet1  = 0;
  LorentzVector*  jet2  = 0;
  LorentzVector*  jet3  = 0;
  Float_t         dPhi;
  Float_t         dR;
  LorentzVector*  dilep = 0;
  UInt_t          type;
  Float_t         met;
  Float_t         trackMet;
  Float_t         pmet;
  Float_t         pTrackMet;
  Float_t         mt;
  Float_t         mt1;
  Float_t         mt2;
  Float_t         dPhiLep1MET;
  Float_t         dPhiLep2MET;
  Float_t         dPhiDiLepMET;
  Float_t         dPhiDiLepJet1;
  Int_t           lq1;
  Int_t           lq2;
  Int_t           lid1;
  Int_t           lid2;
  Int_t           lid3;
  Int_t           processId;
  Float_t         jetLowBtag;
  UInt_t          nSoftMuons;
  Float_t         jet1Btag;
  Float_t         jet2Btag;
  Int_t 	  lep1McId;
  Int_t 	  lep2McId;
  Int_t 	  lep1MotherMcId;
  Int_t 	  lep2MotherMcId;
  Float_t         bdt = 0.0;
  Float_t         bdtd = 0.0;
  Float_t         nn = 0.0;
  Float_t         knn = 0.0;
  Float_t         bdtg = 0.0;
  Float_t         higgsPt = -999;

  //***********************************************************************************************
  //Signal
  //***********************************************************************************************

  signal->SetBranchAddress( "cuts"         , &cuts         );
  signal->SetBranchAddress( "dstype"       , &dstype       );
  signal->SetBranchAddress( "nvtx"         , &nvtx         );
  signal->SetBranchAddress( "npu"          , &npu          );
  signal->SetBranchAddress( "njets"        , &njets        );
  signal->SetBranchAddress( "run"          , &run          );
  signal->SetBranchAddress( "event"        , &event        );
  signal->SetBranchAddress( "scale1fb"     , &scale1fb     );
  signal->SetBranchAddress( "lep1"         , &lep1         );
  signal->SetBranchAddress( "lep2"         , &lep2         );
  signal->SetBranchAddress( "jet1"         , &jet1         );
  signal->SetBranchAddress( "jet2"         , &jet2         );
  signal->SetBranchAddress( "jet3"         , &jet3         );
  signal->SetBranchAddress( "dPhi"         , &dPhi         );
  signal->SetBranchAddress( "dR"           , &dR           );
  signal->SetBranchAddress( "dilep"        , &dilep        );
  signal->SetBranchAddress( "type"         , &type         );
  signal->SetBranchAddress( "met"          , &met         );
  signal->SetBranchAddress( "trackMet"     , &trackMet    );
  signal->SetBranchAddress( "pmet"         , &pmet         );
  signal->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  signal->SetBranchAddress( "met"          , &met          );
  signal->SetBranchAddress( "mt"           , &mt           );
  signal->SetBranchAddress( "mt1"          , &mt1          );
  signal->SetBranchAddress( "mt2"          , &mt2          );
  signal->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  signal->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  signal->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  signal->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  signal->SetBranchAddress( "lq1"          , &lq1          );
  signal->SetBranchAddress( "lq2"          , &lq2          );
  signal->SetBranchAddress( "lid1"         , &lid1         );
  signal->SetBranchAddress( "lid2"         , &lid2         );
  signal->SetBranchAddress( "lid3"         , &lid3         );
  signal->SetBranchAddress( "processId"    , &processId    );
  signal->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  signal->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  signal->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  signal->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  signal->SetBranchAddress( "higgsPt"      , &higgsPt      );

  for (UInt_t i=0; i<signal->GetEntries(); i++) {

    if (i%100000 == 0 )
      printf("--- reading event %5d of %5d\n",i,signal->GetEntries());
    signal->GetEntry(i);

    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair

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

    if (processId == 10010) {
      MinMet_HWW160->Fill(min(met,trackMet),myWeight);
      MinProjectedMet_HWW160->Fill(min(pmet,pTrackMet), myWeight);

      if ((type == SmurfTree::ee || type == SmurfTree::mm) && jet1->pt() > 15
          && TMath::Min(pmet,pTrackMet) > 40 ) {
        DPhiDileptonJet_HWW160->Fill(dPhiDiLepJet1*180.0/TMath::Pi(),myWeight);
      }
    }
    NJets_HWW160->Fill(njets, myWeight);
    
  }



  //***********************************************************************************************
  //Background
  //***********************************************************************************************

  background->SetBranchAddress( "cuts"         , &cuts         );
  background->SetBranchAddress( "dstype"       , &dstype       );
  background->SetBranchAddress( "nvtx"         , &nvtx         );
  background->SetBranchAddress( "npu"          , &npu          );
  background->SetBranchAddress( "njets"        , &njets        );
  background->SetBranchAddress( "run"          , &run          );
  background->SetBranchAddress( "event"        , &event        );
  background->SetBranchAddress( "scale1fb"     , &scale1fb     );
  background->SetBranchAddress( "lep1"         , &lep1         );
  background->SetBranchAddress( "lep2"         , &lep2         );
  background->SetBranchAddress( "jet1"         , &jet1         );
  background->SetBranchAddress( "jet2"         , &jet2         );
  background->SetBranchAddress( "jet3"         , &jet3         );
  background->SetBranchAddress( "dPhi"         , &dPhi         );
  background->SetBranchAddress( "dR"           , &dR           );
  background->SetBranchAddress( "dilep"        , &dilep        );
  background->SetBranchAddress( "type"         , &type         );
  background->SetBranchAddress( "met"          , &met         );
  background->SetBranchAddress( "trackMet"     , &trackMet    );
  background->SetBranchAddress( "pmet"         , &pmet         );
  background->SetBranchAddress( "pTrackMet"    , &pTrackMet    );
  background->SetBranchAddress( "met"          , &met          );
  background->SetBranchAddress( "mt"           , &mt           );
  background->SetBranchAddress( "mt1"          , &mt1          );
  background->SetBranchAddress( "mt2"          , &mt2          );
  background->SetBranchAddress( "dPhiLep1MET"  , &dPhiLep1MET  );
  background->SetBranchAddress( "dPhiLep2MET"  , &dPhiLep2MET  );
  background->SetBranchAddress( "dPhiDiLepMET" , &dPhiDiLepMET );
  background->SetBranchAddress( "dPhiDiLepJet1", &dPhiDiLepJet1);
  background->SetBranchAddress( "lq1"          , &lq1          );
  background->SetBranchAddress( "lq2"          , &lq2          );
  background->SetBranchAddress( "lid1"         , &lid1         );
  background->SetBranchAddress( "lid2"         , &lid2         );
  background->SetBranchAddress( "lid3"         , &lid3         );
  background->SetBranchAddress( "processId"    , &processId    );
  background->SetBranchAddress( "jetLowBtag"   , &jetLowBtag   );
  background->SetBranchAddress( "nSoftMuons"   , &nSoftMuons   );
  background->SetBranchAddress( "jet1Btag"     , &jet1Btag     );
  background->SetBranchAddress( "jet2Btag"     , &jet2Btag     );
  background->SetBranchAddress( "higgsPt"      , &higgsPt      );

  for (UInt_t i=0; i<background->GetEntries(); i++) {

    if (i%100000 == 0 )
      printf("--- reading event %5d of %5d\n",i,background->GetEntries());
    background->GetEntry(i);

    if (!(((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection) 
          && ((cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection))) continue;
    if( lq1*lq2 > 0                 					 ) continue; // cut on opposite-sign leptons
    if( dilep->mass() <= 12.0       					 ) continue; // cut on low dilepton mass
    if( lep1->pt() <= 20	    					 ) continue; // cut on leading lepton pt
    if( lep2->pt() <= 10	    					 ) continue; // cut on trailing lepton pt
    if(fabs(dilep->mass()-91.1876) <= 15 &&
      (type == SmurfTree::mm || 
       type == SmurfTree::ee)                                            ) continue; // cut on Z veto for ee/mm lepton-pair

    double add = 1.0;
    add = add*nPUScaleFactor(fhDPUS4,npu);

    add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
        					      fabs(lep2->eta()), lep2->pt(), 
        					      TMath::Abs( lid1), TMath::Abs(lid2));

    double myWeight = scaleFactorLum * scale1fb * add;

    MinMet_Z->Fill(min(met,trackMet),myWeight);
    MinProjectedMet_Z->Fill(min(pmet,pTrackMet), myWeight);
    if ((type == SmurfTree::ee || type == SmurfTree::mm) && jet1->pt() > 15
      && TMath::Min(pmet,pTrackMet) > 40 ) {
      DPhiDileptonJet_Z->Fill(dPhiDiLepJet1*180.0/TMath::Pi(),myWeight);
    }      

  }


  //********************************************************
  // Normalize Hists
  //********************************************************
  NormalizeHist(MinMet_HWW160);
  NormalizeHist(MinProjectedMet_HWW160);
  NormalizeHist(DPhiDileptonJet_HWW160);
  NormalizeHist(MinMet_Z);
  NormalizeHist(MinProjectedMet_Z);
  NormalizeHist(DPhiDileptonJet_Z);


  //********************************************************
  // Save Plots
  //********************************************************
  TFile* file = new TFile("ThesisPlots.SignalCharacteristics.root", "UPDATE");
  file->WriteTObject(MinMet_HWW160, MinMet_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(MinProjectedMet_HWW160, MinProjectedMet_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(DPhiDileptonJet_HWW160, DPhiDileptonJet_HWW160->GetName(), "WriteDelete");
  file->WriteTObject(MinMet_Z, MinMet_Z->GetName(), "WriteDelete");
  file->WriteTObject(MinProjectedMet_Z, MinProjectedMet_Z->GetName(), "WriteDelete");
  file->WriteTObject(DPhiDileptonJet_Z, DPhiDileptonJet_Z->GetName(), "WriteDelete");
  file->Close();


  return;

}





void PlotMetPlots() {

  TFile *file = new TFile("ThesisPlots.SignalCharacteristics.root", "READ");

  TH1D *MinMet_HWW160 = (TH1D*)file->Get("MinMet_HWW160");
  TH1D *MinProjectedMet_HWW160 = (TH1D*)file->Get("MinProjectedMet_HWW160");
  TH1D *DPhiDileptonJet_HWW160 = (TH1D*)file->Get("DPhiDileptonJet_HWW160");
  TH1D *MinMet_Z = (TH1D*)file->Get("MinMet_Z");
  TH1D *MinProjectedMet_Z = (TH1D*)file->Get("MinProjectedMet_Z");
  TH1D *DPhiDileptonJet_Z = (TH1D*)file->Get("DPhiDileptonJet_Z");



//   //********************************************************
//   // Draw
//   //********************************************************
  TCanvas *cv = 0;
  TLegend *legend = 0;


  cv = new TCanvas("cv", "cv", 800, 600);
  cv->SetLogy();
  legend = new TLegend(0.5,0.75,0.95,0.90);   
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(MinMet_HWW160, "H->WW (m_{H} = 160 GeV)" , "L");
  legend->AddEntry(MinMet_Z, "Z/#gamma^{*} #rightarrow #mu#mu" , "L");

  MinMet_HWW160->SetMaximum(10*max(MinMet_HWW160->GetMaximum(),MinMet_Z->GetMaximum()));
  MinMet_HWW160->SetMinimum(1e-5);
  MinMet_HWW160->GetYaxis()->SetTitleOffset(1.4);
  MinMet_HWW160->GetXaxis()->SetTitleOffset(1.05);
  MinMet_HWW160->SetLineColor(kBlack);
  MinMet_HWW160->Draw("hist");
  MinMet_Z->SetLineColor(kRed);
  MinMet_Z->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Selection_HWW160VsZ_MinMet.gif");
  cv->SaveAs("Selection_HWW160VsZ_MinMet.eps");

  cv = new TCanvas("cv", "cv", 800, 600);
  cv->SetLogy();
  legend = new TLegend(0.5,0.75,0.95,0.90);   
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(MinProjectedMet_HWW160, "H->WW (m_{H} = 160 GeV)" , "L");
  legend->AddEntry(MinProjectedMet_Z, "Z/#gamma^{*} #rightarrow #mu#mu" , "L");

  MinProjectedMet_HWW160->SetMaximum(10*max(MinProjectedMet_HWW160->GetMaximum(),MinProjectedMet_Z->GetMaximum()));
  MinProjectedMet_HWW160->SetMinimum(1e-5);
  MinProjectedMet_HWW160->GetYaxis()->SetTitleOffset(1.4);
  MinProjectedMet_HWW160->GetXaxis()->SetTitleOffset(1.05);
  MinProjectedMet_HWW160->SetLineColor(kBlack);
  MinProjectedMet_HWW160->Draw("hist");
  MinProjectedMet_Z->SetLineColor(kRed);
  MinProjectedMet_Z->Draw("hist,same");
  legend->Draw();
  cv->SaveAs("Selection_HWW160VsZ_MinProjectedMet.gif");
  cv->SaveAs("Selection_HWW160VsZ_MinProjectedMet.eps");

  cv = new TCanvas("cv", "cv", 800, 600);
  legend = new TLegend(0.5,0.75,0.95,0.90);   
  legend->SetTextSize(0.04);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(DPhiDileptonJet_HWW160, "H->WW (m_{H} = 160 GeV)" , "L");
  legend->AddEntry(DPhiDileptonJet_Z, "Z/#gamma^{*} #rightarrow #mu#mu" , "L");

  DPhiDileptonJet_HWW160->SetMaximum(1.2*max(DPhiDileptonJet_HWW160->GetMaximum(),DPhiDileptonJet_Z->GetMaximum()));
  DPhiDileptonJet_HWW160->SetMinimum(1e-5);
  DPhiDileptonJet_HWW160->GetYaxis()->SetTitleOffset(1.4);
  DPhiDileptonJet_HWW160->GetXaxis()->SetTitleOffset(1.05);
  DPhiDileptonJet_HWW160->SetLineColor(kBlack);
  DPhiDileptonJet_HWW160->Draw("hist");
  DPhiDileptonJet_Z->SetLineColor(kRed);
  DPhiDileptonJet_Z->Draw("hist,same");
  legend->Draw();

  TLine *cutline = new TLine(165,0.2,165,0.0);
  cutline->SetLineColor(kBlue);
  cutline->SetLineWidth(4);
  cutline->Draw();

  TArrow *cutarrow = new TArrow(165,0.15,150,0.15,0.03,"|>");
  cutarrow->SetLineColor(kBlue);
  cutarrow->SetLineWidth(4);
  cutarrow->SetArrowSize(0.03);
  cutarrow->SetAngle(30);
  cutarrow->Draw();


  TPaveLabel *label = new TPaveLabel(120,0.16,160,0.195,"Cut : #Delta#phi < 165");
  label->SetBorderSize(0);
  label->SetTextSize(0.6);
  label->SetFillStyle(0);
  label->Draw();

  cv->SaveAs("Selection_HWW160VsZ_DPhiDileptonJet.gif");
  cv->SaveAs("Selection_HWW160VsZ_DPhiDileptonJet.eps");


}


void SelectionMETPlots(Int_t Step = 0) {

  if (Step == 0) {
    PlotProjectedMetSelection( 160 , "",
                               "/data/smurf/data/LP2011/mitf_noweights/hww160.root",
                               "/data/smurf/data/LP2011/mitf_noweights/dymm_nometcut.root"
    );
  } 

  if (Step == 1) {
    PlotMetPlots();
  }

}
