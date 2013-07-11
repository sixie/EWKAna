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
#include "EWKAna/Thesis/Selection/StandardPlot.C"


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
void PlotHiggsSelection
(
  int     mH  ,
  string  Label
  )
{

  string label = Label;
  if (Label != "") label = "_"+Label;

  TChain *chsignal = new TChain("tree");
  chsignal->Add("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/ntuples_160train_0jets_hww160.root");
//   chsignal->Add("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/ntuples_120train_0jets_hww120.root");
  TTree *signal     = (TTree*) chsignal;

  TChain *chbackground = new TChain("tree");
  chbackground->Add("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/mitf-alljets/ntuples_160train_0jets_backgroundA_skim2.root");
  TTree *background = (TTree*) chbackground;


  double scaleFactorLum = 1.545;

  int channel = -1;
  if     (mH == 115) channel = 0;
  else if(mH == 120) channel = 1;
  else if(mH == 130) channel = 2;
  else if(mH == 140) channel = 3;
  else if(mH == 150) channel = 4;
  else if(mH == 160) channel = 5;
  else if(mH == 170) channel = 6;
  else if(mH == 180) channel = 7;
  else if(mH == 190) channel = 8;
  else if(mH == 200) channel = 9;
  else if(mH == 210) channel = 10;
  else if(mH == 220) channel = 11;
  else if(mH == 230) channel = 12;
  else if(mH == 250) channel = 13;
  else if(mH == 300) channel = 14;
  else if(mH == 350) channel = 15;
  else if(mH == 400) channel = 16;
  else if(mH == 450) channel = 17;
  else if(mH == 500) channel = 18;
  else if(mH == 550) channel = 19;
  else if(mH == 600) channel = 20;

  if(channel == -1) return;

  float dilmass_cut = 10000;
   
  if     ( mH == 115 ) dilmass_cut =  70.0;
  else if( mH == 120 ) dilmass_cut =  70.0;
  else if( mH == 130 ) dilmass_cut =  80.0;
  else if( mH == 140 ) dilmass_cut =  90.0;
  else if( mH == 150 ) dilmass_cut = 100.0;
  else if( mH == 160 ) dilmass_cut = 100.0;
  else if( mH == 165 ) dilmass_cut = 100.0;
  else if( mH == 170 ) dilmass_cut = 100.0;
  else if( mH == 180 ) dilmass_cut = 110.0;
  else if( mH == 190 ) dilmass_cut = 120.0;
  else if( mH == 200 ) dilmass_cut = 130.0;
  else if( mH == 210 ) dilmass_cut = 140.0;
  else if( mH == 220 ) dilmass_cut = 150.0;
  else                 dilmass_cut = mH;

  double wwScaleFactor0jCut[21] = {1.07,1.07,1.07,1.05,0.98,0.99,0.98,0.98,0.99,0.98,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor0jMVA[21] = {1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,1.07,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor1jCut[21] = {1.11,1.11,1.12,1.11,1.09,1.09,1.10,1.12,1.12,1.14,
                                   1,1,1,1,1,1,1,1,1,1,1};
  double wwScaleFactor1jMVA[21] = {1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,1.11,
                                   1,1,1,1,1,1,1,1,1,1,1};

  double zjScaleFactor[3][21] = {
  {1.675, 1.536, 1.013, 1.126, 0.520, 0.755, 1.265, 2.132, 0.795, 0.460, 1.000, 1.000, 1.000, 0.159, 0.788, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {0.381, 0.658, 0.438, 0.374, 0.847, 0.811, 1.498, 0.698, 1.416, 1.730, 1.000, 1.000, 1.000, 1.295, 1.278, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000},
  {1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000, 1.000}
  };
  double zjScaleFactorE[3][21] = {
  {2.893, 2.909, 2.307, 1.932, 5.153, 1.952, 2.691, 2.102, 3.029, 6.022, 1.000, 1.000, 7.083, 5.063, 0.609, 0.609, 0.609, 0.609, 0.609, 0.609, 0.609},
  {2.262, 2.243, 3.837, 2.743, 1.861, 1.259, 1.033, 1.388, 1.623, 1.451, 1.000, 1.000, 1.000, 0.319, 1.270, 0.498, 0.498, 0.498, 0.498, 0.498, 0.498},
  {0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378, 0.378}
  };
  double zjScaleFactorWWE[3] = {
   0.609,
   0.498,
   0.378
  };
  double ZXS_E[3] = {0.0, 0.0, 0.0};
  double DYXS[3]  = {0.0, 0.0, 0.0};
  double VVXS[3]  = {0.0, 0.0, 0.0};


  TFile *fLeptonEffFile = TFile::Open("/data/smurf/data/Run2011_Spring11_SmurfV6_42X/lepton_eff/efficiency_results_v6.root");
  TH2D *fhDEffMu = (TH2D*)(fLeptonEffFile->Get("h2_results_muon_selection"));
  TH2D *fhDEffEl = (TH2D*)(fLeptonEffFile->Get("h2_results_electron_selection"));
  fhDEffMu->SetDirectory(0);
  fhDEffEl->SetDirectory(0);
  fLeptonEffFile->Close();
  delete fLeptonEffFile;
  LeptonScaleLookup trigLookup("/data/smurf/data/EPS/auxiliar/efficiency_results_v6.root");

  TFile *fLeptonFRFileM = TFile::Open("/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root");
  TH2D *fhDFRMu = (TH2D*)(fLeptonFRFileM->Get("MuonFakeRate_M2_ptThreshold15_PtEta"));
  assert(fhDFRMu);
  fhDFRMu->SetDirectory(0);
  fLeptonFRFileM->Close();
  delete fLeptonFRFileM;

  TFile *fLeptonFRFileE = TFile::Open("/data/smurf/data/LP2011/auxiliar/FakeRates_SmurfV6.V4HasNod0Cut.root");
  TH2D *fhDFREl = (TH2D*)(fLeptonFRFileE->Get("ElectronFakeRate_V4_ptThreshold35_PtEta"));
  assert(fhDFREl);
  fhDFREl->SetDirectory(0);
  fLeptonFRFileE->Close();
  delete fLeptonFRFileE;


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


  //----------------------------------------------------------------------------

  double cutMassHigh[21]      = { 40, 40, 45, 45, 50, 50, 50, 60, 80, 90,110,120,130,150,200,250,300,350,400,450,500};
  double cutPtMaxLow[21]      = { 20, 20, 25, 25, 27, 30, 34, 36, 38, 40, 44, 48, 52, 55, 70, 80, 90,110,120,130,140};
  double cutPtMinLow[21]      = { 10, 10, 10, 15, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25, 25};
  double cutDeltaphilHigh[21] = {115,115, 90, 90, 90, 60, 60, 70, 90,100,110,120,130,140,175,175,175,175,175,175,175};
  double cutMTLow[21]         = { 70, 70, 75, 80, 80, 90,110,120,120,120,120,120,120,120,120,120,120,120,120,120,120};
  double cutMTHigh[21]        = {110,120,125,130,150,160,170,180,190,200,210,220,230,250,300,350,400,450,500,550,600};



  //***********************************************************************************************
  //Define Histograms & Yields
  //***********************************************************************************************

  //----------------------------------------------------------------------------
  const int nHist = 8;
  int    nBinHis[nHist]  = {    200,     200,            200,        180,  200,          200, 200, 200};
  double minHis[nHist]   = {      0,       0,              0,          0,    0,            0,   0,   0};
  double maxHis[nHist]   = {    100,     100,            200,        180,  200,          100, 7.0, 900};
  string HistName[nHist] = {"PtMax", "PtMin", "DileptonMass", "DeltaPhi", "Mt", "DileptonPt", 
                            "DeltaEtaVBFJets", "DijetMassVBFJets"};
  string HistAxisLabel[nHist] = {"Leading Lepton p_{T} [GeV/c]", "Trailing Lepton p_{T} [GeV/c]",
                                 "M_{ll} [GeV/c^{2}]", "#Delta#phi(l,l)",
                                 "m_{T}^{ll E_{T}^{miss}}",
                                 "p_{T}^{ll} [GeV/c]",
                                 "#Delta#eta(j,j)",
                                 "M_{jj} [GeV/c^{2}]" };

  const int nBkgProcesses = 5;
  string bkgProcessLabel[nBkgProcesses] = {"WW","WZZZ","Top","DYeemm","WJets"};
  //----------------------------------------------------------------------------
  TH1D* signalDistributions[nHist];
  TH1D* bkgDistributions[nHist][nBkgProcesses];
  for(int i=0; i<nHist; i++) {
    signalDistributions[i] = new TH1D( ("signalDistributions_" + HistName[i]).c_str(), 
                                       (";"+HistAxisLabel[i]+";Number of Events").c_str(), 
                                       nBinHis[i], minHis[i], maxHis[i]);
    signalDistributions[i]->Sumw2();

    for(int j=0; j<nBkgProcesses; j++){
      bkgDistributions[i][j] = new TH1D( ("bkgDistributions_" + HistName[i]+"_"+bkgProcessLabel[j]).c_str(), 
                                      (";"+HistAxisLabel[i]+";Number of Events").c_str(), 
                                      nBinHis[i], minHis[i], maxHis[i]);
      bkgDistributions[i][j]->Sumw2();
    }
  }



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
    if( (cuts & SmurfTree::TopTag) == SmurfTree::TopTag                  ) continue; // cut on btagging
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

    if (njets == 0) {
      signalDistributions[0]->Fill(lep1->pt(), myWeight);
      signalDistributions[1]->Fill(lep2->pt(), myWeight);
      signalDistributions[2]->Fill(dilep->mass(), myWeight);
      signalDistributions[3]->Fill(dPhi*180.0/TMath::Pi(), myWeight);
      signalDistributions[4]->Fill(mt, myWeight);
      signalDistributions[5]->Fill(dilep->pt(), myWeight);
    }

    //VBF
    if (processId == 10001) {
      if (jet1->pt() > 30 && jet2->pt() > 30) {
        signalDistributions[6]->Fill(TMath::Abs(jet1->eta()-jet2->eta()), myWeight);
        signalDistributions[7]->Fill((*jet1+*jet2).M(), myWeight);
      }
    }


  }



  //***********************************************************************************************
  //Background
  //***********************************************************************************************
  Double_t all = 0;
  Double_t pass = 0;

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

    bool MinPreselCut = ((cuts & SmurfTree::Lep1FullSelection) == SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) ||
      ((cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection && (cuts & SmurfTree::Lep2FullSelection) == SmurfTree::Lep2FullSelection) ||
      dstype != SmurfTree::data;
    if( MinPreselCut == false                                            ) continue; // cut on MinPreselCut
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
    if( (cuts & SmurfTree::TopTag) == SmurfTree::TopTag                  ) continue; // cut on btagging
    bool dPhiDiLepJet1Cut = jet1->pt() <= 15. ||
                           (dPhiDiLepJet1*180.0/TMath::Pi() < 165. || type == SmurfTree::em || type == SmurfTree::me);
    if( dPhiDiLepJet1Cut == false                                        ) continue; // cut on dPhiDiLepJet1


    int BkgType = 0;
    if(dstype == SmurfTree::qqww  ) BkgType = 0;
    else if(dstype == SmurfTree::ggww  ) BkgType = 0;
    else if(dstype == SmurfTree::wz    ) BkgType = 1;
    else if(dstype == SmurfTree::zz    ) BkgType = 1;
    else if(dstype == SmurfTree::tw    ) BkgType = 2;
    else if(dstype == SmurfTree::ttbar ) BkgType = 2;
    else if(dstype == SmurfTree::dyee  ) BkgType = 3;
    else if(dstype == SmurfTree::dymm  ) BkgType = 3;
    else if(dstype == SmurfTree::dytt  ) BkgType = 3;
    else if(dstype == SmurfTree::wjets ) BkgType = 4;
    else if(dstype == SmurfTree::data  ) BkgType = 4;
    else if(dstype == SmurfTree::wgamma) BkgType = 4;



    double myWeight = 1.0;
    double add      = 1.0;
    int nFake = 0;
    if(dstype == SmurfTree::data ){
      if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    } else {
      if(((cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2)  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2)  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4) && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection) nFake++;
      if(((cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4) && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection) nFake++;
    }
    if(nFake > 1){
      myWeight = 0.0;
    }
    else if(nFake == 1){
      if(dstype == SmurfTree::data){
        add = add*fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
        							      (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
        							      (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);
        BkgType = 4;
        myWeight = add;
       }
      else if(TMath::Abs(lep1McId)*TMath::Abs(lep2McId) > 0 || dstype == SmurfTree::wgamma){
    	add = 1.0;
    	add = add*fakeRate(lep1->pt(), lep1->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep1LooseMuV2)  == SmurfTree::Lep1LooseMuV2  && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection, 
    										      (cuts & SmurfTree::Lep1LooseEleV4) == SmurfTree::Lep1LooseEleV4 && (cuts & SmurfTree::Lep1FullSelection) != SmurfTree::Lep1FullSelection);
    	add = add*fakeRate(lep2->pt(), lep2->eta(), fhDFRMu, fhDFREl, (cuts & SmurfTree::Lep2LooseMuV2)  == SmurfTree::Lep2LooseMuV2  && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection,
    										      (cuts & SmurfTree::Lep2LooseEleV4) == SmurfTree::Lep2LooseEleV4 && (cuts & SmurfTree::Lep2FullSelection) != SmurfTree::Lep2FullSelection);

        add = add*nPUScaleFactor(fhDPUS4,npu);

    	add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
    	add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);

    	add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							  fabs(lep2->eta()), lep2->pt(), 
    							  TMath::Abs( lid1), TMath::Abs(lid2));
        BkgType = 4;
        myWeight	       = -1.0 * scale1fb*scaleFactorLum*add;
      }
      else {
        myWeight = 0.0;
      }
    }
    else if(dstype == SmurfTree::data) myWeight = 0.0;
    else if(dstype != SmurfTree::data){
      add = 1.0;
      add = add*nPUScaleFactor(fhDPUS4,npu);

      add = add*leptonEfficiency(lep1->pt(), lep1->eta(), fhDEffMu, fhDEffEl, lid1);
      add = add*leptonEfficiency(lep2->pt(), lep2->eta(), fhDEffMu, fhDEffEl, lid2);
      add = add*trigLookup.GetExpectedTriggerEfficiency(fabs(lep1->eta()), lep1->pt() , 
    							fabs(lep2->eta()), lep2->pt(), 
    							TMath::Abs( lid1), TMath::Abs(lid2));
      if(BkgType == 3  && (type   == SmurfTree::mm   || type   == SmurfTree::ee)
                      && (dstype == SmurfTree::dyee || dstype == SmurfTree::dymm)) {
    	if(njets == 0) add=add*3.02; 
    	if(njets == 1) add=add*2.81; 
    	if(njets >= 2) add=add*4.84; 
      }
      if(BkgType == 2) {
    	if(njets == 0) add=add*1.51;
    	if(njets == 1) add=add*1.16; 
    	if(njets >= 2) add=add*1.00; 
      }
      if(BkgType == 0){     
        if(njets == 0) add=add*wwScaleFactor0jMVA[channel]; 
        else	       add=add*wwScaleFactor1jMVA[channel]; 
      }
      // CAREFUL HERE, no data-driven corrections, just Higgs k-factors
      // add = 1.0;
      myWeight = scale1fb*scaleFactorLum*add;
    }

    if(myWeight == 0) continue;

    if (njets == 0) {
      bkgDistributions[0][BkgType]->Fill(lep1->pt(), myWeight);
      bkgDistributions[1][BkgType]->Fill(lep2->pt(), myWeight);
      bkgDistributions[2][BkgType]->Fill(dilep->mass(), myWeight);
      bkgDistributions[3][BkgType]->Fill(dPhi*180.0/TMath::Pi(), myWeight);
      bkgDistributions[4][BkgType]->Fill(mt, myWeight);
      bkgDistributions[5][BkgType]->Fill(dilep->pt(), myWeight);
    }

    //VBF
    if (jet1->pt() > 30 && jet2->pt() > 30) {
      all += myWeight;
      if (TMath::Abs(jet1->eta()-jet2->eta()) > 3.5 && (*jet1+*jet2).M() > 450) pass += myWeight;
      bkgDistributions[6][BkgType]->Fill(TMath::Abs(jet1->eta()-jet2->eta()), myWeight);
      bkgDistributions[7][BkgType]->Fill((*jet1+*jet2).M(), myWeight);
    }

  }

  cout << pass << " / " << all << " = " << pass/all << endl;

  //********************************************************
  // Save Plots
  //********************************************************
  TFile* file = new TFile("ThesisPlots.SignalCharacteristics.root", "UPDATE");
  for(int i=0; i<nHist; i++) {
    file->WriteTObject(signalDistributions[i], signalDistributions[i]->GetName(), "WriteDelete");
    for(int j=0; j<nBkgProcesses; j++){
      file->WriteTObject(bkgDistributions[i][j], bkgDistributions[i][j]->GetName(), "WriteDelete");
    }
  }
  file->Close();

  return;

}


void MakePlot(int nsel = 0, int ReBin = 1, char XTitle[300] = "N_{jets}", char units[300] = "", string histname = "", char outputName[300] = "njets",
              bool isLogY = false, int MassH = 160, double lumi = 1.55, Double_t lowCutValue = -999, Double_t highCutValue = -999) {

// void MakePlot(string histname, int ReBin = 1, ) {

//     gInterpreter->ExecuteMacro("EWKAna/Thesis/Selection/GoodStyle.C");
    gROOT->LoadMacro("EWKAna/Thesis/Selection/StandardPlot.C");

    TFile* file = new TFile("ThesisPlots.SignalCharacteristics.root", "read");

    StandardPlot myPlot;
    myPlot.setLumi(lumi);
    myPlot.setLabel(XTitle);
    myPlot.addLabel("");
    myPlot.setUnits(units);

    TH1F* hWW    = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WW").c_str());
    TH1F* hVV    = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WZZZ").c_str());
    TH1F* hTop   = (TH1F*)file->Get(("bkgDistributions_" + histname + "_Top").c_str());
    TH1F* hZJets = (TH1F*)file->Get(("bkgDistributions_" + histname + "_DYeemm").c_str());
    TH1F* hWJets = (TH1F*)file->Get(("bkgDistributions_" + histname + "_WJets").c_str());
    
    assert(hWW);
    assert(hVV);
    assert(hTop);
    assert(hZJets);
    assert(hWJets);

    if(nsel == 0 || nsel == 1){
      if(hWW->GetSumOfWeights(   ) > 0) myPlot.setMCHist(iWW,    (TH1F*)hWW   ->Clone("hWW"));
      if(hVV->GetSumOfWeights()    > 0) myPlot.setMCHist(iVV,    (TH1F*)hVV   ->Clone("hVV")); 
      if(hTop->GetSumOfWeights()   > 0) myPlot.setMCHist(iTop,   (TH1F*)hTop  ->Clone("hTop"));
      if(hZJets->GetSumOfWeights() > 0) myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
      if(hWJets->GetSumOfWeights() > 0) myPlot.setMCHist(iWJets, (TH1F*)hWJets->Clone("hWJets"));
    }
    else if(nsel == 2) {
      //myPlot.setMCHist(iWW,    (TH1F*)hWW   ->Clone("hWW"));
      //myPlot.setMCHist(iZJets, (TH1F*)hZJets->Clone("hZJets"));
      myPlot.setMCHist(iZZ,    (TH1F*)hTop  ->Clone("hTop"));
      myPlot.setMCHist(iWZ,    (TH1F*)hVV   ->Clone("hVV")); 
      myPlot.setMCHist(iFakes, (TH1F*)hWJets->Clone("hWJets"));
    }

    TH1F* hHWW   = (TH1F*)file->Get(("signalDistributions_" + histname).c_str());
    assert(hHWW);



    if(nsel != 1){
      myPlot.setMCHist(iHWW, (TH1F*)hHWW->Clone("hHWW"));
      myPlot._mass = MassH;
    }

    if (lowCutValue != -999) myPlot.setLowCutValue(lowCutValue);
    if (highCutValue != -999) myPlot.setHighCutValue(highCutValue);
    myPlot.Draw(outputName, isLogY, ReBin );  // Can pass a rebin 

    return;

}






void SelectionHiggsCutsPlots(Int_t Step = 1) {

  if (Step == 0) {
    PlotHiggsSelection( 160 , "");
  } 

  if (Step == 1) {
    MakePlot(0, 4, "p_{T}^{Leading Lepton}","GeV/c","PtMax", "Selection_HWW160Selection_PtMax",0,160,1.55, 30, -999 );
    MakePlot(0, 4, "p_{T}^{Trailing Lepton}","GeV/c","PtMin", "Selection_HWW160Selection_PtMin",0,160,1.55, 25, -999 );
    MakePlot(0, 4, "M_{ll}","GeV/c^{2}","DileptonMass", "Selection_HWW160Selection_DileptonMass",0,160,1.55, -999, 50 );
    MakePlot(0, 4, "#Delta#phi(l,l)","degrees","DeltaPhi", "Selection_HWW160Selection_DeltaPhi",0,160,1.55 , -999, 60);
    MakePlot(0, 4, "m_{T}^{ll E_{T}^{miss}}","GeV/c^{2}","Mt", "Selection_HWW160Selection_Mt",0,160,1.55 , 90, 160);
    MakePlot(0, 4, "p_{T}^{ll}","GeV/c","DileptonPt", "Selection_HWW160Selection_DileptonPt",0,160,1.55);
    MakePlot(0, 4, "#Delta#eta(j,j)","","DeltaEtaVBFJets", "Selection_HWW160Selection_DeltaEtaVBFJets",1,160,1.55 , 3.5, -999);
    MakePlot(0, 4, "M_{jj}","GeV/c^{2}","DijetMassVBFJets", "Selection_HWW160Selection_DijetMassVBFJets",1,160,1.55 , 450, -999);
    
  }

}
