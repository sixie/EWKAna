//root -l EWKAna/Hww/Selection/PlotFakePredictions.C+

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "TLegend.h"
#endif

void PlotFakePredictions() {

    TFile *fileSS = new TFile("WWSelectionPlotsFakePrediction_SS.root", "READ");
    TFile *fileOS = new TFile("WWSelectionPlotsFakePrediction.root", "READ");
    
    TH1F *LeptonPtMin_OS = (TH1F*)fileOS->Get("hLeptonPtMin");
    TH1F *LeptonPtMin_SS = (TH1F*)fileSS->Get("hLeptonPtMin");

    TLegend *legend = new TLegend(0.73,0.75,0.93,0.90);
    legend->SetTextSize(0.03);
    legend->SetBorderSize(1);
    
    legend->Clear();
    legend->AddEntry(LeptonPtMin_OS, "OppSign", "L");
    legend->AddEntry(LeptonPtMin_SS, "SameSign", "LP");

    TCanvas *cv = new TCanvas("cv","cv",800,600);
    LeptonPtMin_OS->SetMaximum(0.4);
    LeptonPtMin_OS->Draw("hist");
    LeptonPtMin_SS->SetLineColor(kRed);
    LeptonPtMin_SS->SetMarkerColor(kRed);
    LeptonPtMin_SS->Draw("same,E1");
    legend->Draw();
    cv->SaveAs("LeptonPtMin_FakePrediction.gif");

//     fileSS->Close();
//     fileOS->Close();



    TH1F *DileptonMass_OS = (TH1F*)fileOS->Get("dileptonMass");
    TH1F *DileptonMass_SS = (TH1F*)fileSS->Get("dileptonMass");
    DileptonMass_OS->SetMaximum(0.25);
    DileptonMass_OS->Draw("hist");
    DileptonMass_SS->SetLineColor(kRed);
    DileptonMass_SS->SetMarkerColor(kRed);
    DileptonMass_SS->Draw("same,E1");
    legend->Draw();
    cv->SaveAs("DileptonMass_FakePrediction.gif");

    TH1F *DeltaPhiLeptons_OS = (TH1F*)fileOS->Get("hDeltaPhiLeptons");
    TH1F *DeltaPhiLeptons_SS = (TH1F*)fileSS->Get("hDeltaPhiLeptons");
    DeltaPhiLeptons_OS->SetMaximum(0.25);
    DeltaPhiLeptons_OS->Draw("hist");
    DeltaPhiLeptons_SS->SetLineColor(kRed);
    DeltaPhiLeptons_SS->SetMarkerColor(kRed);
    DeltaPhiLeptons_SS->Draw("same,E1");
    legend->Draw();
    cv->SaveAs("DeltaPhiLeptons_FakePrediction.gif");


}
