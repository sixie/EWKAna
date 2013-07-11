#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TLorentzVector.h>

#include "EWKAna/Ntupler/interface/TGenInfo.hh"
#include "EWKAna/Ntupler/interface/TGenInfoExt.hh"
#endif

//#define __ZMM__
#define __ZEE__

void ResBosNtupler(const TString  resbosfile,   // ResBos output file
                   const Int_t    qid,          // quark ID (1 = down, 2 = up)
		   const Double_t xsec,         // cross section from ResBos
		   const Double_t nevents,      // number of events
		   const TString  outfilename)  // output file name
{
  gBenchmark->Start("ResBosNtupler");

  TTree::SetMaxTreeSize(kMaxLong64);
  
  // Don't write TObject part of the objects  
#ifdef __ZMM__
  mithep::TGenInfo::Class()->IgnoreTObjectStreamer();
  mithep::TGenInfo *gen = new mithep::TGenInfo();
#endif  //__ZMM__
#ifdef __ZEE__
  mithep::TGenInfoExt::Class()->IgnoreTObjectStreamer();
  mithep::TGenInfoExt *gen = new mithep::TGenInfoExt();
#endif  //__ZEE__
  
  //
  // Initialize output file, data trees, and structs
  // 
  TFile* outfile = new TFile(outfilename, "RECREATE");
  TTree *outEventTree = new TTree("Events","Events"); 
  outEventTree->Branch("Gen", &gen);

  
  TLorentzVector boson, lep1, lep2, pho, dilep;
  
  // 
  // parse ResBos output file
  //  
  ifstream ifs;
  ifs.open(resbosfile.Data()); 
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {  
    Int_t evtnum;
    Double_t aqcd;
    stringstream ss1(line);
    ss1 >> evtnum >> aqcd;
    if(evtnum==0) continue;
    
    getline(ifs,line);
    stringstream ss2(line);
    Double_t var1, var2, var3;
    ss2 >> var1 >> var2 >> var3;
    
    UInt_t npho=0;
    Bool_t evtend = kFALSE;
    while(!evtend) {
      getline(ifs,line);
      Int_t pid, foo1, foo2;
      Double_t px, py, pz, e;
      stringstream ss(line);
      ss >> pid;
      if(pid==0) {
        evtend = kTRUE;
      } else {
        ss >> px >> py >> pz >> e >> foo1 >> foo2;
      }
      if(pid==90) {  // Z boson kinematics
        boson.SetPxPyPzE(px,py,pz,e);
      
      } else if((pid==12) || (pid==14)) {  // electron / muon kinematics
        lep1.SetPxPyPzE(px,py,pz,e);
      
      } else if((pid==-12) || (pid==-14)) {  // positron / anti-muon kinematics
        lep2.SetPxPyPzE(px,py,pz,e);
       
      } else if(pid==10) {  // FSR photon kinematics
        npho++;
	pho.SetPxPyPzE(px,py,pz,e);
      }
    }
    dilep = lep1 + lep2;
  
    gen->nGenPart = 0;
    gen->nGenCh   = 0;
    gen->npho     = npho;
    gen->id_1     = qid;
    gen->id_2     = -qid;
    gen->x_1      = 0;
    gen->x_2      = 0;
    gen->weight   = xsec/nevents;
    gen->vmass    = boson.M();
    gen->vpt      = boson.Pt();
    gen->vy       = boson.Rapidity();
    gen->vphi     = boson.Phi();
    gen->mass     = dilep.M();
    gen->pt       = dilep.Pt();
    gen->y        = dilep.Rapidity();
    gen->phi      = dilep.Phi();
    gen->pt_1     = lep1.Pt();
    gen->eta_1    = lep1.Eta();
    gen->phi_1    = lep1.Phi();
    gen->pt_2     = lep2.Pt();
    gen->eta_2    = lep2.Eta(); 
    gen->phi_2    = lep2.Phi();
    gen->phopt    = (npho) ? pho.Pt()  : -999;
    gen->phoeta   = (npho) ? pho.Eta() : -999;
    gen->phophi   = (npho) ? pho.Phi() : -999;
    gen->decx     = 0;
    gen->decy     = 0; 
    gen->decz     = 0;
#ifdef __ZEE__
    gen->scEt_1   = -1;
    gen->scEta_1  = -999;
    gen->scEt_2   = -1;
    gen->scEta_2  = -999;
    gen->scMass   = -1;
#endif  //__ZEE__    
    outEventTree->Fill();
  }
  ifs.close();

  outfile->Write();
  outfile->Close();
  
  delete gen;
  delete outfile;  
  
  std::cout << outfilename << " created!" << std::endl;
  
  gBenchmark->Stop("ResBosNtupler");
}
