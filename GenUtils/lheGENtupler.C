#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TLorentzVector.h>
#endif

#include "GenStructDefs.hh"

TGenData data;

//====================================================================================================

void lheGENtupler()
{
  gBenchmark->Start("lheGENtupler");
  
//   const char* outfname = "/data/blue/sixie/Madgraph/WW/ww_scaleUp_gentuple.root";
  const char* outfname = "ww_gentuple.cteq6m.root";

  vector<string> fnamev;
  fnamev.push_back("/data/blue/sixie/MCGenerators/Madgraph/Wgammastar/MODLLLnoZ_ee1.lhe");

  TFile* outfile = new TFile(outfname,"RECREATE");

  TTree::SetMaxTreeSize(kMaxLong64);
  
  TTree* tree = new TTree("Events","Events");  
  tree->Branch("Events",&data.npho,
"npho/i:scalePdf/F:weight:vmass:vpt:vy:vphi:mass:pt:y:phi:mt:pt_1:eta_1:phi_1:pt_2:eta_2:phi_2:q1/I:q2:phopt/F:phoeta:phophi:x_1:x_2:id_1/I:id_2:acqd/F:aqed:jetpt:jeteta:jetphi:jetid/I"); 
    
  const Int_t bosonID =-24;
  const Int_t lepID1  = 13;
  const Int_t lepID2  =-14;
  
  //--------------------------------------------------------------------------------------------------
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
  
    const char* lhefile = fnamev[ifile].c_str();
    ifstream ifs;
    ifs.open(lhefile);
    assert(ifs.is_open()); 
  
    string line; 
        
    cout << "Processing file: " << fnamev[ifile] << endl;
    
    do {
      getline(ifs,line);
    } while(line.compare("<init>")!=0);
    
    getline(ifs,line);
    stringstream ss1(line);
    int idbmup1, idbmup2; 
    double ebmup1, ebmup2; 
    int pdfgup1, pdfgup2, pdfsup1, pdfsup2, idwtup, nprup;
    ss1 >> idbmup1 >> idbmup2 >> ebmup1 >> ebmup2 >> pdfgup1 >> pdfgup2 >> pdfsup1 >> pdfsup2 >> idwtup >> nprup;
    
    getline(ifs,line);
    stringstream ss2(line);
    double xsecup, xerrup, xmaxup;
    int lprup;
    ss2 >> xsecup >> xerrup >> xmaxup >> lprup;
    getline(ifs,line);
    ss2 >> xsecup >> xerrup >> xmaxup >> lprup;
//     getline(ifs,line);
//     ss2 >> xsecup >> xerrup >> xmaxup >> lprup;
    
    getline(ifs,line);
    assert(line.compare("</init>")==0);
    
    while(getline(ifs,line)) {

      if(line.compare("<event>")==0) {
        getline(ifs,line);
	stringstream ss3(line);
	int nup, idprup;
	double xwgtup, scalup, aqedup, aqcdup;
	ss3 >> nup >> idprup >> xwgtup >> scalup >> aqedup >> aqcdup;
	data.scalePdf = scalup;
	data.weight   = xwgtup;
	data.aqcd     = aqcdup;
	data.aqed     = aqedup;
			
	int idup, istup, mothup1, mothup2, icolup1, icolup2; 
	double  pup1, pup2, pup3, pup4, pup5, vtimup, spinup;
	
	// quark 1 info
        getline(ifs,line);
        stringstream ssq1(line);
	ssq1 >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
	data.x_1  = fabs(pup3/ebmup1);
	data.id_1 = idup;
	
	// quark 2 info
        getline(ifs,line);
        stringstream ssq2(line);
	ssq2 >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
	data.x_2  = fabs(pup3/ebmup2);
	data.id_2 = idup;
	
	TLorentzVector vec1;
	TLorentzVector vec2;
	
        double leadpt = 0;
        data.npho   = 0;
        data.phopt  = -999;
        data.phoeta = -999;
        data.phophi = -999;
	
	getline(ifs,line);	
	while(line.compare("</event>")!=0) {
	  stringstream ss(line);
	  ss >> idup >> istup >> mothup1 >> mothup2 >> icolup1 >> icolup2 >> pup1 >> pup2 >> pup3 >> pup4 >> pup5 >> vtimup >> spinup;
	  
	  if(idup==bosonID) {  // boson info
	    TLorentzVector boson;
	    boson.SetPxPyPzE(pup1,pup2,pup3,pup4);
	    data.vmass = boson.M();
	    data.vpt   = boson.Pt();
	    data.vy    = boson.Rapidity();
	    data.vphi  = boson.Phi();
	  
	  } else if(idup==lepID1) {  // lepton 1 info
	    vec1.SetPxPyPzE(pup1,pup2,pup3,pup4);
	    data.q_1   = -1;
            data.pt_1  = vec1.Pt();
            data.eta_1 = vec1.Eta();
            data.phi_1 = vec1.Phi(); 
	    
	  } else if(idup==lepID2) {  // lepton 2 info
	    vec2.SetPxPyPzE(pup1,pup2,pup3,pup4);
	    data.q_2   = 1;
            data.pt_2  = vec2.Pt();
            data.eta_2 = vec2.Eta();
            data.phi_2 = vec2.Phi(); 
	  
	  } else if(idup==22) {  // photon info
	    data.npho++;
	    TLorentzVector vecpho;
	    vecpho.SetPxPyPzE(pup1,pup2,pup3,pup4);
            if( leadpt < vecpho.Pt() ) {        
	      leadpt = vecpho.Pt();
	      data.phopt  = vecpho.Pt();
	      data.phoeta = vecpho.Eta();
	      data.phophi = vecpho.Phi();
            }
	  } else {
	    TLorentzVector jet;
	    jet.SetPxPyPzE(pup1,pup2,pup3,pup4);
	    data.jetid  = idup;
	    data.jetpt  = jet.Pt();
	    data.jeteta = jet.Rapidity();
	    data.jetphi = jet.Phi();
	  }
	  
	  getline(ifs,line);
	}
	
	// dilepton info
	TLorentzVector dilep = vec1 + vec2;
	data.mass = dilep.M();
	data.pt   = dilep.Pt();
	data.y    = dilep.Rapidity();
	data.phi  = dilep.Phi();
	
	double et1 = sqrt(vec1.Perp2() + vec1.M2());
        double et2 = sqrt(vec2.Perp2() + vec2.M2());
        data.mt = sqrt( (et1+et2)*(et1+et2) - dilep.Pt()*dilep.Pt() );
	
	tree->Fill();
      }  
    }
    ifs.close();
  }
    
  tree->Print();
  outfile->Write();  
  outfile->Close();
    
  gBenchmark->Show("lheGENtupler");
}
