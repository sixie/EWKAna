
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#endif

//=== MAIN MACRO =================================================================================================

void makeTriggerEfficiency(string DataFilename = "Data_EleWPEffTP/basic2_76_106/eff.root",
                   string outputDir = "Data_EleWPEffTP/basic2_76_106/", 
                   string histName = "h2_results_electron_selection",
                   string Label = "SmurfV6")
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  TString outfname = (outputDir + "/eff_table.txt").c_str();
  TFile datafile2(DataFilename.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TH2F *hDataEff2=0, *hDataErrl2=0, *hDataErrh2=0;
  
  hDataEff2  = (TH2F*)datafile2.Get("hEffEtaPt");
  hDataErrl2 = (TH2F*)datafile2.Get("hErrlEtaPt");
  hDataErrh2 = (TH2F*)datafile2.Get("hErrhEtaPt");
  

  //--------------------------------------------------------------------------------------------------------------
  // Update root file histograms
  //==============================================================================================================   
  const Int_t nx = hDataEff2->GetNbinsX();
  const Int_t ny = hDataEff2->GetNbinsY();

  TFile *outputFile = new TFile(("efficiency_results"+label+".root").c_str(), "UPDATE");

  //Do Binning
  Double_t *ptbins = new Double_t[ny+1];
  Double_t *etabins = new Double_t[nx+1];
  for(Int_t iy=1; iy<=ny; iy++) {
      ptbins[iy-1] = hDataEff2->GetYaxis()->GetBinLowEdge(iy);
  }
  for(Int_t ix=1; ix<=nx; ix++) {
    etabins[ix-1]= hDataEff2->GetXaxis()->GetBinLowEdge(ix);
  }
  ptbins[ny] = 60;
  etabins[nx] = 2.5;
    
  TH2F *h2_results_selection = new TH2F(histName.c_str(),"",ny,ptbins,nx,etabins);
  for(Int_t iy=0; iy<=ny+2; iy++) {
    for(Int_t ix=0; ix<=nx+2; ix++) {
      h2_results_selection->SetCellContent(iy,ix, 1.0);
      h2_results_selection->SetCellError(iy,ix, 0.0);
    }
  }

  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Produce Text file table
  //==============================================================================================================   

  ofstream txtfile;
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    Trigger Efficiency         ";  
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      txtfile << "[" << setw(4) << hDataEff2->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hDataEff2->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hDataEff2->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hDataEff2->GetXaxis()->GetBinLowEdge(ix+1) << "]";      
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      
      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      txtfile << endl;
      txtfile.flags(flags);

      cout << "Set " << iy << " " << ix << " " << dataeff << endl;
      h2_results_selection->SetCellContent(iy,ix, dataeff);
      h2_results_selection->SetCellError(iy,ix, (dataerrl + dataerrh)/2);
      //fill overflow bins with the same values as last bin
      if (ix == nx) {
      cout << "Set " << iy << " " << ix+1 << " " << dataeff << endl;
        h2_results_selection->SetCellContent(iy,nx+1, dataeff);
        h2_results_selection->SetCellError(iy,nx+1, (dataerrl + dataerrh)/2);
      }
    }
    if (iy == ny) {
      for(Int_t ix=1; ix<=nx; ix++) {
       Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
        Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
        Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
        
        cout << "Set " << iy+1 << " " << ix << " " << dataeff << endl;
        h2_results_selection->SetCellContent(iy+1,ix, dataeff);
        h2_results_selection->SetCellError(iy+1,ix, (dataerrl + dataerrh)/2);
        //fill overflow bins with the same values as last bin
        if (ix == nx) {
          cout << "Set " << iy+1 << " " << ix+1 << " " << dataeff << endl;
          h2_results_selection->SetCellContent(iy+1,nx+1, dataeff);
          h2_results_selection->SetCellError(iy+1,nx+1, (dataerrl + dataerrh)/2);
        }
      }
    }

    txtfile << endl;
  }
  txtfile.close();
  
  cout << outfname << " created!" << endl;
  



//   //--------------------------------------------------------------------------------------------------------------
//   // Create TEX table
//   //==============================================================================================================   

//   ofstream texfile;
//   texfile.open((outputDir + "/eff_table.tex").c_str());
//   assert(texfile.is_open());

//   texfile << " \\begin{table}[!ht]" << endl;
//   texfile << " \\begin{center} " << endl;
//   texfile << " \\begin{tabular}{|c|c|}" << endl;
//   texfile << " \\hline\n";


//   texfile << " $p_{T}$ / $\\eta$ bin    &  Trigger Efficiency   \\\\  ";
//   texfile << " \\hline           ";
//   texfile << endl;
//   for(Int_t iy=1; iy<=ny; iy++) {
//     for(Int_t ix=1; ix<=nx; ix++) {

//       string binLabel = Form("$%5.1f < p_{T} \\le %5.1f$ , $%5.1f  \\le |\\eta| < %5.1f$", 
//                              hDataEff2->GetYaxis()->GetBinLowEdge(iy), hDataEff2->GetYaxis()->GetBinLowEdge(iy+1),
//                              hDataEff2->GetXaxis()->GetBinLowEdge(ix), hDataEff2->GetXaxis()->GetBinLowEdge(ix+1));
//       if (iy == ny) {
//         binLabel = Form("$%5.1f < p_{T} $ , $%5.1f  \\le |\\eta| < %5.1f$", 
//                         hDataEff2->GetYaxis()->GetBinLowEdge(iy), 
//                         hDataEff2->GetXaxis()->GetBinLowEdge(ix), hDataEff2->GetXaxis()->GetBinLowEdge(ix+1));
//       }
      
//       texfile << binLabel;
//       texfile << "   &   ";

//       ios_base::fmtflags flags = texfile.flags();
//       texfile.precision(4);
      
//       Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
//       Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
//       Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
//       texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
//       texfile << "   \\\\   ";
//       texfile << endl;
//       texfile << "\\hline";
//       texfile << endl;

//     }
//   }
//   texfile << "\\end{tabular}" << endl;
//   texfile << "\\caption{CAPTION.}" << endl;
//   texfile << "\\label{tab:eff_xxx_trigger}" << endl;
//   texfile << "\\end{center}" << endl;
//   texfile << "\\end{table}" << endl;
//   texfile.close();
//   cout << outputDir + "/eff_table.tex" << " created!" << endl;



  //--------------------------------------------------------------------------------------------------------------
  // Create TEX table
  //==============================================================================================================   

  ofstream texfile;
  texfile.open((outputDir + "/eff_table.tex").c_str());
  assert(texfile.is_open());

  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center} " << endl;
  texfile << " \\begin{tabular}{|c|";
  for(Int_t ix=1; ix<=nx; ix++) {
    texfile << "c|";
  }
  texfile << "}" << endl;  
  texfile << " \\hline\n";

  texfile << " Measurement &  ";

  for(Int_t ix=1; ix<=nx; ix++) {
    string tmp = Form("$%5.1f  \\le |\\eta| < %5.1f$",
                      hDataEff2->GetXaxis()->GetBinLowEdge(ix), 
                      hDataEff2->GetXaxis()->GetBinLowEdge(ix+1));
    texfile << tmp ;
    if (ix < nx) texfile << " & ";
    else texfile << " \\\\ " << endl;
  }
  texfile << " \\hline           " << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    
    string tmp;
    if (iy < ny) tmp = Form("$%5.1f < p_{T} \\le %5.1f$",
                            hDataEff2->GetYaxis()->GetBinLowEdge(iy), 
                            hDataEff2->GetYaxis()->GetBinLowEdge(iy+1));
    else tmp = Form("$%5.1f < p_{T} $",
                    hDataEff2->GetYaxis()->GetBinLowEdge(iy));
    texfile << tmp << " & ";

    for(Int_t ix=1; ix<=nx; ix++) {      

      ios_base::fmtflags flags = texfile.flags();
      texfile.precision(4);
      
      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      if (ix < nx) texfile << "  &  ";
      else texfile << " \\\\ " << endl << "\\hline" << endl;
    }
  }
  texfile << "\\end{tabular}" << endl;
  texfile << "\\caption{CAPTION.}" << endl;
  texfile << "\\label{tab:eff_xxx_trigger}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;
  texfile.close();
  cout << outputDir + "/eff_table.tex" << " created!" << endl;






  //--------------------------------------------------------------------------------------------------------------
  // Create TWIKI table
  //==============================================================================================================   

  ofstream twikifile;
  twikifile.open((outputDir + "/eff_table.twiki").c_str());
  assert(twikifile.is_open());
  

  twikifile << "| *pT bin* | *eta bin* | *Trigger Efficiency* | \n";
  
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      string binLabel = Form("| *%4.1f < pT \\le %4.1f* | *%3.1f  <= eta < %3.1f* | ", 
                             hDataEff2->GetYaxis()->GetBinLowEdge(iy), hDataEff2->GetYaxis()->GetBinLowEdge(iy+1),
                             hDataEff2->GetXaxis()->GetBinLowEdge(ix), hDataEff2->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("| *%4.1f < pT* | *%3.1f  <= eta < %3.1f* | ", 
                        hDataEff2->GetYaxis()->GetBinLowEdge(iy), 
                        hDataEff2->GetXaxis()->GetBinLowEdge(ix), hDataEff2->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      twikifile << binLabel ;

      ios_base::fmtflags flags = twikifile.flags();
      twikifile.precision(4);
      
      Double_t dataeff  = hDataEff2->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl2->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh2->GetCellContent(ix,iy);
      twikifile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      twikifile << " |\n";

    }
  }
  twikifile.close();
  cout << outputDir + "/eff_table.twiki" << " created!" << endl;


  outputFile->WriteTObject(h2_results_selection, h2_results_selection->GetName(), "WriteDelete");
  outputFile->Close();

  gBenchmark->Show("printMuonWPEff"); 
}
