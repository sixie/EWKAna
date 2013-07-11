//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; b<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}

void PlotPileupTargets() {

  TFile *f = new TFile("/data/smurf/sixie/Pileup/PUTarget.Run2011A.160404-173692.root");
  TH1F* Target_Run2011A = (TH1F*)f->Get("pileup");
  f = new TFile("/data/smurf/sixie/Pileup/PUTarget.Full2011.160404-180252.root");
  TH1F* Target_Full2011 = (TH1F*)f->Get("pileup");
  f = new TFile("/data/smurf/sixie/Pileup/PUTarget.Run2011B.175832-180252.root");
  TH1F* Target_Run2011B = (TH1F*)f->Get("pileup");

  NormalizeHist(Target_Run2011A);
  NormalizeHist(Target_Run2011B);
  NormalizeHist(Target_Full2011);


  TCanvas *cv = new TCanvas("cv","cv", 800,600);

  TLegend *tmpLegend = new TLegend(0.4,0.75,0.93,0.90);   
  tmpLegend->SetTextSize(0.04);
  tmpLegend->SetBorderSize(0);
  tmpLegend->SetFillStyle(0);

  tmpLegend->AddEntry(Target_Run2011A, "First Half of 2011 Data" , "L");
  tmpLegend->AddEntry(Target_Run2011B, "Second Half of 2011 Data" , "L");
  tmpLegend->AddEntry(Target_Full2011, "Full 2011 Data" , "L");

  Target_Run2011A->SetTitle("");
  Target_Run2011A->GetXaxis()->SetTitle("Number of Pileup Events");
  Target_Run2011A->GetYaxis()->SetTitle("Fraction of Events");
  Target_Run2011A->GetYaxis()->SetTitleOffset(1.4);
  Target_Run2011A->GetXaxis()->SetRangeUser(0,35);
  Target_Run2011A->SetLineColor(kRed);
  Target_Run2011B->SetLineColor(kBlack);
  Target_Full2011->SetLineColor(kBlue);
  Target_Run2011A->SetLineWidth(2);
  Target_Run2011B->SetLineWidth(2);
  Target_Full2011->SetLineWidth(2);

  Target_Run2011A->Draw("hist");
  Target_Run2011B->Draw("hist,same");
  Target_Full2011->Draw("hist,same");

  tmpLegend->Draw();

  cv->SaveAs("PUTargetComparison.gif");
  cv->SaveAs("PUTargetComparison.eps");

  }
