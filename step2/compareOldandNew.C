#include "TFile.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TLatex.h"

void compareOldandNew()
{
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  TLatex lat;  

  int centBounds[6] = {60,80,100,120,140,200};
  for(int i = 0 ; i<5; i++){
  const char * histName = "pPbPbp_FF";
  TFile * fold = TFile::Open("FragmentationFunctionsUE2_nominal.root","read");
  TFile * fnew = TFile::Open("FragmentationFunctionsUE2.root","read");
  TH1D * num = (TH1D*)fnew->Get(Form("%s_%d_%d",histName,centBounds[i],centBounds[i+1]));
  TH1D * den = (TH1D*)fold->Get(Form("%s_%d_%d",histName,centBounds[i],centBounds[i+1]));
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(no L2 residual)/(Nominal)");
  num->GetYaxis()->SetRangeUser(0.7,1.3);
  num->GetXaxis()->SetTitle("p_{T}");
  num->GetXaxis()->SetRangeUser(0.5,200);
  num->Print("All");
  num->Draw();
  
  lat.DrawLatex(1,1.2,Form("%d < p_{T}^{jet} < %d",centBounds[i],centBounds[i+1]));

  c1->SetLogx();
  c1->SaveAs(Form("diffPlots/noL2Residual_%d_%d.png",centBounds[i],centBounds[i+1]));
  c1->SaveAs(Form("diffPlots/noL2Residual_%d_%d.pdf",centBounds[i],centBounds[i+1]));
  fnew->Close();
  fold->Close();
  }  

  /*const char * histName = "PbPbTrackSpectrum_0_5";
  TFile * fnew = TFile::Open("Spectra_Jun9_noChi2Cut.root","read");
  TFile * fold = TFile::Open("Spectra_Jun9_withChi2Cut.root","read");
  TH1D * num = (TH1D*)fnew->Get(histName);
  TH1D * den = (TH1D*)fold->Get(histName);
  //den->Add(num);
  num->Divide(den);
  //num->Scale(2);
  num->GetYaxis()->SetTitle("(No Chi2)/(With Chi2)");
  num->GetYaxis()->SetRangeUser(0.7,1.5);
  num->GetXaxis()->SetRangeUser(0.7,350);
  num->Print("All");
  num->Draw();
  c1->SetLogx();
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.png");
  c1->SaveAs("plots/comparisonPlots/Chi2CutTest_PbPb_0_5.pdf");
  */

  //c1->SaveAs("plots/comparisonPlots/ppChargeFraction2.C");
}
