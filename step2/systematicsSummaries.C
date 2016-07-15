#include "TH1D.h"
#include "TFile.h"
#include "TAttFill.h"
#include "TLegend.h"
#include "TLine.h"
#include "TBox.h"
#include "TColor.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TAttPad.h"
#include "TAttLine.h"
#include "TAttMarker.h" 
#include "TAxis.h"
#include "TAttAxis.h"
#include "TAttText.h"
#include <iostream>

//const int FF_Bins = 5;
//double FF_Bound[FF_Bins+1] = {60,80,100,120,140,200};

const int errorSources = 9;
const char * sources[errorSources] = {"TotUP","TotDOWN","JESTotUP","JESTotDOWN","JERTot","NoChargeCut_MCDiff","pPbPbpDiff","Ratio","ChargeCutDiff"};  
const char * prettySources[errorSources] = 
  {"Total Upwards Error",
   "Total Downwards Error",
   "Upwards JES Error",
   "Downwards JES Error",
   "JER Error",
   "MC Difference Error",
   "pPb vs Pbp Error",
   "UE Subtraction Error",
   "Jet Charge Cut Error"};

void makeSummaries(int v = 0, int UEtype=3)
{ 
  TLatex * tlat = new TLatex(0.6,0.6,"test"); 

  TFile * inf2 = TFile::Open(Form("SystematicsUE%d.root",UEtype),"read");
  TH1D * sys[25];
  for(int i = 0; i<25; i++)
  {
    if(v==6 && i<10)
    {
      sys[i] = (TH1D*)inf2->Get(Form("pp2_%s%d",sources[0],i%5));
      sys[i]->SetDirectory(0);
      sys[i]->Reset();
      continue;
    } 
    if(i<5)        sys[i] = (TH1D*)inf2->Get(Form("pp2_%s%d",sources[v],i%5));
    else if(i<10)  sys[i] = (TH1D*)inf2->Get(Form("pp7_%s%d",sources[v],i%5));
    else if(i<15)  sys[i] = (TH1D*)inf2->Get(Form("pPb5_%s%d",sources[v],i%5));
    else if(i<20)  sys[i] = (TH1D*)inf2->Get(Form("interp_%s%d",sources[v],i%5));
    else           sys[i] = (TH1D*)inf2->Get(Form("FFratio_%s%d",sources[v],i%5));
    sys[i]->Scale(100);
    sys[i]->SetDirectory(0); 
  }
  inf2->Close();

  double max = 0;
  double min = 100000;
  for(int i=0; i<25; i++)
  {    
    for(int j=0; j<sys[i]->GetSize(); j++) 
    {
      if(sys[i]->GetBinContent(j)<0) sys[i]->SetBinContent(j,-sys[i]->GetBinContent(j));
      //sys[i]->SetBinError(j,0);
    }
    for(int j=1; j<sys[i]->GetSize()-1; j++) 
    {
      if(max<sys[i]->GetBinContent(j)) max = sys[i]->GetBinContent(j);
      if(min>sys[i]->GetBinContent(j) && sys[i]->GetBinContent(j)!=0) min = sys[i]->GetBinContent(j);
    }
  }

  TCanvas * c = new TCanvas("sysSummaries","sysSummaries",1200,1200);
  c->SetLeftMargin(0.23);
  c->Divide(5,5,0,0);

  for(int i = 0; i<25; i++)
  {
    c->cd(i+1)->SetLogx();
    c->cd(i+1)->SetLogy();

    sys[i]->SetMaximum(1.5*max);
    sys[i]->SetMinimum(0.5*min);
   
    /*if(v==6 && (i<15||i>19))
    {
      sys[i]->SetMaximum(1.5*max);
      sys[i]->SetMinimum(100000*min);
    }*/
 
    sys[i]->SetLineWidth(1);
    sys[i]->SetLineColor(1);
    sys[i]->SetMarkerStyle(8);
    sys[i]->SetMarkerSize(0.7);
    sys[i]->SetMarkerColor(1);

    if(i%5==0)
    {
      sys[i]->GetYaxis()->SetTitle("|Relative Syst. Error| (%)");
      sys[i]->GetYaxis()->SetTitleOffset(1.6);
      sys[i]->GetYaxis()->SetTitleSize(0.07);
      sys[i]->GetYaxis()->SetLabelOffset(0.003);
      sys[i]->GetYaxis()->SetLabelSize(0.06);
      sys[i]->GetYaxis()->SetNdivisions(505,true);
      tlat->SetTextSize(0.07);
      if(i==20)
      {
        sys[i]->GetYaxis()->SetLabelSize(0.054);
        sys[i]->GetYaxis()->SetTitleSize(0.062);
        sys[i]->GetYaxis()->SetTitleOffset(1.8);
        
        tlat->SetTextSize(0.057);
      }
    } 
    else 
    {
      sys[i]->GetYaxis()->SetTitleSize(0);
      sys[i]->GetYaxis()->SetLabelSize(0);
      tlat->SetTextSize(0.07);
    } 
    if(i>4)
    {
      sys[i]->GetXaxis()->SetTitle("p_{T}^{track} (GeV/c)");
      sys[i]->GetXaxis()->CenterTitle();
      sys[i]->GetXaxis()->SetLabelOffset(0.0001);
      sys[i]->GetXaxis()->SetLabelSize(0.07);
      sys[i]->GetXaxis()->SetTitleOffset(1.0);
      sys[i]->GetXaxis()->SetTitleSize(0.07);
      if(i==20)
      {
        sys[i]->GetXaxis()->SetLabelSize(0.058);
        sys[i]->GetXaxis()->SetTitleSize(0.055);
        sys[i]->GetXaxis()->SetTitleOffset(1.24);
        sys[i]->GetXaxis()->SetLabelOffset(0.012);
      }
    }
    
    sys[i]->Draw("hist p");
    tlat->SetNDC();
    double xmin = 0.05;
    if(i%5==0) xmin = 0.3;
    if(i<5) tlat->DrawLatex(xmin,0.9,Form("%d < p_{T}^{jet} < %d GeV/c",(int)FF_Bound[i],(int)FF_Bound[i+1]));
    if(i==1) tlat->DrawLatex(xmin,0.6,"ak3PF jets");
    if(i==1) tlat->DrawLatex(xmin,0.7,"|#eta_{CM}^{jet}|<1.5");
    if(i==1) tlat->DrawLatex(xmin,0.8,"CMS Preliminary");
    if(i==0) tlat->DrawLatex(xmin,0.8,"2.76 TeV pp FF");
    if(i==5) tlat->DrawLatex(xmin,0.8,"7 TeV pp FF");
    if(i==10) tlat->DrawLatex(xmin,0.8,"5.02 TeV pPb FF");
    if(i==15) tlat->DrawLatex(xmin,0.8,"pp Interpolation");
    if(i==20) tlat->DrawLatex(xmin,0.8,"FF ratio");
    if(i==2) tlat->DrawLatex(xmin,0.8,prettySources[v]);
  }

  c->SaveAs(Form("plots/systematicsSummaries_UE%d_%s.png",UEtype,sources[v]));
  c->SaveAs(Form("plots/systematicsSummaries_UE%d_%s.pdf",UEtype,sources[v]));
}

void systematicSummaries()
{
  for(int v=0; v<errorSources ; v++) 
  {
    if(v!=7) makeSummaries(v,3);
  }
  for(int v=0; v<errorSources ; v++) makeSummaries(v,0);
}
