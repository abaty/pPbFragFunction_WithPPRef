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

const int FF_Bins = 5;
double FF_Bound[FF_Bins+1] = {60,80,100,120,140,200};

TBox* makeBox(double x1, double y1, double x2, double y2)
{
  TBox * box = new TBox(x1,y1,x2,y2);
  box->SetLineColor(1);
  box->SetLineWidth(1);
  box->SetLineStyle(1);
  box->SetFillStyle(1001);
  box->SetFillColor(kOrange);
  return box;
}

void interpRatio(int UEtype, int turnOffPoints = 0)
{  
  TLatex * tlat = new TLatex(0.6,0.6,"test");
  tlat->SetTextFont(42);

  TFile * inf1 = TFile::Open(Form("../FragmentationFunctionsUE%d.root",UEtype),"read");
  TH1D * hist[15];
  for(int i = 0; i<15; i++)
  {
    if(i<5)       hist[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    else if(i<10) hist[i] = (TH1D*)inf1->Get(Form("ppref5TeV_data_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    else          hist[i] = (TH1D*)inf1->Get(Form("pPbPbp_FF_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    hist[i]->SetDirectory(0);
  }
  inf1->Close();

  TFile * inf2 = TFile::Open(Form("../SystematicsUE%d.root",UEtype),"read");
  TH1D * sysU[15];
  for(int i = 0; i<15; i++)
  {
    if(i<5)
    {     
      sysU[i] = (TH1D*)inf2->Get(Form("pPbSyst_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    }   
    else if(i<10) 
    {
      sysU[i] = (TH1D*)inf2->Get(Form("ppSyst_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    }    
    else 
    { 
      sysU[i] = (TH1D*)inf2->Get(Form("ratSyst_%d_%d",(int)FF_Bound[i%5],(int)FF_Bound[i%5+1]));
    }
    sysU[i]->SetDirectory(0);
  }
  inf2->Close();

  TBox* boxArray[15][30];
  for(int i=0; i<15; i++)
  {
    for(int j=1; j<hist[i]->GetSize()-1; j++)
    {
      boxArray[i][j] = makeBox(hist[i]->GetXaxis()->GetBinLowEdge(j),hist[i]->GetBinContent(j)*(1-sysU[i]->GetBinContent(j)),hist[i]->GetXaxis()->GetBinUpEdge(j),hist[i]->GetBinContent(j)*(1+sysU[i]->GetBinContent(j)));
    }
  }

  TCanvas * c = new TCanvas("rawFFs","rawFFs",1200,600);
  c->SetLeftMargin(0.23);
  c->Divide(5,2,0,0);
  
  TCanvas * c1 = new TCanvas("rawFFs1","rawFFs1",1200,350);
  c1->SetLeftMargin(0.23);
  c1->Divide(5,1,0,0);

  for(int i = 0; i<15; i++)
  { 
    if(i<10) c->cd(i+1)->SetLogx();
    if(i<10) c->cd(i+1)->SetLogy();
    if(i>=10) c1->cd(i-10+1)->SetLogx();

    hist[i]->SetMaximum(9);
    hist[i]->SetMinimum(0.0000012);
    
    if(i>9)
    {
      hist[i]->SetMaximum(2.9);
      hist[i]->SetMinimum(0.2);
    }

    hist[i]->SetLineWidth(1);
    hist[i]->SetLineColor(1);
    if(i<5)
    { 
      hist[i]->SetMarkerStyle(8);
      hist[i]->SetMarkerSize(0.7);
    }
    else if(i<10)
    {
      hist[i]->SetMarkerStyle(24);
      hist[i]->SetMarkerSize(0.7);
    }
    else
    {
      hist[i]->SetMarkerStyle(8);
      hist[i]->SetMarkerSize(0.7);
    }
    hist[i]->SetMarkerColor(1);

    if(i%5==0)
    {
      hist[i]->GetYaxis()->SetTitle("#frac{1}{N_{jet}} #frac{dN^{trk}}{dp_{T}^{trk}} (GeV^{-1}c)");
      tlat->SetTextSize(0.068);
      if(i==10)
      {
        hist[i]->GetYaxis()->SetTitle("FF_{pPb}/FF_{pp}^{interpolation}");
        hist[i]->GetYaxis()->CenterTitle();
      tlat->SetTextSize(0.068);
      }  
      hist[i]->GetYaxis()->SetTitleOffset(1.6);
      hist[i]->GetYaxis()->SetTitleSize(0.07);
      hist[i]->GetYaxis()->SetLabelOffset(0.003);
      hist[i]->GetYaxis()->SetLabelSize(0.05);

    } 
    else 
    {
      hist[i]->GetYaxis()->SetTitleSize(0);
      hist[i]->GetYaxis()->SetLabelSize(0);

      tlat->SetTextSize(0.085);
    } 
    if(i>4)
    {
      hist[i]->GetXaxis()->SetTitle("p_{T}^{track} (GeV/c)");
      hist[i]->GetXaxis()->CenterTitle();
      hist[i]->GetXaxis()->SetLabelOffset(0.0001);
      hist[i]->GetXaxis()->SetLabelSize(0.07);
      hist[i]->GetXaxis()->SetTitleOffset(1.0);
      hist[i]->GetXaxis()->SetTitleSize(0.07);
      if(i==10 || i==5)
      {
        hist[i]->GetXaxis()->SetLabelSize(0.058);
        hist[i]->GetXaxis()->SetTitleSize(0.058);
        hist[i]->GetXaxis()->SetTitleOffset(1.2);
        hist[i]->GetXaxis()->SetLabelOffset(0.01);
      }
    }

    
    hist[i]->Draw();
    for(int j=1; j<hist[i]->GetSize()-1; j++)
    { 
      if(turnOffPoints && i==13 && j==21)
      {
        hist[i]->SetBinContent(j,0);
        hist[i]->SetBinError(j,0);
        continue;
      } 
      if(turnOffPoints && i==8 && j==21)
      {
        hist[i]->SetBinContent(j,0);
        hist[i]->SetBinError(j,0);
        continue;
      }
      if(boxArray[i][j]->GetY2()!=0) boxArray[i][j]->Draw("l");
    }
    hist[i]->Draw("same");
    if(i<5) tlat->DrawLatex(0.7,0.000008,Form("%d < p_{T}^{jet} < %d GeV/c",(int)FF_Bound[i],(int)FF_Bound[i+1]));
    if(i>9) tlat->DrawLatex(0.65,0.35,Form("%d < p_{T}^{jet} < %d GeV/c",(int)FF_Bound[i-10],(int)FF_Bound[i+1-10]));
    if(i==0) tlat->DrawLatex(1,2,"30.9 nb^{-1} (5.02 TeV pPb)");
    if(i==5) tlat->DrawLatex(1,2,"5.3 pb^{-1} (2.76 TeV pp)");
    if(i==6) tlat->DrawLatex(1.2,2,"2.5 fb^{-1} (7 TeV pp)");
    //if(i==5) tlat->DrawLatex(1.05,0.00002,"#int_{pPb}Ldt = 30.9 nb^{-1}");
    //if(i==6) tlat->DrawLatex(0.8,0.00002,"#int_{2.76 TeV pp}Ldt = 5.3 pb^{-1}");
    //if(i==7) tlat->DrawLatex(1,0.00002,"#int_{7 TeV pp}Ldt = 2.5 fb^{-1}");
    if(i==1) tlat->DrawLatex(1.2,2,"anti-k_{t} particle flow jets");
    if(i==2) tlat->DrawLatex(20,2,"R=0.3");
    if(i==3) tlat->DrawLatex(20,2,"|#eta_{CM}^{jet}|<1.5");
    if(i==4)
    {
      tlat->SetTextFont(61);
      tlat->SetTextSize(1.0/0.76*tlat->GetTextSize());
      tlat->DrawLatex(36,1.8,"CMS");
      tlat->SetTextFont(52);
      tlat->SetTextSize(0.76*tlat->GetTextSize());
      tlat->DrawLatex(14,0.7,"Preliminary");
      tlat->SetTextFont(42);
    }
    if(i==7 && UEtype==3)
    {
      tlat->DrawLatex(2,2,"No Underlying Event");
      tlat->DrawLatex(12,0.7,"Subtraction");
    }
    if(i==11 && UEtype==3)
    {
      tlat->DrawLatex(1,1.9,"No Underlying Event");
      tlat->DrawLatex(2,1.7,"Subtraction");
    }

    if(i==10) tlat->DrawLatex(0.65,2.35,"30.9 nb^{-1} (5.02 TeV pPb)");
    if(i==10) tlat->DrawLatex(0.65,2.6,"5.3 pb^{-1} (2.76 TeV pp)");
    if(i==10) tlat->DrawLatex(0.65,2.1,"2.5 fb^{-1} (7 TeV pp)");
    if(i==11) tlat->DrawLatex(1,2.6,"anti-k_{t} particle flow jets");
    if(i==12) tlat->DrawLatex(1,2.6,"R=0.3");
    if(i==13) tlat->DrawLatex(1,2.6,"|#eta_{CM}^{jet}|<1.5");
    if(i==14)
    {
      tlat->SetTextFont(61);
      tlat->SetTextSize(1.0/0.76*tlat->GetTextSize());
      tlat->DrawLatex(36,2.6,"CMS");
      tlat->SetTextFont(52);
      tlat->SetTextSize(0.76*tlat->GetTextSize());
      tlat->DrawLatex(14,2.4,"Preliminary");
      tlat->SetTextFont(42);
    }


    if(i>9)
    {
      TLine * l = new TLine(0.5,1, hist[i]->GetBinLowEdge(40),1);
      l->SetLineWidth(1);
      l->SetLineStyle(3);
      l->SetLineColor(1);
      l->Draw("same");
    }

    TLegend * leg;
    TLegend * leg2;
    if(i==0)
    {
      leg = new TLegend(0.28,0.23,0.7,0.5);
      leg->AddEntry(hist[0],"5.02 TeV pPb","p");
      //leg->AddEntry(hist[10],"Fragmentation Function Ratio","p");
      leg->SetTextSize(tlat->GetTextSize());
      leg->SetBorderSize(0);
      leg->SetTextFont(tlat->GetTextFont());
      leg->Draw();
    }
    if(i==5)
    {
      leg2 = new TLegend(0.28,0.4,0.7,0.5);
      leg2->AddEntry(hist[5],"pp interpolation","p");
      //leg->AddEntry(hist[10],"Fragmentation Function Ratio","p");
      leg2->SetTextSize(tlat->GetTextSize());
      leg2->SetBorderSize(0);
      leg2->SetTextFont(tlat->GetTextFont());
      leg2->Draw();
    }
  }
  c->SaveAs(Form("plots/prettyInterp_UE%d.png",UEtype));
  c->SaveAs(Form("plots/prettyInterp_UE%d.pdf",UEtype));
  c1->SaveAs(Form("plots/prettyInterpRatio_UE%d.png",UEtype));
  c1->SaveAs(Form("plots/prettyInterpRatio_UE%d.pdf",UEtype));
  c1->SaveAs(Form("plots/prettyInterpRatio_UE%d.C",UEtype));
  
  for(int i = 10; i<15; i++)
  {
    double yMaximum = 0;
    double yMinimum = 0.60000001;
    if(UEtype == 0) yMaximum = 1.39999;
    if(UEtype == 3) yMaximum = 1.59999;
    hist[i]->GetYaxis()->SetRangeUser(yMinimum,yMaximum);
 
  hist[i]->Draw();
  for(int j=1; j<hist[i]->GetSize()-1; j++)
  { 
    if(turnOffPoints && i==13 && j==21)
    {
      hist[i]->SetBinContent(j,0);
      hist[i]->SetBinError(j,0);
      continue;
    } 
    if(turnOffPoints && i==8 && j==21)
    {
      hist[i]->SetBinContent(j,0);
      hist[i]->SetBinError(j,0);
      continue;
    }
    if(boxArray[i][j]->GetY2()>yMaximum) boxArray[i][j]->SetY2(yMaximum);
    if(boxArray[i][j]->GetY1()<yMinimum) boxArray[i][j]->SetY1(yMinimum);
    if(boxArray[i][j]->GetY2()!=0) boxArray[i][j]->Draw("l");
  }
  hist[i]->Draw("same");
    //tlat->DrawLatex(0.1,0.7,Form("%d < p_{T}^{jet} < %d GeV/c",(int)FF_Bound[i-10],(int)FF_Bound[i+1-10]));
    if(i>9)
    {
      TLine * l = new TLine(0.5,1, hist[i]->GetBinLowEdge(40),1);
      l->SetLineWidth(1);
      l->SetLineStyle(3);
      l->SetLineColor(1);
      l->Draw("same");
    }
  }
  for(int i=10; i<15; i++)
  {
    c1->cd(i-9);
    tlat->DrawLatex(0.6,0.65,Form("%d < p_{T}^{jet} < %d GeV/c",(int)FF_Bound[i-10],(int)FF_Bound[i+1-10]));
  }

  c1->Update();
  c1->SaveAs(Form("plots/prettyInterpRatioZoom_UE%d.png",UEtype));
  c1->SaveAs(Form("plots/prettyInterpRatioZoom_UE%d.pdf",UEtype));

}
