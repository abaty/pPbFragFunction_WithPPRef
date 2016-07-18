#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "TPad.h"
#include "TText.h"
#include "TAttFill.h"
#include "TLine.h"
#include "TColor.h"
#include "TMathText.h"
#include "TMath.h"
#include "TStyle.h"
#include "loadHistograms.h"
#include "makePlots.h"
//#include "interpolationErrors.h"
//#include "systematics.C"
//#include "systematicsSummaries.C"
#include <iostream>
#include <string>

//TODO MakePlots, systematics, systematicssummairies 

//forward declarations
TH1D* getFF_pp(double jetPt_low, double jetPt_high, const char* histTitle, int mode = 0, bool isXi = false);
void getSpectra(int mode);
void Interpolate();

TH1D * jet;
TH2D * trk;
TH2D * trk_xi;
TH2D * trkUE;
TH2D * trkUE_xi;
TH1D * jet_pPb;

//trk pt bins

const int trkBins=22;
const double yAxis[trkBins+1] = {0.5, 0.63, 0.77,  1.03,1.38, 1.84, 2.46, 3.29,  4.40, 5.88,  7.87,  10.52, 14.06,  18.8, 25.13,  33.58,  44.89,  60, 80, 100, 120, 140, 200};

//main execution starts here
void makeFF(int v, int UEtype=3)
{ 
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  //initializing histograms for analysis
  loadHistos(v,UEtype);

  for(int i = 0; i < 2*FF_Bins; i++)
  {
    //pPb direction
    pPb5TeV_data[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_data_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),2,i>=FF_Bins);

    //Pbpdirection
    Pbp5TeV_data[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_data_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),3,i>=FF_Bins);

    ppref5TeV_data[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("ppref5TeV_data_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),4,i>=FF_Bins);
    
    //Reco MC
    pPb5TeV_recoMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_recoMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),6,i>=FF_Bins);
    Pbp5TeV_recoMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_recoMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),51,i>=FF_Bins);
    pPb5Pbp5TeV_recoMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_recoMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),59,i>=FF_Bins);
    pp5TeV_recoMC[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_recoMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),11,i>=FF_Bins);
    
    //Gen MC
    pPb5TeV_genMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_genMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),9,i>=FF_Bins);    
    Pbp5TeV_genMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_genMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),52,i>=FF_Bins);    
    pPb5Pbp5TeV_genMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_genMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),60,i>=FF_Bins);    
    pp5TeV_genMC[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_genMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),10,i>=FF_Bins); 

    //Reco Gen Combinations
    pPb5TeV_rJgTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_rJgTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),35,i>=FF_Bins);
    Pbp5TeV_rJgTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_rJgTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),53,i>=FF_Bins);
    pPb5Pbp5TeV_rJgTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_rJgTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),61,i>=FF_Bins);
    pp5TeV_rJgTMC[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_rJgTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),36,i>=FF_Bins);
    pPb5TeV_gJrTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_gJrTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),39,i>=FF_Bins);
    Pbp5TeV_gJrTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_gJrTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),54,i>=FF_Bins);
    pPb5Pbp5TeV_gJrTMC[i]= getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_gJrTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),62,i>=FF_Bins);
    pp5TeV_gJrTMC[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_gJrTMC_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),40,i>=FF_Bins); 



    //full data set
    pPb5Pbp5TeV_fulldata[i]=getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_fulldata_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),28,i>=FF_Bins);

    //2 and 7 tev pp FF's with no reweighting

    //interpolations
  }
 
  for(int i = 0; i < 2*FF_Bins; i++)
  {
    std::string isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    pPb_FF[i] = (TH1D*) pPb5TeV_data[i]->Clone(Form("pPb_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF[i]->Divide(ppref5TeV_data[i]);
    Pbp_FF[i] = (TH1D*) Pbp5TeV_data[i]->Clone(Form("Pbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF[i]->Divide(ppref5TeV_data[i]);
    pPbPbp_FF[i] = (TH1D*) pPb5Pbp5TeV_fulldata[i]->Clone(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF[i]->Divide(ppref5TeV_data[i]); 

    pPb_FF_recoMC[i] = (TH1D*) pPb5TeV_recoMC[i]->Clone(Form("pPb_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_recoMC[i]->Divide(pp5TeV_recoMC[i]);
    pPb_FF_genMC[i] = (TH1D*) pPb5TeV_genMC[i]->Clone(Form("pPb_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_genMC[i]->Divide(pp5TeV_genMC[i]);
    pPb_FF_rJgTMC[i] = (TH1D*) pPb5TeV_rJgTMC[i]->Clone(Form("pPb_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_rJgTMC[i]->Divide(pp5TeV_rJgTMC[i]);
    pPb_FF_gJrTMC[i] = (TH1D*) pPb5TeV_gJrTMC[i]->Clone(Form("pPb_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_gJrTMC[i]->Divide(pp5TeV_gJrTMC[i]);

    Pbp_FF_recoMC[i] = (TH1D*) Pbp5TeV_recoMC[i]->Clone(Form("Pbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_recoMC[i]->Divide(pp5TeV_recoMC[i]);
    Pbp_FF_genMC[i] = (TH1D*) Pbp5TeV_genMC[i]->Clone(Form("Pbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_genMC[i]->Divide(pp5TeV_genMC[i]);
    Pbp_FF_rJgTMC[i] = (TH1D*) Pbp5TeV_rJgTMC[i]->Clone(Form("Pbp_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_rJgTMC[i]->Divide(pp5TeV_rJgTMC[i]);
    Pbp_FF_gJrTMC[i] = (TH1D*) Pbp5TeV_gJrTMC[i]->Clone(Form("Pbp_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_gJrTMC[i]->Divide(pp5TeV_gJrTMC[i]);

    pPbPbp_FF_recoMC[i] = (TH1D*) pPb5Pbp5TeV_recoMC[i]->Clone(Form("pPbPbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_recoMC[i]->Divide(pp5TeV_recoMC[i]);
    pPbPbp_FF_genMC[i] = (TH1D*) pPb5Pbp5TeV_genMC[i]->Clone(Form("pPbPbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_genMC[i]->Divide(pp5TeV_genMC[i]);
    pPbPbp_FF_rJgTMC[i] = (TH1D*) pPb5Pbp5TeV_rJgTMC[i]->Clone(Form("pPbPbp_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_rJgTMC[i]->Divide(pp5TeV_rJgTMC[i]);
    pPbPbp_FF_gJrTMC[i] = (TH1D*) pPb5Pbp5TeV_gJrTMC[i]->Clone(Form("pPbPbp_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_gJrTMC[i]->Divide(pp5TeV_gJrTMC[i]);

    pp5_FF_recoMC[i] = (TH1D*) pp5TeV_recoMC[i]->Clone(Form("pp5_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pp5_FF_recoMC[i]->Divide(pp5TeV_recoMC[i]);
    pp5_FF_genMC[i] = (TH1D*) pp5TeV_genMC[i]->Clone(Form("pp5_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pp5_FF_genMC[i]->Divide(pp5TeV_genMC[i]);
  }
 
  TFile * outfile = new TFile(Form("FragmentationFunctions%sUE%d.root",variationTag[v],UEtype),"recreate");
  for(int i = 0; i < 2*FF_Bins; i++)
  {
    pPb5TeV_data[i]->Write();
    Pbp5TeV_data[i]->Write();
    ppref5TeV_data[i]->Write();
    pPb5Pbp5TeV_fulldata[i]->Write();

    pPb5TeV_recoMC[i]->Write();
    Pbp5TeV_recoMC[i]->Write();
    pPb5Pbp5TeV_recoMC[i]->Write();
    pp5TeV_recoMC[i]->Write();

    pPb5TeV_rJgTMC[i]->Write();
    Pbp5TeV_rJgTMC[i]->Write();
    pPb5Pbp5TeV_rJgTMC[i]->Write();
    pp5TeV_rJgTMC[i]->Write();
    pPb5TeV_gJrTMC[i]->Write();
    Pbp5TeV_gJrTMC[i]->Write();
    pPb5Pbp5TeV_gJrTMC[i]->Write();
    pp5TeV_gJrTMC[i]->Write();
 

    pPb5TeV_genMC[i]->Write();
    Pbp5TeV_genMC[i]->Write();
    pPb5Pbp5TeV_genMC[i]->Write();
    pp5TeV_genMC[i]->Write(); 


    pPb_FF[i]->Write();
    Pbp_FF[i]->Write();  
    pPbPbp_FF[i]->Write(); 
    pPb_FF_recoMC[i]->Write();
    pPb_FF_genMC[i]->Write();
    pPb_FF_rJgTMC[i]->Write();
    pPb_FF_gJrTMC[i]->Write();

    Pbp_FF_recoMC[i]->Write();
    Pbp_FF_genMC[i]->Write();
    Pbp_FF_rJgTMC[i]->Write();
    Pbp_FF_gJrTMC[i]->Write();

    pPbPbp_FF_recoMC[i]->Write();
    pPbPbp_FF_genMC[i]->Write();
    pPbPbp_FF_rJgTMC[i]->Write();
    pPbPbp_FF_gJrTMC[i]->Write();
  }
  //handing it over to a plotting macro
  makePlots(variationTag[v],UEtype); 
  outfile->Close(); 
}

TH1D* getFF_pp(double jetPt_low, double jetPt_high, const char* histTitle, int mode, bool isXi)
{
  TH1D * FF;
  getSpectra(mode);
  
  if(!isXi)
  {
    FF = new TH1D(histTitle,";p_{T trk};#frac{1}{N_{jet}} #frac{dN_{trk}}{dp_{t trk}}",trkBins ,yAxis);

    //looping over bins in Trk pt
    for(int t = 1; t < trk->GetYaxis()->FindBin(jetPt_high); t++)
    {
      double sum   = 0;
      double error = 0;
      double norm  = 0;        

      //looping over jet pt
      for(int j = jet->FindBin(jetPt_low); j < jet->FindBin(jetPt_high); j++)
      {
        //reweighting factor to 5 TeV pPb
        //for (mode 2) 5TeV pPb jet_pPb = jet so this is by definition one
        double w = jet_pPb->GetBinContent(j)/jet->GetBinContent(j);

        //adding up contributions to the FF from each track and jet bin
        sum   += w*(trk->GetBinContent(j,t) - trkUE->GetBinContent(j,t));
        error += w*w*(TMath::Power(trk->GetBinError(j,t),2) + TMath::Power(trkUE->GetBinError(j,t),2));
        norm  += jet_pPb->GetBinContent(j);
      }

      if(norm!=0)
      {   
        sum   = sum/(norm*FF->GetBinWidth(t));
        error = TMath::Power(error,0.5)/(norm*FF->GetBinWidth(t));
        FF->SetBinContent(t,sum);
        FF->SetBinError(t,error);
      }
    }
  }
  else 
  {
    FF = new TH1D(Form("%s_xi",histTitle),";#xi;#frac{1}{N_{jet}} #frac{dN_{trk}}{d#xi}",24,0.0,6.0);

    //looping over bins in Trk pt
    for(int t = 1; t < 25; t++)
    {
      double sum   = 0;
      double error = 0;
      double norm  = 0;        

      //looping over jet pt
      for(int j = jet->FindBin(jetPt_low); j < jet->FindBin(jetPt_high); j++)
      {
        //reweighting factor to 5 TeV pPb
        //for (mode 2) 5TeV pPb jet_pPb = jet so this is by definition one
        double w = jet_pPb->GetBinContent(j)/jet->GetBinContent(j);

        //adding up contributions to the FF from each track and jet bin
        sum   += w*(trk_xi->GetBinContent(j,t) - trkUE_xi->GetBinContent(j,t));
        error += w*w*(TMath::Power(trk_xi->GetBinError(j,t),2) + TMath::Power(trkUE_xi->GetBinError(j,t),2));
        norm  += jet_pPb->GetBinContent(j);
      }

      if(norm!=0)
      {   
        sum   = sum/(norm*FF->GetBinWidth(t));
        error = TMath::Power(error,0.5)/(norm*FF->GetBinWidth(t));
        FF->SetBinContent(t,sum);
        FF->SetBinError(t,error);
      }
    } 
  }
  return FF;
}


void Interpolate()
{
  for(int v = 0; v<variations; v++) 
  {
    if((v<1 || v>4) && v!=7 && v!=8 && v!=10 && v!=11 && (v<14 || v>23) && v!=29 && (v<31 || v>33))
    {
      makeFF(v,0);
      if(v!=26) makeFF(v,3);
      if(v==0)  makeFF(v,2);
    }
  }

  //SYSTEMATICS HERE
  //systematics(0);
  //systematics(3);
  //systematicSummaries();
}

//tells the getFF_pp exactly which spectra to use in which jetPt range
//Numbers here in the if statements are hard-coded for triggers atm
//to do: combine the spectra before so I can remove this step
void getSpectra(int mode)
{

  if(mode == 2)
  {
    jet     = pPb5_0_jet;
    jet_pPb = pPb5_0_jet;
    trk     = pPb5_0_track;
    trkUE   = pPb5_0_trackUE;
    trk_xi     = pPb5_0_track_xi;
    trkUE_xi   = pPb5_0_trackUE_xi;
  }
  if(mode == 3)
  {
    jet     = Pbp5_0_jet;
    jet_pPb = Pbp5_0_jet;
    trk     = Pbp5_0_track;
    trkUE   = Pbp5_0_trackUE;
    trk_xi     = Pbp5_0_track_xi;
    trkUE_xi   = Pbp5_0_trackUE_xi;
  }
  if(mode == 4)
  {
    jet     = ppref5_0_jet;
    jet_pPb = ppref5_0_jet;
    trk     = ppref5_0_track;
    trkUE   = ppref5_0_trackUE;
    trk_xi     = ppref5_0_track_xi;
    trkUE_xi   = ppref5_0_trackUE_xi;
  }

  if(mode == 6)
  {
    jet     = pPb5_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pPb5_1_track_reco;
    trkUE   = pPb5_1_trackUE_reco;
    trk_xi     = pPb5_1_track_xi_reco;
    trkUE_xi   = pPb5_1_trackUE_xi_reco;
  }

  if(mode == 9)
  {
    jet     = pPb5_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pPb5_1_track_gen;
    trkUE   = pPb5_1_trackUE_gen;
    trk_xi     = pPb5_1_track_xi_gen;
    trkUE_xi   = pPb5_1_trackUE_xi_gen;
  }
  
  if(mode == 10)
  {
    jet     = pp5_1_jet_gen;
    jet_pPb = pp5_1_jet_gen;
    trk     = pp5_1_track_gen;
    trkUE   = pp5_1_trackUE_gen;
    trk_xi     = pp5_1_track_xi_gen;
    trkUE_xi   = pp5_1_trackUE_xi_gen;
  }

  if(mode == 11)
  {
    jet     = pp5_1_jet_reco;
    jet_pPb = pp5_1_jet_reco;
    trk     = pp5_1_track_reco;
    trkUE   = pp5_1_trackUE_reco;
    trk_xi     = pp5_1_track_xi_reco;
    trkUE_xi   = pp5_1_trackUE_xi_reco;
  }
  if(mode==28)
  {
    jet     = pPb5Pbp5_0_jet;
    jet_pPb = pPb5Pbp5_0_jet;
    trk     = pPb5Pbp5_0_track;
    trkUE   = pPb5Pbp5_0_trackUE;
    trk_xi     = pPb5Pbp5_0_track_xi;
    trkUE_xi   = pPb5Pbp5_0_trackUE_xi;
  }
  if(mode == 35)
  {
    jet     = pPb5_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pPb5_1_track_rJgT;
    trkUE   = pPb5_1_trackUE_rJgT;
    trk_xi     = pPb5_1_track_xi_rJgT;
    trkUE_xi   = pPb5_1_trackUE_xi_rJgT;
  } 
  if(mode == 36)
  {
    jet     = pp5_1_jet_reco;
    jet_pPb = pp5_1_jet_reco;
    trk     = pp5_1_track_rJgT;
    trkUE   = pp5_1_trackUE_rJgT;
    trk_xi     = pp5_1_track_xi_rJgT;
    trkUE_xi   = pp5_1_trackUE_xi_rJgT;
  }
  if(mode == 39)
  {
    jet     = pPb5_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pPb5_1_track_gJrT;
    trkUE   = pPb5_1_trackUE_gJrT;
    trk_xi     = pPb5_1_track_xi_gJrT;
    trkUE_xi   = pPb5_1_trackUE_xi_gJrT;
  } 
  if(mode == 40)
  {
    jet     = pp5_1_jet_gen;
    jet_pPb = pp5_1_jet_gen;
    trk     = pp5_1_track_gJrT;
    trkUE   = pp5_1_trackUE_gJrT;
    trk_xi     = pp5_1_track_xi_gJrT;
    trkUE_xi   = pp5_1_trackUE_xi_gJrT;
  } 
  if(mode == 51)
  {
    jet     = Pbp5_1_jet_reco;
    jet_pPb = Pbp5_1_jet_reco;
    trk     = Pbp5_1_track_reco;
    trkUE   = Pbp5_1_trackUE_reco;
    trk_xi     = Pbp5_1_track_xi_reco;
    trkUE_xi   = Pbp5_1_trackUE_xi_reco;
  }
  if(mode == 52)
  {
    jet     = Pbp5_1_jet_gen;
    jet_pPb = Pbp5_1_jet_gen;
    trk     = Pbp5_1_track_gen;
    trkUE   = Pbp5_1_trackUE_gen;
    trk_xi     = Pbp5_1_track_xi_gen;
    trkUE_xi   = Pbp5_1_trackUE_xi_gen;
  }
  if(mode == 53)
  {
    jet     = Pbp5_1_jet_reco;
    jet_pPb = Pbp5_1_jet_reco;
    trk     = Pbp5_1_track_rJgT;
    trkUE   = Pbp5_1_trackUE_rJgT;
    trk_xi     = Pbp5_1_track_xi_rJgT;
    trkUE_xi   = Pbp5_1_trackUE_xi_rJgT;
  } 
  if(mode == 54)
  {
    jet     = Pbp5_1_jet_gen;
    jet_pPb = Pbp5_1_jet_gen;
    trk     = Pbp5_1_track_gJrT;
    trkUE   = Pbp5_1_trackUE_gJrT;
    trk_xi     = Pbp5_1_track_xi_gJrT;
    trkUE_xi   = Pbp5_1_trackUE_xi_gJrT;
  } 
  if(mode == 59)
  {
    jet     = pPb5Pbp5_1_jet_reco;
    jet_pPb = pPb5Pbp5_1_jet_reco;
    trk     = pPb5Pbp5_1_track_reco;
    trkUE   = pPb5Pbp5_1_trackUE_reco;
    trk_xi     = pPb5Pbp5_1_track_xi_reco;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_reco;
  }
  if(mode == 60)
  {
    jet     = pPb5Pbp5_1_jet_gen;
    jet_pPb = pPb5Pbp5_1_jet_gen;
    trk     = pPb5Pbp5_1_track_gen;
    trkUE   = pPb5Pbp5_1_trackUE_gen;
    trk_xi     = pPb5Pbp5_1_track_xi_gen;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_gen;
  }
  if(mode == 61)
  {
    jet     = pPb5Pbp5_1_jet_reco;
    jet_pPb = pPb5Pbp5_1_jet_reco;
    trk     = pPb5Pbp5_1_track_rJgT;
    trkUE   = pPb5Pbp5_1_trackUE_rJgT;
    trk_xi     = pPb5Pbp5_1_track_xi_rJgT;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_rJgT;
  } 
  if(mode == 62)
  {
    jet     = pPb5Pbp5_1_jet_gen;
    jet_pPb = pPb5Pbp5_1_jet_gen;
    trk     = pPb5Pbp5_1_track_gJrT;
    trkUE   = pPb5Pbp5_1_trackUE_gJrT;
    trk_xi     = pPb5Pbp5_1_track_xi_gJrT;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_gJrT;
  } 
  return;
}

