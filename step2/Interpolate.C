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
#include "interpolationErrors.h"
#include "systematics.C"
#include "systematicsSummaries.C"
#include <iostream>
#include <string>

//forward declarations
TH1D* getFF_pp(double jetPt_low, double jetPt_high, const char* histTitle, int mode = 0, bool isXi = false);
TH1D** getInterpolation(double jetPt_low, double jetPt_high, const char* histTitle, int ReweightMode, TH1D * hist2, TH1D * hist7, int doGenGluFrac=0, bool isXi = false);
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

    //reco MC Q/G study
    pPb5TeV_recoMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_recoMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),16,i>=FF_Bins);
    pPb5TeV_recoMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_recoMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),17,i>=FF_Bins);
    Pbp5TeV_recoMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_recoMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),55,i>=FF_Bins);
    Pbp5TeV_recoMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_recoMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),56,i>=FF_Bins);
    pPb5Pbp5TeV_recoMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_recoMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),63,i>=FF_Bins);
    pPb5Pbp5TeV_recoMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_recoMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),64,i>=FF_Bins);
    pp5TeV_recoMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_recoMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),18,i>=FF_Bins);
    pp5TeV_recoMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_recoMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),19,i>=FF_Bins);

    //gen MC Q/G study
    pPb5TeV_genMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_genMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),24,i>=FF_Bins);
    pPb5TeV_genMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5TeV_genMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),25,i>=FF_Bins);
    Pbp5TeV_genMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_genMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),57,i>=FF_Bins);
    Pbp5TeV_genMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("Pbp5TeV_genMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),58,i>=FF_Bins);
    pPb5Pbp5TeV_genMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_genMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),65,i>=FF_Bins);
    pPb5Pbp5TeV_genMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pPb5Pbp5TeV_genMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),66,i>=FF_Bins);
    pp5TeV_genMC_Q[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_genMC_Q_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),26,i>=FF_Bins);
    pp5TeV_genMC_G[i] = getFF_pp(FF_Bound[i%FF_Bins],FF_Bound[i%FF_Bins+1],Form("pp5TeV_genMC_G_%d_%d",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1]),27,i>=FF_Bins);

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
    pPb_FF[i]->Divide(pPb5TeV_data_interp[i][0]);
    Pbp_FF[i] = (TH1D*) Pbp5TeV_data[i]->Clone(Form("Pbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF[i]->Divide(Pbp5TeV_data_interp[i][0]);
    pPbPbp_FF[i] = (TH1D*) pPb5Pbp5TeV_fulldata[i]->Clone(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF[i]->Divide(pPb5Pb5TeV_data_interp[i][0]); 
    pPbPbp_FF_genGluFrac[i] = (TH1D*) pPb5Pbp5TeV_fulldata[i]->Clone(Form("pPbPbp_FF_genGluFrac_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_genGluFrac[i]->Divide(pPb5Pb5TeV_data_interp_genGluFrac[i][0]); 

    pPb_FF_recoMC[i] = (TH1D*) pPb5TeV_recoMC[i]->Clone(Form("pPb_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_recoMC[i]->Divide(pPb5TeV_recoMC_interp[i][0]);
    pPb_FF_genMC[i] = (TH1D*) pPb5TeV_genMC[i]->Clone(Form("pPb_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_genMC[i]->Divide(pPb5TeV_genMC_interp[i][0]);
    pPb_FF_rJgTMC[i] = (TH1D*) pPb5TeV_rJgTMC[i]->Clone(Form("pPb_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_rJgTMC[i]->Divide(pPb5TeV_rJgTMC_interp[i][0]);
    pPb_FF_gJrTMC[i] = (TH1D*) pPb5TeV_gJrTMC[i]->Clone(Form("pPb_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPb_FF_gJrTMC[i]->Divide(pPb5TeV_gJrTMC_interp[i][0]);

    Pbp_FF_recoMC[i] = (TH1D*) Pbp5TeV_recoMC[i]->Clone(Form("Pbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_recoMC[i]->Divide(Pbp5TeV_recoMC_interp[i][0]);
    Pbp_FF_genMC[i] = (TH1D*) Pbp5TeV_genMC[i]->Clone(Form("Pbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_genMC[i]->Divide(Pbp5TeV_genMC_interp[i][0]);
    Pbp_FF_rJgTMC[i] = (TH1D*) Pbp5TeV_rJgTMC[i]->Clone(Form("Pbp_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_rJgTMC[i]->Divide(Pbp5TeV_rJgTMC_interp[i][0]);
    Pbp_FF_gJrTMC[i] = (TH1D*) Pbp5TeV_gJrTMC[i]->Clone(Form("Pbp_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    Pbp_FF_gJrTMC[i]->Divide(Pbp5TeV_gJrTMC_interp[i][0]);

    pPbPbp_FF_recoMC[i] = (TH1D*) pPb5Pbp5TeV_recoMC[i]->Clone(Form("pPbPbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_recoMC[i]->Divide(pPb5Pbp5TeV_recoMC_interp[i][0]);
    pPbPbp_FF_genMC[i] = (TH1D*) pPb5Pbp5TeV_genMC[i]->Clone(Form("pPbPbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_genMC[i]->Divide(pPb5Pbp5TeV_genMC_interp[i][0]);
    pPbPbp_FF_rJgTMC[i] = (TH1D*) pPb5Pbp5TeV_rJgTMC[i]->Clone(Form("pPbPbp_FF_rJgTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_rJgTMC[i]->Divide(pPb5Pbp5TeV_rJgTMC_interp[i][0]);
    pPbPbp_FF_gJrTMC[i] = (TH1D*) pPb5Pbp5TeV_gJrTMC[i]->Clone(Form("pPbPbp_FF_gJrTMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbPbp_FF_gJrTMC[i]->Divide(pPb5Pbp5TeV_gJrTMC_interp[i][0]);

    pp5_FF_recoMC[i] = (TH1D*) pp5TeV_recoMC[i]->Clone(Form("pp5_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pp5_FF_recoMC[i]->Divide(pPb5Pbp5TeV_recoMC_interp[i][0]);
    pp5_FF_genMC[i] = (TH1D*) pp5TeV_genMC[i]->Clone(Form("pp5_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pp5_FF_genMC[i]->Divide(pPb5Pbp5TeV_genMC_interp[i][0]);
  }
 
  TFile * outfile = new TFile(Form("FragmentationFunctions%sUE%d.root",variationTag[v],UEtype),"recreate");
  for(int i = 0; i < 2*FF_Bins; i++)
  {
    pPb5TeV_data[i]->Write();
    Pbp5TeV_data[i]->Write();
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
 
    pPb5TeV_recoMC_Q[i]->Write();
    pPb5TeV_recoMC_G[i]->Write();
    Pbp5TeV_recoMC_Q[i]->Write();
    Pbp5TeV_recoMC_G[i]->Write();
    pPb5Pbp5TeV_recoMC_Q[i]->Write();
    pPb5Pbp5TeV_recoMC_G[i]->Write();
    pp5TeV_recoMC_Q[i]->Write();
    pp5TeV_recoMC_G[i]->Write();

    pPb5TeV_genMC[i]->Write();
    Pbp5TeV_genMC[i]->Write();
    pPb5Pbp5TeV_genMC[i]->Write();
    pp5TeV_genMC[i]->Write(); 
    pPb5TeV_genMC_Q[i]->Write();
    pPb5TeV_genMC_G[i]->Write();
    Pbp5TeV_genMC_Q[i]->Write();
    Pbp5TeV_genMC_G[i]->Write();
    pPb5Pbp5TeV_genMC_Q[i]->Write();
    pPb5Pbp5TeV_genMC_G[i]->Write();
    pp5TeV_genMC_Q[i]->Write();
    pp5TeV_genMC_G[i]->Write();


    for(int indx = 0; indx<3; indx++)
    {
      pPb5TeV_data_interp[i][indx]->Write();
      Pbp5TeV_data_interp[i][indx]->Write();
      pPb5Pb5TeV_data_interp[i][indx]->Write();
      pPb5Pb5TeV_data_interp_genGluFrac[i][indx]->Write();
      pPb5TeV_recoMC_interp[i][indx]->Write();
      pPb5TeV_genMC_interp[i][indx]->Write();
      pPb5TeV_rJgTMC_interp[i][indx]->Write();
      pPb5TeV_gJrTMC_interp[i][indx]->Write();
      Pbp5TeV_recoMC_interp[i][indx]->Write();
      Pbp5TeV_genMC_interp[i][indx]->Write();
      Pbp5TeV_rJgTMC_interp[i][indx]->Write();
      Pbp5TeV_gJrTMC_interp[i][indx]->Write();
      pPb5Pbp5TeV_recoMC_interp[i][indx]->Write();
      pPb5Pbp5TeV_genMC_interp[i][indx]->Write();
      pPb5Pbp5TeV_rJgTMC_interp[i][indx]->Write();
      pPb5Pbp5TeV_gJrTMC_interp[i][indx]->Write();
    } 

    pPb_FF[i]->Write();
    Pbp_FF[i]->Write();  
    pPbPbp_FF[i]->Write(); 
    pPbPbp_FF_genGluFrac[i]->Write(); 
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

TH1D** getInterpolation(double jetPt_low, double jetPt_high, const char* histTitle, int ReweightMode, TH1D * hist2, TH1D * hist7, int doGenGluFrac, bool isXi)
{
  TH1D * tmp_interp;
  TH1D * tmp_Q_interp;
  TH1D * tmp_G_interp;
  
  //setting jet spectrum to 5TeV for use in the interpolation weighting (only "jet" is used)
  getSpectra(ReweightMode);

  int binLimit;
  if(!isXi)
  {
    tmp_interp = new TH1D(Form("%s_%d_%d",histTitle,(int)jetPt_low,(int)jetPt_high),";p_{T trk};#frac{1}{N_{jet}} #frac{dN_{trk}}{dp_{t trk}}",trkBins,yAxis);
    tmp_Q_interp = new TH1D(Form("%s_Q_%d_%d",histTitle,(int)jetPt_low,(int)jetPt_high),";p_{T trk};#frac{1}{N_{jet}} #frac{dN_{trk}}{dp_{t trk}}",trkBins,yAxis);
    tmp_G_interp = new TH1D(Form("%s_G_%d_%d",histTitle,(int)jetPt_low,(int)jetPt_high),";p_{T trk};#frac{1}{N_{jet}} #frac{dN_{trk}}{dp_{t trk}}",trkBins,yAxis);    
    binLimit = trkBins+1;
  }
  else
  {
    tmp_interp = new TH1D(Form("%s_%d_%d_xi",histTitle,(int)jetPt_low,(int)jetPt_high),";#xi;#frac{1}{N_{jet}} #frac{dN_{trk}}{d#xi}",24,0,6);
    tmp_Q_interp = new TH1D(Form("%s_Q_%d_%d_xi",histTitle,(int)jetPt_low,(int)jetPt_high),";#xi;#frac{1}{N_{jet}} #frac{dN_{trk}}{d#xi}",24,0,6);
    tmp_G_interp = new TH1D(Form("%s_G_%d_%d_xi",histTitle,(int)jetPt_low,(int)jetPt_high),";#xi;#frac{1}{N_{jet}} #frac{dN_{trk}}{d#xi}",24,0,6);
    binLimit = 25;
  }
  for(int t = 1; t < binLimit; t++)
  {
    //number of gluon jets for 2/5/7 MC (proportional to the jet fraction because jet spectra are all reweighted to 5TeV)
    double glu2 = 0;
    double glu5 = 0;
    double glu7 = 0;
    double totalJets = 0;    
    double glu2Err = 0;
    double glu5Err = 0;
    double glu7Err = 0;

    for(int j = jet->FindBin(jetPt_low); j < jet->FindBin(jetPt_high); j++)
    {
      double nJet = jet->GetBinContent(j);
      double nJetErr = jet->GetBinError(j); 

      if(doGenGluFrac==0)
      {
        glu2 += gluon_2tev_reco->GetBinContent(j)*nJet;
        glu5 += gluon_5tev_reco->GetBinContent(j)*nJet;
        glu7 += gluon_7tev_reco->GetBinContent(j)*nJet;
        glu2Err += TMath::Power(gluon_2tev_reco->GetBinError(j)*nJet,2)+ TMath::Power(gluon_2tev_reco->GetBinContent(j)*nJetErr,2);
        glu5Err += TMath::Power(gluon_5tev_reco->GetBinError(j)*nJet,2)+ TMath::Power(gluon_5tev_reco->GetBinContent(j)*nJetErr,2);
        glu7Err += TMath::Power(gluon_7tev_reco->GetBinError(j)*nJet,2)+ TMath::Power(gluon_7tev_reco->GetBinContent(j)*nJetErr,2);
      }
      else
      {
        glu2 += gluon_2tev_gen->GetBinContent(j)*nJet;
        glu5 += gluon_5tev_gen->GetBinContent(j)*nJet;
        glu7 += gluon_7tev_gen->GetBinContent(j)*nJet;
        glu2Err += TMath::Power(gluon_2tev_gen->GetBinError(j)*nJet,2)+ TMath::Power(gluon_2tev_gen->GetBinContent(j)*nJetErr,2);
        glu5Err += TMath::Power(gluon_5tev_gen->GetBinError(j)*nJet,2)+ TMath::Power(gluon_5tev_gen->GetBinContent(j)*nJetErr,2);
        glu7Err += TMath::Power(gluon_7tev_gen->GetBinError(j)*nJet,2)+ TMath::Power(gluon_7tev_gen->GetBinContent(j)*nJetErr,2);
      }
      totalJets += nJet;
    }
    glu2Err = TMath::Power(glu2Err,0.5);
    glu5Err = TMath::Power(glu5Err,0.5);
    glu7Err = TMath::Power(glu7Err,0.5);   
 
    double average = ((glu5 - glu7)*hist2->GetBinContent(t)+(glu2-glu5)*hist7->GetBinContent(t))/(glu2-glu7);
    double error = getInterpolationError(glu2,glu2Err,glu5,glu5Err,glu7,glu7Err,hist2->GetBinContent(t),hist2->GetBinError(t),hist7->GetBinContent(t),hist7->GetBinError(t));
    //double error = TMath::Power(TMath::Power((glu5 - glu7)*pp2TeV_data[i]->GetBinError(t)/(glu2-glu7),2)+TMath::Power((glu2-glu5)*pp7TeV_data[i]->GetBinError(t)/(glu2-glu7),2),0.5);
    double averageQ = (glu2*hist7->GetBinContent(t)-glu7*hist2->GetBinContent(t))/(glu2-glu7);
    double errorQ = TMath::Power(TMath::Power(glu2*hist7->GetBinError(t)/(glu2-glu7),2)+TMath::Power(glu7*hist2->GetBinError(t)/(glu2-glu7),2),0.5);
    double averageG = ((totalJets-glu7)*hist2->GetBinContent(t)-(totalJets-glu2)*hist7->GetBinContent(t))/(glu2-glu7);
    double errorG = TMath::Power(TMath::Power((totalJets-glu7)*hist2->GetBinError(t)/(glu2-glu7),2)+TMath::Power((totalJets-glu2)*hist7->GetBinError(t)/(glu2-glu7),2),0.5);     
 
    tmp_interp->SetBinContent(t, average);
    tmp_interp->SetBinError(t, error);
    tmp_Q_interp->SetBinContent(t, averageQ);
    tmp_Q_interp->SetBinError(t, errorQ);
    tmp_G_interp->SetBinContent(t, averageG);
    tmp_G_interp->SetBinError(t, errorG);
  }

  TH1D** outputArray = new TH1D*[3];
  outputArray[0]=tmp_interp;
  outputArray[1]=tmp_Q_interp;
  outputArray[2]=tmp_G_interp;
  return outputArray;
}

void Interpolate()
{
  for(int v = 0; v<variations; v++) 
  {
    if(v<14 || v>19)
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
  if(mode == 0)
  {
    jet     = pp2_0_jet;
    jet_pPb = pPb5_0_jet;
    trk     = pp2_0_track;
    trkUE   = pp2_0_trackUE;
    trk_xi     = pp2_0_track_xi;
    trkUE_xi   = pp2_0_trackUE_xi;
  }
  
  if(mode == 1)
  {
    jet     = pp7_0_jet;
    jet_pPb = pPb5_0_jet;
    trk     = pp7_0_track;
    trkUE   = pp7_0_trackUE;
    trk_xi    = pp7_0_track_xi;
    trkUE_xi  = pp7_0_trackUE_xi;
  }

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
    jet     = pp2_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pp2_1_track_reco;
    trkUE   = pp2_1_trackUE_reco;
    trk_xi     = pp2_1_track_xi_reco;
    trkUE_xi   = pp2_1_trackUE_xi_reco;
  }

  if(mode == 5)
  {
    jet     = pp7_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pp7_1_track_reco;
    trkUE   = pp7_1_trackUE_reco;
    trk_xi     = pp7_1_track_xi_reco;
    trkUE_xi   = pp7_1_trackUE_xi_reco;
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
  if(mode == 7)
  {
    jet     = pp2_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pp2_1_track_gen;
    trkUE   = pp2_1_trackUE_gen;
    trk_xi     = pp2_1_track_xi_gen;
    trkUE_xi   = pp2_1_trackUE_xi_gen;
  }

  if(mode == 8)
  {
    jet     = pp7_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pp7_1_track_gen;
    trkUE   = pp7_1_trackUE_gen;
    trk_xi     = pp7_1_track_xi_gen;
    trkUE_xi   = pp7_1_trackUE_xi_gen;
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
  
  if(mode == 12)
  {
    jet     = pp2_1_jet_reco_Q;
    jet_pPb = pPb5_1_jet_reco_Q;
    trk     = pp2_1_track_reco_Q;
    trkUE   = pp2_1_trackUE_reco_Q;
    trk_xi     = pp2_1_track_xi_reco_Q;
    trkUE_xi   = pp2_1_trackUE_xi_reco_Q;
  }
  if(mode == 13)
  {
    jet     = pp2_1_jet_reco_G;
    jet_pPb = pPb5_1_jet_reco_G;
    trk     = pp2_1_track_reco_G;
    trkUE   = pp2_1_trackUE_reco_G;
    trk_xi     = pp2_1_track_xi_reco_G;
    trkUE_xi   = pp2_1_trackUE_xi_reco_G;
  }
  if(mode == 14)
  {
    jet     = pp7_1_jet_reco_Q;
    jet_pPb = pPb5_1_jet_reco_Q;
    trk     = pp7_1_track_reco_Q;
    trkUE   = pp7_1_trackUE_reco_Q;
    trk_xi     = pp7_1_track_xi_reco_Q;
    trkUE_xi   = pp7_1_trackUE_xi_reco_Q;
  }
  if(mode == 15)
  {
    jet     = pp7_1_jet_reco_G;
    jet_pPb = pPb5_1_jet_reco_G;
    trk     = pp7_1_track_reco_G;
    trkUE   = pp7_1_trackUE_reco_G;
    trk_xi     = pp7_1_track_xi_reco_G;
    trkUE_xi   = pp7_1_trackUE_xi_reco_G;
  }
  if(mode == 16)
  {
    jet     = pPb5_1_jet_reco_Q;
    jet_pPb = pPb5_1_jet_reco_Q;
    trk     = pPb5_1_track_reco_Q;
    trkUE   = pPb5_1_trackUE_reco_Q;
    trk_xi     = pPb5_1_track_xi_reco_Q;
    trkUE_xi   = pPb5_1_trackUE_xi_reco_Q;
  }
  if(mode == 17)
  {
    jet     = pPb5_1_jet_reco_G;
    jet_pPb = pPb5_1_jet_reco_G;
    trk     = pPb5_1_track_reco_G;
    trkUE   = pPb5_1_trackUE_reco_G;
    trk_xi     = pPb5_1_track_xi_reco_G;
    trkUE_xi   = pPb5_1_trackUE_xi_reco_G;
  }
  if(mode == 18)
  {
    jet     = pp5_1_jet_reco_Q;
    jet_pPb = pp5_1_jet_reco_Q;
    trk     = pp5_1_track_reco_Q;
    trkUE   = pp5_1_trackUE_reco_Q;
    trk_xi     = pp5_1_track_xi_reco_Q;
    trkUE_xi   = pp5_1_trackUE_xi_reco_Q;
  }
  if(mode == 19)
  {
    jet     = pp5_1_jet_reco_G;
    jet_pPb = pp5_1_jet_reco_G;
    trk     = pp5_1_track_reco_G;
    trkUE   = pp5_1_trackUE_reco_G;
    trk_xi     = pp5_1_track_xi_reco_G;
    trkUE_xi   = pp5_1_trackUE_xi_reco_G;
  }

  if(mode == 20)
  {
    jet     = pp2_1_jet_gen_Q;
    jet_pPb = pPb5_1_jet_gen_Q;
    trk     = pp2_1_track_gen_Q;
    trkUE   = pp2_1_trackUE_gen_Q;
    trk_xi     = pp2_1_track_xi_gen_Q;
    trkUE_xi   = pp2_1_trackUE_xi_gen_Q;
  }
  if(mode == 21)
  {
    jet     = pp2_1_jet_gen_G;
    jet_pPb = pPb5_1_jet_gen_G;
    trk     = pp2_1_track_gen_G;
    trkUE   = pp2_1_trackUE_gen_G;
    trk_xi     = pp2_1_track_xi_gen_G;
    trkUE_xi   = pp2_1_trackUE_xi_gen_G;
  }
  if(mode == 22)
  {
    jet     = pp7_1_jet_gen_Q;
    jet_pPb = pPb5_1_jet_gen_Q;
    trk     = pp7_1_track_gen_Q;
    trkUE   = pp7_1_trackUE_gen_Q;
    trk_xi     = pp7_1_track_xi_gen_Q;
    trkUE_xi   = pp7_1_trackUE_xi_gen_Q;
  }
  if(mode == 23)
  {
    jet     = pp7_1_jet_gen_G;
    jet_pPb = pPb5_1_jet_gen_G;
    trk     = pp7_1_track_gen_G;
    trkUE   = pp7_1_trackUE_gen_G;
    trk_xi     = pp7_1_track_xi_gen_G;
    trkUE_xi   = pp7_1_trackUE_xi_gen_G;
  }
  if(mode == 24)
  {
    jet     = pPb5_1_jet_gen_Q;
    jet_pPb = pPb5_1_jet_gen_Q;
    trk     = pPb5_1_track_gen_Q;
    trkUE   = pPb5_1_trackUE_gen_Q;
    trk_xi     = pPb5_1_track_xi_gen_Q;
    trkUE_xi   = pPb5_1_trackUE_xi_gen_Q;
  }
  if(mode == 25)
  {
    jet     = pPb5_1_jet_gen_G;
    jet_pPb = pPb5_1_jet_gen_G;
    trk     = pPb5_1_track_gen_G;
    trkUE   = pPb5_1_trackUE_gen_G;
    trk_xi     = pPb5_1_track_xi_gen_G;
    trkUE_xi   = pPb5_1_trackUE_xi_gen_G;
  }
  if(mode == 26)
  {
    jet     = pp5_1_jet_gen_Q;
    jet_pPb = pp5_1_jet_gen_Q;
    trk     = pp5_1_track_gen_Q;
    trkUE   = pp5_1_trackUE_gen_Q;
    trk_xi     = pp5_1_track_xi_gen_Q;
    trkUE_xi   = pp5_1_trackUE_xi_gen_Q;
  }
  if(mode == 27)
  {
    jet     = pp5_1_jet_reco_G;
    jet_pPb = pp5_1_jet_reco_G;
    trk     = pp5_1_track_reco_G;
    trkUE   = pp5_1_trackUE_reco_G;
    trk_xi     = pp5_1_track_xi_reco_G;
    trkUE_xi   = pp5_1_trackUE_xi_reco_G;
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
  if(mode==29)
  {
    jet     = pp2_0_jet;
    jet_pPb = pPb5Pbp5_0_jet;
    trk     = pp2_0_track;
    trkUE   = pp2_0_trackUE;
    trk_xi     = pp2_0_track_xi;
    trkUE_xi   = pp2_0_trackUE_xi;
  }
  if(mode==30)
  {
    jet     = pp7_0_jet;
    jet_pPb = pPb5Pbp5_0_jet;
    trk     = pp7_0_track;
    trkUE   = pp7_0_trackUE;
    trk_xi     = pp7_0_track_xi;
    trkUE_xi   = pp7_0_trackUE_xi;
  }
  if(mode == 31)
  {
    jet     = pp2_0_jet;
    jet_pPb = Pbp5_0_jet;
    trk     = pp2_0_track;
    trkUE   = pp2_0_trackUE;
    trk_xi     = pp2_0_track_xi;
    trkUE_xi   = pp2_0_trackUE_xi;
  }

  if(mode == 32)
  {
    jet     = pp7_0_jet;
    jet_pPb = Pbp5_0_jet;
    trk     = pp7_0_track;
    trkUE   = pp7_0_trackUE;
    trk_xi     = pp7_0_track_xi;
    trkUE_xi   = pp7_0_trackUE_xi;
  }
  if(mode == 33)
  {
    jet     = pp2_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pp2_1_track_rJgT;
    trkUE   = pp2_1_trackUE_rJgT;
    trk_xi     = pp2_1_track_xi_rJgT;
    trkUE_xi   = pp2_1_trackUE_xi_rJgT;
  } 
  if(mode == 34)
  {
    jet     = pp7_1_jet_reco;
    jet_pPb = pPb5_1_jet_reco;
    trk     = pp7_1_track_rJgT;
    trkUE   = pp7_1_trackUE_rJgT;
    trk_xi     = pp7_1_track_xi_rJgT;
    trkUE_xi   = pp7_1_trackUE_xi_rJgT;
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
  if(mode == 37)
  {
    jet     = pp2_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pp2_1_track_gJrT;
    trkUE   = pp2_1_trackUE_gJrT;
    trk_xi     = pp2_1_track_xi_gJrT;
    trkUE_xi   = pp2_1_trackUE_xi_gJrT;
  } 
  if(mode == 38)
  {
    jet     = pp7_1_jet_gen;
    jet_pPb = pPb5_1_jet_gen;
    trk     = pp7_1_track_gJrT;
    trkUE   = pp7_1_trackUE_gJrT;
    trk_xi     = pp7_1_track_xi_gJrT;
    trkUE_xi  = pp7_1_trackUE_xi_gJrT;
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
  if(mode == 41)
  {
    jet     = pp2_0_jet;
    jet_pPb = pp2_0_jet;
    trk     = pp2_0_track;
    trkUE   = pp2_0_trackUE;
    trk_xi     = pp2_0_track_xi;
    trkUE_xi   = pp2_0_trackUE_xi;
  }
  if(mode == 42)
  {
    jet     = pp7_0_jet;
    jet_pPb = pp7_0_jet;
    trk     = pp7_0_track;
    trkUE   = pp7_0_trackUE;
    trk_xi     = pp7_0_track_xi;
    trkUE_xi   = pp7_0_trackUE_xi;
  }
  if(mode == 43)
  {
    jet     = pp2_1_jet_reco;
    jet_pPb = pp2_1_jet_reco;
    trk     = pp2_1_track_reco;
    trkUE   = pp2_1_trackUE_reco;
    trk_xi     = pp2_1_track_xi_reco;
    trkUE_xi   = pp2_1_trackUE_xi_reco;
  }
  if(mode == 44)
  {
    jet     = pp2_1_jet_reco;
    jet_pPb = pp2_1_jet_reco;
    trk     = pp2_1_track_rJgT;
    trkUE   = pp2_1_trackUE_rJgT;
    trk_xi     = pp2_1_track_xi_rJgT;
    trkUE_xi   = pp2_1_trackUE_xi_rJgT;
  }
  if(mode == 45)
  {
    jet     = pp2_1_jet_gen;
    jet_pPb = pp2_1_jet_gen;
    trk     = pp2_1_track_gJrT;
    trkUE   = pp2_1_trackUE_gJrT;
    trk_xi     = pp2_1_track_xi_gJrT;
    trkUE_xi   = pp2_1_trackUE_xi_gJrT;
  }
  if(mode == 46)
  {
    jet     = pp2_1_jet_gen;
    jet_pPb = pp2_1_jet_gen;
    trk     = pp2_1_track_gen;
    trkUE   = pp2_1_trackUE_gen;
    trk_xi     = pp2_1_track_xi_gen;
    trkUE_xi   = pp2_1_trackUE_xi_gen;
  }
  if(mode == 47)
  {
    jet     = pp7_1_jet_reco;
    jet_pPb = pp7_1_jet_reco;
    trk     = pp7_1_track_reco;
    trkUE   = pp7_1_trackUE_reco;
    trk_xi     = pp7_1_track_xi_reco;
    trkUE_xi   = pp7_1_trackUE_xi_reco;
  }
  if(mode == 48)
  {
    jet     = pp7_1_jet_reco;
    jet_pPb = pp7_1_jet_reco;
    trk     = pp7_1_track_rJgT;
    trkUE   = pp7_1_trackUE_rJgT;
    trk_xi     = pp7_1_track_xi_rJgT;
    trkUE_xi   = pp7_1_trackUE_xi_rJgT;
  }
  if(mode == 49)
  {
    jet     = pp7_1_jet_gen;
    jet_pPb = pp7_1_jet_gen;
    trk     = pp7_1_track_gJrT;
    trkUE   = pp7_1_trackUE_gJrT;
    trk_xi     = pp7_1_track_xi_gJrT;
    trkUE_xi   = pp7_1_trackUE_xi_gJrT;
  }
  if(mode == 50)
  {
    jet     = pp7_1_jet_gen;
    jet_pPb = pp7_1_jet_gen;
    trk     = pp7_1_track_gen;
    trkUE   = pp7_1_trackUE_gen;
    trk_xi     = pp7_1_track_xi_gen;
    trkUE_xi   = pp7_1_trackUE_xi_gen;
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
  if(mode == 55)
  {
    jet     = Pbp5_1_jet_reco_Q;
    jet_pPb = Pbp5_1_jet_reco_Q;
    trk     = Pbp5_1_track_reco_Q;
    trkUE   = Pbp5_1_trackUE_reco_Q;
    trk_xi     = Pbp5_1_track_xi_reco_Q;
    trkUE_xi   = Pbp5_1_trackUE_xi_reco_Q;
  }
  if(mode == 56)
  {
    jet     = Pbp5_1_jet_reco_G;
    jet_pPb = Pbp5_1_jet_reco_G;
    trk     = Pbp5_1_track_reco_G;
    trkUE   = Pbp5_1_trackUE_reco_G;
    trk_xi     = Pbp5_1_track_xi_reco_G;
    trkUE_xi   = Pbp5_1_trackUE_xi_reco_G;
  }
  if(mode == 57)
  {
    jet     = Pbp5_1_jet_gen_Q;
    jet_pPb = Pbp5_1_jet_gen_Q;
    trk     = Pbp5_1_track_gen_Q;
    trkUE   = Pbp5_1_trackUE_gen_Q;
    trk_xi     = Pbp5_1_track_xi_gen_Q;
    trkUE_xi   = Pbp5_1_trackUE_xi_gen_Q;
  }
  if(mode == 58)
  {
    jet     = Pbp5_1_jet_gen_G;
    jet_pPb = Pbp5_1_jet_gen_G;
    trk     = Pbp5_1_track_gen_G;
    trkUE   = Pbp5_1_trackUE_gen_G;
    trk_xi     = Pbp5_1_track_xi_gen_G;
    trkUE_xi   = Pbp5_1_trackUE_xi_gen_G;
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
  if(mode == 63)
  {
    jet     = pPb5Pbp5_1_jet_reco_Q;
    jet_pPb = pPb5Pbp5_1_jet_reco_Q;
    trk     = pPb5Pbp5_1_track_reco_Q;
    trkUE   = pPb5Pbp5_1_trackUE_reco_Q;
    trk_xi     = pPb5Pbp5_1_track_xi_reco_Q;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_reco_Q;
  }
  if(mode == 64)
  {
    jet     = pPb5Pbp5_1_jet_reco_G;
    jet_pPb = pPb5Pbp5_1_jet_reco_G;
    trk     = pPb5Pbp5_1_track_reco_G;
    trkUE   = pPb5Pbp5_1_trackUE_reco_G;
    trk_xi     = pPb5Pbp5_1_track_xi_reco_G;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_reco_G;
  }
  if(mode == 65)
  {
    jet     = pPb5Pbp5_1_jet_gen_Q;
    jet_pPb = pPb5Pbp5_1_jet_gen_Q;
    trk     = pPb5Pbp5_1_track_gen_Q;
    trkUE   = pPb5Pbp5_1_trackUE_gen_Q;
    trk_xi     = pPb5Pbp5_1_track_xi_gen_Q;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_gen_Q;
  }
  if(mode == 66)
  {
    jet     = pPb5Pbp5_1_jet_gen_G;
    jet_pPb = pPb5Pbp5_1_jet_gen_G;
    trk     = pPb5Pbp5_1_track_gen_G;
    trkUE   = pPb5Pbp5_1_trackUE_gen_G;
    trk_xi     = pPb5Pbp5_1_track_xi_gen_G;
    trkUE_xi   = pPb5Pbp5_1_trackUE_xi_gen_G;
  }
  return;
}

