#include "TMath.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TRandom.h"

bool isJER_Initialized = false;
TH1D * JER_pp2;
TH1D * JER_pp7;
TH1D * JER_pPb5;
TF1 * gaussJER;

void initializeJER()
{
  TFile * JERFile = TFile::Open("jetJER.root","read");
  JER_pp2 = (TH1D*) JERFile->Get("pp2_JER"); 
  JER_pp7 = (TH1D*) JERFile->Get("pp7_JER");
  JER_pPb5 = (TH1D*) JERFile->Get("pPb5_JER");
 
  gaussJER = new TF1("gausJER","gaus(0)",-20,20);
  gaussJER->SetParameters(1,0,1);

  return;
}

double getJERCorrected(const char * mode, double jetpt, double percentage = 0.02)
{
  if(!isJER_Initialized)
  {
    gRandom->SetSeed(0);
    initializeJER();
    isJER_Initialized = true;
  }

  /*double res_pp2 = JER_pp2->GetBinContent(JER_pp2->FindBin(jetpt));
  double res_pp7 = JER_pp7->GetBinContent(JER_pp7->FindBin(jetpt));
  double res_pPb5 = JER_pPb5->GetBinContent(JER_pPb5->FindBin(jetpt));

  double max = TMath::Max(TMath::Max(res_pp2,res_pp7),res_pPb5);
  double sigma2 = 0;
  if(strcmp("pp2",mode)==0) sigma2 = jetpt*jetpt*(max-res_pp2);
  else if(strcmp("pp7",mode)==0) sigma2 = jetpt*jetpt*(max-res_pp7);
  else sigma2 = jetpt*jetpt*(max-res_pPb5);

  return jetpt+TMath::Power(sigma2,0.5)*gaussJER->GetRandom();*/
  return jetpt+jetpt*percentage*gaussJER->GetRandom();
}

double getPPDataSmearFactor(float jtpt)
{
  float pPb1 = 2.065;
  float pPb2 = 0.597;
  float pp1 = 0.692;
  float pp2 = 0.3746;
  return (pPb1*TMath::Power(jtpt,-pPb2)-pp1*TMath::Power(jtpt,-pp2)>0)?TMath::Power(TMath::Power(pPb1*TMath::Power(jtpt,-pPb2),2)-TMath::Power(pp1*TMath::Power(jtpt,-pp2),2),0.5):0;
}

double getGenMCSmearFactor(const char * mode, float jtpt)
{
  float pPb1 = 2.065;
  float pPb2 = 0.597;
  float pp1 = 0.692;
  float pp2 = 0.3746;
  if(strcmp(mode,"ppref5")==0){
    return pp1*TMath::Power(jtpt,-pp2);
  }else{
    return pPb1*TMath::Power(jtpt,-pPb2);
  }
}

double getFragJECFactor(const char * mode, float leadingPart, float jtpt)
{
  float z = leadingPart/jtpt;
  float Pbp0 = 0.933;
  float Pbp1 = 0.383;
  float Pbp2 = -0.35;
  float pPb0 = 0.938;
  float pPb1 = 0.374;
  float pPb2 = -0.286;
  float pp0 = 0.9261;
  float pp1 = 0.3784;
  float pp2 = -0.2702;
  if(z>1.2) return jtpt;
  if(strcmp(mode,"ppref5")==0){
    return jtpt/(pp0+pp1*z+pp2*z*z);
  }else if(strcmp(mode,"pPb5")==0){
    return jtpt/(pPb0+pPb1*z+pPb2*z*z);
  }else if(strcmp(mode,"Pbp5")==0){
    return jtpt/(Pbp0+Pbp1*z+Pbp2*z*z);
  }else{
    return jtpt;
  }
}
