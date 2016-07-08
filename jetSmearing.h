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
