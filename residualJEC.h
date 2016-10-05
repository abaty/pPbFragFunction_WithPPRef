#include "TF1.h"

bool isJEC_Initialized = false;
TF1 * residualJEC_pt;
TF1 * residualJEC_MCsmear;
TF1 * gaus;
TH1D * residualJEC_eta;

void initializeResidualJEC(const char* mode)
{
  TFile * residualFile = TFile::Open("residualJEC.root","read");
  if(strcmp(mode,"pp5")!=0)
  {
    residualJEC_pt = (TF1*)residualFile->Get(Form("%sJEC_pt",mode));
    residualJEC_MCsmear = (TF1*)residualFile->Get(Form("%sJEC_MCsmear",mode));
    residualJEC_eta = (TH1D*)residualFile->Get(Form("%sJEC_eta",mode));
  }
  else//if pp5 (outdated)
  {
    residualJEC_pt = (TF1*)residualFile->Get("pPb5JEC_pt");
    residualJEC_MCsmear = (TF1*)residualFile->Get("pPb5JEC_MCsmear");
    residualJEC_eta = (TH1D*)residualFile->Get("pPb5JEC_eta");
  }
  gaus = new TF1("gaus","gaus(0)",-20,20);
  gaus->SetParameters(1,0,1);
  return;
}

double getCorrectedJetPt(const char * mode, bool isMC, double jetpt, double jeteta)
{
  if(!isJEC_Initialized)
  {
    initializeResidualJEC(mode);
    isJEC_Initialized = true;
  }

  double correction = 1;
  correction = residualJEC_eta->GetBinContent(residualJEC_eta->FindBin(jeteta));
  correction = correction*residualJEC_pt->Eval(jetpt); 

  if(strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0) correction = correction*0.995;//for UE activity residual that Yue Shi didn't account for (Doga did though)
 
  if(isMC)
  {
    correction = 1;
    //correction = 1+residualJEC_MCsmear->Eval(jetpt)*gaus->GetRandom();
  }
  return correction*jetpt;
}
