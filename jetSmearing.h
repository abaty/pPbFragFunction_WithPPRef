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

TH1D * ppAdditionalSmearingTopPb;

void initializeJER()
{
  /*TFile * JERFile = TFile::Open("jetJER.root","read");
  JER_pp2 = (TH1D*) JERFile->Get("pp2_JER"); 
  JER_pp7 = (TH1D*) JERFile->Get("pp7_JER");
  JER_pPb5 = (TH1D*) JERFile->Get("pPb5_JER");
  */ 

  gaussJER = new TF1("gausJER","gaus(0)",-20,20);
  gaussJER->SetParameters(1,0,1);

  
  ppAdditionalSmearingTopPb = new TH1D("ppAdditionalSmearingTopPb","ppAdditionalSmearingTopPb",60,50,300);
   ppAdditionalSmearingTopPb->SetBinContent(1,0.06696667);
   ppAdditionalSmearingTopPb->SetBinContent(2,0.072717);
   ppAdditionalSmearingTopPb->SetBinContent(3,0.07693734);
   ppAdditionalSmearingTopPb->SetBinContent(4,0.07707438);
   ppAdditionalSmearingTopPb->SetBinContent(5,0.07214285);
   ppAdditionalSmearingTopPb->SetBinContent(6,0.06624958);
   ppAdditionalSmearingTopPb->SetBinContent(7,0.06231565);
   ppAdditionalSmearingTopPb->SetBinContent(8,0.05955913);
   ppAdditionalSmearingTopPb->SetBinContent(9,0.05919462);
   ppAdditionalSmearingTopPb->SetBinContent(10,0.05453558);
   ppAdditionalSmearingTopPb->SetBinContent(11,0.05111116);
   ppAdditionalSmearingTopPb->SetBinContent(12,0.04068843);
   ppAdditionalSmearingTopPb->SetBinContent(13,0.03708558);
   ppAdditionalSmearingTopPb->SetBinContent(14,0.03192078);
   ppAdditionalSmearingTopPb->SetBinContent(15,0.03388508);
   ppAdditionalSmearingTopPb->SetBinContent(16,0.03450703);
   ppAdditionalSmearingTopPb->SetBinContent(17,0.03315393);
   ppAdditionalSmearingTopPb->SetBinContent(18,0.03027993);
   ppAdditionalSmearingTopPb->SetBinContent(19,0.01909295);
   ppAdditionalSmearingTopPb->SetBinContent(20,0.01564647);
   ppAdditionalSmearingTopPb->SetBinContent(21,0.01408416);
   ppAdditionalSmearingTopPb->SetBinContent(22,0.01763726);
   ppAdditionalSmearingTopPb->SetBinContent(23,0.01468521);
   ppAdditionalSmearingTopPb->SetBinContent(24,0.01166856);
   ppAdditionalSmearingTopPb->SetBinContent(25,0.01014812);
   ppAdditionalSmearingTopPb->SetBinContent(26,0.01132335);
   ppAdditionalSmearingTopPb->SetBinContent(27,0.01188625);
   ppAdditionalSmearingTopPb->SetBinContent(28,0.01344728);
   ppAdditionalSmearingTopPb->SetBinContent(29,0.01448107);
   ppAdditionalSmearingTopPb->SetBinContent(30,0.01732975);
   ppAdditionalSmearingTopPb->SetBinContent(31,0.01910914);
   ppAdditionalSmearingTopPb->SetBinContent(32,0.02213025);
   ppAdditionalSmearingTopPb->SetBinContent(33,0.02487918);
   ppAdditionalSmearingTopPb->SetBinContent(34,0.02779003);
   ppAdditionalSmearingTopPb->SetBinContent(35,0.02905112);
   ppAdditionalSmearingTopPb->SetBinContent(36,0.0305647);
   ppAdditionalSmearingTopPb->SetBinContent(37,0.03106759);
   ppAdditionalSmearingTopPb->SetBinContent(38,0.03289593);
   ppAdditionalSmearingTopPb->SetBinContent(39,0.03307773);
   ppAdditionalSmearingTopPb->SetBinContent(40,0.03184452);
   ppAdditionalSmearingTopPb->SetBinContent(41,0.02916415);
   ppAdditionalSmearingTopPb->SetBinContent(42,0.02385031);
   ppAdditionalSmearingTopPb->SetBinContent(43,0.007724378);
   ppAdditionalSmearingTopPb->SetBinContent(44,-0.01007987);
   ppAdditionalSmearingTopPb->SetBinContent(45,-0.02867151);
   ppAdditionalSmearingTopPb->SetBinContent(46,-0.0366714);
   ppAdditionalSmearingTopPb->SetBinContent(47,-0.04397122);
   ppAdditionalSmearingTopPb->SetBinContent(48,-0.04695244);
   ppAdditionalSmearingTopPb->SetBinContent(49,-0.05122457);
   ppAdditionalSmearingTopPb->SetBinContent(50,-0.05443871);
   ppAdditionalSmearingTopPb->SetBinContent(51,-0.05793289);
   ppAdditionalSmearingTopPb->SetBinContent(52,-0.05100638);
   ppAdditionalSmearingTopPb->SetBinContent(53,-0.04600485);
   ppAdditionalSmearingTopPb->SetBinContent(54,-0.04301372);
   ppAdditionalSmearingTopPb->SetBinContent(55,-0.04499522);
   ppAdditionalSmearingTopPb->SetBinContent(56,-0.04151413);
   ppAdditionalSmearingTopPb->SetBinContent(57,-0.03423494);
   ppAdditionalSmearingTopPb->SetBinContent(58,-0.02655439);
   ppAdditionalSmearingTopPb->SetBinContent(59,-0.02186047);
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
  if(!isJER_Initialized)
  {
    gRandom->SetSeed(0);
    initializeJER();
    isJER_Initialized = true;
  }
  if(jtpt<50 || jtpt>220) return 0;
  return ppAdditionalSmearingTopPb->GetBinContent(ppAdditionalSmearingTopPb->FindBin(jtpt));
  /*
  float pp1 = 0.554;
  float pp2 = 0.327;
  float pPb1 = 0.919;
  float pPb2 = 0.429;
  float Pbp1 = 0.565;
  float Pbp2 = 0.337;
  return (((pPb1*TMath::Power(jtpt,-pPb2)*0.71+Pbp1*TMath::Power(jtpt,-Pbp2)*1.29)/2.0-pp1*TMath::Power(jtpt,-pp2)>0)?(pPb1*TMath::Power(jtpt,-pPb2)*0.71+Pbp1*TMath::Power(jtpt,-Pbp2)*1.29)/2.0-pp1*TMath::Power(jtpt,-pp2):0);
  */
}

double getGenMCSmearFactor(const char * mode, float jtpt)
{
  if(jtpt<50) return 0;
  float pp1 = 0.704;
  float pp2 = 0.381;
  float pPb1 = 1.42;
  float pPb2 = 0.521;
  float Pbp1 = 1.21;
  float Pbp2 = 0.501;
  if(strcmp(mode,"ppref5")==0){
    return pp1*TMath::Power(jtpt,-pp2);
  }else if(strcmp(mode,"pPb5")==0){
    return pPb1*TMath::Power(jtpt,-pPb2);
  }else if(strcmp(mode,"Pbp5")==0){
    return Pbp1*TMath::Power(jtpt,-Pbp2);
  }else{
    return 0;
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
