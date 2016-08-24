#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TF1.h"
#include <iostream>

void FitGauss(TH1D* hist_p, Bool_t isPbPb, float& mean, float& meanErr, float& res, float& resErr)
{
  if(hist_p->Integral() == 0) return;
  if(hist_p->GetEntries() == 0) return;

  TF1* f1_p = new TF1("f1_p", "gaus", hist_p->GetXaxis()->GetXmin(), hist_p->GetXaxis()->GetXmax());

  hist_p->Fit("f1_p", "Q M");

  mean = f1_p->GetParameter(1);
  res = f1_p->GetParameter(2);

  meanErr = f1_p->GetParError(1);
  resErr = f1_p->GetParError(2);

  //  if(TMath::Abs(mean - 1.0) < 0.01) return;
  if(f1_p->GetProb() > .01) return;
  if(!isPbPb) return;

  for(Int_t fitIter = 0; fitIter < 1; fitIter++){
    Float_t max = -1;
    Int_t maxBin = -1;

    for(Int_t iter = 0; iter < hist_p->GetNbinsX(); iter++){
      if(hist_p->GetBinContent(iter+1) > max){
	max = hist_p->GetBinContent(iter+1);
	maxBin = iter+1;
      }
    }
    
    Int_t nBins = 1;

    Float_t fitLow = -1;
    Float_t fitHi = -1;
    
    for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
      Int_t tempBin = maxBin - iter;
      if(tempBin < 1) tempBin = 1;
      
      nBins += 2;
      
      if(hist_p->Integral(tempBin, maxBin+iter) > (.95-.05*fitIter)*hist_p->Integral() || tempBin == 1){
	fitLow = hist_p->GetBinCenter(tempBin);
	fitHi = hist_p->GetBinCenter(maxBin+iter);
	break;
      }
    }
    
    if(nBins >= 7){
      hist_p->Fit("f1_p", "Q M", "", fitLow, fitHi);
      
      mean = f1_p->GetParameter(1);
      res = f1_p->GetParameter(2);
      
      meanErr = f1_p->GetParError(1);
      resErr = f1_p->GetParError(2);
      
      //  if(TMath::Abs(mean - 1.0) < 0.01) return;
      if(f1_p->GetProb() > .01) return;
      return;
    }
    nBins = 1;
    
    Int_t meanBin = hist_p->FindBin(hist_p->GetMean());
    fitLow = -1;
    fitHi = -1;
    
    for(Int_t iter = 1; iter < hist_p->GetNbinsX(); iter++){
      Int_t tempBin = meanBin - iter;
      if(tempBin < 1) tempBin = 1;
      
      nBins += 2;
      
      if(hist_p->Integral(tempBin, meanBin+iter) > (.95-.05*fitIter)*hist_p->Integral() || tempBin == 1){
	fitLow = hist_p->GetBinCenter(tempBin);
	fitHi = hist_p->GetBinCenter(meanBin+iter);
	break;
      }
    } 
    
    if(nBins >= 7){
      hist_p->Fit("f1_p", "Q M", "", fitLow, fitHi);
      
      mean = f1_p->GetParameter(1);
      res = f1_p->GetParameter(2);
      
      meanErr = f1_p->GetParError(1);
      resErr = f1_p->GetParError(2);
    }
  }

  delete f1_p;
  return;
}

void checkJESJER(){
  TFile * out = TFile::Open("checkJESJER.root","recreate");
  TFile * f = TFile::Open("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia6_Dijet170_pp5TeV/0.root","read");
  TTree * t = (TTree*)f->Get("ak3PFJetAnalyzer/t");

  int nref;
  int ngen;
  float rawpt[1000];
  float chargedSum[1000];
  float chargedMax[1000];
  float jtpt[1000];
  float jteta[1000];
  float jtphi[1000];
  float refpt[1000];
  float refeta[1000];
  float refphi[1000];
  float genpt[1000];
  float geneta[1000];
  float genphi[1000];

  t->SetBranchAddress("nref",&nref);
  t->SetBranchAddress("ngen",&ngen);
  t->SetBranchAddress("rawpt",&rawpt);
  t->SetBranchAddress("chargedSum",&chargedSum);
  t->SetBranchAddress("chargedMax",&chargedMax);
  t->SetBranchAddress("jtpt",&jtpt);
  t->SetBranchAddress("jteta",&jteta);
  t->SetBranchAddress("jtphi",&jtphi);
  t->SetBranchAddress("refpt",&refpt);
  t->SetBranchAddress("refeta",&refeta);
  t->SetBranchAddress("refphi",&refphi);
  t->SetBranchAddress("genpt",&genpt);
  t->SetBranchAddress("geneta",&geneta);
  t->SetBranchAddress("genphi",&genphi);
 
  TH1D * res[5];
  for(int i = 0; i<5; i++){
    res[i] = new TH1D(Form("res%d",i),"",30,0.6,1.4);
  }
 
  for(int i = 0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    for(int j = 0; j<nref; j++){
      if(TMath::Abs(jteta[j])>1.5) continue;
      if(rawpt[j]<30) continue;
      if(chargedSum[j]/rawpt[j]<0.05 || chargedSum[j]/rawpt[j]>0.95) continue;

      int bin = -1;
      /*if(jtpt[j]>=60) bin=0;
      if(jtpt[j]>=80) bin=1;
      if(jtpt[j]>=100) bin=2;
      if(jtpt[j]>=120) bin=3;
      if(jtpt[j]>=140) bin=4;
      if(jtpt[j]<60 || jtpt[j]>200) bin=-1;*/
      if(refpt[j]>=60) bin=0;
      if(refpt[j]>=80) bin=1;
      if(refpt[j]>=100) bin=2;
      if(refpt[j]>=120) bin=3;
      if(refpt[j]>=140) bin=4;
      if(refpt[j]<60 || refpt[j]>200) bin=-1;
      if(bin == -1) continue;
      res[bin]->Fill(jtpt[j]/refpt[j]); 
    }
  }
  
  const int jtptBin = 5;
  double jtptBins[jtptBin+1] = {60,80,100,120,140,200};
  TH1D * meanPP = new TH1D("meanPP","",jtptBin,jtptBins);
  TH1D * resoPP = new TH1D("resoPP","",jtptBin,jtptBins);
  for(int j = 0; j<5; j++){
    float mean=0, meanErr=0, reso=0, resoErr=0;      
    FitGauss(res[j],0,mean,meanErr,reso,resoErr);
    std::cout << mean << " " << meanErr << " " << reso << " " << resoErr << std::endl;
    meanPP->SetBinContent(j+1,mean);
    meanPP->SetBinError(j+1,meanErr);
    resoPP->SetBinContent(j+1,reso);
    resoPP->SetBinError(j+1,resoErr);
  }
  meanPP->GetYaxis()->SetRangeUser(0.9,1.1);
  resoPP->GetYaxis()->SetRangeUser(-0.01,0.2);
  out->Write();
}
