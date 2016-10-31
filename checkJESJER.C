#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TF1.h"
#include "TLegend.h"
#include "TProfile.h"
//#include "residualJEC.h"
//#include "jetSmearing.h"
//#include "getJEC_2nd.h"
//#include "getJEC_1st.h"
//#include "getJEC_L2L3res.h"
//#include "getJEC_SystError.h"

#include "L2L3ResidualWFits.h"
#include "MCTruthResidual.h"

#include "jetSmearing.h"
#include <cmath>
#include <iostream>

double getdR2(double jet_eta, double jet_phi, double track_eta, double track_phi)
{
  return TMath::Power(jet_eta-track_eta,2)+TMath::Power(acos(cos(jet_phi-track_phi)),2);
}

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

void checkJESJER(bool doMoreBins = false, bool doFragJEC = true, bool doFineBins = false, bool dojtpt = false, bool showMore = false){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();

  TH1D * res[60][4];
  TH1D * res_hard[60][4];
  TH1D * meanPP[4];
  TH1D * resoPP[4];
  TH1D * meanPP_hard[4];
  TH1D * resoPP_hard[4];

  int pthats[7] = {30,50,80,120,170,220,280};
  double crossSection5[7]  = {3.378E-02,3.778E-03,4.412E-04,6.147E-05,1.018E-05,2.477E-06,6.160E-07};

  TFile * out = TFile::Open("checkJESJER.root","recreate");
  TH1D * pthatP[3];
  pthatP[0]  = new TH1D("pthatp0","",100,0,200);
  
  TH2D * h_trackVsJEC[3];
  TH2D * h_trackVsJEC_weights[3];
  TH2D * h_trackVsJEC_z[3];
  TH2D * h_trackVsJEC_weights_z[3];
  TH1D * h_trackVsJEC_1z[3][3];
  TH1D * h_trackVsJEC_weights_1z[3][3];
  TProfile * h_fakeVsChargeSum[3];
  TH1D * h_neutralSum[3];
  const int ntrackBins=23;
  const double axis[ntrackBins] = {0.5, 0.63, 0.77,  1.03,1.38, 1.84, 2.46, 3.29,  4.40, 5.88,  7.87,  10.52, 14.06,  18.8, 25.13,  33.58,  44.89,  60, 80, 100, 120, 140, 200};
  const int ntrackBins_z=11;
  const double axis_z[ntrackBins_z] = {0.0,0.05,0.1,0.15,0.2,0.3,0.4,0.6,0.8,1,1.2};
  const int trkVsJEC_jetBins = 9;
  const double trkVsJEC_jetBinsAxis[trkVsJEC_jetBins+1] = {50,60,80,100,120,140,160,180,200,220};
  for(int i = 0; i<3; i++){
    h_trackVsJEC[i] = new TH2D(Form("h_trackVsJEC%d",i),"",ntrackBins-1,axis,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_trackVsJEC_weights[i] = new TH2D(Form("h_trackVsJEC_weights%d",i),"",ntrackBins-1,axis,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_trackVsJEC_z[i] = new TH2D(Form("h_trackVsJEC%d_z",i),"",ntrackBins_z-1,axis_z,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_trackVsJEC_weights_z[i] = new TH2D(Form("h_trackVsJEC_weights%d_z",i),"",ntrackBins_z-1,axis_z,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    for(int k = 0; k<3; k++){
      h_trackVsJEC_1z[i][k] = new TH1D(Form("h_trackVsJEC%d_%d_1z",i,k),"",ntrackBins_z-1,axis_z);
      h_trackVsJEC_weights_1z[i][k] = new TH1D(Form("h_trackVsJEC_weights%d_%d_1z",i,k),"",ntrackBins_z-1,axis_z);
    }
    h_fakeVsChargeSum[i] = new TProfile(Form("h_fakeVsChargeSum_%d",i),"",25,0,1.25);
    h_neutralSum[i] = new TH1D(Form("h_neutralSum_%d",i),"",25,0,1.25);
  }
  
  for(int file = 0; file<4; file++){
  if(file<3){//file==3 skip to plotting step
  std::cout << file << std::endl;
  for(int file2= 0; file2<5; file2++){
  std::cout << file2 << std::endl;
//for pp
  int radius = 3;
  double etacutcorr = 3;
  TString mode;
  if(file==0) mode = "pp5";
  if(file==1) mode = "pPb5";
  if(file==2) mode = "Pbp5";
  if(file==3) mode = "pp5";
  std::cout << mode << std::endl;
  L2ResidualJES * L2JES = new L2ResidualJES(radius,((int)etacutcorr),mode);
  L2ResidualJER * L2JER = new L2ResidualJER(mode);
  L3ResidualJES * L3JES = new L3ResidualJES(mode);
  MCTruthResidual * MCTruth = new MCTruthResidual(mode);
 

  TFile * f;
  if(file2==0 && (!doFineBins || file==0 || file==4)) continue;
  //if(file2>=4) continue;
  //if(file==0) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/mergedForests/PYTHIA_QCD_TuneCUETP8M1_cfi_GEN_SIM_5020GeV_ppSignal/Pythia8_Dijet%d_pp_TuneCUETP8M1_5020GeV_FOREST_758_PrivMC/0.root",pthats[file2]),"read");//pythia 8
  if(file==0) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia6_Dijet%d_pp5TeV/0.root",pthats[file2]),"read");//pythia 6
//for pPb
  if(doFragJEC){
    if(file==1) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%d/HiForest_v77_merged01/pt%d_HP04_prod16_v77_merged_forest_0.root",pthats[file2],pthats[file2]),"read");
    if(file==2 && file2<4) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/Pbp/HP05/prod24/Hijing_Pythia_pt%d/HiForest_v84_merged02/pt%d_HP05_prod24_v84_merged_forest_0.root",pthats[file2],pthats[file2]),"read");
  }else{
    if(file==1) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForestpPb/MC/HiForest_pPbDir_QCD%d_pythia6_HiWinter13.root",pthats[file2]),"read");
    if(file==2) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForestpPb/MC/HiForest_PbpDir_QCD%d_pythia6_HiWinter13.root",pthats[file2]),"read");
  }
  if(doFragJEC && file==2 && file2>=4) continue;
  TTree * t = (TTree*)f->Get("ak3PFJetAnalyzer/t");
  TTree * hi;
  if(file!=0) hi = (TTree*)f->Get("hiEvtAnalyzer/HiTree");
  TTree * trk;
  if(doFragJEC) trk = (TTree*)f->Get("ppTrack/trackTree");

  int nref;
  int ngen;
  float rawpt[1000];
  float chargedSum[1000];
  float neutralSum[1000];
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
  float pthat;

  float hiHF4;

  int nTrk;
  float trkPt[3000];
  float trkEta[3000];
  float trkPhi[3000];
  float trkPtError[3000];
  float trkDz1[3000];
  float trkDzError1[3000];
  float trkDxy1[3000];
  float trkDxyError1[3000];
  bool highPurity[3000];

  t->SetBranchAddress("nref",&nref);
  t->SetBranchAddress("ngen",&ngen);
  t->SetBranchAddress("pthat",&pthat);
  t->SetBranchAddress("rawpt",&rawpt);
  t->SetBranchAddress("chargedSum",&chargedSum);
  t->SetBranchAddress("neutralSum",&neutralSum);
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
 
  if(doFragJEC){ 
  trk->SetBranchAddress("nTrk",&nTrk);
  trk->SetBranchAddress("trkPt",&trkPt);
  trk->SetBranchAddress("trkEta",&trkEta);
  trk->SetBranchAddress("trkPhi",&trkPhi);
  /*trk->SetBranchAddress("trkPtError",&trkPtError);
  trk->SetBranchAddress("trkDz1",&trkDz1);
  trk->SetBranchAddress("trkDzError1",&trkDzError1);
  trk->SetBranchAddress("trkDxy1",&trkDxy1);
  trk->SetBranchAddress("trkDxyError1",&trkDxyError1);*/
  trk->SetBranchAddress("highPurity",&highPurity);
  }

  if(file==1) hi->SetBranchAddress("hiHFminusEta4",&hiHF4);
  if(file==2) hi->SetBranchAddress("hiHFplusEta4",&hiHF4);
  if(file==1 || file==2) t->AddFriend(hi);
  if(doFragJEC) t->AddFriend(trk);
 
  for(int i = 0; i< (doFineBins?60:(doMoreBins?9:5)) ; i++){
    std::cout << i << std::endl;
    res[i][file] = new TH1D(Form("res%d%d",i,file),"",30,0.6,1.4);
    res_hard[i][file] = new TH1D(Form("res_hard%d%d",i,file),"",30,0.6,1.4);
    if(file==0){
      res[i][3] = new TH1D(Form("res%d%d",i,3),"",30,0.6,1.4);
      res_hard[i][3] = new TH1D(Form("res_hard%d%d",i,3),"",30,0.6,1.4);
    }
  }
 
  for(int i = 0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(pthat>pthats[file2+1]) continue;
    if(file==0) pthatP[0]->Fill(pthat,crossSection5[file2]/t->GetEntries());
    if(i%1000==0) std::cout << i << "/" << t->GetEntries() << std::endl;
    for(int j = 0; j<nref; j++){
      float boost = 0;
      //if(pthat > 2*pthats[file2]) continue;
      if(file==1) boost = 0.4654094531;
      if(file==2) boost = -0.4654094531;
      boost = -boost;//pPb MC has proton in + direction, pPb data has it in minus, and both are reversed for the 'Pbp definition'
      if(TMath::Abs(jteta[j]+boost)>1.5) continue;
      if(jtpt[j]>140 && jtpt[j]<200 && refpt[j]<=0 && rawpt[j]>30) std::cout << jtpt[j] << " " << refpt[j] << " " << chargedMax[j]/rawpt[j] << " " << chargedSum[j]/rawpt[j] << std::endl;
      if(rawpt[j]<30) continue;
      //if(chargedSum[j]/rawpt[j]<0.05 || chargedSum[j]/rawpt[j]>0.95) continue;
      if(chargedSum[j]/rawpt[j]<0.05) continue;

  double rawptT = rawpt[j];
  double jtptT = jtpt[j];
  double jtetaT = jteta[j];
  double correctedPt = jtptT;
  bool doMC = true;
  /*rawptT=100;
  jtptT=100;
  jtetaT=1;
  correctedPt = jtptT;*/

  //if(i%50000==0 && j==0) cout << "rawpt = " << rawpt[j] << " MCpt = " << refpt[j] << endl;
  //JEC part 1
  if(mode == "pPb5" || mode == "Pbp5"){
    correctedPt = MCTruth->getJEC_1st(rawptT,correctedPt,jtetaT); 
  }
  correctedPt = MCTruth->getResidualCorr(correctedPt,jtetaT);
  if((!(mode == "pPb5" || mode == "Pbp5")) && doMC) correctedPt = correctedPt*1.008;
  if(i%50000==0 && j==0) cout << "after MC correction jtpt = " << correctedPt << endl;
  //JEC part 2
  if(!doMC){
   correctedPt = L2JES->getCorrectedPt(correctedPt, jtetaT);
   correctedPt = L3JES->getCorrectedPt(correctedPt);
    if(i%50000==0 && j==0) cout << "after data driven correction jtpt = " << correctedPt << endl;
  }else{   
    correctedPt = L2JER->getSmearedPt(correctedPt, jtetaT);
  }
  
  //fragmentation stuff
  float leadingHadronPt = 0;
  if(doFragJEC){
  for(int t=0; t<nTrk; t++){
    if((TMath::Abs(trkEta[t]-jteta[j])>0.3) || (TMath::Abs(trkPhi[t]-jtphi[j])>0.3)) continue;
    if(trkPt[t]<leadingHadronPt) continue;
    if(trkPt[t] <= 1 || !highPurity[t] || TMath::Abs(trkEta[t])>2.4 ) continue;
    //if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.1)) continue;            
    //else if(strcmp(mode,"ppref5")==0 && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.3 || ((pfEcal[t]+pfHcal[t])/TMath::CosH(trkEta[t])<0.5*trkPt[t] && trkPt[t]>20))) continue;            
    if(getdR2(jteta[j]+boost,jtphi[j],trkEta[t]+boost,trkPhi[t]) < 0.3*0.3){
      leadingHadronPt = trkPt[t];
    }
  }

  //FF experiment:
  
  /*if(file==0 || file==4) jtpt[j] = getFragJECFactor("ppref5",leadingHadronPt,jtpt[j]);
  if(file==1) jtpt[j] = getFragJECFactor("pPb5",leadingHadronPt,jtpt[j]);
  if(file==2) jtpt[j] = getFragJECFactor("Pbp5",leadingHadronPt,jtpt[j]);
    */
  
  if(doMC && refpt[j]>20){
    float weight = crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4);
    h_trackVsJEC[file]->Fill(leadingHadronPt,refpt[j],jtpt[j]/refpt[j]*weight);
    h_trackVsJEC_weights[file]->Fill(leadingHadronPt,refpt[j],weight);
    h_trackVsJEC_z[file]->Fill(leadingHadronPt/refpt[j],refpt[j],jtpt[j]/refpt[j]*weight);
    h_trackVsJEC_weights_z[file]->Fill(leadingHadronPt/refpt[j],refpt[j],weight);
    if(refpt[j]<100 && refpt[j]>=60){
      h_trackVsJEC_1z[file][0]->Fill(leadingHadronPt/refpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][0]->Fill(leadingHadronPt/refpt[j],weight);
    }
    if(refpt[j]<140 && refpt[j]>=100){
      h_trackVsJEC_1z[file][1]->Fill(leadingHadronPt/refpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][1]->Fill(leadingHadronPt/refpt[j],weight);
    }
    if(refpt[j]<200 && refpt[j]>=140){
      h_trackVsJEC_1z[file][2]->Fill(leadingHadronPt/refpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][2]->Fill(leadingHadronPt/refpt[j],weight);
    }
    /*if(jtpt[j]<100 && jtpt[j]>=60){
      h_trackVsJEC_1z[file][0]->Fill(leadingHadronPt/jtpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][0]->Fill(leadingHadronPt/jtpt[j],weight);
    }
    if(jtpt[j]<140 && jtpt[j]>=100){
      h_trackVsJEC_1z[file][1]->Fill(leadingHadronPt/jtpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][1]->Fill(leadingHadronPt/jtpt[j],weight);
    }
    if(jtpt[j]<200 && jtpt[j]>=140){
      h_trackVsJEC_1z[file][2]->Fill(leadingHadronPt/jtpt[j],jtpt[j]/refpt[j]*weight);
      h_trackVsJEC_weights_1z[file][2]->Fill(leadingHadronPt/jtpt[j],weight);
    }*/
  }
  }
  //FF experiment:
  /*if(file==0 || file==4) jtpt[j] = getFragJECFactor("ppref5",leadingHadronPt,jtpt[j]);
  if(file==1) jtpt[j] = getFragJECFactor("pPb5",leadingHadronPt,jtpt[j]);
  if(file==2) jtpt[j] = getFragJECFactor("Pbp5",leadingHadronPt,jtpt[j]);*/
  
  jtpt[j] = correctedPt;

  //smearing pp to match pPb
  float smearedJtPt = jtpt[j];
  if(file==0) smearedJtPt = getJERCorrected("pp5",jtpt[j],getPPDataSmearFactor(jtpt[j]));
  
  /*if(jtpt[j]>60 && jtpt[j]<80){
    h_fakeVsChargeSum[file]->Fill(leadingHadronPt/jtpt[j],((refpt[j]<0)?1:0),crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4));
    h_neutralSum[file]->Fill(neutralSum[j]/rawpt[j],crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4));
  }*/

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
      if(refpt[j]>=200) bin=5;
      if(refpt[j]>=260) bin=6;
      if(refpt[j]>=320) bin=7;
      if(refpt[j]>=380) bin=8;
      if(refpt[j]<60 || refpt[j]>440 || (!doMoreBins && refpt[j]>=200)) bin=-1;

      if(doFineBins){
        for(int b = 0; b<60; b++){
          if(refpt[j]>=50+250.0/60.0*b) bin=b;
          if(refpt[j]<50 || refpt[j]>=300) bin=-1;
        }
      }

     if(dojtpt){
      if(jtpt[j]>=60) bin=0;
      if(jtpt[j]>=80) bin=1;
      if(jtpt[j]>=100) bin=2;
      if(jtpt[j]>=120) bin=3;
      if(jtpt[j]>=140) bin=4;
      if(jtpt[j]>=200) bin=5;
      if(jtpt[j]>=260) bin=6;
      if(jtpt[j]>=320) bin=7;
      if(jtpt[j]>=380) bin=8;
      if(jtpt[j]<60 || jtpt[j]>440 || (!doMoreBins && jtpt[j]>=200)) bin=-1;

      if(doFineBins){
        for(int b = 0; b<60; b++){
          if(jtpt[j]>=50+250.0/60.0*b) bin=b;
          if(jtpt[j]<50 || jtpt[j]>=300) bin=-1;
          }
        }
      }

      if(bin == -1) continue;
      if(file==1 || file==2){
        res[bin][file]->Fill(jtpt[j]/refpt[j],crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4));
      }else{
        res[bin][file]->Fill(jtpt[j]/refpt[j],crossSection5[file2]/t->GetEntries());
        res[bin][3]->Fill(smearedJtPt/refpt[j],crossSection5[file2]/t->GetEntries());
      }

        if(chargedMax[j]/jtpt[j]<0.5) continue;
        else{
          if(file==1 || file==2) res_hard[bin][file]->Fill(jtpt[j]/refpt[j],crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4)); 
          else{
            res_hard[bin][0]->Fill(smearedJtPt/refpt[j],crossSection5[file2]/t->GetEntries()); 
            res_hard[bin][3]->Fill(smearedJtPt/refpt[j],crossSection5[file2]/t->GetEntries()); 
          }
        }
      }
      //std::cout << "test" << std::endl;
    }
  res[0][file]->Print("All");
  if(file==0) res[0][3]->Print("All");
  }//end inner file loop
  }
  
  const int jtptBin = 5;
  double jtptBins[jtptBin+1] = {60,80,100,120,140,200};
  const int jtptBinMore = 9;
  double jtptBinsMore[jtptBinMore+1] = {60,80,100,120,140,200,260,320,380,440};
  const int jtptBinFine = 60;
  if(!doFineBins){
    meanPP[file] = new TH1D(Form("meanPP%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
    resoPP[file] = new TH1D(Form("resoPP%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
    meanPP_hard[file] = new TH1D(Form("meanPP_hard%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
    resoPP_hard[file] = new TH1D(Form("resoPP_hard%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
  }else{
    meanPP[file] = new TH1D(Form("meanPP%d",file),"",jtptBinFine,50,300);
    resoPP[file] = new TH1D(Form("resoPP%d",file),"",jtptBinFine,50,300);
    meanPP_hard[file] = new TH1D(Form("meanPP_hard%d",file),"",jtptBinFine,50,300);
    resoPP_hard[file] = new TH1D(Form("resoPP_hard%d",file),"",jtptBinFine,50,300);
  }
  for(int j = 0; j<(doFineBins?60:(doMoreBins?9:5)); j++){
    float mean=0, meanErr=0, reso=0, resoErr=0;      
    FitGauss(res[j][file],0,mean,meanErr,reso,resoErr);
    std::cout << mean << " " << meanErr << " " << reso << " " << resoErr << std::endl;
    meanPP[file]->SetBinContent(j+1,mean);
    meanPP[file]->SetBinError(j+1,meanErr);
    resoPP[file]->SetBinContent(j+1,reso);
    resoPP[file]->SetBinError(j+1,resoErr);
  }
  for(int j = 0; j<(doFineBins?60:(doMoreBins?9:5)); j++){
    float mean=0, meanErr=0, reso=0, resoErr=0;      
    FitGauss(res_hard[j][file],0,mean,meanErr,reso,resoErr);
    std::cout << mean << " " << meanErr << " " << reso << " " << resoErr << std::endl;
    meanPP_hard[file]->SetBinContent(j+1,mean);
    meanPP_hard[file]->SetBinError(j+1,meanErr);
    resoPP_hard[file]->SetBinContent(j+1,reso);
    resoPP_hard[file]->SetBinError(j+1,resoErr);
  }
  meanPP[file]->GetYaxis()->SetRangeUser(0.9,1.1);
  resoPP[file]->GetYaxis()->SetRangeUser(-0.01,0.2);
  meanPP_hard[file]->GetYaxis()->SetRangeUser(0.9,1.1);
  resoPP_hard[file]->GetYaxis()->SetRangeUser(-0.01,0.2);

    if(file<3){
      h_trackVsJEC[file]->Divide(h_trackVsJEC_weights[file]);
      h_trackVsJEC[file]->Print("All");
      h_trackVsJEC_z[file]->Divide(h_trackVsJEC_weights_z[file]);
      h_trackVsJEC_z[file]->Print("All");
      h_trackVsJEC_1z[file][0]->Divide(h_trackVsJEC_weights_1z[file][0]);
      h_trackVsJEC_1z[file][0]->Print("All");
      h_trackVsJEC_1z[file][1]->Divide(h_trackVsJEC_weights_1z[file][1]);
      h_trackVsJEC_1z[file][1]->Print("All");
      h_trackVsJEC_1z[file][2]->Divide(h_trackVsJEC_weights_1z[file][2]);
      h_trackVsJEC_1z[file][2]->Print("All");
    }
  }//end file loop
  out->Write();

  TCanvas * c1 = new TCanvas("c1","c1",600,600);
  meanPP[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  meanPP[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(dojtpt) meanPP[0]->GetXaxis()->SetTitle("p_{T}^{reco}");
  meanPP[0]->GetXaxis()->SetRangeUser(60,200);
  if(doMoreBins) meanPP[0]->GetXaxis()->SetRangeUser(60,440);
  if(doFineBins && showMore) meanPP[0]->GetXaxis()->SetRangeUser(50,300);
  meanPP[0]->Draw();
  meanPP[1]->SetMarkerColor(kRed);
  meanPP[1]->SetMarkerStyle(24);
  meanPP[1]->SetLineColor(kRed);
  meanPP[1]->Draw("same p");
  meanPP[2]->SetMarkerColor(kBlue);
  meanPP[2]->SetMarkerStyle(24);
  meanPP[2]->SetLineColor(kBlue);
  meanPP[2]->Draw("same p");
  meanPP[3]->SetMarkerColor(kGreen);
  meanPP[3]->SetLineColor(kGreen);
  meanPP[3]->SetMarkerStyle(24);
  meanPP[3]->Draw("same p");
  TLegend * leg1 = new TLegend(0.5,0.7,0.8,0.9);
  leg1->AddEntry(meanPP[0],"pp","p");
  leg1->AddEntry(meanPP[1],"pPb","p");
  leg1->AddEntry(meanPP[2],"Pbp","p");
  leg1->AddEntry(meanPP[3],"Smeared pp","p");
  leg1->SetBorderSize(0);
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/mean_moreBins.png"); 
  resoPP[0]->GetYaxis()->SetTitle("#sigma_{JER}");
  resoPP[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(dojtpt) resoPP[0]->GetXaxis()->SetTitle("p_{T}^{reco}");
  resoPP[0]->GetXaxis()->SetRangeUser(60,200);
  if(doFineBins && showMore) resoPP[0]->GetXaxis()->SetRangeUser(50,300);
  if(doMoreBins) resoPP[0]->GetXaxis()->SetRangeUser(60,440);
  resoPP[0]->Draw();
  resoPP[1]->SetMarkerColor(kRed);
  resoPP[1]->SetMarkerStyle(24);
  resoPP[1]->SetLineColor(kRed);
  resoPP[1]->Draw("same p");
  resoPP[2]->SetMarkerColor(kBlue);
  resoPP[2]->SetMarkerStyle(24);
  resoPP[2]->SetLineColor(kBlue);
  resoPP[2]->Draw("same p");
  resoPP[3]->SetMarkerColor(kGreen);
  resoPP[3]->SetLineColor(kGreen);
  resoPP[3]->SetMarkerStyle(24);
  resoPP[3]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/reso_moreBins.png"); 
  meanPP_hard[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  meanPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(dojtpt) meanPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{reco}");
  if(doMoreBins) meanPP_hard[0]->GetXaxis()->SetRangeUser(60,440);
  if(doFineBins && showMore) meanPP_hard[0]->GetXaxis()->SetRangeUser(50,300);
  meanPP_hard[0]->Draw();
  meanPP_hard[1]->SetMarkerColor(kRed);
  meanPP_hard[1]->Draw("same p");
  meanPP_hard[2]->SetMarkerColor(kBlue);
  meanPP_hard[2]->Draw("same p");
  meanPP_hard[3]->SetMarkerColor(kBlack);
  meanPP_hard[3]->SetMarkerStyle(24);
  meanPP_hard[3]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/mean_hard_0p5_moreBins.png"); 
  resoPP_hard[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  resoPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(dojtpt) resoPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{reco}");
  if(doMoreBins) resoPP_hard[0]->GetXaxis()->SetRangeUser(60,440);
  if(doFineBins && showMore) resoPP_hard[0]->GetXaxis()->SetRangeUser(50,300);
  resoPP_hard[0]->Draw();
  resoPP_hard[1]->SetMarkerColor(kRed);
  resoPP_hard[1]->Draw("same p");
  resoPP_hard[2]->SetMarkerColor(kBlue);
  resoPP_hard[2]->Draw("same p");
  resoPP_hard[3]->SetMarkerColor(kBlack);
  resoPP_hard[3]->SetMarkerStyle(24);
  resoPP_hard[3]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/reso_hard_0p5_moreBins.png"); 

  TF1 * ppFit = new TF1("ppFit","[0]*TMath::Power(x,-[1])",50,250);
  ppFit->SetParameters(0.15,0);
  resoPP[0]->Fit("ppFit","EMRI");  
  TF1 * pPbFit = new TF1("pPbFit","[0]*TMath::Power(x,-[1])",50,250);
  pPbFit->SetParameters(0.15,0);
  TH1D * resoPPForFit2 = (TH1D*)resoPP[1]->Clone("resoPPforFit2");
  resoPPForFit2->Fit("pPbFit","EMRI same");  
  TF1 * PbpFit = new TF1("PbpFit","[0]*TMath::Power(x,-[1])",50,250);
  PbpFit->SetParameters(0.15,0);
  TH1D * resoPPForFit = (TH1D*)resoPP[2]->Clone("resoPPforFit");
  resoPPForFit->Fit("PbpFit","EMRI same");  

  TH1D * addSmear = (TH1D*)resoPPForFit->Clone("additionalSmearing"); 
  addSmear->Reset();
  if(doFineBins){
  for(int r = 1; r<60; r++){
    float additionalSmearingPP;
    additionalSmearingPP = TMath::Power(TMath::Abs(TMath::Power((resoPP[1]->GetBinContent(r)*0.71+resoPP[2]->GetBinContent(r)*1.29)/2.0,2)-TMath::Power(resoPP[0]->GetBinContent(r),2)),0.5);
    if(!((resoPP[1]->GetBinContent(r)*0.71+resoPP[2]->GetBinContent(r)*1.29)/2.0-resoPP[0]->GetBinContent(r) > 0)) additionalSmearingPP = -additionalSmearingPP;
    addSmear->SetBinContent(r,additionalSmearingPP);
    addSmear->SetBinError(r,0);
    std::cout << "Additional smearing for bin " << r << ": " << additionalSmearingPP <<std::endl;
  }
  
  float previousBin[65] = {0};
  for(int r = 1; r<60; r++){
    if(r==1){
      previousBin[r] = addSmear->GetBinContent(r);
      addSmear->SetBinContent(r,(2*previousBin[r]+addSmear->GetBinContent(r+1))/3.0);
    }else if(r>1 && r<59){
      previousBin[r] = addSmear->GetBinContent(r);
      addSmear->SetBinContent(r,(previousBin[r-1]+previousBin[r]+addSmear->GetBinContent(r+1))/3.0);
    }else{
      addSmear->SetBinContent(r,(previousBin[r-1]+addSmear->GetBinContent(r)*2)/3.0);
    } 
  }
  }
  addSmear->SetDirectory(out);
  out->Write();

  c1->SetLogx(0);
  c1->SetLogy(0);
  TLegend * leg3 = new TLegend(0.3,0.15,0.6,0.4);
  leg3->SetBorderSize(0);
  for(int file=0; file<3; file++){
    h_trackVsJEC_1z[file][0]->Draw();
    h_trackVsJEC_1z[file][0]->GetXaxis()->SetTitle("Leading ch. particle Z (p_{T}^{trk}/p_{T}^{gen jet})");
    h_trackVsJEC_1z[file][0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
    h_trackVsJEC_1z[file][1]->SetMarkerColor(kBlue);
    h_trackVsJEC_1z[file][1]->SetLineColor(kBlue);
    h_trackVsJEC_1z[file][1]->Draw("same");
    h_trackVsJEC_1z[file][2]->SetMarkerColor(kRed);
    h_trackVsJEC_1z[file][2]->SetLineColor(kRed);
    h_trackVsJEC_1z[file][2]->Draw("same");

    if(file==0){
      leg3->AddEntry(h_trackVsJEC_1z[file][0],"60 < p_{T}^{jet} < 100 GeV","p");
      leg3->AddEntry(h_trackVsJEC_1z[file][1],"100 < p_{T}^{jet} < 140 GeV","p");
      leg3->AddEntry(h_trackVsJEC_1z[file][2],"140 < p_{T}^{jet} < 200 GeV","p");
    }
    leg3->Draw("same");
    c1->SaveAs(Form("checkJESJERPlots/JEC_vs_Z_%d.png",file)); 
    c1->SaveAs(Form("checkJESJERPlots/JEC_vs_Z_%d.pdf",file)); 
  } 

  /*delete *L2JES; 
  delete *L2JER; 
  delete *L3JES; 
  delete *MCTruth;*/
}
