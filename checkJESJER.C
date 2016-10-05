#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TF1.h"
#include "TLegend.h"

//#include "residualJEC.h"
//#include "jetSmearing.h"
//#include "getJEC_2nd.h"
//#include "getJEC_1st.h"
//#include "getJEC_L2L3res.h"
//#include "getJEC_SystError.h"

#include "L2L3ResidualWFits.h"
#include "MCTruthResidual.h"

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

void checkJESJER(bool doMoreBins = false){
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  

  TH1D * res[9][3];
  TH1D * res_hard[9][3];
  TH1D * meanPP[3];
  TH1D * resoPP[3];
  TH1D * meanPP_hard[3];
  TH1D * resoPP_hard[3];

  int pthats[6] = {50,80,120,170,220,280};
  double crossSection5[6]  = {3.778E-03,4.412E-04,6.147E-05,1.018E-05,2.477E-06,6.160E-07};

  TFile * out = TFile::Open("checkJESJER.root","recreate");
  for(int file = 0; file<3; file++){
  std::cout << file << std::endl;
  for(int file2= 0; file2<6; file2++){
  std::cout << file2 << std::endl;
//for pp
  int radius = 3;
  double etacutcorr = 3;
  TString mode;
  if(file==0) mode = "pp5";
  if(file==1) mode = "pPb5";
  if(file==2) mode = "Pbp5";
  std::cout << mode << std::endl;
  L2ResidualJES * L2JES = new L2ResidualJES(radius,((int)etacutcorr),mode);
  L2ResidualJER * L2JER = new L2ResidualJER(mode);
  L3ResidualJES * L3JES = new L3ResidualJES(mode);
  MCTruthResidual * MCTruth = new MCTruthResidual(mode);
 
  TFile * f;
  if(file==0) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/abaty/mergedForests/Pythia6_Dijet%d_pp5TeV/0.root",pthats[file2]),"read");
//for pPb
  if(file==1) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/pPb/HP04/prod16/Hijing_Pythia_pt%d/HiForest_v77_merged01/pt%d_HP04_prod16_v77_merged_forest_0.root",pthats[file2],pthats[file2]),"read");
  if(file==2 && file2<4) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/Pbp/HP05/prod24/Hijing_Pythia_pt%d/HiForest_v84_merged02/pt%d_HP05_prod24_v84_merged_forest_0.root",pthats[file2],pthats[file2]),"read");
  //if(file==1) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForestpPb/MC/HiForest_pPbDir_QCD%d_pythia6_HiWinter13.root",pthats[file2]),"read");
  //if(file==2) f = TFile::Open(Form("/mnt/hadoop/cms/store/user/dgulhan/mergedForestpPb/MC/HiForest_PbpDir_QCD%d_pythia6_HiWinter13.root",pthats[file2]),"read");
  if(file==2 && file2>=4) continue;
  TTree * t = (TTree*)f->Get("ak3PFJetAnalyzer/t");
  TTree * hi;
  if(file!=0) hi = (TTree*)f->Get("hiEvtAnalyzer/HiTree");

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
  float pthat;

  float hiHF4;

  t->SetBranchAddress("nref",&nref);
  t->SetBranchAddress("ngen",&ngen);
  t->SetBranchAddress("pthat",&pthat);
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

  if(file==1) hi->SetBranchAddress("hiHFminusEta4",&hiHF4);
  if(file==2) hi->SetBranchAddress("hiHFplusEta4",&hiHF4);
  if(file==1 || file==2) t->AddFriend(hi);
 
  for(int i = 0; i< (doMoreBins?9:5) ; i++){
    std::cout << i << std::endl;
    res[i][file] = new TH1D(Form("res%d%d",i,file),"",30,0.6,1.4);
    res_hard[i][file] = new TH1D(Form("res_hard%d%d",i,file),"",30,0.6,1.4);
  }
 
  for(int i = 0; i<t->GetEntries(); i++){
    t->GetEntry(i);
    if(i%1000==0) std::cout << i << "/" << t->GetEntries() << std::endl;
    for(int j = 0; j<nref; j++){
      float boost = 0;
      if(pthat > 2*pthats[file2]) continue;
      if(file==1) boost = 0.4654094531;
      if(file==2) boost = -0.4654094531;
      if(TMath::Abs(jteta[j]+boost)>1.5) continue;
      if(jtpt[j]>140 && jtpt[j]<200 && refpt[j]<=0 && rawpt[j]>30) std::cout << jtpt[j] << " " << refpt[j] << " " << chargedMax[j]/rawpt[j] << " " << chargedSum[j]/rawpt[j] << std::endl;
      if(rawpt[j]<30) continue;
      if(chargedSum[j]/rawpt[j]<0.05 || chargedSum[j]/rawpt[j]>0.95) continue;

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
  if(mode == "pPb5" || mode == "Pbp5"){
    correctedPt = MCTruth->getJEC_1st(rawptT,correctedPt,jtetaT); 
  }
  correctedPt = MCTruth->getResidualCorr(correctedPt,jtetaT);
  if(i%50000==0 && j==0) cout << "after MC correction jtpt = " << correctedPt << endl;
  if(!doMC){
   correctedPt = L2JES->getCorrectedPt(correctedPt, jtetaT);
   correctedPt = L3JES->getCorrectedPt(correctedPt);
  if(i%50000==0 && j==0) cout << "after data driven correction jtpt = " << correctedPt << endl;
  }else{   
   //double correctedPtOld = correctedPt;
   correctedPt = L2JER->getSmearedPt(correctedPt, jtetaT);
   //cout << "smearing factor:" << correctedPt/correctedPtOld << endl;
  //cout << "after additional smearing for MC jtpt = " << correctedPt << endl;
  }

  jtpt[j] = correctedPt;
  //std::cout << correctedPt << std::endl;
      /*
      if(file==0){
        jtpt[j] = 1.008*jtpt[j];
      }
      if(file==1){
        jtpt[j] = getJEC_1st("pPb5",rawpt[j],jtpt[j],jteta[j]); 
        jtpt[j] = getJEC_2nd(jtpt[j],jteta[j],"pPb5");
        jtpt[j] = getCorrectedJetPt("pPb5",true,jtpt[j],jteta[j]);
      }
      if(file==2){
        jtpt[j] = getJEC_1st("Pbp5",rawpt[j],jtpt[j],jteta[j]); 
        jtpt[j] = getJEC_2nd(jtpt[j],jteta[j],"Pbp5");
        jtpt[j] = getCorrectedJetPt("Pbp5",true,jtpt[j],jteta[j]);
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
      if(bin == -1) continue;
      //std::cout << "test" << std::endl;
      //std::cout << jtpt[j]/refpt[j] << " " << crossSection5[file2]/t->GetEntries() << std::endl;
      //std::cout << bin << " " << file << std::endl;
      res[bin][file]->Fill(jtpt[j]/refpt[j],crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4));
      //std::cout << "test" << std::endl;

      if(chargedMax[j]/rawpt[j]<0.5) continue;
      else res_hard[bin][file]->Fill(jtpt[j]/refpt[j],crossSection5[file2]/t->GetEntries()*MCTruth->getHFPbWeight(hiHF4)); 
      //std::cout << "test" << std::endl;
    }
  }
  res[0][file]->Print("All");
  }//end inner file loop
  
  const int jtptBin = 5;
  double jtptBins[jtptBin+1] = {60,80,100,120,140,200};
  const int jtptBinMore = 9;
  double jtptBinsMore[jtptBinMore+1] = {60,80,100,120,140,200,260,320,380,440};
  meanPP[file] = new TH1D(Form("meanPP%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
  resoPP[file] = new TH1D(Form("resoPP%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
  meanPP_hard[file] = new TH1D(Form("meanPP_hard%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
  resoPP_hard[file] = new TH1D(Form("resoPP_hard%d",file),"",doMoreBins?jtptBinMore:jtptBin,doMoreBins?jtptBinsMore:jtptBins);
  for(int j = 0; j<(doMoreBins?9:5); j++){
    float mean=0, meanErr=0, reso=0, resoErr=0;      
    FitGauss(res[j][file],0,mean,meanErr,reso,resoErr);
    std::cout << mean << " " << meanErr << " " << reso << " " << resoErr << std::endl;
    meanPP[file]->SetBinContent(j+1,mean);
    meanPP[file]->SetBinError(j+1,meanErr);
    resoPP[file]->SetBinContent(j+1,reso);
    resoPP[file]->SetBinError(j+1,resoErr);
  }
  for(int j = 0; j<(doMoreBins?9:5); j++){
    float mean=0, meanErr=0, reso=0, resoErr=0;      
    FitGauss(res_hard[j][file],0,mean,meanErr,reso,resoErr);
    std::cout << mean << " " << meanErr << " " << reso << " " << resoErr << std::endl;
    meanPP_hard[file]->SetBinContent(j+1,mean);
    meanPP_hard[file]->SetBinError(j+1,meanErr);
    resoPP_hard[file]->SetBinContent(j+1,reso);
    resoPP_hard[file]->SetBinError(j+1,resoErr);
  }
  meanPP[file]->GetYaxis()->SetRangeUser(0.95,1.05);
  resoPP[file]->GetYaxis()->SetRangeUser(-0.01,0.2);
  meanPP_hard[file]->GetYaxis()->SetRangeUser(0.95,1.15);
  resoPP_hard[file]->GetYaxis()->SetRangeUser(-0.01,0.2);
  }//end file loop
  out->Write();

  TCanvas * c1 = new TCanvas("c1","c1",600,600);
  meanPP[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  meanPP[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  meanPP[0]->GetXaxis()->SetRangeUser(60,200);
  if(doMoreBins) meanPP[0]->GetXaxis()->SetRangeUser(60,440);
  meanPP[0]->Draw();
  meanPP[1]->SetMarkerColor(kRed);
  meanPP[1]->Draw("same p");
  meanPP[2]->SetMarkerColor(kBlue);
  meanPP[2]->Draw("same p");
  TLegend * leg1 = new TLegend(0.6,0.7,0.9,0.9);
  leg1->AddEntry(meanPP[0],"pp","p");
  leg1->AddEntry(meanPP[1],"pPb","p");
  leg1->AddEntry(meanPP[2],"Pbp","p");
  leg1->SetBorderSize(0);
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/mean_moreBins.png"); 
  resoPP[0]->GetYaxis()->SetTitle("#sigma_{JER}");
  resoPP[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  resoPP[0]->GetXaxis()->SetRangeUser(60,200);
  if(doMoreBins) resoPP[0]->GetXaxis()->SetRangeUser(60,440);
  resoPP[0]->Draw();
  resoPP[1]->SetMarkerColor(kRed);
  resoPP[1]->Draw("same p");
  resoPP[2]->SetMarkerColor(kBlue);
  resoPP[2]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/reso_moreBins.png"); 
  meanPP_hard[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  meanPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(doMoreBins) meanPP_hard[0]->GetXaxis()->SetRangeUser(60,440);
  meanPP_hard[0]->Draw();
  meanPP_hard[1]->SetMarkerColor(kRed);
  meanPP_hard[1]->Draw("same p");
  meanPP_hard[2]->SetMarkerColor(kBlue);
  meanPP_hard[2]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/mean_hard_0p5_moreBins.png"); 
  resoPP_hard[0]->GetYaxis()->SetTitle("p_{T}^{reco}/p_{T}^{gen}");
  resoPP_hard[0]->GetXaxis()->SetTitle("p_{T}^{gen}");
  if(doMoreBins) resoPP_hard[0]->GetXaxis()->SetRangeUser(60,440);
  resoPP_hard[0]->Draw();
  resoPP_hard[1]->SetMarkerColor(kRed);
  resoPP_hard[1]->Draw("same p");
  resoPP_hard[2]->SetMarkerColor(kBlue);
  resoPP_hard[2]->Draw("same p");
  leg1->Draw("same");
  c1->SaveAs("checkJESJERPlots/reso_hard_0p5_moreBins.png"); 

  /*delete *L2JES; 
  delete *L2JER; 
  delete *L3JES; 
  delete *MCTruth;*/
}
