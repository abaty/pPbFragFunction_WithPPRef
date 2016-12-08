#include "TH2D.h"
#include "TH1D.h"
#include "TMath.h"
#include "TRandom.h"
#include "string.h"
#include "TDatime.h"
#include <iostream>
#include <fstream>
#include <istream>
#include <cmath>
#include "factorizedPtCorr.h"
#include "SpectraFiles.h"
#include "jetSmearing.h"
/*#include "residualJEC.h"
#include "getJEC_2nd.h"
#include "getJEC_1st.h"
#include "getJEC_L2L3res.h"*/
#include "MCTruthResidual.h"
#include "L2L3ResidualWFits.h"
#include "getJEC_SystError.h"
#include "getTrkCorr.h"

//CAREFUL WITH TRACKING CORRECTIONS >200 GeV?


const double pPbRapidity = 0.4654094531;
//const double pPbRapidity = 0;
const int nJetBins = 400;
const int ntrackBins=25;
const double axis[ntrackBins] = {0.5, 0.63, 0.77,  1.03,1.38, 1.84, 2.46, 3.29,  4.40, 5.88,  7.87,  10.52, 14.06,  18.8, 25.13,  33.58,  44.89,  60, 80, 100, 120, 140, 200,300,400};


//calculates dr^2 to avoid the slow TMath() Sqrt function
double getdR2(double jet_eta, double jet_phi, double track_eta, double track_phi)
{
  return TMath::Power(jet_eta-track_eta,2)+TMath::Power(acos(cos(jet_phi-track_phi)),2);
}

//calculation of Xi
double getXi(double jet_pt, double jet_eta, double jet_phi, double track_pt, double track_eta, double track_phi)
{
  double xi = -2;
  xi = TMath::Log((jet_pt/track_pt)*TMath::Power(TMath::CosH(jet_eta),2)/(TMath::Cos(track_phi-jet_phi) + TMath::SinH(jet_eta)*TMath::SinH(track_eta)));
  return xi;
}

//modes are pp2,pp7,pPb5,Pbp5,pp5
void Spectra(const char* inputJets, const char* inputMB, const char* mode = "pp2", const char* trigger = "jet80", bool isMC = 0,int  jobNum = -1, int typeUE = 3, double jetEtaMin = 0, double jetEtaMax = 1.5)
{
  if(strcmp(mode,"pp2") && strcmp(mode,"pp7") && strcmp(mode,"pPb5") && strcmp(mode,"Pbp5") && strcmp(mode,"pp5") && strcmp(mode,"ppref5"))
  {
    std::cout << "Invalid mode, modes are pp2, pp7, pPb5, Pbp5" << std::endl;
    return;
  }

  bool doMidRapidity = false;
  double coneSize = 0.3;

  std::string trkSkimVars;
  trkSkimVars=   "trkPt:trkEta:trkPhi:jtpt:jteta:jtphi:rawpt:chargedSum:corr:weight";
  TNtuple * highZskim  = new TNtuple("highZskim","",trkSkimVars.data()); 

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2(); 
  TDatime * dateTime = new TDatime();
  TRandom * rand = new TRandom(dateTime->GetTime());
  
  //Tracking initialization
  bool ispPbStyleCorr = true;
  TrkCorrObj* trkCorrObj;
  if(!ispPbStyleCorr && strcmp(mode,"ppref5")==0) trkCorrObj = new TrkCorrObj("TrkCorr_July22_Iterative_pp_eta2p4/");
  //if(ispPbStyleCorr && strcmp(mode,"ppref5")==0) trkCorrObj = new TrkCorrObj("TrkCorr_July22_Iterative_pp_pPbFFCuts/");
  if(ispPbStyleCorr && strcmp(mode,"ppref5")==0) trkCorrObj = new TrkCorrObj("TrkCorr_Aug19_Iterative_pp_pPbFFCuts/");
  const sampleType sType = kPPDATA;
  InitCorrFiles(sType,mode);
  InitCorrHists(sType);

  //JEC initialization
  TString Jtmode;
  if(strcmp(mode,"ppref5")==0) Jtmode = "pp5";
  if(((strcmp(mode,"pPb5")==0) || (strcmp(mode,"pp5")==0)) && isMC==0) Jtmode = "Pbp5";
  if(strcmp(mode,"Pbp5")==0 && isMC==0) Jtmode = "pPb5";
  if(((strcmp(mode,"pPb5")==0) || (strcmp(mode,"pp5")==0)) && isMC==1) Jtmode = "pPb5";
  if(strcmp(mode,"Pbp5")==0 && isMC==1) Jtmode = "Pbp5";
  L2ResidualJES * L2JES = new L2ResidualJES(3,3,Jtmode);
  L2ResidualJER * L2JER = new L2ResidualJER(Jtmode);
  L3ResidualJES * L3JES = new L3ResidualJES(Jtmode);
  MCTruthResidual * MCTruth = new MCTruthResidual(Jtmode);

  //tracking variables calc
  const int ny=20;
  double trkx[ny+1];
  double inix=log(0.5)/log(10);
  double delta=(log(200)-log(0.5))/(ny*log(10));
  for(int ix=0; ix<ny+1;ix++)
  {
    trkx[ix]=pow(10,inix);
    inix+=delta;
  }
  TH1D * totalGenTGenJ_Tracking = new TH1D("totalGenTGenJ_Tracking","",ny,trkx);
  TH1D * totalRecoTGenJ_Tracking = new TH1D("totalRecoTGenJ_Tracking","",ny,trkx);
  TH1D * totalGenTRecoJ_Tracking = new TH1D("totalGenTRecoJ_Tracking","",ny,trkx);
  TH1D * totalRecoTRecoJ_Tracking = new TH1D("totalRecoTRecoJ_Tracking","",ny,trkx);
  TH1D * totalGenTRecoJ_Tracking_eta = new TH1D("totalGenTRecoJ_Tracking_eta","",10,-2.4,2.4);
  TH1D * totalRecoTRecoJ_Tracking_eta = new TH1D("totalRecoTRecoJ_Tracking_eta","",10,-2.4,2.4);
  TH1D * totalGenTRecoJ_Tracking_largeEta = new TH1D("totalGenTRecoJ_Tracking_largeEta","",ny,trkx);
  TH1D * totalRecoTRecoJ_Tracking_largeEta = new TH1D("totalRecoTRecoJ_Tracking_largeEta","",ny,trkx);
  TH1D * totalGenTRecoJ_Tracking_lowEta = new TH1D("totalGenTRecoJ_Tracking_lowEta","",ny,trkx);
  TH1D * totalRecoTRecoJ_Tracking_lowEta = new TH1D("totalRecoTRecoJ_Tracking_lowEta","",ny,trkx);
  TH1D * totalGenTRecoJ_Tracking_largePt = new TH1D("totalGenTRecoJ_Tracking_largePt","",10,-2.4,2.4);
  TH1D * totalRecoTRecoJ_Tracking_largePt = new TH1D("totalRecoTRecoJ_Tracking_largePt","",10,-2.4,2.4);

  //for testing code in interactive mode only
  if(jobNum == -1)
  {
    std::cout << "here" << std::endl;
    getInputFile("/mnt/hadoop/cms/store/user/abaty/FF_forests/skims/pPb5/data/pPb5jet80_0_20150227_0.root",0);
    //getInputFile("/mnt/hadoop/cms/store/user/abaty/FF_forests/skims/ppref5/data/ppref5jet80_0_20160712_6.root",0);
    if(typeUE==2){ 
      getInputFileMix("/mnt/hadoop/cms/store/user/abaty/FF_forests/skims/pPb5/data/pPb5MB_0_20150411_205.root",0);
      //getInputFileMix("/mnt/hadoop/cms/store/user/abaty/FF_forests/skims/ppref5/data/ppref5MB_0_20160718_53.root",0);
    }
  }
  else{
    getInputFile(inputJets,isMC);
    if(typeUE==2) getInputFileMix(inputMB,isMC);
  }

  //different nonzero variations are used for systematics checks, variation 0 is for the basic calculation
  for(int v = 0; v<variations; v++)
  {
    setJetPtRange(mode,trigger,(int)(v==29),(strcmp(mode,"ppref5")==0 &&  isMC)?1:0);
  
    if((strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp5")==0 || strcmp(mode,"ppref5")==0) && !(v==0 || v==1 || v==2 || v==5 || v==6 || v==9 || v==13 || v==26 || v==30 || v==31 || v==34 || v==35 || v==36 || v==37 || v==38)) continue;
    if(typeUE!=0 && v==26) continue;
    if(typeUE==1 && v!=0) continue;
    //if(typeUE==2 && (v!=0 && v!=31)) continue;

    //reco
    h_jet = new TH1D("h_jet","",nJetBins,0,1000); 
    h_track = new TH2D("h_track","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE = new TH2D("h_trackUE","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi = new TH2D("h_track_xi","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi = new TH2D("h_trackUE_xi","",nJetBins,0,1000,24,0,6.0);
  
    //gen/reco combinations  
    h_track_rJgT = new TH2D("h_track_rJgT","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_rJgT = new TH2D("h_trackUE_rJgT","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_rJgT = new TH2D("h_track_xi_rJgT","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_rJgT = new TH2D("h_trackUE_xi_rJgT","",nJetBins,0,1000,24,0,6.0);
  
    h_track_gJrT = new TH2D("h_track_gJrT","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_gJrT = new TH2D("h_trackUE_gJrT","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_gJrT = new TH2D("h_track_xi_gJrT","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_gJrT = new TH2D("h_trackUE_xi_gJrT","",nJetBins,0,1000,24,0,6.0);
  
    //gen
    h_jet_gen = new TH1D("h_jet_gen","",nJetBins,0,1000);
    h_track_gen = new TH2D("h_track_gen","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_gen = new TH2D("h_trackUE_gen","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_gen = new TH2D("h_track_xi_gen","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_gen = new TH2D("h_trackUE_xi_gen","",nJetBins,0,1000,24,0,6.0);   
  
    //Quark Gluons checks 
    h_jet_Q = new TH1D("h_jet_Q","",nJetBins,0,1000);
    h_track_Q = new TH2D("h_track_Q","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_Q = new TH2D("h_trackUE_Q","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_Q = new TH2D("h_track_xi_Q","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_Q = new TH2D("h_trackUE_xi_Q","",nJetBins,0,1000,24,0,6.0);
    h_jet_G = new TH1D("h_jet_G","",nJetBins,0,1000);
    h_track_G = new TH2D("h_track_G","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_G = new TH2D("h_trackUE_G","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_G = new TH2D("h_track_xi_G","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_G = new TH2D("h_trackUE_xi_G","",nJetBins,0,1000,24,0,6.0);
    h_jet_gen_Q = new TH1D("h_jet_gen_Q","",nJetBins,0,1000);
    h_track_gen_Q = new TH2D("h_track_gen_Q","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_gen_Q = new TH2D("h_trackUE_gen_Q","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_gen_Q = new TH2D("h_track_xi_gen_Q","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_gen_Q = new TH2D("h_trackUE_xi_gen_Q","",nJetBins,0,1000,24,0,6.0);
    h_jet_gen_G = new TH1D("h_jet_gen_G","",nJetBins,0,1000);
    h_track_gen_G = new TH2D("h_track_gen_G","",nJetBins,0,1000,ntrackBins-1,axis);
    h_trackUE_gen_G = new TH2D("h_trackUE_gen_G","",nJetBins,0,1000,ntrackBins-1,axis);
    h_track_xi_gen_G = new TH2D("h_track_xi_gen_G","",nJetBins,0,1000,24,0,6.0);
    h_trackUE_xi_gen_G = new TH2D("h_trackUE_xi_gen_G","",nJetBins,0,1000,24,0,6.0); 

    const int trkVsJEC_jetBins = 9;
    const double trkVsJEC_jetBinsAxis[trkVsJEC_jetBins+1] = {50,60,80,100,120,140,160,180,200,220};
    h_trackVsJEC = new TH2D("h_trackVsJEC","",ntrackBins-1,axis,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_trackVsJEC_weights = new TH2D("h_trackVsJEC_weights","",ntrackBins-1,axis,trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_JEC = new TProfile("h_JEC","",trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    h_JEC_weights = new TProfile("h_JEC_weights","",trkVsJEC_jetBins,trkVsJEC_jetBinsAxis);
    TProfile * h_L2ResidualCorr = new TProfile("h_L2ResidualCorr","",160,50,210);
    TProfile * h_L2ResidualCorr_eta = new TProfile("h_L2ResidualCorr_eta","",40,-2,2);

    //trigger combination (for xchecks)
    TH1D * nLeadJet40 = new TH1D("nLeadJet40","nLeadJet40",100,0,1000);
    TH1D * nLeadJet60 = new TH1D("nLeadJet60","nLeadJet60",100,0,1000);
    TH1D * nLeadJet80 = new TH1D("nLeadJet80","nLeadJet80",100,0,1000);
  
    //boosting
    double boost = 0;
    if(strcmp(mode,"pPb5") == 0 || strcmp(mode,"pp5") == 0) boost = pPbRapidity;
    else if(strcmp(mode,"Pbp5") == 0) boost = -pPbRapidity;
    if(isMC) boost = -boost;//pPb MC has proton in + direction, pPb data has it in minus, and both are reversed for the 'Pbp definition'
    std::cout << mode << " mode is specified; using a boost of: " << boost << std::endl;
  
    //variables for mixing
    int startMixEvt = 0;
    int lastMixEvt = 1;
    //mixing variable for pPb systems
    //not updated for pp reference
    /*float mixHFProxyP[300000] = {0};
    float mixHFProxyM[300000] = {0};
    if(typeUE==2 && (strcmp(mode,"pPb5") == 0 || strcmp(mode,"Pbp5") == 0))
    {
      for(int i = 0; i<trackMix->GetEntries();i++)
      {
        if(i%10000==0) std::cout << i << "/" << trackMix->GetEntries() << " Calculating sideband multiplicities for MB matching..." << std::endl;
        trackMix->GetEntry(i);
        for(int t = 0; t<nTrkMix; t++)
        {
          if(strcmp(mode,"pPb5") == 0 || strcmp(mode,"pp5") == 0)
          {
            if(trkEtaMix[t]>2 && trkEtaMix[t]<2.4 && trkPtMix[t]>0.5 && trkPtMix[t]<3 && highPurityMix[t]) mixHFProxy[i]++;
          }
          else if(strcmp(mode,"Pbp5") == 0)
          {
            if(trkEtaMix[t]<-2 && trkEtaMix[t]>-2.4 && trkPtMix[t]>0.5 && trkPtMix[t]<3 && highPurityMix[t]) mixHFProxy[i]++;
          }
        }
      }
    }*/

    //if(strcmp(mode, "Pbp5")==0) startMixEvt = 6743253;
    if(typeUE==2) lastMixEvt = evtMix->GetEntries();
  
    int nEntry = track->GetEntries();
    //nEntry = 10;
    for(int i=0; i<nEntry; i++)
    {
      weight = 1;
      getInputEntry(i);
      //fix to pthat weights if using first production of ppref with Pythia8 (comment otherwise)
      /*if(isMC && strcmp(mode,"ppref5")==0){
        if(TMath::Abs(weight-5.334e-07)<1e-9) weight = weight*4.069E-03/5.335E-01;
        else if(TMath::Abs(weight-1.695e-08)<1e-10)  weight = weight*4.959E-04/3.378E-02;
        else if(TMath::Abs(weight-3.795e-09)<1e-11)  weight = weight*7.096E-5/3.778E-03;
        else if(TMath::Abs(weight-4.936e-10)<1e-12)  weight = weight*1.223E-05/4.412E-04;
      }*/
      //if(weight<1e-10) continue;//cutting out upper pthat samples just to keep gen similar

      if(v==31 && nVtx>1) continue;
      
      /*if(v==32 && nVtx<5) continue;
      if(v==33 && nVtx>=5 && nVtx<7) continue;*/
      
      //vz cut
      if(TMath::Abs(vz)>15) continue;
      totalEvts->Fill(1,weight);     
  
      if(i%10000 == 0) std::cout << i << "/" << nEntry << std::endl;
      //finding a MB event to mix with if needed 
      if(typeUE==2)
      {
        //calculating mutliplicity mixing variable for pPb
        /*float HFProxy = 0;
        if(strcmp(mode,"pPb5") == 0 || strcmp(mode,"pp5") == 0 || strcmp(mode,"Pbp5") == 0)
        {
          track->GetEntry(i);
          for(int t = 0; t<nTrk; t++)
          {
            if(strcmp(mode,"pPb5") == 0 || strcmp(mode,"pp5") == 0)
            {
              if(trkEta[t]>2 && trkEta[t]<2.4 && trkPt[t]>0.5 && trkPt[t]<3 && highPurity[t]) HFProxy++;
            }
            else if(strcmp(mode,"Pbp5") == 0)
            {
              if(trkEta[t]<-2 && trkEta[t]>-2.4 && trkPt[t]>0.5 && trkPt[t]<3 && highPurity[t]) HFProxy++;
            }
          }
        }*/
        
        int loopIter=0;
        int maxIter = trackMix->GetEntries(); 
        while(true)
        {
          //preventing infinite loop
          loopIter++; 
          if(loopIter == maxIter)
          {
            std::cout << "error finding matching MB event, using random MB event" << std::endl;
            break;
          }            
   
          //finding matching event
          lastMixEvt++;
          if(lastMixEvt>startMixEvt+maxIter) lastMixEvt = startMixEvt;
          evtMix->GetEntry(lastMixEvt);  
          if((strcmp(mode,"pPb5")==0 || strcmp(mode,"pp5")==0) && (TMath::Abs(hiHFplusEta4-hiHFplusEta4Mix)<2 || (hiHFplusEta4>20 && hiHFplusEta4Mix>18)) && (TMath::Abs(hiHFminusEta4-hiHFminusEta4Mix)<2 || (hiHFminusEta4>13 && hiHFminusEta4Mix>11))) break;//require the both HFs to be 'close' to each other, but open up the matching space for very high HF values
          else if(strcmp(mode,"Pbp5")==0 && (TMath::Abs(hiHFplusEta4-hiHFplusEta4Mix)<2 || (hiHFplusEta4>13 && hiHFplusEta4Mix>11)) && (TMath::Abs(hiHFminusEta4-hiHFminusEta4Mix)<2 || (hiHFminusEta4>20 && hiHFminusEta4Mix>18))) break;//same thing with sides swapped
          else if(strcmp(mode,"ppref5")==0) break;
        }
      }

      //pPb5 and Pbp5 activity weightings
      if((strcmp(mode,"pPb5")==0) && isMC) weight = weight*MCTruth->getHFPbWeight(hiHFminusEta4);
      if((strcmp(mode,"Pbp5")==0) && isMC) weight = weight*MCTruth->getHFPbWeight(hiHFplusEta4);
      //trigger combining weights (magic numbers are average prescales)
      if((strcmp(mode,"ppref5")==0) && (strcmp(trigger,"jet60")==0)  && !isMC) weight = weight*6.283512;
      if((strcmp(mode,"ppref5")==0) && (strcmp(trigger,"jet40")==0)  && !isMC) weight = weight*28.90433;
      if((strcmp(mode,"pPb5")==0) && (strcmp(trigger,"jet40")==0)  && !isMC) weight = weight*28.79457;
      if((strcmp(mode,"Pbp5")==0) && (strcmp(trigger,"jet40")==0)  && !isMC) weight = weight*28.79457;
 
      //starting jet loop Reco
      bool filledTriggerObject = false;
      for(int j=0; j<nref; j++)
      {
        totalJetsHist->Fill(1,weight);
        if(rawpt[j]<30) continue;
        totalJetsRawPtCut->Fill(1,weight);
        //if((chargedSum[j]/rawpt[j]>0.95 || chargedSum[j]/rawpt[j]<0.05) && v!=30) continue;
        if((chargedSum[j]/rawpt[j]<0.05) && v!=30) continue;
        totalJetsChargeSumCut->Fill(1,weight);
        if(TMath::Abs(jteta[j]+boost) < jetEtaMin || TMath::Abs(jteta[j]+boost) > (doMidRapidity?0.25:jetEtaMax)) continue;
        if(typeUE==1 && (TMath::Abs(jteta[j]+boost))<coneSize) continue;
        totalJetsEtaCutHist->Fill(1,weight);
        if(strcmp(trigger,"jet40")==0 && filledTriggerObject==false){ nLeadJet40->Fill(jtpt[j],1); filledTriggerObject=true;}
        if(strcmp(trigger,"jet60")==0 && filledTriggerObject==false){ nLeadJet60->Fill(jtpt[j],1); filledTriggerObject=true;}
        if(strcmp(trigger,"jet80")==0 && filledTriggerObject==false){ nLeadJet80->Fill(jtpt[j],1); filledTriggerObject=true;}
  
      //JEC and its corrections applied
        float appliedL2Corr = 0;
        if(!(strcmp(mode,"ppref5")==0)){
          jtpt[j] = MCTruth->getJEC_1st(rawpt[j],jtpt[j],jteta[j]); 
          jtpt[j] = MCTruth->getResidualCorr(jtpt[j],jteta[j]);
          if(!isMC){
            appliedL2Corr = L2JES->getCorrectedPt(jtpt[j], jteta[j])/jtpt[j]-1;
            jtpt[j] = L2JES->getCorrectedPt(jtpt[j], jteta[j]);
            jtpt[j] = L3JES->getCorrectedPt(jtpt[j]);// this is 1 right now
          }else{
            jtpt[j] = L2JER->getSmearedPt(jtpt[j], jteta[j]);
          } 
          //if(v==1 || v==3 || v==5)jtpt[j] = jtpt[j]+getJEC_SystError(mode,rawpt[j],jtpt[j],jteta[j],false);
          if(v==1) jtpt[j] = jtpt[j]*1.015;
          if(v==3 || v==5) jtpt[j] = jtpt[j]*1.025;
          //if(v==20 || v==22 || v==24)jtpt[j] = jtpt[j]+getJEC_SystError(mode,rawpt[j],jtpt[j],jteta[j],true);
          //if(v==2 || v==4 || v==6)jtpt[j] = jtpt[j]-getJEC_SystError(mode,rawpt[j],jtpt[j],jteta[j],false);
          if(v==2) jtpt[j] = jtpt[j]*0.985;
          if(v==4 || v==6)jtpt[j] = jtpt[j]*0.975;
          //if(v==21 || v==23 || v==25)jtpt[j] = jtpt[j]-getJEC_SystError(mode,rawpt[j],jtpt[j],jteta[j],true);
          if(v==7 || v==8 || v==9)jtpt[j] = getJERCorrected(mode,jtpt[j],0.05);
        }
        else{
          jtpt[j] = MCTruth->getResidualCorr(jtpt[j],jteta[j]);
          if(!isMC){
            appliedL2Corr = L2JES->getCorrectedPt(jtpt[j], jteta[j])/jtpt[j]-1;
            jtpt[j] = L2JES->getCorrectedPt(jtpt[j], jteta[j]);
            jtpt[j] = L3JES->getCorrectedPt(jtpt[j]);
          }else{
            //constant shift of 0.8% in MC (not in data because is taken care of by L3Residual)
            jtpt[j] = jtpt[j]*1.008;
            jtpt[j] = L2JER->getSmearedPt(jtpt[j], jteta[j]);
          }
        }
 
        //JEC thing here
        /*
        float leadingHadronPt = 0;
        for(int t=0; t<nTrk; t++){
          if((TMath::Abs(trkEta[t]-jteta[j])>coneSize) || (TMath::Abs(trkPhi[t]-jtphi[j])>coneSize)) continue;
          if(trkPt[t]<leadingHadronPt) continue;
          if(trkPt[t] <= 1 || !highPurity[t] || TMath::Abs(trkEta[t])>2.4 ) continue;
          //if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.1)) continue;            
          //else if(strcmp(mode,"ppref5")==0 && !ispPbStyleCorr &&  (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.3 || ((pfEcal[t]+pfHcal[t])/TMath::CosH(trkEta[t])<0.5*trkPt[t] && trkPt[t]>20))) continue;            
          if(getdR2(jteta[j]+boost,jtphi[j],trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize){
            leadingHadronPt = trkPt[t];
          }
        }
        if(isMC && refpt[j]>20){
          h_trackVsJEC->Fill(leadingHadronPt,refpt[j],jtpt[j]/refpt[j]*weight);
          h_trackVsJEC_weights->Fill(leadingHadronPt,refpt[j],weight);
        }
        if(strcmp(mode,"ppref5")==0) jtpt[j] = getFragJECFactor("ppref5",leadingHadronPt,jtpt[j]);
        if(strcmp(mode,"pPb5")==0) jtpt[j] = getFragJECFactor("pPb5",leadingHadronPt,jtpt[j]);
        if(strcmp(mode,"Pbp5")==0) jtpt[j] = getFragJECFactor("Pbp5",leadingHadronPt,jtpt[j]);*/
        //end FF correction
        

        if(strcmp(mode,"ppref5")==0){
          //smearing to match pPb resolution
          jtpt[j] = getJERCorrected("ppref5",jtpt[j],getPPDataSmearFactor(jtpt[j])); 

          if(v==34) jtpt[j] = jtpt[j]*1.025;
          if(v==37) jtpt[j] = jtpt[j]*1.015;
          if(v==35) jtpt[j] = jtpt[j]*0.975;
          if(v==38) jtpt[j] = jtpt[j]*0.985;
          if(v==36) jtpt[j] = getJERCorrected("pPb5",jtpt[j],0.05); 
        }
        //checking fake jet contribution
        //if(isMC && refpt[j]<0) continue;

        //adding 10 GeV buffer here to the lowest bin
        totalJetsPtCutHist->Fill(1,weight);    
        if(jtpt[j]<((lowJetPtBound==60)?50:lowJetPtBound) || jtpt[j]>=((upJetPtBound==1000)?1020:upJetPtBound)) continue;      
        /*if(isMC && refpt[j]>20){
          h_JEC->Fill(refpt[j],jtpt[j]/refpt[j],weight);
          h_JEC_weights->Fill(refpt[j],weight);
        }*/
  
      //quark or gluon only contributions for MC
        bool isQ = false;
        bool isG = false;
        if(isMC && TMath::Abs(refparton_flavor[j])<901 && TMath::Abs(refparton_flavor[j])!=21) isQ=true;
        if(isMC && TMath::Abs(refparton_flavor[j])==21) isG=true;
        if(isQ) h_jet_Q->Fill(jtpt[j],weight);
        if(isG) h_jet_G->Fill(jtpt[j],weight);
        h_jet->Fill(jtpt[j],weight);
   

        if(jtpt[j]<lowJetPtBound  || jtpt[j]>=upJetPtBound) continue;
        if(v==0) h_L2ResidualCorr->Fill(jtpt[j],appliedL2Corr);
        if(v==0) h_L2ResidualCorr_eta->Fill(jteta[j],appliedL2Corr);

        for(int t=0; t<nTrk; t++)
        {
          if((trkCharge[t]!=1 && v==27) || (trkCharge[t]!=-1 && v==28)) continue;     
          if(trkPt[t] <= 0.5 || trkPt[t] > 400 || !highPurity[t] || TMath::Abs(trkEta[t])>2.4 ) continue;
          if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.1)) continue;            
          else if(strcmp(mode,"ppref5")==0 && !ispPbStyleCorr && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.3 || ((pfEcal[t]+pfHcal[t])/TMath::CosH(trkEta[t])<0.5*trkPt[t] && trkPt[t]>20))) continue;            
          //calculating r_min for tracking correction
          double r_min = 9;
          double rmin_pt = 60;
          for(int j2 = 0; j2<nref; j2++)
          {
            if(TMath::Abs(jteta[j2])>2 || TMath::Abs(jtpt[j2]) < 60 || jtpt[j2]>1000) continue;
            double r_min_temp = TMath::Power(getdR2(jteta[j2],jtphi[j2],trkEta[t],trkPhi[t]),0.5);
            if(r_min_temp < r_min){ r_min = r_min_temp;
            rmin_pt = jtpt[j2];}
          }
          r_min = rmin_pt;
   
          //Filling track spectrum in jet cone
          if(getdR2(jteta[j]+boost,jtphi[j],trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
          {
            double trkCorr = 1;
            //std::cout << "before track" ;
            if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
            else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,jtpt[j]);     
            //std::cout << " after track" << std::endl;

            if(v==13) trkCorr=1; 
            if(std::isfinite(trkCorr))
            {  
              totalRecoTRecoJ_Tracking->Fill(trkPt[t],weight*trkCorr);
              totalRecoTRecoJ_Tracking_eta->Fill(trkEta[t],weight*trkCorr);
              if(TMath::Abs(trkEta[t])>1.9)totalRecoTRecoJ_Tracking_largeEta->Fill(trkPt[t],weight*trkCorr);
              if(TMath::Abs(trkEta[t])<1.9)totalRecoTRecoJ_Tracking_lowEta->Fill(trkPt[t],weight*trkCorr);
              if(trkPt[t]>10)totalRecoTRecoJ_Tracking_largePt->Fill(trkEta[t],weight*trkCorr);
              if(v==0 && typeUE==2 &&((jtpt[j]<140 && trkPt[t]>40 && (jtpt[j]-trkPt[t])<40) || trkPt[t]>140)){
                float skimEntry[] = {trkPt[t],trkEta[t],trkPhi[t],jtpt[j],jteta[j],jtphi[j],rawpt[j],chargedSum[j],(float)trkCorr,weight};
                highZskim->Fill(skimEntry);
              }

              h_track->Fill(jtpt[j],trkPt[t],trkCorr*weight);
              h_track_xi->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              if(isQ)
              {
                h_track_Q->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_track_xi_Q->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
              if(isG)
              {
                h_track_G->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_track_xi_G->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
            }
          }
       
          //Phi rotated UE subtraction
          //returns either +-1 to rotate clockwise or ccw randomly
          double rotationDirection = 2*(int)rand->Integer(2)-1;
          if(v==26) rotationDirection = rotationDirection*2.0/3.0;
  
          if(typeUE==0 && getdR2(jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
          {
            double trkCorr = 1;
            if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
            else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,jtpt[j]);     
            
            if(v==13) trkCorr=1;
            if(std::isfinite(trkCorr))
            {
              h_trackUE->Fill(jtpt[j],trkPt[t],trkCorr*weight);  
              h_trackUE_xi->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              if(isQ)
              {
                h_trackUE_Q->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_trackUE_xi_Q->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
              if(isG)
              {
                h_trackUE_G->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_trackUE_xi_G->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
            }
          }
  
          //Eta Reflected UE subtraction
          if(typeUE==1 && getdR2(-1*(jteta[j]+boost),jtphi[j],trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
          {
            double trkCorr = 1;
            if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
            else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,jtpt[j]);     
            
            if(v==13) trkCorr=1;
            if(std::isfinite(trkCorr))
            {
              h_trackUE->Fill(jtpt[j],trkPt[t],trkCorr*weight); 
              h_trackUE_xi->Fill(jtpt[j],getXi(jtpt[j],-1*(jteta[j]+boost),jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              if(isG)
              {
                h_trackUE_G->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_trackUE_xi_G->Fill(jtpt[j],getXi(jtpt[j],-1*(jteta[j]+boost),jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
              if(isQ)
              {
                h_trackUE_Q->Fill(jtpt[j],trkPt[t],trkCorr*weight);
                h_trackUE_xi_Q->Fill(jtpt[j],getXi(jtpt[j],-1*(jteta[j]+boost),jtphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              } 
            }
          }
        }
  
        //UE subtraction w/ MB mixing
        if(typeUE==2)
        {        
          bool noJet = false;
          while(!noJet){
            noJet = true;
            getInputEntryMix(lastMixEvt);
            //veto on jets in the MB event
            for(int j2 = 0; j2<nrefMix; j2++)
            {
              if(TMath::Abs(jtetaMix[j2])>2 || TMath::Abs(jtptMix[j2]) > 20) continue;
              double r_min_temp = TMath::Power(getdR2(jtetaMix[j2],jtphiMix[j2],jteta[j],jtphi[j]),0.5);
              if(r_min_temp < coneSize) noJet=false;
            }                     
            if(!noJet){
              lastMixEvt++;                                                     
              if(lastMixEvt>startMixEvt+(trackMix->GetEntries())) lastMixEvt = startMixEvt;
            }
          }

          for(int t = 0; t < nTrkMix; t++)
          {
            if((trkChargeMix[t]!=1 && v==27) || (trkChargeMix[t]!=-1 && v==28)) continue;     
            if(trkPtMix[t] <= 0.5 || trkPtMix[t] > 400 || !highPurityMix[t] || TMath::Abs(trkEtaMix[t])>2.4 ) continue;
            if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1Mix[t]/trkDxyError1Mix[t]) > 3 || TMath::Abs(trkDz1Mix[t]/trkDzError1Mix[t]) > 3 || trkPtErrorMix[t]/trkPtMix[t] > 0.1)) continue;            
            else if(strcmp(mode,"ppref5")==0 && !ispPbStyleCorr && (TMath::Abs(trkDxy1Mix[t]/trkDxyError1Mix[t]) > 3 || TMath::Abs(trkDz1Mix[t]/trkDzError1Mix[t]) > 3 || trkPtErrorMix[t]/trkPtMix[t] > 0.3 || ((pfEcalMix[t]+pfHcalMix[t])/TMath::CosH(trkEtaMix[t])<0.5*trkPtMix[t] && trkPtMix[t]>20))) continue;            
  
            //calculating r_min for tracking correction
            double r_min = 9;
            double rmin_pt = 60;
            for(int j2 = 0; j2<nrefMix; j2++)
            {
              if(TMath::Abs(jtetaMix[j2])>2 || TMath::Abs(jtptMix[j2]) < 60 || jtptMix[j2]>1000) continue;
              double r_min_temp = TMath::Power(getdR2(jtetaMix[j2],jtphiMix[j2],trkEtaMix[t],trkPhiMix[t]),0.5);
              if(r_min_temp < r_min){ r_min = r_min_temp;
              rmin_pt = jtptMix[j2];}
            }                                                                          
            r_min = rmin_pt;
              
            //Filling track spectrum in jet cone
            if(getdR2(jteta[j]+boost,jtphi[j],trkEtaMix[t]+boost,trkPhiMix[t]) < coneSize*coneSize)
            {
              double trkCorr = 1;
              if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPtMix[t], sType), 1, trkPtMix[t], trkPhiMix[t], trkEtaMix[t], r_min, sType);     
              else trkCorr = trkCorrObj->getTrkCorr(trkPtMix[t],trkEtaMix[t],trkPhiMix[t],1,r_min,60);//just use lowest jet pt in table (60 GeV, should make no difference)     
 
              if(v==13) trkCorr=1;
              if(std::isfinite(trkCorr))
              {
                h_trackUE->Fill(jtpt[j],trkPtMix[t],trkCorr*weight);
                h_trackUE_xi->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPtMix[t],trkEtaMix[t]+boost,trkPhiMix[t]),trkCorr*weight);
                if(isG)
                {
                  h_trackUE_G->Fill(jtpt[j],trkPtMix[t],trkCorr*weight);
                  h_trackUE_xi_G->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPtMix[t],trkEtaMix[t]+boost,trkPhiMix[t]),trkCorr*weight);
                }
                if(isQ)
                {
                  h_trackUE_Q->Fill(jtpt[j],trkPtMix[t],trkCorr*weight);
                  h_trackUE_xi_Q->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],trkPtMix[t],trkEtaMix[t]+boost,trkPhiMix[t]),trkCorr*weight);
                }
              }
            }
          }
        }
 
        //reco jet gen particle 
        if(isMC)
        {
          for(int t=0; t<nParticle; t++)
          { 
            if(pPt[t] <= 0.5 || pPt[t] > 400 || TMath::Abs(pEta[t])>2.4 ) continue;
            //Filling track spectrum in jet cone
            if(getdR2(jteta[j]+boost,jtphi[j],pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
            {
              totalGenTRecoJ_Tracking->Fill(pPt[t],weight);
              totalGenTRecoJ_Tracking_eta->Fill(pEta[t],weight);
              if(TMath::Abs(pEta[t])>1.9)totalGenTRecoJ_Tracking_largeEta->Fill(pPt[t],weight);
              if(TMath::Abs(pEta[t])<1.9)totalGenTRecoJ_Tracking_lowEta->Fill(pPt[t],weight);
              if(pPt[t]>10)totalGenTRecoJ_Tracking_largePt->Fill(pEta[t],weight);
              h_track_rJgT->Fill(jtpt[j],pPt[t],weight);
              h_track_xi_rJgT->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
   
              //Phi rotated UE subtraction
              //returns either +-1 to rotate clockwise or ccw randomly
              double rotationDirection = 2*(int)rand->Integer(2)-1;
              if(v==26) rotationDirection = rotationDirection*2.0/3.0;

              if(typeUE==0 && getdR2(jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
              {
                h_trackUE_rJgT->Fill(jtpt[j],pPt[t],weight);  
                h_trackUE_xi_rJgT->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j]+rotationDirection*TMath::PiOver2(),pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
  
              //Eta Reflected UE subtraction
              if(typeUE==1 && getdR2(-1*(jteta[j]+boost),jtphi[j],pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
              {
                h_trackUE_rJgT->Fill(jtpt[j],pPt[t],weight); 
                h_trackUE_xi_rJgT->Fill(jtpt[j],getXi(jtpt[j],-1*(jteta[j]+boost),jtphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
            }
          }
          //UE subtraction w/ MB mixing
          if(typeUE==2)
          {        
            getInputEntryMix(lastMixEvt);
            for(int t = 0; t < nParticleMix; t++)
            {
              if(pPtMix[t] <= 0.5 || pPtMix[t] > 400 || TMath::Abs(pEtaMix[t])>2.4 ) continue;
            
              //Filling track spectrum in jet cone
              if(getdR2(jteta[j]+boost,jtphi[j],pEtaMix[t]+boost,pPhiMix[t]) < coneSize*coneSize)
              {
                h_trackUE_rJgT->Fill(jtpt[j],pPtMix[t],weight);
                h_trackUE_xi_rJgT->Fill(jtpt[j],getXi(jtpt[j],jteta[j]+boost,jtphi[j],pPtMix[t],pEtaMix[t]+boost,pPhiMix[t]),weight);
              }
            }
          }
        }
      }
  
  
      //starting jet loop Gen
      if(isMC)
      {  
        for(int j=0; j<ngen; j++)
        {  

          //Potential Problem 1

          if(v==1) genpt[j] = genpt[j]*1.015; 
          if(v==3 || v==5)genpt[j] = genpt[j]*1.025;
          //if(v==14 || v==16 || v==18)genpt[j] = genpt[j]*1.01;
          if(v==2) genpt[j] = genpt[j]*0.985; 
          if(v==4 || v==6)genpt[j] = genpt[j]*0.975; 
          //if(v==15 || v==17 || v==19)genpt[j] = genpt[j]*0.99;
          if(v==7 || v==8 || v==9)genpt[j] = getJERCorrected(mode,genpt[j],0.05);
          if(v==34) genpt[j] = genpt[j]*1.025;
          if(v==37) genpt[j] = genpt[j]*1.015;
          if(v==35) genpt[j] = genpt[j]*0.975;
          if(v==38) genpt[j] = genpt[j]*0.985;
          if(v==36) genpt[j] = getJERCorrected("pPb5",genpt[j],0.05); 
          //smearing gen to match reco resolution (for ratio comparisons to be fair)
          //(smears to pPb resolution, first ppref5 is a dummy argument)
          genpt[j] = getJERCorrected("ppref5",genpt[j],getGenMCSmearFactor("pPb5",genpt[j])); 

          if(TMath::Abs(geneta[j]+boost) < jetEtaMin || TMath::Abs(geneta[j]+boost) > (doMidRapidity?0.25:jetEtaMax) || genpt[j]<lowJetPtBound || genpt[j]>=upJetPtBound) continue;
          if(typeUE==1 && (TMath::Abs(geneta[j]+boost))<coneSize) continue;
          
          //getting jet flavor (a bit convoluted because genMatchedID is not filled in forest correctly)
          bool isQ=false;
          bool isG=false;
          for(int j2=0; j2<nref; j2++)
          {
             if(TMath::Abs(refpt[j2] - genpt[j])==0 && TMath::Abs(refeta[j2] - geneta[j])==0 && refparton_flavor[j2]==21) isG = true;
             if(TMath::Abs(refpt[j2] - genpt[j])==0 && TMath::Abs(refeta[j2] - geneta[j])==0 && refparton_flavor[j2]!=21 && TMath::Abs(refparton_flavor[j2])<901) isQ = true; 
             if(isQ || isG) break;
          }
            
          if(isQ) h_jet_gen_Q->Fill(genpt[j],weight);
          if(isG) h_jet_gen_G->Fill(genpt[j],weight);
          h_jet_gen->Fill(genpt[j], weight);  
       
          for(int t=0; t<nParticle; t++)
          { 
            if(pPt[t] <= 0.5 || pPt[t] > 400 || TMath::Abs(pEta[t])>2.4 ) continue;
  
            //Filling track spectrum in jet cone
            if(getdR2(geneta[j]+boost,genphi[j],pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
            {
              totalGenTGenJ_Tracking->Fill(pPt[t],weight);
              h_track_gen->Fill(genpt[j],pPt[t],weight);
              h_track_xi_gen->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
   
              if(isQ)
              {
                h_track_gen_Q->Fill(genpt[j], pPt[t],weight);
                h_track_xi_gen_Q->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
              if(isG)
              {
                h_track_gen_G->Fill(genpt[j],pPt[t],weight);
                h_track_xi_gen_G->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
            }
       
            //Phi rotated UE subtraction
            //returns either +-1 to rotate clockwise or ccw randomly
            double rotationDirection = 2*(int)rand->Integer(2)-1; 
            if(v==26) rotationDirection = rotationDirection*2.0/3.0;

            if(typeUE==0 && getdR2(geneta[j]+boost,genphi[j]+rotationDirection*TMath::PiOver2(),pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
            {
              h_trackUE_gen->Fill(genpt[j],pPt[t],weight);  
              h_trackUE_xi_gen->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j]+rotationDirection*TMath::PiOver2(),pPt[t],pEta[t]+boost,pPhi[t]),weight);
   
              if(isQ)
              {
                h_trackUE_gen_Q->Fill(genpt[j], pPt[t],weight);
                h_trackUE_xi_gen_Q->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j]+rotationDirection*TMath::PiOver2(),pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
              if(isG)
              {
                h_trackUE_gen_G->Fill(genpt[j],pPt[t],weight);
                h_trackUE_xi_gen_G->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j]+rotationDirection*TMath::PiOver2(),pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
            }
  
            //Eta Reflected UE subtraction
            if(typeUE==1 && getdR2(-1*(geneta[j]+boost),genphi[j],pEta[t]+boost,pPhi[t]) < coneSize*coneSize)
            {
              h_trackUE_gen->Fill(genpt[j],pPt[t],weight); 
              h_trackUE_xi_gen->Fill(genpt[j],getXi(genpt[j],-1*(geneta[j]+boost),genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
  
              if(isQ)
              {
                h_trackUE_gen_Q->Fill(genpt[j],pPt[t],weight);
                h_trackUE_xi_gen_Q->Fill(genpt[j],getXi(genpt[j],-1*(geneta[j]+boost),genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
              if(isG)
              {
                h_trackUE_gen_G->Fill(genpt[j],pPt[t],weight);
                h_trackUE_xi_gen_G->Fill(genpt[j],getXi(genpt[j],-1*(geneta[j]+boost),genphi[j],pPt[t],pEta[t]+boost,pPhi[t]),weight);
              }
            }
          }
  
          //UE subtraction w/ MB mixing
          if(typeUE==2)
          {        
            getInputEntryMix(lastMixEvt);
            for(int t = 0; t < nParticleMix; t++)
            {
              if(pPtMix[t] <= 0.5 || pPtMix[t] > 400 || TMath::Abs(pEtaMix[t])>2.4 ) continue;
              
              //Filling track spectrum in jet cone
              if(getdR2(geneta[j]+boost,genphi[j],pEtaMix[t]+boost,pPhiMix[t]) < coneSize*coneSize)
              {
                h_trackUE_gen->Fill(genpt[j],pPtMix[t],weight);
                h_trackUE_xi_gen->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPtMix[t],pEtaMix[t]+boost,pPhiMix[t]),weight);
                
                if(isQ)
                {
                  h_trackUE_gen_Q->Fill(genpt[j],pPtMix[t],weight);
                  h_trackUE_xi_gen_Q->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPtMix[t],pEtaMix[t]+boost,pPhiMix[t]),weight);
                }
                if(isG)
                {
                  h_trackUE_gen_G->Fill(genpt[j],pPtMix[t],weight);
                  h_trackUE_xi_gen_G->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],pPtMix[t],pEtaMix[t]+boost,pPhiMix[t]),weight);
                }
              }
            }
          }
          
          //gen jet reco track
          for(int t=0; t<nTrk; t++)
          {         
            if((trkCharge[t]!=1 && v==27) || (trkCharge[t]!=-1 && v==28)) continue;     
            if(trkPt[t] <= 0.5 || trkPt[t] > 400 || !highPurity[t] || TMath::Abs(trkEta[t])>2.4 ) continue;
            if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.1)) continue;            
            else if(strcmp(mode,"ppref5")==0 && !ispPbStyleCorr && (TMath::Abs(trkDxy1[t]/trkDxyError1[t]) > 3 || TMath::Abs(trkDz1[t]/trkDzError1[t]) > 3 || trkPtError[t]/trkPt[t] > 0.3 || ((pfEcal[t]+pfHcal[t])/TMath::CosH(trkEta[t])<0.5*trkPt[t] && trkPt[t]>20))) continue;            
            //calculating r_min for tracking correction  
            double r_min = 9;
            double rmin_pt = 60;
            for(int j2 = 0; j2<nref; j2++)
            {
              if(TMath::Abs(jteta[j2])>2 || TMath::Abs(jtpt[j2]) < 60 || jtpt[j2]>1000) continue;
              double r_min_temp = TMath::Power(getdR2(jteta[j2],jtphi[j2],trkEta[t],trkPhi[t]),0.5);
              if(r_min_temp < r_min){ r_min = r_min_temp;
              rmin_pt = jtpt[j2];}
            }
            r_min = rmin_pt;
   
            //Filling track spectrum in jet cone
            if(getdR2(geneta[j]+boost,genphi[j],trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
            {
              double trkCorr = 1;
              if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
              else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,genpt[j]);     
            
              if(v==13) trkCorr=1;          
              if(std::isfinite(trkCorr))
              {
                totalRecoTGenJ_Tracking->Fill(trkPt[t],weight*trkCorr);
                h_track_gJrT->Fill(genpt[j],trkPt[t],trkCorr*weight);
                h_track_xi_gJrT->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
            }
       
            //Phi rotated UE subtraction
            //returns either +-1 to rotate clockwise or ccw randomly
            double rotationDirection = 2*(int)rand->Integer(2)-1;
            if(v==26) rotationDirection = rotationDirection*2.0/3.0;

            if(typeUE==0 && getdR2(geneta[j]+boost,genphi[j]+rotationDirection*TMath::PiOver2(),trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
            {
              double trkCorr = 1;
              if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
              else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,genpt[j]);     
              
              if(v==13) trkCorr=1;
              if(std::isfinite(trkCorr))
              {
                h_trackUE_gJrT->Fill(genpt[j],trkPt[t],trkCorr*weight);  
                h_trackUE_xi_gJrT->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j]+rotationDirection*TMath::Pi(),trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight);
              }
            }
  
            //Eta Reflected UE subtraction
            if(typeUE==1 && getdR2(-1*(geneta[j]+boost),genphi[j],trkEta[t]+boost,trkPhi[t]) < coneSize*coneSize)
            {
              double trkCorr = 1;
              if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPt[t], sType), 1, trkPt[t], trkPhi[t], trkEta[t], r_min, sType);     
              else trkCorr = trkCorrObj->getTrkCorr(trkPt[t],trkEta[t],trkPhi[t],1,r_min,genpt[j]);     
              
              if(v==13) trkCorr=1;
              if(std::isfinite(trkCorr))
              {
                h_trackUE_gJrT->Fill(genpt[j],trkPt[t],trkCorr*weight); 
                h_trackUE_xi_gJrT->Fill(genpt[j],getXi(genpt[j],-1*(geneta[j]+boost),genphi[j],trkPt[t],trkEta[t]+boost,trkPhi[t]),trkCorr*weight); 
              }
            }
          }
  
          //UE subtraction w/ MB mixing
          if(typeUE==2)
          {        
            getInputEntryMix(lastMixEvt);
            for(int t = 0; t < nTrkMix; t++)
            {
              if((trkChargeMix[t]!=1 && v==27) || (trkChargeMix[t]!=-1 && v==28)) continue;     
              if(trkPtMix[t] <= 0.5 || trkPtMix[t] > 400 || !highPurityMix[t] || TMath::Abs(trkEtaMix[t])>2.4 ) continue;
              if((!(strcmp(mode,"ppref5")==0) || ispPbStyleCorr) && (TMath::Abs(trkDxy1Mix[t]/trkDxyError1Mix[t]) > 3 || TMath::Abs(trkDz1Mix[t]/trkDzError1Mix[t]) > 3 || trkPtErrorMix[t]/trkPtMix[t] > 0.1)) continue;            
              else if(strcmp(mode,"ppref5")==0 && !ispPbStyleCorr && (TMath::Abs(trkDxy1Mix[t]/trkDxyError1Mix[t]) > 3 || TMath::Abs(trkDz1Mix[t]/trkDzError1Mix[t]) > 3 || trkPtErrorMix[t]/trkPtMix[t] > 0.3 || ((pfEcalMix[t]+pfHcalMix[t])/TMath::CosH(trkEtaMix[t])<0.5*trkPtMix[t] && trkPtMix[t]>20))) continue;            
  
              //calculating r_min for tracking correction                                                   
              double r_min = 9;
              double rmin_pt = 60;
              for(int j2 = 0; j2<nrefMix; j2++)
              {
                if(TMath::Abs(jtetaMix[j2])>2 || TMath::Abs(jtptMix[j2]) < 60 || jtptMix[j2]>1000) continue;
                double r_min_temp = TMath::Power(getdR2(jtetaMix[j2],jtphiMix[j2],trkEtaMix[t],trkPhiMix[t]),0.5);
                if(r_min_temp < r_min){ r_min = r_min_temp;
                rmin_pt = jtptMix[j2];}
              }                                                                          
              r_min = rmin_pt;
              
              //Filling track spectrum in jet cone
              if(getdR2(geneta[j]+boost,genphi[j],trkEtaMix[t]+boost,trkPhiMix[t]) < coneSize*coneSize)
              {
                double trkCorr = 1;
                if(!(strcmp(mode,"ppref5")==0)) trkCorr = factorizedPtCorr(getPtBin(trkPtMix[t], sType), 1, trkPtMix[t], trkPhiMix[t], trkEtaMix[t], r_min, sType);     
                else trkCorr = trkCorrObj->getTrkCorr(trkPtMix[t],trkEtaMix[t],trkPhiMix[t],1,r_min,genpt[j]);     
                
                if(v==13) trkCorr=1;
                if(std::isfinite(trkCorr))
                {
                  h_trackUE_gJrT->Fill(genpt[j],trkPtMix[t],trkCorr*weight);
                  h_trackUE_xi_gJrT->Fill(genpt[j],getXi(genpt[j],geneta[j]+boost,genphi[j],trkPtMix[t],trkEtaMix[t]+boost,trkPhiMix[t]),trkCorr*weight);
                }
              }
            }
          }
        }
      }
    }
    TFile * outf = new TFile(Form("spectra%s%s_%d_%d_UE%d_%d_%d.root",mode,trigger,jobNum,(int)isMC,(int)typeUE,(int)(10*jetEtaMin),(int)(10*jetEtaMax)),"update");
    h_jet->SetDirectory(0);
    h_track->SetDirectory(0);
    h_trackUE->SetDirectory(0);
    h_track_xi->SetDirectory(0);
    h_trackUE_xi->SetDirectory(0);
    h_track_rJgT->SetDirectory(0);
    h_trackUE_rJgT->SetDirectory(0);
    h_track_xi_rJgT->SetDirectory(0);
    h_trackUE_xi_rJgT->SetDirectory(0);
    h_track_gJrT->SetDirectory(0);
    h_trackUE_gJrT->SetDirectory(0);
    h_track_xi_gJrT->SetDirectory(0);
    h_trackUE_xi_gJrT->SetDirectory(0);
    h_jet_gen->SetDirectory(0);
    h_track_gen->SetDirectory(0);
    h_trackUE_gen->SetDirectory(0);
    h_track_xi_gen->SetDirectory(0);
    h_trackUE_xi_gen->SetDirectory(0);
    h_jet_Q->SetDirectory(0);
    h_track_Q->SetDirectory(0);
    h_trackUE_Q->SetDirectory(0);
    h_track_xi_Q->SetDirectory(0);
    h_trackUE_xi_Q->SetDirectory(0);
    h_jet_G->SetDirectory(0);
    h_track_G->SetDirectory(0);
    h_trackUE_G->SetDirectory(0);
    h_track_xi_G->SetDirectory(0);
    h_trackUE_xi_G->SetDirectory(0);
    h_jet_gen_Q->SetDirectory(0);
    h_track_gen_Q->SetDirectory(0);
    h_trackUE_gen_Q->SetDirectory(0);
    h_track_xi_gen_Q->SetDirectory(0);
    h_trackUE_xi_gen_Q->SetDirectory(0);
    h_jet_gen_G->SetDirectory(0);
    h_track_gen_G->SetDirectory(0);
    h_trackUE_gen_G->SetDirectory(0);
    h_track_xi_gen_G->SetDirectory(0);
    h_trackUE_xi_gen_G->SetDirectory(0);
    h_JEC->SetDirectory(0);  
    h_JEC_weights->SetDirectory(0);  
    h_trackVsJEC->SetDirectory(0);  
    h_trackVsJEC_weights->SetDirectory(0);  

    if(typeUE==2 && v==0){
      //highZskim->SetDirectory(0);
      //highZskim->Write();
      h_L2ResidualCorr->SetDirectory(0);
      h_L2ResidualCorr->Write(); 
      h_L2ResidualCorr_eta->SetDirectory(0);
      h_L2ResidualCorr_eta->Write();
      nLeadJet40->SetDirectory(0); 
      nLeadJet60->SetDirectory(0); 
      nLeadJet80->SetDirectory(0); 
      nLeadJet40->Write(); 
      nLeadJet60->Write(); 
      nLeadJet80->Write(); 
    }

    h_jet->Write(Form("%s_reco_jet%s",mode,variationTag[v]));
    h_track->Write(Form("%s_reco_track%s",mode,variationTag[v]));
    h_trackUE->Write(Form("%s_reco_trackUE%s",mode,variationTag[v]));
    h_track_xi->Write(Form("%s_reco_track_xi%s",mode,variationTag[v]));
    h_trackUE_xi->Write(Form("%s_reco_trackUE_xi%s",mode,variationTag[v]));  
    h_track_rJgT->Write(Form("%s_rJgT_track%s",mode,variationTag[v]));
    h_trackUE_rJgT->Write(Form("%s_rJgT_trackUE%s",mode,variationTag[v]));
    h_track_xi_rJgT->Write(Form("%s_rJgT_track_xi%s",mode,variationTag[v]));
    h_trackUE_xi_rJgT->Write(Form("%s_rJgT_trackUE_xi%s",mode,variationTag[v]));
    h_track_gJrT->Write(Form("%s_gJrT_track%s",mode,variationTag[v]));
    h_trackUE_gJrT->Write(Form("%s_gJrT_trackUE%s",mode,variationTag[v]));
    h_track_xi_gJrT->Write(Form("%s_gJrT_track_xi%s",mode,variationTag[v]));
    h_trackUE_xi_gJrT->Write(Form("%s_gJrT_trackUE_xi%s",mode,variationTag[v]));
    h_jet_gen->Write(Form("%s_gen_jet%s",mode,variationTag[v]));
    h_track_gen->Write(Form("%s_gen_track%s",mode,variationTag[v]));
    h_trackUE_gen->Write(Form("%s_gen_trackUE%s",mode,variationTag[v]));
    h_track_xi_gen->Write(Form("%s_gen_track_xi%s",mode,variationTag[v]));
    h_trackUE_xi_gen->Write(Form("%s_gen_trackUE_xi%s",mode,variationTag[v]));
    h_jet_Q->Write(Form("%s_reco_jet_Q%s",mode,variationTag[v]));
    h_track_Q->Write(Form("%s_reco_track_Q%s",mode,variationTag[v]));
    h_trackUE_Q->Write(Form("%s_reco_trackUE_Q%s",mode,variationTag[v]));
    h_track_xi_Q->Write(Form("%s_reco_track_xi_Q%s",mode,variationTag[v]));
    h_trackUE_xi_Q->Write(Form("%s_reco_trackUE_xi_Q%s",mode,variationTag[v]));
    h_jet_G->Write(Form("%s_reco_jet_G%s",mode,variationTag[v]));
    h_track_G->Write(Form("%s_reco_track_G%s",mode,variationTag[v]));
    h_trackUE_G->Write(Form("%s_reco_trackUE_G%s",mode,variationTag[v]));
    h_track_xi_G->Write(Form("%s_reco_track_xi_G%s",mode,variationTag[v]));
    h_trackUE_xi_G->Write(Form("%s_reco_trackUE_xi_G%s",mode,variationTag[v]));
    h_jet_gen_Q->Write(Form("%s_gen_jet_Q%s",mode,variationTag[v]));
    h_track_gen_Q->Write(Form("%s_gen_track_Q%s",mode,variationTag[v]));
    h_trackUE_gen_Q->Write(Form("%s_gen_trackUE_Q%s",mode,variationTag[v]));
    h_track_xi_gen_Q->Write(Form("%s_gen_track_xi_Q%s",mode,variationTag[v]));
    h_trackUE_xi_gen_Q->Write(Form("%s_gen_trackUE_xi_Q%s",mode,variationTag[v]));
    h_jet_gen_G->Write(Form("%s_gen_jet_G%s",mode,variationTag[v]));
    h_track_gen_G->Write(Form("%s_gen_track_G%s",mode,variationTag[v]));
    h_trackUE_gen_G->Write(Form("%s_gen_trackUE_G%s",mode,variationTag[v]));
    h_track_xi_gen_G->Write(Form("%s_gen_track_xi_G%s",mode,variationTag[v]));
    h_trackUE_xi_gen_G->Write(Form("%s_gen_trackUE_xi_G%s",mode,variationTag[v]));
    h_trackVsJEC->Write(Form("%s_trackVsJEC_%s",mode,variationTag[v]));  
    h_trackVsJEC_weights->Write(Form("%s_trackVsJEC_weights_%s",mode,variationTag[v]));  
    h_JEC->Write(Form("%s_JEC_%s",mode,variationTag[v]));  
    h_JEC_weights->Write(Form("%s_JEC_weights_%s",mode,variationTag[v]));  
    if(v==0)
    {
      totalEvts->Write();
      totalJetsHist->Write();
      totalJetsRawPtCut->Write();
      totalJetsChargeSumCut->Write();
      totalJetsEtaCutHist->Write();
      totalJetsPtCutHist->Write();
    }
    if(v==0 || v==31 || v==32)
    {
      totalGenTRecoJ_Tracking->Write(Form("totalGenTRecoJ_Tracking%s",variationTag[v]));
      totalRecoTRecoJ_Tracking->Write(Form("totalRecoTRecoJ_Tracking%s",variationTag[v]));
      totalGenTGenJ_Tracking->Write(Form("totalGenTGenJ_Tracking%s",variationTag[v]));
      totalRecoTGenJ_Tracking->Write(Form("totalRecoTGenJ_Tracking%s",variationTag[v]));
      totalGenTRecoJ_Tracking_eta->Write(Form("totalGenTRecoJ_Tracking_eta%s",variationTag[v]));
      totalRecoTRecoJ_Tracking_eta->Write(Form("totalRecoTRecoJ_Tracking_eta%s",variationTag[v]));
      totalRecoTRecoJ_Tracking_largeEta->Write(Form("totalRecoTRecoJ_Tracking_largeEta%s",variationTag[v]));
      totalGenTRecoJ_Tracking_largeEta->Write(Form("totalGenTRecoJ_Tracking_largeEta%s",variationTag[v]));
      totalRecoTRecoJ_Tracking_lowEta->Write(Form("totalRecoTRecoJ_Tracking_lowEta%s",variationTag[v]));
      totalGenTRecoJ_Tracking_lowEta->Write(Form("totalGenTRecoJ_Tracking_lowEta%s",variationTag[v]));
      totalRecoTRecoJ_Tracking_largePt->Write(Form("totalRecoTRecoJ_Tracking_largePt%s",variationTag[v]));
      totalGenTRecoJ_Tracking_largePt->Write(Form("totalGenTRecoJ_Tracking_largePt%s",variationTag[v]));
    }
    totalGenTRecoJ_Tracking->Reset();
    totalRecoTRecoJ_Tracking->Reset();
    totalGenTGenJ_Tracking->Reset();
    totalRecoTGenJ_Tracking->Reset();
    totalRecoTRecoJ_Tracking_eta->Reset();
    totalGenTRecoJ_Tracking_eta->Reset();
    totalRecoTRecoJ_Tracking_largeEta->Reset();
    totalGenTRecoJ_Tracking_largeEta->Reset();
    totalRecoTRecoJ_Tracking_lowEta->Reset();
    totalGenTRecoJ_Tracking_lowEta->Reset();
    totalRecoTRecoJ_Tracking_largePt->Reset();
    totalGenTRecoJ_Tracking_largePt->Reset();
  
    outf->Close();
  }
}

int main(int argc, const char* argv[])
{
  if(argc != 5)
  {
    std::cout << "Usage: Spectra <fileListJets> <fileListMB> <job> <UEMode>" << std::endl;
    return 1;
  }
  
  std::string fList = argv[1];
  std::string fListMB = argv[2];
  int job = std::atoi(argv[3]);
  int UEMode = std::atoi(argv[4]);

  /*if(job==0) job=27;
  else if(job==1) job=26;
  else if(job==2) job=27;*/
 
  std::string buffer, bufferMB;
  std::vector<std::string> listOfFilesJets;
  std::vector<std::string> listOfFilesMB;  

  std::ifstream inFile(fList.data());
  std::ifstream inFileMB(fListMB.data());
  std::cout << fList << std::endl;
  std::cout << fListMB << std::endl;
  std::cout << inFile.is_open() << std::endl;
  std::cout << inFileMB.is_open() << std::endl;

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return 1;
  }
  else if (!inFileMB.is_open())
  {
    std::cout << "Error opening MB file. Exiting." <<std::endl;
    return 1;
  }
  else
  {
    while(true)
    {
      inFile >> buffer;
      inFileMB >> bufferMB;
      if(inFile.eof()) break;
      if(inFileMB.eof()) break;
      listOfFilesJets.push_back(buffer);
      listOfFilesMB.push_back(bufferMB);
    }
  }

  std::cout << "FileListJets Loaded" << std::endl;
  std::cout << "FileJobJets: " << listOfFilesJets[job] << std::endl;
  std::cout << "FileListMB Loaded" << std::endl;
  std::cout << "FileJobMB: " << listOfFilesMB[job] << std::endl;

  std::string  parsedMode;
  std::string  parsedTrigger;
  if(listOfFilesJets[job].find("pPb5") != std::string::npos) parsedMode = "pPb5"; 
  if(listOfFilesJets[job].find("Pbp5") != std::string::npos) parsedMode = "Pbp5";
  if(listOfFilesJets[job].find("pp5") != std::string::npos) parsedMode = "pp5";
  if(listOfFilesJets[job].find("ppref5") != std::string::npos) parsedMode = "ppref5";

  if(listOfFilesJets[job].find("jet80_") != std::string::npos) parsedTrigger = "jet80";
  if(listOfFilesJets[job].find("jet40_") != std::string::npos) parsedTrigger = "jet40";
  if(listOfFilesJets[job].find("jet60_") != std::string::npos) parsedTrigger = "jet60";

  int MCStatus = 0;
  if(listOfFilesJets[job].find("/MC/") != std::string::npos) MCStatus = 1;

  std::cout << "Results of parsing input files for mode and trigger used:" << std::endl;
  std::cout << "Mode: " <<  parsedMode.data() << "  Trigger: " << parsedTrigger.data() << "  Is MC:" << MCStatus << std::endl;
  
  Spectra(listOfFilesJets[job].data(),listOfFilesMB[job].data(),parsedMode.data(),parsedTrigger.data(),MCStatus,job,UEMode,0,1.5);
  return 0;
}

