#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>

const char* filePath = "../mergedHists/Oct25_systChecks/";
const int variations = 37;
const char * variationTag[variations]= {"","_pp2JESUP3","_pp2JESDOWN3","_pp7JESUP3","_pp7JESDOWN3","_pPb5JESUP3","_pPb5JESDOWN3","_pp2JER5","_pp7JER5","_pPb5JER5","_pp2JER2","_pp7JER2","_pPb5JER2","_NoTrackCorr","_pp2JESUP1","_pp2JESDOWN1","_pp7JESUP1","_pp7JESDOWN1","_pPb5JESUP1","_pPb5JESDOWN1","_pp2JESUP2","_pp2JESDOWN2","_pp7JESUP2","_pp7JESDOWN2","_pPb5JESUP2","_pPb5JESDOWN2","_60DegreeCone","_ChargePlus","_ChargeMinus","_XtScaled","_NoChargeCut","_lowPU","_highPU","_midPU","_pprefJESUP2p5","_pprefJESDOWN2p5","_ppref5JER5"};

// jet pt boundaries
const int FF_Bins = 5;
double FF_Bound[FF_Bins+1] = {60,80,100,120,140,200};

//interpolation histos

//plotting histos
TH1D * pPb5TeV_data[2*FF_Bins];
TH1D * Pbp5TeV_data[2*FF_Bins];
TH1D * pPb5Pbp5TeV_fulldata[2*FF_Bins];
TH1D * ppref5TeV_data[2*FF_Bins];

TH1D * pPb5TeV_recoMC[2*FF_Bins];
TH1D * Pbp5TeV_recoMC[2*FF_Bins];
TH1D * pPb5Pbp5TeV_recoMC[2*FF_Bins];
TH1D * pp5TeV_recoMC[2*FF_Bins];

TH1D * pPb5TeV_genMC[2*FF_Bins];
TH1D * Pbp5TeV_genMC[2*FF_Bins];
TH1D * pPb5Pbp5TeV_genMC[2*FF_Bins];
TH1D * pp5TeV_genMC[2*FF_Bins];

TH1D * pPb5TeV_rJgTMC[2*FF_Bins];
TH1D * Pbp5TeV_rJgTMC[2*FF_Bins];
TH1D * pPb5Pbp5TeV_rJgTMC[2*FF_Bins];
TH1D * pp5TeV_rJgTMC[2*FF_Bins];
TH1D * pPb5TeV_gJrTMC[2*FF_Bins];
TH1D * Pbp5TeV_gJrTMC[2*FF_Bins];
TH1D * pPb5Pbp5TeV_gJrTMC[2*FF_Bins];
TH1D * pp5TeV_gJrTMC[2*FF_Bins];

TH1D * pPb_FF[2*FF_Bins];
TH1D * Pbp_FF[2*FF_Bins];
TH1D * pPb_FF_recoMC[2*FF_Bins];
TH1D * pPb_FF_genMC[2*FF_Bins];
TH1D * pPb_FF_rJgTMC[2*FF_Bins];
TH1D * pPb_FF_gJrTMC[2*FF_Bins];
TH1D * Pbp_FF_recoMC[2*FF_Bins];
TH1D * Pbp_FF_genMC[2*FF_Bins];
TH1D * Pbp_FF_rJgTMC[2*FF_Bins];
TH1D * Pbp_FF_gJrTMC[2*FF_Bins];
TH1D * pPbPbp_FF_recoMC[2*FF_Bins];
TH1D * pPbPbp_FF_genMC[2*FF_Bins];
TH1D * pPbPbp_FF_rJgTMC[2*FF_Bins];
TH1D * pPbPbp_FF_gJrTMC[2*FF_Bins];
TH1D * pp5_FF_recoMC[2*FF_Bins];
TH1D * pp5_FF_genMC[2*FF_Bins];
TH1D * pPbPbp_FF[2*FF_Bins];


TH1D * pPb5_0_jet;
TH2D * pPb5_0_track;
TH2D * pPb5_0_trackUE;
TH2D * pPb5_0_track_xi;
TH2D * pPb5_0_trackUE_xi;

TH1D * Pbp5_0_jet;
TH2D * Pbp5_0_track;
TH2D * Pbp5_0_trackUE;
TH2D * Pbp5_0_track_xi;
TH2D * Pbp5_0_trackUE_xi;

TH1D * pPb5Pbp5_0_jet;
TH2D * pPb5Pbp5_0_track;
TH2D * pPb5Pbp5_0_trackUE;
TH2D * pPb5Pbp5_0_track_xi;
TH2D * pPb5Pbp5_0_trackUE_xi;

TH1D * ppref5_0_jet;
TH2D * ppref5_0_track;
TH2D * ppref5_0_trackUE;
TH2D * ppref5_0_track_xi;
TH2D * ppref5_0_trackUE_xi;

TH1D * pPb5_1_jet_reco;
TH2D * pPb5_1_track_reco;
TH2D * pPb5_1_trackUE_reco;
TH2D * pPb5_1_track_xi_reco;
TH2D * pPb5_1_trackUE_xi_reco;


TH2D * pPb5_1_track_rJgT;
TH2D * pPb5_1_trackUE_rJgT;
TH2D * pPb5_1_track_xi_rJgT;
TH2D * pPb5_1_trackUE_xi_rJgT;
TH2D * pPb5_1_track_gJrT;
TH2D * pPb5_1_trackUE_gJrT;
TH2D * pPb5_1_track_xi_gJrT;
TH2D * pPb5_1_trackUE_xi_gJrT;

TH1D * Pbp5_1_jet_reco;
TH2D * Pbp5_1_track_reco;
TH2D * Pbp5_1_trackUE_reco;
TH2D * Pbp5_1_track_xi_reco;
TH2D * Pbp5_1_trackUE_xi_reco;

TH2D * Pbp5_1_track_rJgT;
TH2D * Pbp5_1_trackUE_rJgT;
TH2D * Pbp5_1_track_xi_rJgT;
TH2D * Pbp5_1_trackUE_xi_rJgT;
TH2D * Pbp5_1_track_gJrT;
TH2D * Pbp5_1_trackUE_gJrT;
TH2D * Pbp5_1_track_xi_gJrT;
TH2D * Pbp5_1_trackUE_xi_gJrT;

TH1D * pPb5Pbp5_1_jet_reco;
TH2D * pPb5Pbp5_1_track_reco;
TH2D * pPb5Pbp5_1_trackUE_reco;
TH2D * pPb5Pbp5_1_track_xi_reco;
TH2D * pPb5Pbp5_1_trackUE_xi_reco;


TH2D * pPb5Pbp5_1_track_rJgT;
TH2D * pPb5Pbp5_1_trackUE_rJgT;
TH2D * pPb5Pbp5_1_track_xi_rJgT;
TH2D * pPb5Pbp5_1_trackUE_xi_rJgT;
TH2D * pPb5Pbp5_1_track_gJrT;
TH2D * pPb5Pbp5_1_trackUE_gJrT;
TH2D * pPb5Pbp5_1_track_xi_gJrT;
TH2D * pPb5Pbp5_1_trackUE_xi_gJrT;


TH1D * pPb5_1_jet_gen;
TH2D * pPb5_1_track_gen;
TH2D * pPb5_1_trackUE_gen;
TH2D * pPb5_1_track_xi_gen;
TH2D * pPb5_1_trackUE_xi_gen;


TH1D * Pbp5_1_jet_gen;
TH2D * Pbp5_1_track_gen;
TH2D * Pbp5_1_trackUE_gen;
TH2D * Pbp5_1_track_xi_gen;
TH2D * Pbp5_1_trackUE_xi_gen;


TH1D * pPb5Pbp5_1_jet_gen;
TH2D * pPb5Pbp5_1_track_gen;
TH2D * pPb5Pbp5_1_trackUE_gen;
TH2D * pPb5Pbp5_1_track_xi_gen;
TH2D * pPb5Pbp5_1_trackUE_xi_gen;


TH1D * pp5_1_jet_reco;
TH2D * pp5_1_track_reco;
TH2D * pp5_1_trackUE_reco;
TH2D * pp5_1_track_xi_reco;
TH2D * pp5_1_trackUE_xi_reco;


TH2D * pp5_1_track_rJgT;
TH2D * pp5_1_trackUE_rJgT;
TH2D * pp5_1_track_xi_rJgT;
TH2D * pp5_1_trackUE_xi_rJgT;
TH2D * pp5_1_track_gJrT;
TH2D * pp5_1_trackUE_gJrT;
TH2D * pp5_1_track_xi_gJrT;
TH2D * pp5_1_trackUE_xi_gJrT;

TH1D * pp5_1_jet_gen;
TH2D * pp5_1_track_gen;
TH2D * pp5_1_trackUE_gen;
TH2D * pp5_1_track_xi_gen;
TH2D * pp5_1_trackUE_xi_gen;



TH1D *Njets_pPbData;
TH1D *Njets_pPbMC;
TH1D *Njets_PbpData;
TH1D *Njets_PbpMC;

void loadHistos(int v, int UEtype)
{
  TFile * spectraFilepPb5 = new TFile(Form("%s/pPb5_UE%d_0_15.root",filePath,UEtype),"read");
  TFile * spectraFilePbp5 = new TFile(Form("%s/Pbp5_UE%d_0_15.root",filePath,UEtype),"read");
  TFile * spectraFileppref5 = new TFile(Form("%s/ppref5_UE%d_0_15.root",filePath,UEtype),"read");
  TFile * MCFilepPb5 = new TFile(Form("%s/pPb5MC_UE%d_0_15.root",filePath,UEtype),"read");
  TFile * MCFilePbp5 = new TFile(Form("%s/Pbp5MC_UE%d_0_15.root",filePath,UEtype),"read");
  //TFile * MCFilepp5 = new TFile(Form("%s/pp5MC_UE%d_0_15.root",filePath,UEtype),"read");
  TFile * MCFilepp5 = new TFile(Form("%s/ppref5MC_UE%d_0_15.root",filePath,UEtype),"read");

  std::string pp2Tag = "", pp7Tag = "", pPb5Tag = "";
  if(v==1 || v==2 || v==7 || v==10 || v==13 || v==14 || v==15 || v==20 || v==21 || v==26 || v==27 || v==28 || v==29 || v==30) pp2Tag = variationTag[v];
  if(v==3 || v==4 || v==8 || v==11 || v==13 || v==16 || v==17 || v==22 || v==23 || v==26 || v==27 || v==28 || v==29 || v==30 || v==31 || v==32 || v==33) pp7Tag = variationTag[v];
  if(v==5 || v==6 || v==9 || v==12 || v==13 || v==18 || v==19 || v==24 || v==25 || v==26 || v==27 || v==28 || v==29 || v==30 || v==31 || v==34 || v==35 || v==36) pPb5Tag = variationTag[v];

  pPb5_0_jet = (TH1D*) spectraFilepPb5->Get(Form("pPb5_reco_jet%s",pPb5Tag.data()));
  pPb5_0_track = (TH2D*) spectraFilepPb5->Get(Form("pPb5_reco_track%s",pPb5Tag.data()));
  pPb5_0_trackUE = (TH2D*) spectraFilepPb5->Get(Form("pPb5_reco_trackUE%s",pPb5Tag.data()));
  pPb5_0_track_xi = (TH2D*) spectraFilepPb5->Get(Form("pPb5_reco_track_xi%s",pPb5Tag.data()));
  pPb5_0_trackUE_xi = (TH2D*) spectraFilepPb5->Get(Form("pPb5_reco_trackUE_xi%s",pPb5Tag.data()));

  Pbp5_0_jet = (TH1D*) spectraFilePbp5->Get(Form("Pbp5_reco_jet%s",pPb5Tag.data()));
  Pbp5_0_track = (TH2D*) spectraFilePbp5->Get(Form("Pbp5_reco_track%s",pPb5Tag.data()));
  Pbp5_0_trackUE = (TH2D*) spectraFilePbp5->Get(Form("Pbp5_reco_trackUE%s",pPb5Tag.data()));
  Pbp5_0_track_xi = (TH2D*) spectraFilePbp5->Get(Form("Pbp5_reco_track_xi%s",pPb5Tag.data()));
  Pbp5_0_trackUE_xi = (TH2D*) spectraFilePbp5->Get(Form("Pbp5_reco_trackUE_xi%s",pPb5Tag.data()));
  
  ppref5_0_jet = (TH1D*) spectraFileppref5->Get(Form("ppref5_reco_jet%s",pPb5Tag.data()));
  ppref5_0_track = (TH2D*) spectraFileppref5->Get(Form("ppref5_reco_track%s",pPb5Tag.data()));
  ppref5_0_trackUE = (TH2D*) spectraFileppref5->Get(Form("ppref5_reco_trackUE%s",pPb5Tag.data()));
  ppref5_0_track_xi = (TH2D*) spectraFileppref5->Get(Form("ppref5_reco_track_xi%s",pPb5Tag.data()));
  ppref5_0_trackUE_xi = (TH2D*) spectraFileppref5->Get(Form("ppref5_reco_trackUE_xi%s",pPb5Tag.data()));

  pPb5Pbp5_0_jet = (TH1D*)pPb5_0_jet->Clone("pPb5Pbp5_reco_jet");
  pPb5Pbp5_0_track = (TH2D*)pPb5_0_track->Clone("pPb5Pbp5_reco_track");
  pPb5Pbp5_0_trackUE = (TH2D*)pPb5_0_trackUE->Clone("pPb5Pbp5_reco_trackUE");
  pPb5Pbp5_0_track_xi = (TH2D*)pPb5_0_track_xi->Clone("pPb5Pbp5_reco_track_xi");
  pPb5Pbp5_0_trackUE_xi = (TH2D*)pPb5_0_trackUE_xi->Clone("pPb5Pbp5_reco_trackUE_xi");
  pPb5Pbp5_0_jet->Add(Pbp5_0_jet);
  pPb5Pbp5_0_track->Add(Pbp5_0_track);
  pPb5Pbp5_0_trackUE->Add(Pbp5_0_trackUE);
  pPb5Pbp5_0_track_xi->Add(Pbp5_0_track_xi);
  pPb5Pbp5_0_trackUE_xi->Add(Pbp5_0_trackUE_xi);

  pPb5_1_jet_reco = (TH1D*) MCFilepPb5->Get(Form("pPb5_reco_jet%s",pPb5Tag.data()));
  pPb5_1_track_reco = (TH2D*) MCFilepPb5->Get(Form("pPb5_reco_track%s",pPb5Tag.data()));
  pPb5_1_trackUE_reco = (TH2D*) MCFilepPb5->Get(Form("pPb5_reco_trackUE%s",pPb5Tag.data()));
  pPb5_1_track_xi_reco = (TH2D*) MCFilepPb5->Get(Form("pPb5_reco_track_xi%s",pPb5Tag.data()));
  pPb5_1_trackUE_xi_reco = (TH2D*)MCFilepPb5->Get(Form("pPb5_reco_trackUE_xi%s",pPb5Tag.data()));


  pPb5_1_track_rJgT = (TH2D*) MCFilepPb5->Get(Form("pPb5_rJgT_track%s",pPb5Tag.data()));
  pPb5_1_trackUE_rJgT = (TH2D*) MCFilepPb5->Get(Form("pPb5_rJgT_trackUE%s",pPb5Tag.data()));
  pPb5_1_track_xi_rJgT = (TH2D*) MCFilepPb5->Get(Form("pPb5_rJgT_track_xi%s",pPb5Tag.data()));
  pPb5_1_trackUE_xi_rJgT = (TH2D*)MCFilepPb5->Get(Form("pPb5_rJgT_trackUE_xi%s",pPb5Tag.data()));
  pPb5_1_track_gJrT = (TH2D*) MCFilepPb5->Get(Form("pPb5_gJrT_track%s",pPb5Tag.data()));
  pPb5_1_trackUE_gJrT = (TH2D*) MCFilepPb5->Get(Form("pPb5_gJrT_trackUE%s",pPb5Tag.data()));
  pPb5_1_track_xi_gJrT = (TH2D*) MCFilepPb5->Get(Form("pPb5_gJrT_track_xi%s",pPb5Tag.data()));
  pPb5_1_trackUE_xi_gJrT = (TH2D*)MCFilepPb5->Get(Form("pPb5_gJrT_trackUE_xi%s",pPb5Tag.data()));

  Pbp5_1_jet_reco = (TH1D*) MCFilePbp5->Get(Form("Pbp5_reco_jet%s",pPb5Tag.data()));
  Pbp5_1_track_reco = (TH2D*) MCFilePbp5->Get(Form("Pbp5_reco_track%s",pPb5Tag.data()));
  Pbp5_1_trackUE_reco = (TH2D*) MCFilePbp5->Get(Form("Pbp5_reco_trackUE%s",pPb5Tag.data()));
  Pbp5_1_track_xi_reco = (TH2D*) MCFilePbp5->Get(Form("Pbp5_reco_track_xi%s",pPb5Tag.data()));
  Pbp5_1_trackUE_xi_reco = (TH2D*)MCFilePbp5->Get(Form("Pbp5_reco_trackUE_xi%s",pPb5Tag.data()));


  Pbp5_1_track_rJgT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_rJgT_track%s",pPb5Tag.data()));
  Pbp5_1_trackUE_rJgT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_rJgT_trackUE%s",pPb5Tag.data()));
  Pbp5_1_track_xi_rJgT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_rJgT_track_xi%s",pPb5Tag.data()));
  Pbp5_1_trackUE_xi_rJgT = (TH2D*)MCFilePbp5->Get(Form("Pbp5_rJgT_trackUE_xi%s",pPb5Tag.data()));
  Pbp5_1_track_gJrT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gJrT_track%s",pPb5Tag.data()));
  Pbp5_1_trackUE_gJrT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gJrT_trackUE%s",pPb5Tag.data()));
  Pbp5_1_track_xi_gJrT = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gJrT_track_xi%s",pPb5Tag.data()));
  Pbp5_1_trackUE_xi_gJrT = (TH2D*)MCFilePbp5->Get(Form("Pbp5_gJrT_trackUE_xi%s",pPb5Tag.data()));

  /*
  pp5_1_jet_reco = (TH1D*) MCFilepp5->Get(Form("pp5_reco_jet%s",pPb5Tag.data()));
  pp5_1_track_reco = (TH2D*) MCFilepp5->Get(Form("pp5_reco_track%s",pPb5Tag.data()));
  pp5_1_trackUE_reco = (TH2D*) MCFilepp5->Get(Form("pp5_reco_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_reco = (TH2D*) MCFilepp5->Get(Form("pp5_reco_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_reco = (TH2D*) MCFilepp5->Get(Form("pp5_reco_trackUE_xi%s",pPb5Tag.data()));


  pp5_1_track_rJgT = (TH2D*) MCFilepp5->Get(Form("pp5_rJgT_track%s",pPb5Tag.data()));
  pp5_1_trackUE_rJgT = (TH2D*) MCFilepp5->Get(Form("pp5_rJgT_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_rJgT = (TH2D*) MCFilepp5->Get(Form("pp5_rJgT_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_rJgT = (TH2D*) MCFilepp5->Get(Form("pp5_rJgT_trackUE_xi%s",pPb5Tag.data()));
  pp5_1_track_gJrT = (TH2D*) MCFilepp5->Get(Form("pp5_gJrT_track%s",pPb5Tag.data()));
  pp5_1_trackUE_gJrT = (TH2D*) MCFilepp5->Get(Form("pp5_gJrT_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_gJrT = (TH2D*) MCFilepp5->Get(Form("pp5_gJrT_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_gJrT = (TH2D*) MCFilepp5->Get(Form("pp5_gJrT_trackUE_xi%s",pPb5Tag.data()));
  */
  pp5_1_jet_reco = (TH1D*) MCFilepp5->Get(Form("ppref5_reco_jet%s",pPb5Tag.data()));
  pp5_1_track_reco = (TH2D*) MCFilepp5->Get(Form("ppref5_reco_track%s",pPb5Tag.data()));
  pp5_1_trackUE_reco = (TH2D*) MCFilepp5->Get(Form("ppref5_reco_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_reco = (TH2D*) MCFilepp5->Get(Form("ppref5_reco_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_reco = (TH2D*) MCFilepp5->Get(Form("ppref5_reco_trackUE_xi%s",pPb5Tag.data()));


  pp5_1_track_rJgT = (TH2D*) MCFilepp5->Get(Form("ppref5_rJgT_track%s",pPb5Tag.data()));
  pp5_1_trackUE_rJgT = (TH2D*) MCFilepp5->Get(Form("ppref5_rJgT_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_rJgT = (TH2D*) MCFilepp5->Get(Form("ppref5_rJgT_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_rJgT = (TH2D*) MCFilepp5->Get(Form("ppref5_rJgT_trackUE_xi%s",pPb5Tag.data()));
  pp5_1_track_gJrT = (TH2D*) MCFilepp5->Get(Form("ppref5_gJrT_track%s",pPb5Tag.data()));
  pp5_1_trackUE_gJrT = (TH2D*) MCFilepp5->Get(Form("ppref5_gJrT_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_gJrT = (TH2D*) MCFilepp5->Get(Form("ppref5_gJrT_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_gJrT = (TH2D*) MCFilepp5->Get(Form("ppref5_gJrT_trackUE_xi%s",pPb5Tag.data()));


  pPb5_1_jet_gen = (TH1D*) MCFilepPb5->Get(Form("pPb5_gen_jet%s",pPb5Tag.data()));
  pPb5_1_track_gen = (TH2D*) MCFilepPb5->Get(Form("pPb5_gen_track%s",pPb5Tag.data()));
  pPb5_1_trackUE_gen = (TH2D*) MCFilepPb5->Get(Form("pPb5_gen_trackUE%s",pPb5Tag.data()));
  pPb5_1_track_xi_gen = (TH2D*) MCFilepPb5->Get(Form("pPb5_gen_track_xi%s",pPb5Tag.data()));
  pPb5_1_trackUE_xi_gen = (TH2D*)MCFilepPb5->Get(Form("pPb5_gen_trackUE_xi%s",pPb5Tag.data()));
 

  Pbp5_1_jet_gen = (TH1D*) MCFilePbp5->Get(Form("Pbp5_gen_jet%s",pPb5Tag.data()));
  Pbp5_1_track_gen = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gen_track%s",pPb5Tag.data()));
  Pbp5_1_trackUE_gen = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gen_trackUE%s",pPb5Tag.data()));
  Pbp5_1_track_xi_gen = (TH2D*) MCFilePbp5->Get(Form("Pbp5_gen_track_xi%s",pPb5Tag.data()));
  Pbp5_1_trackUE_xi_gen = (TH2D*)MCFilePbp5->Get(Form("Pbp5_gen_trackUE_xi%s",pPb5Tag.data()));
 
  /*
  pp5_1_jet_gen = (TH1D*) MCFilepp5->Get(Form("pp5_gen_jet%s",pPb5Tag.data()));
  pp5_1_track_gen = (TH2D*) MCFilepp5->Get(Form("pp5_gen_track%s",pPb5Tag.data()));
  pp5_1_trackUE_gen = (TH2D*) MCFilepp5->Get(Form("pp5_gen_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_gen = (TH2D*) MCFilepp5->Get(Form("pp5_gen_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_gen = (TH2D*) MCFilepp5->Get(Form("pp5_gen_trackUE_xi%s",pPb5Tag.data()));
  */
  pp5_1_jet_gen = (TH1D*) MCFilepp5->Get(Form("ppref5_gen_jet%s",pPb5Tag.data()));
  pp5_1_track_gen = (TH2D*) MCFilepp5->Get(Form("ppref5_gen_track%s",pPb5Tag.data()));
  pp5_1_trackUE_gen = (TH2D*) MCFilepp5->Get(Form("ppref5_gen_trackUE%s",pPb5Tag.data()));
  pp5_1_track_xi_gen = (TH2D*) MCFilepp5->Get(Form("ppref5_gen_track_xi%s",pPb5Tag.data()));
  pp5_1_trackUE_xi_gen = (TH2D*) MCFilepp5->Get(Form("ppref5_gen_trackUE_xi%s",pPb5Tag.data()));

  
//reweighting combined MC to match data in terms of pPb vs Pbp fraction
  Njets_pPbData = (TH1D*) spectraFilepPb5->Get("totalJetsPtCutHist"); 
  Njets_pPbData->SetName("totalJetsPtCutHistpPbdata"); 
  Njets_pPbMC = (TH1D*)   MCFilepPb5->Get("totalJetsPtCutHist"); 
  Njets_pPbMC->SetName("totalJetsPtCutHistpPbMC"); 
  Njets_PbpData = (TH1D*) spectraFilePbp5->Get("totalJetsPtCutHist");
  Njets_PbpData->SetName("totalJetsPtCutHistPbpdata"); 
  Njets_PbpMC = (TH1D*)   MCFilePbp5->Get("totalJetsPtCutHist");
  Njets_PbpMC->SetName("totalJetsPtCutHistPbpMC"); 

  double pPbMCFracCorr = Njets_pPbData->GetBinContent(2)*(Njets_pPbMC->GetBinContent(2)+Njets_PbpMC->GetBinContent(2))/(Njets_pPbMC->GetBinContent(2)*(Njets_pPbData->GetBinContent(2)+Njets_PbpData->GetBinContent(2)));
  double PbpMCFracCorr = Njets_PbpData->GetBinContent(2)*(Njets_pPbMC->GetBinContent(2)+Njets_PbpMC->GetBinContent(2))/(Njets_PbpMC->GetBinContent(2)*(Njets_pPbData->GetBinContent(2)+Njets_PbpData->GetBinContent(2)));
  std::cout << "pPb MC correction factor for combining w/ Pbp:" << pPbMCFracCorr << std::endl;
  std::cout << "Pbp MC correction factor for combining w/ pPb:" << PbpMCFracCorr << std::endl;

  TH1D* pPb = (TH1D*) pPb5_1_jet_reco->Clone("pPbRecoMC_jet");
  TH1D* Pbp = (TH1D*) Pbp5_1_jet_reco->Clone("PbpRecoMC_jet");
  pPb->Scale(pPbMCFracCorr);
  Pbp->Scale(PbpMCFracCorr);
  pPb->Add(Pbp);
  pPb5Pbp5_1_jet_reco = (TH1D*) pPb->Clone(Form("pPb5Pbp5_reco_jet%s",pPb5Tag.data()));

  TH2D* pPb2 = (TH2D*) pPb5_1_track_reco->Clone("pPbRecoMC_track");
  TH2D* Pbp2 = (TH2D*) Pbp5_1_track_reco->Clone("PbpRecoMC_track");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_reco = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_reco_track%s",pPb5Tag.data()));
 
  pPb2 = (TH2D*) pPb5_1_trackUE_reco->Clone("pPbRecoMC_trackUE");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_reco->Clone("PbpRecoMC_trackUE");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_reco = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_reco_trackUE%s",pPb5Tag.data()));
  
  pPb2 = (TH2D*) pPb5_1_track_xi_reco->Clone("pPbRecoMC_track_xi");
  Pbp2 = (TH2D*) Pbp5_1_track_xi_reco->Clone("PbpRecoMC_track_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_xi_reco = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_reco_track_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_trackUE_xi_reco->Clone("pPbRecoMC_trackUE_xi");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_xi_reco->Clone("PbpRecoMC_trackUE_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_xi_reco = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_reco_trackUE_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_track_rJgT->Clone("pPbrJgTMC_track");
  Pbp2 = (TH2D*) Pbp5_1_track_rJgT->Clone("PbprJgTMC_track");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_rJgT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_rJgT_track%s",pPb5Tag.data()));
 
  pPb2 = (TH2D*) pPb5_1_trackUE_rJgT->Clone("pPbrJgTMC_trackUE");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_rJgT->Clone("PbprJgTMC_trackUE");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_rJgT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_rJgT_trackUE%s",pPb5Tag.data()));
  
  pPb2 = (TH2D*) pPb5_1_track_xi_rJgT->Clone("pPbrJgTMC_track_xi");
  Pbp2 = (TH2D*) Pbp5_1_track_xi_rJgT->Clone("PbprJgTMC_track_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_xi_rJgT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_rJgT_track_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_trackUE_xi_rJgT->Clone("pPbrJgTMC_trackUE_xi");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_xi_rJgT->Clone("PbprJgTMC_trackUE_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_xi_rJgT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_rJgT_trackUE_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_track_gJrT->Clone("pPbgJrTMC_track");
  Pbp2 = (TH2D*) Pbp5_1_track_gJrT->Clone("PbpgJrTMC_track");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_gJrT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gJrT_track%s",pPb5Tag.data()));
 
  pPb2 = (TH2D*) pPb5_1_trackUE_gJrT->Clone("pPbgJrTMC_trackUE");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_gJrT->Clone("PbpgJrTMC_trackUE");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_gJrT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gJrT_trackUE%s",pPb5Tag.data()));
  
  pPb2 = (TH2D*) pPb5_1_track_xi_gJrT->Clone("pPbgJrTMC_track_xi");
  Pbp2 = (TH2D*) Pbp5_1_track_xi_gJrT->Clone("PbpgJrTMC_track_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_xi_gJrT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gJrT_track_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_trackUE_xi_gJrT->Clone("pPbgJrTMC_trackUE_xi");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_xi_gJrT->Clone("PbpgJrTMC_trackUE_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_xi_gJrT = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gJrT_trackUE_xi%s",pPb5Tag.data()));

  pPb = (TH1D*) pPb5_1_jet_gen->Clone("pPbGenMC_jet");
  Pbp = (TH1D*) Pbp5_1_jet_gen->Clone("PbpGenMC_jet");
  pPb->Scale(pPbMCFracCorr);
  Pbp->Scale(PbpMCFracCorr);
  pPb->Add(Pbp);
  pPb5Pbp5_1_jet_gen = (TH1D*) pPb->Clone(Form("pPb5Pbp5_gen_jet%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_track_gen->Clone("pPbGenMC_track");
  Pbp2 = (TH2D*) Pbp5_1_track_gen->Clone("PbpGenMC_track");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_gen = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gen_track%s",pPb5Tag.data()));
 
  pPb2 = (TH2D*) pPb5_1_trackUE_gen->Clone("pPbGenMC_trackUE");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_gen->Clone("PbpGenMC_trackUE");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_gen = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gen_trackUE%s",pPb5Tag.data()));
  
  pPb2 = (TH2D*) pPb5_1_track_xi_gen->Clone("pPbGenMC_track_xi");
  Pbp2 = (TH2D*) Pbp5_1_track_xi_gen->Clone("PbpGenMC_track_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_track_xi_gen = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gen_track_xi%s",pPb5Tag.data()));

  pPb2 = (TH2D*) pPb5_1_trackUE_xi_gen->Clone("pPbGenMC_trackUE_xi");
  Pbp2 = (TH2D*) Pbp5_1_trackUE_xi_gen->Clone("PbpGenMC_trackUE_xi");
  pPb2->Scale(pPbMCFracCorr);
  Pbp2->Scale(PbpMCFracCorr);
  pPb2->Add(Pbp2);
  pPb5Pbp5_1_trackUE_xi_gen = (TH2D*) pPb2->Clone(Form("pPb5Pbp5_gen_trackUE_xi%s",pPb5Tag.data()));
}
