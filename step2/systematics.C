#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TPad.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TBox.h"
#include "TF1.h"
#include "TMath.h"

/*const int variations = 20;
const char * variationTag[variations]= {"","_pp2JESUP4","_pp2JESDOWN4","_pp7JESUP4","_pp7JESDOWN4","_pPb5JESUP4","_pPb5JESDOWN4","_pp2JER10","_pp7JER10","_pPb5JER10","_pp2JER2","_pp7JER2","_pPb5JER2","_NoTrackCorr","_pp2JESUP1","_pp2JESDOWN1","_pp7JESUP1","_pp7JESDOWN1","_pPb5JESUP1","_pPb5JESDOWN1"};

const int FF_Bins = 5;
double FF_Bound[FF_Bins+1] = {60,80,100,120,140,200};*/


TH1D** getRatio(const char * mode = "pp2", int v=0, int UEtype1=3, int UEtype2=3)
{
  TH1D** ratioArray = new TH1D*[FF_Bins*2];

  TFile * inf1 = TFile::Open(Form("FragmentationFunctions%sUE%d.root",variationTag[v],UEtype1),"read");
  TFile * inf2 = TFile::Open(Form("FragmentationFunctionsUE%d.root",UEtype2),"read");
  for(int i =0; i<FF_Bins*2; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    if(strcmp(mode,"pPb5")==0)
    {
      ratioArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      ratioArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      ratioArray[i]->SetName(Form("%s%s_Ratio%d%s",mode,variationTag[v],i,isXi.data()));
      ratioArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pp2")==0)
    {
      ratioArray[i] = (TH1D*)inf1->Get(Form("pp2TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      ratioArray[i]->Divide((TH1D*)inf2->Get(Form("pp2TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      ratioArray[i]->SetName(Form("%s%s_Ratio%d%s",mode,variationTag[v],i,isXi.data()));
      ratioArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pp7")==0)
    {
      ratioArray[i] = (TH1D*)inf1->Get(Form("pp7TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      ratioArray[i]->Divide((TH1D*)inf2->Get(Form("pp7TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      ratioArray[i]->SetName(Form("%s%s_Ratio%d%s",mode,variationTag[v],i,isXi.data()));
      ratioArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"interp")==0)
    {
      ratioArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      ratioArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      ratioArray[i]->SetName(Form("%s%s_Ratio%d%s",mode,variationTag[v],i,isXi.data()));
      ratioArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"FFratio")==0)
    {
      ratioArray[i] = (TH1D*)inf1->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      ratioArray[i]->Divide((TH1D*)inf2->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      ratioArray[i]->SetName(Form("%s%s_Ratio%d%s",mode,variationTag[v],i,isXi.data()));
      ratioArray[i]->SetDirectory(0);
    }
  }  
  inf1->Close();
  inf2->Close();
  return ratioArray;
}

TH1D** getpPbPbpDiff(const char * mode = "pPb5", int v=0, int UEtype=3)
{
  TH1D** diffArray = new TH1D*[FF_Bins*2];

  TFile * inf1 = TFile::Open(Form("FragmentationFunctions%sUE%d.root",variationTag[v],UEtype),"read");
  TFile * inf2 = TFile::Open(Form("FragmentationFunctionsUE%d.root",UEtype),"read");
  for(int i =0; i<FF_Bins*2; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    if(strcmp(mode,"pPb5")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("Pbp5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->Scale(0.5);
      diffArray[i]->SetName(Form("%s%s_pPbPbpDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"FFratio")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("Pbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->Scale(0.5);
      diffArray[i]->SetName(Form("%s%s_pPbPbpDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"interp")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->Scale(0.5);
      diffArray[i]->SetName(Form("%s%s_pPbPbpDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
  }  
  inf1->Close();
  inf2->Close();
  return diffArray;
}

TH1D** getMCDiff(const char * mode = "pPb5", int v=0,int UEtype=3)
{
  TH1D** diffArray = new TH1D*[FF_Bins*2];

  TFile * inf1 = TFile::Open(Form("FragmentationFunctions%sUE%d.root",variationTag[v],UEtype),"read");
  TFile * inf2 = TFile::Open(Form("FragmentationFunctionsUE%d.root",UEtype),"read");
  for(int i =0; i<FF_Bins*2; i++)
  {   
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    if(strcmp(mode,"pp2")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pp2TeV_reco_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("pp2TeV_gen_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pp2TeV_reco_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s%s_MCDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pp7")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pp7TeV_reco_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("pp7TeV_gen_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pp7TeV_reco_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s%s_MCDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pPb5")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("pPb5Pbp5TeV_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s%s_MCDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"interp")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_recoMC_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("pPb5Pbp5TeV_genMC_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_recoMC_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s%s_MCDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"FFratio")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPbPbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf1->Get(Form("pPbPbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPbPbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s%s_MCDiff%d%s",mode,variationTag[v],i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
  }  
  inf1->Close();
  inf2->Close();
  return diffArray;
}

TH1D** getChargeCutDiff(const char * mode = "pPb5", int UEtype=3)
{
  TH1D** diffArray = new TH1D*[FF_Bins*2];

  TFile * inf1 = TFile::Open(Form("FragmentationFunctions%sUE%d.root",variationTag[30],UEtype),"read");
  TFile * inf2 = TFile::Open(Form("FragmentationFunctionsUE%d.root",UEtype),"read");
  for(int i =0; i<FF_Bins*2; i++)
  {   
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    if(strcmp(mode,"pp2")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pp2TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf2->Get(Form("pp2TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pp2TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s_ChargeCutDiff%d%s",mode,i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pp7")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pp7TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf2->Get(Form("pp7TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pp7TeV_NoReweight_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s_ChargeCutDiff%d%s",mode,i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"pPb5")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s_ChargeCutDiff%d%s",mode,i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"interp")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPb5Pbp5TeV_data_interp_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s_ChargeCutDiff%d%s",mode,i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
    if(strcmp(mode,"FFratio")==0)
    {
      diffArray[i] = (TH1D*)inf1->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
      diffArray[i]->Add((TH1D*)inf2->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())),-1);
      diffArray[i]->Divide((TH1D*)inf2->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())));
      diffArray[i]->SetName(Form("%s_ChargeCutDiff%d%s",mode,i,isXi.data()));
      diffArray[i]->SetDirectory(0);
    }
  }  
  inf1->Close();
  inf2->Close();
  return diffArray;
}

void FFSystematics(const char * mode, int UEtype)
{
  //getting sources of error as ratios
   TH1D** JESUP;
   TH1D** JESDOWN;
   TH1D** JER;
   TH1D** pPbPbpDiff;
   TH1D** MCDiff;
   TH1D** ChargeCutDiff;
   TH1D** UEDiff;
   //estimating 5% tracking errors
   double track = 0.039;

  if(strcmp(mode,"pPb5")==0)
  {
    JESUP = getRatio(mode,5,UEtype,UEtype);
    JESDOWN = getRatio(mode,6,UEtype,UEtype);
    JER = getRatio(mode,9,UEtype,UEtype);
    pPbPbpDiff = getpPbPbpDiff(mode,0,UEtype);
    MCDiff = getMCDiff(mode,30,UEtype);
    ChargeCutDiff = getChargeCutDiff(mode,UEtype);
//    if(UEtype==0) UEDiff = getRatio(mode,26,UEtype);
    if(UEtype==0) UEDiff = getRatio(mode,0,2,UEtype);
  }
  if(strcmp(mode,"pp2")==0)
  {
    JESUP = getRatio(mode,1,UEtype,UEtype);
    JESDOWN = getRatio(mode,2,UEtype,UEtype);
    JER = getRatio(mode,7,UEtype,UEtype);
    MCDiff = getMCDiff(mode,30,UEtype); 
    ChargeCutDiff = getChargeCutDiff(mode,UEtype);
    //if(UEtype==0) UEDiff = getRatio(mode,26,UEtype);
    if(UEtype==0) UEDiff = getRatio(mode,0,2,UEtype);
  }
  if(strcmp(mode,"pp7")==0)
  {
    JESUP = getRatio(mode,3,UEtype,UEtype);
    JESDOWN = getRatio(mode,4,UEtype,UEtype);
    JER = getRatio(mode,8,UEtype,UEtype);
    MCDiff = getMCDiff(mode,30,UEtype);
    ChargeCutDiff = getChargeCutDiff(mode,UEtype);
    //if(UEtype==0) UEDiff = getRatio(mode,26,UEtype);
    if(UEtype==0) UEDiff = getRatio(mode,0,2,UEtype); 
  }

  TFile * output;
  if(strcmp("pPb5",mode)==0) output = new TFile(Form("SystematicsUE%d.root",UEtype),"recreate");
  else output = new TFile(Form("SystematicsUE%d.root",UEtype),"update");

  //doing asymmetric JES errors
  TH1D* JESTotUP[FF_Bins*2];
  TH1D* JESTotDOWN[FF_Bins*2];
  TH1D* JERTot[FF_Bins*2];
  for(int i = 0; i<FF_Bins*2; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    JESTotUP[i] = (TH1D*)JESUP[i]->Clone(Form("%s_JESTotUP%d%s",mode,i,isXi.data()));
    JESTotDOWN[i] = (TH1D*)JESDOWN[i]->Clone(Form("%s_JESTotDOWN%d%s",mode,i,isXi.data()));
    JERTot[i] = (TH1D*)JER[i]->Clone(Form("%s_JERTot%d%s",mode,i,isXi.data()));

    for(int j = 0; j<JESTotUP[0]->GetSize(); j++)
    {
      if(JESTotUP[i]->GetBinContent(j)<JESDOWN[i]->GetBinContent(j)) JESTotUP[i]->SetBinContent(j,JESDOWN[i]->GetBinContent(j));
      JESTotUP[i]->SetBinContent(j,JESTotUP[i]->GetBinContent(j)-1);
      if(JESTotUP[i]->GetBinContent(j)<0.001) JESTotUP[i]->SetBinContent(j,0.001);

      if(JESTotDOWN[i]->GetBinContent(j)>JESUP[i]->GetBinContent(j)) JESTotDOWN[i]->SetBinContent(j,JESUP[i]->GetBinContent(j));
      JESTotDOWN[i]->SetBinContent(j,JESTotDOWN[i]->GetBinContent(j)-1);
      if(JESTotDOWN[i]->GetBinContent(j)>-0.001) JESTotDOWN[i]->SetBinContent(j,-0.001);

      if(JERTot[i]->GetBinContent(j)!=0)  JERTot[i]->SetBinContent(j,TMath::Abs(JERTot[i]->GetBinContent(j)-1));
    }
  }

  //adding all sources together in quadrature
  TH1D* TotUP[FF_Bins*2];
  TH1D* TotDOWN[FF_Bins*2];
  for(int i = 0; i<2*FF_Bins; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    TotUP[i] = (TH1D*)JESTotUP[i]->Clone(Form("%s_TotUP%d%s",mode,i,isXi.data())); 
    TotDOWN[i] = (TH1D*)JESTotDOWN[i]->Clone(Form("%s_TotDOWN%d%s",mode,i,isXi.data()));
  
    for(int j = 0; j<TotUP[0]->GetSize(); j++)
    {
      //if empty bins
      if(JER[i]->GetBinContent(j)==0 || JESTotUP[i]->GetBinContent(j)==-1 || JESTotDOWN[i]->GetBinContent(j)==-1)
      {
        TotUP[i]->SetBinContent(j,0);
        TotDOWN[i]->SetBinContent(j,0);
        TotUP[i]->SetBinError(j,0);
        TotDOWN[i]->SetBinError(j,0);
        JESTotUP[i]->SetBinContent(j,0);
        JESTotDOWN[i]->SetBinContent(j,0);
        JESTotUP[i]->SetBinError(j,0);
        JESTotDOWN[i]->SetBinError(j,0);
      }
      else 
      {
        TotUP[i]->SetBinContent(j,TMath::Power(TotUP[i]->GetBinContent(j),2)+TMath::Power(JER[i]->GetBinContent(j)-1,2)+track*track);
        TotDOWN[i]->SetBinContent(j,TMath::Power(TotDOWN[i]->GetBinContent(j),2)+TMath::Power(JER[i]->GetBinContent(j)-1,2)+track*track);
        if(strcmp(mode,"pPb5")==0)
        {
          TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(pPbPbpDiff[i]->GetBinContent(j),2));
          TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(pPbPbpDiff[i]->GetBinContent(j),2));
        }
        if(MCDiff[i]->GetBinContent(j)<0) TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(MCDiff[i]->GetBinContent(j),2));
        else if(MCDiff[i]->GetBinContent(j)>0) TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(MCDiff[i]->GetBinContent(j),2));
        if(ChargeCutDiff[i]->GetBinContent(j)>0) TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(ChargeCutDiff[i]->GetBinContent(j),2));
        else if(ChargeCutDiff[i]->GetBinContent(j)<0) TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(ChargeCutDiff[i]->GetBinContent(j),2));

        if(UEtype==0)
        {        
          TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+ TMath::Power(UEDiff[i]->GetBinContent(j)-1,2));
          TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+ TMath::Power(UEDiff[i]->GetBinContent(j)-1,2));
          UEDiff[i]->SetBinContent(j,TMath::Abs(UEDiff[i]->GetBinContent(j)-1));
        }      
 
        TotUP[i]->SetBinContent(j,TMath::Power(TotUP[i]->GetBinContent(j),0.5));
        TotDOWN[i]->SetBinContent(j,TMath::Power(TotDOWN[i]->GetBinContent(j),0.5));
      }
    } 
    JESUP[i]->Write();
    JESDOWN[i]->Write();
    JER[i]->Write();
    if(strcmp(mode,"pPb5")==0) pPbPbpDiff[i]->Write();
    MCDiff[i]->Write();
    ChargeCutDiff[i]->Write();
    if(UEtype==0) UEDiff[i]->Write();
    JESTotUP[i]->Write();
    JESTotDOWN[i]->Write();
    JERTot[i]->Write();
    TotUP[i]->Write();
    TotDOWN[i]->Write();
  }
  output->Close();
}

void Interpolation_and_Ratio_Systematics(const char * mode = "interp", int UEtype=3)
{
  //getting sources of error as ratios
  TH1D** pp2JESUP;
  TH1D** pp2JESDOWN;
  TH1D** pp2JER; 
  TH1D** pp7JESUP;
  TH1D** pp7JESDOWN;
  TH1D** pp7JER;
  TH1D** pPb5JESUP;
  TH1D** pPb5JESDOWN;
  TH1D** pPb5JER;
  TH1D** pPbPbpDiff;
  TH1D** MCDiff;
  TH1D** ChargeCutDiff;
  TH1D** UEDiff;

  //estimating 5% tracking errors
  double track = 0.039;
  double pdf_uncert = 0.005;
 
  if(strcmp(mode,"interp")==0) 
  {
    pp2JESUP = getRatio(mode,1,UEtype,UEtype);
    pp2JESDOWN = getRatio(mode,2,UEtype,UEtype);
    pp2JER = getRatio(mode,7,UEtype,UEtype);
    pp7JESUP = getRatio(mode,3,UEtype,UEtype);
    pp7JESDOWN = getRatio(mode,4,UEtype,UEtype);
    pp7JER = getRatio(mode,8,UEtype,UEtype);
    pPb5JESUP = getRatio(mode,5,UEtype,UEtype);
    pPb5JESDOWN = getRatio(mode,6,UEtype,UEtype);
    pPb5JER = getRatio(mode,9,UEtype,UEtype);
    pPbPbpDiff = getpPbPbpDiff(mode,0,UEtype);
    MCDiff = getMCDiff(mode,30,UEtype);
    ChargeCutDiff = getChargeCutDiff(mode,UEtype);
    //if(UEtype==0) UEDiff = getRatio(mode,26,UEtype);
    if(UEtype==0) UEDiff = getRatio(mode,0,2,UEtype);
  }
  else if(strcmp(mode,"FFratio")==0) 
  {
    track = 0.05515432893; //3.9% times sqrt(2)
    pp2JESUP = getRatio(mode,20,UEtype,UEtype);
    pp2JESDOWN = getRatio(mode,21,UEtype,UEtype);
    pp2JER = getRatio(mode,7,UEtype,UEtype);
    pp7JESUP = getRatio(mode,22,UEtype,UEtype);
    pp7JESDOWN = getRatio(mode,23,UEtype,UEtype);
    pp7JER = getRatio(mode,8,UEtype,UEtype);
    pPb5JESUP = getRatio(mode,24,UEtype,UEtype);
    pPb5JESDOWN = getRatio(mode,25,UEtype,UEtype);
    pPb5JER = getRatio(mode,9,UEtype,UEtype); 
    pPbPbpDiff = getpPbPbpDiff(mode,0,UEtype);
    MCDiff = getMCDiff(mode,30,UEtype);
    ChargeCutDiff = getChargeCutDiff(mode,UEtype);
    //if(UEtype==0) UEDiff = getRatio(mode,26,UEtype);
    if(UEtype==0) UEDiff = getRatio(mode,0,2,UEtype);
  }

  TFile * output = new TFile(Form("SystematicsUE%d.root",UEtype),"update");

  TH1D* pp2JESTotUP[FF_Bins*2];
  TH1D* pp2JESTotDOWN[FF_Bins*2];
  TH1D* pp7JESTotUP[FF_Bins*2];
  TH1D* pp7JESTotDOWN[FF_Bins*2];
  TH1D* pPb5JESTotUP[FF_Bins*2];
  TH1D* pPb5JESTotDOWN[FF_Bins*2];

  for(int i = 0; i<2*FF_Bins; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    pp2JESTotUP[i] = (TH1D*)pp2JESUP[i]->Clone(Form("%s_pp2JESTotUP%d%s",mode,i,isXi.data()));
    pp2JESTotDOWN[i] = (TH1D*)pp2JESDOWN[i]->Clone(Form("%s_pp2JESTotDOWN%d%s",mode,i,isXi.data()));
    pp7JESTotUP[i] = (TH1D*)pp7JESUP[i]->Clone(Form("%s_pp7JESTotUP%d%s",mode,i,isXi.data()));
    pp7JESTotDOWN[i] = (TH1D*)pp7JESDOWN[i]->Clone(Form("%s_pp7JESTotDOWN%d%s",mode,i,isXi.data()));
    pPb5JESTotUP[i] = (TH1D*)pPb5JESUP[i]->Clone(Form("%s_pPb5JESTotUP%d%s",mode,i,isXi.data()));
    pPb5JESTotDOWN[i] = (TH1D*)pPb5JESDOWN[i]->Clone(Form("%s_pPb5JESTotDOWN%d%s",mode,i,isXi.data()));

    for(int j = 0; j<pp2JESTotUP[0]->GetSize(); j++)
    {
      if(pp2JESTotUP[i]->GetBinContent(j)<pp2JESDOWN[i]->GetBinContent(j)) pp2JESTotUP[i]->SetBinContent(j,pp2JESDOWN[i]->GetBinContent(j));
      pp2JESTotUP[i]->SetBinContent(j,pp2JESTotUP[i]->GetBinContent(j)-1);
      if(pp2JESTotUP[i]->GetBinContent(j)<0.001) pp2JESTotUP[i]->SetBinContent(j,0.001);

      if(pp2JESTotDOWN[i]->GetBinContent(j)>pp2JESUP[i]->GetBinContent(j)) pp2JESTotDOWN[i]->SetBinContent(j,pp2JESUP[i]->GetBinContent(j));
      pp2JESTotDOWN[i]->SetBinContent(j,pp2JESTotDOWN[i]->GetBinContent(j)-1);
      if(pp2JESTotDOWN[i]->GetBinContent(j)>-0.001) pp2JESTotDOWN[i]->SetBinContent(j,-0.001);
 
      if(pp7JESTotUP[i]->GetBinContent(j)<pp7JESDOWN[i]->GetBinContent(j)) pp7JESTotUP[i]->SetBinContent(j,pp7JESDOWN[i]->GetBinContent(j));
      pp7JESTotUP[i]->SetBinContent(j,pp7JESTotUP[i]->GetBinContent(j)-1);
      if(pp7JESTotUP[i]->GetBinContent(j)<0.001) pp7JESTotUP[i]->SetBinContent(j,0.001);

      if(pp7JESTotDOWN[i]->GetBinContent(j)>pp7JESUP[i]->GetBinContent(j)) pp7JESTotDOWN[i]->SetBinContent(j,pp7JESUP[i]->GetBinContent(j));
      pp7JESTotDOWN[i]->SetBinContent(j,pp7JESTotDOWN[i]->GetBinContent(j)-1);   
      if(pp7JESTotDOWN[i]->GetBinContent(j)>-0.001) pp7JESTotDOWN[i]->SetBinContent(j,-0.001);

      if(pPb5JESTotUP[i]->GetBinContent(j)<pPb5JESDOWN[i]->GetBinContent(j)) pPb5JESTotUP[i]->SetBinContent(j,pPb5JESDOWN[i]->GetBinContent(j));
      pPb5JESTotUP[i]->SetBinContent(j,pPb5JESTotUP[i]->GetBinContent(j)-1);
      if(pPb5JESTotUP[i]->GetBinContent(j)<0.001) pPb5JESTotUP[i]->SetBinContent(j,0.001);

      if(pPb5JESTotDOWN[i]->GetBinContent(j)>pPb5JESUP[i]->GetBinContent(j)) pPb5JESTotDOWN[i]->SetBinContent(j,pPb5JESUP[i]->GetBinContent(j));
      pPb5JESTotDOWN[i]->SetBinContent(j,pPb5JESTotDOWN[i]->GetBinContent(j)-1);
      if(pPb5JESTotDOWN[i]->GetBinContent(j)>-0.001) pPb5JESTotDOWN[i]->SetBinContent(j,-0.001);
    }
  }

  //adding all sources in quadrature
  TH1D* TotUP[FF_Bins*2];
  TH1D* TotDOWN[FF_Bins*2];
  TH1D* JESTotUP[FF_Bins*2];
  TH1D* JESTotDOWN[FF_Bins*2];
  TH1D* JERTot[FF_Bins*2];

  for(int i = 0; i<2*FF_Bins; i++)
  {
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";
    TotUP[i] = (TH1D*)pp2JESTotUP[i]->Clone(Form("%s_TotUP%d%s",mode,i,isXi.data())); 
    TotDOWN[i] = (TH1D*)pp2JESTotDOWN[i]->Clone(Form("%s_TotDOWN%d%s",mode,i,isXi.data()));
    JESTotUP[i] = (TH1D*)pp2JESTotUP[i]->Clone(Form("%s_JESTotUP%d%s",mode,i,isXi.data()));
    JESTotDOWN[i] = (TH1D*)pp2JESTotDOWN[i]->Clone(Form("%s_JESTotDOWN%d%s",mode,i,isXi.data()));
    JERTot[i] = (TH1D*)pp2JER[i]->Clone(Form("%s_JERTot%d%s",mode,i,isXi.data())); 
  
    
    for(int j = 0; j<TotUP[0]->GetSize(); j++)
    {
      //if empty bins
      if(pp2JER[i]->GetBinContent(j)==0 || pp2JESTotUP[i]->GetBinContent(j)==-1 || pp2JESTotDOWN[i]->GetBinContent(j)==-1 || pp7JER[i]->GetBinContent(j)==0 || pp7JESTotUP[i]->GetBinContent(j)==-1 || pp7JESTotDOWN[i]->GetBinContent(j)==-1 || pPb5JER[i]->GetBinContent(j)==0 || pPb5JESTotUP[i]->GetBinContent(j)==-1 || pPb5JESTotDOWN[i]->GetBinContent(j)==-1)
      {
        TotUP[i]->SetBinContent(j,0);
        TotDOWN[i]->SetBinContent(j,0);
        TotUP[i]->SetBinError(j,0);
        TotDOWN[i]->SetBinError(j,0);
        JESTotUP[i]->SetBinContent(j,0);
        JESTotDOWN[i]->SetBinContent(j,0);
        JESTotUP[i]->SetBinError(j,0);
        JESTotDOWN[i]->SetBinError(j,0);
        JERTot[i]->SetBinContent(j,0);
        JERTot[i]->SetBinContent(j,0);
        JERTot[i]->SetBinError(j,0);
        JERTot[i]->SetBinError(j,0);
      }
      else 
      {
        JESTotUP[i]->SetBinContent(j,TMath::Power(TMath::Power(pp2JESTotUP[i]->GetBinContent(j),2)+TMath::Power(pp7JESTotUP[i]->GetBinContent(j),2)+TMath::Power(pPb5JESTotUP[i]->GetBinContent(j),2),0.5));
        JESTotDOWN[i]->SetBinContent(j,TMath::Power(TMath::Power(pp2JESTotDOWN[i]->GetBinContent(j),2)+TMath::Power(pp7JESTotDOWN[i]->GetBinContent(j),2)+TMath::Power(pPb5JESTotDOWN[i]->GetBinContent(j),2),0.5));    
        JERTot[i]->SetBinContent(j,TMath::Power(TMath::Power(pp2JER[i]->GetBinContent(j)-1,2)+TMath::Power(pp7JER[i]->GetBinContent(j)-1,2)+TMath::Power(pPb5JER[i]->GetBinContent(j)-1,2),0.5)); 

        TotUP[i]->SetBinContent(j,TMath::Power(JESTotUP[i]->GetBinContent(j),2)+TMath::Power(JERTot[i]->GetBinContent(j),2)+track*track+pdf_uncert*pdf_uncert);     
        TotDOWN[i]->SetBinContent(j,TMath::Power(JESTotDOWN[i]->GetBinContent(j),2)+TMath::Power(JERTot[i]->GetBinContent(j),2)+track*track+pdf_uncert*pdf_uncert);
        TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(pPbPbpDiff[i]->GetBinContent(j),2));
        TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(pPbPbpDiff[i]->GetBinContent(j),2));
        if(MCDiff[i]->GetBinContent(j)<0) TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(MCDiff[i]->GetBinContent(j),2));
        else if(MCDiff[i]->GetBinContent(j)>0) TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(MCDiff[i]->GetBinContent(j),2));
        if(ChargeCutDiff[i]->GetBinContent(j)>0) TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+TMath::Power(ChargeCutDiff[i]->GetBinContent(j),2));
        else if(ChargeCutDiff[i]->GetBinContent(j)<0) TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+TMath::Power(ChargeCutDiff[i]->GetBinContent(j),2));
    
        if(UEtype == 0)
        {     
          TotUP[i]->SetBinContent(j,TotUP[i]->GetBinContent(j)+ TMath::Power(UEDiff[i]->GetBinContent(j)-1,2));
          TotDOWN[i]->SetBinContent(j,TotDOWN[i]->GetBinContent(j)+ TMath::Power(UEDiff[i]->GetBinContent(j)-1,2));
          UEDiff[i]->SetBinContent(j,TMath::Abs(UEDiff[i]->GetBinContent(j)-1));
        }

        TotUP[i]->SetBinContent(j,TMath::Power(TotUP[i]->GetBinContent(j),0.5));
        TotDOWN[i]->SetBinContent(j,TMath::Power(TotDOWN[i]->GetBinContent(j),0.5));
      }
    }
    pp2JESUP[i]->Write();
    pp2JESDOWN[i]->Write();
    pp2JER[i]->Write(); 
    pp7JESUP[i]->Write();
    pp7JESDOWN[i]->Write();
    pp7JER[i]->Write();
    pPb5JESUP[i]->Write();
    pPb5JESDOWN[i]->Write();
    pPb5JER[i]->Write();
    pPbPbpDiff[i]->Write(); 
    MCDiff[i]->Write();
    ChargeCutDiff[i]->Write();
    if(UEtype==0) UEDiff[i]->Write();
    pp2JESTotUP[i]->Write();
    pp2JESTotDOWN[i]->Write(); 
    pp7JESTotUP[i]->Write();
    pp7JESTotDOWN[i]->Write();
    pPb5JESTotUP[i]->Write();
    pPb5JESTotDOWN[i]->Write();  
    JESTotUP[i]->Write();
    JESTotDOWN[i]->Write();
    JERTot[i]->Write();
    TotUP[i]->Write();
    TotDOWN[i]->Write();
  }
  output->Close(); 
}

void systematics(int UEtype=3)
{
  TH1::SetDefaultSumw2();
  FFSystematics("pPb5",UEtype);
  FFSystematics("pp2",UEtype);
  FFSystematics("pp7",UEtype);
  Interpolation_and_Ratio_Systematics("interp",UEtype);
  Interpolation_and_Ratio_Systematics("FFratio",UEtype);
}
