#include "TFile.h"
#include "TH1D.h"
#include "TLine.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TColor.h"
#include "TAxis.h"
#include "TPad.h"
#include "TAttMarker.h"
#include "TAttLine.h"
#include "TBox.h"
#include "TF1.h"
#include "TMath.h"

const int FF_Bins = 5;
double FF_Bound[FF_Bins+1] = {60,80,100,120,140,200};

void quad(TH1D* h, int j, double percent){
  h->SetBinContent(j,TMath::Power(h->GetBinContent(j)*h->GetBinContent(j)+percent*percent,0.5));
  return;
}

double quad(double a, double b){
  return TMath::Power(a*a+b*b,0.5);
}

void systematics(int UEtype=2, int alternativeUE = 0)
{
  TH1::SetDefaultSumw2();
  //histos needed for computation
  TH1D* pPbSys[FF_Bins*2];
  TH1D* ppSys[FF_Bins*2];
  TH1D* ratSys[FF_Bins*2];
  TH1D* pPbSys_part[FF_Bins*2][10];
  TH1D* ppSys_part[FF_Bins*2][10];
  TH1D* ratSys_part[FF_Bins*2][10];

  TH1D* pPbFF[FF_Bins*2];
  TH1D* ppFF[FF_Bins*2];
  TH1D* rat[FF_Bins*2];
  
  TH1D* pPbGoing[FF_Bins*2];
  TH1D* PbpGoing[FF_Bins*2];
  
  TH1D* pPbFF_UE[FF_Bins*2];
  TH1D* ppFF_UE[FF_Bins*2];
  TH1D* rat_UE[FF_Bins*2];
  TH1D* pPbFF_UE1[FF_Bins*2];
  TH1D* ppFF_UE1[FF_Bins*2];
  TH1D* rat_UE1[FF_Bins*2];
  
  TH1D* pPbFF_lowPU[FF_Bins*2];
  TH1D* ppFF_lowPU[FF_Bins*2];
  TH1D* rat_lowPU[FF_Bins*2];
  
  TH1D* pPbFF_noChgCut[FF_Bins*2];
  TH1D* ppFF_noChgCut[FF_Bins*2];
  TH1D* rat_noChgCut[FF_Bins*2];

  //JES
  TH1D* ppFF_ppJESU[FF_Bins*2];
  TH1D* rat_ppJESU[FF_Bins*2];
  TH1D* ppFF_ppJESD[FF_Bins*2];
  TH1D* rat_ppJESD[FF_Bins*2];
  TH1D* pPbFF_pPbJESU[FF_Bins*2];
  TH1D* rat_pPbJESU[FF_Bins*2];
  TH1D* pPbFF_pPbJESD[FF_Bins*2];
  TH1D* rat_pPbJESD[FF_Bins*2];
  //smaller values for ratio
  TH1D* ppFF_ppJESU_small[FF_Bins*2];
  TH1D* rat_ppJESU_small[FF_Bins*2];
  TH1D* ppFF_ppJESD_small[FF_Bins*2];
  TH1D* rat_ppJESD_small[FF_Bins*2];
  TH1D* pPbFF_pPbJESU_small[FF_Bins*2];
  TH1D* rat_pPbJESU_small[FF_Bins*2];
  TH1D* pPbFF_pPbJESD_small[FF_Bins*2];
  TH1D* rat_pPbJESD_small[FF_Bins*2];

  //JER
  TH1D* ppFF_ppJER[FF_Bins*2];
  TH1D* rat_ppJER[FF_Bins*2];
  TH1D* pPbFF_pPbJER[FF_Bins*2];
  TH1D* rat_pPbJER[FF_Bins*2];

  //MC nonclsures
  TH1D* ppFF_RnonC[FF_Bins*2];
  TH1D* pPbFF_RnonC[FF_Bins*2];
  TH1D* rat_RnonC[FF_Bins*2];
  TH1D* ppFF_GnonC[FF_Bins*2];
  TH1D* pPbFF_GnonC[FF_Bins*2];
  TH1D* rat_GnonC[FF_Bins*2];

  //files needed for computation
  TFile * inf = TFile::Open(Form("FragmentationFunctionsUE%d.root",UEtype),"read");
  TFile * inf_UE = TFile::Open(Form("FragmentationFunctionsUE%d.root",alternativeUE),"read");
  TFile * inf_UE1 = TFile::Open(Form("FragmentationFunctionsUE%d.root",1),"read");
  TFile * inf_lowPU = TFile::Open(Form("FragmentationFunctions_lowPUUE%d.root",UEtype),"read");
  TFile * inf_noChgCut = TFile::Open(Form("FragmentationFunctions_NoChargeCutUE%d.root",UEtype),"read");
  TFile * inf_ppJESU = TFile::Open(Form("FragmentationFunctions_pprefJESUP2p5UE%d.root",UEtype),"read");
  TFile * inf_ppJESD = TFile::Open(Form("FragmentationFunctions_pprefJESDOWN2p5UE%d.root",UEtype),"read");
  TFile * inf_pPbJESU = TFile::Open(Form("FragmentationFunctions_pPb5JESUP3UE%d.root",UEtype),"read");
  TFile * inf_pPbJESD = TFile::Open(Form("FragmentationFunctions_pPb5JESDOWN3UE%d.root",UEtype),"read");
  TFile * inf_ppJESU_small = TFile::Open(Form("FragmentationFunctions_pprefJESUP1p5UE%d.root",UEtype),"read");
  TFile * inf_ppJESD_small = TFile::Open(Form("FragmentationFunctions_pprefJESDOWN1p5UE%d.root",UEtype),"read");
  TFile * inf_pPbJESU_small = TFile::Open(Form("FragmentationFunctions_pPb5JESUP1p5UE%d.root",UEtype),"read");
  TFile * inf_pPbJESD_small = TFile::Open(Form("FragmentationFunctions_pPb5JESDOWN1p5UE%d.root",UEtype),"read");
  TFile * inf_ppJER = TFile::Open(Form("FragmentationFunctions_ppref5JER5UE%d.root",UEtype),"read");
  TFile * inf_pPbJER = TFile::Open(Form("FragmentationFunctions_pPb5JER5UE%d.root",UEtype),"read");
  TFile * inf_nonC = TFile::Open(Form("FragmentationFunctionsUE%d.root",3),"read");
 
  //loading everything in
  for(int i =0; i<FF_Bins*2; i++)
  {
    std::cout << i << std::endl;   
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";

    //get Result
    pPbFF[i] = (TH1D*)inf->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF[i] = (TH1D*)inf->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat[i] = (TH1D*)inf->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
   
    //beam reversal check
    pPbGoing[i] = (TH1D*)inf->Get(Form("pPb5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    PbpGoing[i] = (TH1D*)inf->Get(Form("Pbp5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
  
    //UE check
    pPbFF_UE[i] = (TH1D*)inf_UE->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_UE[i] = (TH1D*)inf_UE->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_UE[i] = (TH1D*)inf_UE->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbFF_UE1[i] = (TH1D*)inf_UE1->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_UE1[i] = (TH1D*)inf_UE1->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_UE1[i] = (TH1D*)inf_UE1->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
 
    //no PU check 
    pPbFF_lowPU[i] = (TH1D*)inf_lowPU->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_lowPU[i] = (TH1D*)inf_lowPU->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_lowPU[i] = (TH1D*)inf_lowPU->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    
    //no charge Cut
    pPbFF_noChgCut[i] = (TH1D*)inf_noChgCut->Get(Form("pPb5Pbp5TeV_fulldata_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_noChgCut[i] = (TH1D*)inf_noChgCut->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_noChgCut[i] = (TH1D*)inf_noChgCut->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 

    //JES   
    ppFF_ppJESU[i] = (TH1D*)inf_ppJESU->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_ppJESU[i] = (TH1D*)inf_ppJESU->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    ppFF_ppJESD[i] = (TH1D*)inf_ppJESD->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_ppJESD[i] = (TH1D*)inf_ppJESD->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_pPbJESU[i] = (TH1D*)inf_pPbJESU->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_pPbJESU[i] = (TH1D*)inf_pPbJESU->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_pPbJESD[i] = (TH1D*)inf_pPbJESD->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_pPbJESD[i] = (TH1D*)inf_pPbJESD->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    ppFF_ppJESU_small[i] = (TH1D*)inf_ppJESU_small->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_ppJESU_small[i] = (TH1D*)inf_ppJESU_small->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    ppFF_ppJESD_small[i] = (TH1D*)inf_ppJESD_small->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_ppJESD_small[i] = (TH1D*)inf_ppJESD_small->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_pPbJESU_small[i] = (TH1D*)inf_pPbJESU_small->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_pPbJESU_small[i] = (TH1D*)inf_pPbJESU_small->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_pPbJESD_small[i] = (TH1D*)inf_pPbJESD_small->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_pPbJESD_small[i] = (TH1D*)inf_pPbJESD_small->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 

    //JER
    ppFF_ppJER[i] = (TH1D*)inf_ppJER->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_ppJER[i] = (TH1D*)inf_ppJER->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_pPbJER[i] = (TH1D*)inf_pPbJER->Get(Form("ppref5TeV_data_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_pPbJER[i] = (TH1D*)inf_pPbJER->Get(Form("pPbPbp_FF_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    
    //MC nonclosure
    pPbFF_RnonC[i] = (TH1D*)inf_nonC->Get(Form("pPb5Pbp5TeV_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_RnonC[i] = (TH1D*)inf_nonC->Get(Form("pp5TeV_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_RnonC[i] = (TH1D*)inf_nonC->Get(Form("pPbPbp_FF_recoMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    pPbFF_GnonC[i] = (TH1D*)inf_nonC->Get(Form("pPb5Pbp5TeV_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppFF_GnonC[i] = (TH1D*)inf_nonC->Get(Form("pp5TeV_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    rat_GnonC[i] = (TH1D*)inf_nonC->Get(Form("pPbPbp_FF_genMC_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data())); 
    //make sure gen divided by reco is correct here (think so)
    pPbFF_GnonC[i]->Divide(pPbFF_RnonC[i]);
    ppFF_GnonC[i]->Divide(ppFF_RnonC[i]);
    rat_GnonC[i]->Divide(rat_RnonC[i]);

 
    //clone results to get systematics binning
    pPbSys[i] = (TH1D*)pPbFF[i]->Clone(Form("pPbSyst_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ppSys[i] = (TH1D*)ppFF[i]->Clone(Form("ppSyst_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    ratSys[i] = (TH1D*)rat[i]->Clone(Form("ratSyst_%d_%d%s",(int)FF_Bound[i%FF_Bins],(int)FF_Bound[i%FF_Bins+1],isXi.data()));
    pPbSys[i]->Reset();
    ppSys[i]->Reset();
    ratSys[i]->Reset();
    pPbSys[i]->SetDirectory(0);
    ppSys[i]->SetDirectory(0);
    ratSys[i]->SetDirectory(0);
  }

  

  //calculation
  for(int i =0; i<FF_Bins*2; i++)
  { 
    const char* contribution[8] = {"dataMC","PU","pPbPbp","UE","jtchrg","JES","JER","MCnonc"};
    for(int j = 0; j<8; j++){
      pPbSys_part[i][j] = (TH1D*) pPbSys[i]->Clone(Form("pPbSys_part%d_%s",i,contribution[j]));
      ppSys_part[i][j] = (TH1D*) ppSys[i]->Clone(Form("ppSys_part%d_%s",i,contribution[j]));
      ratSys_part[i][j] = (TH1D*) ratSys[i]->Clone(Form("ratSys_part%d_%s",i,contribution[j]));
    }
  
    std::string  isXi = "";
    if(i>=FF_Bins) isXi = "_xi";

    for(int j = 1; j<pPbSys[i]->GetSize()-1; j++){
      //for data MC differences in tracking
      quad(pPbSys[i],j, ((pPbFF[i]->GetBinContent(j)!=0)? 0.039 : 0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0)? 0.039 : 0));
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)? quad(0.039,0.039) : 0)); //(cancellation here?)
         quad(pPbSys_part[i][0],j, ((pPbFF[i]->GetBinContent(j)!=0)? 0.039 : 0));
         quad(ppSys_part[i][0],j, ((ppFF[i]->GetBinContent(j)!=0)? 0.039 : 0));
         quad(ratSys_part[i][0],j,((rat[i]->GetBinContent(j)!=0)? quad(0.039,0.039) : 0)); //(cancellation here?)

      //for pileup check
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs(1-pPbFF_lowPU[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j))):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Min(0.03,TMath::Abs(1-ppFF_lowPU[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j)))  :0));
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Min(0.03,TMath::Abs(1-rat_lowPU[i]->GetBinContent(j)/rat[i]->GetBinContent(j)))    :0));
         quad(pPbSys_part[i][1],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs(1-pPbFF_lowPU[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j))):0));
         quad(ppSys_part[i][1],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Min(0.03,TMath::Abs(1-ppFF_lowPU[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j)))  :0));
         quad(ratSys_part[i][1],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Min(0.03,TMath::Abs(1-rat_lowPU[i]->GetBinContent(j)/rat[i]->GetBinContent(j)))    :0));

      //beam reversal (half the difference of directions), pPb propagates directly to ratio
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs((PbpGoing[i]->GetBinContent(j)-pPbGoing[i]->GetBinContent(j))/(2*pPbFF[i]->GetBinContent(j)))):0));
      quad(ratSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs((PbpGoing[i]->GetBinContent(j)-pPbGoing[i]->GetBinContent(j))/(2*pPbFF[i]->GetBinContent(j)))):0));
         quad(pPbSys_part[i][2],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs((PbpGoing[i]->GetBinContent(j)-pPbGoing[i]->GetBinContent(j))/(2*pPbFF[i]->GetBinContent(j)))):0));
         quad(ratSys_part[i][2],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Min(0.03,TMath::Abs((PbpGoing[i]->GetBinContent(j)-pPbGoing[i]->GetBinContent(j))/(2*pPbFF[i]->GetBinContent(j)))):0));
      
      //Alternative UE check (difference between results)
      bool inThisBin = ((i<5)?(pPbFF[i]->GetBinLowEdge(j)<40):(pPbFF[i]->GetBinLowEdge(j)>1));
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0 && inThisBin)?  TMath::Max(TMath::Abs(1-pPbFF_UE1[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)) ,TMath::Abs(1-pPbFF_UE[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j))):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0 && inThisBin) ?  TMath::Max(TMath::Abs(1-ppFF_UE1[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j)) ,TMath::Abs(1-ppFF_UE[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j))):0));
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0 && inThisBin)  ?  TMath::Max(TMath::Abs(1-rat_UE1[i]->GetBinContent(j)/rat[i]->GetBinContent(j)) ,TMath::Abs(1-rat_UE[i]->GetBinContent(j)/rat[i]->GetBinContent(j))):0));
         quad(pPbSys_part[i][3],j,((pPbFF[i]->GetBinContent(j)!=0 && inThisBin)?  TMath::Max(TMath::Abs(1-pPbFF_UE1[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)) ,TMath::Abs(1-pPbFF_UE[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j))):0));
         quad(ppSys_part[i][3],j, ((ppFF[i]->GetBinContent(j)!=0 && inThisBin) ?  TMath::Max(TMath::Abs(1-ppFF_UE1[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j)) ,TMath::Abs(1-ppFF_UE[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j)))  :0));
         quad(ratSys_part[i][3],j,((rat[i]->GetBinContent(j)!=0 && inThisBin)  ?  TMath::Max(TMath::Abs(1-rat_UE1[i]->GetBinContent(j)/rat[i]->GetBinContent(j)) ,TMath::Abs(1-rat_UE[i]->GetBinContent(j)/rat[i]->GetBinContent(j))):0));

      //no jet charge cut (difference between results)
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_noChgCut[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_noChgCut[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j))  :0));
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_noChgCut[i]->GetBinContent(j)/rat[i]->GetBinContent(j)):0));
         quad(pPbSys_part[i][4],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_noChgCut[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)):0));
         quad(ppSys_part[i][4],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_noChgCut[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j))  :0));
         quad(ratSys_part[i][4],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_noChgCut[i]->GetBinContent(j)/rat[i]->GetBinContent(j)):0));

      //JES for FFs (max of the difference between up and down)
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Max(TMath::Abs(pPbFF_pPbJESU[i]->GetBinContent(j)-pPbFF[i]->GetBinContent(j)),TMath::Abs(pPbFF_pPbJESD[i]->GetBinContent(j)-pPbFF[i]->GetBinContent(j)))/pPbFF[i]->GetBinContent(j):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Max(TMath::Abs(ppFF_ppJESU[i]->GetBinContent(j)-ppFF[i]->GetBinContent(j)),TMath::Abs(ppFF_ppJESD[i]->GetBinContent(j)-ppFF[i]->GetBinContent(j)))/ppFF[i]->GetBinContent(j)  :0));
         quad(pPbSys_part[i][5],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Max(TMath::Abs(pPbFF_pPbJESU[i]->GetBinContent(j)-pPbFF[i]->GetBinContent(j)),TMath::Abs(pPbFF_pPbJESD[i]->GetBinContent(j)-pPbFF[i]->GetBinContent(j)))/pPbFF[i]->GetBinContent(j):0));
         quad(ppSys_part[i][5],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Max(TMath::Abs(ppFF_ppJESU[i]->GetBinContent(j)-ppFF[i]->GetBinContent(j)),TMath::Abs(ppFF_ppJESD[i]->GetBinContent(j)-ppFF[i]->GetBinContent(j)))/ppFF[i]->GetBinContent(j)  :0));
      //same thing for ratios (cancellation here?)
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)?  TMath::Max( TMath::Max(TMath::Abs(rat_ppJESU_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_ppJESD_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j) ,TMath::Max(TMath::Abs(rat_pPbJESU_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_pPbJESD_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j)):0));
      //quad(ratSys[i],j, ((rat[i]->GetBinContent(j)!=0) ?  TMath::Max(TMath::Abs(rat_ppJESU[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_ppJESD[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j)  :0));
         quad(ratSys_part[i][5],j,((rat[i]->GetBinContent(j)!=0)?   TMath::Max( TMath::Max(TMath::Abs(rat_ppJESU_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_ppJESD_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j) ,TMath::Max(TMath::Abs(rat_pPbJESU_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_pPbJESD_small[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j))  :0));
         //quad(ratSys_part[i][5],j, ((rat[i]->GetBinContent(j)!=0) ?  TMath::Max(TMath::Abs(rat_ppJESU[i]->GetBinContent(j)-rat[i]->GetBinContent(j)),TMath::Abs(rat_ppJESD[i]->GetBinContent(j)-rat[i]->GetBinContent(j)))/rat[i]->GetBinContent(j)  :0));
     
      //JER (difference from nominal) 
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_pPbJER[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_ppJER[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j))  :0));
         quad(pPbSys_part[i][6],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_pPbJER[i]->GetBinContent(j)/pPbFF[i]->GetBinContent(j)):0));
         quad(ppSys_part[i][6],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_ppJER[i]->GetBinContent(j)/ppFF[i]->GetBinContent(j))  :0));
      //JER ratio (adding both contributions difference from nominal) 
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)  ? TMath::Max( TMath::Abs(1-rat_ppJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j)) , TMath::Abs(1-rat_pPbJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j))):0));
      //quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_ppJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j)):0));
         quad(ratSys_part[i][6],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Max( TMath::Abs(1-rat_ppJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j)) , TMath::Abs(1-rat_pPbJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j))):0));
         //quad(ratSys_part[i][6],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_ppJER[i]->GetBinContent(j)/rat[i]->GetBinContent(j)):0));

      //nonclosure in MC
      quad(pPbSys[i],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_GnonC[i]->GetBinContent(j)):0));
      quad(ppSys[i],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_GnonC[i]->GetBinContent(j)) :0));
      quad(ratSys[i],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_GnonC[i]->GetBinContent(j)):0));
         quad(pPbSys_part[i][7],j,((pPbFF[i]->GetBinContent(j)!=0)?  TMath::Abs(1-pPbFF_GnonC[i]->GetBinContent(j)):0));
         quad(ppSys_part[i][7],j, ((ppFF[i]->GetBinContent(j)!=0) ?  TMath::Abs(1-ppFF_GnonC[i]->GetBinContent(j)) :0));
         quad(ratSys_part[i][7],j,((rat[i]->GetBinContent(j)!=0)  ?  TMath::Abs(1-rat_GnonC[i]->GetBinContent(j)):0));
    }
    if(i==0){ pPbSys[i]->Print("All"); ppSys[i]->Print("All"); ratSys[i]->Print("All");}
  }

  TFile * output = TFile::Open(Form("SystematicsUE%d.root",UEtype),"recreate");
  for(int i = 0; i<FF_Bins*2; i++){
    pPbSys[i]->SetDirectory(output);
    ppSys[i]->SetDirectory(output);
    ratSys[i]->SetDirectory(output);
    for(int j = 0; j<8; j++){
      if(j>7) continue;
      pPbSys_part[i][j]->SetDirectory(output);
      ratSys_part[i][j]->SetDirectory(output);
      if(j==2) continue;
      ppSys_part[i][j]->SetDirectory(output);
    }
  }

  TCanvas * c1;
  TLegend * leg = new TLegend(0.3,0.4,0.7,0.9);
  c1 = new TCanvas("c1","c1",1200,800);
  float upper[5] = {80,100,120,140,200};
  for(int i = 0; i<FF_Bins*2; i++){
    c1->Clear();
    leg->Clear();
    leg->SetBorderSize(0);
    pPbSys[i]->Draw("c");
    if(i<5){
      pPbSys[i]->GetXaxis()->SetRangeUser(0.5,upper[i]);
      pPbSys[i]->GetXaxis()->SetTitle("p_{T}");
      c1->SetLogx();
    }else{
      pPbSys[i]->GetXaxis()->SetRangeUser(0,5);
      pPbSys[i]->GetXaxis()->SetTitle("#xi");
      c1->SetLogx(0);
    }
    pPbSys[i]->GetYaxis()->SetRangeUser(0,1);
    pPbSys[i]->GetYaxis()->SetTitle("Uncertainty (%)");
    pPbSys_part[i][0]->SetLineColor(kRed); 
    pPbSys_part[i][0]->Draw("c same"); 
    pPbSys_part[i][1]->SetLineColor(kBlue); 
    pPbSys_part[i][1]->Draw("c same"); 
    pPbSys_part[i][2]->SetLineColor(kGreen); 
    pPbSys_part[i][2]->Draw("c same"); 
    pPbSys_part[i][3]->SetLineColor(kOrange); 
    pPbSys_part[i][3]->Draw("c same"); 
    pPbSys_part[i][4]->SetLineColor(kYellow); 
    pPbSys_part[i][4]->Draw("c same"); 
    pPbSys_part[i][5]->SetLineColor(kGray); 
    pPbSys_part[i][5]->Draw("c same"); 
    pPbSys_part[i][6]->SetLineColor(kCyan); 
    pPbSys_part[i][6]->Draw("c same"); 
    pPbSys_part[i][7]->SetLineColor(kMagenta); 
    pPbSys_part[i][7]->Draw("c same"); 
    leg->AddEntry((TObject*)0,"pPb FF", "");
    leg->AddEntry(pPbSys[i],"Total Uncert.", "l");
    leg->AddEntry(pPbSys_part[i][0],"Tracking data vs. MC", "l");
    leg->AddEntry(pPbSys_part[i][1],"Tracking pileup", "l");
    leg->AddEntry(pPbSys_part[i][2],"Beam reversal (pPb/Pbp)", "l");
    leg->AddEntry(pPbSys_part[i][3],"UE subtraction", "l");
    leg->AddEntry(pPbSys_part[i][4],"Jet charge cut", "l");
    leg->AddEntry(pPbSys_part[i][5],"JES", "l");
    leg->AddEntry(pPbSys_part[i][6],"JER", "l");
    leg->AddEntry(pPbSys_part[i][7],"data/MC nonclosure", "l");
    leg->Draw("same");
    c1->SaveAs(Form("sysPlots/pPbsys_%d.png",i));
    c1->SaveAs(Form("sysPlots/pPbsys_%d.pdf",i));
    c1->SaveAs(Form("sysPlots/pPbsys_%d.C",i));

    c1->Clear();
    leg->Clear();
    leg->SetBorderSize(0);
    ppSys[i]->Draw("c");
    if(i<5){
      ppSys[i]->GetXaxis()->SetRangeUser(0.5,upper[i]);
      ppSys[i]->GetXaxis()->SetTitle("p_{T}");
      c1->SetLogx();
    }else{
      ppSys[i]->GetXaxis()->SetRangeUser(0,5);
      ppSys[i]->GetXaxis()->SetTitle("#xi");
      c1->SetLogx(0);
    }
    ppSys[i]->GetYaxis()->SetRangeUser(0,1);
    ppSys[i]->GetYaxis()->SetTitle("Uncertainty (%)");
    ppSys_part[i][0]->SetLineColor(kRed); 
    ppSys_part[i][0]->Draw("c same"); 
    ppSys_part[i][1]->SetLineColor(kBlue); 
    ppSys_part[i][1]->Draw("c same"); 
    //ppSys_part[i][2]->SetLineColor(kGreen); 
    //ppSys_part[i][2]->Draw("c same"); 
    ppSys_part[i][3]->SetLineColor(kOrange); 
    ppSys_part[i][3]->Draw("c same"); 
    ppSys_part[i][4]->SetLineColor(kYellow); 
    ppSys_part[i][4]->Draw("c same"); 
    ppSys_part[i][5]->SetLineColor(kGray); 
    ppSys_part[i][5]->Draw("c same"); 
    ppSys_part[i][6]->SetLineColor(kCyan); 
    ppSys_part[i][6]->Draw("c same"); 
    ppSys_part[i][7]->SetLineColor(kMagenta); 
    ppSys_part[i][7]->Draw("c same"); 
    leg->AddEntry((TObject*)0,"pp FF", "");
    leg->AddEntry(ppSys[i],"Total Uncert.", "l");
    leg->AddEntry(ppSys_part[i][0],"Tracking data vs. MC", "l");
    leg->AddEntry(ppSys_part[i][1],"Tracking pileup", "l");
    leg->AddEntry(ppSys_part[i][3],"UE subtraction", "l");
    leg->AddEntry(ppSys_part[i][4],"Jet charge cut", "l");
    leg->AddEntry(ppSys_part[i][5],"JES", "l");
    leg->AddEntry(ppSys_part[i][6],"JER", "l");
    leg->AddEntry(ppSys_part[i][7],"data/MC nonclosure", "l");
    leg->Draw("same");
    c1->SaveAs(Form("sysPlots/ppsys_%d.png",i));
    c1->SaveAs(Form("sysPlots/ppsys_%d.pdf",i));
    c1->SaveAs(Form("sysPlots/ppsys_%d.C",i));
    
    c1->Clear();
    leg->Clear();
    leg->SetBorderSize(0);
    ratSys[i]->Draw("c");
    if(i<5){
      ratSys[i]->GetXaxis()->SetRangeUser(0.5,upper[i]);
      ratSys[i]->GetXaxis()->SetTitle("p_{T}");
      c1->SetLogx();
    }else{
      ratSys[i]->GetXaxis()->SetRangeUser(0,5);
      ratSys[i]->GetXaxis()->SetTitle("#xi");
      c1->SetLogx(0);
    }
    ratSys[i]->GetYaxis()->SetRangeUser(0,1);
    ratSys[i]->GetYaxis()->SetTitle("Uncertainty (%)");
    ratSys_part[i][0]->SetLineColor(kRed); 
    ratSys_part[i][0]->Draw("c same"); 
    ratSys_part[i][1]->SetLineColor(kBlue); 
    ratSys_part[i][1]->Draw("c same"); 
    ratSys_part[i][2]->SetLineColor(kGreen); 
    ratSys_part[i][2]->Draw("c same"); 
    ratSys_part[i][3]->SetLineColor(kOrange); 
    ratSys_part[i][3]->Draw("c same"); 
    ratSys_part[i][4]->SetLineColor(kYellow); 
    ratSys_part[i][4]->Draw("c same"); 
    ratSys_part[i][5]->SetLineColor(kGray); 
    ratSys_part[i][5]->Draw("c same"); 
    ratSys_part[i][6]->SetLineColor(kCyan); 
    ratSys_part[i][6]->Draw("c same"); 
    ratSys_part[i][7]->SetLineColor(kMagenta); 
    ratSys_part[i][7]->Draw("c same"); 
    leg->AddEntry((TObject*)0,"FF Ratio", "");
    leg->AddEntry(ratSys[i],"Total Uncert.", "l");
    leg->AddEntry(ratSys_part[i][0],"Tracking data vs. MC", "l");
    leg->AddEntry(ratSys_part[i][1],"Tracking pileup", "l");
    leg->AddEntry(ratSys_part[i][2],"Beam reversal (pPb/Pbp)", "l");
    leg->AddEntry(ratSys_part[i][3],"UE subtraction", "l");
    leg->AddEntry(ratSys_part[i][4],"Jet charge cut", "l");
    leg->AddEntry(ratSys_part[i][5],"JES data vs. MC", "l");
    leg->AddEntry(ratSys_part[i][6],"JER", "l");
    leg->AddEntry(ratSys_part[i][7],"MC nonclosure", "l");
    leg->Draw("same");
    c1->SaveAs(Form("sysPlots/ratsys_%d.png",i));
    c1->SaveAs(Form("sysPlots/ratsys_%d.pdf",i));
    c1->SaveAs(Form("sysPlots/ratsys_%d.C",i));
 
  }
  delete c1;

  output->Write();
  output->Close();

  inf->Close();
  inf_UE->Close();
  inf_lowPU->Close();
  inf_noChgCut->Close();
  inf_ppJESU->Close(); 
  inf_ppJESD->Close();
  inf_pPbJESU->Close();
  inf_pPbJESD->Close();
  inf_ppJER->Close();
  inf_pPbJER->Close();
  inf_nonC->Close();
}
