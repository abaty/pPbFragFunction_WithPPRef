//=============================================                                                     
// Author: Chris McGinn                                                                           
//                                                                                                  
// DiJet Analysis Class (MC)                                                                 
//                                                                                            
//=============================================                                                     

#include <TFile.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <iostream>
#include "cfmDiJetIniSkim.h"

//Current # of correction histograms

const Int_t nHistPbPb = 29;
const Int_t nHistPP = 4;

//HI Vs Calo

TFile* VsCaloFile_p[nHistPbPb];
TProfile* VsCalocent_p[nHistPbPb];
TProfile2D* VsCalophiEta_p[nHistPbPb];
TProfile* VsCalopt_p[nHistPbPb];
TProfile* VsCalodelR_p[nHistPbPb];

TFile* FakeVsCaloFile_p[nHistPbPb];
TProfile* FakeVsCalocent_p[nHistPbPb];
TProfile2D* FakeVsCalophiEta_p[nHistPbPb];
TProfile* FakeVsCalopt_p[nHistPbPb];
TProfile* FakeVsCalodelR_p[nHistPbPb];

Int_t posArrPbPb[7] = {0, 5, 10, 15, 20, 25, 28};

//PP Vs Calo

TFile* CaloFile_p[nHistPP];
TProfile2D* CalophiEta_p[nHistPP];
TProfile* Calopt_p[nHistPP];
TProfile* CalodelR_p[nHistPP];

TFile* FakeCaloFile_p[nHistPP];
TProfile2D* FakeCalophiEta_p[nHistPP];
TProfile* FakeCalopt_p[nHistPP];
TProfile* FakeCalodelR_p[nHistPP];

TFile* SecondCaloFile_p;
TH2D* SecondCaloptEta_p;

TFile* MultCaloFile_p;
TH2D* MultCaloptEta_p;

void InitCorrFiles(sampleType sType = kHIDATA, const char * mode = "pp2")
{
  //File names w/ various binnings, ordered by pt and then centrality. Each Jet Algorithm gets a file array

  if(sType == kHIDATA || sType == kHIMC){
    VsCaloFile_p[0] = new TFile("eff_pt50_55_cent0_10.root", "READ");
    VsCaloFile_p[1] = new TFile("eff_pt50_55_cent10_20.root", "READ");
    VsCaloFile_p[2] = new TFile("eff_pt50_55_cent20_30.root", "READ");
    VsCaloFile_p[3] = new TFile("eff_pt50_55_cent30_50.root", "READ");
    VsCaloFile_p[4] = new TFile("eff_pt50_55_cent50_100.root", "READ");
    
    VsCaloFile_p[5] = new TFile("eff_pt55_65_cent0_10.root", "READ");
    VsCaloFile_p[6] = new TFile("eff_pt55_65_cent10_20.root", "READ");
    VsCaloFile_p[7] = new TFile("eff_pt55_65_cent20_30.root", "READ");
    VsCaloFile_p[8] = new TFile("eff_pt55_65_cent30_50.root", "READ");
    VsCaloFile_p[9] = new TFile("eff_pt55_65_cent50_100.root", "READ");
    
    VsCaloFile_p[10] = new TFile("eff_pt65_80_cent0_10.root", "READ");
    VsCaloFile_p[11] = new TFile("eff_pt65_80_cent10_20.root", "READ");
    VsCaloFile_p[12] = new TFile("eff_pt65_80_cent20_30.root", "READ");
    VsCaloFile_p[13] = new TFile("eff_pt65_80_cent30_50.root", "READ");
    VsCaloFile_p[14] = new TFile("eff_pt65_80_cent50_100.root", "READ");
    
    VsCaloFile_p[15] = new TFile("eff_pt80_100_cent0_10.root", "READ");
    VsCaloFile_p[16] = new TFile("eff_pt80_100_cent10_20.root", "READ");
    VsCaloFile_p[17] = new TFile("eff_pt80_100_cent20_30.root", "READ");
    VsCaloFile_p[18] = new TFile("eff_pt80_100_cent30_50.root", "READ");
    VsCaloFile_p[19] = new TFile("eff_pt80_100_cent50_100.root", "READ");
    
    VsCaloFile_p[20] = new TFile("eff_pt100_300_cent0_10.root", "READ");
    VsCaloFile_p[21] = new TFile("eff_pt100_300_cent10_20.root", "READ");
    VsCaloFile_p[22] = new TFile("eff_pt100_300_cent20_30.root", "READ");
    VsCaloFile_p[23] = new TFile("eff_pt100_300_cent30_50.root", "READ");
    VsCaloFile_p[24] = new TFile("eff_pt100_300_cent50_100.root", "READ");

    VsCaloFile_p[25] = new TFile("eff_pt300_800_cent0_10.root", "READ");
    VsCaloFile_p[26] = new TFile("eff_pt300_800_cent10_20.root", "READ");
    VsCaloFile_p[27] = new TFile("eff_pt300_800_cent20_100.root", "READ");
    VsCaloFile_p[28] = new TFile("eff_pt800_30000_cent0_100.root", "READ");
    
    //Fakes VsCalo
    
    FakeVsCaloFile_p[0] = new TFile("fake_pt50_55_cent0_10.root", "READ");
    FakeVsCaloFile_p[1] = new TFile("fake_pt50_55_cent10_20.root", "READ");
    FakeVsCaloFile_p[2] = new TFile("fake_pt50_55_cent20_30.root", "READ");
    FakeVsCaloFile_p[3] = new TFile("fake_pt50_55_cent30_50.root", "READ");
    FakeVsCaloFile_p[4] = new TFile("fake_pt50_55_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[5] = new TFile("fake_pt55_65_cent0_10.root", "READ");
    FakeVsCaloFile_p[6] = new TFile("fake_pt55_65_cent10_20.root", "READ");
    FakeVsCaloFile_p[7] = new TFile("fake_pt55_65_cent20_30.root", "READ");
    FakeVsCaloFile_p[8] = new TFile("fake_pt55_65_cent30_50.root", "READ");
    FakeVsCaloFile_p[9] = new TFile("fake_pt55_65_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[10] = new TFile("fake_pt65_80_cent0_10.root", "READ");
    FakeVsCaloFile_p[11] = new TFile("fake_pt65_80_cent10_20.root", "READ");
    FakeVsCaloFile_p[12] = new TFile("fake_pt65_80_cent20_30.root", "READ");
    FakeVsCaloFile_p[13] = new TFile("fake_pt65_80_cent30_50.root", "READ");
    FakeVsCaloFile_p[14] = new TFile("fake_pt65_80_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[15] = new TFile("fake_pt80_100_cent0_10.root", "READ");
    FakeVsCaloFile_p[16] = new TFile("fake_pt80_100_cent10_20.root", "READ");
    FakeVsCaloFile_p[17] = new TFile("fake_pt80_100_cent20_30.root", "READ");
    FakeVsCaloFile_p[18] = new TFile("fake_pt80_100_cent30_50.root", "READ");
    FakeVsCaloFile_p[19] = new TFile("fake_pt80_100_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[20] = new TFile("fake_pt100_300_cent0_10.root", "READ");
    FakeVsCaloFile_p[21] = new TFile("fake_pt100_300_cent10_20.root", "READ");
    FakeVsCaloFile_p[22] = new TFile("fake_pt100_300_cent20_30.root", "READ");
    FakeVsCaloFile_p[23] = new TFile("fake_pt100_300_cent30_50.root", "READ");
    FakeVsCaloFile_p[24] = new TFile("fake_pt100_300_cent50_100.root", "READ");
    
    FakeVsCaloFile_p[25] = new TFile("fake_pt300_800_cent0_10.root", "READ");
    FakeVsCaloFile_p[26] = new TFile("fake_pt300_800_cent10_20.root", "READ");
    FakeVsCaloFile_p[27] = new TFile("fake_pt300_800_cent20_100.root", "READ");
    FakeVsCaloFile_p[28] = new TFile("fake_pt800_30000_cent0_100.root", "READ");
  }
  else if(sType == kPPDATA || sType == kPPMC){
    /*
    CaloFile_p[0] = new TFile("eff_pt0_1_ak3Calo_dogenjet0.root", "READ");
    CaloFile_p[1] = new TFile("eff_pt1_3_ak3Calo_dogenjet0.root", "READ");
    CaloFile_p[2] = new TFile("eff_pt3_8_ak3Calo_dogenjet0.root", "READ");
    CaloFile_p[3] = new TFile("eff_pt8_300_ak3Calo_dogenjet0.root", "READ");

    FakeCaloFile_p[0] = new TFile("fake_pt0_1_ak3Calo_dogenjet0.root", "READ");
    FakeCaloFile_p[1] = new TFile("fake_pt1_3_ak3Calo_dogenjet0.root", "READ");
    FakeCaloFile_p[2] = new TFile("fake_pt3_8_ak3Calo_dogenjet0.root", "READ");
    FakeCaloFile_p[3] = new TFile("fake_pt8_300_ak3Calo_dogenjet0.root", "READ");*/
    if(strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp5")==0 || strcmp(mode,"pp7")==0)
    {
      CaloFile_p[0] = new TFile("eff_pt0_1_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      CaloFile_p[1] = new TFile("eff_pt1_3_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      CaloFile_p[2] = new TFile("eff_pt3_8_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      CaloFile_p[3] = new TFile("eff_pt8_300_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");

      FakeCaloFile_p[0] = new TFile("fake_pt0_1_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      FakeCaloFile_p[1] = new TFile("fake_pt1_3_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      FakeCaloFile_p[2] = new TFile("fake_pt3_8_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
      FakeCaloFile_p[3] = new TFile("fake_pt8_300_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1_pPb.root", "READ");
    }
    else
    {
      CaloFile_p[0] = new TFile("eff_pt0_1_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      CaloFile_p[1] = new TFile("eff_pt1_3_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      CaloFile_p[2] = new TFile("eff_pt3_8_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      CaloFile_p[3] = new TFile("eff_pt8_300_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");

      FakeCaloFile_p[0] = new TFile("fake_pt0_1_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      FakeCaloFile_p[1] = new TFile("fake_pt1_3_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      FakeCaloFile_p[2] = new TFile("fake_pt3_8_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
      FakeCaloFile_p[3] = new TFile("fake_pt8_300_accept4pt4jtpt4_ak3PF_dogenjet0_doJtPt1.root", "READ");
    }

    SecondCaloFile_p = new TFile("secondary_pp.root", "READ");
    MultCaloFile_p   = new TFile("multiplereco_pp.root", "READ");
  }

  return;
}


void InitCorrHists(sampleType sType = kHIDATA)
{
  if(sType == kHIDATA || sType == kHIMC){
    for(Int_t hIter = 0; hIter < nHistPbPb; hIter++){
      VsCalocent_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_cent");
      VsCalophiEta_p[hIter] = (TProfile2D*)VsCaloFile_p[hIter]->Get("p_eff_acceptance");
      VsCalopt_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_pt");
      VsCalodelR_p[hIter] = (TProfile*)VsCaloFile_p[hIter]->Get("p_eff_rmin");
      
      FakeVsCalocent_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_cent");
      FakeVsCalophiEta_p[hIter] = (TProfile2D*)FakeVsCaloFile_p[hIter]->Get("p_fake_acceptance");
      FakeVsCalopt_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_pt");
      FakeVsCalodelR_p[hIter] = (TProfile*)FakeVsCaloFile_p[hIter]->Get("p_fake_rmin");
    }    
  }
  else if(sType == kPPDATA || sType == kPPMC){
    for(Int_t hIter = 0; hIter < nHistPP; hIter++){
      CalophiEta_p[hIter] = (TProfile2D*)CaloFile_p[hIter]->Get("p_eff_acceptance");
      Calopt_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_pt");
      CalodelR_p[hIter] = (TProfile*)CaloFile_p[hIter]->Get("p_eff_jtpt");

      FakeCalophiEta_p[hIter] = (TProfile2D*)FakeCaloFile_p[hIter]->Get("p_fake_acceptance");
      FakeCalopt_p[hIter] = (TProfile*)FakeCaloFile_p[hIter]->Get("p_fake_pt");
      FakeCalodelR_p[hIter] = (TProfile*)FakeCaloFile_p[hIter]->Get("p_fake_jtpt");
    }

    SecondCaloptEta_p = (TH2D*)SecondCaloFile_p->Get("hpt_eta");
    MultCaloptEta_p = (TH2D*)MultCaloFile_p->Get("hpt_eta"); 
  }

  return;
}


void InitPosArrPbPb(Int_t hiBin)
{
  if(hiBin < 0 || hiBin > 199){
    std::cout << "InitPosArrPbPb: hiBin out of range; check input" << std::endl;
    return;
  }

  Int_t startArrPbPb[7] = {0, 5, 10, 15, 20, 25, 28};
  Int_t hiBinDiv[5] = {20, 40, 60, 100, 200};

  for(Int_t hiBinIter = 0; hiBinIter < 5; hiBinIter++){
    if(hiBin < hiBinDiv[hiBinIter]){
      for(Int_t arrIter = 0; arrIter < 6; arrIter++){
	posArrPbPb[arrIter] = startArrPbPb[arrIter] + hiBinIter;
      }
      break;
    }
  }

  if(hiBin >= 60)
    posArrPbPb[5] = 27;

  return;
}


Int_t getPtBin(Float_t pt, sampleType sType = kHIDATA)
{
  if(pt < .50 || pt > 1000000){
    std::cout << "getPtBin: pt outside of acceptable range; check input" << std::endl;
    return -1;
  }

  if(sType == kHIDATA || sType == kHIMC){    
    Float_t ptArr[8] = {.50, .55, .65, .80, 1.00, 3.00, 8.00, 1000000};
    
    for(Int_t ptPosIter = 0; ptPosIter < 7; ptPosIter++){
      if(pt > ptArr[ptPosIter] && pt <= ptArr[ptPosIter + 1]){
	return posArrPbPb[ptPosIter];
      }
    }
  }
  else if(sType == kPPDATA || sType == kPPMC){
    Float_t ptArr[5] = {.50, 1.00, 3.00, 8.00, 1000000};

    for(Int_t ptPosIter = 0; ptPosIter < 4; ptPosIter++){
      if(pt > ptArr[ptPosIter] && pt <= ptArr[ptPosIter + 1])
	return ptPosIter;
    }
  }

  std::cout << "getPtBin: Shouldn't reach this point; check input" << std::endl;
  return -1;
}


Float_t getEffCorr(Int_t corrBin, Int_t hiBin, Float_t pt, Float_t phi, Float_t eta, Float_t rmin, sampleType sType = kHIDATA)
{
  Float_t effCorr = 1;

  if(sType == kHIDATA || sType == kHIMC){
    effCorr = effCorr*(VsCalocent_p[corrBin]->GetBinContent(VsCalocent_p[corrBin]->FindBin(hiBin)));
    effCorr = effCorr*(VsCalophiEta_p[corrBin]->GetBinContent(VsCalophiEta_p[corrBin]->FindBin(phi, eta)));
    effCorr = effCorr*(VsCalopt_p[corrBin]->GetBinContent(VsCalopt_p[corrBin]->FindBin(pt)));
    effCorr = effCorr*(VsCalodelR_p[corrBin]->GetBinContent(VsCalodelR_p[corrBin]->FindBin(rmin)));
  }
  else if(sType == kPPDATA || kPPMC){
    effCorr = effCorr*(CalophiEta_p[corrBin]->GetBinContent(CalophiEta_p[corrBin]->FindBin(phi, eta)));
    effCorr = effCorr*(Calopt_p[corrBin]->GetBinContent(Calopt_p[corrBin]->FindBin(pt)));
    effCorr = effCorr*(CalodelR_p[corrBin]->GetBinContent(CalodelR_p[corrBin]->FindBin(rmin)));
  }
  return effCorr;
}


Float_t getFakeCorr(Int_t corrBin, Int_t hiBin, Float_t pt, Float_t phi, Float_t eta, Float_t rmin, sampleType sType = kHIDATA)
{
  Float_t fakeCorr = 0;

  if(sType == kHIDATA || sType == kHIMC){
    fakeCorr = fakeCorr + (FakeVsCalocent_p[corrBin]->GetBinContent(FakeVsCalocent_p[corrBin]->FindBin(hiBin)));
    fakeCorr = fakeCorr + (FakeVsCalophiEta_p[corrBin]->GetBinContent(FakeVsCalophiEta_p[corrBin]->FindBin(phi, eta)));
    fakeCorr = fakeCorr + (FakeVsCalopt_p[corrBin]->GetBinContent(FakeVsCalopt_p[corrBin]->FindBin(pt)));
    fakeCorr = fakeCorr + (FakeVsCalodelR_p[corrBin]->GetBinContent(FakeVsCalodelR_p[corrBin]->FindBin(rmin)));
  }
  else if(sType == kPPDATA || sType == kPPMC){
    fakeCorr = fakeCorr + (FakeCalophiEta_p[corrBin]->GetBinContent(FakeCalophiEta_p[corrBin]->FindBin(phi, eta)));
    fakeCorr = fakeCorr + (FakeCalopt_p[corrBin]->GetBinContent(FakeCalopt_p[corrBin]->FindBin(pt)));
    fakeCorr = fakeCorr + (FakeCalodelR_p[corrBin]->GetBinContent(FakeCalodelR_p[corrBin]->FindBin(rmin)));
  }

  return fakeCorr;
}

//Feed variables and the histograms along w/ the appropriate rmincut (currently rmin only defined to 3 for PuPF and 5 for calo)

Float_t factorizedPtCorr(Int_t corrBin, Int_t hiBin, Float_t pt, Float_t phi, Float_t eta, Float_t rmin, sampleType sType = kHIDATA)
{
  if(hiBin < 0 || hiBin > 199){
    std::cout << "factorizedPtCorr: hiBin outside of acceptable range; check input" << std::endl;
    return 1;
  }
  else if(pt < 0.5){
    std::cout << "factorizedPtCorr: pt outside of acceptable range; check input" << std::endl;
  }

  Float_t corrFactor = 1;

  corrFactor = (1 - getFakeCorr(corrBin, hiBin, pt, phi, eta, rmin, sType))/(getEffCorr(corrBin, hiBin, pt, phi, eta, rmin, sType));

  if(sType == kPPDATA || sType == kPPMC){
    corrFactor = corrFactor*(1 - SecondCaloptEta_p->GetBinContent(SecondCaloptEta_p->FindBin(pt, eta)))/(1+MultCaloptEta_p->GetBinContent(MultCaloptEta_p->FindBin(pt,eta)));
  }
  return corrFactor;
}


/*
Example use:

In a directory with correction files open

  InitCorrFiles(sType); // Initializes files
  InitCorrHists(sType); // Initializes histograms

  // Start looping over events

  for(Int_t jentry = 0; jentry < nentries; jentry++){
     tree->GetEntry(jentry);

     InitPosArrPbPb(hiBin); // Initializes array with correct histogram position stored based on hiBin

     // Now loop over tracks and grab corrections
     
     for(Int_t trkIter = 0; trkIter < nTrk; trkIter++){
       // call function to get correct array position (first arg) withing function grabbing correction
       Float_t trkCorr = factorizedPtCorr(getPtBin(trkPt[trkIter], sType), hiBin, trkPt[trkIter], trkPhi[trkIter], trkEta[trkIter], trkRMin[trkIter], sType);
     }
  }

*/
