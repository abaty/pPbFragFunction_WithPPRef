//=============================================                                 
// Author: Chris McGinn                                                         
//                                                                              
// DiJet Analysis Skim Class (MC)                                                    
//                                                                              
//=============================================  
#ifndef cfmDiJetAnaSkim_h
#define cfmDiJetAnaSkim_h

#include "/net/hisrv0001/home/cfmcginn/DijetMPt/CMSSW_5_3_12_patch3/src/DijetInitialSkim/cfmDiJetIniSkim.h"

TTree* trackTreeAna_p;
TTree* jetTreeAna_p;
TTree* genTreeAna_p;

//Track Tree Variables

//Tracks proj. onto Alg (enum ordered above, w/ corrected in back 5), All, Cone, and NotCone

//2nd Dim 0 == 0_1, 1 == 1_2 ... 5 == F


Float_t rAlgImbProjA_[6][6];
Float_t rAlgImbProjAC_[6][6];
Float_t rAlgImbProjANC_[6][6];

//DelR ProjA

Float_t rAlgImbProjAR_[6][6][10];

//Jet Tree Variables

Int_t run_;
Int_t evt_;
Int_t lumi_;
Int_t hiBin_;

Float_t pthat_;

Float_t hiEvtPlane_;
Float_t psin_;

//Set Bool

//Event Set Bool array, [0] == PuPF, [1] == PuCalo, .etc according to enum

Bool_t eventSet_[3];
Float_t centWeight_[3];
Float_t centWeight_2pi3_[3];
Float_t centWeight_120_5pi6_[3];
Float_t centWeight_hatAll_0_[3];
Float_t centWeight_hatAll_1_[3];

Float_t fullWeight_[3];

Bool_t isQuarkJet_[3];
Bool_t isGluonJet_[3];

Float_t pthatWeight_;

//Jet Set, Array by algorithm, according to enum above, 2nd, 0 = lead, 1 =sublead etc.

Float_t AlgJtPt_[3][4];
Float_t AlgJtPhi_[3][4];
Float_t AlgJtEta_[3][4];
Float_t AlgJtTrkMax_[3][4];

Float_t AlgJtMult_[3][2];

Float_t AlgJtAvePhi_[3];
Float_t AlgJtDelPhi_[3];
Float_t AlgJtAsymm_[3];

Float_t AlgRefPt_[3][4];
Float_t AlgRefPhi_[3][4];
Float_t AlgRefEta_[3][4];

Float_t AlgRefAvePhi_[3];
Float_t AlgRefDelPhi_[3];
Float_t AlgRefAsymm_[3];

//Gen Tree Variables

//Gen. proj. onto Jets, ordered by algorithm according to enum, PuPF == [0], PuCalo == [1], etc.

Float_t gAlgImbProjA_[3][6];
Float_t gAlgImbProjAC_[3][6];
Float_t gAlgImbProjANC_[3][6];

//truth delRs
//ProjA delR

Float_t gAlgImbProjAR_[3][6][10];

void SetAnaBranches(Bool_t montecarlo = false, sampleType sType = kHIDATA)
{
  //Track Tree Branches

  std::cout << "Branches Set" << std::endl;

  //Tracks proj. onto Alg, ordered according to enum above, All, Cone, and NotCone

  trackTreeAna_p->Branch("rAlgImbProjA", &rAlgImbProjA_, "rAlgImbProjA[6][6]/F");
  trackTreeAna_p->Branch("rAlgImbProjAC", &rAlgImbProjAC_, "rAlgImbProjAC[6][6]/F");
  trackTreeAna_p->Branch("rAlgImbProjANC", &rAlgImbProjANC_, "rAlgImbProjANC[6][6]/F");

  //ProjA DelRs

  trackTreeAna_p->Branch("rAlgImbProjAR", &rAlgImbProjAR_, "rAlgImbProjAR[6][6][10]/F");

  //Jet Tree Branches

  jetTreeAna_p->Branch("run", &run_, "run/I");
  jetTreeAna_p->Branch("evt", &evt_, "evt/I");
  jetTreeAna_p->Branch("lumi", &lumi_, "lumi/I");

  if(sType == kHIDATA || sType == kHIMC){
    jetTreeAna_p->Branch("hiBin", &hiBin_, "hiBin/I");
    jetTreeAna_p->Branch("hiEvtPlane", &hiEvtPlane_, "hiEvtPlane/F");
    jetTreeAna_p->Branch("psin", &psin_, "psin/F");
  }

  jetTreeAna_p->Branch("eventSet", &eventSet_, "eventSet[3]/O");

  if(sType == kHIMC){
    jetTreeAna_p->Branch("centWeight", &centWeight_, "centWeight[3]/F");
    jetTreeAna_p->Branch("centWeight_2pi3", &centWeight_2pi3_, "centWeight_2pi3[3]/F");
    jetTreeAna_p->Branch("centWeight_120_5pi6", &centWeight_120_5pi6_, "centWeight_120_5pi6[3]/F");
    jetTreeAna_p->Branch("centWeight_hatAll_0", &centWeight_hatAll_0_, "centWeight_hatAll_0[3]/F");
    jetTreeAna_p->Branch("centWeight_hatAll_1", &centWeight_hatAll_1_, "centWeight_hatAll_1[3]/F");
  }    

  jetTreeAna_p->Branch("fullWeight", &fullWeight_, "fullWeight[3]/F");

  jetTreeAna_p->Branch("AlgJtPt", &AlgJtPt_, "AlgJtPt[3][4]/F");
  jetTreeAna_p->Branch("AlgJtPhi", &AlgJtPhi_, "AlgJtPhi[3][4]/F");
  jetTreeAna_p->Branch("AlgJtEta", &AlgJtEta_, "AlgJtEta[3][4]/F");
  jetTreeAna_p->Branch("AlgJtTrkMax", &AlgJtTrkMax_, "AlgJtTrkMax[3][4]/F");

  jetTreeAna_p->Branch("AlgJtMult", &AlgJtMult_, "AlgJtMult[3][2]/F");

  jetTreeAna_p->Branch("AlgJtAvePhi", &AlgJtAvePhi_, "AlgJtAvePhi[3]/F");
  jetTreeAna_p->Branch("AlgJtDelPhi", &AlgJtDelPhi_, "AlgJtDelPhi[3]/F");
  jetTreeAna_p->Branch("AlgJtAsymm", &AlgJtAsymm_, "AlgJtAsymm[3]/F");

  if(montecarlo){
    jetTreeAna_p->Branch("isQuarkJet", &isQuarkJet_, "isQuarkJet[3]/O");
    jetTreeAna_p->Branch("isGluonJet", &isGluonJet_, "isGluonJet[3]/O");

    jetTreeAna_p->Branch("pthatWeight", &pthatWeight_, "pthatWeight/F");
    jetTreeAna_p->Branch("pthat", &pthat_, "pthat/F");

    //refpt for jets immediately above
    jetTreeAna_p->Branch("AlgRefPt", &AlgRefPt_, "AlgRefPt[3][4]/F");
    jetTreeAna_p->Branch("AlgRefPhi", &AlgRefPhi_, "AlgRefPhi[3][4]/F");
    jetTreeAna_p->Branch("AlgRefEta", &AlgRefEta_, "AlgRefEta[3][4]/F");

    jetTreeAna_p->Branch("AlgRefAvePhi", &AlgRefAvePhi_, "AlgRefAvePhi[3]/F");
    jetTreeAna_p->Branch("AlgRefDelPhi", &AlgRefDelPhi_, "AlgRefDelPhi[3]/F");
    jetTreeAna_p->Branch("AlgRefAsymm", &AlgRefAsymm_, "AlgRefAsymm[3]/F");


    //Gen. proj. onto jetAlg, array ordered according to enum

    genTreeAna_p->Branch("gAlgImbProjA", &gAlgImbProjA_, "gAlgImbProjA[3][6]/F");
    genTreeAna_p->Branch("gAlgImbProjAC", &gAlgImbProjAC_, "gAlgImbProjAC[3][6]/F");
    genTreeAna_p->Branch("gAlgImbProjANC", &gAlgImbProjANC_, "gAlgImbProjANC[3][6]/F");

    //Truth Del Rs, Proj
    //Proj A's

    genTreeAna_p->Branch("gAlgImbProjAR", &gAlgImbProjAR_, "gAlgImbProjAR[3][6][10]/F");
  }
}


void InitDiJetAnaSkim(Bool_t montecarlo = false, sampleType sType = kHIDATA)
{
  std::cout << "Init DiJet AnaSkim" << std::endl;

  trackTreeAna_p = new TTree("trackTreeAna", "trackTreeAna");
  jetTreeAna_p = new TTree("jetTreeAna", "jetTreeAna");

  if(montecarlo)
    genTreeAna_p = new TTree("genTreeAna", "genTreeAna");

  SetAnaBranches(montecarlo, sType);
}


void InitJetVar(Bool_t montecarlo = false, sampleType sType = kHIDATA)
{
  for(Int_t initIter = 0; initIter < 3; initIter++){
    eventSet_[initIter] = false;

    if(sType == kHIMC){
      centWeight_[initIter] = 1;
      centWeight_2pi3_[initIter] = 1;
      centWeight_120_5pi6_[initIter] = 1;
      centWeight_hatAll_0_[initIter] = 1;
      centWeight_hatAll_1_[initIter] = 1;
    }

    fullWeight_[initIter] = 1;

    for(Int_t initIter2 = 0; initIter2 < 4; initIter2++){
      AlgJtPt_[initIter][initIter2] = -10;
      AlgJtPhi_[initIter][initIter2] = -10;
      AlgJtEta_[initIter][initIter2] = -10;
      AlgJtTrkMax_[initIter][initIter2] = -10;

      if(initIter2 < 2)
	AlgJtMult_[initIter][initIter2] = 0;
    }

    AlgJtAvePhi_[initIter] = -10;
    AlgJtDelPhi_[initIter] = -10;
    AlgJtAsymm_[initIter] = -10;

    if(montecarlo){
      isQuarkJet_[initIter] = false;
      isGluonJet_[initIter] = false;

      pthatWeight_ = -10;

      for(Int_t initIter2 = 0; initIter2 < 4; initIter2++){
	AlgRefPt_[initIter][initIter2] = -10;
	AlgRefPhi_[initIter][initIter2] = -10;
	AlgRefEta_[initIter][initIter2] = -10;
      }

      AlgRefAvePhi_[initIter] = -10;
      AlgRefDelPhi_[initIter] = -10;
      AlgRefAsymm_[initIter] = -10;
    }
  }
}

void InitProjPerp(Bool_t montecarlo = false)
{
  //Tracks proj. onto Alg, ordered according to enum above, corr in the back 5, All, Cone, and NotCone

  for(Int_t initIter = 0; initIter < 6; initIter++){
    for(Int_t initIter2 = 0; initIter2 < 6; initIter2++){
      rAlgImbProjA_[initIter][initIter2] = 0;
      rAlgImbProjAC_[initIter][initIter2] = 0;
      rAlgImbProjANC_[initIter][initIter2] = 0;

      for(Int_t initIter3 = 0; initIter3 < 10; initIter3++){
	rAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
      }
    }
  }

  if(montecarlo){
    //Gen. proj. onto Truth
    for(Int_t initIter = 0; initIter < 3; initIter++){
      for(Int_t initIter2 = 0; initIter2 < 6; initIter2++){
	gAlgImbProjA_[initIter][initIter2] = 0;
	gAlgImbProjAC_[initIter][initIter2] = 0;
	gAlgImbProjANC_[initIter][initIter2] = 0;

	for(Int_t initIter3 = 0; initIter3 < 10; initIter3++){
	  gAlgImbProjAR_[initIter][initIter2][initIter3] = 0;
	}
      }
    }    
  }
}


#endif
