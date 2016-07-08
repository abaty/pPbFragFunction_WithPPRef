#include <cmath>

const int variations = 34;
const char * variationTag[variations]= {"","_pp2JESUP3","_pp2JESDOWN3","_pp7JESUP3","_pp7JESDOWN3","_pPb5JESUP3","_pPb5JESDOWN3","_pp2JER5","_pp7JER5","_pPb5JER5","_pp2JER2","_pp7JER2","_pPb5JER2","_NoTrackCorr","_pp2JESUP1","_pp2JESDOWN1","_pp7JESUP1","_pp7JESDOWN1","_pPb5JESUP1","_pPb5JESDOWN1","_pp2JESUP2","_pp2JESDOWN2","_pp7JESUP2","_pp7JESDOWN2","_pPb5JESUP2","_pPb5JESDOWN2","_60DegreeCone","_ChargePlus","_ChargeMinus","_XtScaled","_NoChargeCut","_lowPU","_highPU","_midPU"};

TH1D * h_jet;
TH2D * h_track;
TH2D * h_trackUE;
TH2D * h_track_xi;
TH2D * h_trackUE_xi;
TH2D * h_track_rJgT;
TH2D * h_trackUE_rJgT;
TH2D * h_track_xi_rJgT;
TH2D * h_trackUE_xi_rJgT;
TH2D * h_track_gJrT;
TH2D * h_trackUE_gJrT;
TH2D * h_track_xi_gJrT;
TH2D * h_trackUE_xi_gJrT;
TH1D * h_jet_gen;
TH2D * h_track_gen;
TH2D * h_trackUE_gen;
TH2D * h_track_xi_gen;
TH2D * h_trackUE_xi_gen;
TH1D * h_jet_gen_Q;
TH2D * h_track_gen_Q;
TH2D * h_trackUE_gen_Q;
TH2D * h_track_xi_gen_Q;
TH2D * h_trackUE_xi_gen_Q;
TH1D * h_jet_gen_G;
TH2D * h_track_gen_G;
TH2D * h_trackUE_gen_G;
TH2D * h_track_xi_gen_G;
TH2D * h_trackUE_xi_gen_G;
TH1D * h_jet_Q;
TH2D * h_track_Q;
TH2D * h_trackUE_Q;
TH2D * h_track_xi_Q;
TH2D * h_trackUE_xi_Q;
TH1D * h_jet_G;
TH2D * h_track_G;
TH2D * h_trackUE_G;
TH2D * h_track_xi_G;
TH2D * h_trackUE_xi_G;

//simple counting Histograms for tabulating jet numbers
TH1D * totalEvts = new TH1D("totalEvts","",2,-0.5,1.5);
TH1D * totalJetsHist = new TH1D("totalJetsHist","",2,-0.5,1.5);
TH1D * totalJetsRawPtCut = new TH1D("totalJetsRawPtCut","",2,-0.5,1.5);
TH1D * totalJetsChargeSumCut = new TH1D("totalJetsChargeSumCut","",2,-0.5,1.5);
TH1D * totalJetsEtaCutHist = new TH1D("totalJetsEtaCutHist","",2,-0.5,1.5);
TH1D * totalJetsPtCutHist = new TH1D("totalJetsPtCutHist","",2,-0.5,1.5);

TFile * inf, *infMix;
TTree * track, *trackMix;
TTree * ak3PF, *ak3PFMix;
TTree * evt, *evtMix;

const int maxTrack = 25000;
int nTrk, nTrkMix;
int nVtx, nVtxMix;
float trkPt[maxTrack], trkPtMix[maxTrack];
float trkEta[maxTrack], trkEtaMix[maxTrack];
float trkPhi[maxTrack], trkPhiMix[maxTrack];
float trkPtError[maxTrack], trkPtErrorMix[maxTrack];
float trkDzError1[maxTrack], trkDzError1Mix[maxTrack];
float trkDz1[maxTrack], trkDz1Mix[maxTrack];
float trkDxyError1[maxTrack], trkDxyError1Mix[maxTrack];
float trkDxy1[maxTrack], trkDxy1Mix[maxTrack];
float trkFake[maxTrack], trkFakeMix[maxTrack];
int trkCharge[maxTrack], trkChargeMix[maxTrack];
bool highPurity[maxTrack], highPurityMix[maxTrack];
//float trkRmin[maxTrack], trkRminMix[maxTrack];
int nParticle, nParticleMix;
float pEta[2*maxTrack], pEtaMix[2*maxTrack];
float pPhi[2*maxTrack], pPhiMix[2*maxTrack];
float pPt[2*maxTrack], pPtMix[2*maxTrack];
float mtrkPt[2*maxTrack], mtrkPtMix[2*maxTrack];


const int maxJet = 500;
int nref, nrefMix;
float rawpt[maxJet], rawptMix[maxJet];
float jtpt[maxJet], jtptMix[maxJet];
float jteta[maxJet],jtetaMix[maxJet];
float jtphi[maxJet],jtphiMix[maxJet];
float chargedSum[maxJet],chargedSumMix[maxJet];
float refpt[maxJet], refptMix[maxJet];
float refeta[maxJet], refetaMix[maxJet];
float refphi[maxJet], refphiMix[maxJet];
int refparton_flavor[maxJet], refparton_flavorMix[maxJet];
float genChargedSum[maxJet], genChargedSumMix[maxJet];
float pthat, pthatMix;
int ngen, ngenMix;
float genpt[maxJet], genptMix[maxJet];
float geneta[maxJet], genetaMix[maxJet];
float genphi[maxJet], genphiMix[maxJet];

int run, runMix;
int event, eventMix;
float vz, vzMix;
float hiHFplus, hiHFplusMix;
float hiHFminus, hiHFminusMix;
float hiHFplusEta4, hiHFplusEta4Mix;
float hiHFminusEta4, hiHFminusEta4Mix;
float hiHFhitPlus, hiHFhitPlusMix;
float hiHFhitMinus, hiHFhitMinusMix;
float weight, weightMix;

double lowJetPtBound;
double upJetPtBound;

void setJetPtRange(const char * mode, const char* trigger, int doXt = 0)
{
  if((strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp2")==0 || strcmp(mode,"pp5")==0) && strcmp(trigger,"jet80")==0)
  {
    lowJetPtBound = 100; 
    upJetPtBound  = 200;
  }
  if((strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp2")==0 || strcmp(mode,"pp5")==0) && strcmp(trigger,"jet40")==0)
  {
    lowJetPtBound = 60;
    upJetPtBound  = 100;
  }
  if(strcmp(mode,"pp7")==0 && strcmp(trigger,"jet30")==0)
  {
    lowJetPtBound = 60;
    upJetPtBound  = 80;
  }
  if(strcmp(mode,"pp7")==0 && strcmp(trigger,"jet60")==0)
  {
    lowJetPtBound = 80;
    upJetPtBound  = 140;
  }
  if(strcmp(mode,"pp7")==0 && strcmp(trigger,"jet110")==0)
  {
    lowJetPtBound = 140;
    upJetPtBound  = 200;
    if(doXt == 1) upJetPtBound = 300;
  }
  return;
}

void getInputFile(const char * name, int isMC)
{
  inf = new TFile(name,"read");  
  
  track = (TTree*)inf->Get("track"); 
  ak3PF = (TTree*)inf->Get("ak3PF");
  evt = (TTree*)inf->Get("evt");

  //Setting addresses
  track->SetBranchAddress("nTrk",&nTrk); 
  track->SetBranchAddress("trkPt",&trkPt);
  track->SetBranchAddress("trkEta",&trkEta);
  track->SetBranchAddress("trkPhi",&trkPhi);
  track->SetBranchAddress("trkPtError",&trkPtError);
  track->SetBranchAddress("trkDzError1",&trkDzError1);
  track->SetBranchAddress("trkDz1",&trkDz1);
  track->SetBranchAddress("trkDxyError1",&trkDxyError1);
  track->SetBranchAddress("trkDxy1",&trkDxy1);
  track->SetBranchAddress("trkCharge",&trkCharge);
  track->SetBranchAddress("highPurity",&highPurity);
  track->SetBranchAddress("nVtx",&nVtx);

  ak3PF->SetBranchAddress("nref",&nref);
  ak3PF->SetBranchAddress("rawpt",&rawpt);
  ak3PF->SetBranchAddress("jtpt",&jtpt);
  ak3PF->SetBranchAddress("jteta",&jteta);
  ak3PF->SetBranchAddress("jtphi",&jtphi);
  ak3PF->SetBranchAddress("chargedSum",&chargedSum);

  evt->SetBranchAddress("run",&run);
  evt->SetBranchAddress("vz",&vz);
  evt->SetBranchAddress("hiHFplus",&hiHFplus);
  evt->SetBranchAddress("hiHFminus",&hiHFminus);
  evt->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
  evt->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
  evt->SetBranchAddress("hiHFhitPlus",&hiHFhitPlus);
  evt->SetBranchAddress("hiHFhitMinus",&hiHFhitMinus);
  evt->SetBranchAddress("evt",&event);

  if(isMC)
  { 
    track->SetBranchAddress("nParticle",&nParticle);
    track->SetBranchAddress("pEta",&pEta);
    track->SetBranchAddress("pPhi",&pPhi);
    track->SetBranchAddress("pPt",&pPt);
    track->SetBranchAddress("trkFake",&trkFake);
    track->SetBranchAddress("mtrkPt",&mtrkPt);   
 
    ak3PF->SetBranchAddress("refpt",&refpt);
    ak3PF->SetBranchAddress("refeta",&refeta);
    ak3PF->SetBranchAddress("refphi",&refphi);
    ak3PF->SetBranchAddress("refparton_flavor",&refparton_flavor);
    ak3PF->SetBranchAddress("genChargedSum",&genChargedSum);
    ak3PF->SetBranchAddress("pthat",&pthat);
    ak3PF->SetBranchAddress("ngen",&ngen);
    ak3PF->SetBranchAddress("genpt",&genpt);
    ak3PF->SetBranchAddress("geneta",&geneta);
    ak3PF->SetBranchAddress("genphi",&genphi); 

    evt->SetBranchAddress("weight",&weight);
  }
  if(!isMC) weight = 1;

  return;
}

void getInputEntry(int i)
{
  track->GetEntry(i);
  ak3PF->GetEntry(i);
  evt->GetEntry(i);
  return;
}

void getInputFileMix(const char * name, int isMC)
{
  infMix = new TFile(name,"read");  
  
  trackMix = (TTree*)infMix->Get("track"); 
  ak3PFMix = (TTree*)infMix->Get("ak3PF");
  evtMix = (TTree*)infMix->Get("evt");


  //Setting addresses
  trackMix->SetBranchAddress("nTrk",&nTrkMix); 
  trackMix->SetBranchAddress("trkPt",&trkPtMix);
  trackMix->SetBranchAddress("trkEta",&trkEtaMix);
  trackMix->SetBranchAddress("trkPhi",&trkPhiMix);
  trackMix->SetBranchAddress("trkPtError",&trkPtErrorMix);
  trackMix->SetBranchAddress("trkDzError1",&trkDzError1Mix);
  trackMix->SetBranchAddress("trkDz1",&trkDz1Mix);
  trackMix->SetBranchAddress("trkDxyError1",&trkDxyError1Mix);
  trackMix->SetBranchAddress("trkDxy1",&trkDxy1Mix);
  trackMix->SetBranchAddress("trkCharge",&trkChargeMix);
  trackMix->SetBranchAddress("highPurity",&highPurityMix);
  trackMix->SetBranchAddress("nVtx",&nVtxMix);

  ak3PFMix->SetBranchAddress("nref",&nrefMix);
  ak3PFMix->SetBranchAddress("rawpt",&rawptMix);
  ak3PFMix->SetBranchAddress("jtpt",&jtptMix);
  ak3PFMix->SetBranchAddress("jteta",&jtetaMix);
  ak3PFMix->SetBranchAddress("jtphi",&jtphiMix);
  ak3PFMix->SetBranchAddress("chargedSum",&chargedSumMix);

  evtMix->SetBranchAddress("run",&runMix);
  evtMix->SetBranchAddress("vz",&vzMix);
  evtMix->SetBranchAddress("hiHFplus",&hiHFplusMix);
  evtMix->SetBranchAddress("hiHFminus",&hiHFminusMix);
  evtMix->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4Mix);
  evtMix->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4Mix);
  evtMix->SetBranchAddress("hiHFhitPlus",&hiHFhitPlusMix);
  evtMix->SetBranchAddress("hiHFhitMinus",&hiHFhitMinusMix);
  evtMix->SetBranchAddress("evt",&eventMix);

  if(isMC)
  {
    trackMix->SetBranchAddress("nParticle",&nParticleMix);
    trackMix->SetBranchAddress("pEta",&pEtaMix);
    trackMix->SetBranchAddress("pPhi",&pPhiMix);
    trackMix->SetBranchAddress("pPt",&pPtMix);
    trackMix->SetBranchAddress("trkFake",&trkFakeMix);
    trackMix->SetBranchAddress("mtrkPt",&mtrkPtMix);
    
    ak3PFMix->SetBranchAddress("refpt",&refptMix);
    ak3PFMix->SetBranchAddress("refeta",&refetaMix);
    ak3PFMix->SetBranchAddress("refphi",&refphiMix);
    ak3PFMix->SetBranchAddress("refparton_flavor",&refparton_flavorMix);
    ak3PFMix->SetBranchAddress("genChargedSum",&genChargedSumMix);
    ak3PFMix->SetBranchAddress("pthat",&pthatMix);
    ak3PFMix->SetBranchAddress("ngen",&ngenMix);
    ak3PFMix->SetBranchAddress("genpt",&genptMix);
    ak3PFMix->SetBranchAddress("geneta",&genetaMix);
    ak3PFMix->SetBranchAddress("genphi",&genphiMix); 

    evtMix->SetBranchAddress("weight",&weightMix);
  }
  if(!isMC) weightMix = 1;
  return;
}

void getInputEntryMix(int i)
{
  evtMix->GetEntry(i);
  trackMix->GetEntry(i);
  ak3PFMix->GetEntry(i);
  return;
}
