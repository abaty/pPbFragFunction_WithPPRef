#include <fstream>
// MC info
const int npPb5MC = 7;
double pp2PthatBounds[8] = {50,80,120,170,220,280,10000,10000};
double pPb5PthatBounds[9]= {50,80,120,170,220,280,370,10000,10000};
double pp7PthatBounds[8] = {30,50,80,120,170,300,470,600};
double pp5PthatBounds[10] = {15,30,50,80,120,170,220,280,1000,1000};
double ppref5PthatBounds[6] = {50,80,120,170,1000,1000};
double crossSection2[7]  = {1.025E-03,9.865E-05,1.129E-05,1.465E-06,2.837E-07,5.323E-08,0};
double crossSection5[8]  = {3.778E-03,4.412E-04,6.147E-05,1.018E-05,2.477E-06,6.160E-07,1.088E-07,0};
double crossSection7[7]  = {5.31E-2, 6.36E-3, 7.84E-4, 1.15E-4, 2.43E-5, 1.17E-6, 7.02E-8};
double crossSectionpp5[9] = {5.335E-01,3.378E-02,3.778E-03,4.412E-04,6.147E-05,1.018E-05,2.477E-06,6.160E-07,0};
//double crossSectionppref5[5] = {4.069E-03,4.959E-04,7.096E-5,1.223E-05,0};//for Pythia8
double crossSectionppref5[5] = {3.778E-03,4.412E-04,6.147E-05,1.018E-05,0};//for Pythia6

TFile * outf;
TTree * track;
TTree * ak3PF;
TTree * evt;

TFile * inf;
TTree * trackIn;
TTree * ak3PFIn;
TTree * akPu3PFIn;
TTree * skimIn;
TTree * evtIn;
TTree * hltIn;

//variables to be written
const int maxTrack = 25000;
int nTrk;
int nVtx;
float trkPt[maxTrack];
float trkEta[maxTrack];
float trkPhi[maxTrack];
float pfEcal[maxTrack];
float pfHcal[maxTrack];
bool  trkFake[maxTrack];
float trkPtError[maxTrack];
float trkDzError1[maxTrack];
float trkDz1[maxTrack];
float trkDxyError1[maxTrack];
float trkDxy1[maxTrack];
int trkCharge[maxTrack];
int highPurity[maxTrack];
//double trkRmin[maxTrack];i
int nParticle;
float pEta[2*maxTrack];
float pPhi[2*maxTrack];
float pPt[2*maxTrack];
float mtrkPt[2*maxTrack];
float mtrkPtError[2*maxTrack];
int mtrkQual[2*maxTrack];
float mtrkChi2[2*maxTrack];
float mtrkNdof[2*maxTrack];
float mtrkDz1[2*maxTrack];
float mtrkDzError1[2*maxTrack];
float mtrkDxy1[2*maxTrack];
float mtrkDxyError1[2*maxTrack];
float mpfEcal[maxTrack];
float mpfHcal[maxTrack];
float zVtx[100];
int maxPtVtx;

const int maxJet = 500;
int nref;
float rawpt[maxJet];
float jtpt[maxJet];
float jteta[maxJet];
float jtphi[maxJet];
float chargedSum[maxJet];
float refpt[maxJet];
float refeta[maxJet];
float refphi[maxJet];
int refparton_flavor[maxJet];
float genChargedSum[maxJet];
float pthat;
int ngen;
float genpt[maxJet];
float geneta[maxJet];
float genphi[maxJet];

int run;
//int lumi;
int event;
float vz;
float hiHFplus;
float hiHFminus;
float hiHFplusEta4;
float hiHFminusEta4;
float hiHFhitPlus;
float hiHFhitMinus;
float weight;

//unwritten variables cut on
int pcollisionEventSelection;
int pPAcollisionEventSelectionPA;
int pPAprimaryVertexFilter;
int pBeamScrapingFilter;
int HLT_PAJet80_NoJetID_v1;
int HLT_PAJet40_NoJetID_v1;
int HLT_PAJet80_NoJetID_v1_Prescl;
int HLT_PAJet40_NoJetID_v1_Prescl;
int HLT_PAZeroBiasPixel_SingleTrack_v1;
int HLT_L1MinimumBiasHF1OR_part5_v1;
int HLT_AK4PFJet40_Eta5p1_v1;
int HLT_AK4PFJet60_Eta5p1_v1;
int HLT_AK4PFJet80_Eta5p1_v1;
int HLT_ZeroBias;
int HLT_Jet30;
int HLT_Jet60;
int HLT_Jet110;
int pHBHENoiseFilter;

int fileSize;

void openOutFile(const char * mode, const char * trigger, int isMC, int date,int iteration)
{
  //outf = new TFile(Form("/mnt/hadoop/cms/store/user/abaty/FF_forests/skims/%s/%s%s_%d_%d_%d.root",mode,mode,trigger,isMC,date,iteration),"recreate");
  outf = new TFile(Form("/export/d00/scratch/abaty/skims/%s%s_%d_%d_%d.root",mode,trigger,isMC,date, iteration),"recreate");
  track = new TTree("track","track");
  ak3PF = new TTree("ak3PF","ak3PF");
  evt   = new TTree("evt","evt");

  track->Branch("nTrk",&nTrk,"nTrk/I");
  track->Branch("trkPt",&trkPt,"trkPt[nTrk]/F");
  track->Branch("trkEta",&trkEta,"trkEta[nTrk]/F");
  track->Branch("trkPhi",&trkPhi,"trkPhi[nTrk]/F");
  track->Branch("pfEcal",&pfEcal,"pfEcal[nTrk]/F");
  track->Branch("pfHcal",&pfHcal,"pfHcal[nTrk]/F");
  track->Branch("trkPtError",&trkPtError,"trkPtError[nTrk]/F");
  track->Branch("trkDzError1",&trkDzError1,"trkDzError1[nTrk]/F");
  track->Branch("trkDz1",&trkDz1,"trkDz1[nTrk]/F");
  track->Branch("trkDxyError1",&trkDxyError1,"trkDxyError1[nTrk]/F");
  track->Branch("trkDxy1",&trkDxy1,"trkDxy1[nTrk]/F");
  track->Branch("trkCharge",&trkCharge,"trkCharge[nTrk]/I");
  track->Branch("highPurity",&highPurity,"highPurity[nTrk]/O");
  track->Branch("nVtx",&nVtx,"nVtx/I");
  //track->Branch("trkRmin",&trkRmin,"trkRmin[nTrk]/F");

  ak3PF->Branch("nref",&nref,"nref/I");
  ak3PF->Branch("rawpt",&rawpt,"rawpt[nref]/F");
  ak3PF->Branch("jtpt",&jtpt,"jtpt[nref]/F");
  ak3PF->Branch("jteta",&jteta,"jteta[nref]/F");
  ak3PF->Branch("jtphi",&jtphi,"jtphi[nref]/F");
  ak3PF->Branch("chargedSum",&chargedSum,"chargedSum[nref]/F");

  evt->Branch("run",&run,"run/I");
  //evt->Branch("lumi",&lumi,"lumi/I");
  evt->Branch("vz",&vz,"vz/F");
  evt->Branch("hiHFplus",&hiHFplus,"hiHFplus/F");
  evt->Branch("hiHFminus",&hiHFminus,"hiHFminus/F"); 
  evt->Branch("hiHFplusEta4",&hiHFplusEta4,"hiHFplusEta4/F");
  evt->Branch("hiHFminusEta4",&hiHFminusEta4,"hiHFminusEta4/F"); 
  evt->Branch("hiHFhitPlus",&hiHFhitPlus,"hiHFhitPlus/F");
  evt->Branch("hiHFhitMinus",&hiHFhitMinus,"hiHFhitMinus/F"); 
  evt->Branch("evt",&event,"evt/I");

  if(isMC)
  {    
    track->Branch("nParticle",&nParticle,"nParticle/I");
    track->Branch("pEta",&pEta,"pEta[nParticle]/F");
    track->Branch("pPhi",&pPhi,"pPhi[nParticle]/F");
    track->Branch("pPt",&pPt,"pPt[nParticle]/F");
    track->Branch("trkFake",&trkFake,"trkFake[nTrk]/O");
    track->Branch("mtrkPt",&mtrkPt,"mtrkPt[nParticle]/F");
    track->Branch("mtrkPtError",&mtrkPtError,"mtrkPtError[nParticle]/F");
    track->Branch("mtrkQual",&mtrkQual,"mtrkQual[nParticle]/I");
    track->Branch("mtrkChi2",&mtrkChi2,"mtrkChi2[nParticle]/F");
    track->Branch("mtrkNdof",&mtrkNdof,"mtrkNdof[nParticle]/F");
    track->Branch("mtrkDz1",&mtrkDz1,"mtrkDz1[nParticle]/F");
    track->Branch("mtrkDzError1",&mtrkDzError1,"mtrkDzError1[nParticle]/F");
    track->Branch("mtrkDxy1",&mtrkDxy1,"mtrkDxy1[nParticle]/F");
    track->Branch("mtrkDxyError1",&mtrkDxyError1,"mtrkDxyError1[nParticle]/F");
    track->Branch("mpfEcal",&mpfEcal,"mpfEcal[nParticle]/F");
    track->Branch("mpfHcal",&mpfHcal,"mpfHcal[nParticle]/F");

    ak3PF->Branch("refpt",&refpt,"refpt[nref]/F");
    ak3PF->Branch("refeta",&refeta,"refeta[nref]/F");
    ak3PF->Branch("refphi",&refphi,"refphi[nref]/F");
    ak3PF->Branch("refparton_flavor",&refparton_flavor,"refparton_flavor[nref]/I");
    ak3PF->Branch("genChargedSum",&genChargedSum,"genChargedSum[nref]/F");
    ak3PF->Branch("pthat",&pthat,"pthat/F");
    ak3PF->Branch("ngen",&ngen,"ngen/I");
    ak3PF->Branch("genpt",&genpt,"genpt[ngen]/F");
    ak3PF->Branch("geneta",&geneta,"geneta[ngen]/F");
    ak3PF->Branch("genphi",&genphi,"genphi[ngen]/F"); 

    evt->Branch("weight",&weight,"weight/F");
  }

  fileSize = 0;
  std::cout << "returning" << std::endl;
  return;
}

int openInFile(const char * name, const char * trigger, const char * mode, int isMC)
{
  //inf = new TFile(name,"read");  
  inf = TFile::Open(name,"read"); 
  if(inf->GetNbytesInfo()==0) return 0;

  trackIn = (TTree*)inf->Get("anaTrack/trackTree");
  if(trackIn == 0) trackIn = (TTree*) inf->Get("ppTrack/trackTree"); 
  ak3PFIn = (TTree*)inf->Get("ak3PFJetAnalyzer/t");
  skimIn = (TTree*)inf->Get("skimanalysis/HltTree");
  evtIn = (TTree*)inf->Get("hiEvtAnalyzer/HiTree");
  hltIn = (TTree*)inf->Get("hltanalysis/HltTree"); 
  akPu3PFIn = (TTree*)inf->Get("akPu3PFJetAnalyzer/t");

 
  //Setting addresses
  trackIn->SetBranchAddress("nTrk",&nTrk); 
  trackIn->SetBranchAddress("trkPt",&trkPt);
  trackIn->SetBranchAddress("trkEta",&trkEta);
  trackIn->SetBranchAddress("trkPhi",&trkPhi);
  trackIn->SetBranchAddress("pfEcal",&pfEcal);
  trackIn->SetBranchAddress("pfHcal",&pfHcal);
  trackIn->SetBranchAddress("trkPtError",&trkPtError);
  trackIn->SetBranchAddress("trkDzError1",&trkDzError1);
  trackIn->SetBranchAddress("trkDz1",&trkDz1);
  trackIn->SetBranchAddress("trkDxyError1",&trkDxyError1);
  trackIn->SetBranchAddress("trkDxy1",&trkDxy1);
  trackIn->SetBranchAddress("trkCharge",&trkCharge);
  trackIn->SetBranchAddress("highPurity",&highPurity);
  trackIn->SetBranchAddress("nVtx",&nVtx);
  trackIn->SetBranchAddress("zVtx",&zVtx);
  trackIn->SetBranchAddress("maxPtVtx",&maxPtVtx); 

  ak3PFIn->SetBranchAddress("nref",&nref);
  ak3PFIn->SetBranchAddress("rawpt",&rawpt);
  ak3PFIn->SetBranchAddress("jtpt",&jtpt);
  ak3PFIn->SetBranchAddress("jteta",&jteta);
  ak3PFIn->SetBranchAddress("jtphi",&jtphi);
  ak3PFIn->SetBranchAddress("chargedSum",&chargedSum);

  if(strcmp(mode,"pp7")!=0 || (!isMC && strcmp(trigger,"MB")!=0))
  {
    evtIn->SetBranchAddress("run",&run);
    //evtIn->SetBranchAddress("lumi",&lumi);
    evtIn->SetBranchAddress("vz",&vz);
    evtIn->SetBranchAddress("hiHFplus",&hiHFplus);
    evtIn->SetBranchAddress("hiHFminus",&hiHFminus);
    evtIn->SetBranchAddress("hiHFplusEta4",&hiHFplusEta4);
    evtIn->SetBranchAddress("hiHFminusEta4",&hiHFminusEta4);
    evtIn->SetBranchAddress("hiHFhitPlus",&hiHFhitPlus);
    evtIn->SetBranchAddress("hiHFhitMinus",&hiHFhitMinus);
    evtIn->SetBranchAddress("evt",&event);
  }

  skimIn->SetBranchAddress("pPAprimaryVertexFilter",&pPAprimaryVertexFilter);
  skimIn->SetBranchAddress("pBeamScrapingFilter",&pBeamScrapingFilter);
  skimIn->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);
  skimIn->SetBranchAddress("pPAcollisionEventSelectionPA",&pPAcollisionEventSelectionPA);
  skimIn->SetBranchAddress("pHBHENoiseFilter",&pHBHENoiseFilter);

  if(strcmp(mode,"pp7")!=0 && strcmp(mode,"ppref5")!=0)
  {
    hltIn->SetBranchAddress("HLT_PAJet80_NoJetID_v1",&HLT_PAJet80_NoJetID_v1);
    hltIn->SetBranchAddress("HLT_PAJet40_NoJetID_v1",&HLT_PAJet40_NoJetID_v1);
    hltIn->SetBranchAddress("HLT_PAJet80_NoJetID_v1_Prescl",&HLT_PAJet80_NoJetID_v1_Prescl);
    hltIn->SetBranchAddress("HLT_PAJet40_NoJetID_v1_Prescl",&HLT_PAJet40_NoJetID_v1_Prescl);
  }
  if(strcmp(mode,"ppref5")==0){
    hltIn->SetBranchAddress("HLT_AK4PFJet40_Eta5p1_v1",&HLT_AK4PFJet40_Eta5p1_v1);
    hltIn->SetBranchAddress("HLT_AK4PFJet60_Eta5p1_v1",&HLT_AK4PFJet60_Eta5p1_v1);
    hltIn->SetBranchAddress("HLT_AK4PFJet80_Eta5p1_v1",&HLT_AK4PFJet80_Eta5p1_v1);
  }
 
  if(strcmp(mode,"pp7")==0)
  {
    for(int vNum = 0; vNum<100; vNum++)
    {
      int status = hltIn->SetBranchAddress(Form("HLT_ZeroBias_v%d",vNum),&HLT_ZeroBias);
      if(status!=-5) break;
    }
    for(int vNum = 0; vNum<100; vNum++)
    {
      int status = hltIn->SetBranchAddress(Form("HLT_Jet30_v%d",vNum),&HLT_Jet30);
      if(status!=-5) break;
    }
    for(int vNum = 0; vNum<100; vNum++)
    {
      int status = hltIn->SetBranchAddress(Form("HLT_Jet60_v%d",vNum),&HLT_Jet60);
      if(status!=-5) break;
    }
    for(int vNum = 0; vNum<100; vNum++)
    {
      int status = hltIn->SetBranchAddress(Form("HLT_Jet110_v%d",vNum),&HLT_Jet110);
      if(status!=-5) break;
    }
  } 

  hltIn->SetBranchAddress("HLT_PAZeroBiasPixel_SingleTrack_v1",&HLT_PAZeroBiasPixel_SingleTrack_v1);
  hltIn->SetBranchAddress("HLT_L1MinimumBiasHF1OR_part5_v1",&HLT_L1MinimumBiasHF1OR_part5_v1);

  if(isMC)
  {
    trackIn->SetBranchAddress("nParticle",&nParticle);
    trackIn->SetBranchAddress("pEta",&pEta);
    trackIn->SetBranchAddress("pPhi",&pPhi);
    trackIn->SetBranchAddress("pPt",&pPt); 
    trackIn->SetBranchAddress("trkFake",&trkFake); 
    trackIn->SetBranchAddress("mtrkPt",&mtrkPt);
    trackIn->SetBranchAddress("mtrkPtError",&mtrkPtError);
    trackIn->SetBranchAddress("mtrkQual",&mtrkQual);
    trackIn->SetBranchAddress("mtrkChi2",&mtrkChi2);
    trackIn->SetBranchAddress("mtrkNdof",&mtrkNdof);
    trackIn->SetBranchAddress("mtrkDz1",&mtrkDz1);
    trackIn->SetBranchAddress("mtrkDzError1",&mtrkDzError1);
    trackIn->SetBranchAddress("mtrkDxy1",&mtrkDxy1);
    trackIn->SetBranchAddress("mtrkDxyError1",&mtrkDxyError1);
    trackIn->SetBranchAddress("mtrkPfEcal",&mpfEcal);
    trackIn->SetBranchAddress("mtrkPfHcal",&mpfHcal);
 
    ak3PFIn->SetBranchAddress("refpt",&refpt);
    ak3PFIn->SetBranchAddress("refeta",&refeta);
    ak3PFIn->SetBranchAddress("refphi",&refphi);
    ak3PFIn->SetBranchAddress("refparton_flavor",&refparton_flavor);
    ak3PFIn->SetBranchAddress("genChargedSum",&genChargedSum);
    ak3PFIn->SetBranchAddress("pthat",&pthat);
    
    if((strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp5")==0) || strcmp(mode,"pp2")==0)
    {
      akPu3PFIn->SetBranchAddress("ngen",&ngen);
      akPu3PFIn->SetBranchAddress("genpt",&genpt);
      akPu3PFIn->SetBranchAddress("geneta",&geneta);
      akPu3PFIn->SetBranchAddress("genphi",&genphi); 
    }
    else if(strcmp(mode,"ppref5")==0)
    {
      ak3PFIn->SetBranchAddress("ngen",&ngen);
      ak3PFIn->SetBranchAddress("genpt",&genpt);
      ak3PFIn->SetBranchAddress("geneta",&geneta);
      ak3PFIn->SetBranchAddress("genphi",&genphi); 
    }
    else
    {
      ak3PFIn->SetBranchAddress("ngen",&ngen);
      ak3PFIn->SetBranchAddress("genpt",&genpt);
      ak3PFIn->SetBranchAddress("geneta",&geneta);
      ak3PFIn->SetBranchAddress("genphi",&genphi);
    }
  }
  /*
  if(!(trackIn->GetEntries() == ak3PFIn->GetEntries() &&
       trackIn->GetEntries() == evtIn->GetEntries() &&
       trackIn->GetEntries() == skimIn->GetEntries() &&
       trackIn->GetEntries() == hltIn->GetEntries())) std::cout << "\n\nWarning, input files have missmatched number of entries!!!\n\n" << std::endl;*/

  return 1;
}
    
int openInFileFast(const char * name, const char * mode, int isMC)
{ 
  inf = TFile::Open(name,"read"); 
  if(inf->GetNbytesInfo()==0) return 0;
  trackIn = (TTree*)inf->Get("anaTrack/trackTree");
  if(trackIn == 0) trackIn = (TTree*) inf->Get("ppTrack/trackTree");
  return 1;
}

void closeOutFile()
{
  outf->cd();
  track->Write("", TObject::kOverwrite); 
  ak3PF->Write("", TObject::kOverwrite);
  evt->Write("", TObject::kOverwrite);  
  outf->Close();
  delete outf;
  return;
}

void closeInFile(int isGoodFile = 1)
{
  inf->cd();
  if(isGoodFile == 1)
  {
    delete trackIn;
    delete evtIn;
    delete ak3PFIn;
    delete akPu3PFIn;
    delete skimIn;
    delete hltIn;
  }
  inf->Close();
  delete inf;
  return;
}

void closeInFileFast(int isGoodFile = 1)
{
  inf->cd();
  if(isGoodFile == 1) delete trackIn;
  inf->Close();
  delete inf;
  return;
}

std::vector<std::string> readInputFileList(const char* fList)
{
  std::string buffer;
  std::vector<std::string> listOfFiles;

  ifstream inFile(fList);
  std::cout << fList << std::endl;
  std::cout << inFile.is_open() << std::endl;

  if(!inFile.is_open())
  {
    std::cout << "Error opening jet file. Exiting." <<std::endl;
    return listOfFiles;
  }
  else
  {
    while(true)
    {
      inFile >> buffer;
      if(inFile.eof()) break;
      listOfFiles.push_back(buffer);
    }
  }
  return listOfFiles;
}
