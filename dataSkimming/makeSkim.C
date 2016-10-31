//written by Austin Baty, 2/27/2015
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TDatime.h"
#include <iostream>
#include "skimUtilities.h"

void makeSkim(const char * mode = "pp2", const char * trigger = "jet80",int isMC = 0)
{
  if(!(strcmp(mode,"pp2")==0 || strcmp(mode,"pp7")==0 || strcmp(mode,"pPb5")==0 || strcmp(mode,"Pbp5")==0 || strcmp(mode,"pp5")==0 || strcmp(mode,"ppref5")==0))
  {
    std::cout << "invalid mode, terminating execution" << std::endl;
    return;
  }
  if(!(strcmp(trigger,"jet80") == 0 || strcmp(trigger,"jet40") == 0 || strcmp(trigger,"jet60") == 0 || strcmp(trigger,"MB") == 0) && !(strcmp(mode,"pp7") == 0))
  {
    std::cout << "invalid trigger, terminating execution" << std::endl;
    return;
  }
  if(!(strcmp(trigger,"jet30") == 0 || strcmp(trigger,"jet60") == 0 || strcmp(trigger,"110") == 0 || strcmp(trigger,"MB") == 0) && strcmp(mode,"pp7") == 0)
  {
    std::cout << "invalid trigger, terminating execution" << std::endl;
    return;
  }


  //event counters to check how many events we use/cut out
  int totalEvents =0, afterRunCut = 0, afterNoiseCut = 0, afterPCollCut = 0, afterHLTCut = 0;
  
  //date execution started for labeling output files
  TDatime * time = new TDatime();
  int date = time->GetDate();

  //some parameters for  files
  //max output file size
  const int maxOutputFileSize = 60000;
  //if(isMC) 

//setting up files 
  std::vector<std::string> fileList;
  if(isMC)
  {
    if(strcmp("pp2",mode)==0) fileList = readInputFileList("pp2MCFiles.txt");
    if(strcmp("pPb5",mode)==0) fileList = readInputFileList("pPbMCFiles.txt");
    if(strcmp("Pbp5",mode)==0) fileList = readInputFileList("PbpMCFiles.txt");
    if(strcmp("pp7",mode)==0) fileList = readInputFileList("pp7MCFiles.txt");
    if(strcmp("pp5",mode)==0) fileList = readInputFileList("pp5MCFiles.txt");
    if(strcmp("ppref5",mode)==0) fileList = readInputFileList("ppref5_Pythia6.txt");
  }
  else 
  {
    //if((strcmp("pPb5",mode)==0) && strcmp(trigger,"jet80")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pPb_5_02TeV_pA2013/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root"); 
    //if((strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet80")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/Pbp_5_02TeV_pA2013/PA2013_HiForest_PromptReco_JSonPbp_JECdb_forestv84_Jet80_100.root");
    if((strcmp("pPb5",mode)==0) && strcmp(trigger,"jet80")==0) fileList = readInputFileList("pPb5_jet_data_files_v2.txt");
    if((strcmp("pPb5",mode)==0) && strcmp(trigger,"jet40")==0) fileList = readInputFileList("pPb5_jet_data_files_v2.txt");
    if((strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet80")==0) fileList = readInputFileList("Pbp5_jet_data_files_v2.txt");
    if((strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet40")==0) fileList = readInputFileList("Pbp5_jet_data_files_v2.txt");
    if((strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0) && strcmp(trigger,"MB")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pPb_5_02TeV_pA2013/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root");

    if(strcmp("ppref5",mode)==0 && (strcmp(trigger,"MB")==0)) fileList = readInputFileList("pp_MinBias6.txt");
    if(strcmp("ppref5",mode)==0 && (strcmp(trigger,"jet40")==0 || strcmp(trigger,"jet60")==0)) fileList = readInputFileList("HighPtLowerJets.txt");
    if(strcmp("ppref5",mode)==0 && strcmp(trigger,"jet80")==0) fileList = readInputFileList("HIPhysicsJet80.txt");
  }
  int nFiles = fileList.size();

  //pp7 MC reweighting
  int totalEntriesForWeighting[7] = {0};

//start of skim here 
  int outFileNum = 0;
  //looping over forests to skim out of
  //change f= at 2 spots to change starting point, as well as skim outFileNum
  for(int f = 0; f<nFiles; f++)
  //for(int f = 0; f<4; f++)
  { 
    std::cout << "opening file" << std::endl;  
    int isGoodFile = openInFile(fileList[f].data(),trigger,mode,isMC);
    std::cout << 1 << std::endl;
    std::cout << fileList[f].data() << std::endl;
    if(f==0) openOutFile(mode,trigger,isMC,date,outFileNum);
    std::cout << 0 << std::endl;
    if( isGoodFile == 0)
    {
      closeInFile(0);
      continue;
    }
    std::cout << 0 << std::endl;

    int nEntries = trackIn->GetEntries();
    //nEntries = 500;
    for(int i = 0; i<nEntries; i++)
    {
      totalEvents++;
      if(i%10000==0) std::cout <<"file: " << f << " event: " << i << "/" << nEntries << std::endl;
      if(strcmp(mode,"pp7")!=0 || (!isMC && strcmp(trigger,"MB")!=0)) evtIn->GetEntry(i);
      skimIn->GetEntry(i);

      //event and run selections (veto the other-going way as well as first 7 runs for misalignment)
      if(isMC == 0)
      {
        if(strcmp(mode,"pPb5") == 0 && (run>=211313)) continue;
        if(strcmp(mode,"Pbp5") == 0 && run<211313) continue;
        afterRunCut++;
        if(strcmp(mode,"ppref5")!=0 && pHBHENoiseFilter == 0) continue;
        afterNoiseCut++;
      }
      if(strcmp(mode,"ppref5")!=0 && pPAcollisionEventSelectionPA == 0 && pcollisionEventSelection == 0) continue; 
      if(strcmp(mode,"ppref5")==0 && (pPAprimaryVertexFilter==0 || pBeamScrapingFilter == 0)) continue;
      afterPCollCut++;
 
      //trigger selection
      hltIn->GetEntry(i);
      if(strcmp(mode,"pp5")!=0 && strcmp(mode,"ppref5")!=0 && !(strcmp(mode,"pp7")==0 && strcmp(trigger,"MB")==0))
      {
        if(strcmp(trigger,"jet80") == 0 && HLT_PAJet80_NoJetID_v1 == 0) continue;
        if(strcmp(trigger,"jet40") == 0 && HLT_PAJet40_NoJetID_v1 == 0) continue;
        if(strcmp(trigger,"MB") == 0 && HLT_PAZeroBiasPixel_SingleTrack_v1 == 0) continue;
        if(strcmp(trigger,"jet30") == 0 && HLT_Jet30 == 0 ) continue;
        if(strcmp(trigger,"jet60") == 0 && HLT_Jet60 == 0 ) continue;
        if(strcmp(trigger,"jet110") == 0 && HLT_Jet110 == 0 ) continue;
        afterHLTCut++;
      }

      //Filling Trees
      ak3PFIn->GetEntry(i);
      if(strcmp(mode,"ppref5")==0)
      {
        if(!isMC){
          if(strcmp(trigger,"jet80") == 0 && HLT_AK4PFJet80_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"jet60") == 0 && HLT_AK4PFJet60_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"jet40") == 0 && HLT_AK4PFJet40_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"MB") == 0 && HLT_L1MinimumBiasHF1OR_part5_v1 == 0) continue;
        }else{
          if(strcmp(trigger,"jet80") == 0 && HLT_AK4PFJet80_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"jet60") == 0 && HLT_AK4PFJet60_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"jet40") == 0 && HLT_AK4PFJet40_Eta5p1_v1 == 0) continue;
          if(strcmp(trigger,"MB") == 0 && HLT_L1MinimumBiasHF1OR_part5_v1 == 0) continue;
        }
        afterHLTCut++;
      }
      
      trackIn->GetEntry(i);   
      if(isMC && (strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0 || strcmp("pp5",mode)==0 || strcmp("pp2",mode)==0)) akPu3PFIn->GetEntry(i);

      if(isMC)
      {
        if(strcmp("ppref5",mode)==0)
        {
          if(pthat > ppref5PthatBounds[f+2]) continue;
          weight = crossSectionppref5[f]/(float)nEntries;
        }
        if(strcmp("pp5",mode)==0 || strcmp("pPbv25",mode)==0)
        {
          if(pthat > pp5PthatBounds[f+2]) continue;
          weight = crossSectionpp5[f]/(float)nEntries;
        }
        if(strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0)
        {
          if(pthat > pPb5PthatBounds[f+2]) continue; 
          weight = crossSection5[f]/(float)nEntries;
        }
      }
      if((strcmp("pp7",mode)==0 && strcmp("MB",trigger)==0 && !isMC) || strcmp("ppref5",mode)==0) vz=zVtx[maxPtVtx];

      if(isMC && strcmp("ppref5",mode)==0){
        nParticle=0;
        genIn->GetEntry(i);
        for(int t = 0; t<mult; t++){
          if(TMath::Abs(eta->at(t))>2.4) continue;
          if(chg->at(t)==0) continue;
          if(pt->at(t)<0.5) continue;
          pPt[nParticle] = (float)pt->at(t);
          pEta[nParticle] = (float)eta->at(t);
          pPhi[nParticle] = (float)phi->at(t);
          //chg[t] = (float)chg.at(nParticle); 
          nParticle++;
        }
      }

      track->Fill();
      ak3PF->Fill();
      evt->Fill();
      fileSize++; 
      //open a new output file if the old one is "full"
      if(fileSize>=maxOutputFileSize)
      {
        closeOutFile();
        std::cout << "Filled file " << outFileNum << ". Opening new file." << std::endl;
        outFileNum++;              
        openOutFile(mode,trigger,isMC,date,outFileNum);   
      }
    }
    //cleanup so we can open another
    closeInFile();  
  }
  closeOutFile();

  //output number of files used
  std::cout << "Total Events scanned: " << totalEvents << std::endl;
  std::cout << "Total Events after Run cuts " << afterRunCut << std::endl;
  std::cout << "Total Events after Noise cuts: " << afterNoiseCut << std::endl;
  std::cout << "Total Events after PColl cut: " << afterPCollCut << std::endl;
  std::cout << "Total Events after HLT cuts: " << afterHLTCut << std::endl;
}

