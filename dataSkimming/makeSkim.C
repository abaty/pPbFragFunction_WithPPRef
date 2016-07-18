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
  const int maxOutputFileSize = 700000;
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
  }
  else 
  {
    //if((strcmp("pPb5",mode)==0) && strcmp(trigger,"jet80")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pPb_5_02TeV_pA2013/PA2013_HiForest_PromptReco_JSonPPb_forestv77.root"); 
    //if((strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet80")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/Pbp_5_02TeV_pA2013/PA2013_HiForest_PromptReco_JSonPbp_JECdb_forestv84_Jet80_100.root");
    if((strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet80")==0) fileList = readInputFileList("pPb5Pbp5_jet80_data_files.txt");
    if((strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0) && strcmp(trigger,"jet40")==0) fileList = readInputFileList("pPb5Pbp5_jet40_data_files.txt");
    if((strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0) && strcmp(trigger,"MB")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pPb_5_02TeV_pA2013/PA2013_HiForest_PromptReco_KrisztianMB_JSonPPb_forestv84.root");
    if(strcmp("pp2",mode)==0 && strcmp(trigger,"jet80")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pp_2_76TeV_pp2013/PP2013_HiForest_PromptReco_JsonPP_Jet80_PPReco_forestv82.root");
    if(strcmp("pp2",mode)==0 && strcmp(trigger,"jet40")==0) fileList.push_back("/mnt/hadoop/cms/store/user/abaty/FF_forests/data/pp_2_76TeV_pp2013/PP2013_HiForest_PromptReco_JSon_Jet40Jet60_ppTrack_forestv84.root");
    if(strcmp("pp2",mode)==0 && strcmp(trigger,"MB")==0) fileList.push_back("/mnt/hadoop/cms/store/user/luck/pp_minbiasSkim_forest_53x_2013-08-15-0155/pp_minbiasSkim_forest_53x_2013-08-15-0155.root");
    if(strcmp("pp7",mode)==0 && strcmp(trigger,"MB")==0) fileList.push_back("/mnt/hadoop/cms/store/user/velicanu/Run2011A/MinimumBias2/FOREST/12Oct2013-v1-merged/pp_7TeV_2011A_12Oct2013-v1.root");
    else if(strcmp("pp7",mode)==0) fileList = readInputFileList("pp7Files.txt");

    if(strcmp("ppref5",mode)==0 && (strcmp(trigger,"MB")==0)) fileList = readInputFileList("pp_MinBias6.txt");
    if(strcmp("ppref5",mode)==0 && (strcmp(trigger,"jet40")==0 || strcmp(trigger,"jet60")==0)) fileList = readInputFileList("HighPtLowerJets.txt");
    if(strcmp("ppref5",mode)==0 && strcmp(trigger,"jet80")==0) fileList = readInputFileList("HIPhysicsJet80.txt");
  }
  int nFiles = fileList.size();

  //pp7 MC reweighting
  int totalEntriesForWeighting[7] = {0};
  if(strcmp(mode, "pp7")==0 && isMC)
  {
    for(int f=0; f<nFiles; f++)
    { 
      if(f%100 == 0) std::cout << "tabulating entries; file " << f << "/" << nFiles <<std::endl;
      int isGoodFile = openInFileFast(fileList[f].data(),mode,isMC);
      if( isGoodFile == 0)
      { 
        closeInFileFast(0);
        continue;
      }

      for(int i = 0; i<7; i++)
      {
        if(fileList[f].find(Form("%dto%d",(int)pp7PthatBounds[i],(int)pp7PthatBounds[i+1])) != std::string::npos) totalEntriesForWeighting[i] += trackIn->GetEntries();
      }
      closeInFileFast();
    }
    for(int i = 0 ; i<7; i++) std::cout << pp7PthatBounds[i] << " has " << totalEntriesForWeighting[i] << " entries." <<std::endl;
  }

//start of skim here 
  int outFileNum = 0;
  //looping over forests to skim out of
  //change f= at 2 spots to change starting point, as well as skim outFileNum
  for(int f = 0; f<nFiles; f++)
  //for(int f = 0; f<4; f++)
  {   
    int isGoodFile = openInFile(fileList[f].data(),trigger,mode,isMC);
     std::cout << 1 << std::endl;
    if(f==0) openOutFile(mode,trigger,isMC,date,outFileNum);
    std::cout << 0 << std::endl;
    if( isGoodFile == 0)
    {
      closeInFile(0);
      continue;
    }
    std::cout << 0 << std::endl;

    int nEntries = trackIn->GetEntries();
    //nEntries = 50;
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
      if(strcmp(mode,"ppref5")==0)
      {
        if(strcmp(trigger,"jet80") == 0 && HLT_AK4PFJet80_Eta5p1_v1 == 0) continue;
        if(strcmp(trigger,"jet60") == 0 && HLT_AK4PFJet60_Eta5p1_v1 == 0) continue;
        if(strcmp(trigger,"jet40") == 0 && HLT_AK4PFJet40_Eta5p1_v1 == 0) continue;
        if(strcmp(trigger,"MB") == 0 && HLT_L1MinimumBiasHF1OR_part5_v1 == 0) continue;
        afterHLTCut++;
      }

      //Filling Trees 
      trackIn->GetEntry(i);   
      ak3PFIn->GetEntry(i);
      if(isMC && (strcmp("pPb5",mode)==0 || strcmp("Pbp5",mode)==0 || strcmp("pp5",mode)==0 || strcmp("pp2",mode)==0)) akPu3PFIn->GetEntry(i);

      if(isMC)
      {
        if(strcmp("pp2",mode)==0)
        {
          if(pthat > pp2PthatBounds[f+2]) continue;
          weight = crossSection2[f]/(float)nEntries;
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
        if(strcmp("pp7",mode)==0)
        {
          for(int k = 0; k<7; k++)
          {
            if(fileList[f].find(Form("%dto%d",(int)pp7PthatBounds[k],(int)pp7PthatBounds[k+1])) != std::string::npos) weight = crossSection7[k]/(float)totalEntriesForWeighting[k];
          }    
        }
      }
      if((strcmp("pp7",mode)==0 && strcmp("MB",trigger)==0 && !isMC) || strcmp("ppref5",mode)==0) vz=zVtx[maxPtVtx];
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

