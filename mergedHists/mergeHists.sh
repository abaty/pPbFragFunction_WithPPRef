#!/bin/bash

now="processed_$(date +"%Y_%m_%d__%H_%M_%S")"

mkdir $now
hadd $now/pPb5_UE0_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_0_UE0_0_15.root
hadd $now/Pbp5_UE0_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_0_UE0_0_15.root
hadd $now/pPb5MC_UE0_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_1_UE0_0_15.root
hadd $now/Pbp5MC_UE0_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_1_UE0_0_15.root
hadd $now/pp5MC_UE0_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapp5*_1_UE0_0_15.root
hadd $now/pPb5_UE3_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_0_UE3_0_15.root
hadd $now/Pbp5_UE3_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_0_UE3_0_15.root
hadd $now/pPb5MC_UE3_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_1_UE3_0_15.root
hadd $now/Pbp5MC_UE3_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_1_UE3_0_15.root
hadd $now/pp5MC_UE3_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapp5*_1_UE3_0_15.root
hadd $now/pPb5_UE2_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_0_UE2_0_15.root
hadd $now/Pbp5_UE2_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_0_UE2_0_15.root
hadd $now/pPb5MC_UE2_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapPb5*_1_UE2_0_15.root
hadd $now/Pbp5MC_UE2_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectraPbp5*_1_UE2_0_15.root
hadd $now/pp5MC_UE2_0_15.root /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectrapp5*_1_UE2_0_15.root
echo 'removing files in 30 seconds'
sleep 15
echo 'removing files in 15 seconds'
sleep 10
echo 'removing files in 5 seconds'
sleep 5
echo 'removing files'
rm /mnt/hadoop/cms/store/user/abaty/temporaryStorage/spectra*.root
