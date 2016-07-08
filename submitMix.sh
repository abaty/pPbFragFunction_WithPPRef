if [ $# -ne 0 ]
then
  echo "Usage: ./psort.sh <trackqual> <file-list> <tag> <nmin> <nmax> <pttrigmin> <pttrigmax> <ptassmin> <ptassmax>"
  exit 1
fi

now="spectra_$(date +"%Y_%m_%d__%H_%M_%S")"
njobs=$((1*($(wc -l < fileLists/masterList.txt))))
#njobs=2

mkdir $now
cp fake*.root $now
cp eff*.root $now
cp secondary*.root $now
cp multiple*.root $now
cp residualJEC.root $now
cp jetJER.root $now
cp fileLists/masterList.txt $now
cp fileLists/masterMBList.txt $now

cp runMix.sh $now/run.sh

#insert blacklist code
blacklist=""
for i in $(cat /net/hisrv0001/home/abaty/condor_blacklist/condor_blacklist.txt); do
  
  blacklist=$blacklist$i" \&\& "
done
blacklist=$blacklist"endoflinetag"
blacklist=$(echo $blacklist | sed "s@ \\\&\\\& endoflinetag@@g")
echo "blacklist: "$blacklist 
cat run.condor | sed "s@blacklist_here@$blacklist@g" > $now/run.condor


cat $now/run.condor | sed "s@log_flag@$now@g" | sed "s@dir_flag@$PWD/$now@g" | sed "s@user_flag@$USER@g" |  sed "s@arglist@@g" | sed "s@transfer_filelist@run.exe@g" | sed "s@njobs@$njobs@g" > $now/run.condor

NAME="Spectra.C"
g++ $NAME $(root-config --cflags --libs) -Werror -Wall -O2 -o "run.exe"
cp run.exe $now
rm run.exe
echo finished compilation
echo
cat $now/run.condor
echo 
echo condor_submit $now/run.condor
echo
# condor_submit $now/run.condor

