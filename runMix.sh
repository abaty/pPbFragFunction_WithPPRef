if [ $# -ne 1 ]
then
  echo "Usage: ./run.sh <condor_iteration>"
  exit 1
fi

sleep $(($1/2))

echo | awk -v i=$(($1)) -v j=$((2)) '{print "./run.exe masterList.txt masterMBList.txt "i" "j}' 
echo | awk -v i=$(($1)) -v j=$((2)) '{print "./run.exe masterList.txt masterMBList.txt "i" "j}' | bash

echo | awk -v tag=$4 -v user=$USER '{print "mv spectra*.root /mnt/hadoop/cms/store/user/"user"/temporaryStorage/"}'
echo | awk -v tag=$4 -v user=$USER '{print "mv spectra*.root /mnt/hadoop/cms/store/user/"user"/temporaryStorage/"}' | bash
rm *.root
echo "job done successfully"
