# Script to get the path of all files in the selected data run subrun list

# Usage: source filelist_from_runsubrunlist.sh <corresponding data sam definition to run subrun list> <the run subrun list>

samdef=$1
runlist=$2

counter=0

for line in $(cat $runlist); do
   
    if [ $counter -eq 0 ]; then
	    run=$line
    elif [ $counter -eq 1  ]; then
	    subrun=$line
	    runsubrun="$run.$subrun"
            #echo "$runsubrun"
	    samweb list-files "defname:$samdef and run_number $runsubrun" | tee -a ../bin/files_$(basename $runlist)
    fi

    counter=$(($counter+1))
    if [ $counter -eq 3 ]; then counter=0; fi

done
