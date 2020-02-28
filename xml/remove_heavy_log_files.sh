#!/bin/bash

# Script to remove heavy log files which have been copied to the uboone data area.

# USAGE: run in the directory in the job where the files.list lives
# source remove_heavy_log_files.sh <file name to remove>

logname=$1

default="larStage3.out"

# Check if an input is given, if not then default to $default
if [ -z $logname ]; then 
    
    echo "No file name given, so removing files with name $default"	
    logname=$default
fi

for i in {1..200}
do
    # tail the log files to see the error codes
    tail -n 1 */$logname >> larstage3_logfile.txt
    # remove the log file
    rm -v */$logname
    # sleep for 10 mins
    echo "ON LOOP: $i, sleeping for 10 mins.."
    sleep 3600
done
