# Python script to call zarko's POT counting script and get the number of HW triggers and POT for each run and subrun
# Based off a script written by Wouter
# Input is a run subrun list

from __future__ import print_function
# import numpy as np
from subprocess import check_output
import sys


# create the output files:
file = open("run_subrun_list_data_bad_removed.txt", "r")

# loop over the run subrun filelist
for line in file:
    
    runsubrun = line[:-1].split()
    
    run = int(runsubrun[0])
    subrun = int(runsubrun[1])

    # get the pot and counts:
    call = "/uboone/app/users/zarko/getDataInfo.py -v3 --format-numi --prescale -r {_run} -s {_subrun} ".format(_run = run, _subrun = subrun)
    output = check_output(call, shell=True)
    #print(output)
    lines = output.split('\n')

    EXT_Trig = 0
    # Get the EXT Triggers
    for x in lines:
        if (len(x) == 0):
            continue
        if (x.split()[0] == "EXT_NUMIwin_FEMBeamTriggerAlgo"):
            EXT_Trig = x.split()[1]

    # Get the other POT stuff
    pot = lines[1].split()[9]
    EA9CNT =lines[1].split()[7]
    print(run, subrun, pot, EA9CNT, EXT_Trig )
