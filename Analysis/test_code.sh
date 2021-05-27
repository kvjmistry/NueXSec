#!/bin/bash

echo "----------------------------------------- > rm /files/trees/*"
rm files/trees/*

echo "----------------------------------------- > make clean"
make clean

echo "----------------------------------------- > make"
make

echo "----------------------------------------- > source run_selection_mcc9_run1_gpvm.sh"
./nuexsec --run 1 --mc /pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root
#source run_selection_mcc9_run1_gpvm.sh
