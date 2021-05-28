#!/bin/bash

echo "----- removing files"
rm files/trees/*

echo "----- compiling"
make clean

make

echo "----- running code"
./nuexsec --run 1 --mc /pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root --gpvm
