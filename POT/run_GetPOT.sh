# e.g. source run_GetPOT.sh /uboone/data/users/kmistry/work/MCC9/searchingfornues/neutrinoselection_filt.root mc


echo root -l 'GetPOT.C("'$1'","'$2'")'
root -l -q -b 'GetPOT.C("'$1'","'$2'")'
