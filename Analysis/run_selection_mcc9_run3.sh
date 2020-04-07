# Example for running over a overlay file with the pandora n-tuples
./nuexsec --run 3 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_beamon_goodbeam.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_beamoff.root --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_dirt_overlay.root --weight 0 2> /dev/null

source merge/merge_run3_files.sh 

./nuexsec --run 3 --hist files/nuexsec_run3_merged.root --weight 0
