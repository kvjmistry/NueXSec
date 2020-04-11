# Example for running over a overlay file with the pandora n-tuples
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root
#./nuexsec --run 1 --mc /uboone/data/users/elenag/PeLEENtupleNuMI/neutrinoselection_filt_run1_NuMI_overlay.root

# Run the selection
./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_overlay.root  --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_dirt_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_beamon_sample0_beam_good_bad.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_beamoff.root 2> /dev/null | tee log/run1.log

# Merge the files
source merge/merge_run1_files.sh 

# Run the histogram plotter
./nuexsec --run 1 --hist files/nuexsec_run1_merged.root