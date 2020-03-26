# Example for running over a overlay file with the pandora n-tuples
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/numi_uboone_overlay_fhc_mcc9_run1_full_v25.root
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root

./nuexsec --run 1 --mc /uboone/data/users/elenag/PeLEENtupleNuMI/neutrinoselection_filt_run1_NuMI_overlay.root --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_dirt_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_beamon.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_beamoff.root

source merge_run1_files.sh 

./nuexsec --run 1 --hist files/nuexsec_run1_merged.root
