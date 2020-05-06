# Example for running over a overlay file with the pandora n-tuples

if [ -z "$1" ]; then
  ./nuexsec --run 3 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run3b_beamon_beamgood.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run3b_beamoff.root --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_dirt_overlay.root 2> /dev/null | tee log/run3.log

  source merge/merge_run3_files.sh files/nuexsec_mc_run3.root files/nuexsec_run3_merged.root

  ./nuexsec --run 3 --hist files/nuexsec_run3_merged.root

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("3","files/trees/nuexsec_selected_tree_mc_run3.root", "files/trees/nuexsec_selected_tree_data_run3.root", "files/trees/nuexsec_selected_tree_ext_run3.root","files/trees/nuexsec_selected_tree_dirt_run3.root", "")'

  # Now run the cross section calculator
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root 

fi



