# Example for running over a overlay file with the pandora n-tuples
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root
#./nuexsec --run 1 --mc /uboone/data/users/elenag/PeLEENtupleNuMI/neutrinoselection_filt_run1_NuMI_overlay.root


if [ -z "$1" ]; then
  # Run the selection
  ./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_overlay.root  --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_dirt_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_beamon_beamgood_morestats.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_beamoff.root  2> /dev/null | tee log/run1.log

  # Merge the files
  source merge/merge_run1_files.sh files/nuexsec_mc_run1.root files/nuexsec_run1_merged.root

  # Run the histogram plotter
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged.root

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "")'

  # Now run the cross section calculator
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root
fi
# ---------------------


# Running detector variation samples

# CV
if [ $1 == "CV" ]; then
  ./nuexsec --run 1 --var /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_detvar/run1/neutrinoselection_filt_run1_overlay_CV.root CV 2> /dev/null | tee log/run1_CV.log

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_CV.root files/nuexsec_run1_CV_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_CV_merged.root --var dummy CV

  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1_CV.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "CV")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_CV.root --var dummy CV

fi

# ----------------------


# BNB_Diffusion
if [ $1 == "BNB_Diffusion" ]; then
  ./nuexsec --run 1 --var /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_detvar/run1/neutrinoselection_filt_run1_overlay_diffusion.root BNB_Diffusion 2> /dev/null | tee log/run1_BNB_Diffusion.log

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_BNB_Diffusion.root files/nuexsec_run1_BNB_Diffusion_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_BNB_Diffusion_merged.root --var dummy BNB_Diffusion

  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1_BNB_Diffusion.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "BNB_Diffusion")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_BNB_Diffusion.root --var dummy BNB_Diffusion

fi




