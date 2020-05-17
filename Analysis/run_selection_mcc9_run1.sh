# Example for running over a overlay file with the pandora n-tuples
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root
#./nuexsec --run 1 --mc /uboone/data/users/elenag/PeLEENtupleNuMI/neutrinoselection_filt_run1_NuMI_overlay.root


if [ -z "$1" ]; then
  # Run the selection
  # ./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run1_overlay.root  --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_dirt_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_beamon_beamgood.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_beamoff.root --printall  2> /dev/null | tee log/run1.log

  # use the BNB ntuples as a cross check for the analysis
  # mc="./nuexsec --run 1 --mc /uboone/data/users/davidc/searchingfornues/v08_00_00_33/cc0pinp/0304/run1/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root --weight 2"
  # data="./nuexsec --run 1 --data /uboone/data/users/davidc/searchingfornues/v08_00_00_33/cc0pinp/0304/run1/data_bnb_mcc9.1_v08_00_00_25_reco2_C1_beam_good_reco2_5e19.root"
  # ext="./nuexsec --run 1 --ext /uboone/data/users/davidc/searchingfornues/v08_00_00_33/cc0pinp/0304/run1/data_extbnb_mcc9.1_v08_00_00_25_reco2_C1_all_reco2.root"
  # dirt="./nuexsec --run 1 --dirt /uboone/data/users/davidc/searchingfornues/v08_00_00_33/cc0pinp/0304/run1/prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root --weight 2"

  # Parallel processing version
  mc="./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_overlay.root --weight 3"
  data="./nuexsec --run 1 --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_beamon_beamgood.root"
  ext="./nuexsec --run 1 --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_beamoff.root"
  dirt="./nuexsec --run 1 --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v3/neutrinoselection_filt_run1_dirt_overlay.root --weight 3"

  eval $mc | tee log/run1_mc.log | sed -e 's/^/[MC] /' &
  eval $data | tee log/run1_data.log | sed -e 's/^/[Data] /' &
  eval $ext | tee log/run1_ext.log | sed -e 's/^/[EXT] /' &
  eval $dirt | tee log/run1_dirt.log | sed -e 's/^/[Dirt] /' &
  
  # wait for the background processes to finish
  wait

  # put the outputs into 1 log file
  declare -a arr=("log/run1_mc.log" "log/run1_data.log" "log/run1_ext.log" "log/run1_dirt.log")
  echo " " > log/run1.log
  for i in "${arr[@]}"; do
    cat "$i" >> log/run1.log 
  done

  # Print the selection
  ./nuexsec --run 1 --printonly --printall | tee -a log/run1.log 
  
  # Merge the files
  source merge/merge_run1_files.sh files/nuexsec_mc_run1.root files/nuexsec_run1_merged.root

  # Run the histogram plotter
  plotter="./nuexsec --run 1 --hist files/nuexsec_run1_merged.root --weight 3"
  eval $plotter &
  wait

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "")'

  # Now run the cross section calculator
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root | tee -a log/run1.log 
fi
# ---------------------


# Running detector variation samples

# CV
if [ "$1" == "CV" ]; then
  ./nuexsec --run 1 --var /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_detvar/run1/neutrinoselection_filt_run1_overlay_CV.root CV 2> /dev/null | tee log/run1_CV.log

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_CV.root files/nuexsec_run1_CV_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_CV_merged.root --var dummy CV

  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1_CV.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "CV")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_CV.root --var dummy CV

fi

# ----------------------


# BNB_Diffusion
if [ "$1" == "BNB_Diffusion" ]; then
  ./nuexsec --run 1 --var /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_detvar/run1/neutrinoselection_filt_run1_overlay_diffusion.root BNB_Diffusion 2> /dev/null | tee log/run1_BNB_Diffusion.log

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_BNB_Diffusion.root files/nuexsec_run1_BNB_Diffusion_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_BNB_Diffusion_merged.root --var dummy BNB_Diffusion

  root -l -b -q 'merge/merge_uneaventrees.C("1","files/trees/nuexsec_selected_tree_mc_run1_BNB_Diffusion.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "BNB_Diffusion")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_BNB_Diffusion.root --var dummy BNB_Diffusion

fi




