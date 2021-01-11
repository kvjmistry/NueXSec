# Example for running over a overlay file with the pandora n-tuples

if [ -z "$1" ]; then

  # # Parallel processing version
  # mc="./nuexsec --run 3 --mc ../ntuples/neutrinoselection_filt_run3b_overlay_weight.root"
  # data="./nuexsec --run 3 --data ../ntuples/neutrinoselection_filt_run3b_beamon_beamgood.root"
  # ext="./nuexsec --run 3 --ext ../ntuples/neutrinoselection_filt_run3b_beamoff.root"
  # dirt="./nuexsec --run 3 --dirt ../ntuples/neutrinoselection_filt_run3b_dirt_overlay.root"

  # eval $mc | tee log/run3_mc.log | sed -e 's/^/[MC] /' &
  # eval $data | tee log/run3_data.log | sed -e 's/^/[Data] /' &
  # eval $ext | tee log/run3_ext.log | sed -e 's/^/[EXT] /' &
  # eval $dirt | tee log/run3_dirt.log | sed -e 's/^/[Dirt] /' &
  
  # # wait for the background processes to finish
  # wait
  
  # # put the outputs into 1 log file
  # declare -a arr=("log/run3_mc.log" "log/run3_data.log" "log/run3_ext.log" "log/run3_dirt.log")
  # echo " " > log/run3.log
  # for i in "${arr[@]}"; do
  #   cat "$i" >> log/run3.log 
  # done

  # # Overwrite the Nue cc events with a higher stats version
  # ./nuexsec --run 3 --mc ../ntuples/neutrinoselection_filt_run3b_overlay_intrinsic.root --intrinsic intrinsic
  
  # Print the selection
  ./nuexsec --run 3 --printonly --printall | tee -a log/run3.log 

  # Merge the files
  source merge/merge_run3_files.sh files/nuexsec_mc_run3.root files/nuexsec_run3_merged.root

  # Run the histogram plotter
  ./nuexsec --run 3 --hist files/nuexsec_run3_merged.root

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("3", true, "files/trees/nuexsec_selected_tree_mc_run3.root", "files/trees/nuexsec_selected_tree_data_run3.root", "files/trees/nuexsec_selected_tree_ext_run3.root","files/trees/nuexsec_selected_tree_dirt_run3.root", "")'

  # Now run the cross section calculator
  # Not running this until happy its working...
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode default | tee -a log/run3.log 

fi

# Running slimmed down version of pelee ntuples with event weights
if [ "$1" == "weight" ]; then

  # new tune
  ./nuexsec --run 3 --var ../ntuples/neutrinoselection_filt_run3b_overlay_weight.root weight
  ./nuexsec --run 3 --var ../ntuples/neutrinoselection_filt_run3b_overlay_intrinsic.root weight --intrinsic intrinsic

  root -l -b -q 'merge/merge_uneaventrees.C("3", true, "files/trees/nuexsec_selected_tree_mc_run3_weight.root", "files/trees/nuexsec_selected_tree_data_run3.root", "files/trees/nuexsec_selected_tree_ext_run3.root","files/trees/nuexsec_selected_tree_dirt_run3.root", "weight")'

  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --xsecmode reweight unisim default
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --xsecmode reweight ppfx default
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --xsecmode reweight genie default
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --xsecmode reweight reint default
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --xsecmode reweight mcstats default

  ./nuexsec --run 3 --sys reweight

  # for running reweighting by cut -- these are slow, so dont run them by default for now
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run3_overlay_weight.root --xsecmode reweight unisim rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run3_overlay_weight.root --xsecmode reweight ppfx rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run3_overlay_weight.root --xsecmode reweight genie rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run3_overlay_weight.root --xsecmode reweight reint rw_cuts --intrinsic intrinsic


fi


