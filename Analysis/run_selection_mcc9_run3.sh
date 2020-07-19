# Example for running over a overlay file with the pandora n-tuples

if [ -z "$1" ]; then
  #./nuexsec --run 3 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_overlay.root --data /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run3b_beamon_beamgood.root --ext /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_v2/neutrinoselection_filt_run3b_beamoff.root --dirt /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run3b_dirt_overlay.root 2> /dev/null | tee log/run3.log

  # Parallel processing version
  mc="./nuexsec --run 3 --mc ../ntuples/neutrinoselection_filt_run3b_overlay.root"
  data="./nuexsec --run 3 --data ../ntuples/neutrinoselection_filt_run3b_beamon_beamgood.root"
  ext="./nuexsec --run 3 --ext ../ntuples/neutrinoselection_filt_run3b_beamoff.root"
  dirt="./nuexsec --run 3 --dirt ../ntuples/neutrinoselection_filt_run3b_dirt_overlay.root"

  eval $mc | tee log/run3_mc.log | sed -e 's/^/[MC] /' &
  eval $data | tee log/run3_data.log | sed -e 's/^/[Data] /' &
  eval $ext | tee log/run3_ext.log | sed -e 's/^/[EXT] /' &
  eval $dirt | tee log/run3_dirt.log | sed -e 's/^/[Dirt] /' &
  
  # wait for the background processes to finish
  wait
  
  # put the outputs into 1 log file
  declare -a arr=("log/run3_mc.log" "log/run3_data.log" "log/run3_ext.log" "log/run3_dirt.log")
  echo " " > log/run3.log
  for i in "${arr[@]}"; do
    cat "$i" >> log/run3.log 
  done
  
  ./nuexsec --run 3 --printonly --printall | tee -a log/run3.log 


  source merge/merge_run3_files.sh files/nuexsec_mc_run3.root files/nuexsec_run3_merged.root

  ./nuexsec --run 3 --hist files/nuexsec_run3_merged.root &
  wait

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("3","files/trees/nuexsec_selected_tree_mc_run3.root", "files/trees/nuexsec_selected_tree_data_run3.root", "files/trees/nuexsec_selected_tree_ext_run3.root","files/trees/nuexsec_selected_tree_dirt_run3.root", "")'

  # Now run the cross section calculator
  # Not running this until happy its working...
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode default | tee -a log/run3.log 

fi



