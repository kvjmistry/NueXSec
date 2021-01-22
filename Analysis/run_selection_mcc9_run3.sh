# Example for running over a overlay file with the pandora n-tuples

if [ -z "$1" ]; then

  # Parallel processing version
  mc="./nuexsec --run 3 --mc ../ntuples/neutrinoselection_filt_run3b_overlay_newtune.root"
  data="./nuexsec --run 3 --data ../ntuples/neutrinoselection_filt_run3b_beamon_beamgood.root"
  ext="./nuexsec --run 3 --ext ../ntuples/neutrinoselection_filt_run3b_beamoff.root"
  dirt="./nuexsec --run 3 --dirt ../ntuples/neutrinoselection_filt_run3b_dirt_overlay_newtune.root"

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

  # Overwrite the Nue cc events with a higher stats version
  ./nuexsec --run 3 --mc ../ntuples/neutrinoselection_filt_run3b_overlay_intrinsic_newtune.root --intrinsic intrinsic
  
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

  # Electron Energy
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel unisim
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel ppfx
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel genie
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel reint
  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel mcstats

  # ./nuexsec --run 3 --sys reweight

  # Electron Angle
  # ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel unisim  --xsecvar elec_ang
  # ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel ppfx    --xsecvar elec_ang
  # ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel genie   --xsecvar elec_ang
  # ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel reint   --xsecvar elec_ang
  # ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode reweight --xseclabel mcstats --xsecvar elec_ang

  # ./nuexsec --run 3 --sys reweight --xsecvar elec_ang

  # -- 

  # for running reweighting by cut -- these are slow, so dont run them by default for now
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --mc ../ntuples/neutrinoselection_filt_run3_overlay.root --xsecmode reweight --xseclabel unisim --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --mc ../ntuples/neutrinoselection_filt_run3_overlay.root --xsecmode reweight --xseclabel ppfx   --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --mc ../ntuples/neutrinoselection_filt_run3_overlay.root --xsecmode reweight --xseclabel genie  --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --mc ../ntuples/neutrinoselection_filt_run3_overlay.root --xsecmode reweight --xseclabel reint  --xsecplot rw_cuts --intrinsic intrinsic


fi

# Running detector variation samples

# To run the det var samples now do source run_selection_mcc9_run3.sh var <variation name>
# The variation are:
# CV, BNB_Diffusion, LYAttenuation, LYRayleigh, LYDown, SCE, Recomb2, WireModX, WireModYZ, WireModThetaXZ, WireModThetaYZ_withSigmaSplines, WireModThetaYZ_withoutSigmaSplines, WireModdEdX

# Could loop these to make it easier to run them all!

if [ "$1" == "var" ]; then
  ./nuexsec --run 3 --var ../ntuples/detvar_newtune/run3/neutrinoselection_filt_run3_overlay_$2.root $2
  
  # Overwrite the true nue information
  ./nuexsec --run 3 --var ../ntuples/detvar_newtune/run3/intrinsic/neutrinoselection_filt_run3_overlay_$2_intrinsic.root $2 --intrinsic intrinsic

  source merge/merge_run3_files.sh files/nuexsec_mc_run3_$2.root files/nuexsec_run3_$2_merged.root

  ./nuexsec --run 3 --hist files/nuexsec_run3_$2_merged.root --var dummy $2

  root -l -b -q 'merge/merge_uneaventrees.C("3", true, "files/trees/nuexsec_selected_tree_mc_run3_'"$2"'.root", "files/trees/nuexsec_selected_tree_data_run3.root", "files/trees/nuexsec_selected_tree_ext_run3.root","files/trees/nuexsec_selected_tree_dirt_run3.root", "'"$2"'")'

  ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3_$2.root --var dummy $2 --xsecmode default --xsecvar elec_E

fi


#----------------------------------------------------------------------------------------------
# Running the detector systematics

# To run the det sys plotter now do source run_selection_mcc9_run3.sh allvar
# This is just running the piece of code above, but for all the det variations at once, so you do not have to do it manually

# This bit of code will run the piece of code above for all the detector variations, calculate deviation
# the detector variations below are going to be taken into account, please make sure this list matches the _util.vec_var_string in Utility.h

declare -a var=(
  CV
  LYRayleigh
  LYAttenuation
  LYDown
  SCE
  Recomb2
  WireModX
  WireModYZ
  WireModThetaXZ
  WireModThetaYZ_withSigmaSplines
  WireModThetaYZ_withoutSigmaSplines
  WireModdEdX
)
if [ "$1" == "allvar" ]; then

  # run the above script for every det variation 
  for i in "${var[@]}"
  do
    source run_selection_mcc9_run3.sh var "$i"
  done

fi


# To generate the file lists with all event weights then use this command
# ./nuexsec --run 3 --xsec files/trees/nuexsec_tree_merged_run3.root --xsecmode txtlist --xseclabel all --xsecvar elec_E


