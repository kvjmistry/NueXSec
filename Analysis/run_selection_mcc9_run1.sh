# Example for running over a overlay file with the pandora n-tuples
#./nuexsec --run 1 --mc /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files/neutrinoselection_filt_run1_overlay.root
#./nuexsec --run 1 --mc /uboone/data/users/elenag/PeLEENtupleNuMI/neutrinoselection_filt_run1_NuMI_overlay.root

if [ ! -d "log" ]; then 
  echo
  echo "Log folder doesnt exist... creating"
  mkdir -p log
fi


if [ -z "$1" ]; then
  # Run the selection
  mc="./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root"
  data="./nuexsec --run 1 --data ../ntuples/neutrinoselection_filt_run1_beamon_beamgood.root"
  ext="./nuexsec --run 1 --ext ../ntuples/neutrinoselection_filt_run1_beamoff.root"
  dirt="./nuexsec --run 1 --dirt ../ntuples/neutrinoselection_filt_run1_dirt_overlay_newtune.root"

  # This runs each of the strings above in parallel to maximise cpu usage
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

  # Overwrite the Nue cc events with a higher stats version
  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root --intrinsic intrinsic

  # Print the selection
  ./nuexsec --run 1 --printonly --printall | tee -a log/run1.log 
  
  # Merge the files
  source merge/merge_run1_files.sh files/nuexsec_mc_run1.root files/nuexsec_run1_merged.root

  # Run the histogram plotter
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged.root
  # ./nuexsec --run 1 --hist files/nuexsec_run1_merged.root --plotsys tot

  # Merge the ttrees to one file
  root -l -b -q 'merge/merge_uneaventrees.C("1", true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "")'

  # Now run the cross section calculator
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode default | tee -a log/run1.log 
fi
# ---------------------

# Running slimmed down version of pelee ntuples with event weights
if [ "$1" == "weight" ]; then

  # old tune
  # ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_weight.root weight
  # ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic.root weight --intrinsic intrinsic

  # new tune
  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root weight
  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root weight --intrinsic intrinsic

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, "files/trees/nuexsec_selected_tree_mc_run1_weight.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "weight")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --xsecmode reweight unisim default
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --xsecmode reweight ppfx default
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --xsecmode reweight genie default
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --xsecmode reweight reint default
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --xsecmode reweight mcstats default

  ./nuexsec --run 1 --sys reweight

  # for running reweighting by cut -- these are slow, so dont run them by default for now
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run1_overlay_weight.root --xsecmode reweight unisim rw_cuts
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run1_overlay_weight.root --xsecmode reweight ppfx rw_cuts
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run1_overlay_weight.root --xsecmode reweight genie rw_cuts
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_weight.root --var dummy weight --mc ../ntuples/neutrinoselection_filt_run1_overlay_weight.root --xsecmode reweight reint rw_cuts


fi


# Running detector variation samples

# To run the det var samples now do source run_selection_mcc9_run1.sh var <variation name>
# The variation are:
# CV, BNB_Diffusion, LYAttenuation, LYRayleigh, LYDown, SCE, Recomb2, WireModX, WireModYZ, WireModThetaXZ, WireModThetaYZ_withSigmaSplines, WireModThetaYZ_withoutSigmaSplines, WireModdEdX

# Could loop these to make it easier to run them all!

if [ "$1" == "var" ]; then
  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_$2.root $2
  
  # Overwrite the true nue information
  # ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_$2_intrinsic.root $2 --intrinsic intrinsic

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_$2.root files/nuexsec_run1_$2_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_$2_merged.root --var dummy $2

  root -l -b -q 'merge/merge_uneaventrees.C("1", false, "files/trees/nuexsec_selected_tree_mc_run1_'"$2"'.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "'"$2"'")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_$2.root --var dummy $2 --xsecmode default

fi


#----------------------------------------------------------------------------------------------
# Running the detector systematics

# To run the det sys plotter now do source run_selection_mcc9_run1.sh allvar
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
    source run_selection_mcc9_run1.sh var "$i"
  done

fi



