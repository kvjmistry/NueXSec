# Script with commands to run the entire selection on the gpvm

if [ ! -d "log" ]; then 
  echo
  echo "Log folder doesnt exist... creating"
  mkdir -p log
fi


if [ -z "$1" ]; then
  # Run the selection
  mc="./nuexsec --run 1 --mc /pnfs/uboone/persistent/users/davidc/searchingfornues/v08_00_00_43/0928/prodgenie_bnb_nu_uboone_overlay_mcc9.1_v08_00_00_26_filter_run1_reco2_reco2.root --gpvm"
  data="./nuexsec --run 1 --data /uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/farsidebands/run1_neutrinoselection_filt_1e_2showers_sideband_skimmed_extended_v47.root --gpvm"
  ext="./nuexsec --run 1 --ext /uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/run1/nslice/data_extbnb_mcc9.1_v08_00_00_25_reco2_C_all_reco2.root --gpvm"
  dirt="./nuexsec --run 1 --dirt /uboone/data/users/davidc/searchingfornues/v08_00_00_43/0702/run1/nslice/prodgenie_bnb_dirt_overlay_mcc9.1_v08_00_00_26_run1_reco2_reco2.root --gpvm"

  # This runs each of the strings above in parallel to maximise cpu usage
  eval $mc | tee log/run1_mc.log | sed -e 's/^/[MC] /' &
  wait
  eval $data | tee log/run1_data.log | sed -e 's/^/[Data] /' &
  wait
  eval $ext | tee log/run1_ext.log | sed -e 's/^/[EXT] /' &
  wait
  eval $dirt | tee log/run1_dirt.log | sed -e 's/^/[Dirt] /' &
  wait

  # put the outputs into 1 log file
  declare -a arr=("log/run1_mc.log" "log/run1_data.log" "log/run1_ext.log" "log/run1_dirt.log")
  echo " " > log/run1.log
  for i in "${arr[@]}"; do
    cat "$i" >> log/run1.log 
  done

  # Print the selection
  ./nuexsec --run 1 --printonly --printall --gpvm | tee -a log/run1.log 
  
  # Merge the files
  source merge/merge_run1_files.sh files/nuexsec_mc_run1.root files/nuexsec_run1_merged.root

  # Run the histogram plotter
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged.root --gpvm
  # ./nuexsec --run 1 --hist files/nuexsec_run1_merged.root --plotsys tot --gpvm

fi
# ---------------------

# Running slimmed down version of pelee ntuples with event weights
if [ "$1" == "weight" ]; then

  # Electron Energy
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel unisim  --xsec_smear mcc8 --gpvm
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsec_smear mcc8 --gpvm
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel genie   --xsec_smear mcc8 --gpvm
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel reint   --xsec_smear mcc8 --gpvm
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel mcstats --xsec_smear mcc8 --gpvm

  ./nuexsec --run 1 --sys reweight --gpvm

  # Electron Angle
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel unisim  --xsecvar elec_ang --xsec_smear mcc8 --gpvm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsecvar elec_ang --xsec_smear mcc8 --gpvm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel genie   --xsecvar elec_ang --xsec_smear mcc8 --gpvm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel reint   --xsecvar elec_ang --xsec_smear mcc8 --gpvm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel mcstats --xsecvar elec_ang --xsec_smear mcc8 --gpvm

  # ./nuexsec --run 1 --sys reweight --xsecvar elec_ang --gpvm

  # -- 

  # for running reweighting by cut -- these are slow, so dont run them by default for now
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel unisim --xsecplot rw_cuts --intrinsic intrinsic --gpvm
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel ppfx   --xsecplot rw_cuts --intrinsic intrinsic --gpvm
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel genie  --xsecplot rw_cuts --intrinsic intrinsic --gpvm
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay.root --xsecmode reweight --xseclabel reint  --xsecplot rw_cuts --intrinsic intrinsic --gpvm


fi


# Running detector variation samples

# To run the det var samples now do source run_selection_mcc9_run1.sh var <variation name>
# The variation are:
# CV, BNB_Diffusion, LYAttenuation, LYRayleigh, LYDown, SCE, Recomb2, WireModX, WireModYZ, WireModThetaXZ, WireModThetaYZ_withSigmaSplines, WireModThetaYZ_withoutSigmaSplines, WireModdEdX

# Could loop these to make it easier to run them all!

if [ "$1" == "var" ]; then
  ./nuexsec --run 1 --var /uboone/data/users/kmistry/work/MCC9/searchingfornues/ntuple_files_detvar_newtune/run1/extra_stats/neutrinoselection_filt_run1_overlay_$2.root $2 --gpvm
  
  source merge/merge_run1_files.sh files/nuexsec_mc_run1_$2.root files/nuexsec_run1_$2_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_$2_merged.root --var dummy $2 --gpvm

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
    source run_selection_mcc9_run1_gpvm.sh var "$i"
  done

fi
