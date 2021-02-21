# Script with commands to run the entire selection

# Set the smearing mode er = Event Rate, MCC8 = Marco
sm=er
# sm=mcc8

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
  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "")'

  # Now run the cross section calculator
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode default --xsec_smear $sm | tee -a log/run1.log 
fi
# ---------------------

# Running slimmed down version of pelee ntuples with event weights
if [ "$1" == "weight" ]; then

  # Electron Energy
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel unisim  --xsec_smear $sm --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsec_smear $sm --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel genie   --xsec_smear $sm --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel reint   --xsec_smear $sm --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel mcstats --xsec_smear $sm --xsecbins standard --xsecvar elec_ang

  # Re-run the detvar xsec
  source run_selection_mcc9_run1.sh allvarxsec

  ./nuexsec --run 1 --sys reweight --xsec_smear wiener --binscaling standard --xsecvar elec_ang
  ./nuexsec --run 1 --sys reweight --xsec_smear $sm --binscaling width --xsecvar elec_ang

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_mec.root       --var dummy mec       --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nogtune.root   --var dummy nogtune   --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang 
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nopi0tune.root --var dummy nopi0tune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_FLUGG.root     --var dummy FLUGG     --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang 
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_tune1.root     --var dummy tune1     --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang

  # Electron Angle
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel unisim  --xsecvar elec_ang --xsec_smear $sm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsecvar elec_ang --xsec_smear $sm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel genie   --xsecvar elec_ang --xsec_smear $sm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel reint   --xsecvar elec_ang --xsec_smear $sm
  # ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel mcstats --xsecvar elec_ang --xsec_smear $sm

  # ./nuexsec --run 1 --sys reweight --xsecvar elec_ang --xsec_smear $sm

  # -- 

  # for running reweighting by cut -- these are slow, so dont run them by default for now
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --xsecmode reweight --xseclabel unisim --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --xsecmode reweight --xseclabel ppfx   --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --xsecmode reweight --xseclabel genie  --xsecplot rw_cuts --intrinsic intrinsic
  #./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --xsecmode reweight --xseclabel reint  --xsecplot rw_cuts --intrinsic intrinsic

fi


# Running detector variation samples

# To run the det var samples now do source run_selection_mcc9_run1.sh var <variation name>
# The variation are:
# CV, BNB_Diffusion, LYAttenuation, LYRayleigh, LYDown, SCE, Recomb2, WireModX, WireModYZ, WireModThetaXZ, WireModThetaYZ_withSigmaSplines, WireModThetaYZ_withoutSigmaSplines, WireModdEdX

# Could loop these to make it easier to run them all!

if [ "$1" == "var" ]; then
  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_$2.root $2
  
  # Overwrite the true nue information
  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_$2_intrinsic.root $2 --intrinsic intrinsic

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_$2.root files/nuexsec_run1_$2_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_$2_merged.root --var dummy $2

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_'"$2"'.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "'"$2"'")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_$2.root --var dummy $2 --xsecmode default --xsecvar elec_E --xsec_smear $sm

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

if [ "$1" == "allvarxsec" ]; then

  # run the above script for every det variation 
  for i in "${var[@]}"
  do
    ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_$i.root --var dummy $i --xsecmode default --xsecvar elec_ang --xsec_smear $sm --xsecbins standard
  done

fi

# Make a new file with a reweighted MEC model
if [ "$1" == "mec" ]; then

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root mec --tunemec

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root mec --intrinsic intrinsic --tunemec 

  # source merge/merge_run1_files.sh files/nuexsec_mc_run1_mec.root files/nuexsec_run1_mec_merged.root

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_mec.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "mec")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_mec.root --var dummy mec --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_mec.root --var dummy mec --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_mec.root --var dummy mec --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Make a new file with a reweighted MEC model
if [ "$1" == "nogtune" ]; then

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root nogtune --weight_tune 0

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root nogtune --intrinsic intrinsic --weight_tune 0

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_nogtune.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "nogtune")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nogtune.root --var dummy nogtune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nogtune.root --var dummy nogtune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nogtune.root --var dummy nogtune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Make a new file with a reweighted no pi0 tune model
if [ "$1" == "nopi0tune" ]; then

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root nopi0tune --weight_pi0 0

  ./nuexsec --run 1 --var ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root nopi0tune --intrinsic intrinsic --weight_pi0 0

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_nopi0tune.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "nopi0tune")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nopi0tune.root --var dummy nopi0tune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nopi0tune.root --var dummy nopi0tune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_nopi0tune.root --var dummy nopi0tune --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Make a new file with a reweighted no pi0 tune model
if [ "$1" == "FLUGG" ]; then

  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_FLUGG.root FLUGG --weight_ppfx 0

  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_FLUGG_intrinsic.root FLUGG --intrinsic intrinsic --weight_ppfx 0

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_FLUGG.root files/nuexsec_run1_FLUGG_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_FLUGG_merged.root --var dummy FLUGG

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_FLUGG.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "FLUGG")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_FLUGG.root --var dummy FLUGG --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_FLUGG.root --var dummy FLUGG --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_FLUGG.root --var dummy FLUGG --xsecmode default --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Make a new file with genie tune1 model
if [ "$1" == "tune1" ]; then

  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_tune1.root tune1 --weight_tune 0

  ./nuexsec --run 1 --var ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_tune1_intrinsic.root tune1 --intrinsic intrinsic --weight_tune 0

  source merge/merge_run1_files.sh files/nuexsec_mc_run1_tune1.root files/nuexsec_run1_tune1_merged.root

  ./nuexsec --run 1 --hist files/nuexsec_run1_tune1_merged.root --var dummy tune1

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, false, "files/trees/nuexsec_selected_tree_mc_run1_tune1.root", "files/trees/nuexsec_selected_tree_data_run1.root", "files/trees/nuexsec_selected_tree_ext_run1.root","files/trees/nuexsec_selected_tree_dirt_run1.root", "tune1")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_tune1.root --var dummy tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_tune1.root --var dummy tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_tune1.root --var dummy tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi


# To generate the file lists with all event weights then use this command
# ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode txtlist --xseclabel all --xsecvar elec_E

# Use the reweighted mec model as fake data
if [ "$1" == "fakemec" ]; then

  # Run the cross sec calculation in fake data mode
  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --fake mec --tunemec 

  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root --fake mec --intrinsic intrinsic --tunemec 

  # Merge the files
  hadd -f -T files/nuexsec_run1_merged_fakemec.root files/nuexsec_data_run1_mec.root files/nuexsec_dirt_run1_fake.root files/nuexsec_ext_run1_fake.root files/nuexsec_mc_run1.root

  # Plot the files
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged_fakemec.root --fake mec

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1_mec.root", "files/trees/nuexsec_selected_tree_ext_run1_fake.root","files/trees/nuexsec_selected_tree_dirt_run1_fake.root", "fakemec")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakemec.root --fake mec --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakemec.root --fake mec --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakemec.root --fake mec --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Use the reweighted no genie tune model as fake data
if [ "$1" == "fakenogtune" ]; then

  # Run the cross sec calculation in fake data mode
  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --fake nogtune --weight_tune 0

  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root --fake nogtune --intrinsic intrinsic --weight_tune 0

  # Merge the files
  hadd -f -T files/nuexsec_run1_merged_fakenogtune.root files/nuexsec_data_run1_nogtune.root files/nuexsec_dirt_run1_fake.root files/nuexsec_ext_run1_fake.root files/nuexsec_mc_run1.root

  # Plot the files
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged_fakenogtune.root --fake nogtune

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1_nogtune.root", "files/trees/nuexsec_selected_tree_ext_run1_fake.root","files/trees/nuexsec_selected_tree_dirt_run1_fake.root", "fakenogtune")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenogtune.root --fake nogtune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenogtune.root --fake nogtune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenogtune.root --fake nogtune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Use the reweighted no pi0 tune model as fake data
if [ "$1" == "fakenopi0tune" ]; then

  # Run the cross sec calculation in fake data mode
  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_newtune.root --fake nopi0tune --weight_pi0 0

  ./nuexsec --run 1 --mc ../ntuples/neutrinoselection_filt_run1_overlay_intrinsic_newtune.root --fake nopi0tune --intrinsic intrinsic --weight_pi0 0

  # Merge the files
  hadd -f -T files/nuexsec_run1_merged_fakenopi0tune.root files/nuexsec_data_run1_nopi0tune.root files/nuexsec_dirt_run1_fake.root files/nuexsec_ext_run1_fake.root files/nuexsec_mc_run1.root

  # Plot the files
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged_fakenopi0tune.root --fake nopi0tune

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1_nopi0tune.root", "files/trees/nuexsec_selected_tree_ext_run1_fake.root","files/trees/nuexsec_selected_tree_dirt_run1_fake.root", "fakenopi0tune")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenopi0tune.root --fake nopi0tune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenopi0tune.root --fake nopi0tune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakenopi0tune.root --fake nopi0tune --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi



# Use the tune 1 model as fake data
if [ "$1" == "faketune1" ]; then

  # Run the cross sec calculation in fake data mode
  ./nuexsec --run 1 --mc ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_tune1.root --fake tune1 --weight_tune 0

  ./nuexsec --run 1 --mc ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_tune1_intrinsic.root --fake tune1 --intrinsic intrinsic --weight_tune 0

  # Merge the files
  hadd -f -T files/nuexsec_run1_merged_faketune1.root files/nuexsec_data_run1_tune1.root files/nuexsec_dirt_run1_fake.root files/nuexsec_ext_run1_fake.root files/nuexsec_mc_run1.root

  # Plot the files
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged_faketune1.root --fake tune1

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1_tune1.root", "files/trees/nuexsec_selected_tree_ext_run1_fake.root","files/trees/nuexsec_selected_tree_dirt_run1_fake.root", "faketune1")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_faketune1.root --fake tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_faketune1.root --fake tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_faketune1.root --fake tune1 --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi

# Use the tune 1 model as fake data
if [ "$1" == "fakeFLUGG" ]; then

  # Run the cross sec calculation in fake data mode
  ./nuexsec --run 1 --mc ../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_FLUGG.root --fake FLUGG  --weight_ppfx 0

  ./nuexsec --run 1 --mc ../ntuples/detvar_newtune/run1/intrinsic/neutrinoselection_filt_run1_overlay_FLUGG_intrinsic.root --fake FLUGG --intrinsic intrinsic --weight_ppfx 0

  # Merge the files
  hadd -f -T files/nuexsec_run1_merged_fakeFLUGG.root files/nuexsec_data_run1_FLUGG.root files/nuexsec_dirt_run1_fake.root files/nuexsec_ext_run1_fake.root files/nuexsec_mc_run1.root

  # Plot the files
  ./nuexsec --run 1 --hist files/nuexsec_run1_merged_fakeFLUGG.root --fake FLUGG

  root -l -b -q 'merge/merge_uneaventrees.C("1", true, true, "files/trees/nuexsec_selected_tree_mc_run1.root", "files/trees/nuexsec_selected_tree_data_run1_FLUGG.root", "files/trees/nuexsec_selected_tree_ext_run1_fake.root","files/trees/nuexsec_selected_tree_dirt_run1_fake.root", "fakeFLUGG")'

  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakeFLUGG.root --fake FLUGG --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_E
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakeFLUGG.root --fake FLUGG --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_ang
  ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1_fakeFLUGG.root --fake FLUGG --xsecmode default  --xsec_smear er --xsecbins standard --xsecvar elec_cang

fi