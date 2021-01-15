# Script to copy all the beamline plots we want to put into the technote into a common folder
# this will make it easier to drag and drop them into the overleaf folder

# Change this here to the path of where your NueXSec folder lives
USER_PATH=/Users/kvjmistry/work/

file_path=$USER_PATH/NueXSec/Analysis/plots/run1/Systematics/
copy_path=$USER_PATH/NueXSec/Analysis/plots/run1/Systematics/Beamline

# Make the directory
mkdir -p $copy_path

# list of cuts minus the last one
vars=(
  Horn_curr
  Horn1_x
  Horn1_y
  Beam_spot
  Horn2_x
  Horn2_y
  Horn_Water
  Beam_shift_x
  Beam_shift_y
  Target_z
  Horn1_refined_descr
  Decay_pipe_Bfield
  Old_Horn_Geometry
)


# Do the copy
for i in ${!vars[*]}; do 
    echo "cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path
    echo "cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path
    echo "cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path"
    cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path
done


copy_path=$USER_PATH/NueXSec/Analysis/plots/run1/Systematics/Genie_Unisim

# Make the directory
mkdir -p $copy_path

# list of cuts minus the last one
vars=(
  RPA
  CCMEC
  AxFFCCQE
  VecFFCCQE
  DecayAngMEC
  ThetaDelta2Npi
  ThetaDelta2NRad
  RPA_CCQE_Reduced
  NormCCCOH
  NormNCCOH
)

# Do the copy
for i in ${!vars[*]}; do 
    echo "cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path
    echo "cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path
    echo "cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path"
    cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path
done

copy_path=$USER_PATH/NueXSec/Analysis/plots/run1/Systematics/DetVar

# Make the directory
mkdir -p $copy_path

# list of cuts minus the last one
vars=(
  LYRayleigh
  LYAttenuation
  SCE
  Recomb2
  WireModX
  WireModYZ
  WireModThetaXZ
  WireModThetaYZ_withSigmaSplines
  WireModThetaYZ_withoutSigmaSplines
  WireModdEdX
)

# Do the copy
for i in ${!vars[*]}; do 
    echo "cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*bkg* $copy_path
    echo "cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path"
    cp $file_path${vars[$i]}/reco_el_E/*mc_xsec* $copy_path
    echo "cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path"
    cp $file_path${vars[$i]}/true_el_E/*eff* $copy_path
done






