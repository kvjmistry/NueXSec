./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_CV_reco2_v08_00_00_38_run3b_reco2_reco2.root CV --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_LYRayleigh_v08_00_00_37_run3b_reco2_reco2.root LYRayleigh --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_LYAttenuation_v08_00_00_38_run3b_reco2_reco2.root LYAttenuation --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_SCE_reco2_v08_00_00_38_run3b_reco2_reco2.root SCE --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_LYDown_v08_00_00_37_v2_run3b_reco2_reco2.root LYDown --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_Recomb2_reco2_v08_00_00_39_run3b_reco2_reco2.root Recomb2 --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_wiremod_ScaleX_v08_00_00_38_run3b_reco2_reco2.root WireModX --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_wiremod_ScaleYZ_v08_00_00_38_run3b_reco2_reco2.root WireModYZ --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_WireModAngleXZ_v08_00_00_38_exe_run3b_reco2_reco2.root WireModThetaXZ --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_WireModAngleYZ_v08_00_00_38_exe_run3b_reco2_reco2.root WireModThetaYZ_withoutSigmaSplines --gpvm
./nuexsec --run 1 --var /uboone/data/users/davidc/searchingfornues/v08_00_00_44/0724/noweights/prodgenie_bnb_nu_overlay_DetVar_wiremod_ScaledEdX_v08_00_00_39_run3b_reco2_reco2.root WireModdEdX --gpvm

source merge/merge_run1_files.sh files/nuexsec_mc_run1_CV.root files/nuexsec_run1_CV_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_LYRayleigh.root files/nuexsec_run1_LYRayleigh_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_LYAttenuation.root files/nuexsec_run1_LYAttenuation_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_SCE.root files/nuexsec_run1_SCE_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_LYDown.root files/nuexsec_run1_LYDown_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_Recomb2.root files/nuexsec_run1_Recomb2_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_WireModX.root files/nuexsec_run1_WireModX_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_WireModYZ.root files/nuexsec_run1_WireModYZ_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_WireModThetaXZ.root files/nuexsec_run1_WireModThetaXZ_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_WireModThetaYZ_withoutSigmaSplines.root files/nuexsec_run1_WireModThetaYZ_withoutSigmaSplines_merged.root
source merge/merge_run1_files.sh files/nuexsec_mc_run1_WireModdEdX.root files/nuexsec_run1_WireModdEdX_merged.root

./nuexsec --run 1 --hist files/nuexsec_run1_CV_merged.root --var dummy CV --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_LYRayleigh_merged.root --var dummy LYRayleigh --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_LYAttenuation_merged.root --var dummy LYAttenuation --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_SCE_merged.root --var dummy SCE --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_LYDown_merged.root --var dummy LYDown --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_Recomb2_merged.root --var dummy Recomb2 --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_WireModX_merged.root --var dummy WireModX --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_WireModYZ_merged.root --var dummy WireModYZ --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_WireModThetaXZ_merged.root --var dummy WireModThetaXZ --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_WireModThetaYZ_withoutSigmaSplines_merged.root --var dummy WireModThetaYZ_withoutSigmaSplines --gpvm
./nuexsec --run 1 --hist files/nuexsec_run1_WireModdEdX_merged.root --var dummy WireModdEdX --gpvm


