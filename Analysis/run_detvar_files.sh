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
	WireModThetaYZ_withoutSigmaSplines
	WireModdEdX
)

for i in "${var[@]}"
do
	source merge/merge_run1_files.sh files/nuexsec_mc_run1_"$i".root files/nuexsec_run1_"$i"_merged.root
	./nuexsec --run 1 --hist files/nuexsec_run1_"$i"_merged.root --var dummy "$i" --gpvm
	echo "============================================================================================================================================================="
done
