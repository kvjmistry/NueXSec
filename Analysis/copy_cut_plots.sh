# Script to copy the cut plots to another directory for making presentations

# Currently works for run1, but will need to do this for other runs

# Make the directory
mkdir -p /uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/plots/run1/cut_plots/

# list of cuts minus the last one
cuts=(
  Unselected
  SoftwareTrig
  Op_Filter_PE
  Slice_ID
  e_candidate
  e_candidate
  e_candidate
  In_FV
  Topo_Score
  Cosmic_IP
  Cluster_Frac
  Shower_Score
  Shower_Contained
  Michel_Rej
  ShrHits
  HitRatio
  Moliere_Avg
  ShrVtxDistance
  dEdx_y
)

# the list of plots to copy
cuts2=(
reco_softwaretrig.pdf
reco_opfilter_beam.pdf
reco_nslice.pdf
reco_shower_multiplicity.pdf
reco_vtx_x_sce.pdf
reco_vtx_y_sce.pdf
reco_vtx_z_sce.pdf
reco_topological_score.pdf
reco_cosmicIP.pdf 
reco_slclustfrac.pdf
reco_shower_score.pdf
reco_n_shower_contained.pdf
reco_shower_energy_tot_cali.pdf
reco_shr_hits_max.pdf
reco_hits_ratio.pdf
reco_shrmoliereavg.pdf
reco_shower_to_vtx_dist.pdf
reco_shr_tkfit_dedx_y.pdf
)

cut_path=/uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/plots/run1/cuts/
copy_path=/uboone/app/users/kmistry/MCC9_uboonecode_v08_00_00_33/srcs/ubana/ubana/NueXSec/Analysis/plots/run1/cut_plots/

# Do the copy
for i in ${!cuts[*]}; do 
    echo "cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path"
    cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path
done




