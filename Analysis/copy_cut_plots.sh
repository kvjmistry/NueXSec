# Script to copy the cut plots to another directory for making presentations
# basically its saving the plots before we apply the cut to it

# Change this here to the path of where your NueXSec folder lives
USER_PATH=/Users/kvjmistry/work/


# Make the directory
mkdir -p $USER_PATH/NueXSec/Analysis/plots/run1/cut_plots/
mkdir -p $USER_PATH/NueXSec/Analysis/plots/run3/cut_plots/

# list of cuts minus the last one
cuts=(
  Unselected
  SoftwareTrig
  Slice_ID
  e_candidate
  e_candidate
  e_candidate
  In_FV
  Contained_Frac
  Contained_Frac
  Topo_Score
  Cosmic_IP
  Shower_Score
  HitRatio
  Moliere_Avg
  Moliere_Avg
  Moliere_Avg
  ShrVtxDist_dEdx_max
)

# the list of plots to copy
cuts2=(
reco_softwaretrig.pdf
reco_nslice_logy.pdf
reco_shower_multiplicity_logy.pdf
reco_vtx_x_sce.pdf
reco_vtx_y_sce.pdf
reco_vtx_z_sce.pdf
reco_contained_fraction_logy.pdf
reco_topological_score.pdf
reco_topological_score_logy.pdf
reco_CosmicIPAll3D.pdf
reco_shower_score.pdf
reco_hits_ratio.pdf
reco_shrmoliereavg.pdf
reco_shower_to_vtx_dist.pdf
reco_shr_tkfit_dedx_max.pdf
reco_shr_tkfit_dedx_max_with_tracks.pdf
reco_shr_tkfit_dedx_max_no_tracks.pdf
)

cut_path=$USER_PATH/NueXSec/Analysis/plots/run1/cuts/
copy_path=$USER_PATH/NueXSec/Analysis/plots/run1/cut_plots/

# Do the copy
for i in ${!cuts[*]}; do 
    echo "cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path"
    cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path
done

cut_path=$USER_PATH/NueXSec/Analysis/plots/run3/cuts/
copy_path=$USER_PATH/NueXSec/Analysis/plots/run3/cut_plots/

# Do the copy
for i in ${!cuts[*]}; do 
    echo "cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path"
    cp $cut_path${cuts[$i]}/${cuts2[$i]} $copy_path
done




