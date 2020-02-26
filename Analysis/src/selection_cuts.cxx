#include "../include/selection_cuts.h"
// -----------------------------------------------------------------------------
bool selection_cuts::slice_id(SliceContainer &SC){
    if (SC.nslice == 1) return true; // pass 
    else return false;               // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::e_candidate(SliceContainer &SC){
    if (SC.n_showers > 0) return true; // pass 
    else return false;                 // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::topo_score(SliceContainer &SC){
    if (SC.topological_score < 0.15) return true; // pass 
    else return false;                            // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::in_fv(SliceContainer &SC){
    
    // These are looser cuts (whats done in the cc inclusive)
    if ( SC.reco_nu_vtx_sce_x >= 10 && SC.reco_nu_vtx_sce_x <= 246.35 &&
        SC.reco_nu_vtx_sce_y >= -106.5 && SC.reco_nu_vtx_sce_y <= 106.5 &&
        ( (SC.reco_nu_vtx_sce_z >= 20 && SC.reco_nu_vtx_sce_z <= 986.8))
        ){
        return true;   // pass
    }  
    else return false; // fail
    
    // These are tighter cuts
    // if ( SC.reco_nu_vtx_sce_x >= 22 && SC.reco_nu_vtx_sce_x <= 234.35 &&
    //     SC.reco_nu_vtx_sce_y >= -75.1 && SC.reco_nu_vtx_sce_y <= 75.1 &&
    //     ( (SC.reco_nu_vtx_sce_z >= 35 && SC.reco_nu_vtx_sce_z <= 665) || (SC.reco_nu_vtx_sce_z >= 785 && SC.reco_nu_vtx_sce_z <= 941.8))
    //     ){
    //     return true; // pass
    // }  
    // else return false; // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::cluster_frac(SliceContainer &SC){
    if (SC.slclustfrac > 0.4) return true; // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shower_score(SliceContainer &SC){
    if (SC.shr_score < 0.3) return true; // pass 
    else return false;                   // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::michel_rej(SliceContainer &SC){
    if (SC.shr_energy_tot_cali > 0.075) return true; // pass 
    else return false;                               // fail
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------