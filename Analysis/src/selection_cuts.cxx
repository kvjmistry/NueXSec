#include "../include/selection_cuts.h"
// -----------------------------------------------------------------------------
bool selection_cuts::swtrig(SliceContainer &SC, int type){
    
    // Common optical is already applied to data
    if (type == _util.k_mc || type == _util.k_dirt){
        if (SC.swtrig > 0) return true;  // pass 
        else return false;               // fail
    }
    else return true;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::opfilt_pe(SliceContainer &SC, int type){
    
    // Common optical is already applied to data
    if (type == _util.k_mc || type == _util.k_dirt){
        if (SC.opfilter_pe_beam > 0) return true; // pass 
        else return false;               // fail
    }
    else return true;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::opfilt_veto(SliceContainer &SC, int type){

    // Common optical is already applied to data
    if (type == _util.k_mc || type == _util.k_dirt){
        if (SC.opfilter_pe_veto < 20) return true; // pass 
        else return false;               // fail
    }
    else return true;
}
// -----------------------------------------------------------------------------
bool selection_cuts::slice_id(SliceContainer &SC){

    if (SC.nslice == 1) return true; // pass 
    else return false;               // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::e_candidate(SliceContainer &SC){
    if (SC.n_showers >= 1) return true; // pass 
    else return false;                 // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::topo_score(SliceContainer &SC){
    if (SC.topological_score > 0.15) return true; // pass 
    else return false;                            // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::in_fv(SliceContainer &SC){
    
    // These are looser cuts (whats done in the cc inclusive)
    if ( SC.reco_nu_vtx_sce_x >= 10 && SC.reco_nu_vtx_sce_x <= 243 &&
        SC.reco_nu_vtx_sce_y >= -106.5 && SC.reco_nu_vtx_sce_y <= 106.5 &&
        ( (SC.reco_nu_vtx_sce_z >= 20 && SC.reco_nu_vtx_sce_z <= 986))
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
bool selection_cuts::dEdx(SliceContainer &SC){
    if (SC.shr_dedx_Y_cali > 1.7 && SC.shr_dedx_Y_cali < 3.2) return true; // pass 
    else return false;                                                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::selected(SliceContainer &SC){
    if (SC.selected == 1 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_hits(SliceContainer &SC){
    if (SC.shr_hits_tot > 220 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_hits_y_plane(SliceContainer &SC){
    if (SC.shr_hits_y_tot > 200 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_distance(SliceContainer &SC){
    
    if (SC.n_tracks == 0) return true; // Dont apply this cut if there is no tracks
    
    if (SC.shr_distance < 10 ) return true; // pass 
    else return false;                      // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_hitratio(SliceContainer &SC){
    if (SC.hits_ratio > 0.7 ) return true; // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_cosmic_IP(SliceContainer &SC){
    if (SC.CosmicIP > 20 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_contained(SliceContainer &SC){
    if (SC.n_showers_contained >= 1 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_moliere_avg(SliceContainer &SC){
    if (SC.shrmoliereavg < 10 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------