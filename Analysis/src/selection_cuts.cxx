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
    
    if (SC.opfilter_pe_beam > 0) return true; // pass 
    else return false;               // fail
    
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
    if (SC.topological_score > 0.1) return true; // pass 
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
bool selection_cuts::dEdx_y(SliceContainer &SC){
    // if (SC.shr_dedx_Y_cali > 1.7 && SC.shr_dedx_Y_cali < 3.2) return true; // pass 
    // else return false;                                                     // fail

    // Kill the background region
    if ( (SC.shr_tkfit_dedx_Y >= 2.6 && SC.shr_tkfit_dedx_Y < 6.8) || (SC.shr_tkfit_dedx_Y < 0.5) ){
        return false;
    }
    // Try and recover the events in the other planes
    else {
        return true;
    }

}
// -----------------------------------------------------------------------------
bool selection_cuts::dEdx_v(SliceContainer &SC){

    // We only want to run this cut if the dedx calculation failed in the collection plane
    if (SC.shr_tkfit_dedx_Y > 0 ) return true;


    // Kill the background region
    if (  (SC.shr_tkfit_dedx_V > 2.7 && SC.shr_tkfit_dedx_V < 6 ) || (SC.shr_tkfit_dedx_V > 0 && SC.shr_tkfit_dedx_V < 1.5) ){
        return false;
    }
    // Try and recover the events in the other planes
    else {
        return true;
    }

}
// -----------------------------------------------------------------------------
bool selection_cuts::dEdx_u(SliceContainer &SC){

    // Signal region
    if (SC.shr_tkfit_dedx_U >= 1.7 && SC.shr_tkfit_dedx_U < 3.2){
        return true;
    }
    // failed the cut
    else {
        return false;
    }

}
// -----------------------------------------------------------------------------
bool selection_cuts::selected(SliceContainer &SC){
    if (SC.selected == 1 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_hits(SliceContainer &SC){
    if (SC.shr_hits_max > 220 ) return true; // pass 
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
    
    // if (SC.shr_distance < 10 ) return true; // pass 
    if (SC.shr_distance < 4) return true; // pass 
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
    if (SC.shrmoliereavg < 8 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_cyl_frac_1cm(SliceContainer &SC){
    if (SC.CylFrac2h_1cm > 0.01 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------