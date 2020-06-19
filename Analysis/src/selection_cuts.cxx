#include "../include/selection_cuts.h"
// -----------------------------------------------------------------------------
void selection_cuts::Initalise(utility _utility){

    std::cout << "Initalising Selection Cuts Class " << std::endl;
    _util = _utility;

}
// -----------------------------------------------------------------------------
bool selection_cuts::swtrig(SliceContainer &SC, int type){
    
    // Common optical is already applied to data
    if (type == _util.k_mc || type == _util.k_dirt){
        
        // Run 1 uses a different software trigger
        if (SC.run_period == "1"){
            if (SC.swtrig_pre == 1) return true;  // pass 
            else return false;                    // fail
        }
        // Keep using standard software trigger for run 3 for now...
        else {
            if (SC.swtrig == 1) return true;  // pass 
            else return false;               // fail
        }
       
    }
    else return true;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::opfilt_pe(SliceContainer &SC, int type){
    
    // This should be turned on if BNB
    // if (type == _util.k_data || type == _util.k_ext) return true;

    if (SC.opfilter_pe_beam > 0) return true; // pass 
    else return false;               // fail
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::opfilt_veto(SliceContainer &SC, int type){

    // This should be turned on if BNB
    // if (type == _util.k_data || type == _util.k_ext) return true;

    if (SC.opfilter_pe_veto < 20) return true; // pass 
    else return false;               // fail
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
    if (SC.topological_score > 0.2) return true; // pass 
    else return false;                            // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::in_fv(SliceContainer &SC){
    
    // Check if the reco vertex is in the FV
    bool is_in_fv = _util.in_fv(SC.reco_nu_vtx_sce_x, SC.reco_nu_vtx_sce_y, SC.reco_nu_vtx_sce_x);
    return is_in_fv;

}
// -----------------------------------------------------------------------------
bool selection_cuts::cluster_frac(SliceContainer &SC){
    if (SC.slclustfrac > 0.4) return true; // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shower_score(SliceContainer &SC){
    if (SC.shr_score < 0.15) return true; // pass 
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
    if ( (SC.shr_tkfit_dedx_Y >= 3.1 && SC.shr_tkfit_dedx_Y < 5.5) || (SC.shr_tkfit_dedx_Y < 0.5) ){
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
bool selection_cuts::shr_hits(SliceContainer &SC){
    if (SC.shr_hits_max > 220 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_distance(SliceContainer &SC){
    
    if (SC.n_tracks == 0) return true; // Dont apply this cut if there is no tracks
    
    // if (SC.shr_distance < 10 ) return true; // pass 
    if (SC.shr_distance < 6) return true; // pass 
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
bool selection_cuts::contained_frac(SliceContainer &SC){
    if (SC.contained_fraction >= 0.75 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool selection_cuts::shr_moliere_avg(SliceContainer &SC){
    if (SC.shrmoliereavg < 7 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
