#include "../include/SelectionCuts.h"
// -----------------------------------------------------------------------------
void SelectionCuts::Initalise(Utility _utility){

    std::cout << "Initalising Selection Cuts Class " << std::endl;
    _util = _utility;

}
// -----------------------------------------------------------------------------
bool SelectionCuts::swtrig(SliceContainer &SC, int type){
    
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
bool SelectionCuts::opfilt_pe(SliceContainer &SC, int type){
    
    // This should be turned on if BNB
    if (type == _util.k_data || type == _util.k_ext) {
        // Placeholder to remove a compilation warning
        // return true;
    }

    

    if (SC.opfilter_pe_beam >= 20) return true; // pass 
    else return false;               // fail
    
}
// -----------------------------------------------------------------------------
bool SelectionCuts::opfilt_veto(SliceContainer &SC, int type){

    // This should be turned on if BNB
    if (type == _util.k_data || type == _util.k_ext){
        // Placeholder to remove a compilation warning
        // return true;
    }

    if (SC.opfilter_pe_veto < 20) return true; // pass 
    else return false;               // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::slice_id(SliceContainer &SC){

    if (SC.nslice == 1) return true; // pass 
    else return false;               // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::e_candidate(SliceContainer &SC){
    if (SC.n_showers >= 1) return true; // pass 
    else return false;                 // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::topo_score(SliceContainer &SC){
    if (SC.topological_score > 0.2) return true; // pass 
    else return false;                            // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::in_fv(SliceContainer &SC){
    
    // Check if the reco vertex is in the FV
    bool is_in_fv = _util.in_fv(SC.reco_nu_vtx_sce_x, SC.reco_nu_vtx_sce_y, SC.reco_nu_vtx_sce_x);
    return is_in_fv;

}
// -----------------------------------------------------------------------------
bool SelectionCuts::cluster_frac(SliceContainer &SC){
    if (SC.slclustfrac > 0.4) return true; // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::shower_score(SliceContainer &SC){
    if (SC.shr_score < 0.15) return true; // pass 
    else return false;                   // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::michel_rej(SliceContainer &SC){
    if (SC.shr_energy_tot_cali > 0.075) return true; // pass 
    else return false;                               // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::dEdx_y(SliceContainer &SC){
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
bool SelectionCuts::dEdx_max_no_tracks(SliceContainer &SC){

    if (SC.n_tracks > 0) return true; // Dont apply this cut if there is no tracks

    double dedx_max = SC.GetdEdxMax();

    // Switch between Y and max plane
    // double dedx_var = SC.shr_tkfit_dedx_Y;
    double dedx_var = dedx_max;
    
    // Kill the background region
    if ( (dedx_var >= 2.7 && dedx_var < 5.7) || (dedx_var < 1.7) ){
        return false;
    }
    // Try and recover the events in the other planes
    else {
        return true;
    }

}
// -----------------------------------------------------------------------------
bool SelectionCuts::dEdx_v(SliceContainer &SC){

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
bool SelectionCuts::dEdx_u(SliceContainer &SC){

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
bool SelectionCuts::shr_hits(SliceContainer &SC){
    if (SC.shr_hits_max > 220 ) return true; // pass 
    else return false;                  // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::shr_distance(SliceContainer &SC){
    
    if (SC.n_tracks == 0) return true; // Dont apply this cut if there is no tracks
    
    if (SC.shr_distance < 6) return true; // pass 
    else return false;                      // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::shr_dist_dEdx_max(SliceContainer &SC){

    double dedx_max = SC.GetdEdxMax();

    // Switch between Y and max plane
    // double dedx_var = SC.shr_tkfit_dedx_Y;
    double dedx_var = dedx_max;

    if (SC.n_tracks == 0) return true; // Dont apply this cut if there is no tracks

    if (dedx_var < 0) return false;
    
    else if (dedx_var >= 0 && dedx_var < 0.5){
        return false;
    }
    
    else if (dedx_var >= 0.5 && dedx_var < 1.75){
        if (SC.shr_distance > 4 ) return false;
        else return true;
    }

    // For only y plane
    // else if (dedx_var >= 1.75 && dedx_var < 2.3){ 
    //     if (SC.shr_distance > 8 ) return false;
    //     else return true;
    // }
    // else if (dedx_var >= 2.3 && dedx_var < 3.5){
    //     if (SC.shr_distance > 3 ) return false;
    //     else return true;
    // }

    // for all planes
    else if (dedx_var >= 1.75 && dedx_var < 2.5){ 
        if (SC.shr_distance > 12 ) return false;
        else return true;
    }
    else if (dedx_var >= 2.5 && dedx_var < 3.5){
        if (SC.shr_distance > 3 ) return false;
        else return true;
    }


    else if (dedx_var >= 3.5 && dedx_var < 4.7){
        if (SC.shr_distance > 0 ) return false;
        else return true;
    }

    else if (dedx_var >= 4.7){
        if (SC.shr_distance > 3 ) return false;
        else return true;
    }
    else{
        std::cout << "Uncaught dEdx values..." << std::endl;
        return false;
    }

}
// -----------------------------------------------------------------------------
bool SelectionCuts::shr_hitratio(SliceContainer &SC){
    if (SC.hits_ratio > 0.5 ) return true; // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::shr_cosmic_IP(SliceContainer &SC){
    if (SC.CosmicIPAll3D > 15 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::contained_frac(SliceContainer &SC){
    if (SC.contained_fraction >= 0.85 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::shr_moliere_avg(SliceContainer &SC){
    if (SC.shrmoliereavg < 7 ) return true;    // pass 
    else return false;                     // fail
}
// -----------------------------------------------------------------------------
bool SelectionCuts::pi_zero_cuts(SliceContainer &SC){
    
    if (SC.pi0_shrscore1 < 0.5  &&
        SC.pi0_shrscore2 < 0.5  &&
        SC.pi0_dot1      > 0.8  &&
        SC.pi0_dot2      > 0.8  &&
        SC.pi0_radlen1   > 3.0  &&
        SC.pi0_radlen2   > 3.0  &&
        SC.pi0_gammadot  < 0.94 &&
        SC.pi0_energy1_Y > 60   &&
        SC.pi0_energy2_Y > 40   &&
        SC.pi0_dedx1_fit_Y > 1.0 ) {
        return true;
    }
    else return false;

}
// -----------------------------------------------------------------------------
bool SelectionCuts::numu_cuts(SliceContainer &SC){
    
    bool passed = false;

    // Loop over particle vector
    for (unsigned int l = 0; l < SC.trk_sce_start_x_v->size(); l++ ){
        
        // Track start, Track end contained
        bool trk_start = _util.in_fv(SC.trk_sce_start_x_v->at(l), SC.trk_sce_start_y_v->at(l), SC.trk_sce_start_z_v->at(l));
        bool trk_end   = _util.in_fv(SC.trk_sce_end_x_v->at(l),   SC.trk_sce_end_y_v->at(l),   SC.trk_sce_end_z_v->at(l));
    
        // PID Score
        bool pid_score = false;
        if (SC.trk_llr_pid_score_v->at(l) > 0.2) pid_score = true;

        // Track Length
        bool track_length = false;
        if (SC.trk_len_v->at(l) > 0.2) track_length = true;

        // Track Score
        bool track_score = false;
        if (SC.trk_score_v->at(l) > 0.8) track_score = true;

        // PFP Generation
        bool pfp_gen = false;
        if (SC.pfp_generation_v->at(l) > 0.2) pfp_gen = true;

        // Track Distance
        bool track_dist = false;
        if (SC.trk_distance_v->at(l) < 4) track_dist = true;

        // MCS Quality 
        bool b_mcs_quality = false;
        double mcs_quality = (SC.trk_mcs_muon_mom_v->at(l) - SC.trk_range_muon_mom_v->at(l)) / SC.trk_range_muon_mom_v->at(l);
        if (mcs_quality > -0.5 && mcs_quality < 0.5) b_mcs_quality = true;


        if (trk_start && trk_end && pid_score && track_length && track_score && pfp_gen && track_dist && b_mcs_quality) {
            passed = true;
            break;
        }
    }
    
    return passed;

}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
