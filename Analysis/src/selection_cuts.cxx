#include "../include/selection_cuts.h"
// -----------------------------------------------------------------------------
void selection_cuts::SetFlashVariables(std::vector<double> largest_flash_v){
        largest_flash_y 	= largest_flash_v.at(0);
        largest_flash_z 	= largest_flash_v.at(1);
        largest_flash_time 	= largest_flash_v.at(2);
        largest_flash_pe 	= largest_flash_v.at(3);
}
// -----------------------------------------------------------------------------
void selection_cuts::SetTPCObjVariables(xsecAna::TPCObjectContainer tpc_obj, 
                                        double mc_nu_vtx_x, double mc_nu_vtx_y, double mc_nu_vtx_z,
                                        std::vector<double> fv_boundary_v, bool has_pi0){
    
    // Check if the true vtx is in the FV
    true_in_tpc = in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);

    // Classify the event  
    tpc_classification = TPCO_Classifier(tpc_obj, true_in_tpc, has_pi0);

    // TPC Obj vars
    tpc_obj_vtx_x        = tpc_obj.pfpVtxX();
    tpc_obj_vtx_y        = tpc_obj.pfpVtxY();
    tpc_obj_vtx_z        = tpc_obj.pfpVtxZ();
    tpc_obj_mode         = tpc_obj.Mode(); 
    n_pfp                = tpc_obj.NumPFParticles();
    leading_shower_index = tpc_classification.second;

}
// -----------------------------------------------------------------------------
std::pair<std::string, int> selection_cuts::TPCO_Classifier(xsecAna::TPCObjectContainer tpc_obj, bool true_in_tpc, bool has_pi0) {
    int part_nue_cc     = 0;
    int part_nue_bar_cc = 0;
    int part_cosmic     = 0;
    int part_nc         = 0;
    int part_nc_pi0     = 0;
    int part_numu_cc    = 0;
    int part_unmatched  = 0;

    const int tpc_obj_mode = tpc_obj.Mode();
    const int n_pfp = tpc_obj.NumPFParticles();
    int most_hits = 0;
    int leading_index = -1;
    std::string leading_origin = "kNothing";

    for(int j = 0; j < n_pfp; j++)
    {
        auto const part = tpc_obj.GetParticle(j);
        const int n_pfp_hits = part.NumPFPHits();
        const int mc_parent_pdg = part.MCParentPdg();
        const int pfp_pdg = part.PFParticlePdgCode();
        if(pfp_pdg == 11)
        {
            if(n_pfp_hits > most_hits)
            {
                leading_index = j;
                most_hits = n_pfp_hits;
            }
        }
        if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == 12)  { part_nue_cc++; }
        if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && mc_parent_pdg == -12) { part_nue_bar_cc++; }
        if(part.CCNC() == 0 && part.Origin() == "kBeamNeutrino" && (mc_parent_pdg == 14 || mc_parent_pdg == -14)) { part_numu_cc++; }
        if(part.CCNC() == 1 && part.Origin() == "kBeamNeutrino")
        {
            if(has_pi0 == true)  {part_nc_pi0++; }
            if(has_pi0 == false) {part_nc++; }
        }
        if(part.Origin() == "kCosmicRay") { part_cosmic++;    }
        if(part.Origin() == "kUnknown"  ) { part_unmatched++; }
    }
    //some tpc objects actually have 0 hits - crazy!
    if(tpc_obj.NumPFPHits() == 0) {return std::make_pair("bad_reco", 0); }

    //currently, any tpc objects which only have a track end up with a leading_index of -1
    //this index will likely cause code to crash if called before the signal definition cuts

    //also some rare cases where nu_pfp = nue, and shower hits = 0 with track hits > 0 - how does this happen? (NC event?)

    //now to catagorise the tpco
    if(part_cosmic > 0)
    {
        if(part_nue_cc  > 0 || part_nue_bar_cc > 0)  { return std::make_pair("nue_cc_mixed",  leading_index); }
        if(part_numu_cc > 0 )                        { return std::make_pair("numu_cc_mixed", leading_index); }
        if(part_nc  > 0 || part_nc_pi0 > 0)          { return std::make_pair("other_mixed",   leading_index); }
        return std::make_pair("cosmic", leading_index);
    }
    //this uses the true neutrino vertex for this specific event
    //not the true vtx per tpc object - maybe this can be fixed in the future...
    //but using the true nu vtx only matters for the pure signal events,
    //where the neutrino vertex IS the true tpc object vertex
    if(part_cosmic == 0)
    {
        if(part_nue_cc      > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }
        if(part_nue_bar_cc  > 0 && true_in_tpc == false) { return std::make_pair("nue_cc_out_fv", leading_index);   }

        if(part_nue_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_cc_qe",     leading_index);   }
        if(part_nue_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_cc_res",    leading_index);   }
        if(part_nue_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_cc_dis",    leading_index);   }
        if(part_nue_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_cc_coh",    leading_index);   }
        if(part_nue_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_cc_mec",    leading_index);   }

        if(part_nue_bar_cc    > 0 && tpc_obj_mode == 0   ) { return std::make_pair("nue_bar_cc_qe",     leading_index);   }
        if(part_nue_bar_cc    > 0 && tpc_obj_mode == 1   ) { return std::make_pair("nue_bar_cc_res",    leading_index);   }
        if(part_nue_bar_cc    > 0 && tpc_obj_mode == 2   ) { return std::make_pair("nue_bar_cc_dis",    leading_index);   }
        if(part_nue_bar_cc    > 0 && tpc_obj_mode == 3   ) { return std::make_pair("nue_bar_cc_coh",    leading_index);   }
        if(part_nue_bar_cc    > 0 && tpc_obj_mode == 10  ) { return std::make_pair("nue_bar_cc_mec",    leading_index);   }

        if(part_numu_cc     > 0 && tpc_obj_mode == 0   ) { return std::make_pair("numu_cc_qe",    leading_index);   }
        if(part_numu_cc     > 0 && tpc_obj_mode == 1   ) { return std::make_pair("numu_cc_res",   leading_index);   }
        if(part_numu_cc     > 0 && tpc_obj_mode == 2   ) { return std::make_pair("numu_cc_dis",   leading_index);   }
        if(part_numu_cc     > 0 && tpc_obj_mode == 3   ) { return std::make_pair("numu_cc_coh",   leading_index);   }
        if(part_numu_cc     > 0 && tpc_obj_mode == 10  ) { return std::make_pair("numu_cc_mec",   leading_index);   }
        if(part_nc          > 0                        ) { return std::make_pair("nc",            leading_index);   }
        if(part_nc_pi0      > 0                        ) { return std::make_pair("nc_pi0",        leading_index);   }
        if(part_unmatched   > 0                        ) { return std::make_pair("unmatched",     leading_index);   }
    }
    //this never happens :)
    std::cout << "HELP HELP HELP END OF TPCO CLASSIFIER AND NO CLASSIFICATION!" << std::endl;
    //return the string for the tpco id
    return std::make_pair("non_match", 0);
}
// *****************************************************************************
// ----------------------- Selection Cuts Functions ----------------------------
// *****************************************************************************

// -----------------------------------------------------------------------------
bool selection_cuts::FlashinTime_FlashPE(double flash_time_start, double flash_time_end, double flash_pe_threshold, std::vector<double> opt_time_v, std::vector<int> opt_pe_v){

    bool in_time          = false;
    bool sufficient_flash = false;
    
    // Loop over the optical list vec
    for (unsigned int j = 0; j < opt_pe_v.size(); j++) {
        
        auto opt_time = opt_time_v.at(j) + 1.0;
        auto opt_pe   = opt_pe_v.at(j);
        
        // See if flash was in time
        in_time = (opt_time >= flash_time_start && opt_time <= flash_time_end) ? true : false;
        
        // See if flash meets the threshold requirements
        sufficient_flash = (opt_pe >= flash_pe_threshold) ? true : false;
        
        // Flash is both in time and over PE threshold
        if(in_time == true && sufficient_flash == true){
            return true;
            break; // once pased we are done, so dont loop any more otherwise we may overwrite this

        }
    }

    return false;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::HasNue(xsecAna::TPCObjectContainer tpc_obj) {

    bool has_nue = false;
    bool has_valid_shower = false;

    // Loop over the PFP
    for (int j = 0; j < n_pfp; j++) {
        auto const part     = tpc_obj.GetParticle(j);
        const int  pfp_pdg  = part.PFParticlePdgCode();
        const int  pfp_hits = part.NumPFPHits();
        
        if(pfp_pdg == 11 && pfp_hits > 0) has_valid_shower = true; 
        
        if(pfp_pdg == 12) has_nue = true; 
    }

    if(has_nue == true && has_valid_shower == true)
        return true; 
    else 
        return false;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::in_fv(double x, double y, double z, std::vector<double> fv_boundary_v){

    const double det_x1 = 0;
    const double det_x2 = 256.35;
    const double det_y1 = -116.5;
    const double det_y2 = 116.5;
    const double det_z1 = 0;
    const double det_z2 = 1036.8;

    const double x1 = fv_boundary_v.at(0);
    const double x2 = fv_boundary_v.at(1);
    const double y1 = fv_boundary_v.at(2);
    const double y2 = fv_boundary_v.at(3);
    const double z1 = fv_boundary_v.at(4);
    const double z2 = fv_boundary_v.at(5);

    if(x <= det_x1 + x1 || x >= det_x2 - x2) {return false; }
    if(y <= det_y1 + y1 || y >= det_y2 - y2) {return false; }
    if(z <= det_z1 + z1 || z >= det_z2 - z2) {return false; }
    return true;
}
// -----------------------------------------------------------------------------
// Flash reco Vertex Distance 
bool selection_cuts::opt_vtx_distance(double tpc_vtx_y, double tpc_vtx_z, double flash_vtx_y, double flash_vtx_z, double tolerance) {
    const double distance = sqrt(pow((tpc_vtx_y - flash_vtx_y), 2) + pow((tpc_vtx_z - flash_vtx_z), 2) );
    
    if(distance <= tolerance) return true;
    return false;
}
// -----------------------------------------------------------------------------
bool selection_cuts::flashRecoVtxDist(std::vector< double > largest_flash_v, double tolerance, const double tpc_vtx_x, const double tpc_vtx_y, const double tpc_vtx_z){
    bool is_close;
    
    // Flash is upstream
    if(tpc_vtx_z < largest_flash_v.at(1)) 
        is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), tolerance);
    
    // Flash is downstream
    if(tpc_vtx_z >= largest_flash_v.at(1)) 
        is_close = opt_vtx_distance(tpc_vtx_y, tpc_vtx_z, largest_flash_v.at(0), largest_flash_v.at(1), (tolerance - 20));
    
    if (is_close == true )
        return true;
    else
        return false;
}
// -----------------------------------------------------------------------------
bool selection_cuts::VtxNuDistance(xsecAna::TPCObjectContainer tpc_obj,int pfp_pdg_type , double tolerance){
    
    const int n_pfp = tpc_obj.NumPFParticles();
    const double tpc_vtx_x = tpc_obj.pfpVtxX();
    const double tpc_vtx_y = tpc_obj.pfpVtxY();
    const double tpc_vtx_z = tpc_obj.pfpVtxZ();

    const int n_tracks = tpc_obj.NPfpTracks();
    if (n_tracks == 0 && pfp_pdg_type == 13 ) return true; 

    for (int j = 0; j < n_pfp; j++) {

        auto const part   = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();

        if (pfp_pdg == pfp_pdg_type) {

            const double pfp_vtx_x = part.pfpVtxX();
            const double pfp_vtx_y = part.pfpVtxY();
            const double pfp_vtx_z = part.pfpVtxZ();

            const double distance = sqrt(pow((tpc_vtx_x - pfp_vtx_x), 2) + pow((tpc_vtx_y - pfp_vtx_y), 2) + pow((tpc_vtx_z - pfp_vtx_z), 2) );

            if (distance <= tolerance) return true; 
        }

    }
    return false;

}
// -----------------------------------------------------------------------------
bool selection_cuts::HitThreshold(xsecAna::TPCObjectContainer tpc_obj, double threshold, bool useCollection){

    const int n_pfp = tpc_obj.NumPFParticles();

    for (int j = 0; j < n_pfp; j++) {

        auto const pfp_obj = tpc_obj.GetParticle(j);
        int  num_pfp_hits  = pfp_obj.NumPFPHits();
        const int  pfp_pdg = pfp_obj.PFParticlePdgCode();

        if (useCollection) num_pfp_hits = pfp_obj.NumPFPHitsW(); // Collection plane hits

        if (pfp_pdg == 11 && num_pfp_hits >= threshold) return true;
    }
    
    return false;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::OpenAngleCut(xsecAna::TPCObjectContainer tpc_obj, double tolerance_open_angle_min, double tolerance_open_angle_max){

    const int n_pfp = tpc_obj.NumPFParticles();

    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    auto const leading_shower       = tpc_obj.GetParticle(leading_index);
    const double leading_open_angle = leading_shower.pfpOpenAngle() * (180 / 3.1415);

    if (leading_open_angle <= tolerance_open_angle_max && leading_open_angle >= tolerance_open_angle_min)
        return true;
    else 
        return false;
}
// -----------------------------------------------------------------------------
bool selection_cuts::dEdxCut( xsecAna::TPCObjectContainer tpc_obj, const double tolerance_dedx_min, const double tolerance_dedx_max){

    const int n_pfp = tpc_obj.NumPFParticles();
    int leading_index = 0;
    int leading_hits  = 0;

    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    } 
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
    leading_dedx = leading_dedx * (196.979 /242.72);

    if (leading_dedx <= tolerance_dedx_max && leading_dedx >= tolerance_dedx_min) return true;
     
    return false;

}
// -----------------------------------------------------------------------------
bool selection_cuts::SecondaryShowersDistCut(xsecAna::TPCObjectContainer tpc_obj, const double dist_tolerance){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_showers = tpc_obj.NPfpShowers();
    
    // This cut does not target events with fewer than 2 showers
    if (n_pfp_showers <= 1) return true; 
    
    const double tpco_vtx_x = tpc_obj.pfpVtxX();
    const double tpco_vtx_y = tpc_obj.pfpVtxY();
    const double tpco_vtx_z = tpc_obj.pfpVtxZ();
    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    for (int j = 0; j < n_pfp; j++) {
        
        if (j == leading_index) continue; // We assume leading shower == electron shower
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const double pfp_vtx_x = part.pfpVtxX();
        const double pfp_vtx_y = part.pfpVtxY();
        const double pfp_vtx_z = part.pfpVtxZ();
        
        const double distance = sqrt(pow((pfp_vtx_x - tpco_vtx_x),2) + pow((pfp_vtx_y - tpco_vtx_y),2) + pow((pfp_vtx_z - tpco_vtx_z),2));
        
        if (pfp_pdg == 11) {
            if (distance > dist_tolerance) return false;
            // if (distance <= dist_tolerance) return true;
                
        }
    }
    
    return true;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::HitLengthRatioCut(const double pfp_hits_length_tolerance, xsecAna::TPCObjectContainer tpc_obj){
    
    const int n_pfp = tpc_obj.NumPFParticles();
    int leading_index = 0;
    int leading_hits  = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const part = tpc_obj.GetParticle(j);
        const int pfp_pdg = part.PFParticlePdgCode();
        const int n_pfp_hits = part.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
    }
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    const int pfp_pdg = leading_shower.PFParticlePdgCode();
    const double pfp_hits = leading_shower.NumPFPHits();
    const double pfp_length = leading_shower.pfpLength();
    const double pfp_hits_length_ratio = (pfp_hits / pfp_length);

    if (pfp_pdg == 11 && pfp_hits_length_ratio > pfp_hits_length_tolerance ) return true;

    return false;
    
}
// -----------------------------------------------------------------------------
bool selection_cuts::LongestTrackLeadingShowerCut(const double ratio_tolerance, xsecAna::TPCObjectContainer tpc_obj){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_tracks = tpc_obj.NPfpTracks();
    
    if (n_pfp_tracks == 0) return true;
    
    int leading_index = 0;
    int leading_hits  = 0;
    double longest_track = 0;
    
    for (int j = 0; j < n_pfp; j++) {
        
        auto const pfp = tpc_obj.GetParticle(j);
        const int pfp_pdg = pfp.PFParticlePdgCode();
        const int n_pfp_hits = pfp.NumPFPHits();
        
        if (pfp_pdg == 11 && n_pfp_hits > leading_hits) {
            leading_hits = n_pfp_hits;
            leading_index = j;
        }
        
        if(pfp_pdg == 13) {
            
            const double trk_length = pfp.pfpLength();
            
            if (trk_length > longest_track) longest_track = trk_length;
            
        }
    
    } //end loop pfparticles
    
    auto const leading_shower = tpc_obj.GetParticle(leading_index);
    const double leading_shower_length = leading_shower.pfpLength();
    const double longest_track_leading_shower_ratio = longest_track / leading_shower_length;

    //if the ratio is too large:
    if (longest_track_leading_shower_ratio > ratio_tolerance) return false;

    return true;

}
// -----------------------------------------------------------------------------
bool selection_cuts::IsContained(std::vector<double> track_start, std::vector<double> track_end, std::vector<double> fv_boundary_v) {
    
    if(in_fv(track_start.at(0), track_start.at(1), track_start.at(2), fv_boundary_v) == true
       && in_fv(track_end.at(0), track_end.at(1), track_end.at(2), fv_boundary_v) == true) {
        return true;
    }
    else 
        return false;
}
bool selection_cuts::ContainedTracksCut(std::vector<double> fv_boundary_v, xsecAna::TPCObjectContainer tpc_obj){

    const int n_pfp = tpc_obj.NumPFParticles();
    const int n_pfp_tracks = tpc_obj.NPfpTracks();
    
    // This is normally enabled, but due to test cut below it is off
    if (n_pfp_tracks == 0) return true;

    for (int j = 0; j < n_pfp; j++) {
        
        auto const pfp = tpc_obj.GetParticle(j);
        const int pfp_pdg = pfp.PFParticlePdgCode();
        
        if (pfp_pdg == 13) {
            
            const double pfp_vtx_x = pfp.pfpVtxX();
            const double pfp_vtx_y = pfp.pfpVtxY();
            const double pfp_vtx_z = pfp.pfpVtxZ();
            const double pfp_dir_x = pfp.pfpDirX();
            const double pfp_dir_y = pfp.pfpDirY();
            const double pfp_dir_z = pfp.pfpDirZ();
            const double trk_length = pfp.pfpLength();
            const double pfp_end_x = (pfp.pfpVtxX() + (trk_length * pfp_dir_x));
            const double pfp_end_y = (pfp.pfpVtxY() + (trk_length * pfp_dir_y));
            const double pfp_end_z = (pfp.pfpVtxZ() + (trk_length * pfp_dir_z));

            std::vector<double> pfp_start_vtx {pfp_vtx_x, pfp_vtx_y, pfp_vtx_z};
            std::vector<double> pfp_end_vtx {pfp_end_x, pfp_end_y, pfp_end_z};

            const bool is_contained = IsContained(pfp_start_vtx, pfp_end_vtx, fv_boundary_v);

            //if not contained
            if(is_contained == false) return false;
        
        } // end is track
    
    } // end loop pfprticles
    return true;	
}
// -----------------------------------------------------------------------------
void selection_cuts::TabulateOrigins(std::vector<double> &tabulated_origins) {
    
    int nue_cc        = 0;
    int nue_cc_qe     = 0;
    int nue_cc_res    = 0;
    int nue_cc_dis    = 0;
    int nue_cc_coh    = 0;
    int nue_cc_mec    = 0;

    int nue_bar_cc_qe     = 0;
    int nue_bar_cc_res    = 0;
    int nue_bar_cc_dis    = 0;
    int nue_bar_cc_coh    = 0;
    int nue_bar_cc_mec    = 0;

    int total_nue_cc_qe     = 0;
    int total_nue_cc_res    = 0;
    int total_nue_cc_dis    = 0;
    int total_nue_cc_coh    = 0;
    int total_nue_cc_mec    = 0;

    int nue_cc_mixed  = 0;
    int nue_cc_out_fv = 0;
    int cosmic        = 0;
    int nc            = 0;
    int numu_cc       = 0;
    int numu_cc_qe    = 0;
    int numu_cc_res   = 0;
    int numu_cc_dis   = 0;
    int numu_cc_coh   = 0;
    int numu_cc_mec   = 0;
    int numu_cc_mixed = 0;
    int nc_pi0        = 0;
    int unmatched     = 0;
    int other_mixed   = 0;
    int total         = 0;
    int signal_tpco_num = -1;
    int only_nue_cc = 0;
    int only_nue_bar_cc = 0;


    std::string tpco_id = tpc_classification.first;

    if(tpco_id == "nue_cc_qe")       {nue_cc_qe++;  }
    if(tpco_id == "nue_cc_res")      {nue_cc_res++; }
    if(tpco_id == "nue_cc_coh")      {nue_cc_coh++; }
    if(tpco_id == "nue_cc_dis")      {nue_cc_dis++; }
    if(tpco_id == "nue_cc_mec")      {nue_cc_mec++; }

    if(tpco_id == "nue_bar_cc_qe")   {nue_bar_cc_qe++;  }
    if(tpco_id == "nue_bar_cc_res")  {nue_bar_cc_res++; }
    if(tpco_id == "nue_bar_cc_coh")  {nue_bar_cc_coh++; }
    if(tpco_id == "nue_bar_cc_dis")  {nue_bar_cc_dis++; }
    if(tpco_id == "nue_bar_cc_mec")  {nue_bar_cc_mec++; }

    if(tpco_id == "nue_cc_out_fv")   {nue_cc_out_fv++; }
    if(tpco_id == "nue_cc_mixed")    {nue_cc_mixed++; }
    if(tpco_id == "nc")              {nc++; }
    if(tpco_id == "numu_cc_qe")      {numu_cc_qe++; }
    if(tpco_id == "numu_cc_res")     {numu_cc_res++; }
    if(tpco_id == "numu_cc_coh")     {numu_cc_coh++; }
    if(tpco_id == "numu_cc_dis")     {numu_cc_dis++; }
    if(tpco_id == "numu_cc_mec")     {numu_cc_mec++; }
    if(tpco_id == "numu_cc_mixed")   {numu_cc_mixed++; }
    if(tpco_id == "nc_pi0")          {nc_pi0++; }
    if(tpco_id == "cosmic")          {cosmic++; }
    if(tpco_id == "other_mixed")     {other_mixed++; }
    if(tpco_id == "unmatched")       {unmatched++; }

    only_nue_cc = nue_cc_qe + nue_cc_res + nue_cc_dis + nue_cc_coh + nue_cc_mec;
    only_nue_bar_cc = nue_bar_cc_qe + nue_bar_cc_res + nue_bar_cc_dis + nue_bar_cc_coh + nue_bar_cc_mec;

    total_nue_cc_qe  = nue_cc_qe  + nue_bar_cc_qe;
    total_nue_cc_res = nue_cc_res + nue_bar_cc_res;
    total_nue_cc_dis = nue_cc_dis + nue_bar_cc_dis;
    total_nue_cc_coh = nue_cc_coh + nue_bar_cc_coh;
    total_nue_cc_mec = nue_cc_mec + nue_bar_cc_mec;


    nue_cc = only_nue_cc + only_nue_bar_cc;
    numu_cc = numu_cc_qe + numu_cc_res + numu_cc_dis + numu_cc_coh + numu_cc_mec;
    total = nue_cc + nue_cc_mixed + nue_cc_out_fv + cosmic + nc + numu_cc + numu_cc_mixed + nc_pi0 + unmatched + other_mixed;

    tabulated_origins.at(0)  += nue_cc;//this is nue_cc + nue_bar_cc
    tabulated_origins.at(1)  += nue_cc_mixed;
    tabulated_origins.at(2)  += cosmic;
    tabulated_origins.at(3)  += nc;
    tabulated_origins.at(4)  += numu_cc;
    tabulated_origins.at(5)  += unmatched;
    tabulated_origins.at(6)  += other_mixed;
    tabulated_origins.at(7)  += total;
    tabulated_origins.at(8)  += signal_tpco_num;
    tabulated_origins.at(9)  += nue_cc_out_fv;
    tabulated_origins.at(10) += nc_pi0;
    tabulated_origins.at(11) += numu_cc_mixed;
    tabulated_origins.at(12) += total_nue_cc_qe;
    tabulated_origins.at(13) += total_nue_cc_res;
    tabulated_origins.at(14) += total_nue_cc_dis;
    tabulated_origins.at(15) += total_nue_cc_coh;
    tabulated_origins.at(16) += total_nue_cc_mec;
    tabulated_origins.at(17) += numu_cc_qe;
    tabulated_origins.at(18) += numu_cc_res;
    tabulated_origins.at(19) += numu_cc_dis;
    tabulated_origins.at(20) += numu_cc_coh;
    tabulated_origins.at(21) += numu_cc_mec;
    tabulated_origins.at(22) += only_nue_cc;
    tabulated_origins.at(23) += only_nue_bar_cc;
}
// -----------------------------------------------------------------------------