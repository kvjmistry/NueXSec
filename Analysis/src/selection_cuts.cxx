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
    std::pair<std::string, int> tpc_classification = TPCO_Classifier(tpc_obj, true_in_tpc, has_pi0);

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
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------