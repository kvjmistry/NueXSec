#include "../include/SliceContainer.h"

// -----------------------------------------------------------------------------
void SliceContainer::Initialise(TTree *tree, TTree *mc_truth_tree, int type, TFile *f_flux_weights, const char * _run_period){

    std::cout << "Initalising Slice Container" << std::endl;

    tree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

    if (type == _util.k_mc){
        mc_truth_tree->SetBranchAddress("mc_nue_cc_counter",      &mc_nue_cc_counter);
        mc_truth_tree->SetBranchAddress("mc_nue_nc_counter",      &mc_nue_nc_counter);
        mc_truth_tree->SetBranchAddress("mc_numu_cc_counter",     &mc_numu_cc_counter);
        mc_truth_tree->SetBranchAddress("mc_numu_nc_counter",     &mc_numu_nc_counter);
        mc_truth_tree->SetBranchAddress("mc_nue_cc_counter_bar",  &mc_nue_cc_counter_bar);
        mc_truth_tree->SetBranchAddress("mc_numu_cc_counter_bar", &mc_numu_cc_counter_bar);
        mc_truth_tree->SetBranchAddress("mc_nue_nc_counter_bar",  &mc_nue_nc_counter_bar);
        mc_truth_tree->SetBranchAddress("mc_numu_nc_counter_bar", &mc_numu_nc_counter_bar);
        mc_truth_tree->SetBranchAddress("fMCNuEnergy", &mc_nu_energy);
        mc_truth_tree->SetBranchAddress("fMCNuMomentum", &mc_nu_momentum);
        mc_truth_tree->SetBranchAddress("fMCNuID", &mc_nu_id);
        mc_truth_tree->SetBranchAddress("fMCNuVtxX", &mc_nu_vtx_x);
        mc_truth_tree->SetBranchAddress("fMCNuVtxY", &mc_nu_vtx_y);
        mc_truth_tree->SetBranchAddress("fMCNuVtxZ", &mc_nu_vtx_z);
        mc_truth_tree->SetBranchAddress("fMCNuDirX", &mc_nu_dir_x);
        mc_truth_tree->SetBranchAddress("fMCNuDirY", &mc_nu_dir_y);
        mc_truth_tree->SetBranchAddress("fMCNuDirZ", &mc_nu_dir_z);
        mc_truth_tree->SetBranchAddress("fMCNumParticles", &mc_nu_num_particles);
        mc_truth_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
        mc_truth_tree->SetBranchAddress("fMCEleDirX", &mc_ele_dir_x);
        mc_truth_tree->SetBranchAddress("fMCEleDirY", &mc_ele_dir_y);
        mc_truth_tree->SetBranchAddress("fMCEleDirZ", &mc_ele_dir_z);
        mc_truth_tree->SetBranchAddress("fMCEleEnergy", &mc_ele_energy);
        mc_truth_tree->SetBranchAddress("fMCEleMomentum", &mc_ele_momentum);
        mc_truth_tree->SetBranchAddress("has_pi0", &has_pi0);
        mc_truth_tree->SetBranchAddress("fMCNuTime", &mc_nu_time);
    }

}
// -----------------------------------------------------------------------------
void SliceContainer::SetTPCObj(xsecAna::TPCObjectContainer tpc_obj, int type){

    // Number of pfp
    n_pfp          = tpc_obj.NumPFParticles();
    n_pfp_tracks   = tpc_obj.NPfpTracks();
    n_pfp_showers  = tpc_obj.NPfpShowers();
    tpc_obj_origin = tpc_obj.Origin();

    tpco_vtx_x = tpc_obj.pfpVtxX();
    tpco_vtx_y = tpc_obj.pfpVtxY();
    tpco_vtx_z = tpc_obj.pfpVtxZ();

    tpco_hits = tpc_obj.NumPFPHits();


    run    = tpc_obj.RunNumber();
    subrun = tpc_obj.SubRunNumber();
    event  = tpc_obj.EventNumber();

    // Get the leading shower index
    GetLeadingShowerIndex();

    // Check if vertex is in the FV
    if (type == _util.k_mc) infv = InFV();

    // Get the Particle vector
    for (int j = 0; j < n_pfp ; j++){

        auto const pfp_obj = tpc_obj.GetParticle(j);

        // PFP vars
        mc_origin    .push_back(pfp_obj.Origin());
        pfp_pdg      .push_back(pfp_obj.PFParticlePdgCode());
        num_pfp_hits .push_back(pfp_obj.NumPFPHits());
        mc_parent_pdg.push_back(pfp_obj.MCParentPdg());
        n_pfp_hits_w .push_back(pfp_obj.NumPFPHitsW());
        
        pfp_length     .push_back(pfp_obj.pfpLength());
        pfp_open_angle .push_back(pfp_obj.pfpOpenAngle());
        leading_dedx   .push_back(pfp_obj.PfpdEdx().at(2)); 

        pfp_vtx_x .push_back(pfp_obj.pfpVtxX());
        pfp_vtx_y .push_back(pfp_obj.pfpVtxY());
        pfp_vtx_z .push_back(pfp_obj.pfpVtxZ());

        pfp_dir_x .push_back(pfp_obj.pfpDirX());
        pfp_dir_y .push_back(pfp_obj.pfpDirY());
        pfp_dir_z .push_back(pfp_obj.pfpDirZ());
        
        mc_Theta  .push_back(pfp_obj.mcTheta());
        mc_Phi    .push_back(pfp_obj.mcPhi());
        mc_Energy .push_back(pfp_obj.mcEnergy());
        mc_pdg    .push_back(pfp_obj.MCPdgCode());

        pfp_end_x .push_back((pfp_obj.pfpVtxX() + (pfp_obj.pfpLength() * pfp_obj.pfpDirX())));
        pfp_end_y .push_back((pfp_obj.pfpVtxY() + (pfp_obj.pfpLength() * pfp_obj.pfpDirY())));
        pfp_end_z .push_back((pfp_obj.pfpVtxZ() + (pfp_obj.pfpLength() * pfp_obj.pfpDirZ())));

        CCNC.push_back(pfp_obj.CCNC());

        dist_pfp_nu_vtx .push_back(sqrt(pow((pfp_obj.pfpVtxX() - tpco_vtx_x),2) + pow((pfp_obj.pfpVtxY() - tpco_vtx_y),2) + pow((pfp_obj.pfpVtxZ() - tpco_vtx_z),2)));

    }

}
// -----------------------------------------------------------------------------
void SliceContainer::Reset(){

    n_pfp = -1;
    longest_track_index = -1;
    leading_shower_index = -1;
    subleading_shower_index = -1;
    n_pfp_tracks = -1;
    n_pfp_showers = -1;

    tpco_vtx_x = -1;
    tpco_vtx_y = -1;
    tpco_vtx_z = -1;
    tpco_hits = -1;

    tpc_obj_origin = "empty";

    infv = true;


    run    = -1;
    subrun = -1;
    event  = -1;

    mc_origin.clear();
    pfp_pdg.clear();
    num_pfp_hits.clear();
    mc_parent_pdg.clear();
    n_pfp_hits_w.clear();
  
    pfp_length.clear();
    pfp_open_angle.clear();
    leading_dedx.clear();
  
    pfp_vtx_x.clear();
    pfp_vtx_y.clear();
    pfp_vtx_z.clear();
    pfp_dir_x.clear();
    pfp_dir_y.clear();
    pfp_dir_z.clear();
  
    mc_pdg.clear();
    mc_Theta.clear();
    mc_Phi.clear();
    mc_Energy.clear();
    
    pfp_end_x.clear();
    pfp_end_y.clear();
    pfp_end_z.clear();

    dist_pfp_nu_vtx.clear();

    CCNC.clear();


}
// -----------------------------------------------------------------------------
void SliceContainer::GetLeadingShowerIndex(){

    int leading_hits  = 0;
    
    // Loop over the particles and get the leading shower index
    for ( int j = 0; j < n_pfp; j++) {

        if (pfp_pdg.at(j) == 11 && num_pfp_hits.at(j) > leading_hits) {
            leading_hits = num_pfp_hits.at(j);
            leading_shower_index = j;
        }
   
    }


    // Do it again and get the subleading shower index
    for ( int j = 0; j < n_pfp; j++) {

        // Skip the leading shower in this case
        if (j == leading_shower_index ) continue;

        if (pfp_pdg.at(j) == 11 && num_pfp_hits.at(j) > leading_hits) {
            leading_hits = num_pfp_hits.at(j);
            subleading_shower_index = j;
        }
   
    }


}
// -----------------------------------------------------------------------------
std::pair<std::string, int> SliceContainer::SliceClassifier(int type){
    
    // MC Specific classsifications
    if (type == _util.k_mc){

        // Some TPC Objects have 0 hits!
        if(tpco_hits == 0) return std::make_pair("unmatched", _util.k_unmatched);

        

        bool part_neutrino, part_cosmic, part_unmatched, part_numu, part_nue, part_cc;


        // loop over the origin vector and see whats in it
        for (unsigned int k =0; k < mc_origin.size();k++){
            if (mc_origin.at(k) == "kBeamNeutrino") part_neutrino = true;
            if (mc_origin.at(k) == "kCosmicRay")    part_cosmic = true;
            if (mc_origin.at(k) == "kUnknown")      part_unmatched = true;

            if (mc_parent_pdg.at(k) == 12 || mc_parent_pdg.at(k) == -12) part_nue = true;
            if (mc_parent_pdg.at(k) == 14 || mc_parent_pdg.at(k) == -14) part_numu = true;

            if (CCNC.at(k) == _util.k_CC) part_cc = true;
            else part_cc = false;
        }

        if (part_unmatched) return std::make_pair("unmatched", _util.k_unmatched);

        // CC Interactions
        if (part_cc){
            
            // Mixed events have cosmic contamination
            if( part_cosmic) {
                if ( part_nue && part_neutrino)        return std::make_pair("nue_cc_mixed", _util.k_nue_cc_mixed );
                else if ( part_numu && part_neutrino)  return std::make_pair("numu_cc",_util.k_numu_cc );
                else return std::make_pair("cosmic", _util.k_cosmic); 
            }

            // Nue CC
            if (part_nue) return std::make_pair("nue_cc", _util.k_cosmic );
            else return std::make_pair("numu_cc",_util.k_numu_cc  );
        
        }
        // NC events
        else {

            // Mixed NC event
            if( part_cosmic ) {
                return std::make_pair("nc_mixed", _util.k_nc_mixed );
            }
            else {
                if (has_pi0) return std::make_pair("nc_pi0", _util.k_nc_pi0);
                else return std::make_pair("nc",_util.k_nc);
            }
        }
        
    }
    // Data
    else if (type == _util.k_data){
        return std::make_pair("data",_util.k_leg_data);
    }
    // EXT
    else if (type == _util.k_ext){
        return std::make_pair("ext",_util.k_leg_ext);
        
    }
    // Dirt
    else if (type == _util.k_dirt){
        return std::make_pair("dirt",_util.k_leg_dirt);
    }
    // What is this type?
    else return std::make_pair("unmatched",_util.k_unmatched);
    
}
// -----------------------------------------------------------------------------
std::pair<std::string, int> SliceContainer::ParticleClassifier(int type){
    
    // MC Specific classsifications
    // if (type == _util.k_mc){

    //     // Electron
    //     if (shr_bkt_pdg == 11 || shr_bkt_pdg == -11){
    //         return std::make_pair("e",_util.k_electron);
    //     }
    //     // Muon
    //     else if (shr_bkt_pdg == 13 || shr_bkt_pdg == -13){
    //         return std::make_pair("muon",_util.k_muon);
    //     }
    //     // Pion
    //     else if (shr_bkt_pdg == 211 || shr_bkt_pdg == -211){
    //         return std::make_pair("e",_util.k_pion);
    //     }
    //     // Photon 
    //     else if (shr_bkt_pdg == 22 ){
    //         return std::make_pair("photon",_util.k_photon);
    //     }
    //     // Proton
    //     else if (shr_bkt_pdg == 2212){
    //         return std::make_pair("p",_util.k_proton);
    //     }
    //     // Neutron
    //     else if (shr_bkt_pdg == 2112){
    //         return std::make_pair("n",_util.k_neutron);
    //     }
    //     // Kaon
    //     else if (shr_bkt_pdg == 321 || shr_bkt_pdg == -321 ){
    //         return std::make_pair("K",_util.k_kaon);
    //     }
    //     // Other stuff is assumed cosmic
    //     else {
    //         return std::make_pair("cosmic",_util.k_part_cosmic);
    //     }


    // }
    // // Data
    // else if (type == _util.k_data){
    //     return std::make_pair("data",_util.k_part_data);
    // }
    // // EXT
    // else if (type == _util.k_ext){
    //     return std::make_pair("ext",_util.k_part_ext);
        
    // }
    // // Dirt
    // else if (type == _util.k_dirt){
    //     return std::make_pair("dirt",_util.k_part_dirt);
    // }
    // // What is this type?
    // else return std::make_pair("unmatched",_util.k_part_unmatched);
    
}
// -----------------------------------------------------------------------------
std::string SliceContainer::SliceInteractionType(int type){

    // // Only do this for mc, otherwise return data type
    // if (type == _util.k_mc || type == _util.k_dirt){
    //     std::string nu;
    //     std::string CCNC;

    //     // Get the nu flavour
    //     if (nu_pdg == 14 ){
    //         nu = "numu_";
    //     }
    //     else if (nu_pdg == -14){
    //         nu = "numu_bar_";
    //     }
    //     else if (nu_pdg == 12){
    //         nu = "nue_";
    //     }
    //     else if (nu_pdg == -12){
    //         nu = "nue_bar_";
    //     }
    //     else {
    //         nu = "unknown_";
    //     }

    //     // The interaction type
    //     if (ccnc == _util.k_CC){
    //         CCNC = "cc_";
    //     }
    //     else CCNC = "nc_";


    //     if (interaction == _util.k_qe) {
    //         return nu + CCNC + "qe";

    //     }
    //     else if (interaction == _util.k_res ) {
    //         return nu + CCNC + "res";

    //     }
    //     else if (interaction == _util.k_dis ) {
    //         return nu + CCNC + "dis";

    //     }
    //     else if (interaction == _util.k_coh) {
    //         return nu + CCNC + "coh";

    //     }
    //     else if (interaction == _util.k_mec) {
    //         return nu + CCNC + "mec";

    //     }
    // }
    // else return "data";


}
// -----------------------------------------------------------------------------
bool SliceContainer::InFV(){
    
    if (mc_nu_vtx_x <= 0 || mc_nu_vtx_x >= 256.35)     return false;
    if (mc_nu_vtx_y <= -116.5 || mc_nu_vtx_y >= 116.5) return false;
    if (mc_nu_vtx_z <= 0 || mc_nu_vtx_z >= 1036.8)     return false;
    
    return true;

}
// -----------------------------------------------------------------------------