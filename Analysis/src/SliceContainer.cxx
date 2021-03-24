#include "../include/SliceContainer.h"

// -----------------------------------------------------------------------------
template<typename T> void SliceContainer::Initialise(T *tree, int type, Utility util){

    std::cout << "Initalising Slice Container" << std::endl;
    _util = util;

    tree->SetBranchAddress("flash_time_v", &flash_time_v);
    tree->SetBranchAddress("has_nue", &has_nue);
    tree->SetBranchAddress("has_valid_shr", &has_valid_shr);
    tree->SetBranchAddress("trk_vtx_dist_v", &trk_vtx_dist_v);
    tree->SetBranchAddress("shr_vtx_dist_v", &shr_vtx_dist_v);
    tree->SetBranchAddress("sec_shr_vtx_dist_v", &sec_shr_vtx_dist_v);
    tree->SetBranchAddress("trk_start_x_v", &trk_start_x_v);
    tree->SetBranchAddress("trk_start_y_v", &trk_start_y_v);
    tree->SetBranchAddress("trk_start_z_v", &trk_start_z_v);
    tree->SetBranchAddress("trk_end_x_v", &trk_end_x_v);
    tree->SetBranchAddress("trk_end_y_v", &trk_end_y_v);
    tree->SetBranchAddress("trk_end_z_v", &trk_end_z_v);
    tree->SetBranchAddress("flash_vtx_dist", &flash_vtx_dist);
    tree->SetBranchAddress("flash_vtx_dist_updated", &flash_vtx_dist_updated);
    tree->SetBranchAddress("flash_vtx_y", &flash_vtx_y);
    tree->SetBranchAddress("flash_vtx_z", &flash_vtx_z);
    tree->SetBranchAddress("flash_vtx_y_updated", &flash_vtx_y_updated);
    tree->SetBranchAddress("flash_vtx_z_updated", &flash_vtx_z_updated);
    tree->SetBranchAddress("tpc_obj_id", &tpc_obj_index);
    tree->SetBranchAddress("number_tpcobj", &number_tpcobj);
    tree->SetBranchAddress("shr_theta", &shr_theta);
    tree->SetBranchAddress("shr_phi", &shr_phi);
    tree->SetBranchAddress("shr_len", &shr_len);
    tree->SetBranchAddress("shr_px", &shr_px);
    tree->SetBranchAddress("shr_py", &shr_py);
    tree->SetBranchAddress("shr_pz", &shr_pz);
    tree->SetBranchAddress("shr_openangle", &shr_openangle);
    tree->SetBranchAddress("shr_start_x", &shr_start_x);
    tree->SetBranchAddress("shr_start_y", &shr_start_y);
    tree->SetBranchAddress("shr_start_z", &shr_start_z);
    tree->SetBranchAddress("shr_dedx_Y", &shr_dedx_Y);
    tree->SetBranchAddress("shr_dedx_V", &shr_dedx_V);
    tree->SetBranchAddress("shr_dedx_U", &shr_dedx_U);
    tree->SetBranchAddress("shr_distance", &shr_distance);
    tree->SetBranchAddress("trk_len", &trk_len);
    tree->SetBranchAddress("shr_hits_tot", &shr_hits_tot);
    tree->SetBranchAddress("shr_hits_y_tot", &shr_hits_y_tot);
    tree->SetBranchAddress("shr_hits_u_tot", &shr_hits_u_tot);
    tree->SetBranchAddress("shr_hits_v_tot", &shr_hits_v_tot);
    tree->SetBranchAddress("nu_pdg", &nu_pdg);
    tree->SetBranchAddress("ccnc",   &ccnc);
    tree->SetBranchAddress("interaction", &interaction);
    tree->SetBranchAddress("nu_e",   &nu_e);
    tree->SetBranchAddress("nu_pt",  &nu_pt);
    tree->SetBranchAddress("true_nu_vtx_t", &true_nu_vtx_t);
    tree->SetBranchAddress("true_nu_vtx_x", &true_nu_vtx_x);
    tree->SetBranchAddress("true_nu_vtx_y", &true_nu_vtx_y);
    tree->SetBranchAddress("true_nu_vtx_z", &true_nu_vtx_z);    
    tree->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x);
    tree->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y);
    tree->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z);
    tree->SetBranchAddress("nelec", &nelec);
    tree->SetBranchAddress("elec_e", &elec_e);
    tree->SetBranchAddress("npi0", &npi0);
    tree->SetBranchAddress("n_pfps", &n_pfps);
    tree->SetBranchAddress("n_tracks", &n_tracks);
    tree->SetBranchAddress("n_showers", &n_showers);
    tree->SetBranchAddress("flash_pe", &flash_pe);
    tree->SetBranchAddress("flash_time", &flash_time);
    tree->SetBranchAddress("failed_intime_cut", &failed_intime_cut);
    tree->SetBranchAddress("failed_flash_cut", &failed_flash_cut);
    tree->SetBranchAddress("nu_purity_from_pfp", &nu_purity_from_pfp);

}
// Explicitly initialise the templates for this function which include TTree and TChain
template void SliceContainer::Initialise<TTree>(TTree *tree, int type, Utility util);
template void SliceContainer::Initialise<TChain>(TChain *tree, int type, Utility util );
// -----------------------------------------------------------------------------
void SliceContainer::SliceClassifier(int type){

    
    // MC Specific classsifications
    if (type == _util.k_mc){
        
        bool is_in_fv = _util.in_fv(true_nu_vtx_x, true_nu_vtx_y, true_nu_vtx_z);

        // Out of Fiducial Volume Event
        if (!is_in_fv) {
            
            // Cosmic events including mixed
            if (nu_purity_from_pfp <= 0.5){
                
                classification = std::make_pair("cosmic",_util.k_cosmic);
                return;
            }
            // Out of Fv events
            else {
                classification = std::make_pair("nu_out_fv",_util.k_nu_out_fv);
                return;
            }
        }
        // In FV event
        else {

            // Charged Current 
            if (ccnc == _util.k_CC){

                // NuMu CC
                if (nu_pdg == 14 || nu_pdg == -14){

                    // Purity is low so return cosmic
                    if (nu_purity_from_pfp <= 0.5) {
                        classification = std::make_pair("cosmic",_util.k_cosmic);
                        return;
                    }
                    
                    if (npi0 > 0) {
                        classification = std::make_pair("numu_cc", _util.k_numu_cc); // has a pi0
                        return;
                    }
                    else {
                        classification = std::make_pair("numu_cc",_util.k_numu_cc);
                        return;
                    }

                }
                // Nue CC
                else if (nu_pdg == 12){
                    
                    // purity > 0.5% so signal
                    if (nu_purity_from_pfp >= 0.5){
                        classification = std::make_pair("nue_cc",       _util.k_nue_cc);    
                        return;
                    }
                    // These events were not picked up by pandora at all
                    else {
                        classification = std::make_pair("cosmic",_util.k_cosmic); 
                        return;
                    }

                }
                else if (nu_pdg == -12){
                    
                    // purity > 0.5% so signal
                    if (nu_purity_from_pfp >= 0.5) {
                        classification = std::make_pair("nuebar_cc",       _util.k_nuebar_cc); 
                        return;

                    }
                    // These events were not picked up by pandora at all
                    else {
                        classification = std::make_pair("cosmic",_util.k_cosmic); 
                        return;
                    }

                }
                // Unknown Neutrino Type
                else {
                    std::cout << "Unknown Neutrino Type..., This will also mess up the efficecy if this occurs!" << std::endl;
                    classification = std::make_pair("unmatched",_util.k_unmatched);
                    return;
                }

            }
            // Neutral Current
            else {

                // Purity is low so return cosmic
                if (nu_purity_from_pfp <= 0.5){
                    classification = std::make_pair("cosmic",_util.k_cosmic);
                    return;
                }

                if (npi0 > 0) {
                    classification = std::make_pair("nc_pi0",_util.k_nc_pi0);
                    return;
                }
                else {
                    classification = std::make_pair("nc",_util.k_nc);
                    return;
                }
            }
        
        } // End if in FV

    }
    // Data
    else if (type == _util.k_data){
        classification = std::make_pair("data",_util.k_leg_data);
        return;
    }
    // EXT
    else if (type == _util.k_ext){
        classification = std::make_pair("ext",_util.k_leg_ext);
        return;
        
    }
    // Dirt
    else if (type == _util.k_dirt){
        classification = std::make_pair("dirt",_util.k_leg_dirt);
        return;
    }
    // What is this type?
    else {
        std::cout << "Got a case we are calling unmatched, this is going to mess up the efficiency in the current way!" << std::endl;
        classification = std::make_pair("unmatched",_util.k_unmatched);
        return;
    }
    
}
// -----------------------------------------------------------------------------
std::string SliceContainer::SliceCategory(){

    if (category == _util.k_pandora_nu_e_other) {
        return "nue_other";

    }
    else if (category == _util.k_pandora_nu_e_cc0pi0p ) {
        return "nu_e_cc0pi0p";

    }
    else if (category == _util.k_pandora_nu_e_cc0pinp ) {
        return "nu_e_cc0pinp";

    }
    else if (category == _util.k_pandora_nu_mu_other) {
        return "nu_mu_other";

    }
    else if (category == _util.k_pandora_nu_mu_pi0 ) {
        return "nu_mu_pi0";

    }
    else if (category == _util.k_pandora_nc) {
        return "nc";

    }
    else if (category == _util.k_pandora_nc_pi0 ) {
        return "nc_pi0";

    }
    else if (category == _util.k_pandora_cosmic) {
        return "cosmic";

    }
    else if (category == _util.k_pandora_outfv) {
        return "outfv";

    }
    else if (category == _util.k_pandora_other) {
        return "other";

    }
    else if (category == _util.k_pandora_data) {
        return "data";
    }
    else {
        std::cout << "Unknown Category type"<< std::endl;
        return "unknown";
    }
}
// -----------------------------------------------------------------------------
void SliceContainer::ParticleClassifier(int type){
    
    // MC Specific classsifications
    if (type == _util.k_mc){

        // Electron
        if (shr_bkt_pdg == 11 || shr_bkt_pdg == -11){
            particle_type = std::make_pair("e",_util.k_electron);
            return;
        }
        // Muon
        else if (shr_bkt_pdg == 13 || shr_bkt_pdg == -13){
            particle_type = std::make_pair("muon",_util.k_muon);
            return;
        }
        // Pion
        else if (shr_bkt_pdg == 211 || shr_bkt_pdg == -211){
            particle_type = std::make_pair("e",_util.k_pion);
            return;
        }
        // Photon 
        else if (shr_bkt_pdg == 22 ){
            particle_type = std::make_pair("photon",_util.k_photon);
            return;
        }
        // Proton
        else if (shr_bkt_pdg == 2212){
            particle_type = std::make_pair("p",_util.k_proton);
            return;
        }
        // Neutron
        else if (shr_bkt_pdg == 2112){
            particle_type = std::make_pair("n",_util.k_neutron);
            return;
        }
        // Kaon
        else if (shr_bkt_pdg == 321 || shr_bkt_pdg == -321 ){
            particle_type = std::make_pair("K",_util.k_kaon);
            return;
        }
        // Other stuff is assumed cosmic
        else {
            particle_type = std::make_pair("cosmic",_util.k_part_cosmic);
            return;
        }


    }
    // Data
    else if (type == _util.k_data){
        particle_type = std::make_pair("data",_util.k_part_data);
        return;
    }
    // EXT
    else if (type == _util.k_ext){
        particle_type = std::make_pair("ext",_util.k_part_ext);
        return;
        
    }
    // Dirt
    else if (type == _util.k_dirt){
        particle_type = std::make_pair("dirt",_util.k_part_dirt);
        return;
    }
    // What is this type?
    else {
        particle_type = std::make_pair("unmatched",_util.k_part_unmatched);
        return;
    }
    
}
// -----------------------------------------------------------------------------
void SliceContainer::SliceInteractionType(int type){

    // Only do this for mc, otherwise return data type
    if (type == _util.k_mc || type == _util.k_dirt){
        std::string nu = "temp";
        std::string CCNC = "temp";

        // Get the nu flavour
        if (nu_pdg == 14 ){
            nu = "numu_";
        }
        else if (nu_pdg == -14){
            nu = "numu_bar_";
        }
        else if (nu_pdg == 12){
            nu = "nue_";
        }
        else if (nu_pdg == -12){
            nu = "nue_bar_";
        }
        else {
            nu = "unknown_";
        }

        // The interaction type
        if (ccnc == _util.k_CC){
            CCNC = "cc_";
        }
        else CCNC = "nc_";


        if (interaction == _util.k_qe) {
            genie_interaction = nu + CCNC + "qe";
            return;

        }
        else if (interaction == _util.k_res ) {
            genie_interaction = nu + CCNC + "res";
            return;

        }
        else if (interaction == _util.k_dis ) {
            genie_interaction = nu + CCNC + "dis";
            return;

        }
        else if (interaction == _util.k_coh) {
            genie_interaction = nu + CCNC + "coh";
            return;

        }
        else if (interaction == _util.k_mec) {
            genie_interaction = nu + CCNC + "mec";
            return;

        }
        else {
            genie_interaction = nu + CCNC + "unknown";
            return;
        }
    }
    else{
        genie_interaction = "data";
        return;
    }



}
// -----------------------------------------------------------------------------
void SliceContainer::Pi0Classifier(int type){

    // Only do this for mc, otherwise return data type
    if (type == _util.k_mc){
        std::string nu = "temp";
        std::string CCNC = "temp";

        // Get the nu flavour
        if (nu_pdg == 14 ||  nu_pdg == -14){
            nu = "numu_";
        }
        else if (nu_pdg == 12){
            nu = "nue_";
        }
        else if (nu_pdg == -12){
            nu = "nue_bar_";
        }
        else {
            nu = "unknown_";
        }

        // The interaction type
        if (ccnc == _util.k_CC){
            CCNC = "cc";
        }
        else CCNC = "nc";
        

        if (npi0 > 0){
            pi0_classification = nu + CCNC + "_pi0";
            return;

        }
        else if (npi0 == 0 ) {
            pi0_classification = nu + CCNC + "";
            return;

        }
        else {
            pi0_classification = nu + CCNC + "unknown";
            return;
        }
    }
    else {
        pi0_classification = "data";
        return;
    }



}
// -----------------------------------------------------------------------------
void SliceContainer::SetPPFXCVWeight(){
    
    float weight = 1.0;

    double xbin{1.0},ybin{1.0};

    if (nu_pdg == 14) {
        xbin =  hist_ratio.at(k_numu)->GetXaxis()->FindBin(nu_e);
        ybin =  hist_ratio.at(k_numu)->GetYaxis()->FindBin(nu_angle);
        weight =  hist_ratio.at(k_numu)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -14) {
        xbin =  hist_ratio.at(k_numubar)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_numubar)->GetYaxis()->FindBin(nu_angle);
        weight = hist_ratio.at(k_numubar)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == 12) {
        xbin = hist_ratio.at(k_nue)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_nue)->GetYaxis()->FindBin(nu_angle);
        weight = hist_ratio.at(k_nue)->GetBinContent(xbin, ybin);
    }
    if (nu_pdg == -12) {
        xbin = hist_ratio.at(k_nuebar)->GetXaxis()->FindBin(nu_e);
        ybin = hist_ratio.at(k_nuebar)->GetYaxis()->FindBin(nu_angle);
        weight = hist_ratio.at(k_nuebar)->GetBinContent(xbin, ybin);
    }

    // Add some catches to remove unphysical weights
    if (std::isinf(weight))      weight = 1.0; 
    if (std::isnan(weight) == 1) weight = 1.0;
    if (weight > 100)            weight = 1.0;

    ppfx_cv =  weight;

}
// -----------------------------------------------------------------------------
double SliceContainer::GetdEdxMax(){

    double dedx_max = -1;

    // We want to also use the dedx when it is defined properly. Sometimes, the plane can have hits but an undefined dedx
    // use the dedx where we get the max number of hits and the dedx > 0
    // int temp_shr_hits_u_tot = shr_hits_u_tot;
    // int temp_shr_hits_v_tot = shr_hits_v_tot;
    // int temp_shr_hits_y_tot = shr_hits_y_tot;

    int temp_shr_hits_u_tot = shr_tkfit_nhits_U; // These variables give a bigger difference in run 1 and run 3
    int temp_shr_hits_v_tot = shr_tkfit_nhits_V;
    int temp_shr_hits_y_tot = shr_tkfit_nhits_Y;

    // If the dedx is undefined, set the hits to zero
    if (shr_tkfit_dedx_U <= 0) temp_shr_hits_u_tot = 0;
    if (shr_tkfit_dedx_V <= 0) temp_shr_hits_v_tot = 0;
    if (shr_tkfit_dedx_Y <= 0) temp_shr_hits_y_tot = 0;


    // Collection plane is the largest
    if (temp_shr_hits_y_tot > temp_shr_hits_u_tot && temp_shr_hits_y_tot > temp_shr_hits_v_tot ){
        dedx_max = shr_tkfit_dedx_Y;
    }
    // V Plane is the largest
    else if (temp_shr_hits_v_tot > temp_shr_hits_u_tot && temp_shr_hits_v_tot > temp_shr_hits_y_tot) {
        dedx_max = shr_tkfit_dedx_V;
        
    }
    // U Plane is the largest
    else if (temp_shr_hits_u_tot > temp_shr_hits_v_tot && temp_shr_hits_u_tot > temp_shr_hits_y_tot){
        dedx_max = shr_tkfit_dedx_U;
        // std::cout << shr_tkfit_dedx_U << " " << shr_tkfit_dedx_V << " " << shr_tkfit_dedx_Y<< "  " << shr_theta*(180/3.14)<< "  " << npi0 <<  std::endl;
        // std::cout << shr_hits_u_tot << " " << shr_hits_v_tot << " " << shr_hits_y_tot<<"\n" <<std::endl;
    }
    // Ok one plane was equal, so need to prioritise planes in preference of y, v, u
    else {

        // If y == any other plane, then y wins
        if (temp_shr_hits_y_tot == temp_shr_hits_u_tot || temp_shr_hits_y_tot == temp_shr_hits_v_tot ){
            dedx_max = shr_tkfit_dedx_Y;
           
        }
        // U == V, ALL Y cases have been used up, so default to v
        else if (temp_shr_hits_u_tot == temp_shr_hits_v_tot ){
            dedx_max = shr_tkfit_dedx_V;
            
        }
        else {
            dedx_max = shr_tkfit_dedx_U;
            
            
        }
    }

    if (dedx_max == -1) {
        std::cout << shr_tkfit_dedx_U << " " << shr_tkfit_dedx_V << " " << shr_tkfit_dedx_Y<< std::endl;
        std::cout << "edge case of dedx comparisons, your logic is flawed!" << std::endl;
    }

    return dedx_max;
    
}
// -----------------------------------------------------------------------------
void SliceContainer::SetCVWeight(double weight){
    cv_weight = weight;
}
// -----------------------------------------------------------------------------
void SliceContainer::SetSignal(){

    if (classification.second == _util.k_nue_cc           || classification.second == _util.k_nuebar_cc ||
        classification.second == _util.k_unmatched_nue    || classification.second == _util.k_cosmic_nue ||
        classification.second == _util.k_unmatched_nuebar || classification.second == _util.k_cosmic_nuebar){
            is_signal = true;
    }
    else 
        is_signal = false;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetTrueElectronThetaPhi(){

    TVector3 vec(elec_px, elec_py, elec_pz);

    elec_theta = vec.Theta() * 180.0/3.14159;
    elec_phi   = vec.Phi()   * 180.0/3.14159;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetNuMIAngularVariables(){


    // --- Effective Angle -- //
    // Calculate the angle between the shower direction and the vector from the target to the nu vtx
    TVector3 shr_dir(shr_px, shr_py, shr_pz); // Shower direction
    shr_dir.Unit();
    
    TVector3 v_targ_uboone(-31387.58422, -3316.402543, -60100.2414);
    TVector3 v_nu_vtx(reco_nu_vtx_sce_x, reco_nu_vtx_sce_y, reco_nu_vtx_sce_z);
    TVector3 v_targ_to_vtx = (-1*v_targ_uboone + v_nu_vtx).Unit(); // -1 because the vector points from uboone to tgt, we need the other way around

    // Set the values
    effective_angle = shr_dir.Angle(v_targ_to_vtx) * 180 / 3.14159;
    cos_effective_angle = std::cos(shr_dir.Angle(v_targ_to_vtx));
    
    // Calculate the effective angle using true variables
    TVector3 nu_dir(true_nu_px, true_nu_py, true_nu_pz); 
    nu_dir.Unit();

    TVector3 v_nu_vtx_true(true_nu_vtx_sce_x, true_nu_vtx_sce_y, true_nu_vtx_sce_z);
    TVector3 v_targ_to_vtx_true = (-1*v_targ_uboone + v_nu_vtx_true).Unit(); // -1 because the vector points from uboone to tgt, we need the other way around

    TVector3 elec_dir(elec_px, elec_py, elec_pz); 
    elec_dir.Unit();
    // true_effective_angle = elec_dir.Angle(v_targ_to_vtx_true) * 180 / 3.14159;
    true_effective_angle = elec_dir.Angle(nu_dir) * 180 / 3.14159; // Use dot product of elec dir to nu dir for this as the true value
    
    // --

    // Momentum of neutrino
    nu_p = std::sqrt(true_nu_px*true_nu_px + true_nu_py*true_nu_py + true_nu_pz*true_nu_pz);

    // Reconstructed Shower momentum
    shr_p = std::sqrt(shr_px*shr_px + shr_py*shr_py + shr_pz*shr_pz);
    
    // Momentum of electron
    elec_mom = std::sqrt(elec_px*elec_px + elec_py*elec_py + elec_pz*elec_pz);

    // --

    TVector3 nu_p_vec(true_nu_px, true_nu_py, true_nu_pz);

    // True nue theta in BNB coordinates (up from beam dir)
    nu_theta = nu_p_vec.Theta() * 180.0/3.14159;
    
    // True nue phi in BNB coordinates (around beam dir)
    nu_phi = nu_p_vec.Phi()   * 180.0/3.14159;
    
    // True nue angle from numi beamline 
    nu_angle = _util.GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "beam"); 

    // True nue angle wrt numi target to uboone vector
    nu_angle_targ = _util.GetNuMIAngle(true_nu_px, true_nu_py, true_nu_pz, "target"); 

    // True electron angle wrt numi target to uboone vector
    elec_ang_targ = _util.GetNuMIAngle(elec_px, elec_py, elec_pz, "target");

    // --- Calculate the dot-product of the proxy for neutrino direction --- //
    // (the vector from the target to the reco nu vtx) and the true nu angle

    reco_true_nu_ang = nu_dir.Angle(v_targ_to_vtx) * 180/3.14159;
    // ---

    // The angle of the reconstructed shower relative to the NuMI target to detector direction
    shr_ang_numi = _util.GetNuMIAngle(shr_px, shr_py, shr_pz, "target");

}
// -----------------------------------------------------------------------------
void SliceContainer::CalibrateShowerEnergy(){

    // Divide the shower energy by 0.83 to calibrate it properly. 
    shr_energy_cali= shr_energy_cali/0.83;

}
// -----------------------------------------------------------------------------
void SliceContainer::SetThresholdEvent(){

    // Below threshold of nue or electron energy
    if (nu_e < _util.energy_threshold || elec_e < _util.elec_threshold){
        if (nu_pdg == 12){
            classification = std::make_pair("thr_nue",_util.k_thr_nue);
        }
        if (nu_pdg == -12){
            classification = std::make_pair("thr_nuebar",_util.k_thr_nuebar);
        }
    }

}
// -----------------------------------------------------------------------------
void SliceContainer::SetFakeData(){
    
    if (_util.isfakedata)
        classification = std::make_pair("data",_util.k_leg_data);
}