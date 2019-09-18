#include "../include/histogram_helper.h"

// -----------------------------------------------------------------------------
histogram_helper::~histogram_helper() { 
    
    // Make sure the file is closed
    f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_helper::MakeDirectory(std::string type){
        
    f_nuexsec->cd();

    TDirectory *top_dir; // e.g MC, Data, EXT
    bool bool_dir;       // Check if directory exists already
    TString type_tstr = type;
   
    // Create the top directory
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, top_dir, type_tstr );
    if (!bool_dir) top_dir = f_nuexsec->mkdir(type.c_str());
    
    // Make the the top dir the current directory
    top_dir->cd();

    // Create subdirectory for cut type
    TDirectory *dir_plot_types[plot_types.size()];
    
    // Create a new subdirectory for each cut
    const Int_t ncuts = cut_dirs.size();
    TDirectory *dir_cut[ncuts];

    // Create a new subdirectory for each classification
    TDirectory *dir_classification[classification_dirs.size()];
    
    // Loop over the plot types ------------------------------------------------
    for (unsigned int k = 0; k < plot_types.size(); k++) {
        
        // Get the directory 
        bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_plot_types[k] ,Form("%s/%s", type.c_str(), plot_types.at(k).c_str()) );

        // Make the directory
        if (!bool_dir) dir_plot_types[k] = top_dir->mkdir(plot_types.at(k).c_str());

        dir_plot_types[k]->cd();

        // If we have stacked histograms, we make plots by cut 
        if (plot_types.at(k) == "Stack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the classifications -------------------------------
                for (unsigned int j = 0; j < classification_dirs.size(); j++){
                
                    // Get the directory 
                    bool_dir = _utility_instance.GetDirectory(f_nuexsec, dir_classification[j] ,Form("%s/%s/%s/%s", type.c_str(), plot_types.at(k).c_str(), cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_classification[j] = dir_cut[i]->mkdir(classification_dirs.at(j).c_str());
                    dir_classification[j]->cd();
                    top_dir->cd();    // change current directory to top

                } // End loop over the classifications -------------------------
                
                top_dir->cd();    // change current directory to top
                
            } // End loop over the cuts ----------------------------------------
       
        }
         
        top_dir->cd();    // change current directory to top
    
    } // End loop over plot types ----------------------------------------------

    top_dir->Write("",TObject::kOverwrite);
   
}
// -----------------------------------------------------------------------------
void histogram_helper::Initialise(){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject("nuexsec.root") ) {
        f_nuexsec = new TFile("nuexsec.root", "UPDATE");
    }

    MakeDirectory("MC");
    MakeDirectory("Data");
    MakeDirectory("EXT");
    MakeDirectory("Dirt");
}
// -----------------------------------------------------------------------------
void histogram_helper::InitHistograms(){
    
    // Flash Histograms
    h_flash_time_v.resize(k_flash_MAX);
    for (unsigned int i=0; i < h_flash_time_v.size();i++){
        h_flash_time_v.at(i) = new TH1D ( Form("h_flash_time_%s", type_prefix.at(i).c_str()) ,"", 80, 0, 20);
    }

    // Reco Vtx X
    h_reco_vtx_x.resize(k_cut_dirs_MAX);
    h_reco_vtx_y.resize(k_cut_dirs_MAX);
    h_reco_vtx_z.resize(k_cut_dirs_MAX);
    
    for (unsigned int i=0; i < cut_dirs.size();i++){

        h_reco_vtx_x.at(i).resize(k_classifications_MAX);
        h_reco_vtx_y.at(i).resize(k_classifications_MAX);
        h_reco_vtx_z.at(i).resize(k_classifications_MAX);

        for (unsigned int j=0; j < classification_dirs.size();j++){
            h_reco_vtx_x.at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 270);
            h_reco_vtx_y.at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            h_reco_vtx_z.at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 1050);
        }
    }

}
// -----------------------------------------------------------------------------
void histogram_helper::FillMCTruth( double mc_nu_energy,  double mc_nu_momentum,  int mc_nu_id, bool in_tpc,
                                    double mc_nu_vtx_x,   double mc_nu_vtx_y,  double mc_nu_vtx_z,
                                    double mc_nu_dir_x,   double mc_nu_dir_y,  double mc_nu_dir_z,
                                    double mc_ele_dir_x,  double mc_ele_dir_y, double mc_ele_dir_z,
                                    double mc_ele_energy, double mc_ele_momentum ) {
    
    double mc_cos_theta     = -999;
    double theta            = -999;
    double mc_ele_cos_theta = -999;
    double mc_ele_theta     = -999;

    // Caclulate theta and cos(theta)
    if (mc_nu_momentum != 0) {
        mc_cos_theta = mc_nu_dir_z;
        theta        = acos(mc_cos_theta) * (180 / 3.1415);
    }
    
    // Caclulate theta and cos(theta) for the electron
    if (mc_ele_momentum != 0) {
        mc_ele_cos_theta = mc_ele_dir_z;
        mc_ele_theta     = acos(mc_ele_cos_theta) * (180/3.1415);
    }
    
    // Calculate Phi
    double phi        = atan2(mc_nu_dir_y, mc_nu_dir_x) * (180/3.1415);
    double mc_ele_phi = atan2(mc_ele_dir_y, mc_ele_dir_x) * (180/3.1415);
   
    if ((mc_nu_id == 1 || mc_nu_id == 5) && in_tpc == true){
        
        // 1D Hists
        h_nue_true_theta->Fill(theta);
        h_nue_true_phi  ->Fill(phi);
        
        // 2D Hists
        h_nue_true_theta_phi   ->Fill(phi, theta );
        h_nue_true_energy_theta->Fill(mc_nu_momentum, theta);
        h_nue_true_energy_phi  ->Fill(mc_nu_momentum, phi);
        
        h_ele_true_energy_theta->Fill(mc_ele_energy, mc_ele_theta);
        h_ele_true_energy_phi  ->Fill(mc_ele_energy, mc_ele_phi);
    }
    

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteMCTruth(std::string type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // Get the truth directory and cd
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s", type.c_str(), "Truth") );
    if (bool_dir) truth_dir->cd();

    // Now write the histograms
    h_nue_true_theta->Write("",TObject::kOverwrite);
    h_nue_true_phi  ->Write("",TObject::kOverwrite);
    
    // 2D Hists
    h_nue_true_theta_phi   ->Write("",TObject::kOverwrite);
    h_nue_true_energy_theta->Write("",TObject::kOverwrite);
    h_nue_true_energy_phi  ->Write("",TObject::kOverwrite);
    
    h_ele_true_energy_theta->Write("",TObject::kOverwrite);
    h_ele_true_energy_phi  ->Write("",TObject::kOverwrite);


}
// -----------------------------------------------------------------------------
void histogram_helper::FillOptical(std::vector<std::vector<double>> optical_list_flash_time_v, int type){

    // Loop over the optical list events
    for (unsigned int i = 0; i < optical_list_flash_time_v.size();i++){

        // Loop over the flashes in each event
        for (unsigned int j = 0; j < optical_list_flash_time_v.at(i).size();j++){
            h_flash_time_v.at(type)->Fill(optical_list_flash_time_v.at(i).at(j));
        }
    }

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteOptical(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // Get the Optical directory and cd
    bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s", type_prefix.at(type).c_str(), "Optical") );
    if (bool_dir) truth_dir->cd();

    // Now write the histograms
    h_flash_time_v.at(type)->Write("",TObject::kOverwrite);
    
   
}
// -----------------------------------------------------------------------------
int histogram_helper::IndexOfClassification(std::string tpco_id){

    // Nue CC
    if (tpco_id == "nue_cc_qe"  || tpco_id == "nue_bar_cc_qe"  ||
        tpco_id == "nue_cc_res" || tpco_id == "nue_bar_cc_res" ||
        tpco_id == "nue_cc_dis" || tpco_id == "nue_bar_cc_dis" || 
        tpco_id == "nue_cc_coh" || tpco_id == "nue_bar_cc_coh" || 
        tpco_id == "nue_cc_mec" || tpco_id == "nue_bar_cc_mec") {
        return k_nue_cc;
    }

    // Nue CC OOFV
    if (tpco_id == "nue_cc_out_fv"){
        return k_nue_cc_out_fv;
    }

    // NuMu CC
    if (tpco_id == "numu_cc_qe" || tpco_id == "numu_cc_res" ||
       tpco_id == "numu_cc_dis" || tpco_id == "numu_cc_coh" || 
       tpco_id == "numu_cc_mec" || tpco_id == "numu_cc_mixed"){
        return k_numu_cc;
    }

    // NC
    if (tpco_id == "nc"){
        return k_nc;
    }
    
    // NC pi0
    if (tpco_id == "nc_pi0"){
        return k_nc_pi0;
    }

    // Nue CC Mixed
    if (tpco_id == "nue_cc_mixed"){
        return k_nue_cc_mixed;
    }

    // Cosmic
    if (tpco_id == "cosmic"){
        return k_cosmic;
    }

    // Other Mixed
    if (tpco_id == "other_mixed"){
        return k_nc_mixed;
    }

    // Unmatched
    if (tpco_id == "unmatched" || tpco_id == "bad_reco"){
        return k_unmatched;
    }

    // Data
    if (tpco_id == "Data"){
        return  k_leg_data;
    }

    // EXT/ In time cosmics
    if (tpco_id == "EXT"){
        return  k_leg_ext;
    }

    std::cout << "Reached end of index of classifcation, this is BAD!!!! " << tpco_id << std::endl;
    return k_classifications_MAX;
}
// -----------------------------------------------------------------------------
void histogram_helper::FillRecoVtx(int classification_index, int cut_index, const xsecAna::TPCObjectContainer &tpc_obj){

    h_reco_vtx_x.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxX());
    h_reco_vtx_y.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxY());
    h_reco_vtx_z.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxZ());

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteRecoVtx(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth

    // loop over the cut directories
    for (unsigned int i = 0; i < cut_dirs.size(); i++){
        
        // loop over the classification directories
        for (unsigned int j = 0; j < classification_dirs.size(); j++){

            // Get the classification directory and cd
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s/%s", type_prefix.at(type).c_str(), "Stack", cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str() ) );
            if (bool_dir) truth_dir->cd();

            // Now write the histogram
            h_reco_vtx_x.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_y.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_z.at(i).at(j)->Write("",TObject::kOverwrite);
        }

    }
    
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------