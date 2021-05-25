#include "../include/HistogramHelper.h"

// -----------------------------------------------------------------------------
HistogramHelper::~HistogramHelper() { 
    
    // Make sure the file is closed
    // f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void HistogramHelper::MakeDirectory(){
        
    f_nuexsec->cd();

    bool bool_dir; // Check if directory exists already


    // Create subdirectory for cut type
    TDirectory *dir_plot_types[_util.plot_types.size()];
    
    // Create a new subdirectory for each cut
    const Int_t ncuts = _util.cut_dirs.size();
    TDirectory *dir_cut[ncuts];

    // Create a new subdirectory for each classification
    TDirectory *dir_classification[_util.classification_dirs.size()];

    TDirectory *dir_particle[_util.particle_types.size()];
    
    // Loop over the plot types ------------------------------------------------
    for (unsigned int k = 0; k < _util.plot_types.size(); k++) {
        
        // Get the directory 
        bool_dir = _util.GetDirectory(f_nuexsec, dir_plot_types[k] , _util.plot_types.at(k).c_str() );

        // Make the directory
        if (!bool_dir) dir_plot_types[k] = f_nuexsec->mkdir(_util.plot_types.at(k).c_str());

        dir_plot_types[k]->cd();

        // If we have stacked histograms, we make plots by cut 
        if (_util.plot_types.at(k) == "Stack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _util.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s", _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(_util.cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the classifications -------------------------------
                for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){
                
                    // Get the directory 
                    bool_dir = _util.GetDirectory(f_nuexsec, dir_classification[j] ,Form("%s/%s/%s", _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_classification[j] = dir_cut[i]->mkdir(_util.classification_dirs.at(j).c_str());
                    dir_classification[j]->cd();
                    f_nuexsec->cd();    // change current directory to top

                } // End loop over the classifications -------------------------
                
                f_nuexsec->cd();    // change current directory to top
                
            } // End loop over the cuts ----------------------------------------
       
        }
        // We also want stacked histograms by particle type
        else if (_util.plot_types.at(k) == "ParticleStack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _util.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s", _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(_util.cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the particle types -------------------------------
                for (unsigned int j = 0; j < _util.particle_types.size(); j++){
                
                    // Get the directory 
                    bool_dir = _util.GetDirectory(f_nuexsec, dir_particle[j] ,Form("%s/%s/%s", _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_particle[j] = dir_cut[i]->mkdir(_util.particle_types.at(j).c_str());
                    dir_particle[j]->cd();
                    f_nuexsec->cd();    // change current directory to top

                } // End loop over the classifications -------------------------
                
                f_nuexsec->cd();    // change current directory to top
                
            } // End loop over the cuts ----------------------------------------
       
        }
         
        f_nuexsec->cd();    // change current directory to top
    
    } // End loop over plot types ----------------------------------------------

    f_nuexsec->Write("",TObject::kOverwrite);
   
}
// -----------------------------------------------------------------------------
void HistogramHelper::Initialise(int type, const char * file_out, Utility util ){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    _util = util;

    std::string file_out_str = file_out;

    std::string file_name;

    if (type == _util.k_mc){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_mc_run%s.root", _util.run_period);
        else file_name = "files/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_nuexsec = new TFile( file_name.c_str(), "UPDATE");
        }
    }
    else if (type == _util.k_data){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_data_run%s.root", _util.run_period);
        else file_name = "files/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

    }
    else if (type == _util.k_ext){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_ext_run%s.root", _util.run_period);
        else file_name = "files/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

    }
    else if (type == _util.k_dirt){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_dirt_run%s.root", _util.run_period);
        else file_name = "files/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

    }
    else {
        std::cout << "Unknown input type!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }

    // Set the type
    _type = type;

    MakeDirectory();
}
// -----------------------------------------------------------------------------
void HistogramHelper::InitHistograms(){
    
    
    // Resize the histogram vector. plot var, cuts, classifications
    TH1D_hists.resize(k_TH1D_MAX);

    for (unsigned int u = 0; u < k_TH1D_MAX; u++){
        TH1D_hists.at(u).resize(_util.k_cuts_MAX);

        for (unsigned int i=0; i < _util.cut_dirs.size();i++){
            TH1D_hists.at(u).at(i).resize(_util.k_classifications_MAX);
        }
    }

    // loop over and create the histograms
    for (unsigned int i=0; i < _util.cut_dirs.size();i++){
    
        for (unsigned int j=0; j < _util.classification_dirs.size();j++){
        
            TH1D_hists.at(k_reco_pi0mass).at(i).at(j) = new TH1D ( Form("h_reco_pi0mass_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 500);

        }
        
    }

    // ----------
    // For pi0 histograms
    // Resize the histogram vector. plot var, cuts, classifications
    TH1D_pi0_hists.resize(k_TH1D_pi0_MAX);

    for (unsigned int u = 0; u < k_TH1D_pi0_MAX; u++){
        TH1D_pi0_hists.at(u).resize(_util.k_classifications_MAX);
    }

    // loop over and create the histograms
    for (unsigned int j=0; j < _util.classification_dirs.size();j++){
        TH1D_pi0_hists.at(k_pi0_mass).at(j) = new TH1D(Form("h_pi0_mass_%s", _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 500);
        TH1D_pi0_hists.at(k_pi0_mass_norm).at(j) = new TH1D(Form("h_pi0_mass_norm_%s", _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 500);
        TH1D_pi0_hists.at(k_pi0_mass_EScale).at(j) = new TH1D(Form("h_pi0_mass_EScale_%s", _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 500);
        TH1D_pi0_hists.at(k_pi0_energy).at(j) = new TH1D(Form("h_pi0_energy_%s", _util.classification_dirs.at(j).c_str()) ,"", 9, 60, 660);
        TH1D_pi0_hists.at(k_pi0_energy_norm).at(j) = new TH1D(Form("h_pi0_energy_norm_%s", _util.classification_dirs.at(j).c_str()) ,"", 9, 60, 660);
        TH1D_pi0_hists.at(k_pi0_energy_EScale).at(j) = new TH1D(Form("h_pi0_energy_EScale_%s", _util.classification_dirs.at(j).c_str()) ,"", 9, 60, 660);
    }    

}
// -----------------------------------------------------------------------------
void HistogramHelper::FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight){

    TH1D_hists.at(k_reco_pi0mass).at(cut_index).at(classification_index)->Fill(SC.pi0_mass_Y, weight);    

}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteReco(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth
    bool break_early{false};

    // Loop over the histogram variables
    for (unsigned int u = 0 ; u < TH1D_hists.size(); u++){
        
        // loop over the cut directories
        for (unsigned int i = 0; i < _util.cut_dirs.size(); i++){

            // loop over the classification directories
            for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){

                // Choose which folder to fill in based on the type (a re-mapping of enums)
                if (type == _util.k_mc && ( j == _util.k_leg_data || j == _util.k_leg_ext || j == _util.k_leg_dirt)){ 
                    break;
                }
                if (type == _util.k_data){ 
                    j = _util.k_leg_data;
                    break_early = true;
                }
                if (type == _util.k_ext){ 
                    j = _util.k_leg_ext;
                    break_early = true;
                }
                if (type == _util.k_dirt){ 
                    j= _util.k_leg_dirt;
                    break_early = true;
                }

                // Get the classification directory and cd
                bool_dir = _util.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s", "Stack", _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str() ) );
                
                if (bool_dir) truth_dir->cd();

                // Skip the write for the backgrounds in intrinsic nue mode
                // This overwrites the signal histograms
                if (std::string(_util.intrinsic_mode) == "intrinsic" && j != _util.k_nue_cc && j != _util.k_nuebar_cc &&
                                                                        j != _util.k_cosmic_nue && j != _util.k_cosmic_nuebar &&
                                                                        j != _util.k_unmatched_nue && j != _util.k_unmatched_nuebar &&
                                                                        j !=_util.k_thr_nue && j != _util.k_thr_nuebar){
                    continue;
                }

                // Now write the histograms
                TH1D_hists.at(u).at(i).at(j)->SetOption("hist,E");
                TH1D_hists.at(u).at(i).at(j)->Write("",TObject::kOverwrite);

                // Try to clear some memory
                // delete TH1D_hists.at(u).at(i).at(j);

                if (break_early) break;
            }

        }
    }

    // -----------------------------------------------------------------------------
    // -----------------------------------------------------------------------------

    // Do the same again for the 2D histograms

    // Loop over the histogram variables
    for (unsigned int u = 0 ; u < TH2D_hists_cuts.size(); u++){
        
        // loop over the cut directories
        for (unsigned int i = 0; i < _util.cut_dirs.size(); i++){

            // loop over the classification directories
            for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){

                // Choose which folder to fill in based on the type (a re-mapping of enums)
                if (type == _util.k_mc && ( j == _util.k_leg_data || j == _util.k_leg_ext || j == _util.k_leg_dirt)){ 
                    break;
                }
                if (type == _util.k_data){ 
                    j = _util.k_leg_data;
                    break_early = true;
                }
                if (type == _util.k_ext){ 
                    j = _util.k_leg_ext;
                    break_early = true;
                }
                if (type == _util.k_dirt){ 
                    j= _util.k_leg_dirt;
                    break_early = true;
                }

                // Get the classification directory and cd
                bool_dir = _util.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s", "Stack", _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str() ) );
                
                if (bool_dir) truth_dir->cd();

                // Skip the write for the backgrounds in intrinsic nue mode
                // This overwrites the signal histograms
                if (std::string(_util.intrinsic_mode) == "intrinsic" && j != _util.k_nue_cc && j != _util.k_nuebar_cc && 
                                                                        j != _util.k_cosmic_nue && j != _util.k_cosmic_nuebar && 
                                                                        j != _util.k_unmatched_nue && j != _util.k_unmatched_nuebar && 
                                                                        j != _util.k_thr_nue && j != _util.k_thr_nuebar  ){
                    continue;
                }

                // Now write the histograms
                TH2D_hists_cuts.at(u).at(i).at(j)->SetOption("colz");
                TH2D_hists_cuts.at(u).at(i).at(j)->Write("",TObject::kOverwrite);

                if (break_early) break;
            }

        }
    }

}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteRecoPar(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth
    bool break_early{false};

    // Loop over the histogram variables
    for (unsigned int u = 0 ; u < TH1D_hists_particle.size(); u++){
        
        // loop over the cut directories
        for (unsigned int i = 0; i < _util.cut_dirs.size(); i++){

            // loop over the particle type directories
            for (unsigned int j = 0; j < _util.particle_types.size(); j++){

                // Choose which folder to fill in based on the type (a re-mapping of enums)
                if (type == _util.k_mc && ( j == _util.k_part_data || j == _util.k_part_ext || j == _util.k_part_dirt)){ 
                    break;
                }
                if (type == _util.k_data){ 
                    j = _util.k_part_data;
                    break_early = true;
                }
                if (type == _util.k_ext){ 
                    j = _util.k_part_ext;
                    break_early = true;
                }
                if (type == _util.k_dirt){ 
                    j= _util.k_part_dirt;
                    break_early = true;
                }

                // Get the classification directory and cd
                bool_dir = _util.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s", "ParticleStack", _util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str() ) );
                
                if (bool_dir) truth_dir->cd();

                // Now write the histograms
                TH1D_hists_particle.at(u).at(i).at(j)->SetOption("hist,E");
                TH1D_hists_particle.at(u).at(i).at(j)->Write("",TObject::kOverwrite);

                if (break_early) break;
            }

        }
    }
}
// -----------------------------------------------------------------------------
void HistogramHelper::FillTEfficiency(int cut_index, std::string classification, SliceContainer SC, double weight){

    // Fill the histogram at the specified cut

    // If start of selection, efficiency denominator includes everything
    if (cut_index == _util.k_unselected){
        if (classification == "nue_cc"        || classification == "nuebar_cc" ||
            classification == "unmatched_nue" || classification == "unmatched_nuebar" ||
            classification == "cosmic_nue"    || classification == "cosmic_nuebar"){
            TEfficiency_hists.at(k_eff_nu_E).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E).at(cut_index)           ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_many_bins).at(cut_index) ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin).at(cut_index)     ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_single_bin).at(cut_index)  ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_flash_time).at(cut_index)    ->Fill(SC.flash_time + 0.055 -0.359, weight);
            TEfficiency_hists.at(k_eff_nu_theta).at(cut_index)         ->Fill(SC.nu_theta, weight);
            TEfficiency_hists.at(k_eff_nu_phi).at(cut_index)           ->Fill(SC.nu_phi, weight);
            TEfficiency_hists.at(k_eff_elec_theta).at(cut_index)       ->Fill(SC.elec_theta, weight);
            TEfficiency_hists.at(k_eff_elec_phi).at(cut_index)         ->Fill(SC.elec_phi, weight);
            TEfficiency_hists.at(k_eff_beta).at(cut_index)             ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_beta_rebin).at(cut_index)       ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta).at(cut_index)      ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
            TEfficiency_hists.at(k_eff_proton_multi).at(cut_index)     ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi).at(cut_index)       ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi).at(cut_index)  ->Fill(SC.nproton + SC.npion, weight);

        }
        
        // Nue only
        if (classification == "nue_cc" || classification == "unmatched_nue" || classification == "cosmic_nue") {
            TEfficiency_hists.at(k_eff_nu_E_nue).at(cut_index)              ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(cut_index)   ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin_nue).at(cut_index)      ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_proton_multi_nue).at(cut_index)      ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi_nue).at(cut_index)        ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi_nue).at(cut_index)   ->Fill(SC.nproton + SC.npion, weight);
            TEfficiency_hists.at(k_eff_beta_rebin_nue).at(cut_index)        ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin_nue).at(cut_index) ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
        }
        
        // Nuebar only
        if (classification == "nuebar_cc" || classification == "unmatched_nuebar" || classification == "cosmic_nuebar"){
            TEfficiency_hists.at(k_eff_nu_E_nuebar).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(cut_index)  ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin_nuebar).at(cut_index)     ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_proton_multi_nuebar).at(cut_index)     ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi_nuebar).at(cut_index)       ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi_nuebar).at(cut_index)  ->Fill(SC.nproton + SC.npion, weight);
            TEfficiency_hists.at(k_eff_beta_rebin_nuebar).at(cut_index)       ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin_nuebar).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
        }
        
    
    }
    // After this, we consider the low purity (cosmics) interactions background -- keep the unreconstructed stuff, but this shouldnt affect the plots all that much
    else {
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar"){
            TEfficiency_hists.at(k_eff_nu_E).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E).at(cut_index)           ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_many_bins).at(cut_index) ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin).at(cut_index)     ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_single_bin).at(cut_index)  ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_flash_time).at(cut_index)    ->Fill(SC.flash_time+ 0.055 -0.359, weight);
            TEfficiency_hists.at(k_eff_nu_theta).at(cut_index)         ->Fill(SC.nu_theta, weight);
            TEfficiency_hists.at(k_eff_nu_phi).at(cut_index)           ->Fill(SC.nu_phi, weight);
            TEfficiency_hists.at(k_eff_elec_theta).at(cut_index)       ->Fill(SC.elec_theta, weight);
            TEfficiency_hists.at(k_eff_elec_phi).at(cut_index)         ->Fill(SC.elec_phi, weight);
            TEfficiency_hists.at(k_eff_beta).at(cut_index)             ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_beta_rebin).at(cut_index)       ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta).at(cut_index)      ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
            TEfficiency_hists.at(k_eff_proton_multi).at(cut_index)     ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi).at(cut_index)       ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi).at(cut_index)  ->Fill(SC.nproton + SC.npion, weight);
        }
    
        // Nue only
        if (classification == "nue_cc" || classification == "unmatched_nue")       {
            TEfficiency_hists.at(k_eff_nu_E_nue).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(cut_index)  ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin_nue).at(cut_index)     ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_proton_multi_nue).at(cut_index)     ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi_nue).at(cut_index)       ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi_nue).at(cut_index)  ->Fill(SC.nproton + SC.npion, weight);
            TEfficiency_hists.at(k_eff_beta_rebin_nue).at(cut_index)       ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin_nue).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
        }
        
        // Nuebar only
        if (classification == "nuebar_cc" || classification == "unmatched_nuebar") {
            TEfficiency_hists.at(k_eff_nu_E_nuebar).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(cut_index)  ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin_nuebar).at(cut_index)     ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_proton_multi_nuebar).at(cut_index)     ->Fill(SC.nproton, weight);
            TEfficiency_hists.at(k_eff_pion_multi_nuebar).at(cut_index)       ->Fill(SC.npion, weight);
            TEfficiency_hists.at(k_eff_charg_par_multi_nuebar).at(cut_index)  ->Fill(SC.nproton + SC.npion, weight);
            TEfficiency_hists.at(k_eff_beta_rebin_nuebar).at(cut_index)       ->Fill(SC.true_effective_angle, weight);
            TEfficiency_hists.at(k_eff_cosine_beta_rebin_nuebar).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
        }
        
    }
    
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteTEfficiency(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"TEff");
    if (bool_dir) dir->cd();
    
    for (unsigned int v = 0; v < TEfficiency_hists.size(); v++){
        for (unsigned int p = 0; p < TEfficiency_hists.at(v).size(); p++){
            TEfficiency_hists.at(v).at(p)->Write("",TObject::kOverwrite);
        }
    }
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteTrue(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"True");
    if (bool_dir) dir->cd();
    
    // TH1D
    for (unsigned int p = 0; p < TH1D_true_hists.at(0).size(); p++){
        TH1D_true_hists.at(0).at(p)->SetOption("hist,E");
        TH1D_true_hists.at(0).at(p)->Write("",TObject::kOverwrite);
        TH1D_true_hists.at(1).at(p)->SetOption("hist,E");
        TH1D_true_hists.at(1).at(p)->Write("",TObject::kOverwrite);
    }
    
    // TH2D
    for (unsigned int p = 0; p < TH2D_true_hists.at(0).size(); p++){
        TH2D_true_hists.at(0).at(p)->SetOption("colz");
        TH2D_true_hists.at(0).at(p)->Write("",TObject::kOverwrite);
        TH2D_true_hists.at(1).at(p)->SetOption("colz");
        TH2D_true_hists.at(1).at(p)->Write("",TObject::kOverwrite);
    }
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteFlash(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"Flash");
    if (bool_dir) dir->cd();
    
    // TH1D
    for (unsigned int p = 0; p < TH1D_flash_hists.size(); p++){
        TH1D_flash_hists.at(p)->SetOption("hist,E");
        TH1D_flash_hists.at(p)->Write("",TObject::kOverwrite);
    }
    
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteInteractions(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"Interaction");
    if (bool_dir) dir->cd();
    
    // Loop over the histogram types
    for (unsigned int hist = 0; hist < TH1D_interaction_hists.size(); hist++){
        
        // Loop over the stages
        for (unsigned int stage = 0; stage < TH1D_interaction_hists.at(0).size(); stage++){

            // loop over the variables
            for (unsigned int p = 0; p < TH1D_interaction_hists.at(0).at(0).size(); p++){
                TH1D_interaction_hists.at(hist).at(stage).at(p)->SetOption("hist,E");
                TH1D_interaction_hists.at(hist).at(stage).at(p)->Write("",TObject::kOverwrite);
            }
        }
    }
    
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::Write_2DSigBkgHists(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"2D");
    if (bool_dir) dir->cd();
    
    // TH2D
    for (unsigned int p = 0; p < TH2D_hists.size(); p++){
        
        if (std::string(_util.intrinsic_mode) == "default" || std::string(_util.intrinsic_mode) == "intrinsic"){
            TH2D_hists.at(p).at(_util.k_signal)->SetOption("box");
            TH2D_hists.at(p).at(_util.k_signal)->SetFillColor(30);
            TH2D_hists.at(p).at(_util.k_signal)->Write("",TObject::kOverwrite);
        }

        if (std::string(_util.intrinsic_mode) == "default"){
            TH2D_hists.at(p).at(_util.k_background)->SetOption("box");
            TH2D_hists.at(p).at(_util.k_background)->SetFillColor(46);
            TH2D_hists.at(p).at(_util.k_background)->Write("",TObject::kOverwrite);
        }
    }
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::FillPiZeroHists(int classification_index, SliceContainer SC, double weight, int pizero_mode){

    if (pizero_mode == 0){
        TH1D_pi0_hists.at(k_pi0_mass).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
        TH1D_pi0_hists.at(k_pi0_energy).at(classification_index)->Fill(SC.pi0_energy1_Y, weight);
    }
    // Norm fix
    else if (pizero_mode == 1){
        TH1D_pi0_hists.at(k_pi0_mass_norm).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
        TH1D_pi0_hists.at(k_pi0_energy_norm).at(classification_index)->Fill(SC.pi0_energy1_Y, weight);
    }
    // Energy dependent
    else {
        TH1D_pi0_hists.at(k_pi0_mass_EScale).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
        TH1D_pi0_hists.at(k_pi0_energy_EScale).at(classification_index)->Fill(SC.pi0_energy1_Y, weight);
    }
    


}
// -----------------------------------------------------------------------------
void HistogramHelper::FillNuMuHists(int classification_index, SliceContainer SC, double weight){

    
    TH1D_numu_hists.at(k_track_theta).at(classification_index)->Fill(SC.trk_theta* 180/3.14159, weight);
    TH1D_numu_hists.at(k_track_phi).at(classification_index)->Fill(SC.trk_phi* 180/3.14159, weight);
    TH1D_numu_hists.at(k_track_cos_theta).at(classification_index)->Fill(std::cos(SC.trk_theta), weight);
    TH1D_numu_hists.at(k_muon_topo_score).at(classification_index)->Fill(SC.topological_score, weight);
    
    


}
// -----------------------------------------------------------------------------
void HistogramHelper::WritePiZero(int type){
    f_nuexsec->cd();

    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"pizero");
    if (bool_dir) dir->cd();

    bool break_early{false};
    
    // Loop over the histogram variables
    for (unsigned int u = 0 ; u < TH1D_pi0_hists.size(); u++){
        
        // loop over the classification directories
        for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){

            // Choose whether to fill MC type classifications or data/ext/dirt (a re-mapping of enums)
            if (type == _util.k_mc && ( j == _util.k_leg_data || j == _util.k_leg_ext || j == _util.k_leg_dirt)){ 
                break;
            }
            if (type == _util.k_data){ 
                j = _util.k_leg_data;
                break_early = true;
            }
            if (type == _util.k_ext){ 
                j = _util.k_leg_ext;
                break_early = true;
            }
            if (type == _util.k_dirt){ 
                j= _util.k_leg_dirt;
                break_early = true;
            }

            // Now write the histograms
            TH1D_pi0_hists.at(u).at(j)->SetOption("hist,E");
            TH1D_pi0_hists.at(u).at(j)->Write("",TObject::kOverwrite);

            if (break_early) break;
        }

        
    }

}
// -----------------------------------------------------------------------------
void HistogramHelper::WriteNuMu(int type){
    f_nuexsec->cd();

    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"numu");
    if (bool_dir) dir->cd();

    bool break_early{false};
    
    // Loop over the histogram variables
    for (unsigned int u = 0 ; u < TH1D_numu_hists.size(); u++){
        
        // loop over the classification directories
        for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){

            // Choose whether to fill MC type classifications or data/ext/dirt (a re-mapping of enums)
            if (type == _util.k_mc && ( j == _util.k_leg_data || j == _util.k_leg_ext || j == _util.k_leg_dirt)){ 
                break;
            }
            if (type == _util.k_data){ 
                j = _util.k_leg_data;
                break_early = true;
            }
            if (type == _util.k_ext){ 
                j = _util.k_leg_ext;
                break_early = true;
            }
            if (type == _util.k_dirt){ 
                j= _util.k_leg_dirt;
                break_early = true;
            }

            // Now write the histograms
            TH1D_numu_hists.at(u).at(j)->SetOption("hist,E");
            TH1D_numu_hists.at(u).at(j)->Write("",TObject::kOverwrite);

            if (break_early) break;
        }
        
    }

}