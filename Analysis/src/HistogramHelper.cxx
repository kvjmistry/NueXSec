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
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // Stacked histograms by particle type

    // Resize the histogram vector. plot var, cuts, classifications
    TH1D_hists_particle.resize(k_TH1D_par_MAX);

    for (unsigned int u = 0; u < k_TH1D_par_MAX; u++){
        TH1D_hists_particle.at(u).resize(_util.k_cuts_MAX);

        for (unsigned int i=0; i < _util.cut_dirs.size();i++){
            TH1D_hists_particle.at(u).at(i).resize(_util.k_particles_MAX);
        }
    }

    // loop over and create the histograms
    for (unsigned int i=0; i < _util.cut_dirs.size();i++){
    
        for (unsigned int j=0; j < _util.particle_types.size();j++){
        
            // dEdx
            TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_max_par).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_max_par_%s_%s",_util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_y_par).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_par_%s_%s", _util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str()) ,"", 40, 0, 10);

            // Track PID score in the event
            TH1D_hists_particle.at(k_reco_trk_pid_score_par).at(i).at(j) = new TH1D ( Form("h_reco_trk_pid_score_par_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, -1, 1);
        }
    }


    // -------------------------------------------------------------------------


    // Intialising true histograms in here
    if (_type == _util.k_mc){
        
        if (_type == _util.k_mc){
            // Initalise the histograms for the TEfficency
            TEfficiency_hists.resize(k_TH1D_eff_MAX);
            for (unsigned int v = 0; v < TEfficiency_hists.size(); v++){
                TEfficiency_hists.at(v).resize(_util.k_cuts_MAX);
            }
            
            for (unsigned int l = 0; l < _util.k_cuts_MAX; l++ ){
                TEfficiency_hists.at(k_eff_nu_E).at(l)                     = new TH1D( Form("h_true_nu_E_%s",                    _util.cut_dirs.at(l).c_str() ), "", 8, 0, 6 );
                TEfficiency_hists.at(k_eff_nu_E_many_bins).at(l)           = new TH1D( Form("h_true_nu_E_many_bins_%s",          _util.cut_dirs.at(l).c_str() ), "", 100, 0, 1 );
                TEfficiency_hists.at(k_eff_elec_E).at(l)                   = new TH1D( Form("h_true_elec_E_%s",                  _util.cut_dirs.at(l).c_str() ), "", 8, 0, 6 );
                TEfficiency_hists.at(k_eff_elec_E_many_bins).at(l)         = new TH1D( Form("h_eff_elec_E_many_bins_%s",         _util.cut_dirs.at(l).c_str() ), "", 100, 0, 1 );
                
                double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
                TEfficiency_hists.at(k_eff_elec_E_rebin).at(l)             = new TH1D( Form("h_true_elec_E_rebin_%s",            _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_elec_E_rebin_nue).at(l)         = new TH1D( Form("h_true_elec_E_rebin_nue_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_elec_E_rebin_nuebar).at(l)      = new TH1D( Form("h_true_elec_E_rebin_nuebar_%s",     _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_elec_E_rebin_pi0).at(l)         = new TH1D( Form("h_true_elec_E_rebin_pi0_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_nu_E_nue).at(l)                 = new TH1D( Form("h_true_nu_E_nue_%s",                _util.cut_dirs.at(l).c_str() ), "", 8, 0, 6 );
                TEfficiency_hists.at(k_eff_nu_E_nuebar).at(l)              = new TH1D( Form("h_true_nu_E_nuebar_%s",             _util.cut_dirs.at(l).c_str() ), "", 8, 0, 6 );
                TEfficiency_hists.at(k_eff_nu_E_single_bin).at(l)          = new TH1D( Form("h_true_nu_E_single_bin_%s",         _util.cut_dirs.at(l).c_str() ), "", 1, 0, 10);
                TEfficiency_hists.at(k_eff_nu_E_nue_single_bin).at(l)      = new TH1D( Form("h_true_nu_E_nue_single_bin_%s",     _util.cut_dirs.at(l).c_str() ), "", 1, 0, 10 );
                TEfficiency_hists.at(k_eff_nu_E_nuebar_single_bin).at(l)   = new TH1D( Form("h_true_nu_E_nuebar_single_bin_%s",  _util.cut_dirs.at(l).c_str() ), "", 1, 0, 10 );
                TEfficiency_hists.at(k_eff_nu_flash_time).at(l)            = new TH1D( Form("h_eff_nu_flash_time_%s",            _util.cut_dirs.at(l).c_str() ), "", 5, 6.0, 15.0);
                TEfficiency_hists.at(k_eff_nu_theta).at(l)                 = new TH1D( Form("h_eff_nu_theta_%s",                 _util.cut_dirs.at(l).c_str() ), "", 13, 20, 100 );
                TEfficiency_hists.at(k_eff_nu_phi).at(l)                   = new TH1D( Form("h_eff_nu_phi_%s",                   _util.cut_dirs.at(l).c_str() ), "", 14, 0, 70 );
                TEfficiency_hists.at(k_eff_elec_theta).at(l)               = new TH1D( Form("h_eff_elec_theta_%s",               _util.cut_dirs.at(l).c_str() ), "", 13, 0, 190 );
                TEfficiency_hists.at(k_eff_elec_phi).at(l)                 = new TH1D( Form("h_eff_elec_phi_%s",                 _util.cut_dirs.at(l).c_str() ), "", 14, -190, 190 );
                
                TEfficiency_hists.at(k_eff_beta).at(l)                     = new TH1D( Form("h_eff_beta_%s",                     _util.cut_dirs.at(l).c_str() ), "", 13, 0, 190 );
                edges = &_util.reco_shr_bins_ang[0]; // Cast to an array 
                TEfficiency_hists.at(k_eff_beta_rebin).at(l)               = new TH1D( Form("h_eff_beta_rebin_%s",               _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_ang.size()-1, edges );
                TEfficiency_hists.at(k_eff_beta_rebin_nue).at(l)           = new TH1D( Form("h_eff_beta_rebin_nue_%s",           _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_ang.size()-1, edges );
                TEfficiency_hists.at(k_eff_beta_rebin_nuebar).at(l)        = new TH1D( Form("h_eff_beta_rebin_nuebar_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_ang.size()-1, edges );
                
                TEfficiency_hists.at(k_eff_cosine_beta).at(l)              = new TH1D( Form("h_eff_cosine_beta_%s",              _util.cut_dirs.at(l).c_str() ), "", 5, 0.6, 1 );
                TEfficiency_hists.at(k_eff_cosine_beta_nue).at(l)          = new TH1D( Form("h_eff_cosine_beta_nue_%s",          _util.cut_dirs.at(l).c_str() ), "", 5, 0.6, 1 );
                TEfficiency_hists.at(k_eff_cosine_beta_nuebar).at(l)       = new TH1D( Form("h_eff_cosine_beta_nuebar_%s",       _util.cut_dirs.at(l).c_str() ), "", 5, 0.6, 1 );
                edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
                TEfficiency_hists.at(k_eff_cosine_beta_rebin).at(l)        = new TH1D( Form("h_eff_cosine_beta_rebin_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                TEfficiency_hists.at(k_eff_cosine_beta_rebin_nue).at(l)    = new TH1D( Form("h_eff_cosine_beta_rebin_nue_%s",    _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                TEfficiency_hists.at(k_eff_cosine_beta_rebin_nuebar).at(l) = new TH1D( Form("h_eff_cosine_beta_rebin_nuebar_%s", _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                TEfficiency_hists.at(k_eff_cosine_beta_rebin_pi0).at(l)    = new TH1D( Form("h_eff_cosine_beta_rebin_pi0_%s",    _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                
                TEfficiency_hists.at(k_eff_proton_multi).at(l)             = new TH1D( Form("h_eff_proton_multi_%s",             _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_proton_multi_nue).at(l)         = new TH1D( Form("h_eff_proton_multi_nue_%s",         _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_proton_multi_nuebar).at(l)      = new TH1D( Form("h_eff_proton_multi_nuebar_%s",      _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_pion_multi).at(l)               = new TH1D( Form("h_eff_pion_multi_%s",               _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_pion_multi_nue).at(l)           = new TH1D( Form("h_eff_pion_multi_nue_%s",           _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_pion_multi_nuebar).at(l)        = new TH1D( Form("h_eff_pion_multi_nuebar_%s",        _util.cut_dirs.at(l).c_str() ), "", 6, 0, 6 );
                TEfficiency_hists.at(k_eff_charg_par_multi).at(l)          = new TH1D( Form("h_eff_charg_par_multi_%s",          _util.cut_dirs.at(l).c_str() ), "", 8, 0, 8 );
                TEfficiency_hists.at(k_eff_charg_par_multi_nue).at(l)      = new TH1D( Form("h_eff_charg_par_multi_nue_%s",      _util.cut_dirs.at(l).c_str() ), "", 8, 0, 8 );
                TEfficiency_hists.at(k_eff_charg_par_multi_nuebar).at(l)   = new TH1D( Form("h_eff_charg_par_multi_nuebar_%s",   _util.cut_dirs.at(l).c_str() ), "", 8, 0, 8 );
            }
            
        }

        // Initalise the True Nue
        TH1D_true_hists.resize(2); // Unselected and Selected
        TH2D_true_hists.resize(2); 
        
        for (unsigned int i = 0; i < TH1D_true_hists.size(); i++){
            TH1D_true_hists.at(i).resize(k_TH1D_true_MAX);
            TH2D_true_hists.at(i).resize(k_TH2D_true_MAX);
        }
        
        for (unsigned int i = 0; i < TH1D_true_hists.size(); i++){
            std::string cut_stage = "unselected";
            if (i == 1) cut_stage = "selected";            
            TH1D_true_hists.at(i).at(k_true_nue_theta) = new TH1D( Form("h_true_nue_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} #theta [deg]; Entries",           14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_nue_phi)   = new TH1D( Form("h_true_nue_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} #phi [deg]; Entries",              14, 0, 40);
            
            TH1D_true_hists.at(i).at(k_true_nue_angle) = new TH1D( Form("h_true_nue_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} Angle [deg]; Entries", 120, 0, 120 );
            TH1D_true_hists.at(i).at(k_true_nue_px)    = new TH1D( Form("h_true_nue_px_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} P_x [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_py)    = new TH1D( Form("h_true_nue_py_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} P_y [GeV/c]; Entries", 14, 0, 1);
            TH1D_true_hists.at(i).at(k_true_nue_pz)    = new TH1D( Form("h_true_nue_pz_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} P_z [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_e)     = new TH1D( Form("h_true_nue_e_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} E [GeV]; Entries",    50, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_p)     = new TH1D( Form("h_true_nue_p_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} P [GeV/c]; Entries",  14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_vtx_x)     = new TH1D( Form("h_true_vtx_x_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx x [cm]; Entries", 40, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y)     = new TH1D( Form("h_true_vtx_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx y [cm]; Entries", 40, -120, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z)     = new TH1D( Form("h_true_vtx_z_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx z [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_vtx_x_sce) = new TH1D( Form("h_true_vtx_x_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx x SCE Corr. [cm]; Entries", 40, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y_sce) = new TH1D( Form("h_true_vtx_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx y SCE Corr. [cm]; Entries", 40, -120, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z_sce) = new TH1D( Form("h_true_vtx_z_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} or #bar{#nu}_{e} Vtx z SCE Corr. [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_nu_ang_targ)     = new TH1D( Form("h_true_nu_ang_targ_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} or #bar{#nu}_{e} Angle wrt NuMI Targ Vec [deg]; Entries",  40, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_ang_targ)   = new TH1D( Form("h_true_elec_ang_targ_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} or e^{+} Angle wrt NuMI Targ Vec [deg]; Entries", 25, 0, 180 );

            TH1D_true_hists.at(i).at(k_true_elec_E)      = new TH1D( Form("h_true_elec_E_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} or e^{+} Energy [GeV]; Entries",                15, 0, 5 );
            TH1D_true_hists.at(i).at(k_true_elec_theta)  = new TH1D( Form("h_true_elec_theta_%s_%s",  _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} or e^{+} #theta [deg]; Entries",            14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_phi)    = new TH1D( Form("h_true_elec_phi_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} or e^{+} #phi [deg]; Entries",              25, -180, 180);
            TH1D_true_hists.at(i).at(k_true_elec_gamma)  = new TH1D( Form("h_true_elec_gamma_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} or e^{+} NuMI #Phi [deg]; Entries",        18, -180, 180);

            TH1D_true_hists.at(i).at(k_reco_true_ang)    = new TH1D( Form("h_reco_true_ang_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), "; Angle (Reco to True) Nu Dir [deg]; Entries",              25, -2, 20);

            TH2D_true_hists.at(i).at(k_true_nue_phi_theta)    = new TH2D( Form("h_true_nue_phi_theta_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} or #bar{#nu}_{e} #phi [deg];True #nu_{e} #theta [deg]",     14, 0, 80, 16, 20, 140 );
            TH2D_true_hists.at(i).at(k_true_nue_energy_theta) = new TH2D( Form("h_true_nue_energy_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} or #bar{#nu}_{e} E [GeV];True #nu_{e} #theta [deg]",           25, 0, 5, 16, 0, 140);
            TH2D_true_hists.at(i).at(k_true_nue_energy_phi)   = new TH2D( Form("h_true_nue_energy_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} or #bar{#nu}_{e} E [GeV];True #nu_{e} #phi [deg]",             25, 0, 5, 14, 0, 75);
            TH2D_true_hists.at(i).at(k_true_nue_energy_angle) = new TH2D( Form("h_true_nue_energy_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} or #bar{#nu}_{e} E [GeV];True #nu_{e} or #bar{#nu}_{e} Angle [deg]",  20, 0, 3, 33, 7, 40);
        
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y)     = new TH2D( Form("h_true_nue_vtx_z_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Vtx Z [cm] ;True #nu_{e} or #bar{#nu}_{e} Vtx Y [cm]", 40, -10, 1050, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y_sce) = new TH2D( Form("h_true_nue_vtx_z_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),";True #nu_{e} or #bar{#nu}_{e} Vtx Z  SCE Corr. [cm];True #nu_{e} or #bar{#nu}_{e} Vtx Y SCE Corr. [cm]", 40, -10, 1050, 20, -10, 120);
        
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E) = new TH2D( Form("h_true_elec_E_reco_elec_E_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e#lower[-0.5]{-} or e^{+} Energy [GeV] ;Reco e#lower[-0.5]{-} or e^{+} Energy [GeV] [GeV]", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins.size()-1, edges);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E)     = new TH2D( Form("h_true_nu_E_reco_nu_E_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Energy [GeV] ;Reco #nu_{e} or #bar{#nu}_{e} Energy [GeV]", 25, 0, 4, 25, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_%s_%s",               _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e#lower[-0.5]{-} or e^{+} Energy [GeV] ;Reco e#lower[-0.5]{-} or e^{+} Energy [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins_nue) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_nue_%s_%s",       _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e#lower[-0.5]{-} Energy [GeV] ;Reco e#lower[-0.5]{-} Energy  [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins_nuebar) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_nuebar_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e^{+} Energy [GeV] ;Reco e^{+} Energy  [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E_extra_bins)     = new TH2D( Form("h_true_nu_E_reco_nu_E_extra_bins_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Energy [GeV] ;Reco #nu_{e} or #bar{#nu}_{e} Energy [GeV]", 50, 0, 4, 50, 0, 4);

            TH2D_true_hists.at(i).at(k_true_nu_vtx_x_reco_nu_vtx_x)  = new TH2D( Form("h_true_nu_vtx_x_reco_nu_vtx_x_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Vtx X [cm] ;Reco #nu_{e} or #bar{#nu}_{e} Vtx X [cm]", 20, -10, 270, 20, -10, 270);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_y_reco_nu_vtx_y)  = new TH2D( Form("h_true_nu_vtx_y_reco_nu_vtx_y_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Vtx Y [cm] ;Reco #nu_{e} or #bar{#nu}_{e} Vtx Y [cm]", 20, -10, 120, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_z_reco_nu_vtx_z)  = new TH2D( Form("h_true_nu_vtx_z_reco_nu_vtx_z_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} or #bar{#nu}_{e} Vtx Z [cm] ;Reco #nu_{e} or #bar{#nu}_{e} Vtx Z [cm]", 40, -10, 1050, 40, -10, 1050);
            TH2D_true_hists.at(i).at(k_true_shr_energy_purity)       = new TH2D( Form("h_true_shr_energy_purity_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e} [GeV] ;Shower Purity", _util.reco_shr_bins.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_energy_completeness) = new TH2D( Form("h_true_shr_energy_completeness_%s_%s",_util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e} [GeV] ;Shower Completeness", _util.reco_shr_bins.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_energy_resolution_reco) = new TH2D( Form("h_true_shr_energy_resolution_reco_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e} [GeV]; Reco - True / Reco", _util.reco_shr_bins.size()-1, edges, 30, -1.2, 1.2);
            TH2D_true_hists.at(i).at(k_true_shr_energy_resolution_true) = new TH2D( Form("h_true_shr_energy_resolution_true_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e} [GeV] ;Reco - True / True", _util.reco_shr_bins.size()-1, edges, 30, -1.2, 1.2);
        
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_rebin) = new TH2D( Form("h_true_elec_E_reco_elec_E_rebin_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e#lower[-0.5]{-} or e^{+} Energy [GeV] ;Reco. e#lower[-0.5]{-} or e^{+} Energy [GeV]", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins.size()-1, edges);


            double* edges_cbeta = &_util.reco_shr_bins_cang[0]; // Cast to an array 
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_purity)          = new TH2D( Form("h_true_shr_cosbeta_purity_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e} ;Shower Purity", _util.reco_shr_bins_cang.size()-1, edges_cbeta, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_completeness)    = new TH2D( Form("h_true_shr_cosbeta_completeness_%s_%s",_util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e} ;Shower Completeness", _util.reco_shr_bins_cang.size()-1, edges_cbeta, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_resolution_reco) = new TH2D( Form("h_true_shr_cosbeta_resolution_reco_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e}; Reco - True", _util.reco_shr_bins_cang.size()-1, edges_cbeta, 30, -1.2, 1.2);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_resolution_true) = new TH2D( Form("h_true_shr_cosbeta_resolution_true_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e} ;Reco - True", _util.reco_shr_bins_cang.size()-1, edges_cbeta, 30, -1.2, 1.2);

            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta)      = new TH2D( Form("h_elec_true_beta_reco_beta_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),            ";True e#lower[-0.5]{-} or e^{+} #beta [deg]; Reco e#lower[-0.5]{-} or e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta_nue)    = new TH2D( Form("h_elec_true_beta_reco_beta_nue_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),      ";True e#lower[-0.5]{-} #beta [deg]; Reco e#lower[-0.5]{-} or e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta_nuebar) = new TH2D( Form("h_elec_true_beta_reco_beta_nuebar_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e^{+} #beta [deg]; Reco e#lower[-0.5]{-} or e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_theta_reco_theta)    = new TH2D( Form("h_elec_true_theta_reco_theta_%s_%s",  _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e#lower[-0.5]{-} or e^{+} #theta [deg];Reco e#lower[-0.5]{-} or e^{+} #theta [deg]", 36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_phi_reco_phi)        = new TH2D( Form("h_elec_true_phi_reco_phi_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e#lower[-0.5]{-} or e^{+} #phi [deg];  Reco e#lower[-0.5]{-} or e^{+} #phi [deg]",   36, -180, 180, 36, -180, 180 );
        
            TH2D_true_hists.at(i).at(k_elec_true_cosbeta_reco_cosbeta_rebin) = new TH2D( Form("h_elec_true_cosbeta_reco_cosbeta_rebin_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e#lower[-0.5]{-} or e^{+} cos#beta ;Reco. e#lower[-0.5]{-} or e^{+} cos#beta", _util.reco_shr_bins_cang.size()-1, edges_cbeta, _util.reco_shr_bins_cang.size()-1, edges_cbeta);

            TH2D_true_hists.at(i).at(k_true_elec_E_true_beta_rebin) = new TH2D( Form("h_true_elec_E_true_beta_rebin_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e#lower[-0.5]{-} or e^{+} Energy [GeV] ;True e#lower[-0.5]{-} or e^{+} cos#beta", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins_cang.size()-1, edges_cbeta);
        }
    }

    // -------------------------------------------------------------------------

    // Flash histograms
    TH1D_flash_hists.resize(k_TH1D_flash_MAX);
    TH1D_flash_hists.at(k_flash_time) = new TH1D( Form("h_flash_time_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time [us]; Entries", 50, 0, 25 );
    TH1D_flash_hists.at(k_flash_time_single_bin) = new TH1D( Form("h_flash_time_single_bin_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time [us]; Entries", 1, 5.6, 15.4 );
    TH1D_flash_hists.at(k_flash_pe)   = new TH1D( Form("h_flash_pe_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity [PE]; Entries", 40, 0, 10000 );

    TH1D_flash_hists.at(k_flash_time_sid1) = new TH1D( Form("h_flash_time_sid1_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Neutrino Candiate [us]; Entries", 50, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid1)   = new TH1D( Form("h_flash_pe_sid1_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity (Neutrino Candiate) [PE]; Entries", 40, 0, 10000 );

    TH1D_flash_hists.at(k_flash_time_sid0) = new TH1D( Form("h_flash_time_sid0_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Non Neutrino Canidate [us]; Entries", 50, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid0)   = new TH1D( Form("h_flash_pe_sid0_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity (Non Neutrino Canidate) [PE]; Entries", 40, 0, 10000 );
    
    // Interaction Histograms
    TH1D_interaction_hists.resize(k_INTERACTION_MAX); // All the interaction histograms
    
    for (unsigned int l=0; l < TH1D_interaction_hists.size(); l++){
        TH1D_interaction_hists.at(l).resize(2); // Unselected and Selected
    }
    
    // Interaction types
    for (unsigned int l=0; l < TH1D_interaction_hists.size(); l++){
        for (unsigned int i = 0; i < TH1D_interaction_hists.at(l).size(); i++){
            TH1D_interaction_hists.at(l).at(i).resize(_util.k_interactions_MAX);
        }
    }
    

    // Loop over the unselected/selected
    for (int i = 0; i < 2; i++){

        std::string stage = "unselected";
        
        if (i == 1)
            stage = "selected";
    
        // Loop over the interaction types
        for (unsigned int p =0 ; p < TH1D_interaction_hists.at(0).at(0).size(); p++ ){
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 

            TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(i).at(p)        = new TH1D( Form("h_true_nue_nuebar_E_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e} Energy [GeV]; Entries",   15, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_nue).at(i).at(p)               = new TH1D( Form("h_true_nue_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} Energy [GeV]; Entries",                  15, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(i).at(p)            = new TH1D( Form("h_true_nuebar_E_%s_%s",              _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{#nu}_{e} Energy [GeV]; Entries",            15, 0, 6 );
            // TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(i).at(p)        = new TH1D( Form("h_true_nue_nuebar_E_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e} Energy [GeV]; Entries",   8, 0, 6 );
            // TH1D_interaction_hists.at(k_int_nu_E_nue).at(i).at(p)               = new TH1D( Form("h_true_nue_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} Energy [GeV]; Entries",                  8, 0, 6 );
            // TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(i).at(p)            = new TH1D( Form("h_true_nuebar_E_%s_%s",              _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{#nu}_{e} Energy [GeV]; Entries",            8, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(i).at(p)        = new TH1D( Form("h_int_nu_E_single_bin_%s_%s",        _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e}; Entries",                               1, 0, 10 );
            TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(i).at(p)    = new TH1D( Form("h_int_nu_E_nue_single_bin_%s_%s",    _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{nu}_{e}; Entries",                          1, 0, 10 );
            TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(i).at(p) = new TH1D( Form("h_int_nu_E_nuebar_single_bin_%s_%s", _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e}; Entries",                1, 0, 10 );
            TH1D_interaction_hists.at(k_int_elec_E).at(i).at(p)                 = new TH1D( Form("h_int_elec_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e#lower[-0.5]{-}/e^{+} Energy [GeV]; Entries",      15, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_nue).at(i).at(p)             = new TH1D( Form("h_int_elec_E_nue_%s_%s",             _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e#lower[-0.5]{-} [GeV]; Entries",                     15, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_nuebar).at(i).at(p)          = new TH1D( Form("h_int_elec_E_nuebar_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e^{+} Energy [GeV]; Entries",                         15, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_rebin).at(i).at(p)           = new TH1D( Form("h_int_elec_E_rebin_%s_%s",           _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} Energy [GeV]; Entries", _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_E_rebin_nue).at(i).at(p)       = new TH1D( Form("h_int_elec_E_rebin_nue_%s_%s",       _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} Energy [GeV]; Entries",         _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar).at(i).at(p)    = new TH1D( Form("h_int_elec_E_rebin_nuebar_%s_%s",    _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e^{+} Energy [GeV]; Entries",                    _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_theta).at(i).at(p)             = new TH1D( Form("h_int_elec_theta_%s_%s",             _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} Energy [GeV]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_elec_phi).at(i).at(p)               = new TH1D( Form("h_int_elec_phi_%s_%s",               _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} Energy [GeV]; Entries", 14, -190, 190 );
            TH1D_interaction_hists.at(k_int_effective_ang).at(i).at(p)          = new TH1D( Form("h_int_effective_ang_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} Energy [GeV]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_beta_nue).at(i).at(p)               = new TH1D( Form("h_int_beta_nue_%s_%s",               _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} #beta [deg]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_beta_nuebar).at(i).at(p)            = new TH1D( Form("h_int_beta_nuebar_%s_%s",            _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} #beta [deg]; Entries", 13, 0, 190 );
            edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
            TH1D_interaction_hists.at(k_int_cosbeta).at(i).at(p)                = new TH1D( Form("h_int_cosbeta_%s_%s",                _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-}/e^{+} cos#beta; Entries", _util.reco_shr_bins_cang.size()-1, edges );
            TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(i).at(p)      = new TH1D( Form("h_int_cosbeta_rebin_nue_%s_%s",      _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} cos#beta; Entries",         _util.reco_shr_bins_cang.size()-1, edges );
            TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(i).at(p)   = new TH1D( Form("h_int_cosbeta_rebin_nuebar_%s_%s",   _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e^{+} cos#beta; Entries",                    _util.reco_shr_bins_cang.size()-1, edges );
        }
    }
        
    // 2D signal and background separation plots
    
    // Resizing the histograms
    TH2D_hists.resize(k_TH2D_reco_MAX);
    
    // Loop over the histograms to resize
    for (unsigned int i = 0; i < TH2D_hists.size(); i++){
        TH2D_hists.at(i).resize(_util.k_sig_bkg_MAX);
    }

    // Loop over the types 
    for (unsigned int k =0 ; k< TH2D_hists.at(0).size(); k++){
        // dEdx vs shower vertex distance 
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(k)              = new TH2D( Form("h_reco_shr_dEdx_shr_dist_%s",            _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist_post).at(k)         = new TH2D( Form("h_reco_shr_dEdx_shr_dist_post_%s",       _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist).at(k)          = new TH2D( Form("h_reco_shr_dEdx_max_shr_dist_%s",        _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (All Planes) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist_post).at(k)     = new TH2D( Form("h_reco_shr_dEdx_max_shr_dist_post_%s",   _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (All Planes) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist_large_dedx).at(k)   = new TH2D( Form("h_reco_shr_dEdx_shr_dist_large_dedx_%s", _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 80, 10, 30, 15, 0, 20);
    
        // dEdx vs Moliere Average
        TH2D_hists.at(k_reco_shr_dEdx_moliere).at(k)          = new TH2D( Form("h_reco_shr_dEdx_moliere_%s", _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dE/dx (Collection Plane) [MeV/cm];Leading Shower Moliere Avg [deg]", 40, 0, 10, 20, 0, 30);
    
        // Moliere Average vs Shower vtx distance
        TH2D_hists.at(k_reco_shr_moliere_shr_dist).at(k)      = new TH2D( Form("h_reco_shr_moliere_shr_dist_%s", _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower Moliere Avg [deg];Leading Shower to Vertex Distance [cm]", 20, 0, 30, 15, 0, 20);
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

    // ----------
    // For numu histograms
    // Resize the histogram vector. plot var, cuts, classifications
    TH1D_numu_hists.resize(k_TH1D_numu_MAX);

    for (unsigned int u = 0; u < k_TH1D_numu_MAX; u++){
        TH1D_numu_hists.at(u).resize(_util.k_classifications_MAX);
    }

    // loop over and create the histograms
    for (unsigned int j=0; j < _util.classification_dirs.size();j++){
        TH1D_numu_hists.at(k_track_theta).at(j)     = new TH1D(Form("h_track_theta_%s", _util.classification_dirs.at(j).c_str()) ,"", 13, 0, 190 );
        TH1D_numu_hists.at(k_track_cos_theta).at(j) = new TH1D(Form("h_track_cos_theta_%s", _util.classification_dirs.at(j).c_str()) ,"", 16, -1, 1);
        TH1D_numu_hists.at(k_track_phi).at(j)       = new TH1D(Form("h_track_phi_%s", _util.classification_dirs.at(j).c_str()) ,"", 14, -190, 190);
        TH1D_numu_hists.at(k_muon_topo_score).at(j) = new TH1D(Form("h_muon_topo_score_%s", _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 1);
        
    }

    // ----------
    // For 2D histograms broken down by plot var, cut, classification

    // Resize the histogram vector. plot var, cuts, classifications
    TH2D_hists_cuts.resize(k_TH2D_cut_MAX);

    for (unsigned int u = 0; u < k_TH2D_cut_MAX; u++){
        TH2D_hists_cuts.at(u).resize(_util.k_cuts_MAX);

        for (unsigned int i=0; i < _util.cut_dirs.size();i++){
            TH2D_hists_cuts.at(u).at(i).resize(_util.k_classifications_MAX);
        }
    }

    // loop over and create the histograms
    for (unsigned int i=0; i < _util.cut_dirs.size();i++){
    
        for (unsigned int j=0; j < _util.classification_dirs.size();j++){
        
            // Shower dEdx vs Shower energy
            TH2D_hists_cuts.at(k_2D_dedx_shower_energy).at(i).at(j) = new TH2D ( Form("h_2D_dedx_shower_energy_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10, 10, 0, 4);

        }
        
    }
    

}
// -----------------------------------------------------------------------------
void HistogramHelper::FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight){

    // Check if the interaction was in the FV
    //bool true_in_fv = _util.in_fv(SC.true_nu_vtx_x, SC.true_nu_vtx_y, SC.true_nu_vtx_z);

    // Get the trkfit dedx max variable
    //double dedx_max = SC.GetdEdxMax();
    
    // Now fill the histograms!
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
            TEfficiency_hists.at(k_eff_nu_E_many_bins).at(cut_index)   ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E).at(cut_index)           ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_many_bins).at(cut_index) ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin).at(cut_index)     ->Fill(SC.elec_e, weight);
            if (SC.npi0 > 0) TEfficiency_hists.at(k_eff_elec_E_rebin_pi0).at(cut_index) ->Fill(SC.elec_e, weight);
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
            if (SC.npi0 > 0)TEfficiency_hists.at(k_eff_cosine_beta_rebin_pi0).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
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
            TEfficiency_hists.at(k_eff_cosine_beta_nue).at(cut_index)       ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
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
            TEfficiency_hists.at(k_eff_cosine_beta_nuebar).at(cut_index)      ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
        }
        
    
    }
    // After this, we consider the low purity (cosmics) interactions background -- keep the unreconstructed stuff, but this shouldnt affect the plots all that much
    else {
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar"){
            TEfficiency_hists.at(k_eff_nu_E).at(cut_index)             ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_nu_E_many_bins).at(cut_index)   ->Fill(SC.nu_e, weight);
            TEfficiency_hists.at(k_eff_elec_E).at(cut_index)           ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_many_bins).at(cut_index) ->Fill(SC.elec_e, weight);
            TEfficiency_hists.at(k_eff_elec_E_rebin).at(cut_index)     ->Fill(SC.elec_e, weight);
            if (SC.npi0 > 0) TEfficiency_hists.at(k_eff_elec_E_rebin_pi0).at(cut_index) ->Fill(SC.elec_e, weight);
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
            if (SC.npi0 > 0)TEfficiency_hists.at(k_eff_cosine_beta_rebin_pi0).at(cut_index)->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
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
            TEfficiency_hists.at(k_eff_cosine_beta_nue).at(cut_index)      ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
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
            TEfficiency_hists.at(k_eff_cosine_beta_nuebar).at(cut_index)      ->Fill(std::cos(SC.true_effective_angle *  3.14159 / 180.0), weight);
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
