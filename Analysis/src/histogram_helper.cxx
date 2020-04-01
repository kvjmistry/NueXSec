#include "../include/histogram_helper.h"

// -----------------------------------------------------------------------------
histogram_helper::~histogram_helper() { 
    
    // Make sure the file is closed
    // f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_helper::MakeDirectory(){
        
    f_nuexsec->cd();

    bool bool_dir; // Check if directory exists already


    // Create subdirectory for cut type
    TDirectory *dir_plot_types[_util.plot_types.size()];
    
    // Create a new subdirectory for each cut
    const Int_t ncuts = _util.cut_dirs.size();
    TDirectory *dir_cut[ncuts];

    // Create a new subdirectory for each classification
    TDirectory *dir_classification[_util.classification_dirs.size()];
    
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
         
        f_nuexsec->cd();    // change current directory to top
    
    } // End loop over plot types ----------------------------------------------

    f_nuexsec->Write("",TObject::kOverwrite);
   
}
// -----------------------------------------------------------------------------
void histogram_helper::Initialise(int type, const char* run_period, const char * file_out, int weight_cfg ){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    std::string file_out_str = file_out;

    std::string file_name;

    if (type == _util.k_mc){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_mc_run%s.root", run_period);
        else file_name = "files/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( file_name.c_str() ) ) {
            f_nuexsec = new TFile( file_name.c_str(), "UPDATE");
        }
    }
    else if (type == _util.k_data){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_data_run%s.root", run_period);
        else file_name = "files/" + file_out_str;

        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

    }
    else if (type == _util.k_ext){

        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_ext_run%s.root", run_period);
        else file_name = "files/" + file_out_str;
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(file_name.c_str()) ) {
            f_nuexsec = new TFile(file_name.c_str(), "UPDATE");
        }

    }
    else if (type == _util.k_dirt){
        
        // If the file name is empty then we use the default file name
        if (file_out_str == "empty") file_name = Form("files/nuexsec_dirt_run%s.root", run_period);
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

    // Set the weight settings
    if (weight_cfg == 0){
        weight_tune = false;
        weight_ppfx = false;
    }
    else if (weight_cfg == 1){
        weight_tune = true;
        weight_ppfx = true;
    }
    else if (weight_cfg == 2){
        weight_tune = true;
        weight_ppfx = false;
    }
    else if (weight_cfg == 3){
        weight_tune = false;
        weight_ppfx = true;
    }
    else {
        std::cout << "Unknown weight setting specified, using defaults" << std::endl;
    }


    MakeDirectory();
}
// -----------------------------------------------------------------------------
void histogram_helper::InitHistograms(){
    
    
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
        
            // Reco Vtx X, Y, Z
            TH1D_hists.at(k_reco_vtx_x).at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 270);

            TH1D_hists.at(k_reco_vtx_y).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            
            TH1D_hists.at(k_reco_vtx_z).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            TH1D_hists.at(k_reco_vtx_x_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 270);

            TH1D_hists.at(k_reco_vtx_y_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            
            TH1D_hists.at(k_reco_vtx_z_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            // dEdx
            TH1D_hists.at(k_reco_dEdx_cali_y_plane).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_cali_y_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_dEdx_y_plane).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_y_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            // Leading Shower Momentum
            TH1D_hists.at(k_reco_leading_mom).at(i).at(j) = new TH1D ( Form("h_reco_leading_mom_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 2);

            // 2D distance shower vertex to reco nu vertex
            TH1D_hists.at(k_reco_shower_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // 2D distance track vertex to reco nu vertex
            TH1D_hists.at(k_reco_track_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_track_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 20);

            // Leading Shower hits in all planes
            TH1D_hists.at(k_reco_leading_shower_hits_all_planes).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_all_planes_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 2000);

            // Leading Shower hits in collection
            TH1D_hists.at(k_reco_leading_shower_hits_collection_plane).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_collection_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 800);

            // Leading Shower opening angle
            TH1D_hists.at(k_reco_leading_shower_open_angle).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_open_angle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 50);

            // Secondary shower to vertex distance (for events with more than 1 shower)
            TH1D_hists.at(k_reco_secondary_shower_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_secondary_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 80);

            // Leading Shower hits per length
            TH1D_hists.at(k_reco_leading_shower_hits_per_length).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_per_length_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // Longest track to leading shower length
            TH1D_hists.at(k_reco_longest_track_leading_shower_length).at(i).at(j) = new TH1D ( Form("h_reco_longest_track_leading_shower_length_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 3);

            // Track Containment
            TH1D_hists.at(k_reco_track_contained).at(i).at(j) = new TH1D ( Form("h_reco_track_contained_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 2, 0, 2);

            // Leading shower phi
            TH1D_hists.at(k_reco_leading_shower_phi).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_phi_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 12, -180, 180);

            // Leading shower theta
            TH1D_hists.at(k_reco_leading_shower_theta).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 12, 0, 180);

            // Leading shower cos theta
            TH1D_hists.at(k_reco_leading_shower_cos_theta).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_cos_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 16, -1, 1);

            // Leading shower multiplicity
            TH1D_hists.at(k_reco_shower_multiplicity).at(i).at(j) = new TH1D ( Form("h_reco_shower_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 1, 6);

            // Leading track multiplicity
            TH1D_hists.at(k_reco_track_multiplicity).at(i).at(j) = new TH1D ( Form("h_reco_track_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 0, 5);

            // Pandora topological score
            TH1D_hists.at(k_reco_topological_score).at(i).at(j) = new TH1D ( Form("h_reco_topological_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 1);

            // Track shower dist
            TH1D_hists.at(k_reco_track_shower_dist).at(i).at(j) = new TH1D (Form("h_reco_track_shower_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 10);
            
            // Track shower angle
            TH1D_hists.at(k_reco_track_shower_angle).at(i).at(j) = new TH1D (Form("h_reco_track_shower_angle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -60, 60);
            
            // Ratio hits from showers to slice
            TH1D_hists.at(k_reco_hits_ratio).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 1);
            
            // Shower score
            TH1D_hists.at(k_reco_shower_score).at(i).at(j) = new TH1D (Form("h_reco_shower_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 0.5);

            // Track score
            TH1D_hists.at(k_reco_track_score).at(i).at(j) = new TH1D (Form("h_reco_track_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0.5, 1);
            
            // Calibrated energy of all the showers
            TH1D_hists.at(k_reco_shower_energy_tot_cali).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_tot_cali_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 3);
            
            // Total number of hits for the leading shower
            TH1D_hists.at(k_reco_shower_hits).at(i).at(j) = new TH1D (Form("h_reco_shower_hits_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 800);
            
            // Total number of hits for the leading shower in the collection plane
            TH1D_hists.at(k_reco_shower_hits_y_plane).at(i).at(j) = new TH1D (Form("h_reco_shower_hits_y_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 800);


        }
        
    }
    // -------------------------------------------------------------------------

    // Intialising true histograms in here
    if (_type == _util.k_mc || _type == _util.k_dirt){
        
        if (_type == _util.k_mc){
            // Initalise the histograms for the TEfficency
            TEfficiency_hists.resize(_util.k_cuts_MAX);

            for (unsigned int l = 0; l < _util.k_cuts_MAX; l++ ){
                TEfficiency_hists.at(l) = new TH1D( Form("h_true_nu_E_%s",_util.cut_dirs.at(l).c_str() ), "", 10, 0, 4 );
            }
        }

        // Initalise the True Nue
        TH1D_true_hists.resize(k_TH1D_true_MAX);
        TH2D_true_hists.resize(k_TH2D_true_MAX);

        TH1D_true_hists.at(k_true_nue_theta) = new TH1D( Form("h_nue_true_theta_%s", _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Theta [degrees]; Entries",           14, 0, 180 );
        TH1D_true_hists.at(k_true_nue_phi)   = new TH1D( Form("h_true_nue_phi_%s",   _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Phi [degrees]; Entries",             14, -180, 180 );
        TH1D_true_hists.at(k_true_nue_angle) = new TH1D( Form("h_true_nue_angle_%s", _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Angle from NuMI [degrees]; Entries", 18, 0, 180 );
        TH1D_true_hists.at(k_true_nue_px)    = new TH1D( Form("h_true_nue_px_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Px [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_py)    = new TH1D( Form("h_true_nue_py_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Px [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_pz)    = new TH1D( Form("h_true_nue_pz_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Px [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_e)     = new TH1D( Form("h_true_nue_e_%s",     _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} E [GeV]; Entries",    25, 0, 5);
        TH1D_true_hists.at(k_true_nue_p)     = new TH1D( Form("h_true_nue_p_%s",     _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} P [GeV/c]; Entries",  14, 0, 5);
        TH1D_true_hists.at(k_true_vtx_x)     = new TH1D( Form("h_true_vtx_x_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx x [cm]; Entries", 20, -10, 270);
        TH1D_true_hists.at(k_true_vtx_y)     = new TH1D( Form("h_true_vtx_y_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx y [cm]; Entries", 20, -10, 120);
        TH1D_true_hists.at(k_true_vtx_z)     = new TH1D( Form("h_true_vtx_z_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx z [cm]; Entries", 40, -10, 1050);
        TH1D_true_hists.at(k_true_vtx_x_sce) = new TH1D( Form("h_true_vtx_x_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx x Space Charge Corr. [cm]; Entries", 20, -10, 270);
        TH1D_true_hists.at(k_true_vtx_y_sce) = new TH1D( Form("h_true_vtx_y_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx y Space Charge Corr. [cm]; Entries", 20, -10, 120);
        TH1D_true_hists.at(k_true_vtx_z_sce) = new TH1D( Form("h_true_vtx_z_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx z Space Charge Corr. [cm]; Entries", 40, -10, 1050);

        TH2D_true_hists.at(k_true_nue_theta_phi)    = new TH2D( Form("h_true_nue_theta_phi_%s", _util.type_prefix.at(_type).c_str()),   ";True #nu_{e} Theta [degrees]; True #nu_{e} Phi [degrees]",    14, 0, 180, 14, -180, 180 );
        TH2D_true_hists.at(k_true_nue_energy_theta) = new TH2D( Form("h_true_nue_energy_theta_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} E [GeV];True #nu_{e} Theta [degrees]",           25, 0, 5, 14, 0, 180);
        TH2D_true_hists.at(k_true_nue_energy_phi)   = new TH2D( Form("h_true_nue_energy_phi_%s", _util.type_prefix.at(_type).c_str()),  ";True #nu_{e} E [GeV];True #nu_{e} Phi [degrees]",             25, 0, 5, 14, -180, 180);
        TH2D_true_hists.at(k_true_nue_energy_angle) = new TH2D( Form("h_true_nue_energy_angle_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} E [GeV];True #nu_{e} Angle from NuMI [degrees]", 25, 0, 5, 18, 0, 180);
    
        TH2D_true_hists.at(k_true_nue_vtx_z_y)     = new TH2D( Form("h_true_nue_vtx_z_y_%s", _util.type_prefix.at(_type).c_str()),    ";True #nu_{e} Vtx Z [cm] ;True #nu_{e} Vtx Y [cm]", 40, -10, 1050, 20, -10, 120);
        TH2D_true_hists.at(k_true_nue_vtx_z_y_sce) = new TH2D( Form("h_true_nue_vtx_z_y_sce_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} Vtx Z  Space Charge Corr. [cm];True #nu_{e} Vtx Y Space Charge Corr. [cm]", 40, -10, 1050, 20, -10, 120);
    }

    // -------------------------------------------------------------------------

    // Flash histograms
    TH1D_flash_hists.resize(k_TH1D_flash_MAX);
    TH1D_flash_hists.at(k_flash_time) = new TH1D( Form("h_flash_time_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe)   = new TH1D( Form("h_flash_pe_%s", _util.type_prefix.at(_type).c_str()),   "; Flash PE [PE]; Entries", 100, 0, 5000 );

    TH1D_flash_hists.at(k_flash_time_sid1) = new TH1D( Form("h_flash_time_sid1_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Neutrino Candiate [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid1)   = new TH1D( Form("h_flash_pe_sid1_%s", _util.type_prefix.at(_type).c_str()),   "; Flash PE Neutrino Candiate [PE]; Entries", 100, 0, 5000 );

    TH1D_flash_hists.at(k_flash_time_sid0) = new TH1D( Form("h_flash_time_sid0_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Non Neutrino Canidate [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid0)   = new TH1D( Form("h_flash_pe_sid0_%s", _util.type_prefix.at(_type).c_str()),   "; Flash PE Non Neutrino Canidate [PE]; Entries", 100, 0, 5000 );
    

}
// -----------------------------------------------------------------------------
void histogram_helper::FillHists(int type, int classification_index, int cut_index, SliceContainer SC, double &_weight){

    // Get the CV weight
    double weight = GetCVWeight(type, SC);
    _weight = weight;

    // Calculate some variables
    double reco_shr_p = std::sqrt(SC.shr_px*SC.shr_px + SC.shr_py*SC.shr_py + SC.shr_pz*SC.shr_pz);

    // Now fill the histograms!
    TH1D_hists.at(k_reco_vtx_x).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x, weight);
    TH1D_hists.at(k_reco_vtx_y).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y, weight);
    TH1D_hists.at(k_reco_vtx_z).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z, weight);

    TH1D_hists.at(k_reco_vtx_x_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_x, weight);
    TH1D_hists.at(k_reco_vtx_y_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_y, weight);
    TH1D_hists.at(k_reco_vtx_z_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_z, weight);

    TH1D_hists.at(k_reco_dEdx_cali_y_plane).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_Y_cali, weight); // Just the collection plane!

    TH1D_hists.at(k_reco_dEdx_y_plane).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_Y, weight); // Just the collection plane!
    
    TH1D_hists.at(k_reco_leading_mom).at(cut_index).at(classification_index)->Fill(reco_shr_p, weight);
    
    TH1D_hists.at(k_reco_shower_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.shr_distance, weight);

    TH1D_hists.at(k_reco_track_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.trk_distance, weight);
    
    TH1D_hists.at(k_reco_leading_shower_hits_all_planes).at(cut_index).at(classification_index)->Fill(SC.shr_hits_tot, weight);
    
    TH1D_hists.at(k_reco_leading_shower_hits_collection_plane).at(cut_index).at(classification_index)->Fill(SC.shr_hits_y_tot, weight);

    // TH1D_hists.at(k_reco_secondary_shower_to_vtx_dist).at(cut_index).at(classification_index)->Fill();
    
    TH1D_hists.at(k_reco_leading_shower_open_angle).at(cut_index).at(classification_index)->Fill(SC.shr_openangle, weight);
    
    // TH1D_hists.at(k_reco_leading_shower_hits_per_length).at(cut_index).at(classification_index)->Fill();
    
    // TH1D_hists.at(k_reco_longest_track_leading_shower_length).at(cut_index).at(classification_index)->Fill();
    
    // TH1D_hists.at(k_reco_track_contained).at(cut_index).at(classification_index)->Fill();
    
    TH1D_hists.at(k_reco_leading_shower_phi).at(cut_index).at(classification_index)->Fill(SC.shr_phi * 180/3.14159, weight);
    
    TH1D_hists.at(k_reco_leading_shower_theta).at(cut_index).at(classification_index)->Fill(SC.shr_theta * 180/3.14159, weight);
    
    TH1D_hists.at(k_reco_leading_shower_cos_theta).at(cut_index).at(classification_index)->Fill(std::cos(SC.shr_theta * 180/3.14159), weight);
    
    TH1D_hists.at(k_reco_shower_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_showers, weight);
    
    TH1D_hists.at(k_reco_track_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_tracks, weight);

    TH1D_hists.at(k_reco_topological_score).at(cut_index).at(classification_index)->Fill(SC.topological_score, weight);

    TH1D_hists.at(k_reco_track_shower_dist).at(cut_index).at(classification_index)->Fill(SC.tksh_distance, weight);

    TH1D_hists.at(k_reco_track_shower_angle).at(cut_index).at(classification_index)->Fill(SC.tksh_angle*180/3.14159, weight);

    TH1D_hists.at(k_reco_hits_ratio).at(cut_index).at(classification_index)->Fill(SC.hits_ratio, weight);

    TH1D_hists.at(k_reco_shower_score).at(cut_index).at(classification_index)->Fill(SC.shr_score, weight);

    TH1D_hists.at(k_reco_track_score).at(cut_index).at(classification_index)->Fill(SC.trk_score, weight);

    TH1D_hists.at(k_reco_shower_energy_tot_cali).at(cut_index).at(classification_index)->Fill(SC.shr_energy_tot_cali, weight);

    TH1D_hists.at(k_reco_shower_hits).at(cut_index).at(classification_index)->Fill(SC.shr_hits_max, weight);

    TH1D_hists.at(k_reco_shower_hits_y_plane).at(cut_index).at(classification_index)->Fill(SC.shr_hits_y_tot, weight);

    // -------------------------------------------------------------------------
    // Fill truth histograms

    // Only do this for MC or dirt
    if ( (_type == _util.k_mc || _type == _util.k_dirt) && cut_index == _util.k_unselected){

        double p = std::sqrt(SC.true_nu_px*SC.true_nu_px + SC.true_nu_py*SC.true_nu_py + SC.true_nu_pz*SC.true_nu_pz);
        
        // True nue in BNB theta coordinates (up from beam dir)
        double nu_theta = acos(SC.true_nu_pz) * 180 / 3.1415;
        
        // True nue in BNB phi coordinates (around beam dir)
        double nu_phi = atan2(SC.true_nu_py, SC.true_nu_px) * 180 / 3.1415;
        
        // True nue angle from numi beamline 
        double nu_angle = _util.GetTheta(SC.true_nu_px, SC.true_nu_py, SC.true_nu_pz); 

        TH1D_true_hists.at(k_true_nue_theta)->Fill(nu_theta, weight);
        TH1D_true_hists.at(k_true_nue_phi)  ->Fill(nu_phi, weight);
        TH1D_true_hists.at(k_true_nue_angle)->Fill(nu_angle, weight);
        TH1D_true_hists.at(k_true_nue_px)   ->Fill(SC.true_nu_px, weight);
        TH1D_true_hists.at(k_true_nue_py)   ->Fill(SC.true_nu_py, weight);
        TH1D_true_hists.at(k_true_nue_pz)   ->Fill(SC.true_nu_pz, weight);
        TH1D_true_hists.at(k_true_nue_e)    ->Fill(SC.nu_e, weight);
        TH1D_true_hists.at(k_true_nue_p)    ->Fill(p, weight);
        TH1D_true_hists.at(k_true_vtx_x)    ->Fill(SC.true_nu_vtx_x, weight);
        TH1D_true_hists.at(k_true_vtx_y)    ->Fill(SC.true_nu_vtx_y, weight);
        TH1D_true_hists.at(k_true_vtx_z)    ->Fill(SC.true_nu_vtx_z, weight);
        TH1D_true_hists.at(k_true_vtx_x_sce)->Fill(SC.true_nu_vtx_sce_x, weight);
        TH1D_true_hists.at(k_true_vtx_y_sce)->Fill(SC.true_nu_vtx_sce_y, weight);
        TH1D_true_hists.at(k_true_vtx_z_sce)->Fill(SC.true_nu_vtx_sce_z, weight);

        TH2D_true_hists.at(k_true_nue_theta_phi)   ->Fill(nu_theta,  nu_phi, weight);
        TH2D_true_hists.at(k_true_nue_energy_theta)->Fill(SC.nu_e, nu_theta, weight);
        TH2D_true_hists.at(k_true_nue_energy_phi)  ->Fill(SC.nu_e, nu_phi, weight);
        TH2D_true_hists.at(k_true_nue_energy_angle)->Fill(SC.nu_e, nu_angle, weight);
        TH2D_true_hists.at(k_true_nue_vtx_z_y)       ->Fill(SC.true_nu_vtx_z,  SC.true_nu_vtx_y, weight);
        TH2D_true_hists.at(k_true_nue_vtx_z_y_sce)   ->Fill(SC.true_nu_vtx_sce_z,  SC.true_nu_vtx_sce_y, weight);
    }

    // -----------------------------------------------------------------------------

    // Only do this for the unselected histograms
    if (cut_index == _util.k_unselected){

        // Flash histograms
        if (type == _util.k_mc){
            TH1D_flash_hists.at(k_flash_time)->Fill(SC.flash_time + 0.055 -0.359, weight); // See numi documentation page to see what these numbers mean
            TH1D_flash_hists.at(k_flash_pe)->Fill(SC.flash_pe, weight); // See numi documentation page to see what these numbers mean
           
            if (SC.nslice == 1){
                TH1D_flash_hists.at(k_flash_time_sid1)->Fill(SC.flash_time + 0.055 -0.359, weight); // See numi documentation page to see what these numbers mean
                TH1D_flash_hists.at(k_flash_pe_sid1)->Fill(SC.flash_pe, weight); // See numi documentation page to see what these numbers mean

            }
            if (SC.nslice == 0){
                TH1D_flash_hists.at(k_flash_time_sid0)->Fill(SC.flash_time + 0.055 -0.359, weight); // See numi documentation page to see what these numbers mean
                TH1D_flash_hists.at(k_flash_pe_sid0)->Fill(SC.flash_pe, weight); // See numi documentation page to see what these numbers mean

            }
        }
        if (type == _util.k_dirt){
            TH1D_flash_hists.at(k_flash_time)->Fill(SC.flash_time + 0.055 -0.359, weight);
            TH1D_flash_hists.at(k_flash_pe)->Fill(SC.flash_pe, weight);
           
            if (SC.nslice == 1){
                TH1D_flash_hists.at(k_flash_time_sid1)->Fill(SC.flash_time + 0.055 -0.359, weight);
                TH1D_flash_hists.at(k_flash_pe_sid1)->Fill(SC.flash_pe, weight);

            }
            if (SC.nslice == 0){
                TH1D_flash_hists.at(k_flash_time_sid0)->Fill(SC.flash_time + 0.055 -0.359, weight);
                TH1D_flash_hists.at(k_flash_pe_sid0)->Fill(SC.flash_pe, weight);
                
            }
        }
        if (type == _util.k_ext){
            TH1D_flash_hists.at(k_flash_time)->Fill(SC.flash_time + -0.359, weight);
            TH1D_flash_hists.at(k_flash_pe)->Fill(SC.flash_pe, weight);
           
            if (SC.nslice == 1){
                TH1D_flash_hists.at(k_flash_time_sid1)->Fill(SC.flash_time + -0.359, weight);
                TH1D_flash_hists.at(k_flash_pe_sid1)->Fill(SC.flash_pe, weight);

            }
            if (SC.nslice == 0){
                TH1D_flash_hists.at(k_flash_time_sid0)->Fill(SC.flash_time + -0.359, weight);
                TH1D_flash_hists.at(k_flash_pe_sid0)->Fill(SC.flash_pe, weight);
                
            }
        }
        if (type == _util.k_data){
            TH1D_flash_hists.at(k_flash_time)->Fill(SC.flash_time, weight);
            TH1D_flash_hists.at(k_flash_pe)->Fill(SC.flash_pe, weight);
            
            if (SC.nslice == 1){
                TH1D_flash_hists.at(k_flash_time_sid1)->Fill(SC.flash_time, weight);
                TH1D_flash_hists.at(k_flash_pe_sid1)->Fill(SC.flash_pe, weight);

            }
            if (SC.nslice == 0){
                TH1D_flash_hists.at(k_flash_time_sid0)->Fill(SC.flash_time, weight);
                TH1D_flash_hists.at(k_flash_pe_sid0)->Fill(SC.flash_pe, weight);
            }
        }
        
    }
    

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteReco(int type){

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

                // Now write the histograms
                TH1D_hists.at(u).at(i).at(j)->SetOption("hist,E");
                TH1D_hists.at(u).at(i).at(j)->Write("",TObject::kOverwrite);

                if (break_early) break;
            }

        }
    }
}
// -----------------------------------------------------------------------------
void histogram_helper::FillTEfficiency(int cut_index, std::string classification, SliceContainer SC, double weight){

    // Fill the histogram at the specified cut
    if (classification == "nue_cc") TEfficiency_hists.at(cut_index)->Fill(SC.nu_e, weight);

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteTEfficiency(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"TEff");
    if (bool_dir) dir->cd();
    
    for (unsigned int p = 0; p < TEfficiency_hists.size(); p++){

        TEfficiency * teff = new TEfficiency(*TEfficiency_hists.at(p), *TEfficiency_hists.at(_util.k_unselected));
        teff->Write( Form("h_true_nu_E_%s",_util.cut_dirs.at(p).c_str()) , TObject::kOverwrite);
        TEfficiency_hists.at(p)->Write("",TObject::kOverwrite);
    }
    
}
// -----------------------------------------------------------------------------
void histogram_helper::WriteTrue(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"True");
    if (bool_dir) dir->cd();
    
    // TH1D
    for (unsigned int p = 0; p < TH1D_true_hists.size(); p++){
        TH1D_true_hists.at(p)->SetOption("hist,E");
        TH1D_true_hists.at(p)->Write("",TObject::kOverwrite);
    }
    
    // TH2D
    for (unsigned int p = 0; p < TH2D_true_hists.size(); p++){
        TH2D_true_hists.at(p)->SetOption("colz");
        TH2D_true_hists.at(p)->Write("",TObject::kOverwrite);
    }
    
}
// -----------------------------------------------------------------------------
void histogram_helper::WriteFlash(){
    
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
double histogram_helper::GetCVWeight(int type, SliceContainer SC){
    
    double weight = 1.0;

    // Get the tune weight
    if (type == _util.k_mc || type == _util.k_dirt){
        
        weight = SC.weightTune; // Here define the weight
        
        // Catch infinate/nan/unreasonably large tune weights
        if (std::isinf(weight))      weight = 1.0; 
        if (std::isnan(weight) == 1) weight = 1.0;
        if (weight > 100)            weight = 1.0;

        // If tune weight turned off, just set weight to 1.0
        if (!weight_tune) weight = 1.0;

    } 
    else weight = 1.0;
    

    // Get the PPFX CV flux correction weight
    if (type == _util.k_mc){
        double weight_flux = SC.GetPPFXCVWeight();
        
        // If we want ppfx weight then add this into the weight
        if (weight_ppfx) weight = weight * weight_flux;

    }

    return weight;

}
// -----------------------------------------------------------------------------
