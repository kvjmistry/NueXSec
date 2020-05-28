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

            TH1D_hists.at(k_reco_vtx_y).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -120, 120);
            
            TH1D_hists.at(k_reco_vtx_z).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            TH1D_hists.at(k_reco_vtx_x_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 270);

            TH1D_hists.at(k_reco_vtx_y_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -120, 120);
            
            TH1D_hists.at(k_reco_vtx_z_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            // dEdx
            TH1D_hists.at(k_reco_dEdx_cali_u_plane).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_cali_u_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_dEdx_cali_v_plane).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_cali_v_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_dEdx_cali_y_plane).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_cali_y_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            // Leading Shower Momentum
            TH1D_hists.at(k_reco_leading_mom).at(i).at(j) = new TH1D ( Form("h_reco_leading_mom_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 2);

            // 2D distance shower vertex to reco nu vertex
            TH1D_hists.at(k_reco_shower_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // 2D distance track vertex to reco nu vertex
            TH1D_hists.at(k_reco_track_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_track_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 20);

            // Leading Shower hits in all planes
            TH1D_hists.at(k_reco_shr_hits_max).at(i).at(j) = new TH1D ( Form("h_reco_shr_hits_max_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 2000);

            // Number of Tracks Contained
            TH1D_hists.at(k_reco_n_track_contained).at(i).at(j) = new TH1D ( Form("h_reco_n_track_contained_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 0, 5);

            // Number of Showers
            TH1D_hists.at(k_reco_n_shower_contained).at(i).at(j) = new TH1D ( Form("h_reco_n_shower_contained_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 0, 5);

            // Leading shower phi
            TH1D_hists.at(k_reco_leading_shower_phi).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_phi_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 14, -190, 190);

            // Leading shower theta
            TH1D_hists.at(k_reco_leading_shower_theta).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 13, 0, 190);

            // Leading shower cos theta
            TH1D_hists.at(k_reco_leading_shower_cos_theta).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_cos_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 16, -1, 1);

            // Leading shower multiplicity
            TH1D_hists.at(k_reco_shower_multiplicity).at(i).at(j) = new TH1D ( Form("h_reco_shower_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 6, 0, 6);

            // Leading track multiplicity
            TH1D_hists.at(k_reco_track_multiplicity).at(i).at(j) = new TH1D ( Form("h_reco_track_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 6, 0, 6);

            // Pandora topological score
            TH1D_hists.at(k_reco_topological_score).at(i).at(j) = new TH1D ( Form("h_reco_topological_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 1);

            // Track shower dist
            TH1D_hists.at(k_reco_track_shower_dist).at(i).at(j) = new TH1D (Form("h_reco_track_shower_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 10);
            
            // Track shower angle
            TH1D_hists.at(k_reco_track_shower_angle).at(i).at(j) = new TH1D (Form("h_reco_track_shower_angle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -180, 180);
            
            // Ratio hits from showers to slice
            TH1D_hists.at(k_reco_hits_ratio).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);
            
            // Shower score
            TH1D_hists.at(k_reco_shower_score).at(i).at(j) = new TH1D (Form("h_reco_shower_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 0.5);

            // Track score
            TH1D_hists.at(k_reco_track_score).at(i).at(j) = new TH1D (Form("h_reco_track_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0.5, 1);
            
            // Calibrated energy of all the showers
            TH1D_hists.at(k_reco_shower_energy_tot_cali).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_tot_cali_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 3);
            
            // Total number of hits for the leading shower
            TH1D_hists.at(k_reco_shr_hits_tot).at(i).at(j) = new TH1D (Form("h_reco_shr_hits_tot_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 1600);
            
            // Total number of hits for the leading shower in the collection plane
            TH1D_hists.at(k_reco_shr_hits_y_tot).at(i).at(j) = new TH1D (Form("h_reco_shr_hits_y_tot_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 1600);

            // Shower Track fit dEdx variables
            TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_u).at(i).at(j)   = new TH1D ( Form("h_reco_shr_trkfit_2cm_dEdx_u_%s_%s",_util.cut_dirs.at(i).c_str(),     _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_v).at(i).at(j)   = new TH1D ( Form("h_reco_shr_trkfit_2cm_dEdx_v_%s_%s",_util.cut_dirs.at(i).c_str(),     _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_y).at(i).at(j)   = new TH1D ( Form("h_reco_shr_trkfit_2cm_dEdx_y_%s_%s",_util.cut_dirs.at(i).c_str(),     _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_u).at(i).at(j)  = new TH1D ( Form("h_reco_shr_trkfit_gap05_dEdx_u_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_v).at(i).at(j)  = new TH1D ( Form("h_reco_shr_trkfit_gap05_dEdx_v_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_y).at(i).at(j)  = new TH1D ( Form("h_reco_shr_trkfit_gap05_dEdx_y_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_u).at(i).at(j) = new TH1D ( Form("h_reco_shr_trkfit_gap10_dEdx_u_%s_%s",_util.cut_dirs.at(i).c_str(),   _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_v).at(i).at(j) = new TH1D ( Form("h_reco_shr_trkfit_gap10_dEdx_v_%s_%s",_util.cut_dirs.at(i).c_str(),   _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_y).at(i).at(j) = new TH1D ( Form("h_reco_shr_trkfit_gap10_dEdx_y_%s_%s",_util.cut_dirs.at(i).c_str(),   _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            // Op Filter
            TH1D_hists.at(k_reco_opfilter_beam).at(i).at(j) = new TH1D ( Form("h_reco_opfilter_beam_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 500);
            TH1D_hists.at(k_reco_opfilter_veto).at(i).at(j) = new TH1D ( Form("h_reco_opfilter_veto_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 75, 0, 200);

            // Software Trigger
            TH1D_hists.at(k_reco_softwaretrig).at(i).at(j) = new TH1D ( Form("h_reco_softwaretrig_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 2, 0, 2);
        
            // Slice ID
            TH1D_hists.at(k_reco_nslice).at(i).at(j) = new TH1D ( Form("h_reco_nslice_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 2, 0, 2);

            // Slice Cluster Fraction
            TH1D_hists.at(k_reco_slclustfrac).at(i).at(j) = new TH1D ( Form("h_reco_slclustfrac_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 15, 0, 1);
        
            // Cosmic Inpact Parameter
            TH1D_hists.at(k_reco_cosmicIP).at(i).at(j)       = new TH1D ( Form("h_reco_cosmicIP_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 75, 0, 250);
            TH1D_hists.at(k_reco_CosmicIPAll3D).at(i).at(j)  = new TH1D ( Form("h_reco_CosmicIPAll3D_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 75, 0, 250);
            TH1D_hists.at(k_reco_CosmicDirAll3D).at(i).at(j) = new TH1D ( Form("h_reco_CosmicDirAll3D_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -1, 1);

            // dEdx with the trackfit variable
            TH1D_hists.at(k_reco_shr_tkfit_dedx_u).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_u_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_v).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_v_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y_good_theta).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_good_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_v_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_v_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_u_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_u_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            
            // Flash time
            TH1D_hists.at(k_reco_flash_time).at(i).at(j) = new TH1D ( Form("h_reco_flash_time_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 25 );

            // Flash PE
            TH1D_hists.at(k_reco_flash_pe).at(i).at(j) = new TH1D ( Form("h_reco_flash_pe_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 100, 0, 5000 );

            // Shower Subcluster in each plane
            TH1D_hists.at(k_reco_shrsubclusters).at(i).at(j)  = new TH1D ( Form("h_reco_shrsubclusters_%s_%s", _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // Shower Moliers Average
            TH1D_hists.at(k_reco_shrmoliereavg).at(i).at(j) = new TH1D ( Form("h_reco_shrmoliereavg_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 30);

            // Shower Moliere RMS
            TH1D_hists.at(k_reco_shrmoliererms).at(i).at(j) = new TH1D ( Form("h_reco_shrmoliererms_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 30000 );

            // Shower Cylindrical Fraction 1cm in the Second half
            TH1D_hists.at(k_reco_CylFrac2h_1cm).at(i).at(j) = new TH1D ( Form("h_reco_CylFrac2h_1cm_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 0.5 );

            // Shower RMS of spacepoint to the shower center at the second half of the Shower
            TH1D_hists.at(k_reco_DeltaRMS2h).at(i).at(j) = new TH1D ( Form("h_reco_DeltaRMS2h_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 20 );

            // Shower Median PCA Calculated in 5cm blocks
            TH1D_hists.at(k_reco_shrPCA1CMed_5cm).at(i).at(j) = new TH1D ( Form("h_reco_shrPCA1CMed_5cm_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 1 );

            // Shower Momentum from MCS
            TH1D_hists.at(k_reco_shrMCSMom).at(i).at(j) = new TH1D ( Form("h_reco_shrMCSMom_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 300 );

            // Distance between the neutrino vertex and (closest?) cosmic trajectory tagged from CRT
            TH1D_hists.at(k_reco_closestNuCosmicDist).at(i).at(j) = new TH1D ( Form("h_reco_closestNuCosmicDist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 400);
            
            // Longest Track Length if there is one
            TH1D_hists.at(k_reco_trk_len).at(i).at(j) = new TH1D ( Form("h_reco_trk_len_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 400);

            // Reconstructed Neutrino Energy
            TH1D_hists.at(k_reco_nu_e).at(i).at(j) = new TH1D ( Form("h_reco_nu_e_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 3);

            // Contained Fraction
            TH1D_hists.at(k_reco_contained_fraction).at(i).at(j) = new TH1D ( Form("h_reco_contained_fraction_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);

            // Run Number
            TH1D_hists.at(k_reco_run_number).at(i).at(j) = new TH1D ( Form("h_reco_run_number_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"; Run Number; Entries", 270, 4500, 18000);
        
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
            TH1D_hists_particle.at(k_reco_dEdx_cali_y_plane_par).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_cali_y_plane_par_%s_%s",_util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_y_par).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_par_%s_%s", _util.cut_dirs.at(i).c_str(), _util.particle_types.at(j).c_str()) ,"", 40, 0, 10);

        }
    }


    // -------------------------------------------------------------------------

    // Intialising true histograms in here
    if (_type == _util.k_mc || _type == _util.k_dirt){
        
        if (_type == _util.k_mc){
            // Initalise the histograms for the TEfficency
            TEfficiency_hists.resize(_util.k_cuts_MAX);

            for (unsigned int l = 0; l < _util.k_cuts_MAX; l++ ){
                TEfficiency_hists.at(l) = new TH1D( Form("h_true_nu_E_%s",_util.cut_dirs.at(l).c_str() ), "", 15, 0, 5 );
            }
        }

        // Initalise the True Nue
        TH1D_true_hists.resize(k_TH1D_true_MAX);
        TH2D_true_hists.resize(k_TH2D_true_MAX);

        TH1D_true_hists.at(k_true_nue_theta) = new TH1D( Form("h_true_nue_theta_%s", _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Theta [degrees]; Entries",           14, 0, 140 );
        TH1D_true_hists.at(k_true_nue_phi)   = new TH1D( Form("h_true_nue_phi_%s",   _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Phi [degrees]; Entries",             14, 0, 100 );
        TH1D_true_hists.at(k_true_nue_angle) = new TH1D( Form("h_true_nue_angle_%s", _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Angle from NuMI [degrees]; Entries", 18, 0, 180 );
        TH1D_true_hists.at(k_true_nue_px)    = new TH1D( Form("h_true_nue_px_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Px [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_py)    = new TH1D( Form("h_true_nue_py_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Py [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_pz)    = new TH1D( Form("h_true_nue_pz_%s",    _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} Pz [GeV/c]; Entries", 14, 0, 5);
        TH1D_true_hists.at(k_true_nue_e)     = new TH1D( Form("h_true_nue_e_%s",     _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} E [GeV]; Entries",    15, 0, 5);
        TH1D_true_hists.at(k_true_nue_p)     = new TH1D( Form("h_true_nue_p_%s",     _util.type_prefix.at(_type).c_str() ), ";True #nu_{e} P [GeV/c]; Entries",  14, 0, 5);
        TH1D_true_hists.at(k_true_vtx_x)     = new TH1D( Form("h_true_vtx_x_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx x [cm]; Entries", 20, -10, 270);
        TH1D_true_hists.at(k_true_vtx_y)     = new TH1D( Form("h_true_vtx_y_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx y [cm]; Entries", 20, -10, 120);
        TH1D_true_hists.at(k_true_vtx_z)     = new TH1D( Form("h_true_vtx_z_%s",     _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx z [cm]; Entries", 40, -10, 1050);
        TH1D_true_hists.at(k_true_vtx_x_sce) = new TH1D( Form("h_true_vtx_x_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx x Space Charge Corr. [cm]; Entries", 20, -10, 270);
        TH1D_true_hists.at(k_true_vtx_y_sce) = new TH1D( Form("h_true_vtx_y_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx y Space Charge Corr. [cm]; Entries", 20, -10, 120);
        TH1D_true_hists.at(k_true_vtx_z_sce) = new TH1D( Form("h_true_vtx_z_sce_%s", _util.type_prefix.at(_type).c_str() ) ,";True #nu_{e} Vtx z Space Charge Corr. [cm]; Entries", 40, -10, 1050);

        TH2D_true_hists.at(k_true_nue_phi_theta)    = new TH2D( Form("h_true_nue_phi_theta_%s", _util.type_prefix.at(_type).c_str()),   ";True #nu_{e} Phi [degrees];True #nu_{e} Theta [degrees]",     14, 0, 100, 16, 0, 140 );
        TH2D_true_hists.at(k_true_nue_energy_theta) = new TH2D( Form("h_true_nue_energy_theta_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} E [GeV];True #nu_{e} Theta [degrees]",           25, 0, 5, 16, 0, 140);
        TH2D_true_hists.at(k_true_nue_energy_phi)   = new TH2D( Form("h_true_nue_energy_phi_%s", _util.type_prefix.at(_type).c_str()),  ";True #nu_{e} E [GeV];True #nu_{e} Phi [degrees]",             25, 0, 5, 14, 0, 100);
        TH2D_true_hists.at(k_true_nue_energy_angle) = new TH2D( Form("h_true_nue_energy_angle_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} E [GeV];True #nu_{e} Angle from NuMI [degrees]", 25, 0, 5, 18, 0, 180);
    
        TH2D_true_hists.at(k_true_nue_vtx_z_y)     = new TH2D( Form("h_true_nue_vtx_z_y_%s", _util.type_prefix.at(_type).c_str()),    ";True #nu_{e} Vtx Z [cm] ;True #nu_{e} Vtx Y [cm]", 40, -10, 1050, 20, -10, 120);
        TH2D_true_hists.at(k_true_nue_vtx_z_y_sce) = new TH2D( Form("h_true_nue_vtx_z_y_sce_%s", _util.type_prefix.at(_type).c_str()),";True #nu_{e} Vtx Z  Space Charge Corr. [cm];True #nu_{e} Vtx Y Space Charge Corr. [cm]", 40, -10, 1050, 20, -10, 120);
    
        TH2D_true_hists.at(k_true_elec_E_reco_elec_E) = new TH2D( Form("h_true_elec_E_reco_elec_E_%s", _util.type_prefix.at(_type).c_str()),    ";True e^{-} Energy [GeV] ;Reco e^{-} Energy [GeV]", 25, 0, 4, 25, 0, 4);
        TH2D_true_hists.at(k_true_nu_vtx_x_reco_nu_vtx_x) = new TH2D( Form("h_true_nu_vtx_x_reco_nu_vtx_x_%s", _util.type_prefix.at(_type).c_str()),    ";True #nu_{e} Vtx X [cm] ;Reco #nu_{e} Vtx X [cm]", 20, -10, 270, 20, -10, 270);
        TH2D_true_hists.at(k_true_nu_vtx_y_reco_nu_vtx_y) = new TH2D( Form("h_true_nu_vtx_y_reco_nu_vtx_y_%s", _util.type_prefix.at(_type).c_str()),    ";True #nu_{e} Vtx Y [cm] ;Reco #nu_{e} Vtx Y [cm]", 20, -10, 120, 20, -10, 120);
        TH2D_true_hists.at(k_true_nu_vtx_z_reco_nu_vtx_z) = new TH2D( Form("h_true_nu_vtx_z_reco_nu_vtx_z_%s", _util.type_prefix.at(_type).c_str()),    ";True #nu_{e} Vtx Z [cm] ;Reco #nu_{e} Vtx Z [cm]", 40, -10, 1050, 40, -10, 1050);
    }

    // -------------------------------------------------------------------------

    // Flash histograms
    TH1D_flash_hists.resize(k_TH1D_flash_MAX);
    TH1D_flash_hists.at(k_flash_time) = new TH1D( Form("h_flash_time_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe)   = new TH1D( Form("h_flash_pe_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity [PE]; Entries", 40, 0, 10000 );

    TH1D_flash_hists.at(k_flash_time_sid1) = new TH1D( Form("h_flash_time_sid1_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Neutrino Candiate [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid1)   = new TH1D( Form("h_flash_pe_sid1_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity (Neutrino Candiate) [PE]; Entries", 40, 0, 10000 );

    TH1D_flash_hists.at(k_flash_time_sid0) = new TH1D( Form("h_flash_time_sid0_%s", _util.type_prefix.at(_type).c_str()), "; Flash Time Non Neutrino Canidate [us]; Entries", 100, 0, 25 );
    TH1D_flash_hists.at(k_flash_pe_sid0)   = new TH1D( Form("h_flash_pe_sid0_%s", _util.type_prefix.at(_type).c_str()),   "; Largest Flash Intensity (Non Neutrino Canidate) [PE]; Entries", 40, 0, 10000 );
    
    // Interaction Histograms
    TH1D_interaction_hists.resize(_util.k_interactions_MAX);
    for (unsigned int p =0 ; p < TH1D_interaction_hists.size(); p++ ){
        TH1D_interaction_hists.at(p) = new TH1D( Form("h_true_nue_E_%s", _util.interaction_types.at(p).c_str() ), "; True #nu_{e} Energy; Entries", 15, 0, 5 );
    }

    // 2D signal and background separation plots
    
    // Resizing the histograms
    TH2D_hists.resize(k_TH2D_reco_MAX);
    
    // Loop over the histograms to resize
    for (unsigned int i = 0; i < TH2D_hists.size(); i++){
        TH2D_hists.at(i).resize(_util.k_sig_bkg_MAX);
    }

    // Loop over the histograms
    for (unsigned int i = 0; i < TH2D_hists.size(); i++){
       
        // Loop over the types 
        for (unsigned int k =0 ; k< TH2D_hists.at(i).size(); k++){
            // dEdx vs shower vertex distance 
            TH2D_hists.at(i).at(k)     = new TH2D( Form("h_reco_shr_dEdx_shr_dist_%s", _util.sig_bkg_prefix.at(k).c_str()),";Collection Plane dEdx (calibrated) [MeV/cm];Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
       
        }
        
    }

}
// -----------------------------------------------------------------------------
void histogram_helper::FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight){

    // Calculate some variables
    double reco_shr_p = std::sqrt(SC.shr_px*SC.shr_px + SC.shr_py*SC.shr_py + SC.shr_pz*SC.shr_pz);

    // Reconstructed Energy of neutrino
    double INTERCEPT = 0.0;
    double SLOPE = 0.83;
    double reco_nu_e = (SC.shr_energy_tot_cali + INTERCEPT) / SLOPE + SC.trk_energy_tot;

    // Now fill the histograms!
    TH1D_hists.at(k_reco_vtx_x).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x, weight);
    TH1D_hists.at(k_reco_vtx_y).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y, weight);
    TH1D_hists.at(k_reco_vtx_z).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z, weight);

    TH1D_hists.at(k_reco_vtx_x_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_x, weight);
    TH1D_hists.at(k_reco_vtx_y_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_y, weight);
    TH1D_hists.at(k_reco_vtx_z_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_z, weight);

    TH1D_hists.at(k_reco_dEdx_cali_u_plane).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_U_cali, weight); // Just the collection plane!
    TH1D_hists.at(k_reco_dEdx_cali_v_plane).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_V_cali, weight); // Just the collection plane!
    TH1D_hists.at(k_reco_dEdx_cali_y_plane).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_Y_cali, weight); // Just the collection plane!
    
    TH1D_hists.at(k_reco_leading_mom).at(cut_index).at(classification_index)->Fill(reco_shr_p, weight);
    
    TH1D_hists.at(k_reco_shower_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.shr_distance, weight);

    TH1D_hists.at(k_reco_track_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.trk_distance, weight);
    
    TH1D_hists.at(k_reco_shr_hits_max).at(cut_index).at(classification_index)->Fill(SC.shr_hits_max, weight);
    
    TH1D_hists.at(k_reco_n_track_contained).at(cut_index).at(classification_index)->Fill(SC.n_tracks_contained, weight);
    
    TH1D_hists.at(k_reco_n_shower_contained).at(cut_index).at(classification_index)->Fill(SC.n_showers_contained, weight);
    
    TH1D_hists.at(k_reco_leading_shower_phi).at(cut_index).at(classification_index)->Fill(SC.shr_phi * 180/3.14159, weight);
    
    TH1D_hists.at(k_reco_leading_shower_theta).at(cut_index).at(classification_index)->Fill(SC.shr_theta * 180/3.14159, weight);
    
    TH1D_hists.at(k_reco_leading_shower_cos_theta).at(cut_index).at(classification_index)->Fill(std::cos(SC.shr_theta), weight);
    
    TH1D_hists.at(k_reco_shower_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_showers, weight);
    
    TH1D_hists.at(k_reco_track_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_tracks, weight);

    TH1D_hists.at(k_reco_topological_score).at(cut_index).at(classification_index)->Fill(SC.topological_score, weight);

    TH1D_hists.at(k_reco_track_shower_dist).at(cut_index).at(classification_index)->Fill(SC.tksh_distance, weight);

    if (SC.n_tracks > 0) TH1D_hists.at(k_reco_track_shower_angle).at(cut_index).at(classification_index)->Fill(SC.tksh_angle*180/3.14159, weight);

    TH1D_hists.at(k_reco_hits_ratio).at(cut_index).at(classification_index)->Fill(SC.hits_ratio, weight);

    TH1D_hists.at(k_reco_shower_score).at(cut_index).at(classification_index)->Fill(SC.shr_score, weight);

    TH1D_hists.at(k_reco_track_score).at(cut_index).at(classification_index)->Fill(SC.trk_score, weight);

    TH1D_hists.at(k_reco_shower_energy_tot_cali).at(cut_index).at(classification_index)->Fill(SC.shr_energy_tot_cali, weight);

    TH1D_hists.at(k_reco_shr_hits_tot).at(cut_index).at(classification_index)->Fill(SC.shr_hits_tot, weight);

    TH1D_hists.at(k_reco_shr_hits_y_tot).at(cut_index).at(classification_index)->Fill(SC.shr_hits_y_tot, weight);

    TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_u).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_2cm_dedx_U, weight);
    TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_v).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_2cm_dedx_V, weight);    
    TH1D_hists.at(k_reco_shr_trkfit_2cm_dEdx_y).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_2cm_dedx_Y, weight);
    
    TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_u).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap05_dedx_U, weight);
    TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_v).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap05_dedx_V, weight);
    TH1D_hists.at(k_reco_shr_trkfit_gap05_dEdx_y).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap05_dedx_Y, weight);
    
    
    TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_u).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap10_dedx_U, weight);
    TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_v).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap10_dedx_V, weight);
    TH1D_hists.at(k_reco_shr_trkfit_gap10_dEdx_y).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_gap10_dedx_Y, weight);


    TH1D_hists.at(k_reco_opfilter_beam).at(cut_index).at(classification_index)->Fill(SC.opfilter_pe_beam, weight);
    TH1D_hists.at(k_reco_opfilter_veto).at(cut_index).at(classification_index)->Fill(SC.opfilter_pe_veto, weight);
    
    TH1D_hists.at(k_reco_softwaretrig).at(cut_index).at(classification_index)->Fill(SC.swtrig, weight);
    
    TH1D_hists.at(k_reco_nslice).at(cut_index).at(classification_index)->Fill(SC.nslice, weight);

    TH1D_hists.at(k_reco_slclustfrac).at(cut_index).at(classification_index)->Fill(SC.slclustfrac, weight);

    TH1D_hists.at(k_reco_cosmicIP).at(cut_index).at(classification_index)->Fill(SC.CosmicIP, weight);
    TH1D_hists.at(k_reco_CosmicIPAll3D).at(cut_index).at(classification_index)->Fill(SC.CosmicIPAll3D, weight);
    TH1D_hists.at(k_reco_CosmicDirAll3D).at(cut_index).at(classification_index)->Fill(SC.CosmicDirAll3D, weight);

    TH1D_hists.at(k_reco_shr_tkfit_dedx_u).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_U, weight);
    TH1D_hists.at(k_reco_shr_tkfit_dedx_v).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_V, weight);
    TH1D_hists.at(k_reco_shr_tkfit_dedx_y).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);

    // Split the dedx into good and bad angles
    if (SC.shr_theta * 180/3.14159 < 80 || SC.shr_theta * 180/3.14159 > 100) {
        TH1D_hists.at(k_reco_shr_tkfit_dedx_y_good_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);
    }
    else {
        TH1D_hists.at(k_reco_shr_tkfit_dedx_y_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);
        TH1D_hists.at(k_reco_shr_tkfit_dedx_v_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_V, weight);
        TH1D_hists.at(k_reco_shr_tkfit_dedx_u_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_U, weight);
    }
     

    TH1D_hists.at(k_reco_shrsubclusters).at(cut_index).at(classification_index)->Fill(SC.shrsubclusters0 + SC.shrsubclusters1 + SC.shrsubclusters2, weight);

    TH1D_hists.at(k_reco_shrmoliereavg).at(cut_index).at(classification_index)->Fill(SC.shrmoliereavg, weight);
    TH1D_hists.at(k_reco_shrmoliererms).at(cut_index).at(classification_index)->Fill(SC.shrmoliererms, weight);

    TH1D_hists.at(k_reco_CylFrac2h_1cm).at(cut_index).at(classification_index)->Fill(SC.CylFrac2h_1cm, weight);
    
    TH1D_hists.at(k_reco_DeltaRMS2h).at(cut_index).at(classification_index)->Fill(SC.DeltaRMS2h, weight);
    
    TH1D_hists.at(k_reco_shrPCA1CMed_5cm).at(cut_index).at(classification_index)->Fill(SC.shrPCA1CMed_5cm, weight);
    
    TH1D_hists.at(k_reco_shrMCSMom).at(cut_index).at(classification_index)->Fill(SC.shrMCSMom, weight);

    TH1D_hists.at(k_reco_closestNuCosmicDist).at(cut_index).at(classification_index)->Fill(SC._closestNuCosmicDist, weight);

    if (SC.n_tracks > 0) TH1D_hists.at(k_reco_trk_len).at(cut_index).at(classification_index)->Fill(SC.trk_len, weight);

    TH1D_hists.at(k_reco_nu_e).at(cut_index).at(classification_index)->Fill(reco_nu_e, weight);

    TH1D_hists.at(k_reco_contained_fraction).at(cut_index).at(classification_index)->Fill(SC.contained_fraction, weight);

    TH1D_hists.at(k_reco_run_number).at(cut_index).at(classification_index)->Fill(SC.run, weight);
    
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // Partice Type Hists
    TH1D_hists_particle.at(k_reco_dEdx_cali_y_plane_par).at(cut_index).at(_par_type)->Fill(SC.shr_dedx_Y_cali, weight);
    TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_y_par).at(cut_index).at(_par_type)->Fill(SC.shr_tkfit_dedx_Y, weight);

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // Flash histograms
    if (type == _util.k_mc || type == _util.k_dirt){
        // We also want to apply the software trigger to the MC, just to be fair
        TH1D_hists.at(k_reco_flash_time).at(cut_index).at(classification_index)->Fill(SC.flash_time + 0.055 -0.359, weight);
    }
    if (type == _util.k_ext){
        TH1D_hists.at(k_reco_flash_time).at(cut_index).at(classification_index)->Fill(SC.flash_time -0.359, weight);
    }
    if (type == _util.k_data){
        TH1D_hists.at(k_reco_flash_time).at(cut_index).at(classification_index)->Fill(SC.flash_time, weight);
    }
    
    // Also apply the SW Trig
    
    if (type == _util.k_mc || type == _util.k_dirt){
        TH1D_hists.at(k_reco_flash_pe).at(cut_index).at(classification_index)->Fill(SC.flash_pe, weight);
    
    }
    else TH1D_hists.at(k_reco_flash_pe).at(cut_index).at(classification_index)->Fill(SC.flash_pe, weight);

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

        if (classification_index == _util.k_nue_cc || classification_index == _util.k_nue_cc_mixed){
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

            TH2D_true_hists.at(k_true_nue_phi_theta)   ->Fill(nu_phi, nu_theta, weight);
            TH2D_true_hists.at(k_true_nue_energy_theta)->Fill(SC.nu_e, nu_theta, weight);
            TH2D_true_hists.at(k_true_nue_energy_phi)  ->Fill(SC.nu_e, nu_phi, weight);
            TH2D_true_hists.at(k_true_nue_energy_angle)->Fill(SC.nu_e, nu_angle, weight);
            TH2D_true_hists.at(k_true_nue_vtx_z_y)       ->Fill(SC.true_nu_vtx_z,  SC.true_nu_vtx_y, weight);
            TH2D_true_hists.at(k_true_nue_vtx_z_y_sce)   ->Fill(SC.true_nu_vtx_sce_z,  SC.true_nu_vtx_sce_y, weight);

            if (_type == _util.k_mc){ 

                // True vs reco histograms
                if (SC.nslice == 1) TH2D_true_hists.at(k_true_elec_E_reco_elec_E)->Fill(SC.elec_e,  SC.shr_energy_tot_cali, weight);
                TH2D_true_hists.at(k_true_nu_vtx_x_reco_nu_vtx_x)->Fill(SC.true_nu_vtx_sce_x,  SC.reco_nu_vtx_sce_x, weight);
                TH2D_true_hists.at(k_true_nu_vtx_y_reco_nu_vtx_y)->Fill(SC.true_nu_vtx_sce_y,  SC.reco_nu_vtx_sce_y, weight);
                TH2D_true_hists.at(k_true_nu_vtx_z_reco_nu_vtx_z)->Fill(SC.true_nu_vtx_sce_z,  SC.reco_nu_vtx_sce_z, weight);


            
                // True nue interaction histograms
                if (interaction == "nue_cc_qe" || interaction == "nue_bar_cc_qe"){
                    TH1D_interaction_hists.at(_util.k_plot_qe)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_res" || interaction == "nue_bar_cc_res"){
                    TH1D_interaction_hists.at(_util.k_plot_res)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_dis" || interaction == "nue_bar_cc_dis"){
                    TH1D_interaction_hists.at(_util.k_plot_dis)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_coh" || interaction == "nue_bar_cc_coh"){
                    TH1D_interaction_hists.at(_util.k_plot_coh)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_mec" || interaction == "nue_bar_cc_mec"){
                    TH1D_interaction_hists.at(_util.k_plot_mec)->Fill(SC.nu_e, weight);
                }
                // NC
                else {
                    TH1D_interaction_hists.at(_util.k_plot_nc)->Fill(SC.nu_e, weight);
                }
            }

        }
       
    }

    // -----------------------------------------------------------------------------

    // Only do this for after the software trigger
    if (cut_index == _util.k_swtrig){

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

    // Fill 2D histograms for signal background rejection -- only for MC, dirt and EXT
    if (_type != _util.k_data){ 
        
        // For the dEdx vs shower distance we make the histogram just before the cut is applied
        if (cut_index == _util.k_shr_distance - 1 ){
            
            // This is the signal
            if (classification_index == _util.k_nue_cc || classification_index == _util.k_nue_cc_mixed){
                TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(_util.k_signal)->Fill(SC.shr_dedx_Y_cali, SC.shr_distance, weight);
            }
            // This is the background
            else {
                TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(_util.k_background)->Fill(SC.shr_dedx_Y_cali, SC.shr_distance, weight);
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

                // Try to clear some memory
                // delete TH1D_hists.at(u).at(i).at(j);

                if (break_early) break;
            }

        }
    }
}
// -----------------------------------------------------------------------------
void histogram_helper::WriteRecoPar(int type){

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
void histogram_helper::FillTEfficiency(int cut_index, std::string classification, SliceContainer SC, double weight){

    // Fill the histogram at the specified cut
    if (classification == "nue_cc" || classification == "nue_cc_mixed") TEfficiency_hists.at(cut_index)->Fill(SC.nu_e, weight);
}
// -----------------------------------------------------------------------------
void histogram_helper::WriteTEfficiency(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"TEff");
    if (bool_dir) dir->cd();
    
    for (unsigned int p = 0; p < TEfficiency_hists.size(); p++){

        // TEfficiency * teff = new TEfficiency(*TEfficiency_hists.at(p), *TEfficiency_hists.at(_util.k_unselected));
        // teff->Write( Form("h_true_nu_E_%s",_util.cut_dirs.at(p).c_str()) , TObject::kOverwrite);
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
void histogram_helper::WriteInteractions(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"Interaction");
    if (bool_dir) dir->cd();
    
    // TH1D
    for (unsigned int p = 0; p < TH1D_interaction_hists.size(); p++){
        TH1D_interaction_hists.at(p)->SetOption("hist,E");
        TH1D_interaction_hists.at(p)->Write("",TObject::kOverwrite);
    }
    
    
}
// -----------------------------------------------------------------------------
void histogram_helper::Write_2DSigBkgHists(){
    
    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(f_nuexsec, dir ,"2D");
    if (bool_dir) dir->cd();
    
    // TH2D
    for (unsigned int p = 0; p < TH2D_hists.size(); p++){
        
        TH2D_hists.at(p).at(_util.k_signal)->SetOption("box");
        TH2D_hists.at(p).at(_util.k_signal)->SetFillColor(30);
        TH2D_hists.at(p).at(_util.k_signal)->Write("",TObject::kOverwrite);

        TH2D_hists.at(p).at(_util.k_background)->SetOption("box");
        TH2D_hists.at(p).at(_util.k_background)->SetFillColor(46);
        TH2D_hists.at(p).at(_util.k_background)->Write("",TObject::kOverwrite);
    }
    
}
// -----------------------------------------------------------------------------
