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
void HistogramHelper::Initialise(int type, const char* run_period, const char * file_out, int weight_cfg, Utility util ){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    _util = util;

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
        
            // Reco Vtx X, Y, Z
            TH1D_hists.at(k_reco_vtx_x).at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 15, -10, 270);

            TH1D_hists.at(k_reco_vtx_y).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, -120, 120);
            
            TH1D_hists.at(k_reco_vtx_z).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, -10, 1050);

            TH1D_hists.at(k_reco_vtx_x_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 15, -10, 270);

            TH1D_hists.at(k_reco_vtx_y_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, -120, 120);
            
            TH1D_hists.at(k_reco_vtx_z_sce).at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_sce_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, -10, 1050);
            
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
            
            // YZ and X Calibrated energy of all the showers
            TH1D_hists.at(k_reco_shower_energy_tot_cali).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_tot_cali_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 4);
            
            // Calibrated energy of just the leading shower
            TH1D_hists.at(k_reco_shower_energy_cali).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_cali_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 4);

            // Set the bins for the reco energy
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
            TH1D_hists.at(k_reco_shower_energy_tot_cali_rebin).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_tot_cali_rebin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", _util.reco_shr_bins.size()-1, edges);
            TH1D_hists.at(k_reco_shower_energy_cali_rebin).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_cali_rebin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", _util.reco_shr_bins.size()-1, edges);

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
            TH1D_hists.at(k_reco_cosmicIP).at(i).at(j)       = new TH1D ( Form("h_reco_cosmicIP_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 250);
            TH1D_hists.at(k_reco_CosmicIPAll3D).at(i).at(j)  = new TH1D ( Form("h_reco_CosmicIPAll3D_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 200);
            TH1D_hists.at(k_reco_CosmicDirAll3D).at(i).at(j) = new TH1D ( Form("h_reco_CosmicDirAll3D_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -1, 1);

            // dEdx with the trackfit variable
            TH1D_hists.at(k_reco_shr_tkfit_dedx_u).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_u_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_v).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_v_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_max).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_max_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_max_with_tracks).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_max_with_tracks_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y_no_tracks).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_no_tracks_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_max_no_tracks).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_max_no_tracks_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            TH1D_hists.at(k_reco_shr_tkfit_dedx_y_good_theta).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_good_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_y_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_y_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_v_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_v_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            TH1D_hists.at(k_reco_shr_tkfit_dedx_u_bad_theta).at(i).at(j)  = new TH1D ( Form("h_reco_shr_tkfit_dedx_u_bad_theta_%s_%s",_util.cut_dirs.at(i).c_str(),  _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);
            
            
            // Flash time
            TH1D_hists.at(k_reco_flash_time).at(i).at(j) = new TH1D ( Form("h_reco_flash_time_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 25 );

            // Flash PE
            TH1D_hists.at(k_reco_flash_pe).at(i).at(j) = new TH1D ( Form("h_reco_flash_pe_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 5000 );

            // Shower Subcluster in each plane
            TH1D_hists.at(k_reco_shrsubclusters).at(i).at(j)  = new TH1D ( Form("h_reco_shrsubclusters_%s_%s", _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // Shower Moliers Average
            TH1D_hists.at(k_reco_shrmoliereavg).at(i).at(j) = new TH1D ( Form("h_reco_shrmoliereavg_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 30);

            // Shower Moliere RMS
            TH1D_hists.at(k_reco_shrmoliererms).at(i).at(j) = new TH1D ( Form("h_reco_shrmoliererms_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 30000 );

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

            // Neutrino Purity
            TH1D_hists.at(k_reco_nu_purity_from_pfp).at(i).at(j) = new TH1D ( Form("h_reco_nu_purity_from_pfp_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,";", 250, 0, 1.1);
        
            // CRT veto
            TH1D_hists.at(k_reco_crtveto).at(i).at(j) = new TH1D ( Form("h_reco_crtveto_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 2, 0, 2);

            // CRT hit PE
            TH1D_hists.at(k_reco_crthitpe).at(i).at(j) = new TH1D ( Form("h_reco_crthitpe_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 50, 0, 500);

            // Shower angle with respect to the numi target direction
            TH1D_hists.at(k_reco_shr_ang_numi).at(i).at(j) = new TH1D ( Form("h_reco_shr_ang_numi_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 14, -190, 190);

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
                TEfficiency_hists.at(k_eff_nu_E).at(l)   = new TH1D( Form("h_true_nu_E_%s",_util.cut_dirs.at(l).c_str() ), "", 15, 0, 5 );
                TEfficiency_hists.at(k_eff_elec_E).at(l) = new TH1D( Form("h_true_elec_E_%s",_util.cut_dirs.at(l).c_str() ), "", 15, 0, 5 );
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
            
            TH1D_true_hists.at(i).at(k_true_nue_theta) = new TH1D( Form("h_true_nue_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} BNB #theta [degrees]; Entries",           14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_nue_phi)   = new TH1D( Form("h_true_nue_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} BNB #phi [degrees]; Entries",              14, 0, 40);

            TH1D_true_hists.at(i).at(k_true_nue_theta_numi) = new TH1D( Form("h_true_nue_theta_numi_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} NuMI #theta [degrees]; Entries",           14, 0, 140 );
            TH1D_true_hists.at(i).at(k_true_nue_phi_numi)   = new TH1D( Form("h_true_nue_phi_numi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} NuMI #phi [degrees]; Entries",             25, -180, 180 );
            
            TH1D_true_hists.at(i).at(k_true_nue_angle) = new TH1D( Form("h_true_nue_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} Angle from NuMI [degrees]; Entries", 120, 0, 120 );
            TH1D_true_hists.at(i).at(k_true_nue_px)    = new TH1D( Form("h_true_nue_px_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} Px [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_py)    = new TH1D( Form("h_true_nue_py_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} Py [GeV/c]; Entries", 14, 0, 1);
            TH1D_true_hists.at(i).at(k_true_nue_pz)    = new TH1D( Form("h_true_nue_pz_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} Pz [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_e)     = new TH1D( Form("h_true_nue_e_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} E [GeV]; Entries",    15, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_p)     = new TH1D( Form("h_true_nue_p_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} P [GeV/c]; Entries",  14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_vtx_x)     = new TH1D( Form("h_true_vtx_x_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx x [cm]; Entries", 20, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y)     = new TH1D( Form("h_true_vtx_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx y [cm]; Entries", 20, -10, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z)     = new TH1D( Form("h_true_vtx_z_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx z [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_vtx_x_sce) = new TH1D( Form("h_true_vtx_x_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx x Space Charge Corr. [cm]; Entries", 20, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y_sce) = new TH1D( Form("h_true_vtx_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx y Space Charge Corr. [cm]; Entries", 20, -10, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z_sce) = new TH1D( Form("h_true_vtx_z_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} Vtx z Space Charge Corr. [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_nu_ang_targ)     = new TH1D( Form("h_true_nu_ang_targ_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} Angle from NuMI Target [degrees]; Entries",  40, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_ang_targ)   = new TH1D( Form("h_true_elec_ang_targ_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True electron Angle from NuMI Target [degrees]; Entries", 25, 0, 180 );

            TH1D_true_hists.at(i).at(k_true_elec_E)     = new TH1D( Form("h_true_elec_E_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True electron energy [GeV]; Entries",                15, 0, 5 );
            TH1D_true_hists.at(i).at(k_true_elec_theta) = new TH1D( Form("h_true_elec_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True electron BNB #theta [degrees]; Entries",            14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_phi)   = new TH1D( Form("h_true_elec_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True electron BNB #phi [degrees]; Entries",              14, 0, 40);


            TH2D_true_hists.at(i).at(k_true_nue_phi_theta)    = new TH2D( Form("h_true_nue_phi_theta_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} Phi [degrees];True #nu_{e} Theta [degrees]",     14, 0, 80, 16, 20, 140 );
            TH2D_true_hists.at(i).at(k_true_nue_energy_theta) = new TH2D( Form("h_true_nue_energy_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),";True #nu_{e} E [GeV];True #nu_{e} Theta [degrees]",           25, 0, 5, 16, 0, 140);
            TH2D_true_hists.at(i).at(k_true_nue_energy_phi)   = new TH2D( Form("h_true_nue_energy_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),  ";True #nu_{e} E [GeV];True #nu_{e} Phi [degrees]",             25, 0, 5, 14, 0, 75);
            TH2D_true_hists.at(i).at(k_true_nue_energy_angle) = new TH2D( Form("h_true_nue_energy_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),";True #nu_{e} E [GeV];True #nu_{e} Angle from NuMI [degrees]", 25, 0, 5, 30, 0, 120);
        
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y)     = new TH2D( Form("h_true_nue_vtx_z_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Vtx Z [cm] ;True #nu_{e} Vtx Y [cm]", 40, -10, 1050, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y_sce) = new TH2D( Form("h_true_nue_vtx_z_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),";True #nu_{e} Vtx Z  Space Charge Corr. [cm];True #nu_{e} Vtx Y Space Charge Corr. [cm]", 40, -10, 1050, 20, -10, 120);
        
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E) = new TH2D( Form("h_true_elec_E_reco_elec_E_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e^{-} Energy [GeV] ;Reco Leading Shower Energy [GeV]", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins.size()-1, edges);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E)     = new TH2D( Form("h_true_nu_E_reco_nu_E_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Energy [GeV] ;Reco #nu_{e} Energy [GeV]", 25, 0, 4, 25, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e^{-} Energy [GeV] ;Reco Leading Shower Energy [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E_extra_bins)     = new TH2D( Form("h_true_nu_E_reco_nu_E_extra_bins_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Energy [GeV] ;Reco #nu_{e} Energy [GeV]", 50, 0, 4, 50, 0, 4);

            TH2D_true_hists.at(i).at(k_true_nu_vtx_x_reco_nu_vtx_x) = new TH2D( Form("h_true_nu_vtx_x_reco_nu_vtx_x_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Vtx X [cm] ;Reco #nu_{e} Vtx X [cm]", 20, -10, 270, 20, -10, 270);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_y_reco_nu_vtx_y) = new TH2D( Form("h_true_nu_vtx_y_reco_nu_vtx_y_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Vtx Y [cm] ;Reco #nu_{e} Vtx Y [cm]", 20, -10, 120, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_z_reco_nu_vtx_z) = new TH2D( Form("h_true_nu_vtx_z_reco_nu_vtx_z_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} Vtx Z [cm] ;Reco #nu_{e} Vtx Z [cm]", 40, -10, 1050, 40, -10, 1050);
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
    TH1D_interaction_hists.resize(2); // Unselected and Selected

    for (unsigned int i = 0; i < TH1D_true_hists.size(); i++){
        TH1D_interaction_hists.at(i).resize(_util.k_interactions_MAX);
    }
    
    for (unsigned int p =0 ; p < TH1D_interaction_hists.at(0).size(); p++ ){
        TH1D_interaction_hists.at(0).at(p) = new TH1D( Form("h_true_nue_E_%s_unselected", _util.interaction_types.at(p).c_str() ), "; True #nu_{e} Energy; Entries", 15, 0, 5 );
        TH1D_interaction_hists.at(1).at(p) = new TH1D( Form("h_true_nue_E_%s_selected",   _util.interaction_types.at(p).c_str() ), "; True #nu_{e} Energy; Entries", 15, 0, 5 );
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
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(k)              = new TH2D( Form("h_reco_shr_dEdx_shr_dist_%s",            _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist_post).at(k)         = new TH2D( Form("h_reco_shr_dEdx_shr_dist_post_%s",       _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist).at(k)          = new TH2D( Form("h_reco_shr_dEdx_max_shr_dist_%s",        _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (All Planes) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist_post).at(k)     = new TH2D( Form("h_reco_shr_dEdx_max_shr_dist_post_%s",   _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (All Planes) [MeV/cm];Leading Shower to Vertex Distance [cm]", 40, 0, 10, 15, 0, 20);
        TH2D_hists.at(k_reco_shr_dEdx_shr_dist_large_dedx).at(k)   = new TH2D( Form("h_reco_shr_dEdx_shr_dist_large_dedx_%s", _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower to Vertex Distance [cm]", 80, 10, 30, 15, 0, 20);
    
        // dEdx vs Moliere Average
        TH2D_hists.at(k_reco_shr_dEdx_moliere).at(k)          = new TH2D( Form("h_reco_shr_dEdx_moliere_%s", _util.sig_bkg_prefix.at(k).c_str()),";Leading Shower dEdx (Collection Plane) [MeV/cm];Leading Shower Moliere Avg [deg]", 40, 0, 10, 20, 0, 30);
    
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
    

}
// -----------------------------------------------------------------------------
void HistogramHelper::FillHists(int type, int classification_index, std::string interaction, int _par_type, int cut_index, SliceContainer SC, double weight){

    // Calculate some variables
    double reco_shr_p = std::sqrt(SC.shr_px*SC.shr_px + SC.shr_py*SC.shr_py + SC.shr_pz*SC.shr_pz);

    // Reconstructed Energy of neutrino
    double reco_nu_e = SC.shr_energy_tot_cali / 0.83 + SC.trk_energy_tot;

    // Check if the interaction was in the FV
    bool true_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z);

    // The angle of the reconstructed shower relative to the NuMI target to detector direction
    double shr_ang_numi = _util.GetNuMIAngle(SC.shr_px, SC.shr_py, SC.shr_pz, "target");

    // Get the trkfit dedx max variable
    double dedx_max = SC.GetdEdxMax();

    // For filling histograms, we lump the cosmic nue and cosmic nuebar into the generic cosmic category (so the plots can be normalised properly)
    if (classification_index == _util.k_cosmic_nue || classification_index == _util.k_cosmic_nuebar) classification_index = _util.k_cosmic;
    
    // Now fill the histograms!
    TH1D_hists.at(k_reco_vtx_x).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x, weight);
    TH1D_hists.at(k_reco_vtx_y).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y, weight);
    TH1D_hists.at(k_reco_vtx_z).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z, weight);

    TH1D_hists.at(k_reco_vtx_x_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_x, weight);
    TH1D_hists.at(k_reco_vtx_y_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_y, weight);
    TH1D_hists.at(k_reco_vtx_z_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_z, weight);
 
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

    TH1D_hists.at(k_reco_shower_energy_tot_cali).at(cut_index).at(classification_index)->Fill(SC.shr_energy_tot_cali/0.83, weight);
    TH1D_hists.at(k_reco_shower_energy_tot_cali_rebin).at(cut_index).at(classification_index)->Fill(SC.shr_energy_tot_cali/0.83, weight);
    TH1D_hists.at(k_reco_shower_energy_cali).at(cut_index).at(classification_index)->Fill(SC.shr_energy_cali/0.83, weight);
    TH1D_hists.at(k_reco_shower_energy_cali_rebin).at(cut_index).at(classification_index)->Fill(SC.shr_energy_cali/0.83, weight);

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
    TH1D_hists.at(k_reco_shr_tkfit_dedx_max).at(cut_index).at(classification_index)->Fill(dedx_max, weight);
    if (SC.n_tracks > 0) TH1D_hists.at(k_reco_shr_tkfit_dedx_max_with_tracks).at(cut_index).at(classification_index)->Fill(dedx_max, weight);
    
    if (SC.n_tracks == 0) TH1D_hists.at(k_reco_shr_tkfit_dedx_y_no_tracks).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);
    if (SC.n_tracks == 0) TH1D_hists.at(k_reco_shr_tkfit_dedx_max_no_tracks).at(cut_index).at(classification_index)->Fill(dedx_max,  weight);

    // Split the dedx into good and bad angles
    
    if (SC.shr_theta * 180/3.14159 > 80 && SC.shr_theta * 180/3.14159 < 100) {
        TH1D_hists.at(k_reco_shr_tkfit_dedx_y_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);
        TH1D_hists.at(k_reco_shr_tkfit_dedx_v_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_V, weight);
        TH1D_hists.at(k_reco_shr_tkfit_dedx_u_bad_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_U, weight);
    }
    else {
        TH1D_hists.at(k_reco_shr_tkfit_dedx_y_good_theta).at(cut_index).at(classification_index)->Fill(SC.shr_tkfit_dedx_Y, weight);
    }
     

    TH1D_hists.at(k_reco_shrsubclusters).at(cut_index).at(classification_index)->Fill(SC.shrsubclusters0 + SC.shrsubclusters1 + SC.shrsubclusters2, weight);

    TH1D_hists.at(k_reco_shrmoliereavg).at(cut_index).at(classification_index)->Fill(SC.shrmoliereavg, weight);
    TH1D_hists.at(k_reco_shrmoliererms).at(cut_index).at(classification_index)->Fill(SC.shrmoliererms, weight);
    
    TH1D_hists.at(k_reco_shrPCA1CMed_5cm).at(cut_index).at(classification_index)->Fill(SC.shrPCA1CMed_5cm, weight);
    
    TH1D_hists.at(k_reco_shrMCSMom).at(cut_index).at(classification_index)->Fill(SC.shrMCSMom, weight);

    TH1D_hists.at(k_reco_closestNuCosmicDist).at(cut_index).at(classification_index)->Fill(SC._closestNuCosmicDist, weight);

    if (SC.n_tracks > 0) TH1D_hists.at(k_reco_trk_len).at(cut_index).at(classification_index)->Fill(SC.trk_len, weight);

    TH1D_hists.at(k_reco_nu_e).at(cut_index).at(classification_index)->Fill(reco_nu_e, weight);

    TH1D_hists.at(k_reco_contained_fraction).at(cut_index).at(classification_index)->Fill(SC.contained_fraction, weight);

    TH1D_hists.at(k_reco_run_number).at(cut_index).at(classification_index)->Fill(SC.run, weight);

    TH1D_hists.at(k_reco_nu_purity_from_pfp).at(cut_index).at(classification_index)->Fill(SC.nu_purity_from_pfp, weight);

    TH1D_hists.at(k_reco_crtveto).at(cut_index).at(classification_index)->Fill(SC.crtveto, weight);

    TH1D_hists.at(k_reco_crthitpe).at(cut_index).at(classification_index)->Fill(SC.crthitpe, weight);

    TH1D_hists.at(k_reco_shr_ang_numi).at(cut_index).at(classification_index)->Fill(shr_ang_numi, weight);
    
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // Partice Type Hists
    TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_max_par).at(cut_index).at(_par_type)->Fill(dedx_max, weight);
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

    

    // Only do this for MC
    if ( (_type == _util.k_mc) && (cut_index == _util.k_unselected || cut_index == _util.k_cuts_MAX-1)){


        int index = 0;
        if (cut_index == _util.k_cuts_MAX-1) index = 1; // If its the last cut index then we want to fill the second lot of histograms

        // Momentum of neutrino
        double p = std::sqrt(SC.true_nu_px*SC.true_nu_px + SC.true_nu_py*SC.true_nu_py + SC.true_nu_pz*SC.true_nu_pz);
        
        // Momentum of electron
        double p_elec = std::sqrt(SC.elec_px*SC.elec_px + SC.elec_py*SC.elec_py + SC.elec_pz*SC.elec_pz);
        
        // True nue theta in BNB coordinates (up from beam dir)
        double nu_theta = acos(SC.true_nu_pz/p) * 180 / 3.1415;
        // std::cout << SC.true_nu_px<< "  "  << SC.true_nu_py <<"  " << SC.true_nu_pz<< "  "  << nu_theta<<   std::endl;
        
        // True nue phi in BNB coordinates (around beam dir)
        double nu_phi = atan2(SC.true_nu_py, SC.true_nu_px) * 180 / 3.1415;

        // True nue theta in BNB coordinates (up from beam dir)
        double elec_theta = acos(SC.elec_pz/p_elec) * 180 / 3.1415;
        
        // True nue phi in BNB coordinates (around beam dir)
        double elec_phi = atan2(SC.elec_py, SC.elec_px) * 180 / 3.1415;
        
        // True nue angle from numi beamline 
        double nu_angle = _util.GetNuMIAngle(SC.true_nu_px, SC.true_nu_py, SC.true_nu_pz, "beam"); 

        // True nue angle from numi target
        double nu_angle_targ = _util.GetNuMIAngle(SC.true_nu_px, SC.true_nu_py, SC.true_nu_pz, "target"); 

        // True electron angle wrt numi beamline
        double elec_ang_targ = _util.GetNuMIAngle(SC.elec_px, SC.elec_py, SC.elec_pz, "target");

        // True nue theta in NuMI coordinates (up from beam dir)
        double nu_theta_numi = _util.GetNuMIAngle(SC.true_nu_px, SC.true_nu_py, SC.true_nu_pz, "numi_theta");
        
        // True nue phi in NuMI coordinates (around beam dir)
        double nu_phi_numi = _util.GetNuMIAngle(SC.true_nu_px, SC.true_nu_py, SC.true_nu_pz, "numi_phi");

        // Also require in FV
        if ( (classification_index == _util.k_nue_cc || classification_index == _util.k_nuebar_cc || classification_index == _util.k_unmatched_nue || classification_index == _util.k_cosmic_nue || classification_index == _util.k_unmatched_nuebar || classification_index == _util.k_cosmic_nuebar) && true_in_fv ){
            TH1D_true_hists.at(index).at(k_true_nue_theta)->Fill(nu_theta, weight);
            TH1D_true_hists.at(index).at(k_true_nue_phi_numi)  ->Fill(nu_phi_numi, weight);
            TH1D_true_hists.at(index).at(k_true_nue_theta_numi)->Fill(nu_theta_numi, weight);
            TH1D_true_hists.at(index).at(k_true_nue_phi)  ->Fill(nu_phi, weight);
            TH1D_true_hists.at(index).at(k_true_nue_angle)->Fill(nu_angle, weight);
            TH1D_true_hists.at(index).at(k_true_nue_px)   ->Fill(SC.true_nu_px, weight);
            TH1D_true_hists.at(index).at(k_true_nue_py)   ->Fill(SC.true_nu_py, weight);
            TH1D_true_hists.at(index).at(k_true_nue_pz)   ->Fill(SC.true_nu_pz, weight);
            TH1D_true_hists.at(index).at(k_true_nue_e)    ->Fill(SC.nu_e, weight);
            TH1D_true_hists.at(index).at(k_true_nue_p)    ->Fill(p, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_x)    ->Fill(SC.true_nu_vtx_x, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_y)    ->Fill(SC.true_nu_vtx_y, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_z)    ->Fill(SC.true_nu_vtx_z, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_x_sce)->Fill(SC.true_nu_vtx_sce_x, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_y_sce)->Fill(SC.true_nu_vtx_sce_y, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_z_sce)->Fill(SC.true_nu_vtx_sce_z, weight);
            TH1D_true_hists.at(index).at(k_true_nu_ang_targ)  ->Fill(nu_angle_targ, weight);
            TH1D_true_hists.at(index).at(k_true_elec_ang_targ)->Fill(elec_ang_targ, weight);

            TH1D_true_hists.at(index).at(k_true_elec_E)    ->Fill(SC.elec_e, weight);
            TH1D_true_hists.at(index).at(k_true_elec_theta)->Fill(elec_theta, weight);
            TH1D_true_hists.at(index).at(k_true_elec_phi)  ->Fill(elec_phi, weight);

            TH2D_true_hists.at(index).at(k_true_nue_phi_theta)   ->Fill(nu_phi, nu_theta, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_theta)->Fill(SC.nu_e, nu_theta, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_phi)  ->Fill(SC.nu_e, nu_phi, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_angle)->Fill(SC.nu_e, nu_angle, weight);
            TH2D_true_hists.at(index).at(k_true_nue_vtx_z_y)       ->Fill(SC.true_nu_vtx_z,  SC.true_nu_vtx_y, weight);
            TH2D_true_hists.at(index).at(k_true_nue_vtx_z_y_sce)   ->Fill(SC.true_nu_vtx_sce_z,  SC.true_nu_vtx_sce_y, weight);
            

            if (_type == _util.k_mc){ 

                // True vs reco histograms
                TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E)->Fill(SC.elec_e,  SC.shr_energy_cali/0.83, weight);
                TH2D_true_hists.at(index).at(k_true_nu_E_reco_nu_E)    ->Fill(SC.nu_e,  reco_nu_e, weight);
                TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E_extra_bins)->Fill(SC.elec_e,  SC.shr_energy_cali/0.83, weight);
                TH2D_true_hists.at(index).at(k_true_nu_E_reco_nu_E_extra_bins)    ->Fill(SC.nu_e,  reco_nu_e, weight);
                TH2D_true_hists.at(index).at(k_true_nu_vtx_x_reco_nu_vtx_x)->Fill(SC.true_nu_vtx_sce_x,  SC.reco_nu_vtx_sce_x, weight);
                TH2D_true_hists.at(index).at(k_true_nu_vtx_y_reco_nu_vtx_y)->Fill(SC.true_nu_vtx_sce_y,  SC.reco_nu_vtx_sce_y, weight);
                TH2D_true_hists.at(index).at(k_true_nu_vtx_z_reco_nu_vtx_z)->Fill(SC.true_nu_vtx_sce_z,  SC.reco_nu_vtx_sce_z, weight);

                // True nue interaction histograms
                if (interaction == "nue_cc_qe" || interaction == "nue_bar_cc_qe"){
                    TH1D_interaction_hists.at(index).at(_util.k_plot_qe)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_res" || interaction == "nue_bar_cc_res"){
                    TH1D_interaction_hists.at(index).at(_util.k_plot_res)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_dis" || interaction == "nue_bar_cc_dis"){
                    TH1D_interaction_hists.at(index).at(_util.k_plot_dis)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_coh" || interaction == "nue_bar_cc_coh"){
                    TH1D_interaction_hists.at(index).at(_util.k_plot_coh)->Fill(SC.nu_e, weight);
                }
                else if (interaction == "nue_cc_mec" || interaction == "nue_bar_cc_mec"){
                    TH1D_interaction_hists.at(index).at(_util.k_plot_mec)->Fill(SC.nu_e, weight);
                }
                // NC
                else {
                    TH1D_interaction_hists.at(index).at(_util.k_plot_nc)->Fill(SC.nu_e, weight);
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
            TH1D_flash_hists.at(k_flash_time_single_bin)->Fill(SC.flash_time + 0.055 -0.359, weight); // See numi documentation page to see what these numbers mean
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
            TH1D_flash_hists.at(k_flash_time_single_bin)->Fill(SC.flash_time + 0.055 -0.359, weight);
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
            TH1D_flash_hists.at(k_flash_time_single_bin)->Fill(SC.flash_time + -0.359, weight);
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
            TH1D_flash_hists.at(k_flash_time_single_bin)->Fill(SC.flash_time, weight);
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

        // For the comparisons of dedx and shr vtx distance with the moliere average we make the histogram just before the cut is applied
        if (cut_index == _util.k_shr_moliere_avg - 1 ){
            
            // This is the signal -- dont include the unmathced cases since they should all be removed by this point
            if (classification_index == _util.k_nue_cc || classification_index == _util.k_nuebar_cc){
                TH2D_hists.at(k_reco_shr_dEdx_moliere).at(_util.k_signal)->Fill(SC.shr_tkfit_dedx_Y, SC.shrmoliereavg, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_moliere_shr_dist).at(_util.k_signal)->Fill(SC.shrmoliereavg, SC.shr_distance, weight);
            }
            // This is the background
            else {
                TH2D_hists.at(k_reco_shr_dEdx_moliere).at(_util.k_background)->Fill(SC.shr_tkfit_dedx_Y, SC.shrmoliereavg, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_moliere_shr_dist).at(_util.k_background)->Fill(SC.shrmoliereavg, SC.shr_distance, weight);
            }

        }
        
        // For the dEdx vs shower distance we make the histogram just before the cut is applied
        if (cut_index == _util.k_vtx_dist_dedx - 1 ){
            
            // This is the signal -- dont include the unmathced cases since they should all be removed by this point
            if (classification_index == _util.k_nue_cc || classification_index == _util.k_nuebar_cc){
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(_util.k_signal)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist).at(_util.k_signal)->Fill(dedx_max, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist_large_dedx).at(_util.k_signal)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
            }
            // This is the background
            else {
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist).at(_util.k_background)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist).at(_util.k_background)->Fill(dedx_max, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist_large_dedx).at(_util.k_background)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
            }

        }

        // For the dEdx vs shower distance we make the histogram after the cut is applied
        if (cut_index == _util.k_vtx_dist_dedx ){
            
            // This is the signal -- dont include the unmathced cases since they should all be removed by this point
            if (classification_index == _util.k_nue_cc || classification_index == _util.k_nuebar_cc){
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist_post).at(_util.k_signal)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist_post).at(_util.k_signal)->Fill(dedx_max, SC.shr_distance, weight);
            }
            // This is the background
            else {
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_shr_dist_post).at(_util.k_background)->Fill(SC.shr_tkfit_dedx_Y, SC.shr_distance, weight);
                if (SC.n_tracks > 0) TH2D_hists.at(k_reco_shr_dEdx_max_shr_dist_post).at(_util.k_background)->Fill(dedx_max, SC.shr_distance, weight);
            }

        }

    }
    

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
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar" || classification == "cosmic_nue" || classification == "cosmic_nuebar") TEfficiency_hists.at(k_eff_nu_E).at(cut_index)->Fill(SC.nu_e, weight);
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar" || classification == "cosmic_nue" || classification == "cosmic_nuebar") TEfficiency_hists.at(k_eff_elec_E).at(cut_index)->Fill(SC.elec_e, weight);
    }
    // After this, we consider the low purity (cosmics) interactions background -- keep the unreconstructed stuff, but this shouldnt affect the plots all that much
    else {
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar") TEfficiency_hists.at(k_eff_nu_E).at(cut_index)->Fill(SC.nu_e, weight);
        if (classification == "nue_cc" || classification == "nuebar_cc" || classification == "unmatched_nue" || classification == "unmatched_nuebar") TEfficiency_hists.at(k_eff_elec_E).at(cut_index)->Fill(SC.elec_e, weight);
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
    
    // TH1D
    for (unsigned int p = 0; p < TH1D_interaction_hists.at(0).size(); p++){
        TH1D_interaction_hists.at(0).at(p)->SetOption("hist,E");
        TH1D_interaction_hists.at(0).at(p)->Write("",TObject::kOverwrite);
        TH1D_interaction_hists.at(1).at(p)->SetOption("hist,E");
        TH1D_interaction_hists.at(1).at(p)->Write("",TObject::kOverwrite);
    }
    
    
}
// -----------------------------------------------------------------------------
void HistogramHelper::Write_2DSigBkgHists(){
    
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
void HistogramHelper::FillPiZeroHists(int classification_index, SliceContainer SC, double weight, int pizero_mode){

    if (pizero_mode == 0){
        TH1D_pi0_hists.at(k_pi0_mass).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
    }
    // Norm fix
    else if (pizero_mode == 1){
        TH1D_pi0_hists.at(k_pi0_mass_norm).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
    }
    // Energy dependent
    else {
        TH1D_pi0_hists.at(k_pi0_mass_EScale).at(classification_index)->Fill(SC.pi0_mass_Y, weight);
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