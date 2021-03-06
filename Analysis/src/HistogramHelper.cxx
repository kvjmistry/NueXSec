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
            
            // Ratio hits from all showers to slice
            TH1D_hists.at(k_reco_hits_ratio).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);
            TH1D_hists.at(k_reco_hits_ratio_th).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_th_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);

            // Ratio hits from leading shower to slice
            TH1D_hists.at(k_reco_hits_ratio_ldg).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_ldg_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);
            
            // Ratio of hits of leading shower to slice with threshold on hits
            TH1D_hists.at(k_reco_hits_ratio_ldg_th).at(i).at(j) = new TH1D (Form("h_reco_hits_ratio_ldg_th_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 21, 0, 1.05);
            
            // Shower score
            TH1D_hists.at(k_reco_shower_score).at(i).at(j) = new TH1D (Form("h_reco_shower_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 0.5);

            // Track score
            TH1D_hists.at(k_reco_track_score).at(i).at(j) = new TH1D (Form("h_reco_track_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0.5, 1);
            
            // Calibrated energy of just the leading shower
            TH1D_hists.at(k_reco_shower_energy_cali).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_cali_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 4);

            // Angle between the vector from target to nu vtx and the shower direction
            TH1D_hists.at(k_reco_effective_angle).at(i).at(j) = new TH1D ( Form("h_reco_effective_angle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 13, 0, 190);

            // Cosine of the Angle between the vector from target to nu vtx and the shower direction
            TH1D_hists.at(k_reco_effective_cosangle).at(i).at(j) = new TH1D ( Form("h_reco_effective_cosangle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 16, -1, 1);

            // Set the bins for the reco energy
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
            TH1D_hists.at(k_reco_shower_energy_cali_rebin).at(i).at(j) = new TH1D (Form("h_reco_shower_energy_cali_rebin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", _util.reco_shr_bins.size()-1, edges);

            edges = &_util.reco_shr_bins_ang[0]; // Cast to an array 
            TH1D_hists.at(k_reco_effective_angle_rebin).at(i).at(j) = new TH1D (Form("h_reco_effective_angle_rebin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", _util.reco_shr_bins_ang.size()-1, edges);

            edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
            TH1D_hists.at(k_reco_effective_cosangle_rebin).at(i).at(j) = new TH1D (Form("h_reco_effective_cosangle_rebin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", _util.reco_shr_bins_cang.size()-1, edges);

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
            TH1D_hists.at(k_reco_shr_tkfit_dedx_max_no_tracks).at(i).at(j) = new TH1D ( Form("h_reco_shr_tkfit_dedx_max_no_tracks_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 10);

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

            // Single Bin
            TH1D_hists.at(k_reco_single_bin).at(i).at(j) = new TH1D ( Form("h_reco_single_bin_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 1, 0, 40);

            // Track PID score in the event
            TH1D_hists.at(k_reco_trk_pid_score).at(i).at(j) = new TH1D ( Form("h_reco_trk_pid_score_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, -1, 1);

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
                TEfficiency_hists.at(k_eff_elec_E).at(l)                   = new TH1D( Form("h_true_elec_E_%s",                  _util.cut_dirs.at(l).c_str() ), "", 8, 0, 6 );
                TEfficiency_hists.at(k_eff_elec_E_many_bins).at(l)         = new TH1D( Form("h_eff_elec_E_many_bins_%s",         _util.cut_dirs.at(l).c_str() ), "", 100, 0, 1 );
                
                double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
                TEfficiency_hists.at(k_eff_elec_E_rebin).at(l)             = new TH1D( Form("h_true_elec_E_rebin_%s",            _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_elec_E_rebin_nue).at(l)         = new TH1D( Form("h_true_elec_E_rebin_nue_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
                TEfficiency_hists.at(k_eff_elec_E_rebin_nuebar).at(l)      = new TH1D( Form("h_true_elec_E_rebin_nuebar_%s",     _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins.size()-1, edges);
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
                
                TEfficiency_hists.at(k_eff_cosine_beta).at(l)              = new TH1D( Form("h_eff_cosine_beta_%s",              _util.cut_dirs.at(l).c_str() ), "", 16, -1, 1 );
                edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
                TEfficiency_hists.at(k_eff_cosine_beta_rebin).at(l)        = new TH1D( Form("h_eff_cosine_beta_rebin_%s",        _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                TEfficiency_hists.at(k_eff_cosine_beta_rebin_nue).at(l)    = new TH1D( Form("h_eff_cosine_beta_rebin_nue_%s",    _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                TEfficiency_hists.at(k_eff_cosine_beta_rebin_nuebar).at(l) = new TH1D( Form("h_eff_cosine_beta_rebin_nuebar_%s", _util.cut_dirs.at(l).c_str() ), "", _util.reco_shr_bins_cang.size()-1, edges );
                
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
            TH1D_true_hists.at(i).at(k_true_nue_theta) = new TH1D( Form("h_true_nue_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} #theta [deg]; Entries",           14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_nue_phi)   = new TH1D( Form("h_true_nue_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} #phi [deg]; Entries",              14, 0, 40);
            
            TH1D_true_hists.at(i).at(k_true_nue_angle) = new TH1D( Form("h_true_nue_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} Angle [deg]; Entries", 120, 0, 120 );
            TH1D_true_hists.at(i).at(k_true_nue_px)    = new TH1D( Form("h_true_nue_px_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} P_x [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_py)    = new TH1D( Form("h_true_nue_py_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} P_y [GeV/c]; Entries", 14, 0, 1);
            TH1D_true_hists.at(i).at(k_true_nue_pz)    = new TH1D( Form("h_true_nue_pz_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} P_z [GeV/c]; Entries", 14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_e)     = new TH1D( Form("h_true_nue_e_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} E [GeV]; Entries",    15, 0, 5);
            TH1D_true_hists.at(i).at(k_true_nue_p)     = new TH1D( Form("h_true_nue_p_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} P [GeV/c]; Entries",  14, 0, 5);
            TH1D_true_hists.at(i).at(k_true_vtx_x)     = new TH1D( Form("h_true_vtx_x_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx x [cm]; Entries", 20, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y)     = new TH1D( Form("h_true_vtx_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx y [cm]; Entries", 20, -10, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z)     = new TH1D( Form("h_true_vtx_z_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx z [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_vtx_x_sce) = new TH1D( Form("h_true_vtx_x_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx x SCE Corr. [cm]; Entries", 20, -10, 270);
            TH1D_true_hists.at(i).at(k_true_vtx_y_sce) = new TH1D( Form("h_true_vtx_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx y SCE Corr. [cm]; Entries", 20, -10, 120);
            TH1D_true_hists.at(i).at(k_true_vtx_z_sce) = new TH1D( Form("h_true_vtx_z_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ) ,";True #nu_{e} + #bar{#nu}_{e} Vtx z SCE Corr. [cm]; Entries", 40, -10, 1050);
            TH1D_true_hists.at(i).at(k_true_nu_ang_targ)     = new TH1D( Form("h_true_nu_ang_targ_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True #nu_{e} + #bar{#nu}_{e} Angle wrt NuMI Targ Vec [deg]; Entries",  40, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_ang_targ)   = new TH1D( Form("h_true_elec_ang_targ_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} + e^{+} Angle wrt NuMI Targ Vec [deg]; Entries", 25, 0, 180 );

            TH1D_true_hists.at(i).at(k_true_elec_E)      = new TH1D( Form("h_true_elec_E_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries",                15, 0, 5 );
            TH1D_true_hists.at(i).at(k_true_elec_theta)  = new TH1D( Form("h_true_elec_theta_%s_%s",  _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} + e^{+} #theta [deg]; Entries",            14, 0, 180 );
            TH1D_true_hists.at(i).at(k_true_elec_phi)    = new TH1D( Form("h_true_elec_phi_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), ";True e#lower[-0.5]{-} + e^{+} #phi [deg]; Entries",              25, -180, 180);

            TH1D_true_hists.at(i).at(k_reco_true_ang)    = new TH1D( Form("h_reco_true_ang_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str() ), "; Angle (Reco to True) Nu Dir [deg]; Entries",              25, -2, 20);

            TH2D_true_hists.at(i).at(k_true_nue_phi_theta)    = new TH2D( Form("h_true_nue_phi_theta_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} + #bar{#nu}_{e} #phi [deg];True #nu_{e} #theta [deg]",     14, 0, 80, 16, 20, 140 );
            TH2D_true_hists.at(i).at(k_true_nue_energy_theta) = new TH2D( Form("h_true_nue_energy_theta_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} + #bar{#nu}_{e} E [GeV];True #nu_{e} #theta [deg]",           25, 0, 5, 16, 0, 140);
            TH2D_true_hists.at(i).at(k_true_nue_energy_phi)   = new TH2D( Form("h_true_nue_energy_phi_%s_%s",   _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} + #bar{#nu}_{e} E [GeV];True #nu_{e} #phi [deg]",             25, 0, 5, 14, 0, 75);
            TH2D_true_hists.at(i).at(k_true_nue_energy_angle) = new TH2D( Form("h_true_nue_energy_angle_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True #nu_{e} + #bar{#nu}_{e} E [GeV];True #nu_{e} Angle from NuMI [deg]", 25, 0, 5, 30, 0, 120);
        
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y)     = new TH2D( Form("h_true_nue_vtx_z_y_%s_%s",     _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Vtx Z [cm] ;True #nu_{e} + #bar{#nu}_{e} Vtx Y [cm]", 40, -10, 1050, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nue_vtx_z_y_sce) = new TH2D( Form("h_true_nue_vtx_z_y_sce_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),";True #nu_{e} + #bar{#nu}_{e} Vtx Z  SCE Corr. [cm];True #nu_{e} + #bar{#nu}_{e} Vtx Y SCE Corr. [cm]", 40, -10, 1050, 20, -10, 120);
        
            double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E) = new TH2D( Form("h_true_elec_E_reco_elec_E_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True e#lower[-0.5]{-} + e^{+} Energy [GeV] ;E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins.size()-1, edges);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E)     = new TH2D( Form("h_true_nu_E_reco_nu_E_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Energy [GeV] ;Reco #nu_{e} + #bar{#nu}_{e} Energy [GeV]", 25, 0, 4, 25, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_%s_%s",               _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e#lower[-0.5]{-} + e^{+} Energy [GeV] ;E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins_nue) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_nue_%s_%s",       _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e#lower[-0.5]{-} Energy [GeV] ;E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_extra_bins_nuebar) = new TH2D( Form("h_true_elec_E_reco_elec_E_extra_bins_nuebar_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()), ";True e^{+} Energy [GeV] ;E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]", 50, 0, 4, 50, 0, 4);
            TH2D_true_hists.at(i).at(k_true_nu_E_reco_nu_E_extra_bins)     = new TH2D( Form("h_true_nu_E_reco_nu_E_extra_bins_%s_%s",             _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Energy [GeV] ;Reco #nu_{e} + #bar{#nu}_{e} Energy [GeV]", 50, 0, 4, 50, 0, 4);

            TH2D_true_hists.at(i).at(k_true_nu_vtx_x_reco_nu_vtx_x)  = new TH2D( Form("h_true_nu_vtx_x_reco_nu_vtx_x_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Vtx X [cm] ;Reco #nu_{e} + #bar{#nu}_{e} Vtx X [cm]", 20, -10, 270, 20, -10, 270);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_y_reco_nu_vtx_y)  = new TH2D( Form("h_true_nu_vtx_y_reco_nu_vtx_y_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Vtx Y [cm] ;Reco #nu_{e} + #bar{#nu}_{e} Vtx Y [cm]", 20, -10, 120, 20, -10, 120);
            TH2D_true_hists.at(i).at(k_true_nu_vtx_z_reco_nu_vtx_z)  = new TH2D( Form("h_true_nu_vtx_z_reco_nu_vtx_z_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";True #nu_{e} + #bar{#nu}_{e} Vtx Z [cm] ;Reco #nu_{e} + #bar{#nu}_{e} Vtx Z [cm]", 40, -10, 1050, 40, -10, 1050);
            TH2D_true_hists.at(i).at(k_true_shr_energy_purity)       = new TH2D( Form("h_true_shr_energy_purity_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV] ;Shower Purity", _util.reco_shr_bins.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_energy_completeness) = new TH2D( Form("h_true_shr_energy_completeness_%s_%s",_util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV] ;Shower Completeness", _util.reco_shr_bins.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_energy_resolution_reco) = new TH2D( Form("h_true_shr_energy_resolution_reco_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV]; Reco - True / Reco", _util.reco_shr_bins.size()-1, edges, 30, -1.2, 1.2);
            TH2D_true_hists.at(i).at(k_true_shr_energy_resolution_true) = new TH2D( Form("h_true_shr_energy_resolution_true_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";E^{reco}_{e#lower[-0.5]{-} + e^{+}} [GeV] ;Reco - True / True", _util.reco_shr_bins.size()-1, edges, 30, -1.2, 1.2);
        
            TH2D_true_hists.at(i).at(k_true_elec_E_reco_elec_E_rebin) = new TH2D( Form("h_true_elec_E_reco_elec_E_rebin_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";Generated Electron Energy [GeV] ;Measured Electron Energy [GeV]", _util.reco_shr_bins.size()-1, edges, _util.reco_shr_bins.size()-1, edges);


            edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_purity)          = new TH2D( Form("h_true_shr_cosbeta_purity_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e#lower[-0.5]{-} + e^{+}} ;Shower Purity", _util.reco_shr_bins_cang.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_completeness)    = new TH2D( Form("h_true_shr_cosbeta_completeness_%s_%s",_util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e#lower[-0.5]{-} + e^{+}} ;Shower Completeness", _util.reco_shr_bins_cang.size()-1, edges, 21, 0, 1.1);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_resolution_reco) = new TH2D( Form("h_true_shr_cosbeta_resolution_reco_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e#lower[-0.5]{-} + e^{+}}; Reco - True / Reco", _util.reco_shr_bins_cang.size()-1, edges, 30, -1.2, 1.2);
            TH2D_true_hists.at(i).at(k_true_shr_cosbeta_resolution_true) = new TH2D( Form("h_true_shr_cosbeta_resolution_true_%s_%s", _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";cos#beta^{reco}_{e#lower[-0.5]{-} + e^{+}} ;Reco - True / True", _util.reco_shr_bins_cang.size()-1, edges, 30, -1.2, 1.2);

            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta)      = new TH2D( Form("h_elec_true_beta_reco_beta_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),            ";True e#lower[-0.5]{-} + e^{+} #beta [deg]; Reco e#lower[-0.5]{-} + e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta_nue)    = new TH2D( Form("h_elec_true_beta_reco_beta_nue_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),      ";True e#lower[-0.5]{-} #beta [deg]; Reco e#lower[-0.5]{-} + e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_beta_reco_beta_nuebar) = new TH2D( Form("h_elec_true_beta_reco_beta_nuebar_%s_%s",    _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e^{+} #beta [deg]; Reco e#lower[-0.5]{-} + e^{+} #beta [deg]",  36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_theta_reco_theta)    = new TH2D( Form("h_elec_true_theta_reco_theta_%s_%s",  _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e#lower[-0.5]{-} + e^{+} #theta [deg];Reco e#lower[-0.5]{-} + e^{+} #theta [deg]", 36, 0, 180, 36, 0, 180 );
            TH2D_true_hists.at(i).at(k_elec_true_phi_reco_phi)        = new TH2D( Form("h_elec_true_phi_reco_phi_%s_%s",      _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),   ";True e#lower[-0.5]{-} + e^{+} #phi [deg];  Reco e#lower[-0.5]{-} + e^{+} #phi [deg]",   36, -180, 180, 36, -180, 180 );
        
            TH2D_true_hists.at(i).at(k_elec_true_cosbeta_reco_cosbeta_rebin) = new TH2D( Form("h_elec_true_cosbeta_reco_cosbeta_rebin_%s_%s",         _util.type_prefix.at(_type).c_str(), cut_stage.c_str()),    ";Generated Electron cos#beta ;Measured Electron cos#beta", _util.reco_shr_bins_cang.size()-1, edges, _util.reco_shr_bins_cang.size()-1, edges);
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

            // TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(i).at(p)        = new TH1D( Form("h_true_nue_nuebar_E_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e} Energy [GeV]; Entries",   15, 0, 6 );
            // TH1D_interaction_hists.at(k_int_nu_E_nue).at(i).at(p)               = new TH1D( Form("h_true_nue_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} Energy [GeV]; Entries",                  15, 0, 6 );
            // TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(i).at(p)            = new TH1D( Form("h_true_nuebar_E_%s_%s",              _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{#nu}_{e} Energy [GeV]; Entries",            15, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(i).at(p)        = new TH1D( Form("h_true_nue_nuebar_E_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e} Energy [GeV]; Entries",   8, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_nue).at(i).at(p)               = new TH1D( Form("h_true_nue_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} Energy [GeV]; Entries",                  8, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(i).at(p)            = new TH1D( Form("h_true_nuebar_E_%s_%s",              _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{#nu}_{e} Energy [GeV]; Entries",            8, 0, 6 );
            TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(i).at(p)        = new TH1D( Form("h_int_nu_E_single_bin_%s_%s",        _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e}; Entries",                               1, 0, 10 );
            TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(i).at(p)    = new TH1D( Form("h_int_nu_E_nue_single_bin_%s_%s",    _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #bar{nu}_{e}; Entries",                          1, 0, 10 );
            TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(i).at(p) = new TH1D( Form("h_int_nu_E_nuebar_single_bin_%s_%s", _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True #nu_{e} + #bar{nu}_{e}; Entries",                1, 0, 10 );
            TH1D_interaction_hists.at(k_int_elec_E).at(i).at(p)                 = new TH1D( Form("h_int_elec_E_%s_%s",                 _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries",      8, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_nue).at(i).at(p)             = new TH1D( Form("h_int_elec_E_nue_%s_%s",             _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e#lower[-0.5]{-} [GeV]; Entries",                     8, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_nuebar).at(i).at(p)          = new TH1D( Form("h_int_elec_E_nuebar_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; e^{+} Energy [GeV]; Entries",                         8, 0, 6 );
            TH1D_interaction_hists.at(k_int_elec_E_rebin).at(i).at(p)           = new TH1D( Form("h_int_elec_E_rebin_%s_%s",           _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries", _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_E_rebin_nue).at(i).at(p)       = new TH1D( Form("h_int_elec_E_rebin_nue_%s_%s",       _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} Energy [GeV]; Entries",         _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar).at(i).at(p)    = new TH1D( Form("h_int_elec_E_rebin_nuebar_%s_%s",    _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e^{+} Energy [GeV]; Entries",                    _util.reco_shr_bins.size()-1, edges );
            TH1D_interaction_hists.at(k_int_elec_theta).at(i).at(p)             = new TH1D( Form("h_int_elec_theta_%s_%s",             _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_elec_phi).at(i).at(p)               = new TH1D( Form("h_int_elec_phi_%s_%s",               _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries", 14, -190, 190 );
            TH1D_interaction_hists.at(k_int_effective_ang).at(i).at(p)          = new TH1D( Form("h_int_effective_ang_%s_%s",          _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} Energy [GeV]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_beta_nue).at(i).at(p)               = new TH1D( Form("h_int_beta_nue_%s_%s",               _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} #beta [deg]; Entries", 13, 0, 190 );
            TH1D_interaction_hists.at(k_int_beta_nuebar).at(i).at(p)            = new TH1D( Form("h_int_beta_nuebar_%s_%s",            _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} #beta [deg]; Entries", 13, 0, 190 );
            edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
            TH1D_interaction_hists.at(k_int_cosbeta).at(i).at(p)                = new TH1D( Form("h_int_cosbeta_%s_%s",                _util.interaction_types.at(p).c_str(), stage.c_str() ), "; True e#lower[-0.5]{-} + e^{+} cos#beta; Entries", _util.reco_shr_bins_cang.size()-1, edges );
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
    bool true_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z);

    // Get the trkfit dedx max variable
    double dedx_max = SC.GetdEdxMax();
    
    // Now fill the histograms!
    TH1D_hists.at(k_reco_vtx_x).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x, weight);
    TH1D_hists.at(k_reco_vtx_y).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y, weight);
    TH1D_hists.at(k_reco_vtx_z).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z, weight);

    TH1D_hists.at(k_reco_vtx_x_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_x, weight);
    TH1D_hists.at(k_reco_vtx_y_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_y, weight);
    TH1D_hists.at(k_reco_vtx_z_sce).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_sce_z, weight);
 
    TH1D_hists.at(k_reco_leading_mom).at(cut_index).at(classification_index)->Fill(SC.shr_p, weight);
    
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

    if (SC.shr_hits_max > 50) TH1D_hists.at(k_reco_hits_ratio_th).at(cut_index).at(classification_index)->Fill(SC.hits_ratio, weight);

    if (SC.slnhits > 0) TH1D_hists.at(k_reco_hits_ratio_ldg).at(cut_index).at(classification_index)->Fill(double(SC.shr_hits_max)/double(SC.slnhits), weight);
    
    if (SC.shr_hits_max > 50 && SC.slnhits > 0) TH1D_hists.at(k_reco_hits_ratio_ldg_th).at(cut_index).at(classification_index)->Fill(double(SC.shr_hits_max)/double(SC.slnhits), weight);

    TH1D_hists.at(k_reco_shower_score).at(cut_index).at(classification_index)->Fill(SC.shr_score, weight);

    TH1D_hists.at(k_reco_track_score).at(cut_index).at(classification_index)->Fill(SC.trk_score, weight);

    TH1D_hists.at(k_reco_shower_energy_cali).at(cut_index).at(classification_index)->Fill(SC.shr_energy_cali, weight);
    TH1D_hists.at(k_reco_shower_energy_cali_rebin).at(cut_index).at(classification_index)->Fill(SC.shr_energy_cali, weight);

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

    TH1D_hists.at(k_reco_nu_e).at(cut_index).at(classification_index)->Fill(SC.reco_e, weight);

    TH1D_hists.at(k_reco_contained_fraction).at(cut_index).at(classification_index)->Fill(SC.contained_fraction, weight);

    TH1D_hists.at(k_reco_run_number).at(cut_index).at(classification_index)->Fill(SC.run, weight);

    TH1D_hists.at(k_reco_nu_purity_from_pfp).at(cut_index).at(classification_index)->Fill(SC.nu_purity_from_pfp, weight);

    TH1D_hists.at(k_reco_crtveto).at(cut_index).at(classification_index)->Fill(SC.crtveto, weight);

    TH1D_hists.at(k_reco_crthitpe).at(cut_index).at(classification_index)->Fill(SC.crthitpe, weight);

    TH1D_hists.at(k_reco_shr_ang_numi).at(cut_index).at(classification_index)->Fill(SC.shr_ang_numi, weight);

    TH1D_hists.at(k_reco_single_bin).at(cut_index).at(classification_index)->Fill(weight);

    TH1D_hists.at(k_reco_effective_angle).at(cut_index).at(classification_index)->Fill(SC.effective_angle, weight);
    TH1D_hists.at(k_reco_effective_angle_rebin).at(cut_index).at(classification_index)->Fill(SC.effective_angle, weight);
    
    TH1D_hists.at(k_reco_effective_cosangle).at(cut_index).at(classification_index)->Fill(SC.cos_effective_angle, weight);
    TH1D_hists.at(k_reco_effective_cosangle_rebin).at(cut_index).at(classification_index)->Fill(SC.cos_effective_angle, weight);
    
    for (unsigned int trk = 0; trk < SC.trk_llr_pid_score_v->size(); trk++){
        TH1D_hists.at(k_reco_trk_pid_score).at(cut_index).at(classification_index)->Fill(SC.trk_llr_pid_score_v->at(trk), weight);
    }

    if (SC.n_showers > 1)
        TH1D_hists.at(k_reco_pi0mass).at(cut_index).at(classification_index)->Fill(SC.pi0_mass_Y, weight);

    
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // Partice Type Hists
    TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_max_par).at(cut_index).at(_par_type)->Fill(dedx_max, weight);
    TH1D_hists_particle.at(k_reco_shr_tkfit_dedx_y_par).at(cut_index).at(_par_type)->Fill(SC.shr_tkfit_dedx_Y, weight);

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    // 2D Histograms
    TH2D_hists_cuts.at(k_2D_dedx_shower_energy).at(cut_index).at(classification_index)->Fill(dedx_max, SC.shr_energy_cali, weight);

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
        
        // Also require in FV
        if ( (classification_index == _util.k_nue_cc || classification_index == _util.k_nuebar_cc ||
              classification_index == _util.k_unmatched_nue || classification_index == _util.k_cosmic_nue || 
              classification_index == _util.k_unmatched_nuebar || classification_index == _util.k_cosmic_nuebar) && true_in_fv ){
            
            // Just True histograms
            TH1D_true_hists.at(index).at(k_true_nue_theta)       ->Fill(SC.nu_theta, weight);
            TH1D_true_hists.at(index).at(k_true_nue_phi)         ->Fill(SC.nu_phi, weight);
            TH1D_true_hists.at(index).at(k_true_nue_angle)       ->Fill(SC.nu_angle, weight);
            TH1D_true_hists.at(index).at(k_true_nue_px)          ->Fill(SC.true_nu_px, weight);
            TH1D_true_hists.at(index).at(k_true_nue_py)          ->Fill(SC.true_nu_py, weight);
            TH1D_true_hists.at(index).at(k_true_nue_pz)          ->Fill(SC.true_nu_pz, weight);
            TH1D_true_hists.at(index).at(k_true_nue_e)           ->Fill(SC.nu_e, weight);
            TH1D_true_hists.at(index).at(k_true_nue_p)           ->Fill(SC.nu_p, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_x)           ->Fill(SC.true_nu_vtx_x, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_y)           ->Fill(SC.true_nu_vtx_y, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_z)           ->Fill(SC.true_nu_vtx_z, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_x_sce)       ->Fill(SC.true_nu_vtx_sce_x, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_y_sce)       ->Fill(SC.true_nu_vtx_sce_y, weight);
            TH1D_true_hists.at(index).at(k_true_vtx_z_sce)       ->Fill(SC.true_nu_vtx_sce_z, weight);
            TH1D_true_hists.at(index).at(k_true_nu_ang_targ)     ->Fill(SC.nu_angle_targ, weight);
            TH1D_true_hists.at(index).at(k_true_elec_ang_targ)   ->Fill(SC.elec_ang_targ, weight);
            TH1D_true_hists.at(index).at(k_true_elec_E)          ->Fill(SC.elec_e, weight);
            TH1D_true_hists.at(index).at(k_true_elec_theta)      ->Fill(SC.elec_theta, weight);
            TH1D_true_hists.at(index).at(k_true_elec_phi)        ->Fill(SC.elec_phi, weight);
            TH1D_true_hists.at(index).at(k_reco_true_ang)        ->Fill(SC.reco_true_nu_ang, weight);
            TH2D_true_hists.at(index).at(k_true_nue_phi_theta)   ->Fill(SC.nu_phi, SC.nu_theta, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_theta)->Fill(SC.nu_e, SC.nu_theta, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_phi)  ->Fill(SC.nu_e, SC.nu_phi, weight);
            TH2D_true_hists.at(index).at(k_true_nue_energy_angle)->Fill(SC.nu_e, SC.nu_angle, weight);
            TH2D_true_hists.at(index).at(k_true_nue_vtx_z_y)     ->Fill(SC.true_nu_vtx_z,  SC.true_nu_vtx_y, weight);
            TH2D_true_hists.at(index).at(k_true_nue_vtx_z_y_sce) ->Fill(SC.true_nu_vtx_sce_z,  SC.true_nu_vtx_sce_y, weight);

            TH2D_true_hists.at(index).at(k_elec_true_beta_reco_beta)  ->Fill(SC.true_effective_angle, SC.effective_angle, weight);
            if (SC.nu_pdg == 12)TH2D_true_hists.at(index).at(k_elec_true_beta_reco_beta_nue)     ->Fill(SC.true_effective_angle, SC.effective_angle, weight);
            if (SC.nu_pdg == -12)TH2D_true_hists.at(index).at(k_elec_true_beta_reco_beta_nuebar)  ->Fill(SC.true_effective_angle, SC.effective_angle, weight);
            TH2D_true_hists.at(index).at(k_elec_true_theta_reco_theta)->Fill(SC.elec_theta, SC.shr_theta  * 180/3.14159, weight);
            TH2D_true_hists.at(index).at(k_elec_true_phi_reco_phi)    ->Fill(SC.elec_phi, SC.shr_phi  * 180/3.14159, weight);

            // True vs reco histograms
            TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E)           ->Fill(SC.elec_e,  SC.shr_energy_cali, weight);
            TH2D_true_hists.at(index).at(k_true_nu_E_reco_nu_E)               ->Fill(SC.nu_e,    SC.reco_e, weight);
            TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E_extra_bins)->Fill(SC.elec_e,  SC.shr_energy_cali, weight);
            if (SC.nu_pdg == 12) TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E_extra_bins_nue)->Fill(SC.elec_e,  SC.shr_energy_cali, weight);
            if (SC.nu_pdg == -12)TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E_extra_bins_nuebar)->Fill(SC.elec_e,  SC.shr_energy_cali, weight);
            TH2D_true_hists.at(index).at(k_true_nu_E_reco_nu_E_extra_bins)    ->Fill(SC.nu_e,    SC.reco_e, weight);
            TH2D_true_hists.at(index).at(k_true_nu_vtx_x_reco_nu_vtx_x)       ->Fill(SC.true_nu_vtx_sce_x,  SC.reco_nu_vtx_sce_x, weight);
            TH2D_true_hists.at(index).at(k_true_nu_vtx_y_reco_nu_vtx_y)       ->Fill(SC.true_nu_vtx_sce_y,  SC.reco_nu_vtx_sce_y, weight);
            TH2D_true_hists.at(index).at(k_true_nu_vtx_z_reco_nu_vtx_z)       ->Fill(SC.true_nu_vtx_sce_z,  SC.reco_nu_vtx_sce_z, weight);

            TH2D_true_hists.at(index).at(k_true_shr_energy_purity)            ->Fill(SC.shr_energy_cali,  SC.shr_bkt_purity, weight);
            TH2D_true_hists.at(index).at(k_true_shr_energy_completeness)      ->Fill(SC.shr_energy_cali,  SC.shr_bkt_completeness, weight);

            TH2D_true_hists.at(index).at(k_true_shr_energy_resolution_reco)      ->Fill( SC.shr_energy_cali, (SC.shr_energy_cali - SC.elec_e) / SC.shr_energy_cali, weight);
            TH2D_true_hists.at(index).at(k_true_shr_energy_resolution_true)      ->Fill( SC.shr_energy_cali, (SC.shr_energy_cali - SC.elec_e) / SC.elec_e, weight);

            
            TH2D_true_hists.at(index).at(k_true_shr_cosbeta_purity)            ->Fill(SC.cos_effective_angle,  SC.shr_bkt_purity, weight);
            TH2D_true_hists.at(index).at(k_true_shr_cosbeta_completeness)      ->Fill(SC.cos_effective_angle,  SC.shr_bkt_completeness, weight);

            double true_cosbeta = std::cos(SC.true_effective_angle * 3.14159 / 180.0);
            TH2D_true_hists.at(index).at(k_true_shr_cosbeta_resolution_reco)      ->Fill( SC.cos_effective_angle, (SC.cos_effective_angle - true_cosbeta) / SC.cos_effective_angle, weight);
            TH2D_true_hists.at(index).at(k_true_shr_cosbeta_resolution_true)      ->Fill( SC.cos_effective_angle, (SC.cos_effective_angle - true_cosbeta) / true_cosbeta, weight);

            TH2D_true_hists.at(index).at(k_true_elec_E_reco_elec_E_rebin)->Fill(SC.elec_e,  SC.shr_energy_cali, weight);
            TH2D_true_hists.at(index).at(k_elec_true_cosbeta_reco_cosbeta_rebin)  ->Fill(true_cosbeta, SC.cos_effective_angle, weight);

            // True nue interaction histograms
            if (interaction == "nue_cc_qe" || interaction == "nue_bar_cc_qe"){
                
                // All 
                TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_qe)->Fill(SC.nu_e, weight);
                TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_qe)->Fill(weight);
                TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_qe)->Fill(SC.elec_theta, weight);
                TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_qe)->Fill(SC.elec_phi, weight);
                TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_qe)->Fill(SC.true_effective_angle, weight);
                TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_qe)->Fill(true_cosbeta, weight);
                
                // Nue
                if (interaction == "nue_cc_qe")     {
                    TH1D_interaction_hists.at(k_int_nu_E_nue).at(index).at(_util.k_plot_qe)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(index).at(_util.k_plot_qe)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nue)   .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nue)     .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nue).at(index).at(_util.k_plot_qe)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(index).at(_util.k_plot_qe)->Fill(true_cosbeta, weight);
                }

                // Nuebar
                if (interaction == "nue_bar_cc_qe") {
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(index).at(_util.k_plot_qe)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(index).at(_util.k_plot_qe)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar)   .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nuebar)   .at(index).at(_util.k_plot_qe)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nuebar).at(index).at(_util.k_plot_qe)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(index).at(_util.k_plot_qe)->Fill(true_cosbeta, weight);
                }
            }
            else if (interaction == "nue_cc_res" || interaction == "nue_bar_cc_res"){
                
                // All
                TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_res)->Fill(SC.nu_e, weight);
                TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_res)->Fill(weight);
                TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_res)->Fill(SC.elec_theta, weight);
                TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_res)->Fill(SC.elec_phi, weight);
                TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_res)->Fill(SC.true_effective_angle, weight);
                TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_res)->Fill(true_cosbeta, weight);
                
                // Nue
                if (interaction == "nue_cc_res"){     
                    TH1D_interaction_hists.at(k_int_nu_E_nue).at(index).at(_util.k_plot_res)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(index).at(_util.k_plot_res)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nue)   .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nue)         .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nue).at(index).at(_util.k_plot_res)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(index).at(_util.k_plot_res)->Fill(true_cosbeta, weight);
                }

                // Nuebar
                if (interaction == "nue_bar_cc_res"){ 
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(index).at(_util.k_plot_res)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(index).at(_util.k_plot_res)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar)   .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nuebar)         .at(index).at(_util.k_plot_res)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nuebar).at(index).at(_util.k_plot_res)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(index).at(_util.k_plot_res)->Fill(true_cosbeta, weight);
                }
            }
            else if (interaction == "nue_cc_dis" || interaction == "nue_bar_cc_dis"){
                
                // All
                TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_dis)->Fill(SC.nu_e, weight);
                TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_dis)->Fill(weight);
                TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_dis)->Fill(SC.elec_theta, weight);
                TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_dis)->Fill(SC.elec_phi, weight);
                TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_dis)->Fill(SC.true_effective_angle, weight);
                TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_dis)->Fill(true_cosbeta, weight);
                
                // Nue
                if (interaction == "nue_cc_dis"){     
                    TH1D_interaction_hists.at(k_int_nu_E_nue).at(index).at(_util.k_plot_dis)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(index).at(_util.k_plot_dis)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nue)   .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nue)         .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nue).at(index).at(_util.k_plot_dis)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(index).at(_util.k_plot_dis)->Fill(true_cosbeta, weight);
                }

                // Nuebar
                if (interaction == "nue_bar_cc_dis"){ 
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(index).at(_util.k_plot_dis)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(index).at(_util.k_plot_dis)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar)   .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nuebar)         .at(index).at(_util.k_plot_dis)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nuebar).at(index).at(_util.k_plot_dis)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(index).at(_util.k_plot_dis)->Fill(true_cosbeta, weight);
                }
            }
            else if (interaction == "nue_cc_coh" || interaction == "nue_bar_cc_coh"){
                
                // All
                TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_coh)->Fill(SC.nu_e, weight);
                TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_coh)->Fill(weight);
                TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_coh)->Fill(SC.elec_theta, weight);
                TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_coh)->Fill(SC.elec_phi, weight);
                TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_coh)->Fill(SC.true_effective_angle, weight);
                TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_coh)->Fill(true_cosbeta, weight);
                
                // Nue
                if (interaction == "nue_cc_coh"){     
                    TH1D_interaction_hists.at(k_int_nu_E_nue).at(index).at(_util.k_plot_coh)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(index).at(_util.k_plot_coh)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nue)   .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nue)         .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nue).at(index).at(_util.k_plot_coh)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(index).at(_util.k_plot_coh)->Fill(true_cosbeta, weight);
                }

                // Nuebar
                if (interaction == "nue_bar_cc_coh"){ 
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(index).at(_util.k_plot_coh)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(index).at(_util.k_plot_coh)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar)   .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nuebar)         .at(index).at(_util.k_plot_coh)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nuebar).at(index).at(_util.k_plot_coh)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(index).at(_util.k_plot_coh)->Fill(true_cosbeta, weight);
                }
            }
            else if (interaction == "nue_cc_mec" || interaction == "nue_bar_cc_mec"){
                
                // All
                TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_mec)->Fill(SC.nu_e, weight);
                TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_mec)->Fill(weight);
                TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_mec)->Fill(SC.elec_theta, weight);
                TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_mec)->Fill(SC.elec_phi, weight);
                TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_mec)->Fill(SC.true_effective_angle, weight);
                TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_mec)->Fill(true_cosbeta, weight);
                
                // Nue
                if (interaction == "nue_cc_mec"){     
                    TH1D_interaction_hists.at(k_int_nu_E_nue).at(index).at(_util.k_plot_mec)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nue_single_bin).at(index).at(_util.k_plot_mec)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nue)   .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nue)         .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nue).at(index).at(_util.k_plot_mec)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nue).at(index).at(_util.k_plot_mec)->Fill(true_cosbeta, weight);
                }

                // Nuebar
                if (interaction == "nue_bar_cc_mec"){ 
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar).at(index).at(_util.k_plot_mec)->Fill(SC.nu_e, weight);
                    TH1D_interaction_hists.at(k_int_nu_E_nuebar_single_bin).at(index).at(_util.k_plot_mec)->Fill(weight);
                    TH1D_interaction_hists.at(k_int_elec_E_rebin_nuebar)   .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_elec_E_nuebar)         .at(index).at(_util.k_plot_mec)->Fill(SC.elec_e, weight);
                    TH1D_interaction_hists.at(k_int_beta_nuebar).at(index).at(_util.k_plot_mec)->Fill(SC.true_effective_angle, weight);
                    TH1D_interaction_hists.at(k_int_cosbeta_rebin_nuebar).at(index).at(_util.k_plot_mec)->Fill(true_cosbeta, weight);
                }
            }
            // CC Tot
            TH1D_interaction_hists.at(k_int_nu_E_nue_nuebar).at(index).at(_util.k_plot_tot)->Fill(SC.nu_e, weight);
            TH1D_interaction_hists.at(k_int_nu_E_single_bin).at(index).at(_util.k_plot_tot)->Fill(weight);
            TH1D_interaction_hists.at(k_int_elec_E)         .at(index).at(_util.k_plot_tot)->Fill(SC.elec_e, weight);
            TH1D_interaction_hists.at(k_int_elec_E_rebin)   .at(index).at(_util.k_plot_tot)->Fill(SC.elec_e, weight);
            TH1D_interaction_hists.at(k_int_elec_theta)     .at(index).at(_util.k_plot_tot)->Fill(SC.elec_theta, weight);
            TH1D_interaction_hists.at(k_int_elec_phi)       .at(index).at(_util.k_plot_tot)->Fill(SC.elec_phi, weight);
            TH1D_interaction_hists.at(k_int_effective_ang)  .at(index).at(_util.k_plot_tot)->Fill(SC.true_effective_angle, weight);
            TH1D_interaction_hists.at(k_int_cosbeta)        .at(index).at(_util.k_plot_tot)->Fill(true_cosbeta, weight);

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