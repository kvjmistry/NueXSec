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
void histogram_helper::Initialise(int type, const char* run_period ){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;

    if (type == _util.k_mc){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject( Form("files/nuexsec_mc_run%s.root", run_period) ) ) {
            f_nuexsec = new TFile(Form("files/nuexsec_mc_run%s.root", run_period), "UPDATE");
        }
    }
    else if (type == _util.k_data){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(Form("files/nuexsec_data_run%s.root", run_period)) ) {
            f_nuexsec = new TFile(Form("files/nuexsec_data_run%s.root", run_period), "UPDATE");
        }

    }
    else if (type == _util.k_ext){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(Form("files/nuexsec_ext_run%s.root", run_period)) ) {
            f_nuexsec = new TFile(Form("files/nuexsec_ext_run%s.root", run_period), "UPDATE");
        }

    }
    else if (type == _util.k_dirt){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject(Form("files/nuexsec_dirt_run%s.root", run_period)) ) {
            f_nuexsec = new TFile(Form("files/nuexsec_dirt_run%s.root", run_period), "UPDATE");
        }

    }
    else {
        std::cout << "Unknown input type!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
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

            // dEdx
            TH1D_hists.at(k_reco_dEdx).at(i).at(j) = new TH1D ( Form("h_reco_dEdx_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            // Leading Shower Momentum
            TH1D_hists.at(k_reco_leading_mom).at(i).at(j) = new TH1D ( Form("h_reco_leading_mom_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 2);

            // 2D distance shower vertex to reco nu vertex
            TH1D_hists.at(k_reco_shower_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // 2D distance track vertex to reco nu vertex
            TH1D_hists.at(k_reco_track_to_vtx_dist).at(i).at(j) = new TH1D ( Form("h_reco_track_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 20);

            // Leading Shower hits in all planes
            TH1D_hists.at(k_reco_leading_shower_hits_all_planes).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_all_planes_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 600);

            // Leading Shower hits in collection
            TH1D_hists.at(k_reco_leading_shower_hits_collection_plane).at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_collection_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 250);

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

        }
        
    }
    // -------------------------------------------------------------------------

}
// -----------------------------------------------------------------------------
void histogram_helper::FillReco(int classification_index, int cut_index, SliceContainer &SC){


    // Now fill the histograms!
    TH1D_hists.at(k_reco_vtx_x).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x);
    TH1D_hists.at(k_reco_vtx_y).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y);
    TH1D_hists.at(k_reco_vtx_z).at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z);

    TH1D_hists.at(k_reco_dEdx).at(cut_index).at(classification_index)->Fill(SC.shr_dedx_Y_cali); // Just the collection plane!
    
    // TH1D_hists.at(k_reco_leading_mom).at(cut_index).at(classification_index)->Fill();
    
    TH1D_hists.at(k_reco_shower_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.shr_distance);

    TH1D_hists.at(k_reco_track_to_vtx_dist).at(cut_index).at(classification_index)->Fill(SC.trk_distance);
    
    // TH1D_hists.at(k_reco_leading_shower_hits_all_planes).at(cut_index).at(classification_index)->Fill();
    
    TH1D_hists.at(k_reco_leading_shower_hits_collection_plane).at(cut_index).at(classification_index)->Fill(SC.shr_hits_y_tot);

    // TH1D_hists.at(k_reco_secondary_shower_to_vtx_dist).at(cut_index).at(classification_index)->Fill();
    
    TH1D_hists.at(k_reco_leading_shower_open_angle).at(cut_index).at(classification_index)->Fill(SC.shr_openangle * 180/3.14159);
    
    // TH1D_hists.at(k_reco_leading_shower_hits_per_length).at(cut_index).at(classification_index)->Fill();
    
    // TH1D_hists.at(k_reco_longest_track_leading_shower_length).at(cut_index).at(classification_index)->Fill();
    
    // TH1D_hists.at(k_reco_track_contained).at(cut_index).at(classification_index)->Fill(); // We fill this once per tpc obj (only 1 track has to pass)
    
    TH1D_hists.at(k_reco_leading_shower_phi).at(cut_index).at(classification_index)->Fill(SC.shr_phi * 180/3.14159);
    
    TH1D_hists.at(k_reco_leading_shower_theta).at(cut_index).at(classification_index)->Fill(SC.shr_theta * 180/3.14159);
    
    TH1D_hists.at(k_reco_leading_shower_cos_theta).at(cut_index).at(classification_index)->Fill(std::cos(SC.shr_theta));
    
    TH1D_hists.at(k_reco_shower_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_showers);
    
    TH1D_hists.at(k_reco_track_multiplicity).at(cut_index).at(classification_index)->Fill(SC.n_tracks);

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
                TH1D_hists.at(u).at(i).at(j)->Write("",TObject::kOverwrite);

                if (break_early) break;
            }

        }
    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------