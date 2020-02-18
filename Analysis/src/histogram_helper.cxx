#include "../include/histogram_helper.h"

// -----------------------------------------------------------------------------
histogram_helper::~histogram_helper() { 
    
    // Make sure the file is closed
    // f_nuexsec->Close();
}
// -----------------------------------------------------------------------------
void histogram_helper::MakeDirectory(std::string type){
        
    f_nuexsec->cd();

    TDirectory *top_dir; // e.g MC, Data, EXT
    bool bool_dir;       // Check if directory exists already
    TString type_tstr = type;
   
    // Create the top directory
    bool_dir = _util.GetDirectory(f_nuexsec, top_dir, type_tstr );
    if (!bool_dir) top_dir = f_nuexsec->mkdir(type.c_str());
    
    // Make the the top dir the current directory
    top_dir->cd();

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
        bool_dir = _util.GetDirectory(f_nuexsec, dir_plot_types[k] ,Form("%s/%s", type.c_str(), _util.plot_types.at(k).c_str()) );

        // Make the directory
        if (!bool_dir) dir_plot_types[k] = top_dir->mkdir(_util.plot_types.at(k).c_str());

        dir_plot_types[k]->cd();

        // If we have stacked histograms, we make plots by cut 
        if (_util.plot_types.at(k) == "Stack"){
            
            // Loop over the cuts ----------------------------------------------
            for (int i = 0; i < ncuts; i++) {
               
                // Get the directory 
                bool_dir = _util.GetDirectory(f_nuexsec, dir_cut[i] ,Form("%s/%s/%s", type.c_str(), _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str()));

                // Make the directory
                if (!bool_dir) dir_cut[i] = dir_plot_types[k]->mkdir(_util.cut_dirs.at(i).c_str());
                dir_cut[i]->cd();
                
                // Loop over the classifications -------------------------------
                for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){
                
                    // Get the directory 
                    bool_dir = _util.GetDirectory(f_nuexsec, dir_classification[j] ,Form("%s/%s/%s/%s", type.c_str(), _util.plot_types.at(k).c_str(), _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()));

                    // Make the directory
                    if (!bool_dir) dir_classification[j] = dir_cut[i]->mkdir(_util.classification_dirs.at(j).c_str());
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
void histogram_helper::Initialise(int type){

    std::cout << "Initalising Histogram Helper, creating TFile and directories..." << std::endl;


    if (type == _util.k_mc){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject("files/nuexsec_mc.root") ) {
            f_nuexsec = new TFile("files/nuexsec_mc.root", "UPDATE");
        }
    }
    else if (type == _util.k_data){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject("files/nuexsec_data.root") ) {
            f_nuexsec = new TFile("files/nuexsec_data.root", "UPDATE");
        }

    }
    else if (type == _util.k_ext){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject("files/nuexsec_ext.root") ) {
            f_nuexsec = new TFile("files/nuexsec_ext.root", "UPDATE");
        }

    }
    else if (type == _util.k_dirt){
        
        // File not already open, open the file
        if (!gROOT->GetListOfFiles()->FindObject("files/nuexsec_dirt.root") ) {
            f_nuexsec = new TFile("files/nuexsec_dirt.root", "UPDATE");
        }

    }
    else {
        std::cout << "Unknown input type!! "<<  __PRETTY_FUNCTION__ << std::endl;
        exit(1);
    }



    MakeDirectory("MC");
    MakeDirectory("Data");
    MakeDirectory("EXT");
    MakeDirectory("Dirt");
}
// -----------------------------------------------------------------------------
void histogram_helper::InitHistograms(){
    
    // -------------------------------------------------------------------------
    // Reco Vtx X, Y, Z
    h_reco_vtx_x.resize(_util.k_cuts_MAX);
    h_reco_vtx_y.resize(_util.k_cuts_MAX);
    h_reco_vtx_z.resize(_util.k_cuts_MAX);
    
    // dEdx
    h_reco_dEdx.resize(_util.k_cuts_MAX);

    // Leading Shower Momentum
    h_reco_leading_mom.resize(_util.k_cuts_MAX);

    // 2D distance shower vertex to reco nu vertex
    h_reco_shower_to_vtx_dist.resize(_util.k_cuts_MAX);

    // 2D distance track vertex to reco nu vertex
    h_reco_track_to_vtx_dist.resize(_util.k_cuts_MAX);

    // Leading Shower hits in all planes
    h_reco_leading_shower_hits_all_planes.resize(_util.k_cuts_MAX);

    // Leading Shower hits in collection
    h_reco_leading_shower_hits_collection_plane.resize(_util.k_cuts_MAX);

    // Leading Shower opening angle
    h_reco_leading_shower_open_angle.resize(_util.k_cuts_MAX);

    // Secondary shower to vertex distance (for events with more than 1 shower)
    h_reco_secondary_shower_to_vtx_dist.resize(_util.k_cuts_MAX);

    // Leading Shower hits per length
    h_reco_leading_shower_hits_per_length.resize(_util.k_cuts_MAX);

    // Longest track to leading shower length
    h_reco_longest_track_leading_shower_length.resize(_util.k_cuts_MAX);

    // Track Containment
    h_reco_track_contained.resize(_util.k_cuts_MAX);

    // Leading shower phi
    h_reco_leading_shower_phi.resize(_util.k_cuts_MAX);

    // Leading shower theta
    h_reco_leading_shower_theta.resize(_util.k_cuts_MAX);

    // Leading shower cos theta
    h_reco_leading_shower_cos_theta.resize(_util.k_cuts_MAX);

    // Leading shower multiplicity
    h_reco_shower_multiplicity.resize(_util.k_cuts_MAX);

    // Leading track multiplicity
    h_reco_track_multiplicity.resize(_util.k_cuts_MAX);


    for (unsigned int i=0; i < _util.cut_dirs.size();i++){

        // Reco Vtx X, Y, Z
        h_reco_vtx_x.at(i).resize(_util.k_classifications_MAX);
        h_reco_vtx_y.at(i).resize(_util.k_classifications_MAX);
        h_reco_vtx_z.at(i).resize(_util.k_classifications_MAX);

        // dEdx
        h_reco_dEdx.at(i).resize(_util.k_classifications_MAX);

        // Leading Shower Momentum
        h_reco_leading_mom.at(i).resize(_util.k_classifications_MAX);

        // 2D distance shower vertex to reco nu vertex
        h_reco_shower_to_vtx_dist.at(i).resize(_util.k_classifications_MAX);

        // 2D distance track vertex to reco nu vertex
        h_reco_track_to_vtx_dist.at(i).resize(_util.k_classifications_MAX);

        // Leading Shower hits in all planes
        h_reco_leading_shower_hits_all_planes.at(i).resize(_util.k_classifications_MAX);

        // Leading Shower hits in collection
        h_reco_leading_shower_hits_collection_plane.at(i).resize(_util.k_classifications_MAX);

        // Leading Shower opening angle
        h_reco_leading_shower_open_angle.at(i).resize(_util.k_classifications_MAX);

        // Secondary shower to vertex distance (for events with more than 1 shower)
        h_reco_secondary_shower_to_vtx_dist.at(i).resize(_util.k_classifications_MAX);

        // Leading Shower hits per length
        h_reco_leading_shower_hits_per_length.at(i).resize(_util.k_classifications_MAX);

        // Longest track to leading shower length
        h_reco_longest_track_leading_shower_length.at(i).resize(_util.k_classifications_MAX);

        // Track Containment
        h_reco_track_contained.at(i).resize(_util.k_classifications_MAX);

        // Leading shower phi
        h_reco_leading_shower_phi.at(i).resize(_util.k_classifications_MAX);

        // Leading shower theta
        h_reco_leading_shower_theta.at(i).resize(_util.k_classifications_MAX);

        // Leading shower cos theta
        h_reco_leading_shower_cos_theta.at(i).resize(_util.k_classifications_MAX);

        // Leading shower multiplicity
        h_reco_shower_multiplicity.at(i).resize(_util.k_classifications_MAX);

        // Leading track multiplicity
        h_reco_track_multiplicity.at(i).resize(_util.k_classifications_MAX);


        for (unsigned int j=0; j < _util.classification_dirs.size();j++){
           
            // Reco Vtx X, Y, Z
            h_reco_vtx_x.at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 270);
            h_reco_vtx_y.at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            h_reco_vtx_z.at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            // dEdx
            h_reco_dEdx.at(i).at(j) = new TH1D ( Form("h_reco_dEdx_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            // Leading Shower Momentum
            h_reco_leading_mom.at(i).at(j) = new TH1D ( Form("h_reco_leading_mom_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 2);

            // 2D distance shower vertex to reco nu vertex
            h_reco_shower_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // 2D distance track vertex to reco nu vertex
            h_reco_track_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_track_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 20);

            // Leading Shower hits in all planes
            h_reco_leading_shower_hits_all_planes.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_all_planes_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 30, 0, 600);

            // Leading Shower hits in collection
            h_reco_leading_shower_hits_collection_plane.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_collection_plane_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 250);

            // Leading Shower opening angle
            h_reco_leading_shower_open_angle.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_open_angle_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 25, 0, 50);

            // Secondary shower to vertex distance (for events with more than 1 shower)
            h_reco_secondary_shower_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_secondary_shower_to_vtx_dist_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 80);

            // Leading Shower hits per length
            h_reco_leading_shower_hits_per_length.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_per_length_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // Longest track to leading shower length
            h_reco_longest_track_leading_shower_length.at(i).at(j) = new TH1D ( Form("h_reco_longest_track_leading_shower_length_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 20, 0, 3);

            // Track Containment
            h_reco_track_contained.at(i).at(j) = new TH1D ( Form("h_reco_track_contained_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 2, 0, 2);

            // Leading shower phi
            h_reco_leading_shower_phi.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_phi_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 12, -180, 180);

            // Leading shower theta
            h_reco_leading_shower_theta.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 12, 0, 180);

            // Leading shower cos theta
            h_reco_leading_shower_cos_theta.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_cos_theta_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 16, -1, 1);

            // Leading shower multiplicity
            h_reco_shower_multiplicity.at(i).at(j) = new TH1D ( Form("h_reco_shower_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 1, 6);

            // Leading track multiplicity
            h_reco_track_multiplicity.at(i).at(j) = new TH1D ( Form("h_reco_track_multiplicity_%s_%s",_util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str()) ,"", 5, 0, 5);
        }
    }
    // -------------------------------------------------------------------------

}
// -----------------------------------------------------------------------------
void histogram_helper::FillReco(int classification_index, int cut_index, SliceContainer &SC){


    // Now fill the histograms!
    h_reco_vtx_x.at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_x);
    h_reco_vtx_y.at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_y);
    h_reco_vtx_z.at(cut_index).at(classification_index)->Fill(SC.reco_nu_vtx_z);

    h_reco_dEdx.at(cut_index).at(classification_index)->Fill(SC.shr_dedx_Y_cali); // Just the collection plane!
    
    // h_reco_leading_mom.at(cut_index).at(classification_index)->Fill();
    
    h_reco_shower_to_vtx_dist.at(cut_index).at(classification_index)->Fill(SC.shr_distance);

    h_reco_track_to_vtx_dist.at(cut_index).at(classification_index)->Fill(SC.trk_distance);
    
    // h_reco_leading_shower_hits_all_planes.at(cut_index).at(classification_index)->Fill();
    
    h_reco_leading_shower_hits_collection_plane.at(cut_index).at(classification_index)->Fill(SC.shr_hits_y_tot);

    // h_reco_secondary_shower_to_vtx_dist.at(cut_index).at(classification_index)->Fill();
    
    h_reco_leading_shower_open_angle.at(cut_index).at(classification_index)->Fill(SC.shr_openangle * 180/3.14159);
    
    // h_reco_leading_shower_hits_per_length.at(cut_index).at(classification_index)->Fill();
    
    // h_reco_longest_track_leading_shower_length.at(cut_index).at(classification_index)->Fill();
    
    // h_reco_track_contained.at(cut_index).at(classification_index)->Fill(); // We fill this once per tpc obj (only 1 track has to pass)
    
    h_reco_leading_shower_phi.at(cut_index).at(classification_index)->Fill(SC.shr_phi * 180/3.14159);
    
    h_reco_leading_shower_theta.at(cut_index).at(classification_index)->Fill(SC.shr_theta * 180/3.14159);
    
    h_reco_leading_shower_cos_theta.at(cut_index).at(classification_index)->Fill(std::cos(SC.shr_theta));
    
    h_reco_shower_multiplicity.at(cut_index).at(classification_index)->Fill(SC.n_showers);
    
    h_reco_track_multiplicity.at(cut_index).at(classification_index)->Fill(SC.n_tracks);

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteReco(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth
    bool break_early{false};

    // loop over the cut directories
    for (unsigned int i = 0; i < _util.cut_dirs.size(); i++){

        // loop over the classification directories
        for (unsigned int j = 0; j < _util.classification_dirs.size(); j++){

            // Choose which folder to fill in based on the type
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
            bool_dir = _util.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s/%s", _util.type_prefix.at(type).c_str(), "Stack", _util.cut_dirs.at(i).c_str(), _util.classification_dirs.at(j).c_str() ) );
            
            if (bool_dir) truth_dir->cd();

            // Now write the histograms
            h_reco_vtx_x.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_y.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_z.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_dEdx .at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_mom.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_shower_to_vtx_dist.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_track_to_vtx_dist.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_hits_all_planes.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_hits_collection_plane.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_open_angle.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_secondary_shower_to_vtx_dist.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_hits_per_length.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_longest_track_leading_shower_length.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_track_contained.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_phi.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_theta.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_leading_shower_cos_theta.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_shower_multiplicity.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_track_multiplicity.at(i).at(j)->Write("",TObject::kOverwrite);

            if (break_early) break;
        }

    }
}
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------