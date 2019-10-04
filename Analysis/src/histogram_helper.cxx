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
    
    // -------------------------------------------------------------------------
    // Flash Histograms
    h_flash_time_v.resize(k_flash_MAX);
    for (unsigned int i=0; i < h_flash_time_v.size();i++){
        h_flash_time_v.at(i) = new TH1D ( Form("h_flash_time_%s", type_prefix.at(i).c_str()) ,"", 80, 0, 20);
    }

    // -------------------------------------------------------------------------
    // Reco Vtx ZY
    h_reco_vtx_zy.resize(k_vertex_MAX);

    for (unsigned int i=0; i < h_reco_vtx_zy.size();i++){
        h_reco_vtx_zy.at(i) = new TH2D ( Form("h_reco_vtx_zy_%s", vertex_strings.at(i).c_str()) ,"", 40, 0, 1038, 40, -117, 117);
    }

    // -------------------------------------------------------------------------
    // Reco Vtx X, Y, Z
    h_reco_vtx_x.resize(k_cuts_MAX);
    h_reco_vtx_y.resize(k_cuts_MAX);
    h_reco_vtx_z.resize(k_cuts_MAX);
    
    // dEdx
    h_reco_dEdx.resize(k_cuts_MAX);

    // Largest flash Z
    h_largest_flash_z.resize(k_cuts_MAX);

    // Leading Shower Momentum
    h_reco_leading_mom.resize(k_cuts_MAX);

    // 2D distance largest flash to reco nu vertex
    h_reco_flash_to_vtx_dist.resize(k_cuts_MAX);

    // 2D distance shower vertex to reco nu vertex
    h_reco_shower_to_vtx_dist.resize(k_cuts_MAX);

    // 2D distance track vertex to reco nu vertex
    h_reco_track_to_vtx_dist.resize(k_cuts_MAX);

    // Leading Shower hits in all planes
    h_reco_leading_shower_hits_all_planes.resize(k_cuts_MAX);

    // Leading Shower hits in collection
    h_reco_leading_shower_hits_collection_plane.resize(k_cuts_MAX);

    // Leading Shower opening angle
    h_reco_leading_shower_open_angle.resize(k_cuts_MAX);

    // Secondary shower to vertex distance (for events with more than 1 shower)
    h_reco_leading_shower_open_angle.resize(k_cuts_MAX);

    // Leading Shower hits per length
    h_reco_leading_shower_hits_per_length.resize(k_cuts_MAX);

    // Longest track to leading shower length
    h_reco_longest_track_eading_shower_length.resize(k_cuts_MAX);

    // Track Containment
    h_reco_track_contained.resize(k_cuts_MAX);

    // Leading shower phi
    h_reco_leading_shower_phi.resize(k_cuts_MAX);

    // Leading shower theta
    h_reco_leading_shower_theta.resize(k_cuts_MAX);

    // Leading shower cos theta
    h_reco_leading_shower_cos_theta.resize(k_cuts_MAX);

    // Leading shower multiplicity
    h_reco_shower_multiplicity.resize(k_cuts_MAX);

    // Leading track multiplicity
    h_reco_track_multiplicity.resize(k_cuts_MAX);


    for (unsigned int i=0; i < cut_dirs.size();i++){

        // Reco Vtx X, Y, Z
        h_reco_vtx_x.at(i).resize(k_classifications_MAX);
        h_reco_vtx_y.at(i).resize(k_classifications_MAX);
        h_reco_vtx_z.at(i).resize(k_classifications_MAX);

        // dEdx
        h_reco_dEdx.at(i).resize(k_classifications_MAX);

        // Largest flash Z
        h_largest_flash_z.resize(k_classifications_MAX));

        // Leading Shower Momentum
        h_reco_leading_mom.resize(k_classifications_MAX));

        // 2D distance largest flash to reco nu vertex
        h_reco_flash_to_vtx_dist.resize(k_classifications_MAX));

        // 2D distance shower vertex to reco nu vertex
        h_reco_shower_to_vtx_dist.resize(k_classifications_MAX));

        // 2D distance track vertex to reco nu vertex
        h_reco_track_to_vtx_dist.resize(k_classifications_MAX));

        // Leading Shower hits in all planes
        h_reco_leading_shower_hits_all_planes.resize(k_classifications_MAX));

        // Leading Shower hits in collection
        h_reco_leading_shower_hits_collection_plane.resize(k_classifications_MAX));

        // Leading Shower opening angle
        h_reco_leading_shower_open_angle.resize(k_classifications_MAX));

        // Secondary shower to vertex distance (for events with more than 1 shower)
        h_reco_leading_shower_open_angle.resize(k_classifications_MAX));

        // Leading Shower hits per length
        h_reco_leading_shower_hits_per_length.resize(k_classifications_MAX));

        // Longest track to leading shower length
        h_reco_longest_track_eading_shower_length.resize(k_classifications_MAX));

        // Track Containment
        h_reco_track_contained.resize(k_classifications_MAX));

        // Leading shower phi
        h_reco_leading_shower_phi.resize(k_classifications_MAX));

        // Leading shower theta
        h_reco_leading_shower_theta.resize(k_classifications_MAX));

        // Leading shower cos theta
        h_reco_leading_shower_cos_theta.resize(k_classifications_MAX));

        // Leading shower multiplicity
        h_reco_shower_multiplicity.resize(k_classifications_MAX));

        // Leading track multiplicity
        h_reco_track_multiplicity.resize(k_classifications_MAX));


        for (unsigned int j=0; j < classification_dirs.size();j++){
           
            // Reco Vtx X, Y, Z
            h_reco_vtx_x.at(i).at(j) = new TH1D ( Form("h_reco_vtx_x_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 270);
            h_reco_vtx_y.at(i).at(j) = new TH1D ( Form("h_reco_vtx_y_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, -10, 120);
            h_reco_vtx_z.at(i).at(j) = new TH1D ( Form("h_reco_vtx_z_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 40, -10, 1050);

            // dEdx
            h_reco_dEdx.at(i).at(j) = new TH1D ( Form("h_reco_dEdx_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 40, 0, 10);

            // Largest flash Z
            h_largest_flash_z.at(i).at(j) = new TH1D ( Form("h_largest_flash_z_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 100, 0, 1040);

            // Leading Shower Momentum
            h_reco_leading_mom.at(i).at(j) = new TH1D ( Form("h_reco_leading_mom_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 2);

            // 2D distance largest flash to reco nu vertex
            h_reco_flash_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_flash_to_vtx_dist._%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 40, 0, 200);

            // 2D distance shower vertex to reco nu vertex
            h_reco_shower_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_shower_to_vtx_dist_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // 2D distance track vertex to reco nu vertex
            h_reco_track_to_vtx_dist.at(i).at(j) = new TH1D ( Form("h_reco_track_to_vtx_dist_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 25, 0, 20);

            // Leading Shower hits in all planes
            h_reco_leading_shower_hits_all_planes.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_all_planes_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 30, 0, 600);

            // Leading Shower hits in collection
            h_reco_leading_shower_hits_collection_plane.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_collection_plane_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 250);

            // Leading Shower opening angle
            h_reco_leading_shower_open_angle.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_open_angle_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 25, 0, 50);

            // Secondary shower to vertex distance (for events with more than 1 shower)
            h_reco_leading_shower_open_angle.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_open_angle_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 80);

            // Leading Shower hits per length
            h_reco_leading_shower_hits_per_length.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_hits_per_length_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 20);

            // Longest track to leading shower length
            h_reco_longest_track_eading_shower_length.at(i).at(j) = new TH1D ( Form("h_reco_longest_track_eading_shower_length_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 20, 0, 3);

            // Track Containment
            h_reco_track_contained.at(i).at(j) = new TH1D ( Form("h_reco_track_contained_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 2, 0, 2);

            // Leading shower phi
            h_reco_leading_shower_phi.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_phi_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 12, -180, 180);

            // Leading shower theta
            h_reco_leading_shower_theta.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_theta_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 12, 0, 180);

            // Leading shower cos theta
            h_reco_leading_shower_cos_theta.at(i).at(j) = new TH1D ( Form("h_reco_leading_shower_cos_theta_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 16, -1, 1);

            // Leading shower multiplicity
            h_reco_shower_multiplicity.at(i).at(j) = new TH1D ( Form("h_reco_shower_multiplicity_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 5, 1, 6);

            // Leading track multiplicity
            h_reco_track_multiplicity.at(i).at(j) = new TH1D ( Form("h_reco_track_multiplicity_%s_%s",cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str()) ,"", 5, 0, 5);
        }
    }
    // -------------------------------------------------------------------------

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
    if (tpco_id == "numu_cc_qe"  || tpco_id == "numu_cc_res" ||
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

    // Dirt
    if (tpco_id == "Dirt"){
        return  k_leg_dirt;
    }

    std::cout << "Reached end of index of classifcation, this is BAD!!!! " << tpco_id << std::endl;
    return k_classifications_MAX;
}
// -----------------------------------------------------------------------------
void histogram_helper::FillReco(int classification_index, int cut_index, const xsecAna::TPCObjectContainer &tpc_obj, int leading_shower_index, std::string type){

    // For now skip cases where there is no leading shower index, need to invesitgate why
    if (leading_shower_index < 0) return;
    
    auto const leading_shower = tpc_obj.GetParticle(leading_shower_index);
    double leading_dedx = leading_shower.PfpdEdx().at(2);//just the collection plane!
    if (type == "MC" || type == "Dirt") leading_dedx = leading_dedx * (196.979 /242.72); // Only calibrate the MC

    h_reco_vtx_x.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxX());
    h_reco_vtx_y.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxY());
    h_reco_vtx_z.at(cut_index).at(classification_index)->Fill(tpc_obj.pfpVtxZ());

    h_reco_dEdx.at(cut_index).at(classification_index)->Fill(leading_dedx); //just the collection plane!

}
// -----------------------------------------------------------------------------
void histogram_helper::FillRecoVtxZY(int classification_index, const xsecAna::TPCObjectContainer &tpc_obj){

    if (classification_index == k_nue_cc)   h_reco_vtx_zy.at(k_vtx_signal)->Fill(tpc_obj.pfpVtxZ(),tpc_obj.pfpVtxY() ); // Signal
    if (classification_index == k_leg_data) h_reco_vtx_zy.at(k_vtx_data)  ->Fill(tpc_obj.pfpVtxZ(),tpc_obj.pfpVtxY() ); // Data
    if (classification_index == k_leg_ext)  h_reco_vtx_zy.at(k_vtx_ext)   ->Fill(tpc_obj.pfpVtxZ(),tpc_obj.pfpVtxY() ); // EXT
    if (classification_index != k_leg_data && classification_index != k_leg_dirt ) h_reco_vtx_zy.at(k_vtx_mc_ext)  ->Fill(tpc_obj.pfpVtxZ(),tpc_obj.pfpVtxY() ); // MC + EXT

}
// -----------------------------------------------------------------------------
void histogram_helper::WriteReco(int type){

    f_nuexsec->cd();

    bool bool_dir;
    TDirectory *truth_dir; // e.g MC/Truth, Data/Truth, EXT/Truth
    bool break_early{false};

    // loop over the cut directories
    for (unsigned int i = 0; i < cut_dirs.size(); i++){

        // loop over the classification directories
        for (unsigned int j = 0; j < classification_dirs.size(); j++){

            // Choose which folder to fill in based on the type
            if (type == k_mc && ( j == k_leg_data || j == k_leg_ext || j == k_leg_dirt)){ 
                break;
            }
            if (type == k_data){ 
                j = k_leg_data;
                break_early = true;
            }
            if (type == k_ext){ 
                j = k_leg_ext;
                break_early = true;
            }
            if (type == k_dirt){ 
                j= k_leg_dirt;
                break_early = true;
            }

            // Get the classification directory and cd
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir ,Form("%s/%s/%s/%s", type_prefix.at(type).c_str(), "Stack", cut_dirs.at(i).c_str(), classification_dirs.at(j).c_str() ) );
            
            if (bool_dir) truth_dir->cd();

            // Now write the histogram
            h_reco_vtx_x.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_y.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_vtx_z.at(i).at(j)->Write("",TObject::kOverwrite);
            h_reco_dEdx.at(i).at(j)->Write("",TObject::kOverwrite);

            if (break_early) break;
        }

    }

    // Now fill the 2D reco vtx histograms -- had to hard code a bit since one of them is MC + EXT!!
    for (unsigned int i = 0; i < vertex_strings.size(); i++){

        if (i == k_vtx_signal || i == k_vtx_mc_ext) {
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir , "MC/Reco" );
            
            if (bool_dir) truth_dir->cd();
            h_reco_vtx_zy.at(i)->Write("",TObject::kOverwrite);
        }
        
        if (i == k_vtx_data) {
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir , "Data/Reco" );
            
            if (bool_dir) truth_dir->cd();
            h_reco_vtx_zy.at(i)->Write("",TObject::kOverwrite);
        }

        if (i == k_vtx_ext) {
            bool_dir = _utility_instance.GetDirectory(f_nuexsec, truth_dir , "EXT/Reco" );
            
            if (bool_dir) truth_dir->cd();
            h_reco_vtx_zy.at(i)->Write("",TObject::kOverwrite);
        }

    }

    
}
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------