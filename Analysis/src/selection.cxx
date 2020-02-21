#include "../include/selection.h"

namespace xsecSelection {
// -----------------------------------------------------------------------------
void selection::Initialise( const char * mc_file,
                            const char * ext_file,
                            const char * data_file,
                            const char * dirt_file,
                            const char * variation_file,
                            const std::vector<double> _config,
                            bool _slim,
                            int num_events){
    
    std::cout << "\nInitialising..." << std::endl;
    if (_slim){
        std::cout << "\033[0;32m-------------------------------" << std::endl;
        std::cout << "     Running in Slim Mode!" << std::endl;
        std::cout << "-------------------------------\033[0m" << std::endl;
        slim = _slim;
    }

    // Set the scale factors
    mc_scale_factor     = Data_POT  / MC_POT;
    dirt_scale_factor   = Data_POT  / Dirt_POT;
    intime_scale_factor = Data_trig / EXT_trig;

    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:   "   << mc_scale_factor     << "\n" <<
    "Dirt Scale factor: "   << dirt_scale_factor   << "\n" <<
    "EXT Scale factor:  "   << intime_scale_factor << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // Set the maximum number of events tp process
    if (num_events > 0 ) max_events = num_events;

    // Resize the selection cuts/histogram helper instance vectors, one instance per type e.g MC, data, ..
    _scuts.resize(_util.k_type_MAX);
    _hhelper.resize(_util.k_type_MAX);

    // Print the input files
    std::cout <<
    "MC   File Path:      " << mc_file        <<"\n" <<
    "Ext  File Path:      " << ext_file       <<"\n" <<
    "Data File Path:      " << data_file      <<"\n" <<
    "Dirt File Path:      " << dirt_file      <<"\n" <<
    "Variation File Path: " << variation_file <<"\n" <<
    std::endl;

    // Now get the files, if file isnt specified then set bool to skip
    bool_use_mc        = _util.GetFile(f_mc,        mc_file);
    bool_use_ext       = _util.GetFile(f_ext,       ext_file);
    bool_use_data      = _util.GetFile(f_data,      data_file);
    bool_use_dirt      = _util.GetFile(f_dirt,      dirt_file);
    bool_use_variation = _util.GetFile(f_variation, variation_file);

    // Configure the externally configurable cut parameters
    std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
    
    // Resize the counter vector
    counter_v.resize(_util.k_COUNTER_MAX);

    for (unsigned int t = 0; t < counter_v.size(); t++){
        counter_v.at(t).resize(_util.k_COUNTER_MAX, 0);
    }

    // Get MC variables --------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nInitialising MC" << std::endl;
        _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the mc slice container
        mc_SC.Initialise(mc_tree);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_mc).Initialise(_util.k_mc);
        if (!_slim) _hhelper.at(_util.k_mc).InitHistograms();
        
        mc_tree_total_entries = mc_tree->GetEntries();
        std::cout << "Total MC Events:         " << mc_tree_total_entries << std::endl;

        // Resize the Passed vector
        mc_passed_v.resize(mc_tree_total_entries);

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of MC Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;
    } // End getting MC variables

    // Initialise Data specific ------------------------------------------------
    if (bool_use_data){
        std::cout << "\nInitialising Data" << std::endl;
        _util.GetTree(f_data, data_tree, "nuselection/NeutrinoSelectionFilter");
        
        // Initialise all the data slice container
        data_SC.Initialise(data_tree);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_data).Initialise(_util.k_data);
        if (!_slim) _hhelper.at(_util.k_data).InitHistograms();
        
        data_tree_total_entries = data_tree->GetEntries();
        std::cout << "Total Data Events:         " << data_tree_total_entries << std::endl;

        // Resize the Passed vector
        data_passed_v.resize(data_tree_total_entries);

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of Data Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of Data variables

    // Initialise EXT specific -------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nInitialising EXT" << std::endl;

        _util.GetTree(f_ext, ext_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        ext_SC.Initialise(ext_tree);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_ext).Initialise(_util.k_ext);
        if (!_slim) _hhelper.at(_util.k_ext).InitHistograms();
        
        ext_tree_total_entries = ext_tree->GetEntries();
        std::cout << "Total MC Events:         " << ext_tree_total_entries << std::endl;

        // Resize the Passed vector
        ext_passed_v.resize(ext_tree_total_entries);

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of EXT Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of ext variables

    // Initialise Dirt specific ------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nInitialising Dirt" << std::endl;

        _util.GetTree(f_dirt, dirt_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        dirt_SC.Initialise(dirt_tree);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_dirt).Initialise(_util.k_dirt);
        if (!_slim) _hhelper.at(_util.k_dirt).InitHistograms();
        
        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total MC Events:         " << dirt_tree_total_entries << std::endl;

        // Resize the Passed vector
        dirt_passed_v.resize(dirt_tree_total_entries);

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of Dirt Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of dirt variables
    
    // Invoke main selection function
    make_selection();

} // END Initialise function
// -----------------------------------------------------------------------------
// Main function for selection
void selection::make_selection(){
    std::cout << "Now Running the selection!"<< std::endl;
    
    // MC ----------------------------------------------------------------------
    if (bool_use_mc){
        std::cout << "Starting Selection over MC" << std::endl;

        // Event loop
        for (int ievent = 0; ievent < mc_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            mc_tree->GetEntry(ievent); // TPC Objects

            // Classify the event
            std::pair<std::string, int> classification = mc_SC.SliceClassifier(_util.k_mc);      // Classification of the event
            std::string interaction                    = mc_SC.SliceInteractionType(_util.k_mc); // Genie interaction type
            std::string category                       = mc_SC.SliceCategory();                  // The pandora group slice category

            // Tabulate the selection
            _util.Tabulate(interaction, classification.first, _util.k_mc, counter_v.at(_util.k_unselected) );

            // Fill the reco histograms
            if (!slim) _hhelper.at(_util.k_mc).FillReco(classification.second, _util.k_unselected, mc_SC);

            // std::cout << "Interaction: " <<  interaction << "   classification: " << classification << "   category: " << category<< std::endl;
            // if (mc_SC.slpdg > 0) std::cout << "run: " << mc_SC.run << "  subrun: " << mc_SC.sub << "  event: " << mc_SC.evt << std::endl;
            // if (mc_SC.slpdg > 0) std::cout << "slpdg: " << mc_SC.slpdg << "  topo score: " << mc_SC.topological_score << std::endl;
            // if (mc_SC.slpdg > 0) std::cout << "Category: " << mc_SC.category << "   interaction: "<< mc_SC.interaction << "   Purity: " << mc_SC.nu_purity_from_pfp <<   std::endl;
            // if (mc_SC.slpdg > 0) std::cout << "ccnc: "<< mc_SC.ccnc << std::endl;
            // if (mc_SC.slpdg < 0) std::cout << "Classification: " << classification  << "  Category: " << category2 << std::endl;
            // if (mc_SC.slpdg < 0) std::cout << "Interaction: " <<  mc_SC.SliceInteractionType(_util.k_mc) << std::endl;

            std::cout << "Shower daughter vector size: " << mc_SC.pfp_shr_daughters_v->size() << "   num showers: " << mc_SC.n_showers << std::endl;
            if (mc_SC.slpdg > 0){
                for (unsigned int k=0 ; k < mc_SC.trk_score_v->size();k++){
                    if (mc_SC.trk_score_v->at(k) < 0.5) std::cout << mc_SC.trk_score_v->at(k) << std::endl;
                }

            }

            

            // Apply the selection cuts 
            bool pass = ApplyCuts(_util.k_mc, ievent, counter_v, mc_passed_v, mc_SC, classification.first, interaction);
            if (!pass) continue;

        } // End Event loop

        std::cout << std::endl;
        std::cout << "Ending Selection over MC" << std::endl;

        // Get the total number of in cryostat nues
        tot_true_cryo_nues = counter_v.at(_util.k_unselected).at(_util.k_count_total_nue_cc_qe)  + 
                             counter_v.at(_util.k_unselected).at(_util.k_count_total_nue_cc_res) + 
                             counter_v.at(_util.k_unselected).at(_util.k_count_total_nue_cc_dis) + 
                             counter_v.at(_util.k_unselected).at(_util.k_count_total_nue_cc_coh) + 
                             counter_v.at(_util.k_unselected).at(_util.k_count_total_nue_cc_mec);

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Total Nue's in the Cryostat: " << tot_true_cryo_nues << std::endl;
        std::cout << "-------------------------------" << std::endl;
    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){
        std::cout << "Starting Selection over Data" << std::endl;

        for (int ievent = 0; ievent < data_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            data_tree->GetEntry(ievent); // TPC Objects

            // Classify the event
            std::pair<std::string, int> classification = data_SC.SliceClassifier(_util.k_data);      // Classification of the event
            std::string interaction                    = data_SC.SliceInteractionType(_util.k_data); // Genie interaction type
            std::string category                       = data_SC.SliceCategory();                    // The pandora group slice category

            // Tabulate the selection
            _util.Tabulate(interaction, classification.first, _util.k_data, counter_v.at(_util.k_unselected) );

            // Fill reconstructed histograms
            if (!slim) _hhelper.at(_util.k_data).FillReco(classification.second, _util.k_unselected, data_SC);
            // std::cout << "Interaction: " <<  interaction << "   classification: " << classification << "   category: " << category<< std::endl;
            // if (data_SC.slpdg > 0) std::cout << "run: " << data_SC.run << "  subrun: " << data_SC.sub << "  event: " << data_SC.evt << std::endl;
            // if (data_SC.slpdg > 0) std::cout << "slpdg: " << data_SC.slpdg << "  topo score: " << data_SC.topological_score << std::endl;
            // if (data_SC.slpdg > 0) std::cout << "Category: " << data_SC.category << std::endl;
        }

        std::cout << std::endl;
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "Starting Selection over EXT" << std::endl;
         
        std::cout << std::endl;
        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){
        std::cout << "Starting Selection over Dirt" << std::endl;
         
        std::cout << std::endl;
        std::cout << "Ending Selection over Dirt" << std::endl;

    }
    // -------------------------------------------------------------------------
    std::cout << "Finished running the selection!"<< std::endl;

    // Print information from the selection
    _util.PrintInfo(counter_v.at(_util.k_unselected), intime_scale_factor, mc_scale_factor, dirt_scale_factor, _util.cut_dirs.at(_util.k_unselected), counter_v.at(_util.k_unselected).at(_util.k_count_nue_cc));

    // Now save all the outputs to file
    if (!slim) SavetoFile();

    return;
} // End Selection
// -----------------------------------------------------------------------------
bool selection::ApplyCuts(int type, int ievent,std::vector<std::vector<int>> &counter_v,
                           std::vector<Passed_Container> &passed_v, SliceContainer SC,
                           std::string classification, std::string interaction){

    // Here we apply the selection cuts ----------------------------------------
    bool pass; // A flag to see if an event passes an event


    // *************************************************************************
    // Pandora Output ----------------------------------------------------------
    // *************************************************************************

    // *************************************************************************
    return true;

}
// -----------------------------------------------------------------------------
void selection::SavetoFile(){

    // Now saving histograms to file
    std::cout << "Now Saving Histograms to file" << std::endl;
    if (bool_use_mc) {
        // _hhelper.WriteMCTruth("MC");
        // _hhelper.WriteOptical(_util.k_mc);
        _hhelper.at(_util.k_mc).WriteReco(_util.k_mc);
    }
    if (bool_use_data) {
        // _hhelper.WriteOptical(_util.k_data);
        _hhelper.at(_util.k_data).WriteReco(_util.k_data);
    }
    if (bool_use_ext) {
        // _hhelper.WriteOptical(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteReco(_util.k_ext);
    }
    if (bool_use_dirt) {
        // _hhelper.WriteMCTruth("Dirt");
        // _hhelper.WriteOptical(_util.k_dirt);
        _hhelper.at(_util.k_data).WriteReco(_util.k_dirt);
    }

} // End save to file
// -----------------------------------------------------------------------------
void selection::MakeHistograms(){
    std::cout << "Creating histograms and making plots" << std::endl;

    // histogram_plotter histogram_plotter_instance;

    // histogram_plotter_instance.Initalise();

    // // Loop over the cuts and plot histograms by plot type
    // for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){
        
    //     // Create a set of strings for creating a dynamic directory
    //     // Directory structure that is created will take the form plots/<cut>/
    //     std::string a = "if [ ! -d \"plots/";
    //     std::string b = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
    //     std::string c = "; fi";
    //     std::string command = a +_util.cut_dirs.at(i) + b + _util.cut_dirs.at(i) + c ;
    //     system(command.c_str()); 

    //     // Reco X
    //     histogram_plotter_instance.MakeStack("h_reco_vtx_x",_util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Reco Vertex X [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_vtx_x.pdf", _util.cut_dirs.at(i).c_str()) );
        
    //     // Reco Y
    //     histogram_plotter_instance.MakeStack("h_reco_vtx_y",_util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Reco Vertex Y [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_vtx_y.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Reco Z
    //     histogram_plotter_instance.MakeStack("h_reco_vtx_z",_util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Reco Vertex Z [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_vtx_z.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // dEdx
    //     histogram_plotter_instance.MakeStack("h_reco_dEdx",_util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Collection Plane dEdx [MeV/cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_dEdx_collection.pdf", _util.cut_dirs.at(i).c_str()) );
   
    //     // Leading Shower Momentum
    //     histogram_plotter_instance.MakeStack("h_reco_leading_mom", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Momentum [MeV/c]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_mom.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // 2D distance largest flash to reco nu vertex
    //     histogram_plotter_instance.MakeStack("h_reco_flash_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Flash to Vertex Distance [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_flash_to_vtx_dist.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // 2D distance shower vertex to reco nu vertex
    //     histogram_plotter_instance.MakeStack("h_reco_shower_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Shower to Vertex Distance [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_shower_to_vtx_dist.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // 2D distance track vertex to reco nu vertex
    //     histogram_plotter_instance.MakeStack("h_reco_trac_util.k_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Track to Vertex Distance [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_trac_util.k_to_vtx_dist.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading Shower hits in all planes
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_all_planes", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Hits All Planes", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_hits_all_planes.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading Shower hits in collection
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_collection_plane", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Hits Collection Plane", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_hits_collection_plane.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading Shower opening angle
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_open_angle", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Open Angle [degrees]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_open_angle.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Secondary shower to vertex distance (for events with more than 1 shower)
    //     histogram_plotter_instance.MakeStack("h_reco_secondary_shower_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Secondary Shower to Vertex Distance (>1 shower) [cm]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_secondary_shower_to_vtx_dist.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading Shower hits per length
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_per_length", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Hits / Length [cm^{-1}]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_hits_per_length.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Longest track to leading shower length
    //     histogram_plotter_instance.MakeStack("h_reco_longest_trac_util.k_leading_shower_length", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Longest Track Length / Leading Shower Length", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_longest_trac_util.k_leading_shower_length.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Track Containment
    //     histogram_plotter_instance.MakeStack("h_reco_trac_util.k_contained", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Contained Tracks", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_trac_util.k_contained.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading shower phi
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_phi", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Phi [degrees]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_phi.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading shower theta
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_theta", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Theta [degrees]", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_theta.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading shower cos theta
    //     histogram_plotter_instance.MakeStack("h_reco_leading_shower_cos_theta", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Leading Shower Cos(#theta)", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_leading_shower_cos_theta.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading shower multiplicity
    //     histogram_plotter_instance.MakeStack("h_reco_shower_multiplicity", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Shower Multiplicty", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_shower_multiplicity.pdf", _util.cut_dirs.at(i).c_str()) );

    //     // Leading track multiplicity
    //     histogram_plotter_instance.MakeStack("h_reco_trac_util.k_multiplicity", _util.cut_dirs.at(i).c_str(),
    //                                        false,  false, "Track Multiplicty", mc_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
    //                                        Form("plots/%s/reco_trac_util.k_multiplicity.pdf", _util.cut_dirs.at(i).c_str()) );
   
    // }
    

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection