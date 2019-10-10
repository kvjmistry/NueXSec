#include "../include/selection.h"
#include <omp.h>

namespace xsecSelection {
// -----------------------------------------------------------------------------
void selection::Initialise( const char * mc_file,
                            const char * ext_file,
                            const char * data_file,
                            const char * dirt_file,
                            const char * variation_file,
                            const std::vector<double> _config,
                            bool _slim){
    
    std::cout << "\nInitialising..." << std::endl;
    if (_slim){
        std::cout << "\033[0;32m-------------------------------" << std::endl;
        std::cout << "     Running in Slim Mode!" << std::endl;
        std::cout << "-------------------------------\033[0m" << std::endl;
        slim = _slim;
    }

    // Initialise the histogram helper
    if (!_slim) histogram_helper_instance.Initialise();
    if (!_slim) histogram_helper_instance.InitHistograms();

    // Resize the selection cuts instance vector, one instance per type e.g MC, data, ..
    selection_cuts_instance.resize(histogram_helper::k_type_MAX);

    // Print the input files
    std::cout <<
    "MC   File Path:      " << mc_file        <<"\n" <<
    "Ext  File Path:      " << ext_file       <<"\n" <<
    "Data File Path:      " << data_file      <<"\n" <<
    "Dirt File Path:      " << dirt_file      <<"\n" <<
    "Variation File Path: " << variation_file <<"\n" <<
    std::endl;

    // Now get the files, if file isnt specified then set bool to skip
    bool_use_mc        = _utility_instance.GetFile(f_mc,        mc_file);
    bool_use_ext       = _utility_instance.GetFile(f_ext,       ext_file);
    bool_use_data      = _utility_instance.GetFile(f_data,      data_file);
    bool_use_dirt      = _utility_instance.GetFile(f_dirt,      dirt_file);
    bool_use_variation = _utility_instance.GetFile(f_variation, variation_file);

    // Configure the externally configurable cut parameters
    std::cout << "\n --- Configuring Parameters --- \n" << std::endl;
    _x1                = _config[0];
    _x2                = _config[1];
    _y1                = _config[2];
    _y2                = _config[3];
    _z1                = _config[4];
    _z2                = _config[5];
    flash_pe_threshold = _config[6];
    flash_time_start   = _config[7];
    flash_time_end     = _config[8];
    tolerance          = _config[9];
    shwr_nue_tolerance = _config[10];
    trk_nue_tolerance  = _config[11];
    shwr_hit_threshold = _config[12];
    shwr_hit_threshold_collection = _config[13];
    tolerance_open_angle_min      = _config[14];
    tolerance_open_angle_max      = _config[15];
    tolerance_dedx_min            = _config[16];
    tolerance_dedx_max            = _config[17];
    dist_tolerance                = _config[18];
    pfp_hits_length_tolerance     = _config[19];
    ratio_tolerance               = _config[20];
    detector_variations           = _config[21];

    std::cout <<
    "_x1                          : " << _config[0] << "\n" << 
    "_x2                          : " << _config[1] << "\n" << 
    "_y1                          : " << _config[2] << "\n" << 
    "_y2                          : " << _config[3] << "\n" << 
    "_z1                          : " << _config[4] << "\n" << 
    "_z2                          : " << _config[5] << "\n" << 
    "flash_pe_threshold           : " << _config[6] << "\n" << 
    "flash_time_start             : " << _config[7] << "\n" << 
    "flash_time_end               : " << _config[8] << "\n" << 
    "tolerance                    : " << _config[9] << "\n" << 
    "shwr_nue_tolerance           : " << _config[10]<< "\n" << 
    "trk_nue_tolerance            : " << _config[11]<< "\n" << 
    "shwr_hit_threshold           : " << _config[12]<< "\n" << 
    "shwr_hit_threshold_collection: " << _config[13] << "\n" <<
    "tolerance_open_angle_min     : " << _config[14] << "\n" <<
    "tolerance_open_angle_max     : " << _config[15] << "\n" <<
    "tolerance_dedx_min           : " << _config[16] << "\n" <<
    "tolerance_dedx_max           : " << _config[17] << "\n" <<
    "dist_tolerance               : " << _config[18] << "\n" <<
    "pfp_hits_length_tolerance    : " << _config[19] << "\n" <<
    "ratio_tolerance              : " << _config[20] << "\n" <<
    "detector_variations          : " << _config[21] << "\n" <<
    "\033[0;31m-------------------------------\033[0m" <<
    std::endl;

    // Set the Fiducial volume
    fv_boundary_v = {_x1, _x2, _y1, _y2, _z1, _z2};

    // Get MC variables --------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nInitialising MC" << std::endl;
        _utility_instance.GetTree(f_mc, mytree, "AnalyzeTPCO/tree");
        _utility_instance.GetTree(f_mc, optree, "AnalyzeTPCO/optical_tree");
        _utility_instance.GetTree(f_mc, mctree, "AnalyzeTPCO/mcparticle_tree");
        _utility_instance.GetTree(f_mc, mctruth_counter_tree, "AnalyzeTPCO/mctruth_counter_tree");

        // Set the branches of the truth counter tree -- function in utility?
        mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter",      &mc_nue_cc_counter);
        mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter",      &mc_nue_nc_counter);
        mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter",     &mc_numu_cc_counter);
        mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter",     &mc_numu_nc_counter);
        mctruth_counter_tree->SetBranchAddress("mc_nue_cc_counter_bar",  &mc_nue_cc_counter_bar);
        mctruth_counter_tree->SetBranchAddress("mc_numu_cc_counter_bar", &mc_numu_cc_counter_bar);
        mctruth_counter_tree->SetBranchAddress("mc_nue_nc_counter_bar",  &mc_nue_nc_counter_bar);
        mctruth_counter_tree->SetBranchAddress("mc_numu_nc_counter_bar", &mc_numu_nc_counter_bar);
        mctruth_counter_tree->SetBranchAddress("fMCNuEnergy",            &mc_nu_energy);
        mctruth_counter_tree->SetBranchAddress("fMCNuMomentum",          &mc_nu_momentum);
        mctruth_counter_tree->SetBranchAddress("fMCNuID",                &mc_nu_id);
        mctruth_counter_tree->SetBranchAddress("fMCNuVtxX",              &mc_nu_vtx_x);
        mctruth_counter_tree->SetBranchAddress("fMCNuVtxY",              &mc_nu_vtx_y);
        mctruth_counter_tree->SetBranchAddress("fMCNuVtxZ",              &mc_nu_vtx_z);
        mctruth_counter_tree->SetBranchAddress("fMCNuDirX",              &mc_nu_dir_x);
        mctruth_counter_tree->SetBranchAddress("fMCNuDirY",              &mc_nu_dir_y);
        mctruth_counter_tree->SetBranchAddress("fMCNuDirZ",              &mc_nu_dir_z);
        mctruth_counter_tree->SetBranchAddress("fMCNumParticles",        &mc_nu_num_particles);
        mctruth_counter_tree->SetBranchAddress("fMCNumChargedParticles", &mc_nu_num_charged_particles);
        mctruth_counter_tree->SetBranchAddress("fMCEleDirX",             &mc_ele_dir_x);
        mctruth_counter_tree->SetBranchAddress("fMCEleDirY",             &mc_ele_dir_y);
        mctruth_counter_tree->SetBranchAddress("fMCEleDirZ",             &mc_ele_dir_z);
        mctruth_counter_tree->SetBranchAddress("fMCEleEnergy",           &mc_ele_energy);
        mctruth_counter_tree->SetBranchAddress("fMCEleMomentum",         &mc_ele_momentum);
        mctruth_counter_tree->SetBranchAddress("has_pi0",                &has_pi0);
        mctruth_counter_tree->SetBranchAddress("fMCNuTime",              &mc_nu_time);


        // Set the branch of the TPC Object Container
        mytree->SetBranchAddress("TpcObjectContainerV", &tpc_object_container_v);

        // Now get all the MC counters for the nuetrinos
        const int total_mc_entries = mctruth_counter_tree->GetEntries();
        true_in_tpc_v.resize(total_mc_entries, false);
        std::cout << "Total MC Entries: " << total_mc_entries << std::endl;
       
        // Loop over the total MC entries
        for (int i = 0; i < total_mc_entries; i++) {
            mctruth_counter_tree->GetEntry(i);
            
            bool true_in_tpc = _utility_instance.in_fv(mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v);
            true_in_tpc_v.at(i) = true_in_tpc;
            
            // Total Nue in Cryostat
            if (mc_nu_id == 1) _mc_nue_cc_counter++; 
            if (mc_nu_id == 3) _mc_nue_nc_counter++; 
            if (mc_nu_id == 5) _mc_nue_cc_counter_bar++; 
            if (mc_nu_id == 7) _mc_nue_nc_counter_bar++; 
            
            // Total Nue in FV
            if (true_in_tpc == true && (mc_nu_id == 1)) total_mc_entries_inFV_nue++;
            if (true_in_tpc == true && (mc_nu_id == 5)) total_mc_entries_inFV_nue_bar++;
            if (true_in_tpc == true && (mc_nu_id == 3)) total_mc_entries_inFV_nue_nc++;
            if (true_in_tpc == true && (mc_nu_id == 7)) total_mc_entries_inFV_nue_nc_bar++;
            
            // Total Numu in Cryostat
            if (mc_nu_id == 2) _mc_numu_cc_counter++;
            if (mc_nu_id == 4) _mc_numu_nc_counter++;
            if (mc_nu_id == 6) _mc_numu_cc_counter_bar++;
            if (mc_nu_id == 8) _mc_numu_nc_counter_bar++;
            
            // Total Numu in FV
            if (true_in_tpc == true && (mc_nu_id == 2)) total_mc_entries_inFV_numu_cc++;
            if (true_in_tpc == true && (mc_nu_id == 4)) total_mc_entries_inFV_numu_nc++;
            if (true_in_tpc == true && (mc_nu_id == 6)) total_mc_entries_inFV_numu_cc_bar++;
            if (true_in_tpc == true && (mc_nu_id == 8)) total_mc_entries_inFV_numu_nc_bar++;
        }

        total_mc_entries_inFV = total_mc_entries_inFV_nue + total_mc_entries_inFV_nue_bar;
    
        std::cout << "MC Nue CC Counter      --- " << _mc_nue_cc_counter << std::endl;
        std::cout << "MC Nue NC Counter      --- " << _mc_nue_nc_counter << std::endl;
        std::cout << "MC Nue CC Counter Bar  --- " << _mc_nue_cc_counter_bar << std::endl;
        std::cout << "MC Nue NC Counter Bar  --- " << _mc_nue_nc_counter_bar << std::endl;
        
        std::cout << "MC Numu CC Counter     --- " << _mc_numu_cc_counter << std::endl;
        std::cout << "MC Numu NC Counter     --- " << _mc_numu_nc_counter << std::endl;
        std::cout << "MC Numu CC Counter Bar --- " << _mc_numu_cc_counter_bar << std::endl;
        std::cout << "MC Numu NC Counter Bar --- " << _mc_numu_nc_counter_bar << std::endl;
        std::cout << "-------------------------------" << std::endl;
    
        // Get the largest flash vector and optical lists for PE and flash time
        _utility_instance.GetLargestFlashVector(optree, flash_time_start, flash_time_end, flash_pe_threshold, mc_largest_flash_v_v, mc_optical_list_pe_v, mc_optical_list_flash_time_v);
        
        // Fill the flash histograms
        if (!_slim) histogram_helper_instance.FillOptical(mc_optical_list_flash_time_v, histogram_helper::k_mc);

        // Get the total number of events in TPC Obj Tree
        tree_total_entries = mytree->GetEntries();
        std::cout << "Total MC Events:           " << tree_total_entries << std::endl;

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of MC Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;
    } // End getting MC variables

    // Initialise Data specific ------------------------------------------------
    if (bool_use_data){
        std::cout << "\nInitialising Data" << std::endl;
        _utility_instance.GetTree(f_data, data_tree,   "AnalyzeTPCO/tree");
        _utility_instance.GetTree(f_data, data_optree, "AnalyzeTPCO/optical_tree");
        
        // Data TPC Object Container
        data_tree->SetBranchAddress("TpcObjectContainerV", &data_tpc_object_container_v);
 
        // Get the largest flash vector and optical lists for PE and flash time
        _utility_instance.GetLargestFlashVector(data_optree, flash_time_start, flash_time_end,flash_pe_threshold, data_largest_flash_v_v, data_optical_list_pe_v, data_optical_list_flash_time_v);

        // Fill the flash histograms
        if (!_slim) histogram_helper_instance.FillOptical(data_optical_list_flash_time_v, histogram_helper::k_data);

        // Get the total number of events in TPC Obj Tree
        data_tree_total_entries = data_tree->GetEntries();
        std::cout << "Total Data Events:         " << data_tree_total_entries << std::endl;

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of Data Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of Data variables

    // Initialise EXT specific -------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nInitialising EXT" << std::endl;
        _utility_instance.GetTree(f_ext, ext_tree,   "AnalyzeTPCO/tree");
        _utility_instance.GetTree(f_ext, ext_optree, "AnalyzeTPCO/optical_tree");
        
        // EXT TPC Object Container
        ext_tree->SetBranchAddress("TpcObjectContainerV", &ext_tpc_object_container_v);
 
        // Get the largest flash vector and optical lists for PE and flash time
        _utility_instance.GetLargestFlashVector(ext_optree, flash_time_start, flash_time_end, flash_pe_threshold, ext_largest_flash_v_v, ext_optical_list_pe_v, ext_optical_list_flash_time_v);

        // Fill the flash histograms
        if (!_slim) histogram_helper_instance.FillOptical(ext_optical_list_flash_time_v, histogram_helper::k_ext);

        // Get the total number of events in TPC Obj Tree
        ext_tree_total_entries = ext_tree->GetEntries();
        std::cout << "Total ext Events:          " << ext_tree_total_entries << std::endl;

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of EXT Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of ext variables

    // Initialise Dirt specific ------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nInitialising Dirt" << std::endl;
        _utility_instance.GetTree(f_dirt, dirt_tree,   "AnalyzeTPCO/tree");
        _utility_instance.GetTree(f_dirt, dirt_optree, "AnalyzeTPCO/optical_tree");
        
        // Dirt TPC Object Container
        dirt_tree->SetBranchAddress("TpcObjectContainerV", &dirt_tpc_object_container_v);
 
        // Get the largest flash vector and optical lists for PE and flash time
        _utility_instance.GetLargestFlashVector(dirt_optree, flash_time_start, flash_time_end, flash_pe_threshold, dirt_largest_flash_v_v, dirt_optical_list_pe_v, dirt_optical_list_flash_time_v);

        // Fill the flash histograms
        if (!_slim) histogram_helper_instance.FillOptical(dirt_optical_list_flash_time_v, histogram_helper::k_dirt);

        // Get the total number of events in TPC Obj Tree
        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total dirt Events:         " << dirt_tree_total_entries << std::endl;

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
    
    // *************************************************************************
    //  ----- Loop over the TPC Objects and Apply the selection cuts -----------
    // *************************************************************************

    // MC ----------------------------------------------------------------------
    if (bool_use_mc){
        std::cout << "Starting Selection over MC" << std::endl;

        // resize the Passed vector
        mc_passed_v.resize(tree_total_entries);
        
        // Resize the Counter Vector
        mc_counter_v.resize(Passed_Container::k_cuts_MAX);
        for (unsigned int i = 0; i < mc_counter_v.size(); i++) mc_counter_v.at(i).resize(24, 0.0 ); // Should remove hardcoded number
        
        // Loop over the Events
        for (int event = 0; event < tree_total_entries; event++){
            // Alert the user
            if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k " << std::flush;
        
            // Get the entry in the tree
            mytree->GetEntry(event);               // TPC Objects
            mctruth_counter_tree->GetEntry(event); // MC Tree

            // Fill MC Truth Histograms
            if (!slim) histogram_helper_instance.FillMCTruth( mc_nu_energy, mc_nu_momentum, mc_nu_id, true_in_tpc_v.at(event), mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, mc_nu_dir_x, mc_nu_dir_y,
                                                   mc_nu_dir_z, mc_ele_dir_x, mc_ele_dir_y, mc_ele_dir_z, mc_ele_energy, mc_ele_momentum);
            
            // The total number of TPC Objects
            int n_tpc_obj = tpc_object_container_v->size();

            // Largest flash information
            std::vector<double> largest_flash_v = mc_largest_flash_v_v.at(event); // Vec with the largest flash

            // Create an instance of the selection cut (also initialises flash info)
            selection_cuts_instance.at(histogram_helper::k_mc).SetFlashVariables(largest_flash_v);
            
            // Loop over the TPC Objects ---------------------------------------
            // (In Pandora Consolidated, there should be 1 TPC Object per event)
            for (int i = 0; i < n_tpc_obj; i++){
                const xsecAna::TPCObjectContainer tpc_obj = tpc_object_container_v->at(i); // Get the TPC Obj

                // Apply the selection cuts 
                bool pass = ApplyCuts(histogram_helper::k_mc, event, tpc_obj, mc_largest_flash_v_v, mc_optical_list_pe_v, mc_optical_list_flash_time_v, mc_counter_v, mc_passed_v);
                if (!pass) continue;

            } // End loop over the TPC Objects

        
        } // End loop over the Events
        std::cout << std::endl;
        std::cout << "Ending Selection over MC" << std::endl;
    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){
        std::cout << "Starting Selection over Data" << std::endl;

        // resize the Passed vector
        data_passed_v.resize(data_tree_total_entries);
        
        // Resize the Counter Vector
        data_counter_v.resize(Passed_Container::k_cuts_MAX);
        for (unsigned int i = 0; i < data_counter_v.size(); i++) data_counter_v.at(i).resize(1, 0 ); 
        
        // Loop over the Events
        for (int event = 0; event < data_tree_total_entries; event++){
            // Alert the user
            if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k " << std::flush;
        
            // Get the entry in the tree
            data_tree->GetEntry(event); // TPC Objects

            // The total number of TPC Objects
            int n_tpc_obj = data_tpc_object_container_v->size();

            // Largest flash information
            std::vector<double> largest_flash_v = data_largest_flash_v_v.at(event); // Vec with the largest flash

            // Create an instance of the selection cut (also initialises flash info)
            selection_cuts_instance.at(histogram_helper::k_data).SetFlashVariables(largest_flash_v);
             
            // Loop over the TPC Objects ---------------------------------------
            // (In Pandora Consolidated, there should be 1 TPC Object per event)
            for (int i = 0; i < n_tpc_obj; i++){
                const xsecAna::TPCObjectContainer tpc_obj = data_tpc_object_container_v->at(i); // Get the TPC Obj

                // Apply the selection cuts 
                bool pass = ApplyCuts(histogram_helper::k_data, event, tpc_obj, data_largest_flash_v_v, data_optical_list_pe_v, data_optical_list_flash_time_v, data_counter_v, data_passed_v);
                if (!pass) continue;

            } // End loop over the TPC Objects

        
        } // End loop over the Events
        std::cout << std::endl;
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "Starting Selection over EXT" << std::endl;
         
        // resize the Passed vector
        ext_passed_v.resize(ext_tree_total_entries);
        
        // Resize the Counter Vector
        ext_counter_v.resize(Passed_Container::k_cuts_MAX);
        for (unsigned int i = 0; i < ext_counter_v.size(); i++) ext_counter_v.at(i).resize(1, 0 ); 
        
        // Loop over the Events
        for (int event = 0; event < ext_tree_total_entries; event++){
            // Alert the user
            if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k " << std::flush;
        
            // Get the entry in the tree
            ext_tree->GetEntry(event);  // TPC Objects

            // The total number of TPC Objects
            int n_tpc_obj = ext_tpc_object_container_v->size();

            // Largest flash information
            std::vector<double> largest_flash_v = ext_largest_flash_v_v.at(event); // Vec with the largest flash

            // Create an instance of the selection cut (also initialises flash info)
            selection_cuts_instance.at(histogram_helper::k_ext).SetFlashVariables(largest_flash_v);
             
            // Loop over the TPC Objects ---------------------------------------
            // (In Pandora Consolidated, there should be 1 TPC Object per event)
            for (int i = 0; i < n_tpc_obj; i++){
                const xsecAna::TPCObjectContainer tpc_obj = ext_tpc_object_container_v->at(i); // Get the TPC Obj

                // Apply the selection cuts 
                bool pass = ApplyCuts(histogram_helper::k_ext, event, tpc_obj, ext_largest_flash_v_v, ext_optical_list_pe_v, ext_optical_list_flash_time_v, ext_counter_v, ext_passed_v);
                if (!pass) continue;

            } // End loop over the TPC Objects

        
        } // End loop over the Events
        std::cout << std::endl;
        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){
        std::cout << "Starting Selection over Dirt" << std::endl;
         
        // resize the Passed vector
        dirt_passed_v.resize(dirt_tree_total_entries);
        
        // Resize the Counter Vector
        dirt_counter_v.resize(Passed_Container::k_cuts_MAX);
        for (unsigned int i = 0; i < dirt_counter_v.size(); i++) dirt_counter_v.at(i).resize(1, 0 ); 
        
        // Loop over the Events
        for (int event = 0; event < dirt_tree_total_entries; event++){
            // Alert the user
            if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k " << std::flush;
        
            // Get the entry in the tree
            dirt_tree->GetEntry(event); // TPC Objects

            // The total number of TPC Objects
            int n_tpc_obj = dirt_tpc_object_container_v->size();

            // Largest flash information
            std::vector<double> largest_flash_v = dirt_largest_flash_v_v.at(event); // Vec with the largest flash

            // Create an instance of the selection cut (also initialises flash info)
            selection_cuts_instance.at(histogram_helper::k_dirt).SetFlashVariables(largest_flash_v);
             
            // Loop over the TPC Objects ---------------------------------------
            // (In Pandora Consolidated, there should be 1 TPC Object per event)
            for (int i = 0; i < n_tpc_obj; i++){
                const xsecAna::TPCObjectContainer tpc_obj = dirt_tpc_object_container_v->at(i); // Get the TPC Obj

                // Apply the selection cuts 
                bool pass = ApplyCuts(histogram_helper::k_dirt, event, tpc_obj, dirt_largest_flash_v_v, dirt_optical_list_pe_v, dirt_optical_list_flash_time_v, dirt_counter_v, dirt_passed_v);
                if (!pass) continue;

            } // End loop over the TPC Objects

        
        } // End loop over the Events
        std::cout << std::endl;
        std::cout << "Ending Selection over Dirt" << std::endl;

    }
    // -------------------------------------------------------------------------
    int num_ext{0}, num_dirt{0};

    // Print the results of the selection -- needs configuring for all the correct 
    // scale factors. For now they are all set to 1.
    if (bool_use_mc) { 
        for (unsigned int k = 0; k < mc_counter_v.size(); k++) {
            if (bool_use_ext)  num_ext  = ext_counter_v .at(k).at(0);
            if (bool_use_dirt) num_dirt = dirt_counter_v.at(k).at(0);
            if (bool_use_mc)   selection_cuts_instance.at(histogram_helper::k_mc)  .PrintInfo(total_mc_entries_inFV, mc_counter_v.at(k), num_ext, intime_scale_factor, data_scale_factor, num_dirt, dirt_scale_factor, cut_names.at(k));
            if (bool_use_data) selection_cuts_instance.at(histogram_helper::k_data).PrintInfoData(data_counter_v.at(k).at(0));
        }
    }

    // If no MC file specified then just print the data
    if (bool_use_data && !bool_use_mc) {
        for (unsigned int k = 0; k < data_counter_v.size(); k++) {
            if (bool_use_data) selection_cuts_instance.at(histogram_helper::k_data).PrintInfoData(data_counter_v.at(k).at(0));
        }
    }

    std::cout << "Finished running the selection!"<< std::endl;
    return;
} // End Selection
// -----------------------------------------------------------------------------
bool selection::ApplyCuts(int type, int event, const xsecAna::TPCObjectContainer &tpc_obj,
                           std::vector<std::vector<double>> &largest_flash_v_v,
                           std::vector<std::vector<int>>    &optical_list_pe_v,
                           std::vector<std::vector<double>> &optical_list_flash_time_v,
                           std::vector<std::vector<double>> &counter_v,
                           std::vector<Passed_Container>    &passed_v){

    std::string type_str;
    // Initalise cut instance with tpc object specifics such as num pfp
    if (type == histogram_helper::k_data){
        type_str = "Data";
        selection_cuts_instance.at(type).SetTPCObjVariables(tpc_obj, type_str);
    } 
    else if (type == histogram_helper::k_ext){ 
        type_str  = "EXT";
        selection_cuts_instance.at(type).SetTPCObjVariables(tpc_obj, type_str);
    }
    else if (type == histogram_helper::k_dirt) {
        type_str  = "Dirt";
        selection_cuts_instance.at(type).SetTPCObjVariables(tpc_obj, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v, has_pi0, type_str);
    }
    else {
        type_str  = "MC";
        selection_cuts_instance.at(type).SetTPCObjVariables(tpc_obj, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v, has_pi0, type_str);
    }

    // Get some variables for filling histograms -- we can put these into a function to initalise the histogram helper and tidy this up later
    int leading_shower_index = selection_cuts_instance.at(type).GetLeadingIndex();

    // Get the classification index for filling stacked histograms
    int classification_index;
    if (!slim) classification_index = histogram_helper_instance.IndexOfClassification(selection_cuts_instance.at(type).tpc_classification.first);
    
    // Here we apply the selection cuts ----------------------------------------
    bool pass; // A flag to see if an event passes an event

    // *************************************************************************
    // Pandora Output ----------------------------------------------------------
    // *************************************************************************
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_pandora_output, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);
    if (!slim) histogram_helper_instance.FillRecoVtxZY(classification_index, tpc_obj);
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_pandora_output), type_str); 

    // *************************************************************************
    // Flash is in time --------------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).FlashinTime(flash_time_start, flash_time_end, optical_list_flash_time_v.at(event), type_str);
    passed_v.at(event).cut_v.at(Passed_Container::k_in_time) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_in_time), type_str); 
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_in_time, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);
    
    // *************************************************************************
    // Flash has more than the required PE -------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).FlashPE(flash_pe_threshold, optical_list_pe_v.at(event));
    passed_v.at(event).cut_v.at(Passed_Container::k_flash_pe) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_flash_pe), type_str); 
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_flash_pe, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);
    
    // *************************************************************************
    // Has a valid Nue ---------------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).HasNue(tpc_obj);
    passed_v.at(event).cut_v.at(Passed_Container::k_has_nue) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_has_nue), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_has_nue, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Is in the FV ------------------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).in_fv(tpc_obj.pfpVtxX(), tpc_obj.pfpVtxY(), tpc_obj.pfpVtxZ(), fv_boundary_v);
    passed_v.at(event).cut_v.at(Passed_Container::k_in_fv) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_in_fv), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_in_fv, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply flash vtx cut -----------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).flashRecoVtxDist(largest_flash_v_v.at(event), tolerance, tpc_obj.pfpVtxY(), tpc_obj.pfpVtxZ());
    passed_v.at(event).cut_v.at(Passed_Container::k_vtx_to_flash) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_vtx_to_flash), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_vtx_to_flash, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply vtx nu distance cut -----------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).VtxNuDistance( tpc_obj, 11, shwr_nue_tolerance);
    passed_v.at(event).cut_v.at(Passed_Container::k_shwr_nue_dist) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_shwr_nue_dist), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_shwr_nue_dist, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply track vtx nu distance cut -----------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).VtxNuDistance( tpc_obj, 13, trk_nue_tolerance);
    passed_v.at(event).cut_v.at(Passed_Container::k_trk_nue_dist) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_trk_nue_dist), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_trk_nue_dist, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Hit threshold cut -------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).HitThreshold(tpc_obj, shwr_hit_threshold, false);
    passed_v.at(event).cut_v.at(Passed_Container::k_shwr_hit_threshold) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_shwr_hit_threshold), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_shwr_hit_threshold, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Hit threshold collection cut --------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).LeadingHitThreshold(tpc_obj, shwr_hit_threshold_collection);
    passed_v.at(event).cut_v.at(Passed_Container::k_shwr_hit_threshold_collection) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_shwr_hit_threshold_collection), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_shwr_hit_threshold_collection, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Open Angle cut ----------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).OpenAngleCut(tpc_obj, tolerance_open_angle_min, tolerance_open_angle_max);
    passed_v.at(event).cut_v.at(Passed_Container::k_shwr_open_angle) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_shwr_open_angle), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_shwr_open_angle, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply dEdx cut ----------------------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).dEdxCut(tpc_obj, tolerance_dedx_min, tolerance_dedx_max, type_str);
    passed_v.at(event).cut_v.at(Passed_Container::k_shwr_dedx) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_shwr_dedx), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_shwr_dedx, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Secondary shower dist cut -----------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).SecondaryShowersDistCut(tpc_obj, dist_tolerance);
    passed_v.at(event).cut_v.at(Passed_Container::k_dist_nue_vtx) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_dist_nue_vtx), type_str); 
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_dist_nue_vtx, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply hit per lengh ratio cut -------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).HitLengthRatioCut( pfp_hits_length_tolerance, tpc_obj);
    passed_v.at(event).cut_v.at(Passed_Container::k_pfp_hits_length) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_pfp_hits_length), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_pfp_hits_length, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Longest Track Leading Shower cut ----------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).LongestTrackLeadingShowerCut(ratio_tolerance, tpc_obj);
    passed_v.at(event).cut_v.at(Passed_Container::k_longest_trk_leading_shwr_length) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_longest_trk_leading_shwr_length), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_longest_trk_leading_shwr_length, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    // Apply Contained Track Cut -----------------------------------------------
    // *************************************************************************
    pass = selection_cuts_instance.at(type).ContainedTracksCut(fv_boundary_v, tpc_obj);
    passed_v.at(event).cut_v.at(Passed_Container::k_trk_contained) = pass;
    if(!pass) return false; // Failed the cut!

    // Set counters for how many passed the cut
    selection_cuts_instance.at(type).TabulateOrigins(counter_v.at(Passed_Container::k_trk_contained), type_str);
    if (!slim) histogram_helper_instance.FillReco(classification_index, histogram_helper::k_trk_contained, tpc_obj, leading_shower_index, largest_flash_v_v.at(event), fv_boundary_v, type_str);

    // *************************************************************************
    return true;

}
// -----------------------------------------------------------------------------
void selection::SavetoFile(){

    // Now saving histograms to file
    std::cout << "Now Saving Histograms to file" << std::endl;
    if (bool_use_mc) {
        histogram_helper_instance.WriteMCTruth("MC");
        histogram_helper_instance.WriteOptical(histogram_helper::k_mc);
        histogram_helper_instance.WriteReco(histogram_helper::k_mc);
    }
    if (bool_use_data) {
        histogram_helper_instance.WriteOptical(histogram_helper::k_data);
        histogram_helper_instance.WriteReco(histogram_helper::k_data);
    }
    if (bool_use_ext) {
        histogram_helper_instance.WriteOptical(histogram_helper::k_ext);
        histogram_helper_instance.WriteReco(histogram_helper::k_ext);
    }
    if (bool_use_dirt) {
        histogram_helper_instance.WriteMCTruth("Dirt");
        histogram_helper_instance.WriteOptical(histogram_helper::k_dirt);
        histogram_helper_instance.WriteReco(histogram_helper::k_dirt);
    }

    

} // End save to file
// -----------------------------------------------------------------------------
void selection::MakeHistograms(){
    std::cout << "Creating histograms and making plots" << std::endl;

    // Loop over the cuts and plot histograms by plot type
    for (unsigned int i = 0 ; i < histogram_helper::k_cuts_MAX; i++){
        
        // Create a set of strings for creating a dynamic directory
        // Directory structure that is created will take the form plots/<cut>/
        std::string a = "if [ ! -d \"plots/";
        std::string b = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
        std::string c = "; fi";
        std::string command = a +histogram_helper_instance.cut_dirs.at(i) + b + histogram_helper_instance.cut_dirs.at(i) + c ;
        system(command.c_str()); 

        // Reco X
        histogram_plotter_instance.MakeStack("h_reco_vtx_x",histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Reco Vertex X [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_vtx_x.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );
        
        // Reco Y
        histogram_plotter_instance.MakeStack("h_reco_vtx_y",histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Reco Vertex Y [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_vtx_y.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Reco Z
        histogram_plotter_instance.MakeStack("h_reco_vtx_z",histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Reco Vertex Z [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_vtx_z.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // dEdx
        histogram_plotter_instance.MakeStack("h_reco_dEdx",histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Collection Plane dEdx [MeV/cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_dEdx_collection.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );
   
        // Leading Shower Momentum
        histogram_plotter_instance.MakeStack("h_reco_leading_mom", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Momentum [MeV/c]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_mom.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // 2D distance largest flash to reco nu vertex
        histogram_plotter_instance.MakeStack("h_reco_flash_to_vtx_dist", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Flash to Vertex Distance [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_flash_to_vtx_dist.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // 2D distance shower vertex to reco nu vertex
        histogram_plotter_instance.MakeStack("h_reco_shower_to_vtx_dist", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Shower to Vertex Distance [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_shower_to_vtx_dist.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // 2D distance track vertex to reco nu vertex
        histogram_plotter_instance.MakeStack("h_reco_track_to_vtx_dist", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Track to Vertex Distance [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_track_to_vtx_dist.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading Shower hits in all planes
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_all_planes", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Hits All Planes", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_hits_all_planes.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading Shower hits in collection
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_collection_plane", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Hits Collection Plane", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_hits_collection_plane.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading Shower opening angle
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_open_angle", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Open Angle [degrees]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_open_angle.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Secondary shower to vertex distance (for events with more than 1 shower)
        histogram_plotter_instance.MakeStack("h_reco_secondary_shower_to_vtx_dist", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Secondary Shower to Vertex Distance (>1 shower) [cm]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_secondary_shower_to_vtx_dist.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading Shower hits per length
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_hits_per_length", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Hits / Length [cm^{-1}]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_hits_per_length.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Longest track to leading shower length
        histogram_plotter_instance.MakeStack("h_reco_longest_track_leading_shower_length", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Longest Track Length / Leading Shower Length", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_longest_track_leading_shower_length.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Track Containment
        histogram_plotter_instance.MakeStack("h_reco_track_contained", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Contained Tracks", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_track_contained.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading shower phi
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_phi", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Phi [degrees]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_phi.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading shower theta
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_theta", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Theta [degrees]", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_theta.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading shower cos theta
        histogram_plotter_instance.MakeStack("h_reco_leading_shower_cos_theta", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Leading Shower Cos(#theta)", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_leading_shower_cos_theta.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading shower multiplicity
        histogram_plotter_instance.MakeStack("h_reco_shower_multiplicity", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Shower Multiplicty", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_shower_multiplicity.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );

        // Leading track multiplicity
        histogram_plotter_instance.MakeStack("h_reco_track_multiplicity", histogram_helper_instance.cut_dirs.at(i).c_str(),
                                           false,  false, "Track Multiplicty", data_scale_factor, 1.0, intime_scale_factor, dirt_scale_factor, 0.8, 0.98, 0.98, 0.50,
                                           Form("plots/%s/reco_track_multiplicity.pdf", histogram_helper_instance.cut_dirs.at(i).c_str()) );
   
    }
    

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection