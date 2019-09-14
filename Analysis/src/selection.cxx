#include "../include/selection.h"

namespace xsecSelection {

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
        mc_counter_v.resize(k_cuts_MAX);
        for (unsigned int i = 0; i < mc_counter_v.size(); i++) mc_counter_v.at(i).resize(24, 0.0 );
        
        // Loop over the Events
        for (int event = 0; event < tree_total_entries; event++){
            // Alert the user
    	    if (event % 100000 == 0) std::cout << "On entry " << event/100000.0 <<"00k" << std::endl;
        
            // Get the entry in the tree
            mytree->GetEntry(event);               // TPC Objects
            mctruth_counter_tree->GetEntry(event); // MC Tree
            
            // The total number of TPC Objects
            int n_tpc_obj = tpc_object_container_v->size();

            // Largest flash information
            std::vector<double> largest_flash_v = mc_largest_flash_v_v.at(event); // Vec with the largest flash

            // Create an instance of the selection cut (also initialises flash info)
            selection_cuts selection_cuts_instance;
            selection_cuts_instance.SetFlashVariables(largest_flash_v);
            
            // Loop over the TPC Objects ---------------------------------------
            // (In Pandora Consolidated, there should be 1 TPC Object per event)
            for (int i = 0; i < n_tpc_obj; i++){
                const xsecAna::TPCObjectContainer tpc_obj = tpc_object_container_v->at(i); // Get the TPC Obj

                // Initalise cut instance with tpc object specifics such as num pfp
                selection_cuts_instance.SetTPCObjVariables(tpc_obj, mc_nu_vtx_x, mc_nu_vtx_y, mc_nu_vtx_z, fv_boundary_v, has_pi0);

                // Here we apply the selection cuts ----------------------------
                bool pass;

                // Flash is in time and has more than the required PE ----------
                pass = selection_cuts_instance.FlashinTime_FlashPE(flash_time_start, flash_time_end, flash_pe_threshold, mc_optical_list_flash_time_v.at(event), mc_optical_list_pe_v.at(event));
                mc_passed_v.at(event).cut_v.at(k_flash_pe_intime) = pass;
                selection_cuts_instance.TabulateOrigins(mc_counter_v.at(k_flash_pe_intime));
                if(!pass) continue; // Failed the cut!
                
                // Has a valid Nue ---------------------------------------------
                pass = selection_cuts_instance.HasNue(tpc_obj);
                mc_passed_v.at(event).cut_v.at(k_has_nue) = pass;
                if(!pass) continue; // Failed the cut!

                // Is in the FV ------------------------------------------------
                pass = selection_cuts_instance.in_fv(tpc_obj.pfpVtxX(), tpc_obj.pfpVtxY(), tpc_obj.pfpVtxZ(), fv_boundary_v);
                mc_passed_v.at(event).cut_v.at(k_in_fv) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply flash vtx cut -----------------------------------------
                pass = selection_cuts_instance.flashRecoVtxDist(mc_largest_flash_v_v.at(event), tolerance, tpc_obj.pfpVtxX(), tpc_obj.pfpVtxY(), tpc_obj.pfpVtxZ());
                mc_passed_v.at(event).cut_v.at(k_vtx_to_flash) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply vtx nu distance cut -----------------------------------
                pass = selection_cuts_instance.VtxNuDistance( tpc_obj, 11, shwr_nue_tolerance);
                mc_passed_v.at(event).cut_v.at(k_shwr_nue_dist) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply track vtx nu distance cut -----------------------------
                pass = selection_cuts_instance.VtxNuDistance( tpc_obj, 13, trk_nue_tolerance);
                mc_passed_v.at(event).cut_v.at(k_trk_nue_dist) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply Hit threshold cut -------------------------------------
                pass = selection_cuts_instance.HitThreshold(tpc_obj, shwr_hit_threshold, false);
                mc_passed_v.at(event).cut_v.at(k_shwr_hit_threshold) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply Hit threshold collection cut --------------------------
                pass = selection_cuts_instance.HitThreshold(tpc_obj, shwr_hit_threshold_collection, true);
                mc_passed_v.at(event).cut_v.at(k_shwr_hit_threshold_collection) = pass;
                if(!pass) continue; // Failed the cut!
            
                // Apply Open Angle cut ----------------------------------------
                pass = selection_cuts_instance.OpenAngleCut(tpc_obj, tolerance_open_angle_min, tolerance_open_angle_max);
                mc_passed_v.at(event).cut_v.at(k_shwr_open_angle) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply dEdx cut ----------------------------------------------
                pass = selection_cuts_instance.dEdxCut(tpc_obj, tolerance_dedx_min, tolerance_dedx_max);
                mc_passed_v.at(event).cut_v.at(k_shwr_dedx) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply Secondary shower dist cut -----------------------------
                pass = selection_cuts_instance.SecondaryShowersDistCut(tpc_obj, dist_tolerance);
                mc_passed_v.at(event).cut_v.at(k_dist_nue_vtx) = pass;
                if(!pass) continue; // Failed the cut! 

                // Apply hit per lengh ratio cut -------------------------------
                pass = selection_cuts_instance.HitLengthRatioCut( pfp_hits_length_tolerance, tpc_obj);
                mc_passed_v.at(event).cut_v.at(k_pfp_hits_length) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply Longest Track Leading Shower cut ----------------------
                pass = selection_cuts_instance.LongestTrackLeadingShowerCut(ratio_tolerance, tpc_obj);
                mc_passed_v.at(event).cut_v.at(k_longest_trk_leading_shwr_length) = pass;
                if(!pass) continue; // Failed the cut!

                // Apply Contained Track Cut -----------------------------------
                pass = selection_cuts_instance.ContainedTracksCut(fv_boundary_v, tpc_obj);
                mc_passed_v.at(event).cut_v.at(k_trk_contained) = pass;
                if(!pass) continue; // Failed the cut!

            } // End loop over the TPC Objects

            
        
        } // End loop over the Events
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(0)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(1)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(2)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(3)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(4)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(5)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(6)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(7)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(8)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(9)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(10)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(11)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(12)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(13)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(14)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(15)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(16)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(17)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(18)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(19)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(20)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(21)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(22)  << std::endl;
        std::cout << mc_counter_v.at(k_flash_pe_intime).at(23)  << std::endl;

        std::cout << "Ending Selection over MC" << std::endl;
    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){

    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){

    }
    // -------------------------------------------------------------------------


    // Now fill hisograms and write to a file
    

    // Plots the histograms and write to a file
    
    
    
    std::cout << "Finished running the selection!"<< std::endl;
    return;
}



} // END NAMESPACE xsecSelection

