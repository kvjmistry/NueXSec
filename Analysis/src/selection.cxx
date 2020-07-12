#include "../include/selection.h"

namespace xsecSelection {
// -----------------------------------------------------------------------------
void selection::Initialise( const char * mc_file,
                            const char * ext_file,
                            const char * data_file,
                            const char * dirt_file,
                            const char * mc_file_out,
                            const char * ext_file_out,
                            const char * data_file_out,
                            const char * dirt_file_out,
                            const char * mc_tree_file_name_out,
                            utility _utility,
                            bool _slim,
                            int num_events,
                            const char * run_period,
                            int _verbose,
                            int weight_cfg){
    
    std::cout << "\nInitialising..." << std::endl;

    _util = _utility;

    // Initialise the selection cuts class
    _scuts.Initalise(_utility);

    // Display slimmed selection
    if (_slim){
        std::cout << "\033[0;32m-------------------------------" << std::endl;
        std::cout << "     Running in Slim Mode!" << std::endl;
        std::cout << "-------------------------------\033[0m" << std::endl;
        slim = _slim;
    }

    std::cout << "\nSetting verbose level to: " << _verbose << std::endl;
    verbose = _verbose;

    std::cout << "\nUsing a weight setting of: " << weight_cfg << std::endl;
    std::cout << "If this is set to 1 (default) then we apply the Genie Tune and PPFX weights to the CV" << std::endl;
    _weight_cfg = weight_cfg;

    // Create the file directory if it does not exist already
    gSystem->Exec("if [ ! -d \"files/trees/\" ]; then echo \"\nfiles folder does not exist... creating\"; mkdir -p files/trees; fi"); 

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run1_Data_trig) / _util.config_v.at(_util.k_Run1_EXT_trig);
        _run_period = 1;
    }
    else if (strcmp(run_period, "3") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run3_Data_trig) / _util.config_v.at(_util.k_Run3_EXT_trig);
        _run_period = 3;
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }


    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:   "   << mc_scale_factor     << "\n" <<
    "Dirt Scale factor: "   << dirt_scale_factor   << "\n" <<
    "EXT Scale factor:  "   << intime_scale_factor << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // Set the maximum number of events tp process
    if (num_events > 0 ) max_events = num_events;

    // Resize the histogram helper instance vectors, one instance per type e.g MC, data, ..
    _hhelper.resize(_util.k_type_MAX);
    
    // Resize the tree helper instance vectors, one instance per type e.g MC, data, ..
    _thelper.resize(_util.k_type_MAX);

    // Print the input files
    std::cout <<
    "Run Period Configured: run" << run_period<<"\n" <<
    "MC   File Path:      " << mc_file        <<"\n" <<
    "Ext  File Path:      " << ext_file       <<"\n" <<
    "Data File Path:      " << data_file      <<"\n" <<
    "Dirt File Path:      " << dirt_file      <<"\n" <<
    std::endl;

    // Now get the files, if file isnt specified then set bool to skip
    bool_use_mc        = _util.GetFile(f_mc,        mc_file);
    bool_use_ext       = _util.GetFile(f_ext,       ext_file);
    bool_use_data      = _util.GetFile(f_data,      data_file);
    bool_use_dirt      = _util.GetFile(f_dirt,      dirt_file);

    // Load in the flux weights file
    std::cout << "Getting the CV flux file..."<< std::endl;
    if (strcmp(run_period, "1") == 0) f_flux_weights = new TFile("../Systematics/f_flux_CV_weights_fhc.root", "READ");
    if (strcmp(run_period, "3") == 0) f_flux_weights = new TFile("../Systematics/f_flux_CV_weights_rhc.root", "READ");

    // Resize the counter vector
    counter_v.resize(_util.k_cuts_MAX);

    for (unsigned int t = 0; t < counter_v.size(); t++){
        counter_v.at(t).resize(_util.k_COUNTER_MAX, 0.0);
    }

    // Get MC variables --------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nInitialising MC" << std::endl;
        _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the mc slice container
        mc_SC.Initialise(mc_tree, _util.k_mc, f_flux_weights, run_period, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_mc).Initialise(_util.k_mc, run_period, mc_tree_file_name_out);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_mc).Initialise(_util.k_mc, run_period, mc_file_out, weight_cfg, _util);
        if (!_slim) _hhelper.at(_util.k_mc).InitHistograms();

        mc_tree_total_entries = mc_tree->GetEntries();
        std::cout << "Total MC Events:         " << mc_tree_total_entries << std::endl;

        // Resize the Passed vector
        mc_passed_v.resize(mc_tree_total_entries);

        for (unsigned int y = 0; y < mc_passed_v.size(); y++ ){
            mc_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "Initialisation of MC Complete!" << std::endl;
    } // End getting MC variables

    // Initialise Data specific ------------------------------------------------
    if (bool_use_data){
        std::cout << "\nInitialising Data" << std::endl;
        _util.GetTree(f_data, data_tree, "nuselection/NeutrinoSelectionFilter");
        
        // Initialise all the data slice container
        data_SC.Initialise(data_tree, _util.k_data, f_flux_weights, run_period, _util);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_data).Initialise(_util.k_data, run_period, data_file_out, weight_cfg, _util);
        if (!_slim) _hhelper.at(_util.k_data).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_data).Initialise(_util.k_data, run_period, "empty");

        data_tree_total_entries = data_tree->GetEntries();
        std::cout << "Total Data Events:         " << data_tree_total_entries << std::endl;

        // Resize the Passed vector
        data_passed_v.resize(data_tree_total_entries);
        
        for (unsigned int y = 0; y < data_passed_v.size(); y++ ){
            data_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "Initialisation of Data Complete!" << std::endl;

    } // End intialisation of Data variables

    // Initialise EXT specific -------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nInitialising EXT" << std::endl;

        _util.GetTree(f_ext, ext_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        ext_SC.Initialise(ext_tree, _util.k_ext, f_flux_weights, run_period, _util);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_ext).Initialise(_util.k_ext, run_period, ext_file_out, weight_cfg, _util);
        if (!_slim) _hhelper.at(_util.k_ext).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_ext).Initialise(_util.k_ext, run_period, "empty");

        ext_tree_total_entries = ext_tree->GetEntries();
        std::cout << "Total EXT Events:        " << ext_tree_total_entries << std::endl;

        // Resize the Passed vector
        ext_passed_v.resize(ext_tree_total_entries);

        for (unsigned int y = 0; y < ext_passed_v.size(); y++ ){
            ext_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "Initialisation of EXT Complete!" << std::endl;

    } // End intialisation of ext variables

    // Initialise Dirt specific ------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nInitialising Dirt" << std::endl;

        _util.GetTree(f_dirt, dirt_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        dirt_SC.Initialise(dirt_tree, _util.k_dirt, f_flux_weights, run_period, _util);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_dirt).Initialise(_util.k_dirt, run_period, dirt_file_out, weight_cfg, _util);
        if (!_slim) _hhelper.at(_util.k_dirt).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_dirt).Initialise(_util.k_dirt, run_period, "empty");

        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total Dirt Events:         " << dirt_tree_total_entries << std::endl;

        // Resize the Passed vector
        dirt_passed_v.resize(dirt_tree_total_entries);
        
        for (unsigned int y = 0; y < dirt_passed_v.size(); y++ ){
            dirt_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "Initialisation of Dirt Complete!" << std::endl;


    } // End intialisation of dirt variables    
    
    // Invoke main selection function
    MakeSelection();

} // END Initialise function
// -----------------------------------------------------------------------------
// Main function for selection
void selection::MakeSelection(){
    std::cout << "\n\033[0;32mNow Running the selection!\033[0m"<< std::endl;
    
    int counter = 0;

    // MC ----------------------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nStarting Selection over MC" << std::endl;

        int numu_pi{0}, numubar_pi{0}, nue_pi{0}, nuebar_pi{0}; // Total number of 1pi events

        // Event loop
        for (int ievent = 0; ievent < mc_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;

            // Get the entry in the tree
            mc_tree->GetEntry(ievent); 

            // Get the 1 pion events for kirsty
            if (mc_SC.isVtxInFiducial){
                if (mc_SC.nu_pdg == 14){
                    if (mc_SC.npion == 1) numu_pi ++;
                }
                if (mc_SC.nu_pdg == -14){
                    if (mc_SC.npion == 1) numubar_pi ++;
                }

                if (mc_SC.nu_pdg == 12){
                    if (mc_SC.npion == 1) nue_pi ++;
                }

                if (mc_SC.nu_pdg == -12){
                    if (mc_SC.npion == 1) nuebar_pi ++;
                }
            }
            

            // std::cout << mc_SC.run << " " << mc_SC.sub<<" " << mc_SC.evt<<  std::endl;

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_mc, mc_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_mc, mc_SC);
            
            // Apply the selection cuts 
            bool pass = ApplyCuts(_util.k_mc, ievent, counter_v, mc_passed_v, mc_SC);
            if (!pass) continue;

        } // End Event loop

        std::cout << "Ending Selection over MC" << std::endl;
        std::cout << "numu_pi: " << numu_pi << std::endl;
        std::cout << "numubar_pi: " << numubar_pi << std::endl;
        std::cout << "nue_pi: "  << nue_pi  << std::endl;
        std::cout << "nuebar_pi: "  << nuebar_pi  << std::endl;

        // Loop again to look at background events that still pass
        // Event loop
        // for (int ievent = 0; ievent < mc_tree_total_entries; ievent++){
            
        //     if (mc_passed_v.at(ievent).cut_v.at(_util.k_cuts_MAX - 1 ) == true ){
            
        //         mc_tree->GetEntry(ievent); 

        //         std::pair<std::string, int> classification = mc_SC.SliceClassifier(_util.k_mc);

        //         // Background events
        //         if (classification.second == _util.k_nue_cc){
        //             if (mc_SC.shr_tkfit_dedx_Y > 6.8){
        //                 std::cout <<  mc_SC.run << " " << mc_SC.sub << " " << mc_SC.evt <<  std::endl;
        //             }

        //         }
        //     }

        // } // End Event loop


    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){
        std::cout << "\nStarting Selection over Data" << std::endl;

        for (int ievent = 0; ievent < data_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            data_tree->GetEntry(ievent);

            // Skip the RHC events contaminated in the FHC files
            if (_run_period == 3 && data_SC.run < 16880 ){
                //continue;
            }

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_data, data_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_data, data_SC);

            bool pass = ApplyCuts(_util.k_data, ievent, counter_v, data_passed_v, data_SC);
            if (!pass) continue;
        }

        // Make the run subrun event file list after the selection 
        if (make_list){

            std::cout << "Making the run subrun list for selected data events..." << std::endl;

            int run, subrun, event;
            std::ofstream run_subrun_file;
            run_subrun_file.open(Form("files/run%i_run_subrun_list_data.txt",_run_period));

            // Loop over the data events and make a run_subrun_event filelist
            for (int ievent = 0; ievent < data_tree_total_entries; ievent++){
            
                if (data_passed_v.at(ievent).cut_v.at(_util.k_cuts_MAX - 1 ) == true ){
                    
                    // Get the entry in the tree
                    data_tree->GetEntry(ievent);

                    run    = data_SC.run;
                    subrun = data_SC.sub;
                    event  = data_SC.evt;

                    double INTERCEPT = 0.0;
                    double SLOPE = 0.83;
                    double reco_nu_e = (data_SC.shr_energy_tot_cali + INTERCEPT) / SLOPE + data_SC.trk_energy_tot;

                    // std::cout << run << " " << subrun << " " << event << " " << reco_nu_e <<  '\n';

                    run_subrun_file << run << " " << subrun << " " << event << '\n';

                }
                
            }

            run_subrun_file.close();
        }
        
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nStarting Selection over EXT" << std::endl;

        for (int ievent = 0; ievent < ext_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            ext_tree->GetEntry(ievent); // TPC Objects

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_ext, ext_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_ext, ext_SC);

            bool pass = ApplyCuts(_util.k_ext, ievent, counter_v, ext_passed_v, ext_SC);
            if (!pass) continue;
        }
         
        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nStarting Selection over Dirt" << std::endl;

        for (int ievent = 0; ievent < dirt_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            dirt_tree->GetEntry(ievent);

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_dirt, dirt_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_dirt, dirt_SC);

            bool pass = ApplyCuts(_util.k_dirt, ievent, counter_v, dirt_passed_v, dirt_SC);
            if (!pass) continue;
        }
         
        std::cout << "Ending Selection over Dirt" << std::endl;

    }
    // -------------------------------------------------------------------------
    std::cout << "Finished running the selection!"<< std::endl;

    // Print information from the selection, loop over the cuts
    for (unsigned int p=0; p < counter_v.size();p++){

        // Fill the counter trees
        if (bool_use_mc)   _thelper.at(_util.k_mc)  .Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_data) _thelper.at(_util.k_data).Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_ext)  _thelper.at(_util.k_ext) .Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_dirt) _thelper.at(_util.k_dirt).Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);

    }
    
    // Now save all the outputs to file
    if (!slim) SavetoFile();

    return;
} // End Selection
// -----------------------------------------------------------------------------
bool selection::ApplyCuts(int type, int ievent,std::vector<std::vector<double>> &counter_v,
                           std::vector<Passed_Container> &passed_v, SliceContainer &SC){

    // Here we apply the selection cuts ----------------------------------------
    bool pass; // A flag to see if an event passes an event

    
    // Classify the event
    std::pair<std::string, int> classification = SC.SliceClassifier(type);      // Classification of the event
    std::string interaction                    = SC.SliceInteractionType(type); // Genie interaction type
    //std::string interaction = "nue_cc_qe";
    std::string category                       = SC.SliceCategory();            // The pandora group slice category
    std::pair<std::string, int> particle_type  = SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower



    // Test code to isolate the low E nues in truth
    // if (type == _util.k_mc && SC.nu_e > 0.5) return false;
    // if (type == _util.k_mc && SC.shr_dedx_Y_cali > 7 && classification.second == _util.k_nue_cc){
        // std::cout << SC.run << " " << SC.sub<<" " << SC.evt<<  std::endl;
    // }

    // *************************************************************************
    // Unselected---------------------------------------------------------------
    // *************************************************************************
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_unselected, counter_v );
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    passed_v.at(ievent).cut_v.at(_util.k_swtrig) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_swtrig, counter_v );

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    passed_v.at(ievent).cut_v.at(_util.k_slice_id) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_slice_id, counter_v );
    
    // *************************************************************************
    // Electron Candidate ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.e_candidate(SC);
    passed_v.at(ievent).cut_v.at(_util.k_e_candidate) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_e_candidate, counter_v );

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    passed_v.at(ievent).cut_v.at(_util.k_in_fv) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_in_fv, counter_v );
    
    // *************************************************************************
    // Topological Score -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.topo_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_topo_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_topo_score, counter_v );

    // *************************************************************************
    // Slice Contained Fraction ------------------------------------------------
    // *************************************************************************
    pass = _scuts.contained_frac(SC);
    passed_v.at(ievent).cut_v.at(_util.k_contained_frac) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_contained_frac, counter_v );

    // *************************************************************************
    // Shower Score ------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shower_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shower_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_shower_score, counter_v );

    // *************************************************************************
    // Michel Rejection --------------------------------------------------------
    // *************************************************************************
    pass = _scuts.michel_rej(SC);
    passed_v.at(ievent).cut_v.at(_util.k_michel_rej) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_michel_rej, counter_v );

    // *************************************************************************
    // Shower Hits -------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_hits(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shr_hits) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_shr_hits, counter_v );

    // *************************************************************************
    // Shower Hit Ratio  -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_hitratio(SC);
    passed_v.at(ievent).cut_v.at(_util.k_hit_ratio) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_hit_ratio, counter_v );

    // *************************************************************************
    // Shower Moliere Average --------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_moliere_avg(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shr_moliere_avg) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_shr_moliere_avg, counter_v );

    // *************************************************************************
    // 2D cut for Shower to Vertex Distance and dEdx ---------------------------
    // *************************************************************************
    pass = _scuts.shr_dist_dEdx_y(SC);
    passed_v.at(ievent).cut_v.at(_util.k_vtx_dist_dedx) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_vtx_dist_dedx, counter_v );

    // *************************************************************************
    // dEdx in y plane for 0 track events --------------------------------------
    // *************************************************************************
    pass = _scuts.dEdx_y_no_tracks(SC);
    passed_v.at(ievent).cut_v.at(_util.k_dEdx_y_no_tracks) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, classification, interaction, particle_type, _util.k_dEdx_y_no_tracks, counter_v );


    // // *************************************************************************
    // // Shower to Vertex Distance --------------------------------------------
    // // *************************************************************************
    // pass = _scuts.shr_distance(SC);
    // passed_v.at(ievent).cut_v.at(_util.k_shr_distance) = pass;
    // if(!pass) return false; // Failed the cut!
    
    // SelectionFill(type, SC, classification, interaction, particle_type, _util.k_shr_distance, counter_v );

    // // *************************************************************************
    // // dEdx in y plane ---------------------------------------------------------
    // // *************************************************************************
    // pass = _scuts.dEdx_y(SC);
    // passed_v.at(ievent).cut_v.at(_util.k_dEdx_y) = pass;
    // // if(!pass) return false; // Failed the cut!
    
    // SelectionFill(type, SC, classification, interaction, particle_type, _util.k_dEdx_y, counter_v );

    // ************************************************************************n*
    return true;

}
// -----------------------------------------------------------------------------
void selection::SavetoFile(){

    // Now saving histograms to file
    std::cout << "Now Saving Histograms to file" << std::endl;
    if (bool_use_mc) {

        _hhelper.at(_util.k_mc).WriteTrue();
        _hhelper.at(_util.k_mc).WriteTEfficiency();
        _hhelper.at(_util.k_mc).WriteReco(_util.k_mc);
        _hhelper.at(_util.k_mc).WriteRecoPar(_util.k_mc);
        _hhelper.at(_util.k_mc).WriteFlash();
        _hhelper.at(_util.k_mc).WriteInteractions();
        _hhelper.at(_util.k_mc).Write_2DSigBkgHists();
        _hhelper.at(_util.k_mc).WritePiZero(_util.k_mc);
        _hhelper.at(_util.k_mc).WriteNuMu(_util.k_mc);

        _thelper.at(_util.k_mc).WriteTree();


    }
    if (bool_use_data) {
        _hhelper.at(_util.k_data).WriteReco(_util.k_data);
        _hhelper.at(_util.k_data).WriteRecoPar(_util.k_data);
        _hhelper.at(_util.k_data).WriteFlash();
        _hhelper.at(_util.k_data).WritePiZero(_util.k_data);
        _hhelper.at(_util.k_data).WriteNuMu(_util.k_data);

        _thelper.at(_util.k_data).WriteTree();

    }
    if (bool_use_ext) {
        _hhelper.at(_util.k_ext).WriteReco(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteRecoPar(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteFlash();
        _hhelper.at(_util.k_ext).Write_2DSigBkgHists();
        _hhelper.at(_util.k_ext).WritePiZero(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteNuMu(_util.k_ext);

        _thelper.at(_util.k_ext).WriteTree();

    }
    if (bool_use_dirt) {
        // _hhelper.at(_util.k_dirt).WriteTrue(); // Only turn this on to inspect. It wont merge since the file names are not uniique yet!!
        _hhelper.at(_util.k_dirt).WriteReco(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteRecoPar(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteFlash();
        _hhelper.at(_util.k_dirt).Write_2DSigBkgHists();
        _hhelper.at(_util.k_dirt).WritePiZero(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteNuMu(_util.k_dirt);

        _thelper.at(_util.k_dirt).WriteTree();

    }

} // End save to file
// -----------------------------------------------------------------------------
void selection::SelectionFill(int type, SliceContainer &SC, std::pair<std::string, int> classification, std::string interaction, std::pair<std::string, int> par_type, int cut_index, std::vector<std::vector<double>> &counter_v){
    
    // Get the CV weight
    double weight = 1.0;
    weight = GetCVWeight(type, SC);

    // Try scaling the pi0
    // 0 == no weighting, 1 == normalisation fix, 2 == energy dependent scaling
    GetPiZeroWeight(weight, 2, SC);

    
    // This is in many places, need to have a way for setting this number by default
    double INTERCEPT = 0.0;
    double SLOPE = 0.83;
    double reco_nu_e = (SC.shr_energy_tot_cali + INTERCEPT) / SLOPE + SC.trk_energy_tot;

    // Fill Histograms
    if (!slim) _hhelper.at(type).FillHists(type, classification.second, interaction, par_type.second, cut_index, SC, weight);

    // Set counters for the cut
    _util.Tabulate(SC.isVtxInFiducial, interaction, classification.first, type, counter_v.at(cut_index), weight );

    // Fill Plots for Efficiency
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(cut_index, classification.first, SC, weight);

    // For the last cut we fill the tree  or the first cut and nue_cc (generated and unselected)
    if ( (cut_index == _util.k_cuts_MAX - 1) || (cut_index == _util.k_unselected && (classification.second == _util.k_nue_cc || classification.second == _util.k_nuebar_cc) ) ){

        // This is a generated event, but unselected
        if (cut_index == _util.k_unselected && (classification.second == _util.k_nue_cc || classification.second == _util.k_nuebar_cc )){
            _thelper.at(type).FillVars(SC, classification, true, weight, reco_nu_e);
        }
        else {
            _thelper.at(type).FillVars(SC, classification, false, weight, reco_nu_e);
        }

    }

    // Fill the dedx ttree before shr dist cut and after cut dedx
    // if (cut_index == _util.k_shr_distance - 1 || cut_index == _util.k_shr_distance || cut_index == _util.k_dEdx_y ){
    //     _thelper.at(type).Fill_dedxVars(SC, classification, _util.cut_dirs.at(cut_index), weight);
    // }

}
// -----------------------------------------------------------------------------
double selection::GetCVWeight(int type, SliceContainer SC){
    

    // Always give weights of 1 to the off beam and data
    if (type == _util.k_data || type == _util.k_ext) return 1.0;

    double weight = 1.0;
    bool weight_tune{true}, weight_ppfx{true};

    // Set the weight settings
    if (_weight_cfg == 0){
        weight_tune = false;
        weight_ppfx = false;
        return weight;
    }
    else if (_weight_cfg == 1){
        weight_tune = true;
        weight_ppfx = true;
    }
    else if (_weight_cfg == 2){
        weight_tune = true;
        weight_ppfx = false;
    }
    else if (_weight_cfg == 3){
        weight_tune = false;
        weight_ppfx = true;
    }
    else {
        std::cout << "Unknown weight setting specified, using defaults" << std::endl;
    }


    // Get the tune weight
    weight = SC.weightSplineTimesTune; // Here define the weight
    
    // Catch infinate/nan/unreasonably large tune weights
    if (std::isinf(weight))      weight = 1.0; 
    if (std::isnan(weight) == 1) weight = 1.0;
    if (weight > 100)            weight = 1.0;
    if (weight < 0)              weight = 1.0;

    // If tune weight turned off, just set weight to 1.0
    if (!weight_tune) weight = 1.0;

    // Get the PPFX CV flux correction weight
    double weight_flux = SC.ppfx_cv;

    if (std::isinf(weight_flux))      weight_flux = 1.0; 
    if (std::isnan(weight_flux) == 1) weight_flux = 1.0;
    if (weight_flux > 100)            weight_flux = 1.0;
    if (weight_flux < 0)              weight_flux = 1.0;

    if (weight_ppfx) weight = weight * weight_flux;

    // std::cout << SC.weightSplineTimesTune << "   "<< SC.ppfx_cv << std::endl;

    return weight;

}
// -----------------------------------------------------------------------------
void selection::ApplyPiZeroSelection(int type, SliceContainer &SC){

    bool pass; // A flag to see if an event passes an event

    // Classify the event
    std::pair<std::string, int> classification = SC.SliceClassifier(type);      // Classification of the event
    std::string interaction                    = SC.SliceInteractionType(type); // Genie interaction type
    std::pair<std::string, int> particle_type  = SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    if(!pass) return; // Failed the cut!
    
    // *************************************************************************
    // Common Optical Filter PE  -----------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_pe(SC, type);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // Common Optical Filter Veto  ---------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_veto(SC, type);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    if(!pass) return; // Failed the cut!
    
    // *************************************************************************
    // Pi0 Selection Cuts ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.pi_zero_cuts(SC);
    if(!pass) return; // Failed the cut!

    // Get the Central Value weight
    double weight = GetCVWeight(type, SC);

    double weight_norm = weight;
    double weight_Escale = weight;

    // Try scaling the pi0
    // 0 == no weighting, 1 == normalisation fix, 2 == energy dependent scaling
    GetPiZeroWeight(weight, 0, SC);
    GetPiZeroWeight(weight_norm, 1, SC);
    GetPiZeroWeight(weight_Escale, 2, SC);

    // Now Fill the histograms
    if (!slim) _hhelper.at(type).FillPiZeroHists(classification.second, SC, weight, 0);
    if (!slim) _hhelper.at(type).FillPiZeroHists(classification.second, SC, weight_norm, 1);
    if (!slim) _hhelper.at(type).FillPiZeroHists(classification.second, SC, weight_Escale, 2);


}
// -----------------------------------------------------------------------------
void selection::GetPiZeroWeight(double &weight, int pizero_mode, SliceContainer &SC){

    // Fix the normalisation
    if (pizero_mode == 1){
        
        if (SC.npi0 > 0) {
            weight = weight * 0.759;
        }

    }
    // Try energy dependent scaling for pi0
    else if (pizero_mode == 2){
        
        if (SC.npi0 > 0) {
            double pi0emax = 0.6;
            if (SC.pi0_e > 0.1 && SC.pi0_e < pi0emax){
                weight = weight * (1 - 0.4 * SC.pi0_e);
            }
            else if (SC.pi0_e > 0.1 && SC.pi0_e >= pi0emax){
                weight = weight * (1 - 0.4 * pi0emax);
            }
            
        }
    }
    else {
        // Dont touch the weight
    }
    

}
// -----------------------------------------------------------------------------
void selection::ApplyNuMuSelection(int type, SliceContainer &SC){

    bool pass; // A flag to see if an event passes an event

    // Classify the event
    std::pair<std::string, int> classification = SC.SliceClassifier(type);      // Classification of the event
    std::string interaction                    = SC.SliceInteractionType(type); // Genie interaction type
    std::pair<std::string, int> particle_type  = SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    if(!pass) return; // Failed the cut!
    
    // *************************************************************************
    // Common Optical Filter PE  -----------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_pe(SC, type);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // Common Optical Filter Veto  ---------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_veto(SC, type);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    if(!pass) return; // Failed the cut!

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    if(!pass) return; // Failed the cut!
        
    // *************************************************************************
    // Topological Score -------------------------------------------------------
    // *************************************************************************
    if (SC.topological_score < 0.5) return;
    
    // *************************************************************************
    // Slice Contained Fraction ------------------------------------------------
    // *************************************************************************
    pass = _scuts.contained_frac(SC);
    if(!pass) return; // Failed the cut!

    // Cuts to kill the shower backgrounds and EXT
    // Skip slices which have a reco shower
    if (SC.n_showers > 0) return;  

    // Require more than 0 tracks
    // if (SC.n_tracks <= 1) return;   

    // *************************************************************************
    // NuMu Selection Cuts ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.numu_cuts(SC);
    if(!pass) return; // Failed the cut!

    // Get the Central Value weight
    double weight = GetCVWeight(type, SC);
    
    // Also apply the pi0 weight
    GetPiZeroWeight(weight, 2, SC);

    if (!slim) _hhelper.at(type).FillNuMuHists(classification.second, SC, weight);

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection
