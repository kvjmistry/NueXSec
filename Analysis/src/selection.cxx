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
                            int num_events,
                            const char * run_period){
    
    std::cout << "\nInitialising..." << std::endl;
    if (_slim){
        std::cout << "\033[0;32m-------------------------------" << std::endl;
        std::cout << "     Running in Slim Mode!" << std::endl;
        std::cout << "-------------------------------\033[0m" << std::endl;
        slim = _slim;
    }

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_MC_POT);
        dirt_scale_factor   = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_Dirt_POT);
        intime_scale_factor = _config.at(_util.k_config_Run1_Data_trig) / _config.at(_util.k_config_Run1_EXT_trig);
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

    // Print the input files
    std::cout <<
    "Run Period Configured: run" << run_period<<"\n" <<
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
    counter_v.resize(_util.k_cuts_MAX);

    for (unsigned int t = 0; t < counter_v.size(); t++){
        counter_v.at(t).resize(_util.k_COUNTER_MAX, 0.0);
    }

    // Get MC variables --------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nInitialising MC" << std::endl;
         _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");
        //_util.GetTree(f_mc, mc_tree, "NeutrinoSelectionFilter");

        // Initialise all the mc slice container
        mc_SC.Initialise(mc_tree, _util.k_mc);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_mc).Initialise(_util.k_mc, run_period);
        if (!_slim) _hhelper.at(_util.k_mc).InitHistograms();
        
        mc_tree_total_entries = mc_tree->GetEntries();
        std::cout << "Total MC Events:         " << mc_tree_total_entries << std::endl;

        // Resize the Passed vector
        mc_passed_v.resize(mc_tree_total_entries);

        for (unsigned int y = 0; y < mc_passed_v.size(); y++ ){
            mc_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of MC Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;
    } // End getting MC variables

    // Initialise Data specific ------------------------------------------------
    if (bool_use_data){
        std::cout << "\nInitialising Data" << std::endl;
        _util.GetTree(f_data, data_tree, "nuselection/NeutrinoSelectionFilter");
        
        // Initialise all the data slice container
        data_SC.Initialise(data_tree, _util.k_data);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_data).Initialise(_util.k_data, run_period);
        if (!_slim) _hhelper.at(_util.k_data).InitHistograms();
        
        data_tree_total_entries = data_tree->GetEntries();
        std::cout << "Total Data Events:         " << data_tree_total_entries << std::endl;

        // Resize the Passed vector
        data_passed_v.resize(data_tree_total_entries);
        
        for (unsigned int y = 0; y < data_passed_v.size(); y++ ){
            data_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of Data Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of Data variables

    // Initialise EXT specific -------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nInitialising EXT" << std::endl;

        _util.GetTree(f_ext, ext_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        ext_SC.Initialise(ext_tree, _util.k_ext);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_ext).Initialise(_util.k_ext, run_period);
        if (!_slim) _hhelper.at(_util.k_ext).InitHistograms();
        
        ext_tree_total_entries = ext_tree->GetEntries();
        std::cout << "Total EXT Events:        " << ext_tree_total_entries << std::endl;

        // Resize the Passed vector
        ext_passed_v.resize(ext_tree_total_entries);

        for (unsigned int y = 0; y < ext_passed_v.size(); y++ ){
            ext_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

        std::cout << "-------------------------------" << std::endl;
        std::cout << "Initialisation of EXT Complete!" << std::endl;
        std::cout << "\033[0;31m-------------------------------\033[0m" << std::endl;

    } // End intialisation of ext variables

    // Initialise Dirt specific ------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nInitialising Dirt" << std::endl;

        _util.GetTree(f_dirt, dirt_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the data slice container
        dirt_SC.Initialise(dirt_tree, _util.k_dirt);

        // Initialise the histogram helper
        if (!_slim) _hhelper.at(_util.k_dirt).Initialise(_util.k_dirt, run_period);
        if (!_slim) _hhelper.at(_util.k_dirt).InitHistograms();
        
        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total Dirt Events:         " << dirt_tree_total_entries << std::endl;

        // Resize the Passed vector
        dirt_passed_v.resize(dirt_tree_total_entries);
        
        for (unsigned int y = 0; y < dirt_passed_v.size(); y++ ){
            dirt_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

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
    
    int counter = 0;

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

            // std::cout << mc_SC.run << " " << mc_SC.sub<<" " << mc_SC.evt<<  std::endl;
            // Get the entry in the tree
            mc_tree->GetEntry(ievent); 

            // Apply the selection cuts 
            bool pass = ApplyCuts(_util.k_mc, ievent, counter_v, mc_passed_v, mc_SC);
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
            data_tree->GetEntry(ievent);

            bool pass = ApplyCuts(_util.k_data, ievent, counter_v, data_passed_v, data_SC);
            if (!pass) continue;
        }

        // Make the run subrun event file list after the selection 
        if (make_list){

            std::cout << "Making the run subrun list for selected data events..." << std::endl;

            int run, subrun, event;
            std::ofstream run_subrun_file;
            run_subrun_file.open("files/run_subrun_list_data.txt");

            // Loop over the data events and make a run_subrun_event filelist
            for (int ievent = 0; ievent < data_tree_total_entries; ievent++){
            
                if (data_passed_v.at(ievent).cut_v.at(_util.k_cuts_MAX - 1 ) == true ){
                    
                    // Get the entry in the tree
                    data_tree->GetEntry(ievent);

                    run    = data_SC.run;
                    subrun = data_SC.sub;
                    event  = data_SC.evt;

                    run_subrun_file << run << " " << subrun << " " << event << '\n';

                }
                
            }

            run_subrun_file.close();
        }
        
        std::cout << std::endl;
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "Starting Selection over EXT" << std::endl;

        for (int ievent = 0; ievent < ext_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            ext_tree->GetEntry(ievent); // TPC Objects

            bool pass = ApplyCuts(_util.k_ext, ievent, counter_v, ext_passed_v, ext_SC);
            if (!pass) continue;
        }
         
        std::cout << std::endl;
        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){
        std::cout << "Starting Selection over Dirt" << std::endl;

        for (int ievent = 0; ievent < dirt_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            dirt_tree->GetEntry(ievent);

            bool pass = ApplyCuts(_util.k_dirt, ievent, counter_v, dirt_passed_v, dirt_SC);
            if (!pass) continue;
        }
         
        std::cout << std::endl;
        std::cout << "Ending Selection over Dirt" << std::endl;

    }
    // -------------------------------------------------------------------------
    std::cout << "Finished running the selection!"<< std::endl;

    // Print information from the selection -- need a loop for all cuts!
    for (unsigned int p=0; p < counter_v.size();p++){
        _util.PrintInfo(counter_v.at(p), intime_scale_factor, mc_scale_factor, dirt_scale_factor, _util.cut_dirs.at(p), counter_v.at(_util.k_unselected).at(_util.k_count_nue_cc));
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
    std::string category                       = SC.SliceCategory();            // The pandora group slice category

    // *************************************************************************
    // Unselected---------------------------------------------------------------
    // *************************************************************************
    
    // Fill reconstructed histograms
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_unselected, SC);
    
    // Tabulate the selection
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_unselected) );

    // Fill the hists for making the efficiency plot
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_unselected, classification.first, SC);
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    passed_v.at(ievent).cut_v.at(_util.k_swtrig) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_swtrig, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_swtrig) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_swtrig, classification.first, SC);

    // *************************************************************************
    // Op Filt PE -- MC Only ---------------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_pe(SC, type);
    passed_v.at(ievent).cut_v.at(_util.k_opfilt_pe) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_opfilt_pe, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_opfilt_pe) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_opfilt_pe, classification.first, SC);

    // *************************************************************************
    // Op Filt Michel Veto -- MC Only ------------------------------------------
    // *************************************************************************
    pass = _scuts.opfilt_veto(SC, type);
    passed_v.at(ievent).cut_v.at(_util.k_opfilt_veto) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_opfilt_veto, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_opfilt_veto) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_opfilt_veto, classification.first, SC);

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    passed_v.at(ievent).cut_v.at(_util.k_slice_id) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_slice_id, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_slice_id) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_slice_id, classification.first, SC);
    
    // *************************************************************************
    // Electron Candidate ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.e_candidate(SC);
    passed_v.at(ievent).cut_v.at(_util.k_e_candidate) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_e_candidate, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_e_candidate) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_e_candidate, classification.first, SC);

    // *************************************************************************
    // Topological Score -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.topo_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_topo_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_topo_score, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_topo_score) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_topo_score, classification.first, SC);

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    passed_v.at(ievent).cut_v.at(_util.k_in_fv) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_in_fv, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_in_fv) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_in_fv, classification.first, SC);

    // *************************************************************************
    // Cluster Fraction --------------------------------------------------------
    // *************************************************************************
    pass = _scuts.cluster_frac(SC);
    passed_v.at(ievent).cut_v.at(_util.k_cluster_frac) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_cluster_frac, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_cluster_frac) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_cluster_frac, classification.first, SC);

    // *************************************************************************
    // Shower Score ------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shower_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shower_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_shower_score, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_shower_score) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_shower_score, classification.first, SC);

    // *************************************************************************
    // Michel Rejection --------------------------------------------------------
    // *************************************************************************
    pass = _scuts.michel_rej(SC);
    passed_v.at(ievent).cut_v.at(_util.k_michel_rej) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_michel_rej, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_michel_rej) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_michel_rej, classification.first, SC);

    // *************************************************************************
    // dEdx --------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.dEdx(SC);
    passed_v.at(ievent).cut_v.at(_util.k_dEdx) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_dEdx, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_dEdx) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_dEdx, classification.first, SC);

    // *************************************************************************
    // Selected ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.selected(SC);
    passed_v.at(ievent).cut_v.at(_util.k_selected) = pass;
    if(!pass) return false; // Failed the cut!
    
    if (!slim) _hhelper.at(type).FillReco(type, classification.second, _util.k_selected, SC);
    _util.Tabulate(interaction, classification.first, type, counter_v.at(_util.k_selected) );
    if (!slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(_util.k_selected, classification.first, SC);
    
    
    // // *************************************************************************
    return true;

}
// -----------------------------------------------------------------------------
void selection::SavetoFile(){

    // Now saving histograms to file
    std::cout << "Now Saving Histograms to file" << std::endl;
    if (bool_use_mc) {
        // _hhelper.WriteMCTruth("MC");
        _hhelper.at(_util.k_mc).WriteTEfficiency();
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
        _hhelper.at(_util.k_dirt).WriteReco(_util.k_dirt);
    }

} // End save to file
// -----------------------------------------------------------------------------
void selection::MakeHistograms(const char * hist_file_name, const char *run_period, const std::vector<double> _config){
    std::cout << "Creating histograms and making plots" << std::endl;
    
    double Data_POT;

    // Set the scale factors
    if (strcmp(run_period, "1") == 0){
        mc_scale_factor     = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_MC_POT) ;
        dirt_scale_factor   = _config.at(_util.k_config_Run1_Data_POT)  / _config.at(_util.k_config_Run1_Dirt_POT);
        intime_scale_factor = _config.at(_util.k_config_Run1_Data_trig) / _config.at(_util.k_config_Run1_EXT_trig);
        Data_POT = _config.at(_util.k_config_Run1_Data_POT); // Define this variable here for easier reading
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

    histogram_plotter _hplot; // One for each type e.g MC, Data, EXT..

    _hplot.Initalise(hist_file_name, run_period, mc_scale_factor, intime_scale_factor, dirt_scale_factor);

    // Loop over the cuts and plot histograms by plot type
    for (unsigned int i = 0 ; i < _util.k_cuts_MAX; i++){
        
        // Create a set of strings for creating a dynamic directory
        // Directory structure that is created will take the form plots/<cut>/
        std::string a = "if [ ! -d \"plots/";
        std::string b = "run" + std::string(run_period) + "/" + _util.cut_dirs.at(i);
        std::string c = "\" ]; then echo \"\nPlots folder does not exist... creating\"; mkdir -p plots/";
        std::string d = "run" + std::string(run_period) + "/" + _util.cut_dirs.at(i);
        std::string e = "; fi";
        std::string command = a + b + c + d + e ;
        system(command.c_str()); 

        // MakeStack(std::string hist_name, std::string cut_name, bool area_norm, bool logy, const char* x_axis_name, double y_scale_factor, 
        //                             const double leg_x1, const double leg_x2, const double leg_y1, const double leg_y2, const char* print_name )

        // Reco X
        _hplot.MakeStack("h_reco_vtx_x",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex X [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_x.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
        
        // Reco Y
         _hplot.MakeStack("h_reco_vtx_y",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex Y [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_y.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Reco Z
        _hplot.MakeStack("h_reco_vtx_z",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex Z [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_z.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Reco X SCE
        _hplot.MakeStack("h_reco_vtx_x_sce",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex X (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_x_sce.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
        
        // Reco Y SCE
         _hplot.MakeStack("h_reco_vtx_y_sce",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex Y (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_y_sce.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Reco Z SCE
        _hplot.MakeStack("h_reco_vtx_z_sce",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Reco Vertex Z (Space Charge Corr) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_vtx_z_sce.pdf", run_period, _util.cut_dirs.at(i).c_str()) );


        // dEdx
        _hplot.MakeStack("h_reco_dEdx_y_plane",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Collection Plane dEdx (uncalibrated) [MeV/cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_dEdx_y_plane.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // dEdx cali
        _hplot.MakeStack("h_reco_dEdx_cali_y_plane",_util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Collection Plane dEdx (calibrated) [MeV/cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_dEdx_cali_y_plane.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
   
        // Leading Shower Momentum
        _hplot.MakeStack("h_reco_leading_mom", _util.cut_dirs.at(i).c_str(),
                         false, false, 1.0, "Leading Shower Momentum [MeV/c]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_mom.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // 2D distance largest flash to reco nu vertex
        _hplot.MakeStack("h_reco_flash_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Flash to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_flash_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // 2D distance shower vertex to reco nu vertex
        _hplot.MakeStack("h_reco_shower_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Shower to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // 2D distance track vertex to reco nu vertex
        _hplot.MakeStack("h_reco_trac_util.k_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Track to Vertex Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_trac_util.k_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading Shower hits in all planes
        _hplot.MakeStack("h_reco_leading_shower_hits_all_planes", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Hits All Planes", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_hits_all_planes.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading Shower hits in collection
        _hplot.MakeStack("h_reco_leading_shower_hits_collection_plane", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Hits Collection Plane", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_hits_collection_plane.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading Shower opening angle
        _hplot.MakeStack("h_reco_leading_shower_open_angle", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Open Angle [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_open_angle.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Secondary shower to vertex distance (for events with more than 1 shower)
        _hplot.MakeStack("h_reco_secondary_shower_to_vtx_dist", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Secondary Shower to Vertex Distance (>1 shower) [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_secondary_shower_to_vtx_dist.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading Shower hits per length
        _hplot.MakeStack("h_reco_leading_shower_hits_per_length", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Hits / Length [cm^{-1}]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_hits_per_length.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Longest track to leading shower length
        _hplot.MakeStack("h_reco_longest_trac_util.k_leading_shower_length", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Longest Track Length / Leading Shower Length", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_longest_track_leading_shower_length.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Track Containment
        _hplot.MakeStack("h_reco_track_contained", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Contained Tracks", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_track_contained.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading shower phi
        _hplot.MakeStack("h_reco_leading_shower_phi", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Phi [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_phi.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading shower theta
        _hplot.MakeStack("h_reco_leading_shower_theta", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Theta [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_theta.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading shower cos theta
        _hplot.MakeStack("h_reco_leading_shower_cos_theta", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Leading Shower Cos(#theta)", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_leading_shower_cos_theta.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading shower multiplicity
        _hplot.MakeStack("h_reco_shower_multiplicity", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Shower Multiplicty", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_multiplicity.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Leading track multiplicity
        _hplot.MakeStack("h_reco_track_multiplicity", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Track Multiplicty", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_track_multiplicity.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Topological Score
        _hplot.MakeStack("h_reco_topological_score", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Topological Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_topological_score.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Track shower dist
        _hplot.MakeStack("h_reco_track_shower_dist", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Longest Track Leading Shower Distance [cm]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_track_shower_dist.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
        
        // Track shower angle
        _hplot.MakeStack("h_reco_track_shower_angle", _util.cut_dirs.at(i).c_str(),
                         false,  true, 1.0, "Longest Track Leading Shower Angle [degrees]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_track_shower_angle.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
   
        // Ratio hits from showers to slice
        _hplot.MakeStack("h_reco_hits_ratio", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Hit Ratio of all Showers and the Slice", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_hits_ratio.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
        
         // Shower score
        _hplot.MakeStack("h_reco_shower_score", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Shower Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_score.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Track score
        _hplot.MakeStack("h_reco_track_score", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Track Score", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_track_score.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Calibrated energy of all the showers
        _hplot.MakeStack("h_reco_shower_energy_tot_cali", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Total Calibrated Energy of all Showers [GeV]", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_energy_tot_cali.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        // Total number of hits for the leading showe
        _hplot.MakeStack("h_reco_shower_hits", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Total Num of hits for the leading Shower", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_hits.pdf", run_period, _util.cut_dirs.at(i).c_str()) );

        
        // Total number of hits for the leading shower in the collection plane
        _hplot.MakeStack("h_reco_shower_hits_y_plane", _util.cut_dirs.at(i).c_str(),
                         false,  false, 1.0, "Total Num of hits for the leading Shower in Collection Plane", 0.8, 0.98, 0.87, 0.32, Data_POT,
                         Form("plots/run%s/%s/reco_shower_hits_y_plane.pdf", run_period, _util.cut_dirs.at(i).c_str()) );
        

    }
    

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection
