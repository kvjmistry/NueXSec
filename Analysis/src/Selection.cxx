#include "../include/Selection.h"

namespace xsecSelection {
// -----------------------------------------------------------------------------
void Selection::Initialise(Utility _utility){
    
    std::cout << "\nInitialising..." << std::endl;

    _util = _utility; 

    // Initialise the Selection cuts class
    _scuts.Initalise(_util);

    // Display slimmed Selection
    if (_util.slim){
        std::cout << "\033[0;32m-------------------------------" << std::endl;
        std::cout << "     Running in Slim Mode!" << std::endl;
        std::cout << "-------------------------------\033[0m" << std::endl;
    }

    // Create the file directory if it does not exist already
    gSystem->Exec("if [ ! -d \"files/trees/\" ]; then echo \"\nfiles folder does not exist... creating\"; mkdir -p files/trees; fi"); 

    // Set the maximum number of events tp process
    if (_util.num_events > 0 ) max_events = _util.num_events;

    // Resize the histogram helper instance vectors, one instance per type e.g MC, data, ..
    _hhelper.resize(_util.k_type_MAX);
    
    // Resize the tree helper instance vectors, one instance per type e.g MC, data, ..
    _thelper.resize(_util.k_type_MAX);

    // Print the input files
    std::cout <<
    "Run Period Configured: run" << _util.run_period<<"\n" <<
    "MC   File Path:      " << _util.mc_file_name        <<"\n" <<
    "EXT  File Path:      " << _util.ext_file_name       <<"\n" <<
    "Data File Path:      " << _util.data_file_name      <<"\n" <<
    "Dirt File Path:      " << _util.dirt_file_name      <<"\n" <<
    std::endl;

    // Now get the files, if file isnt specified then set bool to skip
    bool_use_mc        = _util.GetFile(f_mc,        _util.mc_file_name);
    bool_use_ext       = _util.GetFile(f_ext,       _util.ext_file_name);
    bool_use_data      = _util.GetFile(f_data,      _util.data_file_name);
    bool_use_dirt      = _util.GetFile(f_dirt,      _util.dirt_file_name);

    // Resize the counter vector
    counter_v.resize(_util.k_cuts_MAX);

    for (unsigned int t = 0; t < counter_v.size(); t++){
        counter_v.at(t).resize(_util.k_COUNTER_MAX, 0.0);
    }

    // Get MC variables --------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nInitialising MC" << std::endl;
        _util.GetTree(f_mc, mc_tree, "NeutrinoSelectionFilter");

        // Initialise all the mc slice container
        mc_SC.Initialise(mc_tree, _util.k_mc, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_mc).Initialise(_util.k_mc, _util.mc_tree_file_name_out, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_mc).Initialise(_util.k_mc, _util.mc_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_mc).InitHistograms();

        mc_tree_total_entries = mc_tree->GetEntries();
        std::cout << "Total MC Events:         " << mc_tree_total_entries << std::endl;

        std::cout << "Initialisation of MC Complete!" << std::endl;
    } // End getting MC variables

    // Initialise Data specific ------------------------------------------------
    if (bool_use_data){
        std::cout << "\nInitialising Data" << std::endl;
        _util.GetTree(f_data, data_tree, "NeutrinoSelectionFilter");
        
        // Initialise all the data slice container
        data_SC.Initialise(data_tree, _util.k_data, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_data).Initialise(_util.k_data, _util.data_tree_file_name_out, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_data).Initialise(_util.k_data, _util.data_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_data).InitHistograms();
        
        data_tree_total_entries = data_tree->GetEntries();
        std::cout << "Total Data Events:         " << data_tree_total_entries << std::endl;

        std::cout << "Initialisation of Data Complete!" << std::endl;

    } // End intialisation of Data variables

    // Initialise EXT specific -------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nInitialising EXT" << std::endl;

        _util.GetTree(f_ext, ext_tree, "NeutrinoSelectionFilter");

        // Initialise all the data slice container
        ext_SC.Initialise(ext_tree, _util.k_ext, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_ext).Initialise(_util.k_ext, _util.ext_tree_file_name_out, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_ext).Initialise(_util.k_ext, _util.ext_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_ext).InitHistograms();
        
        ext_tree_total_entries = ext_tree->GetEntries();
        std::cout << "Total EXT Events:        " << ext_tree_total_entries << std::endl;

        std::cout << "Initialisation of EXT Complete!" << std::endl;

    } // End intialisation of ext variables

    // Initialise Dirt specific ------------------------------------------------
    if (bool_use_dirt){
        std::cout << "\nInitialising Dirt" << std::endl;

        _util.GetTree(f_dirt, dirt_tree, "NeutrinoSelectionFilter");

        // Initialise all the data slice container
        dirt_SC.Initialise(dirt_tree, _util.k_dirt, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_dirt).Initialise(_util.k_dirt, _util.dirt_tree_file_name_out, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_dirt).Initialise(_util.k_dirt, _util.dirt_file_name_out ,_util);
        if (!_util.slim) _hhelper.at(_util.k_dirt).InitHistograms();
        
        

        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total Dirt Events:         " << dirt_tree_total_entries << std::endl;

        std::cout << "Initialisation of Dirt Complete!" << std::endl;


    } // End intialisation of dirt variables    
    
    // Invoke main Selection function
    MakeSelection();

} // END Initialise function
// -----------------------------------------------------------------------------
// Main function for Selection
void Selection::MakeSelection(){
    std::cout << "\n\033[0;32mNow Running the Selection!\033[0m"<< std::endl;
    
    // MC ----------------------------------------------------------------------
    if (bool_use_mc){
        std::cout << "\nStarting Selection over MC" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_mc;
        run_subrun_file_mc.open(Form("files/txt/run%s_run_subrun_list_mc.txt",_util.run_period));

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

            // if (ievent != 578)
            //     continue;

            // std::cout << ievent << std::endl;

            // std::cout << mc_SC.run << " " << mc_SC.sub<<" " << mc_SC.evt<<  std::endl;

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_mc, mc_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_mc, mc_SC);
            
            // Apply the Selection cuts 
            bool pass = ApplyCuts(_util.k_mc, counter_v, mc_SC);

            // Only fill passed events for fake data
            // if (_util.isfakedata){
                
            //     // if (pass)
            //         // _thelper.at(_util.k_mc).FillVars(mc_SC, pass);
            // }
            // // Fill the output tree if the event passed or it was signal
            // else {
            //     // if (pass || mc_SC.is_signal)
            //         // _thelper.at(_util.k_mc).FillVars(mc_SC, pass);
            // }
            
            // // If the event passed the selection then save the run subrun event to file
            // if (pass) run_subrun_file_mc << mc_SC.run << " " << mc_SC.sub << " " << mc_SC.evt << '\n';

        } // End Event loop

        run_subrun_file_mc.close();
        std::cout << "Ending Selection over MC" << std::endl;

    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){
        std::cout << "\nStarting Selection over Data" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_data;
        run_subrun_file_data.open(Form("files/txt/run%s_run_subrun_list_data.txt",_util.run_period));

        for (int ievent = 0; ievent < data_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            data_tree->GetEntry(ievent);

            // Skip the events with different sw trigger configured
            if (std::string(_util.run_period) == "3" && data_SC.run > 16880 ){
            // if (_run_period == 3 && data_SC.run < 16880 ){
                // continue;
            }

            // Look at different regions of run 1
            if (std::string(_util.run_period) == "1" ){

                // 4p6
                // if ((data_SC.run > 6035 && data_SC.run < 6284) || data_SC.run > 6510 ){
                //     continue;
                // }
                // 6p6
                // if ((data_SC.run > 6285 && data_SC.run < 6510) || data_SC.run < 6035 ){
                //     continue;
                // }

            }

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_data, data_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_data, data_SC);

            bool pass = ApplyCuts(_util.k_data, counter_v, data_SC);

            // Fill the output tree if the event passed the selection
            // if (pass) _thelper.at(_util.k_data).FillVars(data_SC, pass);
            
            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_data << data_SC.run << " " << data_SC.sub << " " << data_SC.evt << '\n';

        }

        run_subrun_file_data.close();
        
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nStarting Selection over EXT" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_ext;
        run_subrun_file_ext.open(Form("files/txt/run%s_run_subrun_list_ext.txt",_util.run_period));

        for (int ievent = 0; ievent < ext_tree_total_entries; ievent++){

            // See if we want to process all the events
            if (max_events > 0){
                if (ievent >= max_events) break;
            }

            // Alert the user
            if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;
        
            // Get the entry in the tree
            ext_tree->GetEntry(ievent); // TPC Objects

            // Skip the RHC events contaminated in the FHC files
            if (std::string(_util.run_period) == "3" && ext_SC.run > 16880 ){
            // if (_run_period == 3 && ext_SC.run < 16880 ){
                // continue;
            }

            // Look at different regions of run 1
            if (std::string(_util.run_period) == "1" ){

                // Set 1
                // if (ext_SC.run >= 6550) continue;

                // Set2
                // if (ext_SC.run < 6550 || ext_SC.run > 7013) continue;

                // Set 3
                // if (ext_SC.run < 7013) continue;
            }

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_ext, ext_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_ext, ext_SC);

            bool pass = ApplyCuts(_util.k_ext, counter_v, ext_SC);
            
            // Fill the output tree if the event passed the selection
            // if (pass) _thelper.at(_util.k_ext).FillVars(ext_SC, pass);

            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_ext << ext_SC.run << " " << ext_SC.sub << " " << ext_SC.evt << '\n';

        }
         
        run_subrun_file_ext.close();

        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){

        std::cout << "\nStarting Selection over Dirt" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_dirt;
        run_subrun_file_dirt.open(Form("files/txt/run%s_run_subrun_list_dirt.txt",_util.run_period));

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

            bool pass = ApplyCuts(_util.k_dirt, counter_v, dirt_SC);

            // Fill the output tree if the event passed the selection
            // if (pass) _thelper.at(_util.k_dirt).FillVars(dirt_SC, pass);
            
            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_dirt << dirt_SC.run << " " << dirt_SC.sub << " " << dirt_SC.evt << '\n';

        }

        run_subrun_file_dirt.close();
         
        std::cout << "Ending Selection over Dirt" << std::endl;

    }
    // -------------------------------------------------------------------------
    std::cout << "Finished running the Selection!"<< std::endl;

    // Save information from the Selection so we can print
    for (unsigned int p=0; p < counter_v.size();p++){

        // Fill the counter trees
        if (bool_use_mc)   _thelper.at(_util.k_mc)  .Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_data) _thelper.at(_util.k_data).Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_ext)  _thelper.at(_util.k_ext) .Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);
        if (bool_use_dirt) _thelper.at(_util.k_dirt).Fill_counters(counter_v.at(p), bool_use_mc, bool_use_ext, bool_use_data, bool_use_dirt);

    }
    
    // Now save all the outputs to file
    if (!_util.slim) SavetoFile();

    return;
} // End Selection
// -----------------------------------------------------------------------------
bool Selection::ApplyCuts(int type,std::vector<std::vector<double>> &counter_v, SliceContainer &SC){

    // Here we apply the Selection cuts ----------------------------------------
    bool pass; // A flag to see if an event passes an event

    if (SC.tpc_obj_index == 0) {
        SC.already_filled.clear(); // set this to false every time when we have a new tpc object
        SC.already_filled.resize(_util.k_cuts_MAX);
        std::fill(SC.already_filled.begin(), SC.already_filled.end(), false);

    }


    // Classify the event -- sets variable in the slice contianer
    SC.SliceClassifier(type);      // Classification of the event

    // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
    SC.SetThresholdEvent();
    
    SC.SliceInteractionType(type); // Genie interaction type
    SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
    SC.Pi0Classifier(type); 

    // Set derived variables in the slice container
    SC.SetSignal();                // Set the event as either signal or other
    SC.SetFakeData();              // Set the classifcation as data if fake data mode
    SC.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
    SC.SetNuMIAngularVariables();  // Set the NuMI angular variables
    SC.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

    // *************************************************************************
    // Unselected---------------------------------------------------------------
    // *************************************************************************
    SelectionFill(type, SC, _util.k_unselected, counter_v );
    
    // Let apply some cuts and see what happens
    // *************************************************************************
    // In Time-----------------------------------------------------------------
    // *************************************************************************
    // if ((classification.first == "nue_cc" || classification.first == "nuebar_cc") && (SC.flash_time > 16.0 || SC.flash_time < 5.5 )) std::cout << SC.flash_time<< std::endl;
    if (SC.failed_intime_cut == 1) return false;
    SelectionFill(type, SC, _util.k_in_time, counter_v );  
    
    // *************************************************************************
    // Flash PE ----------------------------------------------------------------
    // *************************************************************************
    if (SC.failed_flash_cut == 1) return false;
    SelectionFill(type, SC, _util.k_flash_pe, counter_v );
    
    // *************************************************************************
    // Reco Nue-----------------------------------------------------------------
    // *************************************************************************
    if (SC.has_valid_shr == 0 || SC.has_nue == 0) return false;
    SelectionFill(type, SC, _util.k_reco_nue, counter_v );

    // *************************************************************************
    // IN FV--------------------------------------------------------------------
    // *************************************************************************
    if (_util.in_fv(SC.reco_nu_vtx_x, SC.reco_nu_vtx_y, SC.reco_nu_vtx_z) == false) return false;
    SelectionFill(type, SC, _util.k_in_fv, counter_v );
    
    
    // *************************************************************************
    // Vertex to Flash ---------------------------------------------------------
    // *************************************************************************
    if (SC.reco_nu_vtx_z > SC.flash_vtx_z){
        if (SC.flash_vtx_dist > 60 ) return false;
    }
    else {
        if (SC.flash_vtx_dist > 80 ) return false;
    }
    SelectionFill(type, SC, _util.k_vtx_to_flash, counter_v );
    
    
    // *************************************************************************
    // Shower to vertex --------------------------------------------------------
    // *************************************************************************
    bool valid_shr = false;
    if (SC.shr_vtx_dist_v->size() == 0) valid_shr= true;
    for (unsigned int u = 0; u < SC.shr_vtx_dist_v->size(); u++){
        if (SC.shr_vtx_dist_v->at(u) < 4) valid_shr = true;
    }
    if (!valid_shr) return valid_shr;
    SelectionFill(type, SC, _util.k_shr_to_vtx, counter_v );
    
    
    // *************************************************************************
    // Track to Vertex----------------------------------------------------------
    // *************************************************************************
    bool valid_track = false;
    if (SC.trk_vtx_dist_v->size() == 0) valid_track= true;
    for (unsigned int u = 0; u < SC.trk_vtx_dist_v->size(); u++){
        if (SC.trk_vtx_dist_v->at(u) < 4) valid_track = true;
    }
    if (!valid_track) return valid_track;
    
    SelectionFill(type, SC, _util.k_track_to_vtx, counter_v );
    
    
    // *************************************************************************
    // Hit Threshold -----------------------------------------------------------
    // *************************************************************************
    if (SC.shr_hits_tot < 200) return false; 
    SelectionFill(type, SC, _util.k_hit_thresh, counter_v );
    
    
    // *************************************************************************
    // Hit Threshold Collection ------------------------------------------------
    // *************************************************************************
    if (SC.shr_hits_y_tot < 80) return false; 
    SelectionFill(type, SC, _util.k_hit_thresh_y, counter_v );
    
    
    // *************************************************************************
    // Open Angle --------------------------------------------------------------
    // *************************************************************************
    if (SC.shr_openangle < 2 || SC.shr_openangle > 15) return false; 
    SelectionFill(type, SC, _util.k_open_angle, counter_v );
    
    // *************************************************************************
    // Shower dEdx -------------------------------------------------------------
    // *************************************************************************
    if (SC.shr_dedx_Y < 1.4 || SC.shr_dedx_Y > 3) return false; 
    SelectionFill(type, SC, _util.k_dedx, counter_v );

    // *************************************************************************
    // Shr Vertex Distance > 1 shower ------------------------------------------
    // *************************************************************************
    for (unsigned int u = 0; u < SC.sec_shr_vtx_dist_v->size(); u++){
        if (SC.sec_shr_vtx_dist_v->at(u) > 22) return false;
    }
    SelectionFill(type, SC, _util.k_shr_vtx_dist_gt_1shr, counter_v );
    
    
    // *************************************************************************
    // Hits per length ---------------------------------------------------------------
    // *************************************************************************
    if (double(SC.shr_hits_tot/SC.shr_len) < 3.0) return false;
    SelectionFill(type, SC, _util.k_hit_per_lenth, counter_v );
    
    
    // *************************************************************************
    // Track Shower Length Ratio -----------------------------------------------
    // *************************************************************************
    double track_shr_ratio = SC.trk_len / SC.shr_len;
    // std::cout << track_shr_ratio << std::endl;
    if ( track_shr_ratio > 1.0 && SC.n_tracks > 0) return false;
    SelectionFill(type, SC, _util.k_trk_shr_lengh, counter_v );
    
    
    // *************************************************************************
    // Track Containment -------------------------------------------------------
    // *************************************************************************
    for (unsigned int u = 0; u < SC.trk_start_x_v->size(); u++){
        if ( _util.in_fv(SC.trk_start_x_v->at(u), SC.trk_start_y_v->at(u), SC.trk_start_z_v->at(u)) == false) return false;
        if ( _util.in_fv(SC.trk_end_x_v->at(u),   SC.trk_end_y_v->at(u),   SC.trk_end_z_v->at(u))   == false) return false;
    
    }
    SelectionFill(type, SC, _util.k_trk_containment, counter_v );   

    if (SC.tpc_obj_index == 0) SC.tpc_obj_counter_prev = 0; // Set back to 0
    else                       SC.tpc_obj_counter_prev++;   // Add to counter
        
    // if (( SC.nu_pdg == 12 || SC.nu_pdg == -12 ) && SC.nu_e < 0.3) std::cout << "\033[0;33m" << SC.nu_e << "\033[0m" << std::endl;

    // ************************************************************************n*
    return true;

}
// -----------------------------------------------------------------------------
void Selection::SavetoFile(){

    // Now saving histograms to file
    std::cout << "Now Saving Histograms to file" << std::endl;
    if (bool_use_mc) {
    
        if (std::string(_util.intrinsic_mode) == "default" || std::string(_util.intrinsic_mode) == "intrinsic"){
            _hhelper.at(_util.k_mc).WriteTEfficiency();
            _hhelper.at(_util.k_mc).WriteTrue();
            _hhelper.at(_util.k_mc).WriteInteractions();
            _hhelper.at(_util.k_mc).Write_2DSigBkgHists();
            
            if (_util.isfakedata){
                if (std::string(_util.intrinsic_mode) == "default"){
                    _hhelper.at(_util.k_mc).WriteReco(_util.k_data);
                }
            }
            else
                _hhelper.at(_util.k_mc).WriteReco(_util.k_mc);
        }
        
        if (std::string(_util.intrinsic_mode) == "default") {
            _hhelper.at(_util.k_mc).WriteFlash();
            
            if (_util.isfakedata){
                if (std::string(_util.intrinsic_mode) == "default"){
                    _hhelper.at(_util.k_mc).WriteRecoPar(_util.k_data);
                    _hhelper.at(_util.k_mc).WritePiZero(_util.k_data);
                    _hhelper.at(_util.k_mc).WriteNuMu(_util.k_data);
                }
            }
            else {
                _hhelper.at(_util.k_mc).WriteRecoPar(_util.k_mc);
                _hhelper.at(_util.k_mc).WritePiZero(_util.k_mc);
                _hhelper.at(_util.k_mc).WriteNuMu(_util.k_mc);
            }
        }

        _thelper.at(_util.k_mc).WriteTree(_util.k_mc);


    }
    if (bool_use_data) {
        _hhelper.at(_util.k_data).WriteReco(_util.k_data);
        _hhelper.at(_util.k_data).WriteRecoPar(_util.k_data);
        _hhelper.at(_util.k_data).WriteFlash();
        _hhelper.at(_util.k_data).WritePiZero(_util.k_data);
        _hhelper.at(_util.k_data).WriteNuMu(_util.k_data);

        _thelper.at(_util.k_data).WriteTree(_util.k_data);

    }
    if (bool_use_ext) {
        _hhelper.at(_util.k_ext).WriteReco(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteRecoPar(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteFlash();
        _hhelper.at(_util.k_ext).Write_2DSigBkgHists();
        _hhelper.at(_util.k_ext).WritePiZero(_util.k_ext);
        _hhelper.at(_util.k_ext).WriteNuMu(_util.k_ext);

        _thelper.at(_util.k_ext).WriteTree(_util.k_ext);

    }
    if (bool_use_dirt) {
        // _hhelper.at(_util.k_dirt).WriteTrue(); // Only turn this on to inspect. It wont merge since the file names are not uniique yet!!
        _hhelper.at(_util.k_dirt).WriteReco(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteRecoPar(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteFlash();
        _hhelper.at(_util.k_dirt).Write_2DSigBkgHists();
        _hhelper.at(_util.k_dirt).WritePiZero(_util.k_dirt);
        _hhelper.at(_util.k_dirt).WriteNuMu(_util.k_dirt);

        _thelper.at(_util.k_dirt).WriteTree(_util.k_dirt);

    }

} // End save to file
// -----------------------------------------------------------------------------
void Selection::SelectionFill(int type, SliceContainer &SC, int cut_index, std::vector<std::vector<double>> &counter_v){
    
    // *************************************************************************
    // Get the CV weight
    // *************************************************************************
    double weight = 1.0;

    // *************************************************************************
    // Fill Histograms
    // *************************************************************************
    // Fill almost all the histograms with this function call
    if (!_util.slim) _hhelper.at(type).FillHists(type, SC.classification.second, SC.genie_interaction, SC.particle_type.second, cut_index, SC, weight); 

    // Fill Plots for Efficiency
    if (!_util.slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(cut_index, SC.classification.first, SC, weight);

    // Fill the dedx ttree before shr dist dedx cut and after cut dedx
    // We can use this tree to play around and optimise the dedx cut. Not essential for the analysis
    // if (cut_index == _util.k_vtx_dist_dedx - 1 || cut_index == _util.k_vtx_dist_dedx){
    //     _thelper.at(type).Fill_dedxVars(SC, SC.classification, _util.cut_dirs.at(cut_index), weight);
    // }
    
    // *************************************************************************
    // Tabulate the selection i.e count everything
    // *************************************************************************
    bool is_in_fv = _util.in_fv(SC.true_nu_vtx_x, SC.true_nu_vtx_y, SC.true_nu_vtx_z); // This variable is only used in the case of MC, so it should be fine 
    bool filled = SC.already_filled.at(cut_index);
    _util.Tabulate(is_in_fv, SC.genie_interaction, SC.classification.first, SC.pi0_classification, type, counter_v.at(cut_index), weight, filled );
    SC.already_filled.at(cut_index) = filled;

}
// -----------------------------------------------------------------------------
void Selection::ApplyPiZeroSelection(int type, SliceContainer &SC){

    bool pass; // A flag to see if an event passes an event

    // Classify the event
    SC.SliceClassifier(type);      // Classification of the event

    // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
    SC.SetThresholdEvent();

    SC.SliceInteractionType(type); // Genie interaction type
    SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
        
    // *************************************************************************
    // Pi0 Selection Cuts ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.pi_zero_cuts(SC);
    if(!pass) return; // Failed the cut!

    bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

    // Get the Central Value weight
    double weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv, SC.interaction);
    double weight_norm = weight;
    double weight_Escale = weight;

    // Try scaling the pi0
    // 0 == no weighting, 1 == normalisation fix, 2 == energy dependent scaling
    _util.GetPiZeroWeight(weight,        0, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);
    _util.GetPiZeroWeight(weight_norm,   1, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);
    _util.GetPiZeroWeight(weight_Escale, 2, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);

    // Now Fill the histograms
    if (!_util.slim) _hhelper.at(type).FillPiZeroHists(SC.classification.second, SC, weight, 0);
    if (!_util.slim) _hhelper.at(type).FillPiZeroHists(SC.classification.second, SC, weight_norm, 1);
    if (!_util.slim) _hhelper.at(type).FillPiZeroHists(SC.classification.second, SC, weight_Escale, 2);


}
// -----------------------------------------------------------------------------
void Selection::ApplyNuMuSelection(int type, SliceContainer &SC){

    bool pass; // A flag to see if an event passes an event

    // Classify the event
    SC.SliceClassifier(type);      // Classification of the event

    // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
    SC.SetThresholdEvent();

    SC.SliceInteractionType(type); // Genie interaction type
    SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
    
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

    bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

    // Get the Central Value weight
    double weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv, SC.interaction);
    
    // Also apply the pi0 weight
    _util.GetPiZeroWeight(weight, _util.pi0_correction, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);

    if (!_util.slim) _hhelper.at(type).FillNuMuHists(SC.classification.second, SC, weight);

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection
