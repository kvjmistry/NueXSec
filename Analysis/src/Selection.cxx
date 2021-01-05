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
        _util.GetTree(f_mc, mc_tree, "nuselection/NeutrinoSelectionFilter");

        // Initialise all the mc slice container
        mc_SC.Initialise(mc_tree, _util.k_mc, _util);

        // Initialise the Tree Helper
        _thelper.at(_util.k_mc).Initialise(_util.k_mc, _util.mc_tree_file_name_out, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_mc).Initialise(_util.k_mc, _util.mc_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_mc).InitHistograms();

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
        data_SC.Initialise(data_tree, _util.k_data, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_data).Initialise(_util.k_data, _util.data_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_data).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_data).Initialise(_util.k_data, "empty", _util);

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
        ext_SC.Initialise(ext_tree, _util.k_ext, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_ext).Initialise(_util.k_ext, _util.ext_file_name_out, _util);
        if (!_util.slim) _hhelper.at(_util.k_ext).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_ext).Initialise(_util.k_ext, "empty", _util);

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
        dirt_SC.Initialise(dirt_tree, _util.k_dirt, _util);

        // Initialise the histogram helper
        if (!_util.slim) _hhelper.at(_util.k_dirt).Initialise(_util.k_dirt, _util.dirt_file_name_out ,_util);
        if (!_util.slim) _hhelper.at(_util.k_dirt).InitHistograms();
        
        // Initialise the Tree Helper
        _thelper.at(_util.k_dirt).Initialise(_util.k_dirt, "empty", _util);

        dirt_tree_total_entries = dirt_tree->GetEntries();
        std::cout << "Total Dirt Events:         " << dirt_tree_total_entries << std::endl;

        // Resize the Passed vector
        dirt_passed_v.resize(dirt_tree_total_entries);
        
        for (unsigned int y = 0; y < dirt_passed_v.size(); y++ ){
            dirt_passed_v.at(y).cut_v.resize(_util.k_cuts_MAX, false);
        }

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
        run_subrun_file_mc.open(Form("files/run%s_run_subrun_list_mc.txt",_util.run_period));

        // Create a file for the selected events and their properties
        std::ofstream evt_dist_sig;
        evt_dist_sig.open(Form("files/run%s_evt_dist_sig.txt",_util.run_period));
        evt_dist_sig << "true_E," << "reco_E," << "w" << "\n";

        // Create a file for the generated events and their properties
        std::ofstream evt_dist_gen;
        evt_dist_gen.open(Form("files/run%s_evt_dist_gen.txt",_util.run_period));
        evt_dist_gen << "true_E,"<< "w" << "\n";

        // Create a file for the generated events and their properties
        std::ofstream evt_dist_bkg;
        evt_dist_bkg.open(Form("files/run%s_evt_dist_bkg.txt",_util.run_period));
        evt_dist_bkg << "reco_E," << "w" << "\n";

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

            // std::cout << mc_SC.run << " " << mc_SC.sub<<" " << mc_SC.evt<<  std::endl;

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_mc, mc_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_mc, mc_SC);
            
            // Apply the Selection cuts 
            bool pass = ApplyCuts(_util.k_mc, ievent, counter_v, mc_passed_v, mc_SC);

            // Fill the output tree if the event passed or it was signal
            if (pass || mc_SC.is_signal) _thelper.at(_util.k_mc).FillVars(mc_SC, pass);
            
            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_mc << mc_SC.run << " " << mc_SC.sub << " " << mc_SC.evt << '\n';

            // if passed ans is a signal the write events to the file
            if (pass && mc_SC.is_signal && mc_SC.cv_weight != 0.0) {
                // evt_dist_sig << mc_SC.elec_e << "," << mc_SC.elec_theta << "," << mc_SC.elec_phi << "," << mc_SC.shr_energy_cali << "," << mc_SC.cv_weight*_util.mc_scale_factor << "\n";
                evt_dist_sig << mc_SC.elec_e << "," << mc_SC.shr_energy_cali << "," << mc_SC.cv_weight<< "\n";
            }

            // Generated events
            if (mc_SC.is_signal && mc_SC.cv_weight != 0.0) {
                // evt_dist_gen << mc_SC.elec_e << "," << mc_SC.elec_theta << "," << mc_SC.elec_phi << "," << mc_SC.cv_weight*_util.mc_scale_factor << "\n";
                evt_dist_gen << mc_SC.elec_e << "," << mc_SC.cv_weight << "\n";
            }

            // if passed and is a background the write events to the file
            if (pass && !mc_SC.is_signal && mc_SC.cv_weight != 0.0) {
                evt_dist_bkg << mc_SC.shr_energy_cali << "," << mc_SC.cv_weight*_util.mc_scale_factor << "\n";
            }


        } // End Event loop

        run_subrun_file_mc.close();
        evt_dist_sig.close();
        evt_dist_gen.close();
        evt_dist_bkg.close();

        std::cout << "Ending Selection over MC" << std::endl;

    }
    // Data --------------------------------------------------------------------
    if (bool_use_data){
        std::cout << "\nStarting Selection over Data" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_data;
        run_subrun_file_data.open(Form("files/run%s_run_subrun_list_data.txt",_util.run_period));

        // Create a file for the selected events
        std::ofstream evt_dist_data;
        evt_dist_data.open(Form("files/run%s_evt_dist_data.txt",_util.run_period));
        evt_dist_data << "reco_E" << "\n";


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
                // if (data_SC.run > 5900) continue;

                // 6p6
                // if (data_SC.run < 6600) continue;
            }

            // Apply Pi0 Selection
            ApplyPiZeroSelection(_util.k_data, data_SC);

            // Apply NuMu Selection
            ApplyNuMuSelection(_util.k_data, data_SC);

            bool pass = ApplyCuts(_util.k_data, ievent, counter_v, data_passed_v, data_SC);

            // Fill the output tree if the event passed the selection
            if (pass) _thelper.at(_util.k_data).FillVars(data_SC, pass);
            
            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_data << data_SC.run << " " << data_SC.sub << " " << data_SC.evt << '\n';

            if (pass) evt_dist_data << data_SC.shr_energy_cali << "\n";
        }

        run_subrun_file_data.close();
        evt_dist_data.close();
        
        std::cout << "Ending Selection over Data" << std::endl;
    }
    // EXT ---------------------------------------------------------------------
    if (bool_use_ext){
        std::cout << "\nStarting Selection over EXT" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_ext;
        run_subrun_file_ext.open(Form("files/run%s_run_subrun_list_ext.txt",_util.run_period));

        // Create a file for the bkg events and their properties
        std::ofstream evt_dist_ext;
        evt_dist_ext.open(Form("files/run%s_evt_dist_ext.txt",_util.run_period));

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

            bool pass = ApplyCuts(_util.k_ext, ievent, counter_v, ext_passed_v, ext_SC);
            
            // Fill the output tree if the event passed the selection
            if (pass) _thelper.at(_util.k_ext).FillVars(ext_SC, pass);

            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_ext << ext_SC.run << " " << ext_SC.sub << " " << ext_SC.evt << '\n';

            // if passed and is a background the write events to the file
            if (pass) {
                evt_dist_ext << ext_SC.shr_energy_cali << "," << ext_SC.cv_weight*_util.ext_scale_factor << "\n";
            }
        }
         
        run_subrun_file_ext.close();
        evt_dist_ext.close();

        std::cout << "Ending Selection over EXT" << std::endl;

    }
    // Dirt --------------------------------------------------------------------
    if (bool_use_dirt){

        std::cout << "\nStarting Selection over Dirt" << std::endl;

        // Create file for saving run, subrun event
        std::ofstream run_subrun_file_dirt;
        run_subrun_file_dirt.open(Form("files/run%s_run_subrun_list_dirt.txt",_util.run_period));

        // Create a file for the bkg events and their properties
        std::ofstream evt_dist_dirt;
        evt_dist_dirt.open(Form("files/run%s_evt_dist_dirt.txt",_util.run_period));

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

            // Fill the output tree if the event passed the selection
            if (pass) _thelper.at(_util.k_dirt).FillVars(dirt_SC, pass);
            
            // If the event passed the selection then save the run subrun event to file
            if (pass) run_subrun_file_dirt << dirt_SC.run << " " << dirt_SC.sub << " " << dirt_SC.evt << '\n';

            // if passed and is a background the write events to the file
            if (pass) {
                evt_dist_dirt << dirt_SC.shr_energy_cali << "," << dirt_SC.cv_weight*_util.dirt_scale_factor << "\n";
            }
        }

        run_subrun_file_dirt.close();
        evt_dist_dirt.close();
         
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
bool Selection::ApplyCuts(int type, int ievent,std::vector<std::vector<double>> &counter_v,
                           std::vector<PassedContainer> &passed_v, SliceContainer &SC){

    // Here we apply the Selection cuts ----------------------------------------
    bool pass; // A flag to see if an event passes an event

    // Classify the event -- sets variable in the slice contianer
    SC.SliceClassifier(type);      // Classification of the event
    SC.SliceInteractionType(type); // Genie interaction type
    SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
    SC.Pi0Classifier(type); 

    // Set derived variables in the slice container
    SC.SetSignal();                // Set the event as either signal or other
    SC.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
    SC.SetNuMIAngularVariables();  // Set the NuMI angular variables
    SC.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

    // *************************************************************************
    // Unselected---------------------------------------------------------------
    // *************************************************************************
    SelectionFill(type, SC, _util.k_unselected, counter_v );
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    passed_v.at(ievent).cut_v.at(_util.k_swtrig) = pass;
    // if (type == _util.k_mc && !pass && classification.first == "nu_out_fv") std::cout<< SC.run << " " << SC.sub << " " << SC.evt << std::endl; // Spit out the run subrun event numbers of failed events for an event display
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_swtrig, counter_v );

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    passed_v.at(ievent).cut_v.at(_util.k_slice_id) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_slice_id, counter_v );
    
    // *************************************************************************
    // Electron Candidate ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.e_candidate(SC);
    passed_v.at(ievent).cut_v.at(_util.k_e_candidate) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_e_candidate, counter_v );

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    passed_v.at(ievent).cut_v.at(_util.k_in_fv) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_in_fv, counter_v );
    
    // *************************************************************************
    // Slice Contained Fraction ------------------------------------------------
    // *************************************************************************
    pass = _scuts.contained_frac(SC);
    passed_v.at(ievent).cut_v.at(_util.k_contained_frac) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_contained_frac, counter_v );

    // *************************************************************************
    // Topological Score -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.topo_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_topo_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_topo_score, counter_v );

    // *************************************************************************
    // Cosmic Impact Parameter -------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_cosmic_IP(SC);
    passed_v.at(ievent).cut_v.at(_util.k_cosmic_ip) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_cosmic_ip, counter_v );

    // *************************************************************************
    // Shower Score ------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shower_score(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shower_score) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_shower_score, counter_v );

    // *************************************************************************
    // Shower Hit Ratio  -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_hitratio(SC);
    passed_v.at(ievent).cut_v.at(_util.k_hit_ratio) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_hit_ratio, counter_v );

    // *************************************************************************
    // Shower Moliere Average --------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_moliere_avg(SC);
    passed_v.at(ievent).cut_v.at(_util.k_shr_moliere_avg) = pass;

    // if (!pass && SC.shrmoliereavg >= 7) std::cout << SC.run << " " << SC.sub << " " << SC.evt << std::endl;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_shr_moliere_avg, counter_v );

    // *************************************************************************
    // 2D cut for Shower to Vertex Distance and dEdx ---------------------------
    // *************************************************************************
    pass = _scuts.shr_dist_dEdx_max(SC);
    passed_v.at(ievent).cut_v.at(_util.k_vtx_dist_dedx) = pass;
    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_vtx_dist_dedx, counter_v );

    // *************************************************************************
    // dEdx in all planes for 0 track events -----------------------------------
    // *************************************************************************
    pass = _scuts.dEdx_max_no_tracks(SC);
    passed_v.at(ievent).cut_v.at(_util.k_dEdx_max_no_tracks) = pass;

    if(!pass) return false; // Failed the cut!
    
    SelectionFill(type, SC, _util.k_dEdx_max_no_tracks, counter_v );

    // Skip unnaturally high shower energies?
    // if (SC.shr_energy_cali > 6.0)
    //     std::cout << "reco shr energy: "  << SC.shr_energy_cali << "  elec E: "<< SC.elec_e << "  nu E: " << SC.nu_e<< std::endl;

    // if (SC.is_signal && SC.nu_e < 0.3) std::cout<<"Low elec E!: " <<SC.elec_e*1000 << " MeV" << "  | E Nu: "<< SC.nu_e*1000 << " MeV" <<  std::endl; 
    // if (type == _util.k_mc && (SC.nu_pdg == -12 || SC.nu_pdg == 12) && SC.nu_e <= 0.125 && SC.ccnc == _util.k_NC) std::cout << "Got nue NC event selected thats below th: " << SC.nu_e << std::endl;

    // **************************************************************************
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
            _hhelper.at(_util.k_mc).WriteReco(_util.k_mc);
            _hhelper.at(_util.k_mc).Write_2DSigBkgHists();
        }
        
        if (std::string(_util.intrinsic_mode) == "default") {
            _hhelper.at(_util.k_mc).WriteRecoPar(_util.k_mc);
            _hhelper.at(_util.k_mc).WriteFlash();
            _hhelper.at(_util.k_mc).WritePiZero(_util.k_mc);
            _hhelper.at(_util.k_mc).WriteNuMu(_util.k_mc);
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
    bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 
    weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv);

    // Try scaling the pi0 -- need to implement this as a configurable option
    // 0 == no weighting, 1 == normalisation fix, 2 == energy dependent scaling
    _util.GetPiZeroWeight(weight, _util.pi0_correction, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);

    // Set the CV weight variable in the slice container
    SC.SetCVWeight(weight);
    
    // *************************************************************************
    // Fill Histograms
    // *************************************************************************
    // Fill almost all the histograms with this function call
    if (!_util.slim) _hhelper.at(type).FillHists(type, SC.classification.second, SC.genie_interaction, SC.particle_type.second, cut_index, SC, weight);

    // Fill Plots for Efficiency
    if (!_util.slim && type == _util.k_mc) _hhelper.at(type).FillTEfficiency(cut_index, SC.classification.first, SC, weight);

    // Fill the dedx ttree before shr dist dedx cut and after cut dedx
    // We can use this tree to play around and optimise the dedx cut. Not essential for the analysis
    if (cut_index == _util.k_vtx_dist_dedx - 1 || cut_index == _util.k_vtx_dist_dedx){
        _thelper.at(type).Fill_dedxVars(SC, SC.classification, _util.cut_dirs.at(cut_index), weight);
    }
    
    // *************************************************************************
    // Tabulate the selection i.e count everything
    // *************************************************************************
    _util.Tabulate(is_in_fv, SC.genie_interaction, SC.classification.first, SC.pi0_classification, type, counter_v.at(cut_index), weight );

}
// -----------------------------------------------------------------------------
void Selection::ApplyPiZeroSelection(int type, SliceContainer &SC){

    bool pass; // A flag to see if an event passes an event

    // Classify the event
    SC.SliceClassifier(type);      // Classification of the event
    SC.SliceInteractionType(type); // Genie interaction type
    SC.ParticleClassifier(type);   // The truth matched particle type of the leading shower
        
    // *************************************************************************
    // Pi0 Selection Cuts ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.pi_zero_cuts(SC);
    if(!pass) return; // Failed the cut!

    bool is_in_fv = _util.in_fv(SC.true_nu_vtx_sce_x, SC.true_nu_vtx_sce_y, SC.true_nu_vtx_sce_z); // This variable is only used in the case of MC, so it should be fine 

    // Get the Central Value weight
    double weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv);

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
    double weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv);
    
    // Also apply the pi0 weight
    _util.GetPiZeroWeight(weight, _util.pi0_correction, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);

    if (!_util.slim) _hhelper.at(type).FillNuMuHists(SC.classification.second, SC, weight);

}
// -----------------------------------------------------------------------------
} // END NAMESPACE xsecSelection
