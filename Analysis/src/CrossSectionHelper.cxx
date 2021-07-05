#include "../include/CrossSectionHelper.h"

// -----------------------------------------------------------------------------
void CrossSectionHelper::Initialise(Utility _utility){

    std::cout << "Initalising Cross Section Helper..." << std::endl;
    _util = _utility;

    // Set the scale factors
    if (strcmp(_util.run_period, "1") == 0){
        mc_flux_scale_factor   = flux_scale_factor * _util.config_v.at(_util.k_Run1_MC_POT);
        data_flux_scale_factor = flux_scale_factor * _util.config_v.at(_util.k_Run1_Data_POT);
    }
    else if (strcmp(_util.run_period, "3") == 0){
        mc_flux_scale_factor   = flux_scale_factor * _util.config_v.at(_util.k_Run3_MC_POT); 
        data_flux_scale_factor = flux_scale_factor * _util.config_v.at(_util.k_Run3_Data_POT);
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }
    
    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Flux Scale factor:   " << mc_flux_scale_factor     << "\n" <<
    "Data Flux Scale factor: " << data_flux_scale_factor
    << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( _util.tree_file_name ) ) {
        f_nuexsec = new TFile( _util.tree_file_name, "READ");
    }

    // Get the Input TTree
    _util.GetTree(f_nuexsec, tree, "tree");
    if (tree == NULL) {
        std::cout << "Error failed to get the standard tree, maybe this is a MC only file, so trying to get the mc_tree"<< std::endl;
        _util.GetTree(f_nuexsec, tree, "mc_tree");
    }

    // Initialise the tree (set a bunch of tree branches)
    InitTree();

    // Initialise the Flux file
    if (std::string(_util.run_period) == "1"){
        
        // Switch the file path depending on whether we are on the gpvm or not
        if (!_util.use_gpvm)
            flux_file_name = "Systematics/output_fhc_uboone_run0.root";
        else
            flux_file_name = "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/sept/output_uboone_fhc_run0_merged.root";


        std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f_flux, flux_file_name);
        if (boolfile == false) gSystem->Exit(0); 
    }
    else if (std::string(_util.run_period) == "3") {
        if (!_util.use_gpvm)
            flux_file_name = "Systematics/output_rhc_uboone_run0.root";
        else
            flux_file_name = "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/sept/output_uboone_rhc_run0_merged.root";
        
        
        std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f_flux, flux_file_name );
        if (boolfile == false) gSystem->Exit(0); 
    }
    else{
        std::cout << "Unknown Run period configured, exiting..." << std::endl;
        return;
    }
    
    // Get the integrated flux for the CV
    integrated_flux = GetIntegratedFluxCV();

    // Turn this on for using the FLUGG flux as the integrated flux
    if (_util.usefluggflux){
        std::cout << "Overwriting default with flugg flux input"<< std::endl;
        integrated_flux = GetIntegratedFluxFLUGG();
    }


    // Now lets open the beamline variation files
    GetBeamlineHists();

    f_nuexsec->cd();

    // Calculate the volume as used in the cuts
    volume = (_util.config_v.at(_util.k_config_x2) - _util.config_v.at(_util.k_config_x1)) * 
             (_util.config_v.at(_util.k_config_y2) - _util.config_v.at(_util.k_config_y1)) * 
             (_util.config_v.at(_util.k_config_z2) - _util.config_v.at(_util.k_config_z1));

    std::cout << "Volume used in cuts: " << volume << std::endl;

    N_target_MC   = (lar_density_mc   * volume * NA * N_nuc) / m_mol; // Now use the same number of targets for MC and data since we fixed the simulated density in MCC9
    std::cout << "Number of Target Nucleons in MC: " << N_target_MC << std::endl;
    
    N_target_Data = (lar_density_data * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in Data: " << N_target_Data << std::endl;
    std::cout << "  "<< std::endl;

    // Create and initialise vector of histograms -- dont do this in the case we want to just write out the file lists
    InitialiseHistograms(std::string(_util.xsecmode));

    // Reweight the events by cut to get the sys uncertainty for certain plots in the selection
    // do this rather than reweight events at the end of the selection
    if (std::string(_util.xsec_rw_mode) == "rw_cuts"){
        LoopEventsbyCut();
        return;
    }

    // Initialise the file lists
    if (std::string(_util.xsecmode) == "txtlist"){
        // Initialise output file lists
        evt_dist_sig.open(Form("files/txt/run%s_evt_dist_sig.txt",_util.run_period));
        
        evt_dist_gen.open(Form("files/txt/run%s_evt_dist_gen.txt",_util.run_period));
        
        evt_dist_bkg.open(Form("files/txt/run%s_evt_dist_bkg.txt",_util.run_period));

        evt_dist_data.open(Form("files/txt/run%s_evt_dist_data.txt",_util.run_period));
        evt_dist_data << _util.vars.at(k_var_recoX) << "\n";

        f_out = TFile::Open("files/trees/krish_test.root", "UPDATE");

        // Also init the TTrees
        event_tree = new TTree("events", "events");
        event_tree->Branch("isData",     &isData);
        event_tree->Branch("isSignal",   &isSignal);
        event_tree->Branch("isSelected", &isSelected);
        event_tree->Branch("xTrue",      &xTrue);
        event_tree->Branch("xReco",      &xReco);
        event_tree->Branch("cv_weight",  &weight);
        event_tree->Branch("weights", "std::vector<float>", &ev_weight);
        meta_tree  = new TTree("metadata", "metadata");
        integrated_flux_tree = integrated_flux*flux_scale_factor;
        meta_tree->Branch("flux",       &integrated_flux_tree);
        meta_tree->Branch("nTargets",   &n_targ);
        meta_tree->Branch("potData",    &potData);
        meta_tree->Branch("potMC",      &potMC);
        meta_tree->Branch("nUniverses", &nUniverses);
        meta_tree->Fill();
    }

    // Load in the detector variation CV only if detector variation
    if (_util.isvariation && std::string(_util.variation) != "CV"){
        LoadDetvarCVHist();
    }
    
    // Load in NuWro file and add true cross section to file
    if (std::string(_util.xsec_rw_mode) == "gen"){
        SaveGenXSec();
        return;
    }

    // Load in NuWro file and add true cross section to file
    if (std::string(_util.xsec_rw_mode) == "pi0"){
        CheckPi0CrossSection();
        return;
    }

    if (std::string(_util.xsec_rw_mode) == "fluxshape"){
        StudyFluxShape();
        return;
    }

    // Now loop over events and caluclate the cross section
    LoopEvents();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoopEvents(){

    
    std::cout << "Total Tree Entries: "<< tree->GetEntries() << std::endl;

    // Loop over the tree entries and weight the events in each universe
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){
        
        tree->GetEntry(ievent); 

        if (ievent % 20000 == 0 && ievent > 0) std::cout << "On entry " << ievent/10000.0 <<"0k " << std::endl;

        // Set the reco and true variables to bin in
        
        // Electron/Shower Energy
        if (std::string(_util.xsec_var) =="elec_E"){
            recoX = shr_energy_cali;
            trueX = elec_e;
        }
        // Electron/Shower beta
        else if (std::string(_util.xsec_var) =="elec_ang"){
            recoX = effective_angle;
            trueX = true_effective_angle;
        }
        // Electron/Shower cos(beta)
        else if (std::string(_util.xsec_var) =="elec_cang"){
            recoX = cos_effective_angle;
            trueX = cos_true_effective_angle;
        }
        else {
            std::cout << "Unsupported parameter...exiting!" << std::endl;
            return;
        }

        if (shr_energy_cali > 8.0 && *classification == "data" ) std::cout << _util.red << "reco shower energy was:  " << shr_energy_cali << "  Consider updating the bins " << run << " " << subrun << " " << event  << _util.reset<<std::endl;

        double cv_weight = weight; // SplinetimesTune * PPFX CV * Pi0 Tune
        double weight_dirt = weight; // Use this for estimating dirt and POT sys
        double weight_ext  = weight;  // Use this for estimating ext

        // Loop over the reweighter labels
        for (unsigned int label = 0; label < reweighter_labels.size(); label++){
            
            // Call switch function
            SwitchReweighterLabel(reweighter_labels.at(label));

            // Ensure this isnt ever zero
            if (vec_universes.size() == 0) {
                for (unsigned int t=0; t < h_cross_sec.at(label).size(); t++){
                    vec_universes.push_back(1.0);
                }
            }
            
            if (ievent == 0) std::cout << "Running over reweighter label: " << reweighter_labels.at(label) << " with " << vec_universes.size() << " universes" << std::endl;

            // Now loop over the universes
            for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){
                
                // Update the CV weight to CV * universe i
                double weight_uni{1.0}; 

                // Set the weight for universe i
                SetUniverseWeight(reweighter_labels.at(label), weight_uni, weight_dirt, weight_ext, weightSplineTimesTune, *classification, cv_weight, uni, nu_pdg, true_energy, numi_ang, npi0, pi0_e, ccnc);

                // Signal event
                if ((*classification == "nue_cc" || *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "unmatched_nuebar") && passed_selection) {
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){
                           
                        FillHists(label, uni, k_xsec_sig, weight_uni, recoX, trueX);
                        FillHists(label, uni, k_xsec_sel, cv_weight, recoX, trueX);  // Selected events (N term) we dont weight

                        // Fill 2D histo of true energy and angle only for the CV
                        if (std::string(_util.xsec_bin_mode) == "e_ang" && reweighter_labels.at(label) == "CV" && uni == 0){
                            h_2D_CV_binning->Fill(elec_e, true_effective_angle, cv_weight);
                        }
                    
                    }
                    // Just save the event weights
                    else {
                        ev_weight.push_back(weight_uni);
                    }

                }

                // Background event
                if ( *classification == "nu_out_fv"  || *classification == "cosmic"       ||
                     *classification == "numu_cc"    || *classification == "numu_cc_pi0"  || *classification == "nc" || 
                     *classification == "nc_pi0"     || ((*classification == "cosmic_nue" || *classification == "cosmic_nuebar" || 
                     *classification == "thr_nue"    || *classification == "thr_nuebar") && passed_selection)  ) {
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){


                        FillHists(label, uni, k_xsec_bkg, weight_uni, recoX, trueX);
                        FillHists(label, uni, k_xsec_sel, cv_weight, recoX, trueX);  // Selected events (N term) we dont weight
                    
                    }
                    // Just save the event weights
                    else {
                        ev_weight.push_back(weight_uni);
                    }
                    
                }

                // Generated event
                if (  *classification == "nue_cc"           || *classification == "nuebar_cc" || 
                      *classification == "unmatched_nue"    || *classification == "cosmic_nue" ||
                      *classification == "unmatched_nuebar" || *classification == "cosmic_nuebar") {
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){
                        
                        FillHists(label, uni, k_xsec_gen, weight_uni, recoX, trueX);
                        FillHists(label, uni, k_xsec_gen_smear, weight_uni, recoX, trueX);
                        FillHists(label, uni, k_xsec_gen_shape, weight_uni, recoX, trueX);
                    
                    }
                    // Just save the event weights
                    else {
                        if (!passed_selection)
                            ev_weight.push_back(weight_uni);
                    }
                }
                
                // Data event
                if (*classification == "data"){

                    if (cv_weight != 1.0 && !_util.isfakedata) std::cout << "Error weight for data is not 1, this means your weighting the data... bad!"<< std::endl;

                    // If fake data, dont fill if we have a generated event!
                    if (_util.isfakedata && !passed_selection)
                        continue;
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){
                    
                        FillHists(label, uni, k_xsec_data, cv_weight, recoX, trueX);

                    }
                    // Just save the event weights
                    else {
                        ev_weight.push_back(1.0);
                    }
                }

                // Off beam event
                if (*classification == "ext"){
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){
                        
                        FillHists(label, uni, k_xsec_ext, weight_ext, recoX, trueX);

                        // Apply additional weight to the ext events to get the N selected number correct
                        FillHists(label, uni, k_xsec_sel, cv_weight*(_util.ext_scale_factor / _util.mc_scale_factor), recoX, trueX);

                    }
                    // Just save the event weights
                    else {
                        ev_weight.push_back(cv_weight*(_util.ext_scale_factor / _util.mc_scale_factor));
                    }
                }

                // Dirt event
                if (*classification == "dirt"){
                    
                    // Fill histograms
                    if (std::string(_util.xsecmode) != "txtlist"){
                        
                        FillHists(label, uni, k_xsec_dirt, weight_dirt, recoX, trueX);

                        // Apply additional weight to the dirt events to get the N selected number correct
                        FillHists(label, uni, k_xsec_sel, cv_weight*(_util.dirt_scale_factor / _util.mc_scale_factor), recoX, trueX);
                    }
                    // Just save the event weights
                    else {
                        ev_weight.push_back(cv_weight*(_util.dirt_scale_factor / _util.mc_scale_factor));
                    }
                }
            } // End loop over uni

        } // End loop over labels

        // Call the function to save the events and weights to a file
        if (std::string(_util.xsecmode) == "txtlist"){
            SaveEvent(*classification, passed_selection, ev_weight, recoX, trueX);
        }

        // Clear vector to store the weights
        ev_weight.clear();
        
    } // End loop over events

    // If we want to just write out the event numbers to a txt file
    if (std::string(_util.xsecmode) == "txtlist"){
        std::cout << "Writing out event list to text file..."<< std::endl;
        f_out->cd();
        event_tree->Write("", TObject::kOverwrite);
        meta_tree ->Write("", TObject::kOverwrite);

        return;
    }


    std::cout << "\n\nFinished Event Loop, now calculating the cross-sections\n\n"<< std::endl;

    // ----

    // Now we have rewieghted the evnts, we want to calculate the cross-sections for each label for each universe

    // This part seems quite slow now we added many variables to plot the cross-section against. Maybe there is a way to speed this up?
    
    // Loop over the reweighter labels
    for (unsigned int label = 0; label < reweighter_labels.size(); label++){

        std::cout << "calculating cross sections for: " << reweighter_labels.at(label) << std::endl;
        
        // Now loop over the universes
        for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){

            if (uni % 100 == 0 && uni > 0) std::cout << "On universe " << uni << std::endl;

            // Loop over the variables
            for (unsigned int var = 0; var < h_cross_sec.at(label).at(uni).size(); var ++){

                // Set the integrated flux
                double temp_integrated_flux = integrated_flux;

                // If we are reweighting by the PPFX Multisims, we need to change the integrated flux too
                if (reweighter_labels.at(label) == "weightsPPFX") 
                    temp_integrated_flux = GetIntegratedFluxHP(uni, "ppfx_ms_UBPPFX");

                // If this is a beamline variation then we use the corresponding beamline flux
                if (CheckBeamline(reweighter_labels.at(label))) {
                    int var_index = GetBeamlineIndex(reweighter_labels.at(label));
                    temp_integrated_flux = beamline_flux.at(var_index);
                }

                // We choose the POT counting uncertainty to be exactly 2% rather than just weighting the events
                // To do this, multiply the flux with a 2% scale factor. Since the flux is divided in the x-sec calc,
                // the scale factors are reverse
                double weight_POT = 1.0;
                if (reweighter_labels.at(label) == "POTup") weight_POT =1.0/1.02;
                if (reweighter_labels.at(label) == "POTdn") weight_POT =1.0/0.98;

                // Calculate the efficiency histogram with smearing of the truth
                if (var == k_var_recoX){
                    Smear(h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_sig), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen),
                          h_smear.at(label).at(uni).at(k_var_recoX), h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_eff));
                }
                // Calculate the efficiency histogram by dividing the sig and gen
                else { 

                    // In the case of PPFX and Beamline, use the wiener method of taking the ratio of event distributions
                    if ( (reweighter_labels.at(label) == "weightsPPFX" || CheckBeamline(reweighter_labels.at(label))) && var == k_var_integrated){
                        h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff)->Divide(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig), h_cross_sec.front().front().at(var).at(k_xsec_gen));
                    }
                    else {
                        h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff)->Divide(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig), h_cross_sec.at(label).at(uni).at(var).at(k_xsec_gen));
                    }
                }

                // For the x-sec in true space, appy a response matrix to the generated events
                if (var == k_var_trueX){
                    
                    // Use standard binning to smear
                    if (std::string(_util.xsec_bin_mode) == "standard" || std::string(_util.xsec_bin_mode) == "ratio"){

                        // For shape uncertainties
                        // TH1D* h_shape_temp = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Clone();
                        // ApplyResponseMatrix(h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape),
                        // h_shape_temp, h_smear.front().front().at(k_var_trueX), false);
                        
                        // For making covariance matrix on response * gen + B
                        if (reweighter_labels.at(label) == "weightsPPFX" || CheckBeamline(reweighter_labels.at(label))){
                            TH2D* h_smear_temp = (TH2D*) h_smear.at(label).at(uni).at(k_var_trueX)->Clone();

                            ApplyResponseMatrix( h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape),
                            h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_smear_temp, true);
                            delete h_smear_temp;

                            // Standard Smearing
                            ApplyResponseMatrix(h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_smear),
                            h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_smear.at(label).at(uni).at(k_var_trueX), true);

                        }
                        else {
                            TH2D* h_smear_temp = (TH2D*) h_smear.at(label).at(uni).at(k_var_trueX)->Clone();
                            ApplyResponseMatrix(h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape),
                            h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_smear_temp, true);
                            delete h_smear_temp;

                            // Standard Smearing
                            ApplyResponseMatrix(h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_smear),
                            h_cross_sec.front().front().at(k_var_trueX).at(k_xsec_gen), h_smear.at(label).at(uni).at(k_var_trueX), true);
                        }

                        TH1D* h_ext     = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_ext)->Clone();
                        TH1D* h_dirt    = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_dirt)->Clone();
                        TH1D* h_ext_CV  = (TH1D*)h_cross_sec.front().front().at(k_var_recoX).at(k_xsec_ext)->Clone();
                        TH1D* h_dirt_CV = (TH1D*)h_cross_sec.front().front().at(k_var_recoX).at(k_xsec_dirt)->Clone();

                        // Scale histograms to add in with right scaling
                        h_ext    ->Scale(_util.ext_scale_factor / _util.mc_scale_factor);
                        h_dirt   ->Scale(_util.dirt_scale_factor / _util.mc_scale_factor);
                        h_ext_CV ->Scale(_util.ext_scale_factor / _util.mc_scale_factor);
                        h_dirt_CV->Scale(_util.dirt_scale_factor / _util.mc_scale_factor);

                        // Add in the reco background
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_bkg));
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_ext);
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_dirt);

                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_cross_sec.front().front().at(k_var_recoX).at(k_xsec_bkg), -1);
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_ext_CV, -1);
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape)->Add(h_dirt_CV, -1);

                        // Store the bkg so we can draw it
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_bkg)->Add(h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_bkg));
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_bkg)->Add(h_ext);
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_bkg)->Add(h_dirt);

                        delete h_ext;
                        delete h_dirt;
                        delete h_ext_CV;
                        delete h_dirt_CV;
                        
                    }
                }

                // MC Cross section -- currently using eventrate
                CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sel), // N Sel
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),  // Eff
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),  // N Bkg
                                _util.mc_scale_factor,
                                weight_POT * integrated_flux * mc_flux_scale_factor, // Flux
                                _util.ext_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),  // N EXT
                                _util.dirt_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dirt), // N Dirt
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec),
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig),  // N Sig
                                 N_target_MC, "MC", var);

                if (var != k_var_trueX){
                    // Data Cross section -- currently using eventrate
                    CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_data), // N Sel
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),   // Eff
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),   // N Bkg
                                    _util.mc_scale_factor,
                                    weight_POT * integrated_flux * data_flux_scale_factor, // Flux
                                    _util.ext_scale_factor,
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),   // N EXT
                                    _util.dirt_scale_factor,
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dirt),  // N Dirt
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dataxsec),
                                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig),   // N Sig
                                    N_target_Data, "Data", var);
                }

                if (var == k_var_trueX){
                    // Calculate the the smeared cross section prediction
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_smear)->Add(h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_smear), 1);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_smear)->Scale(1.0 / (integrated_flux * mc_flux_scale_factor * N_target_MC));
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_smear)->Scale(1.0e39);

                    // Calculate the the cross section shape prediction
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_shape)->Add(h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_gen_shape), 1);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_shape)->Scale(1.0 / (weight_POT * integrated_flux * mc_flux_scale_factor * N_target_MC));
                    // h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_shape)->Scale(integrated_flux / temp_integrated_flux);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_shape)->Scale(1.0e39);

                    // Normalise the data event rate for calculations and an easier way to propagate the satistical uncertainties
                    TH1D* h_bkg     = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_bkg)->Clone();
                    TH1D* h_ext     = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_ext)->Clone();
                    TH1D* h_dirt    = (TH1D*)h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_dirt)->Clone();
                    h_bkg    ->Scale(_util.mc_scale_factor);
                    h_ext    ->Scale(_util.ext_scale_factor);
                    h_dirt   ->Scale(_util.dirt_scale_factor);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec)->Add(h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_data), 1);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec)->Add(h_bkg, -1);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec)->Add(h_ext, -1);
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec)->Add(h_dirt, -1);
                    delete h_bkg; delete h_ext; delete h_dirt;
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec)->Scale(1.0 / (weight_POT * integrated_flux * data_flux_scale_factor * N_target_Data));

                    // If the crosssec is a function of angle then scale again by a factor of 100
                    if (std::string(_util.xsec_var) == "elec_ang"){
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_smear)->Scale(100);
                        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_mcxsec_shape)->Scale(100);
                    }

                    // Scale the true bkg down
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_bkg)->Scale(1.0 / (weight_POT * integrated_flux * mc_flux_scale_factor * N_target_MC));
                    h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_bkg)->Scale(1.0e39);

                }

                // Scale the histograms to avoid working with really small numbers
                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec)->Scale(1.0e39);
                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dataxsec)->Scale(1.0e39);

                // If the crosssec is a function of angle then scale again by a factor of 100
                if (std::string(_util.xsec_var) == "elec_ang" && var != k_var_integrated){
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec)->Scale(100);
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dataxsec)->Scale(100);
                }

                // If using the ratio option then calculate the cross section ratio
                if (std::string(_util.xsec_bin_mode) == "ratio" && var != k_var_integrated){
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dataxsec)->Scale(1.0/h_cross_sec.at(label).at(uni).at(k_var_integrated).at(k_xsec_dataxsec)->GetBinContent(1));
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec)->Scale(1.0/h_cross_sec.at(label).at(uni).at(k_var_integrated).at(k_xsec_mcxsec)->GetBinContent(1));
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec_smear)->Scale((1.0/h_cross_sec.at(label).at(uni).at(k_var_integrated).at(k_xsec_mcxsec)->GetBinContent(1)));
                    h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec_shape)->Scale(1.0/h_cross_sec.at(label).at(uni).at(k_var_integrated).at(k_xsec_mcxsec)->GetBinContent(1));
                }

                // if (var == k_var_trueX){
                //     UnregularizedUnfold(h_cross_sec.at(label).at(uni).at(k_var_recoX).at(k_xsec_dataxsec), h_cross_sec.at(label).at(uni).at(k_var_trueX).at(k_xsec_dataxsec),  h_smear.at(label).at(uni).at(k_var_trueX));
                // }


                // if (var == 0) std::cout << reweighter_labels.at(label) << ": " << _util.red << h_cross_sec.at(label).at(uni).at(k_var_integrated).at(k_xsec_mcxsec)  ->Integral() << _util.reset<< " x10^-39 cm2/nucleon" << std::endl;

            } // End loop over the _util.vars

            
        
        } // End loop over universes
    
    } // End loop over labels

    // Print the CV Results for the Flux Integrated Measurement
    // Label 0 should always be the CV with 1 universe (which counts from 0)

    std::cout << "\n" <<
    "Selected MC:     " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sel) ->Integral() << "\n" << 
    "Signal MC:       " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sig) ->Integral() << "\n" << 
    "Background MC:   " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_bkg) ->Integral() << "\n" << 
    "Generated MC:    " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_gen) ->Integral() << "\n" << 
    "Efficiency:      " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_eff) ->Integral() << "\n" << 
    "EXT MC:          " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_ext) ->Integral()* (_util.ext_scale_factor / _util.mc_scale_factor) << "\n" << 
    "Dirt MC:         " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dirt)->Integral()* (_util.dirt_scale_factor / _util.mc_scale_factor) << "\n\n" << 
    
    "Selected Data:   " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_data)->Integral() << "\n" << 
    "Signal Data:     " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sig) ->Integral()* _util.mc_scale_factor << "\n" << 
    "Background Data: " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_bkg) ->Integral()* _util.mc_scale_factor << "\n" << 
    "Generated Data:  " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_gen) ->Integral()* _util.mc_scale_factor << "\n" << 
    "EXT Data:        " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_ext) ->Integral()* _util.ext_scale_factor  << "\n" << 
    "Dirt Data:       " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dirt)->Integral()* _util.dirt_scale_factor << "\n"
    << std::endl;

    std::cout <<
    "GENIE X-Sec:           " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sig) ->Integral() / (h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_eff)->Integral() *integrated_flux * mc_flux_scale_factor * N_target_MC) * 1.0e39 << " x10^-39 cm2/nucleon" << "\n" << 
    "MC CC Cross Section:   " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_mcxsec)  ->Integral() << " x10^-39 cm2/nucleon" << "\n" << 
    "Data CC Cross Section: " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dataxsec)->Integral() << " x10^-39 cm2/nucleon      \n"
    << std::endl;

    // Write the histograms to file for inspection
    WriteHists();


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoopEventsbyCut(){

    // Create the class instances we need
    SelectionCuts  _scuts;
    _scuts.Initalise(_util);
    
    SliceContainer SC;

    int treeNumber = -1; // TTree number in TChain

    // Create the TChain and add the MC and intrinsic events
    TChain *in_chain = new TChain("nuselection/NeutrinoSelectionFilter");

    // For now we need to keep the intrinsic file name in this format unless we add more options for configuration
    if (std::string(_util.intrinsic_mode) == "intrinsic")
        in_chain->Add(Form("../ntuples/neutrinoselection_filt_run%s_overlay_intrinsic_newtune.root", _util.run_period)); 
    
    // Add the Main MC file
    in_chain->Add(_util.mc_file_name);

    // Initialise the TChain and Slice Container class
    SC.Initialise(in_chain, _util.k_mc, _util);

    // Event loop
    for (int ievent = 0; ievent < in_chain->GetEntries(); ++ievent){

        // See if we want to process all the events
        if (_util.num_events > 0){
            if (ievent >= _util.num_events) break;
        }

        // Alert the user
        if (ievent % 100000 == 0) std::cout << "On entry " << ievent/100000.0 <<"00k " << std::endl;

        // Get the entry in the tree
        in_chain->GetEntry(ievent);

        // Print the TChain number
        if(treeNumber != in_chain->GetTreeNumber()) {
            treeNumber = in_chain->GetTreeNumber();
            std::cout << "Moving to tree number " << treeNumber << "." << std::endl;
            
            // if there is more than one tree and its the second file, make sure the intrinsic mode is turned off
            // otherwise we will apply the intrinsic weights to the bkg, this is not what we want
            if (treeNumber == 1){
                _util.TurnoffIntrinsicMode();
            }

            std::cout << "Intrinsic Mode is: " << std::string(_util.intrinsic_mode) << std::endl;

        }

        // Apply the selection cuts and fill histograms for each universe
        ApplyCuts(_util.k_mc, SC, _scuts, treeNumber);

    }

    // Create the direcotry structure in the output file if it doesnt already exist and then write the histograms to file
    WriteHists();

}
// -----------------------------------------------------------------------------
bool CrossSectionHelper::ApplyCuts(int type, SliceContainer &SC, SelectionCuts _scuts, int treeNum){

    bool pass = true;

    SC.ReClassifyPileUps(type);

    // Set derived variables in the slice container
    // Classify the event
    SC.SliceClassifier(type);      // Classification of the event

    // If we have a signal event that is below threshold, then set its category to thr_nue or thr_nuebar
    SC.SetThresholdEvent(type);

    SC.SetSignal();                // Set the event as either signal or other
    SC.SetTrueElectronThetaPhi();  // Set the true electron theta and phi variables
    SC.SetNuMIAngularVariables();  // Set the NuMI angular variables
    SC.CalibrateShowerEnergy();    // Divide the shower energy by 0.83 so it is done in one place

    // Skip signal events in the standard MC file so we dont double count the signal events
    if (treeNum == 1 && SC.is_signal){
        return false;
    }

    // *************************************************************************
    // Unselected---------------------------------------------------------------
    // *************************************************************************
    FillCutHists(type, SC, SC.classification, _util.k_unselected );
    
    // *************************************************************************
    // Software Trigger -- MC Only  --------------------------------------------
    // *************************************************************************
    pass = _scuts.swtrig(SC, type);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_swtrig );

    // *************************************************************************
    // Slice ID ----------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.slice_id(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_slice_id );
    
    // *************************************************************************
    // Electron Candidate ------------------------------------------------------
    // *************************************************************************
    pass = _scuts.e_candidate(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_e_candidate );

    // *************************************************************************
    // In FV -------------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.in_fv(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_in_fv );
    
    // *************************************************************************
    // Slice Contained Fraction ------------------------------------------------
    // *************************************************************************
    pass = _scuts.contained_frac(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_contained_frac );

    // *************************************************************************
    // Topological Score -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.topo_score(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_topo_score );

    // *************************************************************************
    // Cosmic Impact Parameter -------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_cosmic_IP(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_cosmic_ip );

    // *************************************************************************
    // Shower Score ------------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shower_score(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_shower_score );

    // *************************************************************************
    // Shower Hit Ratio  -------------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_hitratio(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_hit_ratio );

    // *************************************************************************
    // Shower Moliere Average --------------------------------------------------
    // *************************************************************************
    pass = _scuts.shr_moliere_avg(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_shr_moliere_avg );

    // *************************************************************************
    // 2D cut for Shower to Vertex Distance and dEdx ---------------------------
    // *************************************************************************
    pass = _scuts.shr_dist_dEdx_max(SC);
    if(!pass) return false; // Failed the cut!
    
    FillCutHists(type, SC, SC.classification, _util.k_vtx_dist_dedx );

    // *************************************************************************
    // dEdx in all planes for 0 track events -----------------------------------
    // *************************************************************************
    pass = _scuts.dEdx_max_no_tracks(SC);
    if(!pass) return false; // Failed the cut!

    // If the backtracked pdg of the leading shower is not an electron then alter classification
    // Turn these off to get the efficiencies at low energies correct
    SC.SetNonLdgShrEvent(type);
    
    FillCutHists(type, SC, SC.classification, _util.k_dEdx_max_no_tracks );

    // **************************************************************************
    return true;

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::FillCutHists(int type, SliceContainer &SC, std::pair<std::string, int> classification, int cut_index){

    // Loop over the reweighter labels
    for (unsigned int label = 0; label < reweighter_labels.size(); label++){

        // Call switch function
        SwitchReweighterLabel(reweighter_labels[label], SC);

        bool is_in_fv = _util.in_fv(SC.true_nu_vtx_x, SC.true_nu_vtx_y, SC.true_nu_vtx_z); // This variable is only used in the case of MC, so it should be fine 

        // Get the CV weight
        double cv_weight = _util.GetCVWeight(type, SC.weightSplineTimesTune, SC.ppfx_cv, SC.nu_e, SC.nu_pdg, is_in_fv, SC.interaction, SC.elec_e);
        _util.GetPiZeroWeight(cv_weight, _util.pi0_correction, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e);
        
        double weight_dirt = cv_weight;
        double weight_ext = cv_weight;

        // This bit of code wont work if the vector size is zero, put this here to catch the error
        if (vec_universes.size() == 0) {
            std::cout << "Vector size is zero" <<  std::endl;
            exit(3);
        }

        // Now loop over the universes
        for (unsigned int uni = 0; uni < vec_universes.size(); uni++){

            // Update the CV weight to CV * universe i
            double weight_uni{cv_weight}; 

            SetUniverseWeight(reweighter_labels[label], weight_uni, weight_dirt, weight_ext, SC.weightSplineTimesTune, classification.first, cv_weight, uni, SC.nu_pdg, SC.nu_e, SC.nu_angle, SC.npi0, SC.pi0_e, SC.ccnc);

            double dedx_max = SC.GetdEdxMax();

            double pi0_tuned_weight = weight_uni;
            _util.GetPiZeroWeight(pi0_tuned_weight, 1, SC.nu_pdg, SC.ccnc, SC.npi0, SC.pi0_e); // 0 == no weighting, 1 == normalisation fix, 2 == energy dependent scaling


            // Now we got the weight for universe i, lets fill the histograms :D
            // Use [] rather than .at() to speed this process up. This can cause errors if indexes go out of bound
            h_cut_v[label][cut_index][_util.k_cut_nslice][uni]                   ->Fill(SC.nslice,                 weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shower_multiplicity][uni]      ->Fill(SC.n_showers,              weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_track_multiplicity][uni]       ->Fill(SC.n_tracks,               weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_topological_score][uni]        ->Fill(SC.topological_score,      weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_vtx_x_sce][uni]                ->Fill(SC.reco_nu_vtx_sce_x,      weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_vtx_y_sce][uni]                ->Fill(SC.reco_nu_vtx_sce_y,      weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_vtx_z_sce][uni]                ->Fill(SC.reco_nu_vtx_sce_z,      weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shower_score][uni]             ->Fill(SC.shr_score,              weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shr_tkfit_dedx_max][uni]       ->Fill(dedx_max, weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shr_tkfit_dedx_max_tune][uni]  ->Fill(dedx_max, pi0_tuned_weight);
            
            if (SC.n_tracks > 0)
                h_cut_v[label][cut_index][_util.k_cut_shr_tkfit_dedx_max_with_tracks][uni]    ->Fill(dedx_max, weight_uni);
            
            if (SC.n_tracks == 0)
                h_cut_v[label][cut_index][_util.k_cut_shr_tkfit_dedx_max_no_tracks][uni]      ->Fill(dedx_max, weight_uni);
            
            h_cut_v[label][cut_index][_util.k_cut_shower_to_vtx_dist][uni]       ->Fill(SC.shr_distance,           weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_hits_ratio][uni]               ->Fill(SC.hits_ratio,             weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_CosmicIPAll3D][uni]            ->Fill(SC.CosmicIPAll3D,          weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_contained_fraction][uni]       ->Fill(SC.contained_fraction,     weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shrmoliereavg][uni]            ->Fill(SC.shrmoliereavg,          weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_leading_shower_theta][uni]     ->Fill(SC.shr_theta * 180/3.14159,weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_leading_shower_phi][uni]       ->Fill(SC.shr_phi * 180/3.14159,  weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shower_energy_cali][uni]       ->Fill(SC.shr_energy_cali,        weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_shower_energy_cali_rebin][uni] ->Fill(SC.shr_energy_cali,        weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_flash_time][uni]               ->Fill(SC.flash_time,             weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_flash_pe][uni]                 ->Fill(SC.flash_pe,               weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_effective_angle][uni]          ->Fill(SC.effective_angle,        weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_effective_cosangle][uni]       ->Fill(SC.cos_effective_angle,    weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_effective_angle_rebin][uni]          ->Fill(SC.effective_angle,        weight_uni);
            h_cut_v[label][cut_index][_util.k_cut_effective_cosangle_rebin][uni]       ->Fill(SC.cos_effective_angle,    weight_uni);
        
        }

    }



}
// -----------------------------------------------------------------------------
void CrossSectionHelper::SetUniverseWeight(std::string label, double &weight_uni, double &weight_dirt, double &weight_ext, double _weightSplineTimesTune,
                                           std::string _classification, double cv_weight, int uni, int _nu_pdg, double _true_energy, double _numi_ang, int _npi0, double _pi0_e, int _ccnc ){

    // Weight equal to universe weight times cv weight
    if (label == "weightsReint" || label == "CV" || label == "xsr_scc_Fv3up" || label == "xsr_scc_Fa3up" || label == "xsr_scc_Fv3dn" || label == "xsr_scc_Fa3dn"){
        
        _util.CheckWeight(vec_universes[uni]);

        weight_uni = cv_weight * vec_universes[uni];
    }
    // Hadron Production weights
    else if (label == "weightsPPFX"){

        // Get weight from ratio of flux histograms
        if (_util.usefluggflux){
            weight_uni = cv_weight * GetHPWeight(uni, "ppfx_ms_UBPPFX", _nu_pdg,  _true_energy, _numi_ang);

            if (std::isnan(weight_uni) == 1 || std::isinf(weight_uni) || weight_uni < 0 || weight_uni > 30 || weight_uni == 1.0) {
                 weight_uni = cv_weight;
            }

        }
        // Get weight from the art-root event
        else {

            _util.CheckWeight(vec_universes[uni]);

            weight_uni = cv_weight * vec_universes[uni];
        }

    }
    // This is a mc stats variation which studies the statisitcal uncertainty on the smearing matrix and efficiency by 
    // varying the signal and generated events
    // dont touch the cosmic contaiminated nues, dont know the expected behaviour for this
    else if (label == "MCStats"){
        
        // Weight the Signal events that factor into the smearing matrix
        if (_classification == "nue_cc" || _classification == "nuebar_cc" || _classification == "unmatched_nue" || _classification == "unmatched_nuebar" || _classification == "cosmic_nue" || _classification == "cosmic_nuebar"){
            
            weight_uni = cv_weight * PoissonRandomNumber(uni);
        
        }
        else
            weight_uni = cv_weight;
        
    }
    // This is a beamline variation
    else if (CheckBeamline(label)){
        weight_uni = cv_weight * GetBeamlineWeight(label, _nu_pdg, _true_energy, _numi_ang);
        // std::cout << GetBeamlineWeight(label, _nu_pdg, _true_energy, _numi_ang) << std::endl;
    }
    // Dirt reweighting
    else if ( label == "Dirtup" || label == "Dirtdn"){
        
        if (label == "Dirtup")
            weight_dirt = cv_weight*2.0; // increase the dirt by 100%
        else
            weight_dirt = cv_weight*0.0; // decrease the dirt by 0%
        
        weight_uni = cv_weight;
    }
    // EXT reweighting
    else if ( label == "EXT"){
        
        weight_ext = cv_weight*1.01; // increase the ext norm by 1%
        weight_uni = cv_weight;
    }
    // POT Counting
    else if ( label == "POTup" || label == "POTdn"){
        weight_uni  = cv_weight;
        weight_dirt = cv_weight;
    }
    // Pi0 Tune -- This undoes the pi0 tune to see what effect it has
    else if ( label == "pi0"){

        // Fix the normalisation
        if (_util.pi0_correction == 1){
            
            // Add an extra 10% to the CC pi0 events to account for stat err
            if (_npi0 > 0) {
                weight_uni = cv_weight * 0.91; 
            }
            else {
                weight_uni = cv_weight;
            }

        }
        // Try energy dependent scaling for pi0
        else if (_util.pi0_correction == 2){
            
            if (_npi0 > 0) {
                double pi0emax = 0.6;
                if (_pi0_e > 0.1 && _pi0_e < pi0emax){
                    weight_uni = cv_weight / (1 - 0.4 * pi0_e);
                }
                else if (_pi0_e > 0.1 && _pi0_e >= pi0emax){
                    weight_uni = cv_weight / (1 - 0.4 * pi0emax);
                }    
                else {
                    weight_uni = cv_weight;
                }            
            }
        }
        else {
            // Dont touch the weight
            weight_uni = cv_weight;
        }
    }
    // If we are using the genie systematics and unisim systematics then we want to undo the genie tune on them so we dont double count
    else {
        // Note we actually dont want to divide out by the spline, but since this is 1 in numi, it doesnt matter!
        // We do this because the interaction systematics are shifted about the genie tune as the CV

        // if (_npi0 > 0) {
        //     weight_uni = cv_weight;
        //     return;
        // }

        // Check the spline times tune weight
        _util.CheckWeight(_weightSplineTimesTune);

        // Check the uiverse weight -- Also need to set weights equal to 1.0 to the tune weight. Othwerwise we will introduce a bias when we divide this universe weight by the tune weight below.
        if (std::isnan(vec_universes[uni]) == 1 || std::isinf(vec_universes[uni]) || vec_universes[uni] < 0 || vec_universes[uni] > 30 || vec_universes[uni] == 1.0) {
            
            // We set the universe to be the spline times tune, so it cancels with the divide below to just return the CV weight
            // i.e. a universe weight of 1
            vec_universes.at(uni)   = _weightSplineTimesTune;
        }

        if (_weightSplineTimesTune == 0) weight_uni = cv_weight; // Special case where the genie tune or we have thresholded events are zero
        else weight_uni = (cv_weight * vec_universes[uni]) / _weightSplineTimesTune;
    }

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::CalcCrossSec(double sel, double gen, double sig, double bkg, double flux, double ext, double dirt, double targ){

    bool DEBUG{false};

    if (DEBUG) {
        std::cout << 
        "DEBUG:       " << "\n" << 
        "Sel:         " << sel  << "\n" << 
        "Bkg:         " << bkg  << "\n" << 
        "Flux:        " << flux << "\n" << 
        "Targets:     " << targ << "\n" << 
        "EXT:         " << ext  << "\n" << 
        "Dirt:        " << dirt << "\n" << 
        "efficiency:\t" << (sig / gen ) << std::endl;
    }

    // (S - B) / (eff * targ * flux)
    return (sel - ( bkg + ext + dirt) ) / ( (sig / gen ) * targ * flux );

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double ext_scale_factor, TH1D* h_ext, double dirt_scale_factor ,TH1D* h_dirt, TH1D* h_xsec, TH1D* h_sig, double targ, std::string mcdata, int _var){


    // I think this is the slow bit -- maybe make copies only once?
    TH1D *h_bkg_clone  = (TH1D*)h_bkg ->Clone("h_bkg_temp");
    TH1D *h_ext_clone  = (TH1D*)h_ext ->Clone("h_ext_temp");
    TH1D *h_dirt_clone = (TH1D*)h_dirt->Clone("h_dirt_temp");

    // Scale the relavent histograms to the MC/Data POT/Triggers
    if (mcdata == "MC"){
        h_ext_clone ->Scale(ext_scale_factor / mc_scale_factor);
        h_dirt_clone->Scale(dirt_scale_factor / mc_scale_factor);
    }
    else if (mcdata == "Data"){
        h_bkg_clone ->Scale(mc_scale_factor);
        h_ext_clone ->Scale(ext_scale_factor);
        h_dirt_clone->Scale(dirt_scale_factor);
    }
    else{
        std::cout << "error wrong mode entering xsec calculation" << std::endl;
    }

    // (N - B) / (eff * targ * flux)
    h_xsec->Sumw2();
    
    if (_var == k_var_trueX){
        h_xsec->Add(h_sig,         1);
    }
    else {
        h_xsec->Add(h_sel,         1);
    }
    
    // Subtract the backgrounds
    if (_var != k_var_trueX){
        h_xsec->Add(h_bkg_clone,  -1);
        h_xsec->Add(h_ext_clone,  -1);
        h_xsec->Add(h_dirt_clone, -1);
    }

    // If using Marco's Method we correct by the efficiency
    if (std::string(_util.xsec_smear_mode) == "mcc8" || _var == k_var_integrated || _var == k_var_trueX){
        h_xsec->Divide(h_eff) ;
    }

    h_xsec->Scale(1.0 / (targ*flux) );

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFluxCV(){

    f_flux->cd();

    TH2D* h_nue, *h_nuebar;

    double xbin;
    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    double xbin_th{0.0};
    double integral_nue{0.0}, integral_nuebar{0.0};
    bool boolhist;

    // Get the nue flux histogram from the file
    boolhist = _util.GetHist(f_flux, h_nue,         "nue/Detsmear/nue_CV_AV_TPC_2D");       
    if (boolhist == false) gSystem->Exit(0); 
    
    // Get the nuebar flux histogram from the file
    boolhist      = _util.GetHist(f_flux, h_nuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_2D");
    if (boolhist == false) gSystem->Exit(0); 

    // Get the POT from the flux histogram
    POT_flux = GetPOT(f_flux);

    // Nue Integrated Flux
    xbin_th = h_nue->GetXaxis()->FindBin(_util.energy_threshold);                     // Find the x bin to integrate from
    integral_nue = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1, 0, h_nue->GetNbinsY()+1); // Integrate over the flux for nue
    std::cout << "\nIntegral Nue Flux: " << flux_scale_factor * integral_nue / POT_flux << " nue / POT / cm2" << std::endl;

    // Nuebar Integrated Flux
    xbin_th   = h_nuebar->GetXaxis()->FindBin(_util.energy_threshold);                                   // Find the x bin to integrate from
    integral_nuebar = h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1, 0, h_nuebar->GetNbinsY()+1); // Integrate over the flux for nue
    std::cout << "\nIntegral Nuebar Flux: " << flux_scale_factor * integral_nuebar / POT_flux << " nuebar / POT / cm2" << "\n" << std::endl;

    // Print the flux scaled to the MC and Data POT
    std::cout << "Integral Flux MC POT: "   << mc_flux_scale_factor * (integral_nue + integral_nuebar) / POT_flux   << " nu / MC POT / cm2"   << "\n" << std::endl;
    std::cout << "Integral Flux Data POT: " << data_flux_scale_factor * (integral_nue + integral_nuebar) / POT_flux << " nu / Data POT / cm2" << "\n" << std::endl;

    // Return the flux per POT
    return (integral_nue + integral_nuebar) / POT_flux;


}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFluxFLUGG(){

    std::cout << "\n\n"<< "Printing Flux values from FLUGG flux file" << std::endl; 

    TFile *f_flugg = TFile::Open("Systematics/MCC8_FHC_flux.root");

    TH1D* h_nue, *h_nuebar;

    double xbin;
    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    double xbin_th{0.0};
    double integral_nue{0.0}, integral_nuebar{0.0};
    bool boolhist;

    // Get the nue flux histogram from the file
    boolhist = _util.GetHist(f_flugg, h_nue,         "nueFluxHisto");       
    if (boolhist == false) gSystem->Exit(0); 
    
    // Get the nuebar flux histogram from the file
    boolhist      = _util.GetHist(f_flugg, h_nuebar, "anueFluxHisto");
    if (boolhist == false) gSystem->Exit(0); 

    // Get the POT from the flux histogram
    POT_flux = 2.4e20;

    // Nue Integrated Flux
    xbin_th = h_nue->GetXaxis()->FindBin(_util.energy_threshold);                     // Find the x bin to integrate from
    integral_nue = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1); // Integrate over the flux for nue
    std::cout << "\nIntegral Nue Flux: " <<  integral_nue / POT_flux << " nue / POT / cm2" << std::endl;

    // Nuebar Integrated Flux
    xbin_th   = h_nuebar->GetXaxis()->FindBin(_util.energy_threshold);                                   // Find the x bin to integrate from
    integral_nuebar = h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1); // Integrate over the flux for nue
    std::cout << "\nIntegral Nuebar Flux: " << integral_nuebar / POT_flux<< " nuebar / POT / cm2" << "\n" << std::endl;

    // Print the flux scaled to the MC and Data POT
    std::cout << "Integral Flux MC POT: "   << mc_flux_scale_factor * integral_nuebar / (POT_flux * flux_scale_factor)   << " nu / MC POT / cm2"   << "\n" << std::endl;
    std::cout << "Integral Flux Data POT: " << data_flux_scale_factor * integral_nuebar / (POT_flux * flux_scale_factor) << " nu / Data POT / cm2" << "\n" << std::endl;

    // Return the flux per POT
    return (integral_nue + integral_nuebar) / (POT_flux * flux_scale_factor);

    f_flugg->Close();

    f_flux->cd();


}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFluxHP(int uni, std::string label){

    f_flux->cd();

    TH2D *h_uni_nue, *h_uni_nuebar;
    double xbin, ybin;
    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    double xbin_th{0.0};
    double integral_nue{0.0}, integral_nuebar{0.0};
    bool boolhist;
    
    // Get the POT from the file
    POT_flux = GetPOT(f_flux);

    // Get Nue histogram
    boolhist = _util.GetHist(f_flux, h_uni_nue, Form("nue/Multisims/nue_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
    if (boolhist == false) gSystem->Exit(0);
    
    // Get Nuebar histogram
    boolhist = _util.GetHist(f_flux, h_uni_nuebar, Form("nuebar/Multisims/nuebar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
    if (boolhist == false) gSystem->Exit(0);

    // Integrate the Nue flux histogram
    xbin_th = h_uni_nue->GetXaxis()->FindBin(_util.energy_threshold);                                    // Find the x bin to integrate from
    integral_nue = h_uni_nue->Integral( xbin_th, h_uni_nue->GetNbinsX()+1, 0, h_uni_nue->GetNbinsY()+1); // Integrate over the flux for nue
    // std::cout << "\nIntegral Nue Flux: " << flux_scale_factor * integral_nue / POT_flux << " nue / POT / GeV / cm2" << std::endl;

    // Integrate the Nuebar flux histogram
    xbin_th   = h_uni_nuebar->GetXaxis()->FindBin(_util.energy_threshold);                                           // Find the x bin to integrate from 
    integral_nuebar = h_uni_nuebar->Integral( xbin_th, h_uni_nuebar->GetNbinsX()+1, 0, h_uni_nuebar->GetNbinsY()+1); // Integrate over the flux for nuebar
    // std::cout << "Integral Nuebar Flux: " << flux_scale_factor * integral_nuebar / POT_flux << " nuebar / POT / GeV / cm2" << "\n" << std::endl;

    // Return the flux per POT
    return (integral_nue + integral_nuebar) / POT_flux;

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetBeamlineWeight(std::string variation, int _nu_pdg, double _true_energy, double _numi_ang){

    f_flux->cd();

    // double weight = {1.0};
    double xbin, ybin;
    double xbin_th{0.0};
    
    // Get the beamline variation index so we can get the right histogram to weight from
    int var_index = GetBeamlineIndex(variation);

    // Nue
    if (_nu_pdg == 12){
        xbin = beamline_hists.at(var_index).at(k_nue)->GetXaxis()->FindBin(_true_energy);
        ybin = beamline_hists.at(var_index).at(k_nue)->GetYaxis()->FindBin(_numi_ang);
        return ( beamline_hists.at(var_index).at(k_nue)->GetBinContent(xbin, ybin));
    }
    // Nuebar
    else if (_nu_pdg == -12){
        xbin = beamline_hists.at(var_index).at(k_nuebar)->GetXaxis()->FindBin(_true_energy);
        ybin = beamline_hists.at(var_index).at(k_nuebar)->GetYaxis()->FindBin(_numi_ang);
        return ( beamline_hists.at(var_index).at(k_nuebar)->GetBinContent(xbin, ybin));
    }
    // Numu
    else if (_nu_pdg == 14){
        xbin = beamline_hists.at(var_index).at(k_numu)->GetXaxis()->FindBin(_true_energy);
        ybin = beamline_hists.at(var_index).at(k_numu)->GetYaxis()->FindBin(_numi_ang);
        return ( beamline_hists.at(var_index).at(k_numu)->GetBinContent(xbin, ybin));
    }
    // NumuBar
    else if (_nu_pdg == -14){
        xbin = beamline_hists.at(var_index).at(k_numubar)->GetXaxis()->FindBin(_true_energy);
        ybin = beamline_hists.at(var_index).at(k_numubar)->GetYaxis()->FindBin(_numi_ang);
        return ( beamline_hists.at(var_index).at(k_numubar)->GetBinContent(xbin, ybin));
    }
    else return 1.0;


}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetHPWeight(int uni, std::string label, int _nu_pdg, double _true_energy, double _numi_ang){

    f_flux->cd();

    TH2D *h_nue, *h_nuebar ,*h_uni_nue, *h_uni_nuebar;
    TH2D *h_numu, *h_numubar, *h_uni_numu, *h_uni_numubar;
    double xbin, ybin;
    double xbinCV, ybinCV;
    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    double xbin_th{0.0};
    double integral_nue{0.0}, integral_nuebar{0.0};
    double ratio{1.0};
    bool boolhist;
    
    // Get the POT from the file
    POT_flux = GetPOT(f_flux);

  
    // Nue
    if (_nu_pdg == 12){

        // Get Nue histogram
        boolhist = _util.GetHist(f_flux, h_uni_nue, Form("nue/Multisims/nue_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        if (boolhist == false) gSystem->Exit(0);
    
        xbin = h_uni_nue->GetXaxis()->FindBin(_true_energy);
        ybin = h_uni_nue->GetYaxis()->FindBin(_numi_ang);
        ratio = h_uni_nue->GetBinContent(xbin, ybin);
        
        // Get the nue flux histogram from the file CV
        boolhist = _util.GetHist(f_flux, h_nue,         "nue/Detsmear/nue_CV_AV_TPC_2D");       
        if (boolhist == false) gSystem->Exit(0); 

        xbinCV = h_nue->GetXaxis()->FindBin(_true_energy);
        ybinCV = h_nue->GetYaxis()->FindBin(_numi_ang);
        ratio /= h_nue->GetBinContent(xbin, ybin);

        return (ratio);
    }
    // Nuebar
    else if (_nu_pdg == -12){
        
        // Get Nuebar histogram
        boolhist = _util.GetHist(f_flux, h_uni_nuebar, Form("nuebar/Multisims/nuebar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        if (boolhist == false) gSystem->Exit(0);
    
        xbin = h_uni_nuebar->GetXaxis()->FindBin(_true_energy);
        ybin = h_uni_nuebar->GetYaxis()->FindBin(_numi_ang);
        ratio = h_uni_nuebar->GetBinContent(xbin, ybin);
        
        // Get the nuebar flux histogram from the file CV
        boolhist = _util.GetHist(f_flux, h_nuebar,         "nuebar/Detsmear/nuebar_CV_AV_TPC_2D");       
        if (boolhist == false) gSystem->Exit(0); 

        xbinCV = h_nuebar->GetXaxis()->FindBin(_true_energy);
        ybinCV = h_nuebar->GetYaxis()->FindBin(_numi_ang);
        ratio /= h_nuebar->GetBinContent(xbin, ybin);

        return (ratio);
    }
    // Numu
    else if (_nu_pdg == 14){
        
        // Get Numu histogram
        boolhist = _util.GetHist(f_flux, h_uni_numu, Form("numu/Multisims/numu_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        if (boolhist == false) gSystem->Exit(0);
    
        xbin = h_uni_numu->GetXaxis()->FindBin(_true_energy);
        ybin = h_uni_numu->GetYaxis()->FindBin(_numi_ang);
        ratio = h_uni_numu->GetBinContent(xbin, ybin);
        
        // Get the numu flux histogram from the file CV
        boolhist = _util.GetHist(f_flux, h_numu,         "numu/Detsmear/numu_CV_AV_TPC_2D");       
        if (boolhist == false) gSystem->Exit(0); 

        xbinCV = h_numu->GetXaxis()->FindBin(_true_energy);
        ybinCV = h_numu->GetYaxis()->FindBin(_numi_ang);
        ratio /= h_numu->GetBinContent(xbin, ybin);

        return (ratio);
    }
    // NumuBar
    else if (_nu_pdg == -14){
        
        // Get Numubar histogram
        boolhist = _util.GetHist(f_flux, h_uni_numubar, Form("numubar/Multisims/numubar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        if (boolhist == false) gSystem->Exit(0);
    
        xbin = h_uni_numubar->GetXaxis()->FindBin(_true_energy);
        ybin = h_uni_numubar->GetYaxis()->FindBin(_numi_ang);
        ratio = h_uni_numubar->GetBinContent(xbin, ybin);
        
        // Get the numubar flux histogram from the file CV
        boolhist = _util.GetHist(f_flux, h_numubar,         "numubar/Detsmear/numubar_CV_AV_TPC_2D");       
        if (boolhist == false) gSystem->Exit(0); 

        xbinCV = h_numubar->GetXaxis()->FindBin(_true_energy);
        ybinCV = h_numubar->GetYaxis()->FindBin(_numi_ang);
        ratio /= h_numubar->GetBinContent(xbin, ybin);

        return (ratio);
    }
    else
        return 1.0;


}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetPOT(TFile* f){
    TTree* TPOT = (TTree*) f->Get("POT");
    if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

    double fPOT{0};
    TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
    TPOT->GetEntry(0);
    double total_entries = TPOT->GetEntries(); // if using hadd, this will not be 1 equal to 1 anymore
    fPOT*=total_entries;
    // std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;

    return fPOT;
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::WriteHists(){

    std::cout << "Now writing histograms to file!"<< std::endl;

    // Now open the output file
    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( Form("files/crosssec_run%s.root", _util.run_period ) )) {

        if (std::string(_util.xsec_bin_mode) == "ratio"){
            fnuexsec_out = new TFile( Form("files/crosssec_run%s_ratio.root", _util.run_period) , "UPDATE");
        }
        else {
            fnuexsec_out = new TFile( Form("files/crosssec_run%s.root", _util.run_period) , "UPDATE");
        }
    }

    fnuexsec_out->cd();

    // If the variation flag has been set then cheekily change the reeighter label name to use a different folder
    if (_util.isvariation && std::string(_util.variation) != "weight"){
        if ( std::string(_util.variation) == "CV"){
            std::cout << "Doing a cheeky swap of the CV name to: " << "detvar_CV" << std::endl;
            reweighter_labels = { "detvar_CV" };
        }
        else {
            std::cout << "Doing a cheeky swap of the CV name to: " << std::string(_util.variation) << std::endl;
            reweighter_labels = { std::string(_util.variation)};
        }
    }

    if (_util.isfakedata){

        std::string fold_name = "fake" + std::string(_util.fakedataname);
        std::cout << "Doing a cheeky swap of the CV name to: " << fold_name << std::endl;
            reweighter_labels = { fold_name };
    }

    // Create subdirectory for each reweighter
    TDirectory *dir_labels[reweighter_labels.size()];

    // Create subdirectory for each variable
    TDirectory *dir_labels_var[_util.vars.size()];

    // We dont want to overwrite these histograms with empty ones
    if (std::string(_util.xsec_rw_mode) != "rw_cuts"){
        
        // Loop over the labels
        for (unsigned int label = 0; label < reweighter_labels.size(); label++) {
            
            // See if the directory already exists
            bool bool_dir = _util.GetDirectory(fnuexsec_out, dir_labels[label], reweighter_labels.at(label).c_str());

            // If it doesnt exist then create it
            if (!bool_dir) dir_labels[label] = fnuexsec_out->mkdir(reweighter_labels.at(label).c_str());

            // Go into the directory
            dir_labels[label]->cd();

            // Loop over the universes
            for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++ ){

                // Loop over the variables
                for (unsigned int var = 0; var < h_cross_sec.at(label).at(uni).size(); var ++){

                    // See if the directory already exists
                    bool bool_dir = _util.GetDirectory(fnuexsec_out, dir_labels_var[var], Form("%s/%s", reweighter_labels.at(label).c_str(), _util.vars.at(var).c_str()));

                    // If it doesnt exist then create it
                    if (!bool_dir) dir_labels_var[var] = dir_labels[label]->mkdir(_util.vars.at(var).c_str());

                    // Go into the directory
                    dir_labels_var[var]->cd();
                
                    // Now write the histograms, 
                    for (unsigned int p = 0; p < h_cross_sec.at(label).at(uni).at(var).size(); p++){

                        // Certain histograms we want to divide out by the bin width
                        if ((var == k_var_recoX || var == k_var_trueX) && p != k_xsec_eff )
                            // h_cross_sec.at(label).at(uni).at(var).at(p)->Scale(1.0, "width");

                        h_cross_sec.at(label).at(uni).at(var).at(p)->SetOption("hist");
                        h_cross_sec.at(label).at(uni).at(var).at(p)->Write("",TObject::kOverwrite);
                    }

                    // Write the smearing matrix to the reco folder
                    if (var == k_var_recoX || var == k_var_trueX){
                        h_smear.at(label).at(uni).at(var)->SetOption("col,text00");
                        h_smear.at(label).at(uni).at(var)->Write(Form("h_run%s_%s_%i_smearing",_util.run_period, reweighter_labels.at(label).c_str(), uni),TObject::kOverwrite);
                    }

                    if (var == k_var_trueX && std::string(_util.xsec_bin_mode) == "e_ang" && reweighter_labels.at(label) == "CV" && uni == 0){
                            h_2D_CV_binning->SetOption("colz");
                            h_2D_CV_binning->Write("",TObject::kOverwrite);
                    }
                
                } // End loop over the variables

                fnuexsec_out->cd();    // change current directory to top
            
            } // End loop over universes

            fnuexsec_out->cd();    // change current directory to top

        }
    }
    else {

        // This is if we want to write the histograms by cut -- stole this cheeky bit of code stucture from Marina, thanks ;) !
        for (unsigned int label = 0; label < reweighter_labels.size(); label++) {
            
            // Loop over the Cuts
            for (unsigned int cut = 0; cut < h_cut_v.at(label).size(); cut++) {

                // Loop over the variables
                for (unsigned int var = 0; var < h_cut_v.at(label).at(cut).size(); var++) {

                    fnuexsec_out->cd();    // change current directory to top

                    if(!fnuexsec_out->GetDirectory(Form("%s/Cuts/%s/%s", reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), _util.vec_hist_name.at(var).c_str() ))) 
                        fnuexsec_out->mkdir(Form("%s/Cuts/%s/%s", reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), _util.vec_hist_name.at(var).c_str() )); // if the directory does not exist, create it
                    
                    fnuexsec_out->cd(Form("%s/Cuts/%s/%s", reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), _util.vec_hist_name.at(var).c_str() )); // open the directory

                    // loop over the universes
                    for (unsigned int uni = 0; uni < h_cut_v.at(label).at(cut).at(var).size(); uni++){
                        h_cut_v.at(label).at(cut).at(var).at(uni)->SetOption("hist");
                        h_cut_v.at(label).at(cut).at(var).at(uni)->Write("",TObject::kOverwrite);
                    }
                }
            }

        }
    }

    fnuexsec_out->cd();    // change current directory to top

    fnuexsec_out->Close();
    
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::InitTree(){
    
    // Set the tree branches
    tree->SetBranchAddress("run",    &run);
    tree->SetBranchAddress("subrun", &subrun);
    tree->SetBranchAddress("event",  &event);
    tree->SetBranchAddress("gen",    &gen);
    tree->SetBranchAddress("passed_selection",    &passed_selection);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("true_energy", &true_energy);
    tree->SetBranchAddress("reco_energy", &reco_energy);
    tree->SetBranchAddress("classification",   &classification);
    tree->SetBranchAddress("shr_energy_cali", &shr_energy_cali);
    tree->SetBranchAddress("elec_e",  &elec_e);
    tree->SetBranchAddress("ppfx_cv",  &ppfx_cv);
    tree->SetBranchAddress("weightSplineTimesTune",  &weightSplineTimesTune);
    tree->SetBranchAddress("numi_ang",  &numi_ang);
    tree->SetBranchAddress("nu_pdg",  &nu_pdg);
    tree->SetBranchAddress("npi0",  &npi0);
    tree->SetBranchAddress("pi0_e",  &pi0_e);
    tree->SetBranchAddress("effective_angle",            &effective_angle);
    tree->SetBranchAddress("cos_effective_angle",        &cos_effective_angle);
    tree->SetBranchAddress("true_effective_angle",       &true_effective_angle);
    tree->SetBranchAddress("cos_true_effective_angle",   &cos_true_effective_angle);
    tree->SetBranchAddress("ccnc",  &ccnc);
    
    tree->SetBranchAddress("weightsGenie",          &weightsGenie);
    tree->SetBranchAddress("weightsReint",          &weightsReint);
    tree->SetBranchAddress("weightsPPFX",           &weightsPPFX);
    tree->SetBranchAddress("knobRPAup",             &knobRPAup);
    tree->SetBranchAddress("knobRPAdn",             &knobRPAdn);
    tree->SetBranchAddress("knobCCMECup",           &knobCCMECup);
    tree->SetBranchAddress("knobCCMECdn",           &knobCCMECdn);
    tree->SetBranchAddress("knobAxFFCCQEup",        &knobAxFFCCQEup);
    tree->SetBranchAddress("knobAxFFCCQEdn",        &knobAxFFCCQEdn);
    tree->SetBranchAddress("knobVecFFCCQEup",       &knobVecFFCCQEup);
    tree->SetBranchAddress("knobVecFFCCQEdn",       &knobVecFFCCQEdn);
    tree->SetBranchAddress("knobDecayAngMECup",     &knobDecayAngMECup);
    tree->SetBranchAddress("knobDecayAngMECdn",     &knobDecayAngMECdn);
    tree->SetBranchAddress("knobThetaDelta2Npiup",  &knobThetaDelta2Npiup);
    tree->SetBranchAddress("knobThetaDelta2Npidn",  &knobThetaDelta2Npidn);
    tree->SetBranchAddress("knobThetaDelta2NRadup", &knobThetaDelta2NRadup);
    tree->SetBranchAddress("knobThetaDelta2NRaddn", &knobThetaDelta2NRaddn);
    tree->SetBranchAddress("knobRPA_CCQE_Reducedup",&knobRPA_CCQE_Reducedup);
    tree->SetBranchAddress("knobRPA_CCQE_Reduceddn",&knobRPA_CCQE_Reduceddn);
    tree->SetBranchAddress("knobNormCCCOHup",       &knobNormCCCOHup);
    tree->SetBranchAddress("knobNormCCCOHdn",       &knobNormCCCOHdn);
    tree->SetBranchAddress("knobNormNCCOHup",       &knobNormNCCOHup);
    tree->SetBranchAddress("knobNormNCCOHdn",       &knobNormNCCOHdn);
    tree->SetBranchAddress("knobxsr_scc_Fv3up",     &knobxsr_scc_Fv3up);
    tree->SetBranchAddress("knobxsr_scc_Fv3dn",     &knobxsr_scc_Fv3dn);
    tree->SetBranchAddress("knobxsr_scc_Fa3up",     &knobxsr_scc_Fa3up);
    tree->SetBranchAddress("knobxsr_scc_Fa3dn",     &knobxsr_scc_Fa3dn);

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::SwitchReweighterLabel(std::string label){

    // Clear it and go again
    vec_universes.clear();

    // Genie All
    if (label == "weightsGenie"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < weightsGenie->size(); j++){
            vec_universes.push_back( (double) weightsGenie->at(j)/1000.0 );
        }

    }
    // Geant Reinteractions
    else if (label == "weightsReint"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < weightsReint->size(); j++){
            vec_universes.push_back( (double) weightsReint->at(j)/1000.0 );
        }

    }
    // PPFX All
    else if (label == "weightsPPFX"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < weightsPPFX->size(); j++){
            vec_universes.push_back( (double) weightsPPFX->at(j)/1000.0);
        }

    }
    // MC Stats
    else if (label == "MCStats"){
        vec_universes.resize(uni_mcstats, 1.0);
    }
    else if (label == "RPAup"){
        vec_universes.push_back(knobRPAup);
    }
    else if (label == "CCMECup"){
        vec_universes.push_back(knobCCMECup);
    }
    else if (label == "AxFFCCQEup"){
        vec_universes.push_back(knobAxFFCCQEup);
    }
    else if (label == "VecFFCCQEup"){
        vec_universes.push_back(knobVecFFCCQEup);
    }
    else if (label == "DecayAngMECup"){
        vec_universes.push_back(knobDecayAngMECup);
    }
    else if (label == "ThetaDelta2Npiup"){
        vec_universes.push_back(knobThetaDelta2Npiup);
    }
    else if (label == "ThetaDelta2NRadup"){
        vec_universes.push_back(knobThetaDelta2NRadup);
    }
    else if (label == "RPA_CCQE_Reducedup"){
        vec_universes.push_back(knobRPA_CCQE_Reducedup);
    }
    else if (label == "NormCCCOHup"){
        vec_universes.push_back(knobNormCCCOHup);
    }
    else if (label == "NormNCCOHup"){
        vec_universes.push_back(knobNormNCCOHup);
    }
    else if (label == "RPAdn"){
        vec_universes.push_back(knobRPAdn);
    }
    else if (label == "CCMECdn"){
        vec_universes.push_back(knobCCMECdn);
    }
    else if (label == "AxFFCCQEdn"){
        vec_universes.push_back(knobAxFFCCQEdn);
    }
    else if (label == "VecFFCCQEdn"){
        vec_universes.push_back(knobVecFFCCQEdn);
    }
    else if (label == "DecayAngMECdn"){
        vec_universes.push_back(knobDecayAngMECdn);
    }
    else if (label == "ThetaDelta2Npidn"){
        vec_universes.push_back(knobThetaDelta2Npidn);
    }
    else if (label == "ThetaDelta2NRaddn"){
        vec_universes.push_back(knobThetaDelta2NRaddn);
    }
    else if (label == "RPA_CCQE_Reduceddn"){
        vec_universes.push_back(knobRPA_CCQE_Reduceddn);
    }
    else if (label == "NormCCCOHdn"){
        vec_universes.push_back(knobNormCCCOHdn);
    }
    else if (label == "NormNCCOHdn"){
        vec_universes.push_back(knobNormNCCOHdn);
    }
    else if (label == "xsr_scc_Fv3up"){
        vec_universes.push_back(knobxsr_scc_Fv3up);
    }
    else if (label == "xsr_scc_Fv3dn"){
        vec_universes.push_back(knobxsr_scc_Fv3dn);
    }
    else if (label == "xsr_scc_Fa3up"){
        vec_universes.push_back(knobxsr_scc_Fa3up);
    }
    else if (label == "xsr_scc_Fa3dn"){
        vec_universes.push_back(knobxsr_scc_Fa3dn);
    }
    // This can be the CV or any beamline variation
    else {
        vec_universes.push_back(1.0);
    }


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::SwitchReweighterLabel(std::string label, SliceContainer &SC){

    // Clear it and go again
    vec_universes.clear();

    // Genie All
    if (label == "weightsGenie"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < SC.weightsGenie->size(); j++){
            vec_universes.push_back( (double) SC.weightsGenie->at(j)/1000.0 );
        }

    }
    // Geant Reinteractions
    else if (label == "weightsReint"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < SC.weightsReint->size(); j++){
            vec_universes.push_back( (double) SC.weightsReint->at(j)/1000.0 );
        }

    }
    // PPFX All
    else if (label == "weightsPPFX"){
        
        // Convert from unsigned short to double and push back -- divide by 1000 to undo previous *1000
        for (unsigned int j = 0; j < SC.weightsPPFX->size(); j++){
            vec_universes.push_back( (double) SC.weightsPPFX->at(j)/1000.0);
        }

    }
    // MC Stats
    else if (label == "MCStats"){
        vec_universes.resize(uni_mcstats, 1.0);
    }
    else if (label == "RPAup"){
        vec_universes.push_back(SC.knobRPAup);
    }
    else if (label == "CCMECup"){
        vec_universes.push_back(SC.knobCCMECup);
    }
    else if (label == "AxFFCCQEup"){
        vec_universes.push_back(SC.knobAxFFCCQEup);
    }
    else if (label == "VecFFCCQEup"){
        vec_universes.push_back(SC.knobVecFFCCQEup);
    }
    else if (label == "DecayAngMECup"){
        vec_universes.push_back(SC.knobDecayAngMECup);
    }
    else if (label == "ThetaDelta2Npiup"){
        vec_universes.push_back(SC.knobThetaDelta2Npiup);
    }
    else if (label == "ThetaDelta2NRadup"){
        vec_universes.push_back(SC.knobThetaDelta2NRadup);
    }
    else if (label == "RPA_CCQE_Reducedup"){
        vec_universes.push_back(SC.knobRPA_CCQE_Reducedup);
    }
    else if (label == "NormCCCOHup"){
        vec_universes.push_back(SC.knobNormCCCOHup);
    }
    else if (label == "NormNCCOHup"){
        vec_universes.push_back(SC.knobNormNCCOHup);
    }
    else if (label == "RPAdn"){
        vec_universes.push_back(SC.knobRPAdn);
    }
    else if (label == "CCMECdn"){
        vec_universes.push_back(SC.knobCCMECdn);
    }
    else if (label == "AxFFCCQEdn"){
        vec_universes.push_back(SC.knobAxFFCCQEdn);
    }
    else if (label == "VecFFCCQEdn"){
        vec_universes.push_back(SC.knobVecFFCCQEdn);
    }
    else if (label == "DecayAngMECdn"){
        vec_universes.push_back(SC.knobDecayAngMECdn);
    }
    else if (label == "ThetaDelta2Npidn"){
        vec_universes.push_back(SC.knobThetaDelta2Npidn);
    }
    else if (label == "ThetaDelta2NRaddn"){
        vec_universes.push_back(SC.knobThetaDelta2NRaddn);
    }
    else if (label == "RPA_CCQE_Reduceddn"){
        vec_universes.push_back(SC.knobRPA_CCQE_Reduceddn);
    }
    else if (label == "NormCCCOHdn"){
        vec_universes.push_back(SC.knobNormCCCOHdn);
    }
    else if (label == "NormNCCOHdn"){
        vec_universes.push_back(SC.knobNormNCCOHdn);
    }
    else if (label == "xsr_scc_Fv3up"){
        vec_universes.push_back(SC.knobxsr_scc_Fv3up);
    }
    else if (label == "xsr_scc_Fv3dn"){
        vec_universes.push_back(SC.knobxsr_scc_Fv3dn);
    }
    else if (label == "xsr_scc_Fa3up"){
        vec_universes.push_back(SC.knobxsr_scc_Fa3up);
    }
    else if (label == "xsr_scc_Fa3dn"){
        vec_universes.push_back(SC.knobxsr_scc_Fa3dn);
    }
    // This can be the CV or any beamline variation
    else {
        vec_universes.push_back(1.0);
    }


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::InitialiseHistograms(std::string run_mode){

    // If mode is CV then set the reweighter labels to just the CV
    if (run_mode == "default"){
        std::cout << "Using the default mode, so resetting the reweighter label vector to just use the CV" << std::endl;
        reweighter_labels.clear();
        reweighter_labels.push_back("CV");
    }

    if (run_mode != "default") {
        std::cout << "XSec Labels to run: " << std::string(_util.xsec_labels) << std::endl;

        // Run Unisims only
        if (std::string(_util.xsec_labels) == "unisim"){
            std::cout << "XSec reweighting mode set to unisim, so only running unisims" << std::endl;
            reweighter_labels.clear();
            reweighter_labels = {
                "CV",    
                "Horn_p2kA",
                "Horn_m2kA",
                "Horn1_x_p3mm",
                "Horn1_x_m3mm",
                "Horn1_y_p3mm",
                "Horn1_y_m3mm",
                "Beam_spot_1_1mm",
                "Beam_spot_1_5mm",
                "Horn2_x_p3mm",
                "Horn2_x_m3mm",
                "Horn2_y_p3mm",
                "Horn2_y_m3mm",
                "Horns_0mm_water",
                "Horns_2mm_water",
                "Beam_shift_x_p1mm",
                "Beam_shift_x_m1mm",
                "Beam_shift_y_p1mm",
                "Beam_shift_y_m1mm",
                "Target_z_p7mm",
                "Target_z_m7mm",
                "Horn1_refined_descr",
                "Decay_pipe_Bfield",
                "Old_Horn_Geometry",
                "RPAup",
                "CCMECup",
                "AxFFCCQEup",
                "VecFFCCQEup",
                "DecayAngMECup",
                "ThetaDelta2Npiup",
                "ThetaDelta2NRadup",
                "RPA_CCQE_Reducedup",
                "NormCCCOHup",
                "NormNCCOHup",
                "RPAdn",
                "CCMECdn",
                "AxFFCCQEdn",
                "VecFFCCQEdn",
                "DecayAngMECdn",
                "ThetaDelta2Npidn",
                "ThetaDelta2NRaddn",
                "RPA_CCQE_Reduceddn",
                "NormCCCOHdn",
                "NormNCCOHdn",
                "xsr_scc_Fv3up",
                "xsr_scc_Fa3up",
                "xsr_scc_Fv3dn",
                "xsr_scc_Fa3dn",
                "Dirtup",
                "Dirtdn",
                "POTup",
                "POTdn", 
                "pi0",
                "EXT"
            };
        }
        // Only run PPFX
        else if (std::string(_util.xsec_labels) == "ppfx"){
            std::cout << "XSec reweighting mode set to ppfx" << std::endl;
            reweighter_labels.clear();
            reweighter_labels = {"CV", "weightsPPFX"};
        }
        // Only run Genie All
        else if (std::string(_util.xsec_labels) == "genie"){
            std::cout << "XSec reweighting mode set to genie all" << std::endl;
            reweighter_labels.clear();
            reweighter_labels = {"CV", "weightsGenie"};
        }
        // only run geant re-interactions
        else if (std::string(_util.xsec_labels) == "reint"){
            std::cout << "XSec reweighting mode set to reint" << std::endl;
            reweighter_labels.clear();
            reweighter_labels = {"CV", "weightsReint"};
        }
        // only run mc stats
        else if (std::string(_util.xsec_labels) == "mcstats"){
            std::cout << "XSec reweighting mode set to MCStats" << std::endl;
            reweighter_labels.clear();
            reweighter_labels = {"CV", "MCStats"};
        }
        // Run EVERYTHING
        else if (std::string(_util.xsec_labels) == "all"){
            std::cout << "Running everything!" << std::endl;
        }
        // Just default to whatever is configured in the header
        else {
            std::cout << "Unknown Cross section mode entered, running default mode" << std::endl;
            reweighter_labels = {"CV"};
        }
    }

    // Set the histogram bins

    // Electron/Shower Energy
    if (std::string(_util.xsec_var) =="elec_E"){
        hist_bins = _util.reco_shr_bins;
        fine_bins = _util.true_shr_bins;
    }
    // Electron/Shower effective angle
    else if (std::string(_util.xsec_var) =="elec_ang"){
        hist_bins = _util.reco_shr_bins_ang;
        fine_bins = _util.true_shr_bins_ang;
    }
    // Electron/Shower effective cangle
    else if (std::string(_util.xsec_var) =="elec_cang"){
        hist_bins = _util.reco_shr_bins_cang;
        fine_bins = _util.true_shr_bins_cang;
    }
    else {
        std::cout << "Unsupported parameter...exiting!" << std::endl;
        return;
    }

    // Set the bins here
    bins.resize(_util.vars.size());
    bins_fine.resize(_util.vars.size());

    // Integrated X-Section Bin definition
    bins.at(k_var_integrated)      = { 0.0, 1.1 };
    bins_fine.at(k_var_integrated) = { 0.0, 1.1 };

    // Reconstructed electron energy Bin definition
    bins.at(k_var_recoX)      = hist_bins;
    bins_fine.at(k_var_recoX) = fine_bins;

    // True electron energy Bin definition
    bins.at(k_var_trueX)      = hist_bins;
    bins_fine.at(k_var_trueX) = fine_bins;

    // Resize to the number of reweighters
    h_cross_sec.resize(reweighter_labels.size());
    h_smear.resize(reweighter_labels.size());
    
    // Resize each reweighter to their number of universes
    for (unsigned int j=0; j < reweighter_labels.size(); j++){

        // Specific resizing -- hardcoded and may break in the future
        if ( reweighter_labels.at(j) == "weightsGenie"){
            std::cout << "Setting Genie All Histogram universe vector to size: " << uni_genie << std::endl;
            h_cross_sec.at(j).resize(uni_genie);
            h_smear.at(j).resize(uni_genie);
            
        }
        // Specific resizing -- hardcoded and may break in the future
        else if ( reweighter_labels.at(j) == "weightsPPFX"){
            std::cout << "Setting PPFX All Histogram universe vector to size: " << uni_ppfx << std::endl;
            h_cross_sec.at(j).resize(uni_ppfx);
            h_smear.at(j).resize(uni_ppfx);
            
        }
        // Specific resizing -- hardcoded and may break in the future
        else if ( reweighter_labels.at(j) == "weightsReint" ){
            std::cout << "Setting Geant Reinteractions Histogram universe vector to size: " << uni_reint << std::endl;
            h_cross_sec.at(j).resize(uni_reint);
            h_smear.at(j).resize(uni_reint);
            
        }
         // Specific resizing -- hardcoded and may break in the future
        else if ( reweighter_labels.at(j) == "MCStats" ){
            std::cout << "Setting MCStats Histogram universe vector to size: " << uni_mcstats << std::endl;
            h_cross_sec.at(j).resize(uni_mcstats);
            h_smear.at(j).resize(uni_mcstats);
            
        }
        // Default size of 1
        else {
            h_cross_sec.at(j).resize(1);
            h_smear.at(j).resize(1);
        }

    }

    // Now resize the universe to each type of variable we want to plot/reweight
    
    // Resize each reweighter to their number of universes
    for (unsigned int label=0; label < reweighter_labels.size() ;label++){

        for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){
            // Resize the histogram vector. plot var, cuts, classifications
            h_cross_sec.at(label).at(uni).resize(_util.vars.size());
            h_smear.at(label).at(uni).resize(_util.vars.size());
        }

    }

    // Now resize the universe to each type of histogram we want to plot/reweight
    
    // Resize each reweighter to their number of universes
    for (unsigned int j=0; j < reweighter_labels.size(); j++){

        // Loop over the universes
        for (unsigned int y=0; y < h_cross_sec.at(j).size(); y++){
            
            // Loop over the variables
            for (unsigned int var = 0; var < h_cross_sec.at(j).at(y).size(); var ++){
                
                // Resize the histogram vector. plot var, cuts, classifications
                h_cross_sec.at(j).at(y).at(var).resize(xsec_types.size());

            }
            
        }

    }

    std::vector<double> temp_bins;
    std::vector<double> temp_bins_fine;


    // Loop over the rewighters
    for (unsigned int label=0; label < reweighter_labels.size(); label++){

        // loop over the universes
        for (unsigned int uni=0; uni < h_cross_sec.at(label).size(); uni++){

            // Define bin labels fro smearing matrix
            int nbins_smear;
            double* edges_smear;
            int nbins_smear_fine;
            double* edges_smear_fine;
            
            // Loop over the variables
            for (unsigned int var = 0; var < h_cross_sec.at(label).at(uni).size(); var ++){

                if (std::string(_util.xsecmode) == "txtlist")
                    continue;

                // Get the number of bins and the right vector
                int const nbins = bins.at(var).size()-1;
                temp_bins.clear();
                temp_bins = bins.at(var);
                double* edges = &temp_bins[0]; // Cast to an array 

                int const nbins_fine = bins_fine.at(var).size()-1;
                temp_bins_fine.clear();
                temp_bins_fine = bins_fine.at(var);
                double* edges_fine = &temp_bins_fine[0]; // Cast to an array 

                if (var == k_var_recoX){
                    nbins_smear = nbins;
                    edges_smear = edges;
                    nbins_smear_fine = nbins_fine;
                    edges_smear_fine = edges_fine;
                }

                // loop over and create the histograms
                for (unsigned int i=0; i < xsec_types.size();i++){    
                    if (i == k_xsec_sel || i == k_xsec_bkg || i == k_xsec_gen || i == k_xsec_gen_smear || i == k_xsec_sig || i == k_xsec_ext || i == k_xsec_dirt || i == k_xsec_data || i == k_xsec_gen_shape){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, _util.vars.at(var).c_str(), xsec_types.at(i).c_str()) ,Form("%s", _util.var_labels_events.at(var).c_str()), nbins, edges);
                    }
                    else if (i == k_xsec_eff){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, _util.vars.at(var).c_str(), xsec_types.at(i).c_str()) ,Form("%s", _util.var_labels_eff.at(var).c_str()), nbins, edges);
                    }
                    else if (i == k_xsec_dataxsec || i == k_xsec_mcxsec || i == k_xsec_mcxsec_smear || i == k_xsec_mcxsec_shape){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, _util.vars.at(var).c_str(), xsec_types.at(i).c_str()) ,Form("%s", _util.var_labels_xsec.at(var).c_str()), nbins, edges);
                    }
                    else {
                        // If this is the case then there is an uncaught case
                        std::cout << "Houston we have a problem!" << std::endl;
                    }
                }

                // We dont care about the integrated smearing, so set it to some arbitary value
                if (var == k_var_integrated){
                    h_smear.at(label).at(uni).at(var)      = new TH2D ( Form("h_run%s_%s_%i_%s_smearing",_util.run_period, reweighter_labels.at(label).c_str(), uni, _util.vars.at(var).c_str()), _util.smear_hist_name.c_str(), 1, 0, 1, 1, 0, 1);
                    
                }
                // Set the reco and true bins
                else {
                    h_smear.at(label).at(uni).at(var)      = new TH2D ( Form("h_run%s_%s_%i_%s_smearing",_util.run_period, reweighter_labels.at(label).c_str(), uni, _util.vars.at(var).c_str()), _util.smear_hist_name.c_str(), nbins_smear, edges_smear, nbins_smear, edges_smear);
                }
                

            } // End loop over the variables

        
        } // End loop over the universes
    
    } // End loop over the labels


    // Create histogram for binning
    if (std::string(_util.xsec_bin_mode) == "e_ang"){
        // Get the number of bins and the right vector
        int const nbins_E = _util.true_shr_bins.size()-1;
        std::vector<double> temp_bins_E = _util.true_shr_bins;
        double* edges_E = &temp_bins_E[0]; // Cast to an array 
        int const nbins_ang = _util.true_shr_bins_ang.size()-1;
        std::vector<double>temp_bins_ang = _util.true_shr_bins_ang;
        double* edges_ang = &temp_bins_ang[0]; // Cast to an array 
        h_2D_CV_binning = new TH2D( "h_2D_CV_binning", ";E^{true}_{e#lower[-0.5]{-} + e^{+}} [GeV]; #beta^{true}_{e#lower[-0.5]{-} + e^{+}} [deg]", nbins_E, edges_E, nbins_ang, edges_ang );
    }


    // Now lets create the hostogram vector for the cuts

    // Resize to the number of reweighters
    h_cut_v.resize(reweighter_labels.size());

    // Resize the cut vec to the number of cuts
    for (unsigned int label = 0; label < h_cut_v.size(); label++){
        h_cut_v.at(label).resize(_util.k_cuts_MAX);
    }

    // reszie to the numner of variables
    for (unsigned int label = 0; label < h_cut_v.size(); label++){
        for (unsigned int cut = 0; cut < h_cut_v.at(label).size(); cut++){
            h_cut_v.at(label).at(cut).resize(_util.k_cut_vars_max);
        }
    }

    // Loop over the labels
    for (unsigned int label = 0; label < h_cut_v.size(); label++){
        
        // Loop over the cuts
        for (unsigned int cut = 0; cut < h_cut_v.at(label).size(); cut++){
            
            // Loop over the variables
            for (unsigned int var=0; var < h_cut_v.at(label).at(cut).size(); var++){

                // Now resize by the universes

                // Specific resizing -- hardcoded and may break in the future
                if ( reweighter_labels.at(label) == "weightsGenie"){
                    h_cut_v.at(label).at(cut).at(var).resize(uni_genie);
                }
                // Specific resizing -- hardcoded and may break in the future
                else if ( reweighter_labels.at(label) == "weightsPPFX"){
                    h_cut_v.at(label).at(cut).at(var).resize(uni_ppfx);
                }
                // Specific resizing -- hardcoded and may break in the future
                else if ( reweighter_labels.at(label) == "weightsReint" ){
                    h_cut_v.at(label).at(cut).at(var).resize(uni_reint);
                }
                // Specific resizing -- hardcoded and may break in the future
                else if ( reweighter_labels.at(label) == "MCStats" ){
                    h_cut_v.at(label).at(cut).at(var).resize(uni_mcstats);
                }
                // Default size of 1
                else {
                    h_cut_v.at(label).at(cut).at(var).resize(1);
                }

            }
        }
    }


    // Loop over the labels
    for (unsigned int label = 0; label < h_cut_v.size(); label++){
        
        // Loop over the cuts
        for (unsigned int cut = 0; cut < h_cut_v.at(label).size(); cut++){
            
            // Loop over the universes
            for (unsigned int uni=0; uni < h_cut_v.at(label).at(cut).at(0).size(); uni++){
                
                // Only initialise if this mode to stop loading too much into memory
                if (std::string(_util.xsec_rw_mode) == "rw_cuts"){
                    h_cut_v.at(label).at(cut).at(_util.k_cut_nslice).at(uni)                         = new TH1D(Form("h_reco_nslice_%s_%s_%i",                         reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 2, 0, 2);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shower_multiplicity).at(uni)            = new TH1D(Form("h_reco_shower_multiplicity_%s_%s_%i",            reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 6, 0, 6);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_track_multiplicity).at(uni)             = new TH1D(Form("h_reco_track_multiplicity_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 6, 0, 6);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_topological_score).at(uni)              = new TH1D(Form("h_reco_topological_score_%s_%s_%i",              reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 30, 0, 1);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_vtx_x_sce).at(uni)                      = new TH1D(Form("h_reco_vtx_x_sce_%s_%s_%i",                      reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 15, -10, 270);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_vtx_y_sce).at(uni)                      = new TH1D(Form("h_reco_vtx_y_sce_%s_%s_%i",                      reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 30, -120, 120);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_vtx_z_sce).at(uni)                      = new TH1D(Form("h_reco_vtx_z_sce_%s_%s_%i",                      reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 30, -10, 1050);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shower_score).at(uni)                   = new TH1D(Form("h_reco_shower_score_%s_%s_%i",                   reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 20, 0, 0.5);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shr_tkfit_dedx_max).at(uni)             = new TH1D(Form("h_reco_shr_tkfit_dedx_max_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 40, 0, 10);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shr_tkfit_dedx_max_tune).at(uni)        = new TH1D(Form("h_reco_shr_tkfit_dedx_max_tune_%s_%s_%i",        reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 20, 0, 10);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shr_tkfit_dedx_max_with_tracks).at(uni) = new TH1D(Form("h_reco_shr_tkfit_dedx_max_with_tracks_%s_%s_%i", reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 40, 0, 10);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shr_tkfit_dedx_max_no_tracks).at(uni)   = new TH1D(Form("h_reco_shr_tkfit_dedx_max_no_tracks_%s_%s_%i",   reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 20, 0, 10);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shower_to_vtx_dist).at(uni)             = new TH1D(Form("h_reco_shower_to_vtx_dist_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 20, 0, 20);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_hits_ratio).at(uni)                     = new TH1D(Form("h_reco_hits_ratio_%s_%s_%i",                     reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 21, 0, 1.05);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_CosmicIPAll3D).at(uni)                  = new TH1D(Form("h_reco_CosmicIPAll3D_%s_%s_%i",                  reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 40, 0, 200);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_contained_fraction).at(uni)             = new TH1D(Form("h_reco_contained_fraction_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 21, 0, 1.05);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shrmoliereavg).at(uni)                  = new TH1D(Form("h_reco_shrmoliereavg_%s_%s_%i",                  reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 30, 0, 30);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_leading_shower_theta).at(uni)           = new TH1D(Form("h_reco_leading_shower_theta_%s_%s_%i",           reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 13, 0, 190);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_leading_shower_phi).at(uni)             = new TH1D(Form("h_reco_leading_shower_phi_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 14, -190, 190);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shower_energy_cali).at(uni)             = new TH1D(Form("h_reco_shower_energy_cali_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 20, 0, 4);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_flash_time).at(uni)                     = new TH1D(Form("h_reco_flash_time_%s_%s_%i",                     reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 50, 0, 25);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_flash_pe).at(uni)                       = new TH1D(Form("h_reco_flash_pe_%s_%s_%i",                       reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 25, 0, 5000);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_effective_angle).at(uni)                = new TH1D(Form("h_reco_effective_angle_%s_%s_%i",                reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 13, 0, 190);
                    h_cut_v.at(label).at(cut).at(_util.k_cut_effective_cosangle).at(uni)             = new TH1D(Form("h_reco_effective_cosangle_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", 16, -1, 1);


                    double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
                    h_cut_v.at(label).at(cut).at(_util.k_cut_shower_energy_cali_rebin).at(uni)  = new TH1D(Form("h_reco_shower_energy_cali_rebin_%s_%s_%i",  reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", _util.reco_shr_bins.size()-1, edges);

                    edges = &_util.reco_shr_bins_ang[0]; // Cast to an array 
                    h_cut_v.at(label).at(cut).at(_util.k_cut_effective_angle_rebin).at(uni)                = new TH1D(Form("h_reco_effective_angle_rebin_%s_%s_%i",                reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", _util.reco_shr_bins_ang.size()-1, edges);
                    
                    edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
                    h_cut_v.at(label).at(cut).at(_util.k_cut_effective_cosangle_rebin).at(uni)             = new TH1D(Form("h_reco_effective_cosangle_rebin_%s_%s_%i",             reweighter_labels.at(label).c_str(), _util.cut_dirs.at(cut).c_str(), uni), "", _util.reco_shr_bins_cang.size()-1, edges);
                }
            }
        }
    }



    std::cout << "Initialisation of cross-section histograms is complete!" << std::endl;


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::FillHists(int label, int uni, int xsec_type, double weight_uni, float _recoX, float _trueX){

    // Integrated
    h_cross_sec.at(label).at(uni).at(k_var_integrated).at(xsec_type)->Fill(1.0, weight_uni);

    // Reco Electron Energy
    h_cross_sec.at(label).at(uni).at(k_var_recoX).at(xsec_type)->Fill(_recoX, weight_uni);

    // True Electron Energy, but dont fill the bkg
    if (xsec_type != k_xsec_bkg){
        h_cross_sec.at(label).at(uni).at(k_var_trueX).at(xsec_type)->Fill(_trueX, weight_uni);
    }
    

    // Smearing Matrix -- only fill for selected signal events
    if (xsec_type == k_xsec_sig){
        h_smear.at(label).at(uni).at(k_var_recoX)->Fill(_trueX, _recoX, weight_uni);
        h_smear.at(label).at(uni).at(k_var_trueX)->Fill(_trueX, _recoX, weight_uni);
    }

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::GetBeamlineHists(){

    std::cout  << "Initialising Beamline Histograms" << std::endl;

    // Resize the vector
    beamline_hists.resize(beamline_map.size());

    for (unsigned int f = 0; f < beamline_map.size(); f++){
        beamline_hists.at(f).resize(k_flav_max);
    }

    // resize the beamline flux vector
    beamline_flux.resize(beamline_map.size(), 0.0);


    // Loop over the beamline files and the ratio of beamline histograms to the cv for each
    // flavour along with the integrated flux. This will allow us to reweight the events
    for (unsigned int f = 0; f < beamline_map.size(); f++){
        TFile *f_temp;
        TH2D *htemp;
        TH2D* h_cv;

        // Get the beamline file
        if (std::string(_util.run_period) == "1"){
            
            // Select the file depending on whether we are on the gpvm or not
            if (!_util.use_gpvm)
                f_temp = TFile::Open(Form("Systematics/beamline/FHC/output_uboone_fhc_run%i.root", beamline_map.at(f).second));
            else
                f_temp = TFile::Open(Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/FHC/output_uboone_fhc_run%i.root", beamline_map.at(f).second));
        
        }
        else if(std::string(_util.run_period) == "3") {
            
            // Select the file depending on whether we are on the gpvm or not
            if (!_util.use_gpvm)
                f_temp = TFile::Open(Form("Systematics/beamline/RHC/output_uboone_rhc_run%i.root", beamline_map.at(f).second));
            else
                f_temp = TFile::Open(Form("/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold_v46/RHC/output_uboone_rhc_run%i.root", beamline_map.at(f).second));
        
        }
        else {
            std::cout << "Unknown run period, exiting...."<< std::endl;
            exit(1);
        }

        if (f_temp == NULL){
            std::cout << "Error in opening beamline file!, EXITING"<< std::endl;
            exit(2);
        }


        // -------------------------------------------------------------------------
        
        // Get the NuMu histogram from the beamline file
        htemp = (TH2D*)f_temp->Get("numu/Detsmear/numu_CV_AV_TPC_2D");
        beamline_hists.at(f).at(k_numu) = (TH2D*) htemp->Clone(Form("%s_%i", beamline_map.at(f).first.c_str(), k_numu));
        beamline_hists.at(f).at(k_numu)->Scale(1.0/GetPOT(f_temp)); 
        beamline_hists.at(f).at(k_numu)->SetDirectory(f_flux);
        
        // Now divide it by the CV
        f_flux->cd();
        _util.GetHist(f_flux, h_cv,         "numu/Detsmear/numu_CV_AV_TPC_2D");      
        beamline_hists.at(f).at(k_numu)->Divide(h_cv);
        f_temp->cd();

        // -------------------------------------------------------------------------
        
        // Get the NuMubar histogram from the beamline file
        htemp = (TH2D*)f_temp->Get("numubar/Detsmear/numubar_CV_AV_TPC_2D");
        beamline_hists.at(f).at(k_numubar) = (TH2D*) htemp->Clone(Form("%s_%i", beamline_map.at(f).first.c_str(), k_numubar));
        beamline_hists.at(f).at(k_numubar)->Scale(1.0/GetPOT(f_temp)); 
        beamline_hists.at(f).at(k_numubar)->SetDirectory(f_flux);

        // Now divide it by the CV
        f_flux->cd();
        _util.GetHist(f_flux, h_cv,         "numubar/Detsmear/numubar_CV_AV_TPC_2D");       
        beamline_hists.at(f).at(k_numubar)->Divide(h_cv);
        f_temp->cd();

        // -------------------------------------------------------------------------

        // Get the Nue histogram from the beamline file
        htemp = (TH2D*)f_temp->Get("nue/Detsmear/nue_CV_AV_TPC_2D");
        beamline_hists.at(f).at(k_nue) = (TH2D*) htemp->Clone(Form("%s_%i", beamline_map.at(f).first.c_str(), k_nue));
        beamline_hists.at(f).at(k_nue)->Scale(1.0/GetPOT(f_temp)); 
        beamline_hists.at(f).at(k_nue)->SetDirectory(f_flux);

        // Get the flux before dividing
        double xbin_th   = beamline_hists.at(f).at(k_nue)->GetXaxis()->FindBin(_util.energy_threshold);              // find the x bin to integrate from 
        beamline_flux.at(f) += beamline_hists.at(f).at(k_nue)->Integral( xbin_th, beamline_hists.at(f).at(k_nue)->GetNbinsX()+1, 0, beamline_hists.at(f).at(k_nue)->GetNbinsY()+1); // Integrate over the flux
    
        // Now divide it by the CV
        f_flux->cd();
        _util.GetHist(f_flux, h_cv,         "nue/Detsmear/nue_CV_AV_TPC_2D");       
        beamline_hists.at(f).at(k_nue)->Divide(h_cv);
        f_temp->cd();

        // -------------------------------------------------------------------------

        // Get the Nuebar histogram from the beamline file
        htemp = (TH2D*)f_temp->Get("nuebar/Detsmear/nuebar_CV_AV_TPC_2D");
        beamline_hists.at(f).at(k_nuebar) = (TH2D*) htemp->Clone(Form("%s_%i", beamline_map.at(f).first.c_str(), k_nuebar));
        beamline_hists.at(f).at(k_nuebar)->Scale(1.0/GetPOT(f_temp)); 
        beamline_hists.at(f).at(k_nuebar)->SetDirectory(f_flux);
        
        // Get the flux before dividing
        xbin_th   = beamline_hists.at(f).at(k_nuebar)->GetXaxis()->FindBin(_util.energy_threshold);              // find the x bin to integrate from 
        beamline_flux.at(f) += beamline_hists.at(f).at(k_nuebar)->Integral( xbin_th, beamline_hists.at(f).at(k_nuebar)->GetNbinsX()+1, 0, beamline_hists.at(f).at(k_nuebar)->GetNbinsY()+1); // Integrate over the flux

        // Now divide it by the CV
        f_flux->cd();
        _util.GetHist(f_flux, h_cv,         "nuebar/Detsmear/nuebar_CV_AV_TPC_2D");       
        beamline_hists.at(f).at(k_nuebar)->Divide(h_cv);
        f_temp->cd();

        // std::cout << beamline_flux.at(f) << std::endl;

        f_temp->Close();

        // -------------------------------------------------------------------------
    }

}
// -----------------------------------------------------------------------------
bool CrossSectionHelper::CheckBeamline(std::string variation){

    if (variation == "Horn_p2kA"           ||
        variation == "Horn_m2kA"           ||
        variation == "Horn1_x_p3mm"        ||
        variation == "Horn1_x_m3mm"        ||
        variation == "Horn1_y_p3mm"        ||
        variation == "Horn1_y_m3mm"        ||
        variation == "Beam_spot_1_1mm"     ||
        variation == "Beam_spot_1_5mm"     ||
        variation == "Horn2_x_p3mm"        ||
        variation == "Horn2_x_m3mm"        ||
        variation == "Horn2_y_p3mm"        ||
        variation == "Horn2_y_m3mm"        ||
        variation == "Horns_0mm_water"     ||
        variation == "Horns_2mm_water"     ||
        variation == "Beam_shift_x_p1mm"   ||
        variation == "Beam_shift_x_m1mm"   ||
        variation == "Beam_shift_y_p1mm"   ||
        variation == "Beam_shift_y_m1mm"   ||
        variation == "Target_z_p7mm"       ||
        variation == "Target_z_m7mm"       ||
        variation == "Horn1_refined_descr" ||
        variation == "Decay_pipe_Bfield"   ||
        variation == "Old_Horn_Geometry") {
        return true;
    }
    else {
        return false;
    }

}
// -----------------------------------------------------------------------------
int CrossSectionHelper::GetBeamlineIndex(std::string variation){

    if       (variation == "Horn_p2kA"          ) return k_Horn_p2kA;
    else if  (variation == "Horn_m2kA"          ) return k_Horn_m2kA;
    else if  (variation == "Horn1_x_p3mm"       ) return k_Horn1_x_p3mm;
    else if  (variation == "Horn1_x_m3mm"       ) return k_Horn1_x_m3mm;
    else if  (variation == "Horn1_y_p3mm"       ) return k_Horn1_y_p3mm;
    else if  (variation == "Horn1_y_m3mm"       ) return k_Horn1_y_m3mm;
    else if  (variation == "Beam_spot_1_1mm"    ) return k_Beam_spot_1_1mm;
    else if  (variation == "Beam_spot_1_5mm"    ) return k_Beam_spot_1_5mm;
    else if  (variation == "Horn2_x_p3mm"       ) return k_Horn2_x_p3mm;
    else if  (variation == "Horn2_x_m3mm"       ) return k_Horn2_x_m3mm;
    else if  (variation == "Horn2_y_p3mm"       ) return k_Horn2_y_p3mm;
    else if  (variation == "Horn2_y_m3mm"       ) return k_Horn2_y_m3mm;
    else if  (variation == "Horns_0mm_water"    ) return k_Horns_0mm_water;
    else if  (variation == "Horns_2mm_water"    ) return k_Horns_2mm_water;
    else if  (variation == "Beam_shift_x_p1mm"  ) return k_Beam_shift_x_p1mm;
    else if  (variation == "Beam_shift_x_m1mm"  ) return k_Beam_shift_x_m1mm;
    else if  (variation == "Beam_shift_y_p1mm"  ) return k_Beam_shift_y_p1mm;
    else if  (variation == "Beam_shift_y_m1mm"  ) return k_Beam_shift_y_m1mm;
    else if  (variation == "Target_z_p7mm"      ) return k_Target_z_p7mm;
    else if  (variation == "Target_z_m7mm"      ) return k_Target_z_m7mm;
    else if  (variation == "Horn1_refined_descr") return k_Horn1_refined_descr;
    else if  (variation == "Decay_pipe_Bfield"  ) return k_Decay_pipe_Bfield;
    else if  (variation == "Old_Horn_Geometry"  ) return k_Old_Horn_Geometry;
    else return 0; 

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::PoissonRandomNumber(int uni){

    // Concatenate the run, subrun, event to a single number so we have a seed
    int seed = ConcatRunSubRunEvent(run, subrun, event);

    // Set the seed of the TRandom 3 based on the run,subrun,event
    // Only set if uni == 0, i.e. a new event
    if (uni == 0) rand->SetSeed(seed);

    // Generate the weight, using a poisson dist with mean 1
    double weight_poisson = rand->Poisson(1);

    // std::cout << weight << std::endl;

    return weight_poisson;

}
// -----------------------------------------------------------------------------
int CrossSectionHelper::ConcatRunSubRunEvent(int run, int subrun, int event){ 
  
    // Convert both the integers to string 
    std::string srun    = std::to_string(run); 
    std::string ssubrun = std::to_string(subrun);
    std::string sevent  = std::to_string(event);

    // Concatenate the subrun and event. Dont add the run because it makes the number too long for storing as an int
    std::string s =  ssubrun + sevent; 

    // std::cout << srun << "  " << ssubrun << "  " << sevent<< std::endl;
  
    // Convert the concatenated string 
    // to integer 
    int c = stoi(s); 
  
    // return the formed integer 
    return c; 
} 
// -----------------------------------------------------------------------------
void CrossSectionHelper::Smear(TH1D* h_sig, TH1D* h_gen, TH2D* h_smear, TH1D* h_eff){

    // Change this bool to print info
    bool debug = false;

    //  ------ First normalise the smearing matrix ------
    // Loop over rows
    for (int row=1; row<h_smear->GetXaxis()->GetNbins()+2; row++) {
        double integral = 0;

        // Loop over columns and get the integral
        for (int col=1; col<h_smear->GetYaxis()->GetNbins()+2; col++){
            integral+=h_smear->GetBinContent(row, col);            
        }

        // Now normalise the column entries by the integral
        for (int col=1; col<h_smear->GetYaxis()->GetNbins()+2; col++){
            
            if (integral == 0){
                h_smear->SetBinContent(row,col, 0.0 );
            }
            else{
                h_smear->SetBinContent(row,col, h_smear->GetBinContent(row, col)/ integral );
            }
            
        }
    } 

    double nbins;

    // Electron/Shower Energy
    if (std::string(_util.xsec_var) =="elec_E"){
        nbins = _util.reco_shr_bins.size()+1;
    }
    // Electron/Shower beta
    else if (std::string(_util.xsec_var) =="elec_ang"){
        nbins = _util.reco_shr_bins_ang.size()+1;
    }
    // Electron/Shower cos(beta)
    else if (std::string(_util.xsec_var) =="elec_cang"){
        nbins = _util.reco_shr_bins_cang.size()+1;
    }

    // Convert the smearing matrix and histograms to a Matrix
    TMatrixD SMatrix(nbins, nbins, h_smear->GetArray(), "F");
    TMatrixD SigMatrix(nbins, 1, h_sig->GetArray(), "F");
    TMatrixD GenMatrix(nbins, 1, h_gen->GetArray(), "F");
    

    // Smear the signal selected events
    TMatrixD SigMatrixSmear(nbins, 1);
    SigMatrixSmear.Mult(SMatrix,SigMatrix);
    
    // Smear the generated events
    TMatrixD GenMatrixSmear(nbins, 1);
    GenMatrixSmear.Mult(SMatrix,GenMatrix);

    if (debug){
        SMatrix.Print();
        SigMatrix.Print();
        SigMatrixSmear.Print();
        GenMatrix.Print();
        GenMatrixSmear.Print();
    }

    // Set the efficiency bins
    for (int bin = 1; bin < h_sig->GetNbinsX()+1; bin++){
        
        if (GenMatrixSmear(bin,0) == 0)
            h_eff->SetBinContent(bin, 0.0);
        else
            h_eff->SetBinContent(bin, SigMatrixSmear(bin,0) / GenMatrixSmear(bin,0));
        
        h_eff->SetBinError(bin, 0);
        
        if (debug)
            std::cout << 100 * SigMatrix(bin,0) / GenMatrix(bin,0) << "  "<< 100 * SigMatrixSmear(bin,0) / GenMatrixSmear(bin,0) << std::endl;
    }
    
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::ApplyResponseMatrix(TH1D* h_gen, TH1D* h_gen_smear, TH1D* h_gen_CV, TH2D* h_smear, bool norm){

    // Change this bool to print info
    bool debug = false;

    // If we want to normalise the response matrix
    if (norm){
        //  ------ First normalise the smearing matrix ------
        for (int col=1; col<h_smear->GetYaxis()->GetNbins()+2; col++){
            double integral = 0;

            // Now normalise the column entries by the number of events in the 1D generated histogram
            for (int row=1; row<h_smear->GetXaxis()->GetNbins()+2; row++) { 
                
                if (h_gen->GetBinContent(row) == 0)
                    h_smear->SetBinContent(row,col, 0.0 );
                else {
                    // h_smear->SetBinContent(row,col, h_gen->GetBinWidth(row) * h_smear->GetBinContent(row, col)/ h_gen->GetBinContent(row) );
                    h_smear->SetBinContent(row,col, h_smear->GetBinContent(row, col)/ h_gen->GetBinContent(row) );
                }
            }
        } 
    }


    // Clear the Bins
    for (int bin = 0; bin < h_gen_smear->GetNbinsX()+2; bin++){
        h_gen_smear->SetBinContent(bin, 0);
    }

    // --- Do the matrix multiplication with the CV gen events --- 
    // Loop over cols
    for (int i=1; i<h_smear->GetYaxis()->GetNbins()+2; i++){
        double integral = 0;

        // Now normalise the column entries by the number of events in the 1D generated histogram
        for (int j=1; j<h_smear->GetXaxis()->GetNbins()+2; j++) { 

            if (debug)
                std::cout <<  "R_" << j << i << " * " << j << "  " << h_smear->GetBinContent(j, i) << " * " << h_gen_CV->GetBinContent(j) << std::endl;
            
            // h_gen_smear->SetBinContent(i, h_gen_smear->GetBinContent(i) + h_smear->GetBinContent(j, i) * (h_gen_CV->GetBinContent(j)/ h_gen_CV->GetBinWidth(j)));
            h_gen_smear->SetBinContent(i, h_gen_smear->GetBinContent(i) + h_smear->GetBinContent(j, i) * (h_gen_CV->GetBinContent(j)));
        }

        if (debug)
            std::cout << std::endl;
    } 

    // Set the bin errors to zero
    for (int bin = 1; bin < h_gen_smear->GetNbinsX()+1; bin++){
        h_gen_smear->SetBinError(bin, 0);
    }
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::SaveEvent(std::string _classification, bool _passed_selection, std::vector<float> _ev_weight, double recoX, double trueX){

    // Get rid of the CV
    _ev_weight.erase(_ev_weight.begin());

    // If the CV weight is zero then lets not fill the event (which is all zeros)
    if (weight == 0)
        return;

    // Signal event
    if ((_classification == "nue_cc" || _classification == "nuebar_cc" || _classification == "unmatched_nue" || _classification == "unmatched_nuebar") && _passed_selection) {
        
        isSignal = true;
        isSelected = true;

        // Initialise the file 
        if (!filled_sig){
            evt_dist_sig << _util.vars.at(k_var_trueX)<< "," << _util.vars.at(k_var_recoX) << ", " <<"w";
        
            for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
                evt_dist_sig << "," << "w_" << _uni;
            }
        
            evt_dist_sig << "\n";
            filled_sig = true;
        }
    
        evt_dist_sig << trueX << "," << recoX << "," << weight*_util.mc_scale_factor;
        
        for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
            evt_dist_sig << "," << _ev_weight.at(_uni)*_util.mc_scale_factor;
        }

        evt_dist_sig << "\n";
    }
    
    // Background event
    if ( _classification == "nu_out_fv"  || _classification == "cosmic"       ||
         _classification == "numu_cc"    || _classification == "numu_cc_pi0"  || _classification == "nc" || 
         _classification == "nc_pi0"     || ((_classification == "cosmic_nue" || _classification == "cosmic_nuebar") && _passed_selection)  ) {

        // Initialise the file 
        if (!filled_bkg){
            evt_dist_bkg << _util.vars.at(k_var_recoX) << "," << "w";
        
            for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
                evt_dist_bkg << "," << "w_" << _uni;
            }
        
            evt_dist_bkg << "\n";
            filled_bkg = true;
        }

        evt_dist_bkg << recoX << "," << weight*_util.mc_scale_factor;
        
        for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
            evt_dist_bkg << "," << _ev_weight.at(_uni)*_util.mc_scale_factor;
        }

        evt_dist_bkg << "\n";

    }
   
    // Generated event
    if (  _classification == "nue_cc"           || _classification == "nuebar_cc" || 
          _classification == "unmatched_nue"    || _classification == "cosmic_nue" ||
          _classification == "unmatched_nuebar" || _classification == "cosmic_nuebar") {
        
        isSignal = true;

        // Initialise the file 
        if (!filled_gen){
            evt_dist_gen << _util.vars.at(k_var_trueX) <<"," <<"w";
        
            for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
                evt_dist_gen << "," << "w_" << _uni;
            }
        
            evt_dist_gen << "\n";

            filled_gen = true;
        }

        evt_dist_gen << trueX << "," << weight*_util.mc_scale_factor;
        
        for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
            evt_dist_gen << "," << _ev_weight.at(_uni)*_util.mc_scale_factor;
        }

        evt_dist_gen << "\n";

    }
    
    // Data
    if (_classification == "data"){
        isData = true;
        isSelected = true;
        evt_dist_data << recoX << "\n";
    }

    // EXT background event
    if (_classification == "ext" ){
        isSelected = true;

        // Initialise the file 
        if (!filled_bkg){
            evt_dist_bkg << _util.vars.at(k_var_recoX) << "," << "w";
        
            for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
                evt_dist_bkg << "," << "w_" << _uni;
            }
        
            evt_dist_bkg << "\n";

            filled_bkg = true;
        }

        evt_dist_bkg << recoX << "," << weight*_util.ext_scale_factor;
        
        for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
            evt_dist_bkg << "," << _ev_weight.at(_uni)*_util.ext_scale_factor;
        }

        evt_dist_bkg << "\n";
    }
    // Dirt background event
    if (_classification == "dirt"){
        isSelected = true;

        // Initialise the file 
        if (!filled_bkg){
            evt_dist_bkg << _util.vars.at(k_var_recoX) << "," << "w";
        
            for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
                evt_dist_bkg << "," << "w_" << _uni;
            }
        
            evt_dist_bkg << "\n";

            filled_bkg = true;
        }

        evt_dist_bkg << recoX << "," << weight*_util.dirt_scale_factor;
        
        for (unsigned int _uni = 0; _uni <_ev_weight.size(); _uni++ ){
            evt_dist_bkg << "," << _ev_weight.at(_uni)*_util.dirt_scale_factor;
        }

        evt_dist_bkg << "\n";
    }

    xTrue = trueX;
    xReco = recoX;

    event_tree->Fill();

    // Reset TTree _util.vars
    isData = false;
    isSignal = false;
    isSelected = false;
    

}
// -----------------------------------------------------------------------------
int CrossSectionHelper::GetBinIndex(){

    
    int xbin = h_2D_CV_binning->GetXaxis()->FindBin(elec_e);
    int ybin = h_2D_CV_binning->GetYaxis()->FindBin(true_effective_angle);

    return xbin*ybin;
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoadDetvarCVHist(){

    std::cout << _util.red << "Loading in detector variation CV histogram to smear" << _util.reset << std::endl;
    
    // Load in the cross section file
    TFile* f = new TFile( Form("files/crosssec_run%s.root", _util.run_period) , "UPDATE");
    
    TH1D* h_temp;
    
    // Get the CV histogram we want to smear
    h_temp = (TH1D*)f->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_bkg", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
    h_detvar_cv_bkg = (TH1D*)h_temp->Clone();
    h_detvar_cv_bkg->SetDirectory(0);

    h_temp = (TH1D*)f->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_ext", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
    h_temp->Scale(_util.ext_scale_factor / _util.mc_scale_factor);
    h_detvar_cv_bkg->Add(h_temp);

    h_temp = (TH1D*)f->Get(Form("detvar_CV/%s/h_run%s_CV_0_%s_dirt", _util.vars.at(k_var_recoX).c_str(), _util.run_period, _util.vars.at(k_var_recoX).c_str()));
    h_temp->Scale(_util.dirt_scale_factor / _util.mc_scale_factor);
    h_detvar_cv_bkg->Add(h_temp);

    delete h_temp;

    f->Close();
}
// -----------------------------------------------------------------------------
void CrossSectionHelper::UnregularizedUnfold(TH1D *h_data_xsec_reco, TH1D* h_data_xsec_true, TH2D* h_response){

    // Getting covariance matrix in TMatrix form
    TMatrix res_inv;
    res_inv.Clear();
    res_inv.ResizeTo(h_response->GetNbinsX(), h_response->GetNbinsY());

    // loop over rows
    for (int i = 0; i < h_response->GetNbinsX(); i++) {

        // loop over columns
        for (int j = 0; j < h_response->GetNbinsY(); j++) {
            res_inv[i][j] = h_response->GetBinContent(i+1, j+1);
        }
    
    }

    // res_inv.Print();

    // Inverting the response matrix
    TMatrix inverse_res_m = res_inv.Invert();

    // inverse_res_m.Print();


    double nbins;

    // Electron/Shower Energy
    if (std::string(_util.xsec_var) =="elec_E"){
        nbins = _util.reco_shr_bins.size();
    }
    // Electron/Shower beta
    else if (std::string(_util.xsec_var) =="elec_ang"){
        nbins = _util.reco_shr_bins_ang.size();
    }
    // Electron/Shower cos(beta)
    else if (std::string(_util.xsec_var) =="elec_cang"){
        nbins = _util.reco_shr_bins_cang.size();
    }

    // Clear the Bins
    for (int bin = 0; bin < h_data_xsec_true->GetNbinsX()+2; bin++){
        h_data_xsec_true->SetBinContent(bin, 0);
    }

    // --- Do the matrix multiplication --- 
    // Loop over cols
    for (int i=0; i<inverse_res_m.GetNcols(); i++){

        // Now normalise the column entries by the number of events in the 1D generated histogram
        for (int j=0; j<inverse_res_m.GetNrows(); j++) {         
            
            h_data_xsec_true->SetBinContent(i+1, h_data_xsec_true->GetBinContent(i+1) + inverse_res_m[j][i] * (h_data_xsec_reco->GetBinContent(j+1)));
        }

        h_data_xsec_true->SetBinError(i+1, h_data_xsec_true->GetBinContent(i+1) * (h_data_xsec_reco->GetBinError(i+1) / h_data_xsec_reco->GetBinContent(i+1)) );

    } 


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::SaveGenXSec(){

    // To run :
    // ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecplot gen --var ../ntuples/genie_v2_intrinsic_events.root geniev2gen
    // ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecplot gen --var ../ntuples/nuwro_intrinsic_events.root nuwrogen

    TFile *f_gen = TFile::Open(_util.mc_file_name, "READ");


    TTree *t_gen = (TTree*)f_gen->Get("tree");

    double* edges = &_util.reco_shr_bins[0]; // Cast to an array 
    TH1D *htemp_e= new TH1D("h_elec_E","", _util.reco_shr_bins.size()-1, edges);

    edges = &_util.reco_shr_bins_cang[0]; // Cast to an array 
    TH1D *htemp_cang= new TH1D("h_elec_cang","", _util.reco_shr_bins_cang.size()-1, edges);

    TH1D *htemp_tot= new TH1D("h_elec_tot","", 1, 0, 1.1);

    t_gen->Draw("elec_e >> h_elec_E", "ppfx_cv*(elec_e > 0.12 && nu_e > 0.06 && infv && ccnc == 0)");
    t_gen->Draw("cosbeta >> h_elec_cang", "ppfx_cv*(elec_e > 0.12 && nu_e > 0.06 && infv && ccnc == 0)");
    t_gen->Draw("0.5 >> h_elec_tot", "ppfx_cv*(elec_e > 0.12 && nu_e > 0.06 && infv && ccnc == 0)");

    TFile *f_out;
    
    if (std::string(_util.xsec_bin_mode) == "ratio"){
        f_out = TFile::Open(Form("files/crosssec_run%s_ratio.root", _util.run_period), "UPDATE");
    }
    else {
        f_out = TFile::Open(Form("files/crosssec_run%s.root", _util.run_period), "UPDATE");
    }

    htemp_e->Scale(1.0 / (integrated_flux * mc_flux_scale_factor * N_target_MC));
    htemp_cang->Scale(1.0 / (integrated_flux * mc_flux_scale_factor * N_target_MC));
    htemp_tot->Scale(1.0 / (integrated_flux * mc_flux_scale_factor * N_target_MC));

    htemp_e->Scale(1.0e39);
    htemp_cang->Scale(1.0e39);
    htemp_tot->Scale(1.0e39);

    if (std::string(_util.xsec_bin_mode) == "ratio"){
        htemp_e->Scale(1.0/htemp_e->Integral());
        htemp_cang->Scale(1.0/htemp_cang->Integral());
    }

    htemp_e->SetOption("hist");
    htemp_cang->SetOption("hist");
    htemp_tot->SetOption("hist");
    
    // Create a directory to save the histograms
    TDirectory *dir;
    
    // Get the directory 
    bool bool_dir = _util.GetDirectory(f_out, dir, _util.variation );

    // Make the directory
    if (!bool_dir){
        dir = f_out->mkdir(_util.variation);
        _util.GetDirectory(f_out, dir, _util.variation );
    }

    dir->cd();

    htemp_e->Write("",TObject::kOverwrite);
    htemp_cang->Write("",TObject::kOverwrite);
    htemp_tot->Write("",TObject::kOverwrite);

    f_out->Close();

    f_gen->cd();
    f_gen->Close();


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::CheckPi0CrossSection(){


    TFile *f_genie = TFile::Open("../ntuples/neutrinoselection_filt_run1_overlay_newtune.root", "READ");
    TTree *t_genie = (TTree*)f_genie->Get("nuselection/NeutrinoSelectionFilter");

    TH1D *h_genie_cc = new TH1D("h_genie_cc","GENIE v3;#pi^{0} energy [GeV]; CC events", 20, 0, 0.6);
    TH1D *h_genie_nc = new TH1D("h_genie_nc","GENIE v3;#pi^{0} energy [GeV]; NC events", 20, 0, 0.6);
    TH1D *h_genie    = new TH1D("h_genie","GENIE v3;#pi^{0} energy [GeV]; CC + NC events", 20, 0, 0.6);

    t_genie->Draw("pi0_e >> h_genie_cc", "npi0 > 0 && ccnc == 0");
    t_genie->Draw("pi0_e >> h_genie_nc", "npi0 > 0 && ccnc == 1");
    t_genie->Draw("pi0_e >> h_genie",    "npi0 > 0");

    TFile *f_nuwro = TFile::Open("../ntuples/detvar_newtune/run1/neutrinoselection_filt_run1_overlay_nuwro.root", "READ");
    TTree *t_nuwro = (TTree*)f_nuwro->Get("nuselection/NeutrinoSelectionFilter");

    TH1D *h_nuwro_cc = new TH1D("h_nuwro_cc","NuWro;#pi^{0} energy [GeV]; CC events", 20, 0, 0.6);
    TH1D *h_nuwro_nc = new TH1D("h_nuwro_nc","NuWro;#pi^{0} energy [GeV]; NC events", 20, 0, 0.6);
    TH1D *h_nuwro    = new TH1D("h_nuwro","NuWro;#pi^{0} energy [GeV]; CC + NCevents", 20, 0, 0.6);

    t_nuwro->Draw("pi0_e >> h_nuwro_cc", "npi0 > 0 && ccnc == 0");
    t_nuwro->Draw("pi0_e >> h_nuwro_nc", "npi0 > 0 && ccnc == 1");
    t_nuwro->Draw("pi0_e >> h_nuwro",    "npi0 > 0");

    h_nuwro_cc->Scale(3.5165631);
    h_nuwro_nc->Scale(3.5165631);
    h_nuwro   ->Scale(3.5165631);

    TFile *f_out = TFile::Open("files/pi0_ratios.root", "RECREATE");

    h_genie_cc->SetOption("hist");
    h_genie_nc->SetOption("hist");
    h_genie->SetOption("hist");
    h_genie_cc->SetLineWidth(3);
    h_genie_nc->SetLineWidth(3);
    h_genie->SetLineWidth(3);
    h_genie_cc->SetLineColor(kRed+2);
    h_genie_nc->SetLineColor(kRed+2);
    h_genie->SetLineColor(kRed+2);

    h_genie_cc->Write();
    h_genie_nc->Write();
    h_genie->Write();

    h_genie_cc->Divide(h_nuwro_cc);
    h_genie_nc->Divide(h_nuwro_nc);
    h_genie->Divide(h_nuwro);

    h_genie_cc->SetTitle("GENIE v3/ NuWro CC #pi^{0};#pi^{0} energy [GeV]; Ratio");
    h_genie_nc->SetTitle("GENIE v3/ NuWro NC #pi^{0};#pi^{0} energy [GeV]; Ratio");
    h_genie->SetTitle("GENIE v3/ NuWro CC+NC #pi^{0};#pi^{0} energy [GeV]; Ratio");

    

    h_nuwro_cc->Write();
    h_nuwro_nc->Write();
    h_nuwro->Write();

    h_genie_cc->Write("h_ratio_cc");
    h_genie_nc->Write("h_ratio_nc");
    h_genie->Write("h_ratio");

    f_out->Close();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::StudyFluxShape(){

    // To Run:
    // ./nuexsec --run 1 --xsec files/trees/nuexsec_tree_merged_run1.root --xsecmode reweight --xseclabel ppfx    --xsec_smear er --xsecbins standard --xsecvar elec_E --xsecplot fluxshape

    gStyle->SetOptStat(0);

    TH1D* h_CV_ppfx = new TH1D("h_cv", "", 30, 0, 6);

    std::vector<TH1D*> h_ppfx_uni(600);

    for (unsigned int u = 0; u < h_ppfx_uni.size(); u++){
        h_ppfx_uni.at(u) = new TH1D("", "", 30, 0, 6);
    }


    std::cout << "Total Tree Entries: "<< tree->GetEntries() << std::endl;

    // Loop over the tree entries and weight the events in each universe
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){
        
        tree->GetEntry(ievent); 

        if (ievent % 20000 == 0 && ievent > 0) std::cout << "On entry " << ievent/10000.0 <<"0k " << std::endl;

        double cv_weight = weight; // SplinetimesTune * PPFX CV * Pi0 Tune
        double weight_dirt = weight; // Use this for estimating dirt and POT sys
        double weight_ext  = weight;  // Use this for estimating ext

        // Loop over the reweighter labels
        for (unsigned int label = 0; label < reweighter_labels.size(); label++){
            
            // Call switch function
            SwitchReweighterLabel(reweighter_labels.at(label));

            // Ensure this isnt ever zero
            if (vec_universes.size() == 0) {
                for (unsigned int t=0; t < h_cross_sec.at(label).size(); t++){
                    vec_universes.push_back(1.0);
                }
            }
            
            if (ievent == 0) std::cout << "Running over reweighter label: " << reweighter_labels.at(label) << " with " << vec_universes.size() << " universes" << std::endl;

            // Now loop over the universes
            for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){
                
                // Update the CV weight to CV * universe i
                double weight_uni{1.0}; 

                // Set the weight for universe i
                SetUniverseWeight(reweighter_labels.at(label), weight_uni, weight_dirt, weight_ext, weightSplineTimesTune, *classification, cv_weight, uni, nu_pdg, true_energy, numi_ang, npi0, pi0_e, ccnc);

                // Generated event
                if (  *classification == "nue_cc"           || *classification == "nuebar_cc" || 
                      *classification == "unmatched_nue"    || *classification == "cosmic_nue" ||
                      *classification == "unmatched_nuebar" || *classification == "cosmic_nuebar") {
                    
                    if (reweighter_labels.at(label) == "CV"){
                        h_CV_ppfx->Fill(true_energy, weight_uni);
                    }
                    else {
                        h_ppfx_uni.at(uni)->Fill(true_energy, weight_uni);
                    }                    
                    
                }
                
            } // End loop over uni

        } // End loop over labels
        
    } // End loop over events


    // Set the integrated flux
    double temp_integrated_flux = integrated_flux;
    h_CV_ppfx->Scale(1.0/(integrated_flux* mc_flux_scale_factor * N_target_MC));
    h_CV_ppfx->Scale(1.0e-39);

    // If we are reweighting by the PPFX Multisims, we need to change the integrated flux too
    for (int uni = 0; uni < 600; uni++){
        temp_integrated_flux = GetIntegratedFluxHP(uni, "ppfx_ms_UBPPFX");
        h_ppfx_uni.at(uni)->Scale(1.0/(temp_integrated_flux* mc_flux_scale_factor * N_target_MC));
        h_ppfx_uni.at(uni)->Scale(1.0e-39);
    }

    // Now make the covariance/correlation matrix
    int n_bins = h_CV_ppfx->GetNbinsX();
    TH2D* cov  = new TH2D("h_cov",   ";Bin i; Bin j",n_bins, 1, n_bins+1, n_bins, 1, n_bins+1 ); // Create the covariance matrix
    _util.CalcCovariance(h_ppfx_uni, h_CV_ppfx, cov);

    TH2D* cor = (TH2D*)cov->Clone();
    _util.CalcCorrelation(h_CV_ppfx, cov, cor);

    TCanvas *c = new TCanvas("c", "c", 500, 500);
    cor->GetXaxis()->CenterLabels(kTRUE);
    cor->GetYaxis()->CenterLabels(kTRUE);
    cor->GetZaxis()->CenterTitle();
    cor->GetZaxis()->SetTitleOffset(1.45);
    cor->SetMarkerColor(kRed+1);
    gStyle->SetPaintTextFormat("0.3f");
    cor->SetTitle("Correlation Matrix");
    cor->GetZaxis()->SetTitle("Correlation");
    if (std::string(_util.xsec_var) == "elec_cang"){
        cor->GetXaxis()->SetNdivisions(cor->GetNbinsX(), 0, 0, kFALSE);
        cor->GetYaxis()->SetNdivisions(cor->GetNbinsY(), 0, 0, kFALSE);
    }
    cor->Draw("colz");
    _util.IncreaseLabelSize(cor, c);
    cor->SetMarkerSize(1.3);
    _util.Draw_ubooneSim(c, 0.30, 0.915, 0.30, 0.915);
    cor->GetZaxis()->SetRangeUser(-1,1);
    c->Print("test_correlation.pdf");
    delete c;


}