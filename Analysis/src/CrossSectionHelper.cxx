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
        std::cout << "Using Flux file name: \033[0;31m" << "Systematics/output_fhc_uboone_run0.root" << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f_flux, "Systematics/output_fhc_uboone_run0.root");
        if (boolfile == false) gSystem->Exit(0); 
    }
    else if (std::string(_util.run_period) == "3") {
        std::cout << "Using Flux file name: \033[0;31m" << "Systematics/output_rhc_uboone_run0.root" << "\033[0m" <<  std::endl;
        bool boolfile = _util.GetFile(f_flux, "Systematics/output_rhc_uboone_run0.root" );
        if (boolfile == false) gSystem->Exit(0); 
    }
    else{
        std::cout << "Unknown Run period configured, exiting..." << std::endl;
        return;
    }
    
    // Get the integrated flux for the CV
    integrated_flux = GetIntegratedFlux(0, "CV", "", "");


    // Now lets open the beamline variation files
    GetBeamlineHists();

    f_nuexsec->cd();

    // Calculate the volume as used in the cuts
    volume = (_util.config_v.at(_util.k_config_x2) - _util.config_v.at(_util.k_config_x1)) * 
             (_util.config_v.at(_util.k_config_y2) - _util.config_v.at(_util.k_config_y1)) * 
             (_util.config_v.at(_util.k_config_z2) - _util.config_v.at(_util.k_config_z1));

    std::cout << "Volume used in cuts: " << volume << std::endl;

    N_target_MC   = (lar_density_mc   * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in MC: " << N_target_MC << std::endl;
    
    N_target_Data = (lar_density_data * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in Data: " << N_target_Data << std::endl;
    std::cout << "  "<< std::endl;


    // Create and initialise vector of histograms
    InitialiseHistograms(std::string(_util.xsecmode));


    // Now loop over events and caluclate the cross section
    LoopEvents();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoopEvents(){

    
    std::cout << "Total Tree Entries: "<< tree->GetEntries() << std::endl;

    int n_gen = 0;

    // Loop over the tree entries and weight the events in each universe
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        if (shr_energy_cali > 4.0 && *classification == "data" ) std::cout << "reco shower energy was:  " << shr_energy_cali << "  Consider updating the bins" <<std::endl;

        double cv_weight = weight;

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
                if (std::isnan(vec_universes.at(uni)) == 1)   vec_universes.at(uni)   = 1.0;
                if (std::isinf(vec_universes.at(uni)))        vec_universes.at(uni)   = 1.0;
                double weight_uni{1.0}; 

                // Weight equal to universe weight times cv weight
                if (reweighter_labels.at(label) == "weightsReint" || reweighter_labels.at(label) == "weightsPPFX" || reweighter_labels.at(label) == "CV" ){
                    
                    _util.CheckWeight(vec_universes.at(uni));

                    weight_uni = cv_weight * vec_universes.at(uni);
                }
                // This is a beamline variation
                else if (CheckBeamline(reweighter_labels.at(label))){
                    weight_uni = cv_weight * GetIntegratedFlux(0, "", "", reweighter_labels.at(label));
                    // std::cout << GetIntegratedFlux(0, "", "", reweighter_labels.at(label)) << std::endl;
                }
                // If we are using the genie systematics and unisim systematics then we want to undo the genie tune on them
                else {
                    // Note we actually dont want to divide out by the spline, but since this is 1 in numi, it doesnt matter!
                    // We do this because the interaction systematics are shifted about the genie tune as the CV

                    // Check the spline times tune weight
                    _util.CheckWeight(weightSplineTimesTune);

                    // Check the uiverse weight
                    if (std::isnan(vec_universes.at(uni)) == 1 || std::isinf(vec_universes.at(uni)) || vec_universes.at(uni) < 0 || vec_universes.at(uni) > 30) {
                        
                        // We set the universe to be the spline times tune, so it cancels with the divide below to just return the CV weight
                        // i.e. a universe weight of 1
                        vec_universes.at(uni)   = weightSplineTimesTune;
                    }

                    if (weightSplineTimesTune == 0) weight_uni = 0.0; // Special case where the genie tune or we have thresholded events are zero
                    else weight_uni = (cv_weight * vec_universes.at(uni)) / weightSplineTimesTune;

                    // std::cout << vec_universes.at(uni) << " " << weight_uni << "   "<< weightSplineTimesTune<< std::endl;

                }


                // Signal event
                if ((*classification == "nue_cc" || *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "unmatched_nuebar") && gen == false) {
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_sig, weight_uni, shr_energy_cali, elec_e);
                    FillHists(label, uni, k_xsec_sel, cv_weight, shr_energy_cali, elec_e);  // Selected events (N term) we dont weight

                }

                // Background event
                if ( *classification == "nu_out_fv"  || *classification == "cosmic"       ||
                     *classification == "numu_cc"    || *classification == "numu_cc_pi0"  || *classification == "nc" || 
                     *classification == "nc_pi0"     || ((*classification == "cosmic_nue" || *classification == "cosmic_nuebar") && !gen)  ) {
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_bkg, weight_uni, shr_energy_cali, elec_e);
                    FillHists(label, uni, k_xsec_sel, cv_weight, shr_energy_cali, elec_e);  // Selected events (N term) we dont weight
                    
                }

                // Generated event
                if ( (*classification == "nue_cc"|| *classification == "nuebar_cc" || *classification == "unmatched_nue" || *classification == "cosmic_nue" || *classification == "unmatched_nuebar" || *classification == "cosmic_nuebar") && gen == true) {
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_gen, weight_uni, shr_energy_cali, elec_e);
                }
                
                // Data event
                if (*classification == "data"){

                    if (cv_weight != 1.0) std::cout << "Error weight for data is not 1, this means your weighting the data... bad!"<< std::endl;
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_data, cv_weight, shr_energy_cali, elec_e);
                }

                // Off beam event
                if (*classification == "ext"){
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_ext, cv_weight, shr_energy_cali, elec_e);
                }

                // Dirt event
                if (*classification == "dirt"){
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_dirt, cv_weight, shr_energy_cali, elec_e);
                }
            } // End loop over uni

        } // End loop over labels
        
    } // End loop over events

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
                if (reweighter_labels.at(label) == "weightsPPFX") temp_integrated_flux = GetIntegratedFlux(uni, "HP", "ppfx_ms_UBPPFX", reweighter_labels.at(label));

                // If this is a beamline variation then we use the corresponding beamline flux
                if (CheckBeamline(reweighter_labels.at(label))) {
                    int var_index = GetBeamlineIndex(reweighter_labels.at(label));
                    temp_integrated_flux = beamline_flux.at(var_index);
                }

                // Calculate the efficiency histogram -- this is incorrect, we need to smear the truth efficiency!
                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff)->Divide(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig), h_cross_sec.at(label).at(uni).at(var).at(k_xsec_gen));

                // MC Cross section -- currently using eventrate
                CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sel), // N Sel
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),  // Eff
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),  // N Bkg
                                _util.mc_scale_factor,
                                temp_integrated_flux * mc_flux_scale_factor,           // Flux
                                _util.ext_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),  // N EXT
                                _util.dirt_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dirt), // N Dirt
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec),
                                N_target_MC, "MC");

                // Data Cross section -- currently using eventrate
                CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_data), // N Sel
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),   // Eff
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),   // N Bkg
                                _util.mc_scale_factor,
                                temp_integrated_flux * data_flux_scale_factor,          // Flux
                                _util.ext_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),   // N EXT
                                _util.dirt_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dirt),  // N Dirt
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dataxsec),
                                N_target_Data, "Data");

            } // End loop over the vars
        
        } // End loop over universes
    
    } // End loop over labels

    // Print the CV Results for the Flux Integrated Measurement
    // Label 0 should always be the CV with 1 universe (which counts from 0)

    std::cout << "\n" <<
    "Selected MC:     " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sel) ->Integral() << "\n" << 
    "Signal MC:       " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sig) ->Integral() << "\n" << 
    "Background MC:   " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_bkg) ->Integral() << "\n" << 
    "Generated MC:    " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_gen) ->Integral() << "\n" << 
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
    "MC Flux Norm. Event Rate: "   << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_mcxsec)  ->Integral() << " cm2" << "\n" << 
    "Data Flux Norm. Event Rate: " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dataxsec)->Integral() << " cm2      \n"
    << std::endl;


    // Write the histograms to file for inspection
    WriteHists();


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
void CrossSectionHelper::CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double ext_scale_factor, TH1D* h_ext, double dirt_scale_factor ,TH1D* h_dirt, TH1D* h_xsec, double targ, std::string mcdata){


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
    h_xsec->Add(h_sel,         1);
    h_xsec->Add(h_bkg_clone,  -1);
    h_xsec->Add(h_ext_clone,  -1);
    h_xsec->Add(h_dirt_clone, -1);
    
    // h_xsec->Divide(h_eff) ; // For flux normalised event rate we dont do anything to the efficiency

    h_xsec->Scale(1.0 / (targ*flux) );

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFlux(int uni, std::string value, std::string label, std::string variation){

    f_flux->cd();

    TH2D* h_nue, *h_nuebar;
    TH2D *h_uni_nue, *h_uni_nuebar;
    double weight = {1.0};
    double xbin, ybin;
    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    double xbin_th{0.0};
    double integral_nue{0.0}, integral_nuebar{0.0};
    bool boolfile, boolhist;
    
    // Get the CV flux from the 5MeV binned histograms
    if (value == "CV"){

        // Get the nue flux histogram from the file
        boolhist = _util.GetHist(f_flux, h_nue,         "nue/Detsmear/nue_CV_AV_TPC_2D");       
        if (boolhist == false) gSystem->Exit(0); 
        
        // Get the nuebar flux histogram from the file
        boolhist      = _util.GetHist(f_flux, h_nuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_2D");
        if (boolhist == false) gSystem->Exit(0); 

        std::cout << "Using Energy Threshold of: \033[0;31m" << energy_threshold * 1000 << " MeV" << "\033[0m" <<  std::endl;

        POT_flux = GetPOT(f_flux);

        double xbin_th = h_nue->GetXaxis()->FindBin(energy_threshold);   // find the x bin to integrate from
        integral_nue = h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1, 0, h_nue->GetNbinsY()+1); // Integrate over the flux for nue

        std::cout << "\nIntegral Nue Flux: " << flux_scale_factor * integral_nue / POT_flux << " nue / POT / GeV / cm2" << std::endl;

        xbin_th   = h_nuebar->GetXaxis()->FindBin(energy_threshold);              // find the x bin to integrate from
        integral_nuebar = h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1, 0, h_nuebar->GetNbinsY()+1); // Integrate over the flux for nue

        std::cout << "\nIntegral Nuebar Flux: " << flux_scale_factor * integral_nuebar / POT_flux << " nuebar / POT / GeV / cm2" << "\n" << std::endl;

        // Return the flux per POT
        return (integral_nue + integral_nuebar) / POT_flux;

    }
    // Will return the integrated nue+nuebar flux for universe i
    else if (value == "HP"){

        // Get the POT from the file
        POT_flux = GetPOT(f_flux);

        boolhist = _util.GetHist(f_flux, h_uni_nue, Form("nue/Multisims/nue_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));

        if (boolhist == false) gSystem->Exit(0);
        
        boolhist = _util.GetHist(f_flux, h_uni_nuebar, Form("nuebar/Multisims/nuebar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        if (boolhist == false) gSystem->Exit(0);

        xbin_th = h_uni_nue->GetXaxis()->FindBin(energy_threshold);   // find the x bin to integrate from
        integral_nue = h_uni_nue->Integral( xbin_th, h_uni_nue->GetNbinsX()+1, 0, h_uni_nue->GetNbinsY()+1); // Integrate over the flux for nue

        // std::cout << xbin_th << std::endl;

        // std::cout << "\nIntegral Nue Flux: " << flux_scale_factor * integral_nue / POT_flux << " nue / POT / GeV / cm2" << std::endl;

        xbin_th   = h_uni_nuebar->GetXaxis()->FindBin(energy_threshold);              // find the x bin to integrate from 
        integral_nuebar = h_uni_nuebar->Integral( xbin_th, h_uni_nuebar->GetNbinsX()+1, 0, h_uni_nuebar->GetNbinsY()+1); // Integrate over the flux for nuebar

        // std::cout << "Integral Nuebar Flux: " << flux_scale_factor * integral_nuebar / POT_flux << " nuebar / POT / GeV / cm2" << "\n" << std::endl;

        // Return the flux per POT
        return (integral_nue + integral_nuebar) / POT_flux;

    }
    
    // Return a weight based on the energy and angle of the event -- used for weighting by beamline variations and ppfx HP types broken down 
    else {

        int var_index = GetBeamlineIndex(variation);

        // std::cout <<true_energy << "  " <<  numi_ang << "  " << nu_pdg<< std::endl;

        // Nue
        if (nu_pdg == 12){
            xbin = beamline_hists.at(var_index).at(k_nue)->GetXaxis()->FindBin(true_energy);
            ybin = beamline_hists.at(var_index).at(k_nue)->GetYaxis()->FindBin(numi_ang);
            return ( beamline_hists.at(var_index).at(k_nue)->GetBinContent(xbin, ybin));
        }
        // Nuebar
        else if (nu_pdg == -12){
            xbin = beamline_hists.at(var_index).at(k_nuebar)->GetXaxis()->FindBin(true_energy);
            ybin = beamline_hists.at(var_index).at(k_nuebar)->GetYaxis()->FindBin(numi_ang);
            return ( beamline_hists.at(var_index).at(k_nuebar)->GetBinContent(xbin, ybin));
        }
        // Numu
        else if (nu_pdg == 14){
            xbin = beamline_hists.at(var_index).at(k_numu)->GetXaxis()->FindBin(true_energy);
            ybin = beamline_hists.at(var_index).at(k_numu)->GetYaxis()->FindBin(numi_ang);
            return ( beamline_hists.at(var_index).at(k_numu)->GetBinContent(xbin, ybin));
        }
        // NumuBar
        else if (nu_pdg == -14){
            xbin = beamline_hists.at(var_index).at(k_numubar)->GetXaxis()->FindBin(true_energy);
            ybin = beamline_hists.at(var_index).at(k_numubar)->GetYaxis()->FindBin(numi_ang);
            return ( beamline_hists.at(var_index).at(k_numubar)->GetBinContent(xbin, ybin));
        }
        else return 1.0;
    }

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
        fnuexsec_out = new TFile( Form("files/crosssec_run%s.root", _util.run_period) , "UPDATE");
    }

    fnuexsec_out->cd();

    // Create subdirectory for each reweighter
    TDirectory *dir_labels[reweighter_labels.size()];

    // Create subdirectory for each variable
    TDirectory *dir_labels_var[vars.size()];

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
                bool bool_dir = _util.GetDirectory(fnuexsec_out, dir_labels_var[var], Form("%s/%s", reweighter_labels.at(label).c_str(), vars.at(var).c_str()));

                // If it doesnt exist then create it
                if (!bool_dir) dir_labels_var[var] = dir_labels[label]->mkdir(vars.at(var).c_str());

                // Go into the directory
                dir_labels_var[var]->cd();
            
                // Now write the histograms, 
                for (unsigned int p = 0; p < h_cross_sec.at(label).at(uni).at(var).size(); p++){
                    h_cross_sec.at(label).at(uni).at(var).at(p)->SetOption("hist");
                    h_cross_sec.at(label).at(uni).at(var).at(p)->Write("",TObject::kOverwrite);
                }
            
            } // End loop over the variables

            fnuexsec_out->cd();    // change current directory to top
        
        } // End loop over universes

        fnuexsec_out->cd();    // change current directory to top

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

    // Set the bins here
    bins.resize(vars.size());

    // Integrated X-Section Bin definition
    bins.at(k_var_integrated) = { 0.0, 1.1 };

    // Reconstructed electron energy Bin definition
    bins.at(k_var_reco_el_E) = _util.reco_shr_bins;

    // True electron energy Bin definition
    bins.at(k_var_true_el_E) = _util.reco_shr_bins;

    // True neutrino energy Bin definition
    // bins.at(k_var_true_nu_E) = { 0.0, 0.25, 0.56, 0.89, 1.61, 3.37, 5.0};

    // Reconstructed neutrino energy Bin definition
    // bins.at(k_var_reco_nu_E) = { 0.0, 0.25, 0.56, 0.89, 1.61, 3.37, 5.0};

    // Resize to the number of reweighters
    h_cross_sec.resize(reweighter_labels.size());

    // Resize each reweighter to their number of universes
    for (unsigned int j=0; j < reweighter_labels.size(); j++){

        // Specific resizing -- hardcoded and may break in the future
        if ( reweighter_labels.at(j) == "weightsGenie"){
            std::cout << "Setting Genie All Histogram universe vector to size: " << uni_genie << std::endl;
            h_cross_sec.at(j).resize(uni_genie);
        }
        // Specific resizing -- hardcoded and may break in the future
        else if ( reweighter_labels.at(j) == "weightsPPFX"){
            std::cout << "Setting PPFX All Histogram universe vector to size: " << uni_ppfx << std::endl;
            h_cross_sec.at(j).resize(uni_ppfx);
        }
        // Specific resizing -- hardcoded and may break in the future
        else if ( reweighter_labels.at(j) == "weightsReint" ){
            std::cout << "Setting Geant Reinteractions Histogram universe vector to size: " << uni_reint << std::endl;
            h_cross_sec.at(j).resize(uni_reint);
        }
        // Default size of 1
        else {
            h_cross_sec.at(j).resize(1);
        }

    }

    // Now resize the universe to each type of variable we want to plot/reweight
    
    // Resize each reweighter to their number of universes
    for (unsigned int label=0; label < reweighter_labels.size() ;label++){

        for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){
            // Resize the histogram vector. plot var, cuts, classifications
            h_cross_sec.at(label).at(uni).resize(vars.size());
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


    // Loop over the rewighters
    for (unsigned int label=0; label < reweighter_labels.size(); label++){

        // loop over the universes
        for (unsigned int uni=0; uni < h_cross_sec.at(label).size(); uni++){
            
            // Loop over the variables
            for (unsigned int var = 0; var < h_cross_sec.at(label).at(uni).size(); var ++){

                // Get the number of bins and the right vector
                int const nbins = bins.at(var).size()-1;
                temp_bins.clear();
                temp_bins = bins.at(var);
                double* edges = &temp_bins[0]; // Cast to an array 


                // loop over and create the histograms
                for (unsigned int i=0; i < xsec_types.size();i++){    
                    if (i == k_xsec_sel || i == k_xsec_bkg || i == k_xsec_gen || i == k_xsec_sig || i == k_xsec_ext || i == k_xsec_dirt || i == k_xsec_data){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,";Reco. Leading Shower Energy [GeV]; Entries", nbins, edges);
                    }
                    else if (i == k_xsec_eff){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,";Reco. Leading Shower Energy [GeV]; Efficiency", nbins, edges);
                    }
                    else if (i == k_xsec_dataxsec || i == k_xsec_mcxsec){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",_util.run_period, reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,Form("%s", var_labels.at(var).c_str()), nbins, edges);
                    }
                    else {
                        // If this is the case then there is an uncaught case
                        std::cout << "Houston we have a problem!" << std::endl;
                    }
                }
            } // End loop over the variables
        
        } // End loop over the universes
    
    } // End loop over the labels

    std::cout << "Initialisation of cross-section histograms is complete!" << std::endl;


}
// -----------------------------------------------------------------------------
void CrossSectionHelper::FillHists(int label, int uni, int xsec_type, double weight_uni, float shr_energy_cali, float elec_e){

    // Integrated
    h_cross_sec.at(label).at(uni).at(k_var_integrated).at(xsec_type)->Fill(1.0, weight_uni);

    // Reco Electron Energy
    h_cross_sec.at(label).at(uni).at(k_var_reco_el_E).at(xsec_type)->Fill(shr_energy_cali, weight_uni);

    // True Electron Energy
    h_cross_sec.at(label).at(uni).at(k_var_true_el_E).at(xsec_type)->Fill(elec_e, weight_uni);

    // True Neutrino Energy
    // h_cross_sec.at(label).at(uni).at(k_var_true_nu_E).at(xsec_type)->Fill(true_energy, weight_uni);

    // Reconstructed Neutrino Energy
    // h_cross_sec.at(label).at(uni).at(k_var_reco_nu_E).at(xsec_type)->Fill(reco_energy, weight_uni);

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::GetBeamlineHists(){

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
            f_temp = TFile::Open(Form("Systematics/beamline/FHC/output_uboone_fhc_run%i.root", beamline_map.at(f).second));
        }
        else if(std::string(_util.run_period) == "3") {
            f_temp = TFile::Open(Form("Systematics/beamline/RHC/output_uboone_rhc_run%i.root", beamline_map.at(f).second));
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
        double xbin_th   = beamline_hists.at(f).at(k_nue)->GetXaxis()->FindBin(energy_threshold);              // find the x bin to integrate from 
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
        xbin_th   = beamline_hists.at(f).at(k_nuebar)->GetXaxis()->FindBin(energy_threshold);              // find the x bin to integrate from 
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
        variation == "Horm1_x_m3mm"        ||
        variation == "Horn1_y_p3mm"        ||
        variation == "Horn1_y_m3mm"        ||
        variation == "Beam_spot_1_1mm"     ||
        variation == "Beam_spot_1_5mm"     ||
        variation == "Horn2_x_p3mm"        ||
        variation == "Horm2_x_m3mm"        ||
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
    else if  (variation == "Horm1_x_m3mm"       ) return k_Horm1_x_m3mm;
    else if  (variation == "Horn1_y_p3mm"       ) return k_Horn1_y_p3mm;
    else if  (variation == "Horn1_y_m3mm"       ) return k_Horn1_y_m3mm;
    else if  (variation == "Beam_spot_1_1mm"    ) return k_Beam_spot_1_1mm;
    else if  (variation == "Beam_spot_1_5mm"    ) return k_Beam_spot_1_5mm;
    else if  (variation == "Horn2_x_p3mm"       ) return k_Horn2_x_p3mm;
    else if  (variation == "Horm2_x_m3mm"       ) return k_Horm2_x_m3mm;
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