#include "../include/CrossSectionHelper.h"

// -----------------------------------------------------------------------------
void CrossSectionHelper::Initialise(const char *_run_period, const char * xsec_file_in, Utility _utility, const char* run_mode){

    std::cout << "Initalising Cross Section Helper..." << std::endl;
    _util = _utility;

    // Set the run period
    run_period = std::string(_run_period);

    // Set the scale factors
    if (strcmp(_run_period, "1") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run1_Data_POT)  / _util.config_v.at(_util.k_Run1_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run1_Data_trig) / _util.config_v.at(_util.k_Run1_EXT_trig);
        mc_flux_scale_factor   = flux_scale_factor * _util.config_v.at(_util.k_Run1_MC_POT);
        data_flux_scale_factor = flux_scale_factor * _util.config_v.at(_util.k_Run1_Data_POT);
    }
    else if (strcmp(_run_period, "3") == 0){
        mc_scale_factor     = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_MC_POT);
        dirt_scale_factor   = _util.config_v.at(_util.k_Run3_Data_POT)  / _util.config_v.at(_util.k_Run3_Dirt_POT);
        intime_scale_factor = _util.config_v.at(_util.k_Run3_Data_trig) / _util.config_v.at(_util.k_Run3_EXT_trig);
        mc_flux_scale_factor   = flux_scale_factor * _util.config_v.at(_util.k_Run3_MC_POT); 
        data_flux_scale_factor = flux_scale_factor * _util.config_v.at(_util.k_Run3_Data_POT);
    }
    else {
        std::cout << "Error Krish... You havent defined the run3b POT numbers yet you donut!" << std::endl;
        exit(1);
    }
    
    std::cout << "\033[0;32m-------------------------------" << std::endl;
    std::cout << "Scale Factors:\n" <<
    "MC Scale factor:        " << mc_scale_factor          << "\n" <<
    "Dirt Scale factor:      " << dirt_scale_factor        << "\n" <<
    "EXT Scale factor:       " << intime_scale_factor      << "\n" <<
    "MC Flux Scale factor:   " << mc_flux_scale_factor     << "\n" <<
    "Data Flux Scale factor: " << data_flux_scale_factor
    << std::endl;
    std::cout << "-------------------------------\033[0m" << std::endl;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( xsec_file_in ) ) {
        f_nuexsec = new TFile( xsec_file_in, "READ");
    }

    // Get the Input TTree
    _util.GetTree(f_nuexsec, tree, "tree");
    if (tree == NULL) {
        std::cout << "Error failed to get the standard tree, maybe this is a MC only file, so trying to get the mc_tree"<< std::endl;
        _util.GetTree(f_nuexsec, tree, "mc_tree");
    }

    // Initialise the tree (set a bunch of tree branches)
    InitTree();

    // Get the integrated flux
    integrated_flux = GetIntegratedFlux();

    // Calculate the volume as used in the cuts
    volume = (_util.config_v.at(_util.k_config_x2) - _util.config_v.at(_util.k_config_x1)) * 
             (_util.config_v.at(_util.k_config_y2) - _util.config_v.at(_util.k_config_y1)) * 
             (_util.config_v.at(_util.k_config_z2) - _util.config_v.at(_util.k_config_z1));

    std::cout << "Volume used in cuts: " << volume << std::endl;

    N_target_MC   = (lar_density_mc   * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in MC: " << N_target_MC << std::endl;
    
    N_target_Data = (lar_density_data * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in Data: " << N_target_Data << std::endl;


    // Create and initialise vector of histograms
    InitialiseHistograms(std::string(run_mode));


    // Now loop over events and caluclate the cross section
    LoopEvents();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoopEvents(){

    // Loop over the tree entries
    std::cout << "Total Tree Entries: "<< tree->GetEntries() << std::endl;

    int n_gen = 0;

    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        if (shr_energy_cali > 5.0 && *classifcation == "data" ) std::cout << "reco shower energy was:  " << shr_energy_cali << "  Consider updating the bins" <<std::endl;

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

                // If we are using the genie systematics and unisim systematics then we want to undo the genie tune on them
                if (reweighter_labels.at(label) == "weightsReint" || reweighter_labels.at(label) == "weightsPPFX" || reweighter_labels.at(label) == "CV" ){
                    weight_uni = cv_weight * vec_universes.at(uni);
                }
                else {
                    // Note we actually dont want to divide out by the spline, but since this is 1 in numi, it doesnt matter!
                    // We do this because the interaction systematics are shifted about the genie tune as the CV

                    if (std::isnan(weightSplineTimesTune) == 1)   weightSplineTimesTune   = 1.0;
                    if (std::isinf(weightSplineTimesTune))        weightSplineTimesTune   = 1.0; 
                    if (weightSplineTimesTune == 0) weightSplineTimesTune = 1.0;
                    weight_uni = cv_weight * vec_universes.at(uni) / weightSplineTimesTune;
                }


                // Signal event
                if ((*classifcation == "nue_cc" || *classifcation == "nuebar_cc" || *classifcation == "unmatched_nue" || *classifcation == "unmatched_nuebar") && gen == false) {
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_sig, weight_uni, shr_energy_cali, elec_e);
                    FillHists(label, uni, k_xsec_sel, weight_uni, shr_energy_cali, elec_e);

                }

                // Background event
                if ( *classifcation == "nu_out_fv"  || *classifcation == "cosmic"      ||
                     *classifcation == "numu_cc"    || *classifcation == "numu_cc_pi0" || *classifcation == "nc" || 
                     *classifcation == "nc_pi0"     || *classifcation == "cosmic_nue" || *classifcation == "cosmic_nuebar"){
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_bkg, weight_uni, shr_energy_cali, elec_e);
                    FillHists(label, uni, k_xsec_sel, weight_uni, shr_energy_cali, elec_e);
                    
                }
                
                // Generated event
                if ( (*classifcation == "nue_cc"|| *classifcation == "nuebar_cc" || *classifcation == "unmatched_nue" || *classifcation == "cosmic_nue" || *classifcation == "unmatched_nuebar" || *classifcation == "cosmic_nuebar") && gen == true) {
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_gen, weight_uni, shr_energy_cali, elec_e);
                }

                // Data event
                if (*classifcation == "data"){

                    if (cv_weight != 1.0) std::cout << "Error weight for data is not 1, this means your weighting the data... bad!"<< std::endl;
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_data, cv_weight, shr_energy_cali, elec_e);
                }

                // Off beam event
                if (*classifcation == "ext"){

                    if (cv_weight != 1.0) std::cout << "Error weight for data is not 1, this means your weighting the data... bad!"<< std::endl;
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_ext, cv_weight, shr_energy_cali, elec_e);
                }

                // Dirt event
                if (*classifcation == "dirt"){
                    
                    // Fill histograms
                    FillHists(label, uni, k_xsec_dirt, cv_weight, shr_energy_cali, elec_e);
                }
            } // End loop over uni

        } // End loop over labels
        
    } // End loop over events

    std::cout << "Finished Event Loop, now calculating the cross-sections"<< std::endl;

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

                // Calculate the efficiency histogram
                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff)->Divide(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sig), h_cross_sec.at(label).at(uni).at(var).at(k_xsec_gen));

                // MC Cross section
                CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_sel), // N Sel
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),  // Eff
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),  // N Bkg
                                mc_scale_factor,
                                integrated_flux * mc_flux_scale_factor,                // Flux
                                intime_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),  // N EXT
                                dirt_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_dirt), // N Dirt
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_mcxsec),
                                N_target_MC, "MC");

                // Data Cross section
                CalcCrossSecHist(h_cross_sec.at(label).at(uni).at(var).at(k_xsec_data), // N Sel
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_eff),   // Eff
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_bkg),   // N Bkg
                                mc_scale_factor,
                                integrated_flux * data_flux_scale_factor,               // Flux
                                intime_scale_factor,
                                h_cross_sec.at(label).at(uni).at(var).at(k_xsec_ext),   // N EXT
                                dirt_scale_factor,
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
    "EXT MC:          " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_ext) ->Integral()* (intime_scale_factor / mc_scale_factor) << "\n" << 
    "Dirt MC:         " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dirt)->Integral()* (dirt_scale_factor / mc_scale_factor) << "\n\n" << 
    
    "Selected Data:   " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_data)->Integral() << "\n" << 
    "Signal Data:     " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_sig) ->Integral()* mc_scale_factor << "\n" << 
    "Background Data: " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_bkg) ->Integral()* mc_scale_factor << "\n" << 
    "Generated Data:  " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_gen) ->Integral()* mc_scale_factor << "\n" << 
    "EXT Data:        " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_ext) ->Integral()* intime_scale_factor  << "\n" << 
    "Dirt Data:       " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dirt)->Integral()* dirt_scale_factor << "\n"
    << std::endl;

    std::cout << 
    "MC XSEC: "   << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_mcxsec)  ->Integral() << " cm2" << "\n" << 
    "Data XSEC: " << h_cross_sec.at(0).at(0).at(k_var_integrated).at(k_xsec_dataxsec)->Integral() << " cm2      \n"
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
void CrossSectionHelper::CalcCrossSecHist(TH1D* h_sel, TH1D* h_eff, TH1D* h_bkg, double mc_scale_factor, double flux, double intime_scale_factor, TH1D* h_ext, double dirt_scale_factor ,TH1D* h_dirt, TH1D* h_xsec, double targ, std::string mcdata){


    // I think this is the slow bit -- maybe make copies only once?
    TH1D *h_bkg_clone  = (TH1D*)h_bkg ->Clone("h_bkg_temp");
    TH1D *h_ext_clone  = (TH1D*)h_ext ->Clone("h_ext_temp");
    TH1D *h_dirt_clone = (TH1D*)h_dirt->Clone("h_dirt_temp");

    // Scale the relavent histograms to the MC/Data POT/Triggers
    if (mcdata == "MC"){
        h_ext_clone ->Scale(intime_scale_factor / mc_scale_factor);
        h_dirt_clone->Scale(dirt_scale_factor / mc_scale_factor);
    }
    else if (mcdata == "Data"){
        h_bkg_clone ->Scale(mc_scale_factor);
        h_ext_clone ->Scale(intime_scale_factor);
        h_dirt_clone->Scale(dirt_scale_factor);
    }
    else{
        std::cout << "error wrong mode entering xsec calculation" << std::endl;
    }

    // (S - B) / (eff * targ * flux)
    h_xsec->Add(h_sel,         1);
    h_xsec->Add(h_bkg_clone,  -1);
    h_xsec->Add(h_ext_clone,  -1);
    h_xsec->Add(h_dirt_clone, -1);
    
    h_xsec->Divide(h_eff) ;
    
    h_xsec->Scale(1.0 / (targ*flux) );

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFlux(){

    TFile * f_flux;

    TH1D* h_nue, *h_nuebar;

    double energy_threshold = 0.00; // Set threshold to integrate the flux from [GeV]
    std::cout << "Using a flux energy threshold of :" << energy_threshold * 1000 << " MeV"<< std::endl;

    // Hardcoded for run 1 right now.. urgh krish you lazy
    bool boolfile;
    
    std::string flux_file_name;

    if (run_period == "1"){
        flux_file_name = "Systematics/output_fhc_uboone_run0.root";
        boolfile = _util.GetFile(f_flux, flux_file_name);
    }
    else {
        flux_file_name = "Systematics/output_rhc_uboone_run0.root";
        boolfile = _util.GetFile(f_flux, flux_file_name );
    }
    std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;

    if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV

    bool boolhist = _util.GetHist(f_flux, h_nue, "nue/Detsmear/nue_CV_AV_TPC_5MeV_bin");     if (boolhist == false) gSystem->Exit(0); // Get the PPFX weighted CV for nue
    boolhist      = _util.GetHist(f_flux, h_nuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin"); if (boolhist == false) gSystem->Exit(0); // Get the PPFX weighted CV for nuebar

    double integral_nue{0};
    double integral_nuebar{0};

    double xbin_th = h_nue->GetXaxis()->FindBin(energy_threshold); // find the x bin to integrate from -- remove the MuDAR Peak where the cross sec is zero
    integral_nue += h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1); // Integrate over the flux for nue

    std::cout << "Integral Nue Flux: " << integral_nue << std::endl;

    xbin_th   = h_nuebar->GetXaxis()->FindBin(energy_threshold); // find the x bin to integrate from -- remove the MuDAR Peak where the cross sec is zero
    integral_nuebar += h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1); // Integrate over the flux for nue

    std::cout << "Integral Nuebar Flux: " << integral_nuebar << std::endl;


    double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
    POT_flux = GetPOT(f_flux);

    f_flux->Close();

    // Return the flux per POT
    return (integral_nue + integral_nuebar) / POT_flux;

}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetFluxUni(int uni, std::string value, std::string label){

    TFile * f_flux;

    TH1D* h_nue, *h_nuebar;
    TH1D* h_uni;
    double weight = {1.0};
    double xbin, ybin;

    double energy_threshold = 0.00; // Set threshold to integrate the flux from [GeV]
    std::cout << "Using a flux energy threshold of :" << energy_threshold * 1000 << " MeV"<< std::endl;

    bool boolfile;
    
    std::string flux_file_name;

    if (run_period == "1"){
        flux_file_name = "Systematics/output_fhc_uboone_run0.root";
        boolfile = _util.GetFile(f_flux, flux_file_name);
    }
    else {
        flux_file_name = "Systematics/output_rhc_uboone_run0.root";
        boolfile = _util.GetFile(f_flux, flux_file_name );
    }
    std::cout << "Using Flux file name: \033[0;31m" << flux_file_name << "\033[0m" <<  std::endl;

    if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV

    bool boolhist;
    
    // Will return the integrated nue+nuebar flux for universe i
    if (value == "flux"){

        boolhist = _util.GetHist(f_flux, h_nue, "nue/Detsmear/nue_CV_AV_TPC_5MeV_bin");
        if (boolhist == false) gSystem->Exit(0);
        
        boolhist      = _util.GetHist(f_flux, h_nuebar, "nuebar/Detsmear/nuebar_CV_AV_TPC_5MeV_bin");
        if (boolhist == false) gSystem->Exit(0);

        double integral_nue{0}, integral_nuebar{0};

        double xbin_th = h_nue->GetXaxis()->FindBin(energy_threshold); // find the x bin to integrate from -- remove the MuDAR Peak where the cross sec is zero
        integral_nue += h_nue->Integral( xbin_th, h_nue->GetNbinsX()+1); // Integrate over the flux for nue


        xbin_th   = h_nuebar->GetXaxis()->FindBin(energy_threshold); // find the x bin to integrate from -- remove the MuDAR Peak where the cross sec is zero
        integral_nuebar += h_nuebar->Integral( xbin_th, h_nuebar->GetNbinsX()+1); // Integrate over the flux for nue

        // std::cout << "Integral Flux: " << integral_nue << std::endl;

        double POT_flux{0.0}; // The POT of the flux file (i.e the POT used in the flux histogram)
        POT_flux = GetPOT(f_flux);

        f_flux->Close();

        // Return the flux per POT
        return (integral_nue + integral_nuebar) / POT_flux;

    }
    
    // Return a weight based on the energy and angle of the event -- used for weighting by beamline variations and ppfx HP types broken down 
    else {

        // Nue
        if (nu_pdg == 12){
            boolhist = _util.GetHist(f_flux, h_uni, Form("nue/Multisims/nue_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        }
        // Nuebar
        else if (nu_pdg == -12){
            boolhist = _util.GetHist(f_flux, h_uni, Form("nuebar/Multisims/nuebar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        }
        // Numu
        else if (nu_pdg == 14){
            boolhist = _util.GetHist(f_flux, h_uni, Form("numu/Multisims/numu_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        }
        // NumuBar
        else if (nu_pdg == -14){
            boolhist = _util.GetHist(f_flux, h_uni, Form("numubar/Multisims/numubar_%s_Uni_%i_AV_TPC_2D", label.c_str(), uni));
        }
        
        if (boolhist == false) gSystem->Exit(0);

        xbin = h_uni->GetXaxis()->FindBin(true_energy);
        ybin = h_uni->GetYaxis()->FindBin(numi_ang);
        weight = h_uni->GetBinContent(xbin, ybin);

        return weight;

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
    if (!gROOT->GetListOfFiles()->FindObject( Form("files/crosssec_run%s.root", run_period.c_str()) ) ) {
        fnuexsec_out = new TFile( Form("files/crosssec_run%s.root", run_period.c_str()) , "UPDATE");
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
    tree->SetBranchAddress("classifcation",   &classifcation);
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
    // This can be the CV
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
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,";Reco Electron Shower Energy [GeV]; Entries", nbins, edges);
                    }
                    else if (i == k_xsec_eff){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,";Reco Electron Shower Energy [GeV]; Efficiency", nbins, edges);
                    }
                    else if (i == k_xsec_dataxsec || i == k_xsec_mcxsec){
                        h_cross_sec.at(label).at(uni).at(var).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni, vars.at(var).c_str(), xsec_types.at(i).c_str()) ,Form("%s", var_labels.at(var).c_str()), nbins, edges);
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