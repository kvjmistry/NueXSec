#include "../include/CrossSectionHelper.h"

// -----------------------------------------------------------------------------
void CrossSectionHelper::Initialise(const char *_run_period, const char * xsec_file_in, utility _utility, const char* run_mode){

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

    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        if (shr_energy_tot_cali > 5.0 && *classifcation == "data" ) std::cout << "reco shower energy was:  " << shr_energy_tot_cali << "  Consider updating the bins" <<std::endl;

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
                double weight_uni = cv_weight * vec_universes.at(uni);


                // Signal event
                if ((*classifcation == "nue_cc" || *classifcation == "nuebar_cc") && gen == false) {
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_sig)->Fill(shr_energy_tot_cali, weight_uni);
                    h_cross_sec.at(label).at(uni).at(k_xsec_sel)->Fill(shr_energy_tot_cali, cv_weight); // Selected events (N term) we dont weight
                }

                // Background event
                if ( *classifcation == "nu_out_fv"  || *classifcation == "cosmic"      ||
                     *classifcation == "numu_cc"    || *classifcation == "numu_cc_pi0" || *classifcation == "nc" || 
                     *classifcation == "nc_pi0"     || *classifcation == "unmatched"){
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_bkg)->Fill(shr_energy_tot_cali, weight_uni);
                    h_cross_sec.at(label).at(uni).at(k_xsec_sel)->Fill(shr_energy_tot_cali, cv_weight); // Selected events (N term) we dont weight
                }
                
                // Generated event
                if ( (*classifcation == "nue_cc"|| *classifcation == "nuebar_cc" ) && gen == true) {
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_gen)->Fill(shr_energy_tot_cali, weight_uni);
                }

                // Data event
                if (*classifcation == "data"){

                    if (cv_weight != 1.0) std::cout << "Error weight for data is not 1, this means your weighting the data... bad!"<< std::endl;
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_data)->Fill(shr_energy_tot_cali, cv_weight);
                }

                // Off beam event
                if (*classifcation == "ext"){
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_ext)->Fill(shr_energy_tot_cali, cv_weight);
                }

                // Dirt event
                if (*classifcation == "dirt"){
                    
                    // Fill histograms
                    h_cross_sec.at(label).at(uni).at(k_xsec_dirt)->Fill(shr_energy_tot_cali, cv_weight);
                }
            } // End loop over uni

        } // End loop over labels
        
    } // End loop over events

    std::cout << "Finished Event Loop"<< std::endl;

    // ----

    // Now we have rewieghted the evnts, we want to calculate the cross-sections for each label for each universe
    
    // Loop over the reweighter labels
    for (unsigned int label = 0; label < reweighter_labels.size(); label++){
        
        // Now loop over the universes
        for (unsigned int uni = 0; uni < h_cross_sec.at(label).size(); uni++){
    
            // loop over the bins in the histogram and calculate the cross section
            for (int bin =1; bin < h_cross_sec.at(label).at(uni).at(k_xsec_sel)->GetNbinsX()+1; bin++){

                double temp_xsec_mc = CalcCrossSec(h_cross_sec.at(label).at(uni).at(k_xsec_sel)->GetBinContent(bin), // N Sel
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_gen)->GetBinContent(bin), // N Gen
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_sig)->GetBinContent(bin), // N Sig
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_bkg)->GetBinContent(bin), // N Bkg
                                                   integrated_flux * mc_flux_scale_factor,                           // Flux
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_ext)->GetBinContent(bin) * (intime_scale_factor / mc_scale_factor), // N EXT
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_dirt)->GetBinContent(bin)* (dirt_scale_factor / mc_scale_factor),   // N Dirt
                                                   N_target_MC);


                double temp_xsec_data = CalcCrossSec(h_cross_sec.at(label).at(uni).at(k_xsec_data)->GetBinContent(bin),                 // N Sel
                                                     h_cross_sec.at(label).at(uni).at(k_xsec_gen)->GetBinContent(bin)* mc_scale_factor, // N Gen
                                                     h_cross_sec.at(label).at(uni).at(k_xsec_sig)->GetBinContent(bin)* mc_scale_factor, // N Sig
                                                     h_cross_sec.at(label).at(uni).at(k_xsec_bkg)->GetBinContent(bin)* mc_scale_factor, // N Bkg
                                                     integrated_flux * data_flux_scale_factor,                                          // Flux
                                                     h_cross_sec.at(label).at(uni).at(k_xsec_ext)->GetBinContent(bin) * intime_scale_factor, // N EXT
                                                     h_cross_sec.at(label).at(uni).at(k_xsec_dirt)->GetBinContent(bin)* dirt_scale_factor, // N Dirt
                                                     N_target_Data);


                // Validate the cross sec, only accept if its > 0...
                if (std::isnan(temp_xsec_mc) == 1)   temp_xsec_mc   = 0.0;
                if (std::isnan(temp_xsec_data) == 1) temp_xsec_data = 0.0;
                if (std::isinf(temp_xsec_mc))        temp_xsec_mc   = 0.0; 
                if (std::isinf(temp_xsec_data))      temp_xsec_data = 0.0; 


                if (temp_xsec_mc > 0)   h_cross_sec.at(label).at(uni).at(k_xsec_mcxsec)  ->SetBinContent(bin, temp_xsec_mc/(1e-39) );
                if (temp_xsec_data > 0) h_cross_sec.at(label).at(uni).at(k_xsec_dataxsec)->SetBinContent(bin, temp_xsec_data/(1e-39));

                // std::cout << "Bin: " << bin << "  MC XSec: " << temp_xsec_mc << std::endl;
                // std::cout << "Bin: " << bin << "  Data XSec: " << temp_xsec_data << std::endl;
            
            } // End loop over bins

            // Calculate the Flux integrated Cross section
            double temp_xsec_mc_int = CalcCrossSec(h_cross_sec.at(label).at(uni).at(k_xsec_sel)->Integral(), // N Sel
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_gen)->Integral(), // N Gen
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_sig)->Integral(), // N Sig
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_bkg)->Integral(), // N Bkg
                                                   integrated_flux * mc_flux_scale_factor,                           // Flux
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_ext)->Integral()  * (intime_scale_factor / mc_scale_factor), // N EXT
                                                   h_cross_sec.at(label).at(uni).at(k_xsec_dirt)->Integral() * (dirt_scale_factor / mc_scale_factor),   // N Dirt
                                                   N_target_MC);


            double temp_xsec_data_int = CalcCrossSec(h_cross_sec.at(label).at(uni).at(k_xsec_data)->Integral(),             // N Sel
                                                    h_cross_sec.at(label).at(uni).at(k_xsec_gen)->Integral()* mc_scale_factor, // N Gen
                                                    h_cross_sec.at(label).at(uni).at(k_xsec_sig)->Integral()* mc_scale_factor, // N Sig
                                                    h_cross_sec.at(label).at(uni).at(k_xsec_bkg)->Integral()* mc_scale_factor, // N Bkg
                                                    integrated_flux * data_flux_scale_factor,                                       // Flux
                                                    h_cross_sec.at(label).at(uni).at(k_xsec_ext)->Integral() * intime_scale_factor, // N EXT
                                                    h_cross_sec.at(label).at(uni).at(k_xsec_dirt)->Integral()* dirt_scale_factor,   // N Dirt
                                                    N_target_Data);

            // Validate the cross sec, only accept if its > 0...
            if (std::isnan(temp_xsec_mc_int) == 1)   temp_xsec_mc_int   = 0.0;
            if (std::isnan(temp_xsec_data_int) == 1) temp_xsec_data_int = 0.0;
            if (std::isinf(temp_xsec_mc_int))        temp_xsec_mc_int   = 0.0; 
            if (std::isinf(temp_xsec_data_int))      temp_xsec_data_int = 0.0; 


            if (temp_xsec_mc_int > 0)   h_cross_sec.at(label).at(uni).at(k_xsec_mcxsec_int)->SetBinContent(1, temp_xsec_mc_int/(1e-39) );
            if (temp_xsec_data_int > 0) h_cross_sec.at(label).at(uni).at(k_xsec_dataxsec_int)->SetBinContent(1, temp_xsec_data_int/(1e-39));
        
        } // End loop over universes
    
    } // End loop over labels

    // Print the CV Results for the Flux Integrated Measurement
    // Label 0 should always be the CV with 1 universe

    std::cout << "\n" <<
    "Selected MC:     " << h_cross_sec.at(0).at(0).at(k_xsec_sel)->Integral() << "\n" << 
    "Signal MC:       " << h_cross_sec.at(0).at(0).at(k_xsec_sig)->Integral() << "\n" << 
    "Background MC:   " << h_cross_sec.at(0).at(0).at(k_xsec_bkg)->Integral() << "\n" << 
    "Generated MC:    " << h_cross_sec.at(0).at(0).at(k_xsec_gen)->Integral() << "\n" << 
    "EXT MC:          " << h_cross_sec.at(0).at(0).at(k_xsec_ext)->Integral()* (intime_scale_factor / mc_scale_factor) << "\n" << 
    "Dirt MC:         " << h_cross_sec.at(0).at(0).at(k_xsec_dirt)->Integral()* (dirt_scale_factor / mc_scale_factor) << "\n\n" << 
    
    "Selected Data:   " << h_cross_sec.at(0).at(0).at(k_xsec_data)->Integral() << "\n" << 
    "Signal Data:     " << h_cross_sec.at(0).at(0).at(k_xsec_sig)->Integral()* mc_scale_factor << "\n" << 
    "Background Data: " << h_cross_sec.at(0).at(0).at(k_xsec_bkg)->Integral()* mc_scale_factor << "\n" << 
    "Generated Data:  " << h_cross_sec.at(0).at(0).at(k_xsec_gen)->Integral()* mc_scale_factor << "\n" << 
    "EXT Data:        " << h_cross_sec.at(0).at(0).at(k_xsec_ext)->Integral()* intime_scale_factor  << "\n" << 
    "Dirt Data:       " << h_cross_sec.at(0).at(0).at(k_xsec_dirt)->Integral()* dirt_scale_factor << "\n"
    << std::endl;

    std::cout << 
    "MC XSEC: "   << h_cross_sec.at(0).at(0).at(k_xsec_mcxsec_int)->Integral() << "e-39 cm2" << "\n" << 
    "Data XSEC: " << h_cross_sec.at(0).at(0).at(k_xsec_dataxsec_int)->Integral() << "e-39 cm2\n"
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
double CrossSectionHelper::GetIntegratedFlux(){

    TFile * f_flux;

    TH1D* h_nue, *h_nuebar;

    double energy_threshold = 0.05; // Set threshold to integrate the flux from [GeV]
    std::cout << "Using a flux energy threshold of :" << energy_threshold * 1000 << " MeV"<< std::endl;

    // Hardcoded for run 1 right now.. urgh krish you lazy
    bool boolfile;
    
    std::string flux_file_name;

    if (run_period == "1"){
        flux_file_name = "../Systematics/output_fhc_uboone_run0.root";
        boolfile = _util.GetFile(f_flux, flux_file_name);
    }
    else {
        flux_file_name = "../Systematics/output_rhc_uboone_run0.root";
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
double CrossSectionHelper::GetPOT(TFile* f){
    TTree* TPOT = (TTree*) f->Get("POT");
    if (TPOT == NULL) std::cout << "Error cant get POT info" << std::endl;

    double fPOT{0};
    TPOT->SetBranchAddress("POT", &fPOT); // Get the POT
    TPOT->GetEntry(0);
    double total_entries = TPOT->GetEntries(); // if using hadd, this will not be 1 equal to 1 anymore
    fPOT*=total_entries;
    std::cout << "TOTAL POT READ IN:\t" << fPOT << std::endl;

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
            
            // Now write the histograms, 
            for (unsigned int p = 0; p < h_cross_sec.at(label).at(uni).size(); p++){
                h_cross_sec.at(label).at(uni).at(p)->SetOption("hist");
                h_cross_sec.at(label).at(uni).at(p)->Write("",TObject::kOverwrite);
            }
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
    tree->SetBranchAddress("reco_energy", &reco_energy);
    tree->SetBranchAddress("classifcation",   &classifcation);
    tree->SetBranchAddress("shr_energy_tot_cali", &shr_energy_tot_cali);
    tree->SetBranchAddress("elec_e",  &elec_e);
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

    // Now resize the universe to each type of histogram we want to plot/reweight
    // Resize each reweighter to their number of universes
    for (unsigned int j=0; j < reweighter_labels.size(); j++){

        for (unsigned int y=0; y < h_cross_sec.at(j).size(); y++){
            // Resize the histogram vector. plot var, cuts, classifications
            h_cross_sec.at(j).at(y).resize(k_TH1D_xsec_MAX);
        }

    }

    // Loop over the rewighters
    for (unsigned int label=0; label < reweighter_labels.size(); label++){

        // loop over the universes
        for (unsigned int uni=0; uni < h_cross_sec.at(label).size(); uni++){
            
            // loop over and create the histograms
            for (unsigned int i=0; i < xsec_types.size();i++){    
                if (i == k_xsec_sel || i == k_xsec_bkg || i == k_xsec_gen || i == k_xsec_sig || i == k_xsec_ext || i == k_xsec_dirt || i == k_xsec_data){
                    h_cross_sec.at(label).at(uni).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni , xsec_types.at(i).c_str()) ,";Reco Electron Shower Energy [GeV]; Entries", nbins, edges);
                }
                else if (i == k_xsec_dataxsec || i == k_xsec_mcxsec){
                    h_cross_sec.at(label).at(uni).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni , xsec_types.at(i).c_str()) ,";Reco Electron Shower Energy [GeV]; #frac{d#sigma_{#nu_{e} + #bar{#nu}_{e}}}{dE^{reco}_{e}} CC Cross-Section [10^{-39} cm^{2}]", nbins, edges);
                }
                else {
                    // Flux Integrated Histograms
                    h_cross_sec.at(label).at(uni).at(i) = new TH1D ( Form("h_run%s_%s_%i_%s",run_period.c_str(), reweighter_labels.at(label).c_str(), uni , xsec_types.at(i).c_str()) ,";; #nu_{e} + #bar{#nu}_{e} CC Cross-Section [10^{-39} cm^{2}]", 1, 0, 1);
                }
            }
        }
    }

    std::cout << "Initialisation of cross-section histograms is complete!" << std::endl;


}