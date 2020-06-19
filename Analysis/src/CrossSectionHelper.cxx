#include "../include/CrossSectionHelper.h"

// -----------------------------------------------------------------------------
void CrossSectionHelper::Initialise(const char *_run_period, const char * xsec_file_in, utility _utility){

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

    // Set the tree branches
    tree->SetBranchAddress("run",    &run);
    tree->SetBranchAddress("subrun", &subrun);
    tree->SetBranchAddress("event",  &event);
    tree->SetBranchAddress("gen",    &gen);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("reco_energy", &reco_energy);
    tree->SetBranchAddress("classifcation",   &classifcation);

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

    // Histogram stuff -- may want to put this in another function
    // Resize the histogram vector. plot var, cuts, classifications
    h_cross_sec.resize(k_TH1D_xsec_MAX);

    // loop over and create the histograms
    for (unsigned int i=0; i < xsec_types.size();i++){    
        h_cross_sec.at(i) = new TH1D ( Form("h_run%s_%s",_run_period, xsec_types.at(i).c_str()) ,"", nbins, edges);
    }
    // Now loop over events and caluclate the cross section
    LoopEvents();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::LoopEvents(){

    double n_data{0.0}, n_dirt{0.0}, n_ext{0.0}, n_selected{0.0}, n_gen{0.0}, n_bkg{0.0}, n_sig{0.0}, n_sel{0.0};

    // Loop over the tree entries
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        if (reco_energy > 6.0 && *classifcation == "data" ) std::cout << "reco energy was:  " << reco_energy << "  Consider updating the bins" <<std::endl;

        // Signal event
        if ((*classifcation == "nue_cc" || *classifcation == "nue_cc_mixed") && gen == false) {
            n_sel+=weight;
            n_sig+=weight;
            h_cross_sec.at(k_xsec_sel)->Fill(reco_energy, weight);
            h_cross_sec.at(k_xsec_sig)->Fill(reco_energy, weight);
        }

        // Background event
        if ( *classifcation == "nu_out_fv"  || *classifcation == "cosmic" ||
           *classifcation == "numu_cc" || *classifcation == "numu_cc_pi0" || *classifcation == "nc" || 
           *classifcation == "nc_pi0" || *classifcation == "unmatched"){
            n_bkg+=weight;
            n_sel+=weight;
            h_cross_sec.at(k_xsec_bkg)->Fill(reco_energy, weight);
            h_cross_sec.at(k_xsec_sel)->Fill(reco_energy, weight);
        }
        
        // Generated event
        if ( (*classifcation == "nue_cc"|| *classifcation == "nue_cc_mixed" ) && gen == true) {
            n_gen+=weight;
            h_cross_sec.at(k_xsec_gen)->Fill(reco_energy, weight);
        }

        // Data event
        if (*classifcation == "data"){
            n_data+=weight;
            h_cross_sec.at(k_xsec_data)->Fill(reco_energy, weight);
        }

        // Off beam event
        if (*classifcation == "ext"){
            n_ext+=weight;
            h_cross_sec.at(k_xsec_ext)->Fill(reco_energy, weight);
        }

        // Dirt event
        if (*classifcation == "dirt"){
            n_dirt+=weight;
            h_cross_sec.at(k_xsec_dirt)->Fill(reco_energy, weight);
        }

    }

    // MC 
    double n_dirt_mc   = n_dirt * (dirt_scale_factor / mc_scale_factor);
    double n_ext_mc    = n_ext  * (intime_scale_factor / mc_scale_factor);
    
    double n_sig_data = n_sig * mc_scale_factor;
    double n_bkg_data = n_bkg * mc_scale_factor;
    double n_gen_data = n_gen * mc_scale_factor;
    double n_dirt_data = n_dirt * dirt_scale_factor;
    double n_ext_data  = n_ext  * intime_scale_factor;

    std::cout << "\n" <<
    "Selected MC:     " << n_sel << "\n" << 
    "Signal MC:       " << n_sig << "\n" << 
    "Background MC:   " << n_bkg << "\n" << 
    "Generated MC:    " << n_gen << "\n" << 
    "EXT MC:          " << n_ext_mc << "\n" << 
    "Dirt MC:         " << n_dirt_mc << "\n\n" << 
    
    "Selected Data:   " << n_data << "\n" << 
    "Signal Data:     " << n_sig_data << "\n" << 
    "Background Data: " << n_bkg_data << "\n" << 
    "Generated Data:  " << n_gen_data << "\n" << 
    "EXT Data:        " << n_ext_data  << "\n" << 
    "Dirt Data:       " << n_dirt_data << "\n"
    << std::endl;

    // MC Cross section
    double mc_xsec   = CalcCrossSec(n_sel, n_gen, n_sig, n_bkg, integrated_flux * mc_flux_scale_factor , n_ext_mc, n_dirt_mc, N_target_MC);

    double data_xsec = CalcCrossSec(n_data, n_gen_data, n_sig_data, n_bkg_data, integrated_flux * data_flux_scale_factor , n_ext_data, n_dirt_data, N_target_Data);

    std::cout << 
    "MC XSEC: "   << mc_xsec << "\n" << 
    "Data XSEC: " << data_xsec
    << std::endl;


    // Now loop over the bins in the histogram and calculate the cross section
    for (int k =1; k < h_cross_sec.at(0)->GetNbinsX()+1; k++){

        double temp_xsec_mc = CalcCrossSec(h_cross_sec.at(k_xsec_sel)->GetBinContent(k),
                                           h_cross_sec.at(k_xsec_gen)->GetBinContent(k),
                                           h_cross_sec.at(k_xsec_sig)->GetBinContent(k),
                                           h_cross_sec.at(k_xsec_bkg)->GetBinContent(k), 
                                           integrated_flux * mc_flux_scale_factor,
                                           h_cross_sec.at(k_xsec_ext)->GetBinContent(k) * (intime_scale_factor / mc_scale_factor),
                                           h_cross_sec.at(k_xsec_dirt)->GetBinContent(k)* (dirt_scale_factor / mc_scale_factor),
                                           N_target_MC);


        double temp_xsec_data = CalcCrossSec(h_cross_sec.at(k_xsec_data)->GetBinContent(k),
                                             h_cross_sec.at(k_xsec_gen)->GetBinContent(k)* mc_scale_factor,
                                             h_cross_sec.at(k_xsec_sig)->GetBinContent(k)* mc_scale_factor,
                                             h_cross_sec.at(k_xsec_bkg)->GetBinContent(k)* mc_scale_factor, 
                                             integrated_flux * data_flux_scale_factor,
                                             h_cross_sec.at(k_xsec_ext)->GetBinContent(k) * intime_scale_factor,
                                             h_cross_sec.at(k_xsec_dirt)->GetBinContent(k)* dirt_scale_factor*0.45,
                                             N_target_Data);

        if (std::isnan(temp_xsec_mc) == 1) temp_xsec_mc = 0.0;
        if (std::isnan(temp_xsec_data) == 1) temp_xsec_data = 0.0;

        h_cross_sec.at(k_xsec_mcxsec)->SetBinContent(k, temp_xsec_mc);
        h_cross_sec.at(k_xsec_dataxsec)->SetBinContent(k, temp_xsec_data);

        std::cout << "Bin: " << k << "  MC XSec: " << temp_xsec_mc << std::endl;
        std::cout << "Bin: " << k << "  Data XSec: " << temp_xsec_data << std::endl;
    }

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

    // Now open the output file
    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( Form("files/crosssec_run%s.root", run_period.c_str()) ) ) {
        fnuexsec_out = new TFile( Form("files/crosssec_run%s.root", run_period.c_str()) , "UPDATE");
    }

    fnuexsec_out->cd();

    TDirectory *dir;
    bool bool_dir = _util.GetDirectory(fnuexsec_out, dir ,"xsec");
    if (bool_dir) dir->cd();
    
    for (unsigned int p = 0; p < h_cross_sec.size(); p++){
        h_cross_sec.at(p)->Write("",TObject::kOverwrite);
    }

    fnuexsec_out->Close();
    
}
// -----------------------------------------------------------------------------
