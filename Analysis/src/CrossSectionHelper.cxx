#include "../include/CrossSectionHelper.h"

// -----------------------------------------------------------------------------
void CrossSectionHelper::Initialise(const char *_run_period, const char * xsec_file_in){

    std::cout << "Initalising Cross Section Helper..." << std::endl;

    // Set the run period
    run_period = _run_period;

    // File not already open, open the file
    if (!gROOT->GetListOfFiles()->FindObject( xsec_file_in ) ) {
        f_nuexsec = new TFile( xsec_file_in, "READ");
    }

    // Get the Input TTree
    _util.GetTree(f_nuexsec, tree, "tree");

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

    N_target_MC   = (lar_density_mc   * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in MC: " << N_target_MC << std::endl;
    
    N_target_Data = (lar_density_data * volume * NA * N_nuc) / m_mol;
    std::cout << "Number of Target Nucleons in Data: " << N_target_Data << std::endl;

    // Now calculate the cross section
    CalcCrossSec();

}
// -----------------------------------------------------------------------------
void CrossSectionHelper::CalcCrossSec(){

    double n_data{0.0}, n_dirt{0.0}, n_ext{0.0}, n_selected{0.0}, n_gen{0.0}, n_bkg{0.0}, n_sig{0.0}, n_sel{0.0};

    // Loop over the tree entries
    for (unsigned int ievent = 0; ievent < tree->GetEntries(); ievent++){

        tree->GetEntry(ievent); 

        // Signal event
        if (*classifcation == "nue_cc" && gen == false) {
            n_sel+=weight;
            n_sig+=weight;
        }

        // Background event
        if (*classifcation == "nue_cc_mixed" || *classifcation == "nu_out_fv"  || *classifcation == "cosmic" ||
           *classifcation == "numu_cc" || *classifcation == "numu_cc_pi0" || *classifcation == "nc" || 
           *classifcation == "nc_pi0" || *classifcation == "unmatched"){
            n_bkg+=weight;
            n_sel+=weight;
        }
        
        // Generated event
        if (*classifcation == "nue_cc" && gen == true) {
            n_gen+=weight;
        }

        // Data event
        if (*classifcation == "data"){
            n_data+=weight;
        }

        // Off beam event
        if (*classifcation == "ext"){
            n_ext+=weight;
        }

        // Dirt event
        if (*classifcation == "dirt"){
            n_dirt+=weight;
        }

    }

    std::cout << 
    "Selected:   " << n_sel << "\n" << 
    "Signal:     " << n_sig << "\n" << 
    "Background: " << n_bkg << "\n" << 
    "Generated:  " << n_gen << "\n" << 
    "EXT:        " << n_ext << "\n" << 
    "Dirt:       " << n_dirt << "\n" << 
    "Data:       " << n_data << "\n" 
    << std::endl;

    
}
// -----------------------------------------------------------------------------
double CrossSectionHelper::GetIntegratedFlux(){

    TFile * f_flux;

    TH1D* h_nue, *h_nuebar;

    double energy_threshold = 0.05; // Set threshold to integrate the flux from [GeV]

    // Hardcoded for run 1 right now.. urgh krish you lazy
    bool boolfile  = _util.GetFile(f_flux , "/uboone/data/users/kmistry/work/PPFX/uboone/beamline_zero_threshold/output_uboone_run0.root"); if (boolfile == false) gSystem->Exit(0); // Most up to date version of CV

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

    return integral_nue + integral_nuebar;

}
// -----------------------------------------------------------------------------