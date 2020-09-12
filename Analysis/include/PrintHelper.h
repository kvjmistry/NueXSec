#ifndef PRINTHELPER_H
#define PRINTHELPER_H

#include "Utility.h"

// Class for printing the selection results
class PrintHelper{

    public:
    // Default constructor
    // PrintHelper(){};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(Utility _utility );
    // -------------------------------------------------------------------------
    // Function to print the selection
    void PrintResults();
    // -------------------------------------------------------------------------
    // Get the histogram files we need to get the efficiency and purity uncertainties
    void GetHists();
    // -------------------------------------------------------------------------
    
    // Destructor 
    // ~PrintHelper(); 

    // The output file
    TFile* f_mc, *f_data, *f_ext, *f_dirt;
    TFile* f_mc_hist, *f_data_hist, *f_ext_hist, *f_dirt_hist; // histogram files

    // Class instances
    Utility _util;

   
    TTree * mc_counter_tree; // MC Counter Tree
    TTree * data_counter_tree; // Data Counter Tree
    TTree * ext_counter_tree; // EXT Counter Tree
    TTree * dirt_counter_tree; // Dirt Counter Tree

    TTree * eff_tree;   // Efficiency and Purity tree

    // Histograms for getting the efficiency along with its uncertainty
    std::vector<std::vector<TH1D*>> TEfficiency_hists; // nue/nuebar -- cut index
    std::vector<std::vector<TH1D*>> TPurity_hists; // classification -- cut index
    std::vector<std::vector<TH1D*>> TPurity_hists_tot; // cut index -- numerator/denominator

    // For storing the errors on the efficiency
    std::vector<std::vector<double>> vec_err;

    // Define vectors for selected and generated events
    std::vector<std::vector<double>> vec_n;
    std::vector<std::vector<double>> vec_N;


    int tree_total_entries{0}; // Should equal number of cuts

    // enum for plots by efficiency, we only care about the single bin efficiencies
    enum TH1D_eff_vars {
        k_eff_nu_E_single_bin,       // True Electron-neutrino energy, single bin
        k_eff_nu_E_nue_single_bin,   // True Electron-neutrino energy single bin
        k_eff_nu_E_nuebar_single_bin,// True anti Electron-neutrino energy single bin
        k_TH1D_eff_MAX
    };

    // Scale factors (everything is scaled to data)
    double mc_scale_factor     = 1.0;
    double ext_scale_factor    = 1.0;
    double dirt_scale_factor   = 1.0;

    double efficiency{0.0}, purity{0.0};

    double tot_true_infv_nues{1.0};

    std::string *cutdata = NULL;
    std::string *cut = NULL;

    // Counters
    double count_nue_cc_qe{0.0};
    double count_nue_cc_res{0.0};
    double count_nue_cc_dis{0.0};
    double count_nue_cc_coh{0.0};
    double count_nue_cc_mec{0.0};
    
    double count_nuebar_cc_qe{0.0};
    double count_nuebar_cc_res{0.0};
    double count_nuebar_cc_dis{0.0};
    double count_nuebar_cc_coh{0.0};
    double count_nuebar_cc_mec{0.0};
    
    double count_nue_cc_infv{0.0};
    double count_nuebar_cc_infv{0.0};
    double count_nue_cc_incryo{0.0};
    double count_nuebar_cc_incryo{0.0};
    
    double count_numu_cc_qe{0.0};
    double count_numu_cc_res{0.0};
    double count_numu_cc_dis{0.0};
    double count_numu_cc_coh{0.0};
    double count_numu_cc_mec{0.0};
    
    double count_numubar_cc_qe{0.0};
    double count_numubar_cc_res{0.0};
    double count_numubar_cc_dis{0.0};
    double count_numubar_cc_coh{0.0};
    double count_numubar_cc_mec{0.0};
    
    double count_numu_cc_infv{0.0};
    double count_numubar_cc_infv{0.0};
    double count_numu_cc_incryo{0.0};
    double count_numubar_cc_incryo{0.0};
    
    double count_pi0_nue_cc_nopi0{0.0};
    double count_pi0_nue_cc_pi0{0.0};
    double count_pi0_nuebar_cc_nopi0{0.0};
    double count_pi0_nuebar_cc_pi0{0.0};
    double count_pi0_numu_cc_nopi0{0.0};
    double count_pi0_numu_cc_pi0{0.0};
    double count_pi0_nc_nopi0{0.0};
    double count_pi0_nc_pi0{0.0};

    double count_nue_cc{0.0};
    double count_nuebar_cc{0.0};
    double count_nu_out_fv{0.0};
    double count_cosmic{0.0};
    double count_numu_cc{0.0};
    double count_numu_cc_pi0{0.0};
    double count_nc{0.0};
    double count_nc_pi0{0.0};
    double count_unmatched{0.0};
    double count_unmatched_nue{0.0};
    double count_cosmic_nue{0.0};
    double count_unmatched_nuebar{0.0};
    double count_cosmic_nuebar{0.0};
    double count_total_mc{0.0};
    double count_data{0.0};
    double count_ext{0.0};
    double count_dirt{0.0};

    // Counters at no cuts applied
    double init_count_nue_cc{0.0};
    double init_count_nuebar_cc{0.0};
    double init_count_nu_out_fv{0.0};
    double init_count_cosmic{0.0};
    double init_count_numu_cc{0.0};
    double init_count_numu_cc_pi0{0.0};
    double init_count_nc{0.0};
    double init_count_nc_pi0{0.0};
    double init_count_unmatched{0.0};
    double init_count_unmatched_nue{0.0};
    double init_count_cosmic_nue{0.0};
    double init_count_unmatched_nuebar{0.0};
    double init_count_cosmic_nuebar{0.0};
    double init_count_ext{0.0};
    double init_count_dirt{0.0};

    // Counters for previous cut
    double prev_count_nue_cc{1.0};
    double prev_count_nuebar_cc{1.0};
    double prev_count_nu_out_fv{1.0};
    double prev_count_cosmic{1.0};
    double prev_count_numu_cc{1.0};
    double prev_count_numu_cc_pi0{1.0};
    double prev_count_nc{1.0};
    double prev_count_nc_pi0{1.0};
    double prev_count_unmatched{1.0};
    double prev_count_unmatched_nue{1.0};
    double prev_count_cosmic_nue{1.0};
    double prev_count_unmatched_nuebar{1.0};
    double prev_count_cosmic_nuebar{1.0};
    double prev_count_ext{1.0};
    double prev_count_dirt{1.0};

    // The efficiency and purity from the previous cut
    double efficiency_last{0.0};
    double purity_last{0.0};
    

    private:

    // Here we create the trees 



}; // End Class PrintHelper

#endif