#ifndef PRINTHELPER_H
#define PRINTHELPER_H

#include "Utility.h"

// Class for printing the selection results
class PrintHelper{

    public:
    // Default constructor
    PrintHelper(){};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char* run_period, const char * mc_file_in, bool _print_mc, bool _print_data, bool _print_ext, bool _print_dirt, Utility _utility );
    // -------------------------------------------------------------------------
    // Function to print the selection
    void PrintResults();
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    
    // Destructor 
    // ~PrintHelper(); 

    // The output file
    TFile* f_mc, *f_data, *f_ext, *f_dirt;

    // Class instances
    Utility _util;

   
    TTree * mc_counter_tree; // MC Counter Tree
    TTree * data_counter_tree; // Data Counter Tree
    TTree * ext_counter_tree; // EXT Counter Tree
    TTree * dirt_counter_tree; // Dirt Counter Tree

    TTree * eff_tree;   // Efficiency and Purity tree

    bool print_mc;
    bool print_data;
    bool print_ext;
    bool print_dirt;

    int tree_total_entries{0}; // Should equal number of cuts

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