#ifndef TREEHELPER_h
#define TREEHELPER_h

#include "SliceContainer.h"

// Class for filling and saving ttree variables to a file
// This will be useful for producing a flat tree for weighting events
class TreeHelper{

    public:
    // Default constructor
    TreeHelper(){};

    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(int type, const char *run_period, const char * file_out, int weight_cfg );
    // -------------------------------------------------------------------------
    // Function to fill the tree vars
    void FillVars(SliceContainer &SC, std::pair<std::string, int> _classification, bool _gen);
    // -------------------------------------------------------------------------
    // Writes the tree to file
    void WriteTree(int type);
    // -------------------------------------------------------------------------
    // Fill the counter tree
    void Fill_counters(std::vector<double> counter_v, std::string cut_name, bool bool_use_mc, bool bool_use_ext, bool bool_use_data, bool bool_use_dirt);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    
    // Destructor 
    // ~TreeHelper(); 

    // The output file
    TFile* f_nuexsec;

    // Class instances
    utility _util;

    int _type{1};

    bool weight_tune = true; // Apply genie tune weight
    bool weight_ppfx = true; // Apply ppfx cv weight

    // Tree variables
    int run{0}, subrun{0}, event{0};
    std::string classifcation; // The classification of the event
    
    // Is the event a true signal event in the FV that was not selected?
    // We still need these for the efficiency
    bool gen{false};           
    
    double weight{0.0};        // This is not going to be integer if we already weight the CV

    TTree * tree;       // Main tree with the selected events
    TTree * counter_tree; // Tree for storing the selection results


    std::string cut;

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
    
    double count_nue_cc{0.0};
    double count_nue_cc_mixed{0.0};
    double count_nu_out_fv{0.0};
    double count_cosmic{0.0};
    double count_numu_cc{0.0};
    double count_numu_cc_pi0{0.0};
    double count_nc{0.0};
    double count_nc_pi0{0.0};
    double count_unmatched{0.0};
    double count_total_mc{0.0};
    double count_data{0.0};
    double count_ext{0.0};
    double count_dirt{0.0};

    

    private:

    // Here we create the trees 



}; // End Class TreeHelper

#endif