#ifndef SYSTEMATICSHELPER_h
#define SYSTEMATICSHELPER_h

#include "utility.h"

// Class for making plots for systematic studies
class SystematicsHelper{

    public:
    // Default constructor
    SystematicsHelper(){};
    
    // The output file
    std::vector<TFile*> f_vars;

    // Class instances
    utility _util;

    // Variables
    int run{0}, subrun{0}, event{0};
    
    TTree * tree;

    std::string run_period;

    std::vector<double> POT_v; // vector of POT for each variation 


    // -------------------------------------------------------------------------
    // Initialiser function
    void Initialise(const char *run_period, utility _utility);
    // -------------------------------------------------------------------------
    // Function to loop over events and calculate the cross section
    void MakeHistograms(); 
    // -------------------------------------------------------------------------
    void GetPOT(const char* run_period);
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------

    enum enum_variations {
        k_CV,
        k_bnb_diffusion,
        k_vars_MAX
    };

    std::string var_string = {
        "CV",
        "BNB_Diffusion"
    };

    std::string var_string_pretty = {
        "CV",
        "BNB Diffusion"
    };


}; // End Class SystematicsHelper

#endif