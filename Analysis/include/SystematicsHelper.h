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
    // Plots the variation comparisons
    void PlotVariations(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name);
    // -------------------------------------------------------------------------
    // Plots the variation comparisons
    void PlotVariationsEXT(std::string hist_name, const char* print_name, std::string cut_name, const char* x_axis_name);
    // -------------------------------------------------------------------------
    void SetVariationProperties(TH1D* h, int index);
    // -------------------------------------------------------------------------
    void SetTPadOptions(TPad * topPad, TPad * bottomPad );
    // -------------------------------------------------------------------------
    void CreateDirectory(std::string folder, std::string run_period);
    // -------------------------------------------------------------------------

    std::string mode{"default"}; // what mode to run this class in


    enum enum_variations {
        k_CV,
        k_bnb_diffusion,
        k_vars_MAX
    };

    enum enum_ext {
        k_BNB,
        k_NuMI,
        k_ext_MAX
    };

    std::vector<std::string> var_string = {
        "CV",
        "BNB_Diffusion"
    };

    std::vector<std::string> var_string_pretty = {
        "CV",
        "BNB Diffusion"
    };


}; // End Class SystematicsHelper

#endif